### 0. PACKAGES AND INPUT

options(repos = c(RSPM = "https://packagemanager.posit.co/cran/latest"))
packages <- c("dplyr","readr", "sf", "terra","leaflet", "randomForest", "blockCV", "predicts", "devtools", "htmlwidgets")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    pak::pkg_install(pkg, build = FALSE)
    library(pkg, character.only = TRUE)
  }
}

# select an example species for demonstration
target_species <- read_lines("target_species.txt")

### 1. DOWNLOAD DATA

# download raster data
system("pip install gdown", intern = TRUE)

shared_link <- "https://drive.google.com/file/d/15s8ApYuYDP5MK-lP0aTKayFcdHuX2Waz/view?usp=sharing"

file_id <- sub(".*?/d/(.*?)/.*", "\\1", shared_link)
download_cmd <- sprintf("gdown https://drive.google.com/uc?id=%s -O raster_data.zip", file_id)

system(download_cmd)

unzip("raster_data.zip", exdir = "content/")

source_url("https://raw.githubusercontent.com/Biodiversity-Meets-Data/SDMExample/refs/heads/master/bmd_functions.R")

# load bat occurrence data and select required columns
presence_full <- read_delim("https://raw.githubusercontent.com/Biodiversity-Meets-Data/SDMExample/refs/heads/master/bats_europe.csv", delim = "\t") %>%
  select(species, decimalLongitude, decimalLatitude)

# load raster data

landmask <- rast("content/raster_data/landseamask_agg10.tif")
landmask_europe <- rast("content/raster_data/landseamask_europe_agg10.tif")

# load raster datasets containing the 19 bioclimatic variables from CHELSA
climate_current <- rast(list.files("content/raster_data/1981_2010/", "_agg10.tif", full.names = TRUE))
names(climate_current) <- str_extract(names(climate_current), "bio[0-9]{1,2}")
climate_future <- rast(list.files("content/raster_data/2071_2100//", "_agg10.tif", full.names = TRUE))
names(climate_future) <- str_extract(names(climate_future), "bio[0-9]{1,2}")


# convert to spatial sf object
presence_full_sf <- st_as_sf(presence_full, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# which species are present in the dataset?
all_species <- unique(presence_full$species)

### 2. SAMPLE BACKGROUND
presence_sf <- presence_full_sf %>%
  filter(species == target_species) %>%
  remove_duplicates_per_cell(landmask)

presence_tg_sf <- presence_full_sf %>%
  filter(species != target_species) %>%
  remove_duplicates_per_cell(landmask)

n_samples <- ifelse(nrow(presence_sf)*5 > nrow(presence_tg_sf), nrow(presence_tg_sf), nrow(presence_sf)*5)

background_points <- presence_tg_sf %>% 
  sample_n(n_samples) %>%
  dplyr::select(-species)

pb <- bind_rows(
  presence_sf %>% dplyr::select() %>% mutate(occ = 1),
  st_sf(background_points) %>% mutate(occ = 0)
)

### 3. EXTRACT PREDICTORS
pb_ext <- terra::extract(climate_current, vect(pb), ID = FALSE)

model_df <- pb %>%
  bind_cols(pb_ext) %>%
  mutate(occ = as.factor(occ)) %>%
  bind_cols(st_coordinates(.) %>% setNames(c("X", "Y"))) %>%
  st_set_geometry(NULL) %>%
  filter(complete.cases(.))

model_sf <- st_as_sf(model_df, coords = c("X", "Y"), crs = 4326)

### 4. SELECT VARIABLES
selected_variables <- select_uncorrelated_features(model_df, 
                                                   variable_names = setdiff(names(climate_current), "bio14"), 
                                                   univariate_metric = "aic", cor_cut = 0.7)
model_sf_select <- model_sf %>% dplyr::select(all_of(c("occ", selected_variables)))

model_df_folds <- folds_spatial(model_sf_select, k = 5) 

### 5. TUNE HYPERPARAMETERS
tune_grid <- tibble(nodesize = round(seq(1, nrow(model_df_folds)*0.05, length.out = 7)))

tunings_results_rf <- tune_model(fit_rf, 
                                 predictors = selected_variables, 
                                 pred_fun = pred_rf,
                                 model_df = model_df_folds, 
                                 tune_grid = tune_grid,
                                 ncores = 7)

selected_run <- which.max(tunings_results_rf$auc)
nodesize <- tune_grid$nodesize[selected_run]


### 6. FIT MODEL
m_rf <- fit_rf(model_sf_select, selected_variables, ntree = 500, nodesize = nodesize)


### 7. PREDICT MODEL
predictors_europe_current <- crop(mask(subset(climate_current, selected_variables), landmask), landmask_europe)
predictors_europe_future <- crop(mask(subset(climate_future, selected_variables), landmask), landmask_europe)

ref_df <- model_df %>% dplyr::select(all_of(selected_variables))

extrapolation <- extrapolation_exdet(predictors_europe_future, ref_df)
extrapolation_nt1 <- extrapolation[["nt1"]] < 0

pred_current_rf <- terra::predict(predictors_europe_current, model = m_rf, fun = pred_rf)

pred_future_rf <- terra::predict(predictors_europe_future, model = m_rf, fun = pred_rf)

pred_diff_rf <- pred_future_rf - pred_current_rf

presence_records <- model_sf %>% filter(occ == 1)


### 8. MAP VISUALIZATION

pal_diff <- colorNumeric(palette = "RdBu", domain = values(pred_diff_rf), reverse = TRUE)
pal_default <- colorNumeric(palette = "viridis", domain = c(0,1))

m <- leaflet() %>%
  addTiles() %>%
  addRasterImage(extrapolation_nt1, colors = c("transparent", "red"), opacity = 1,
                 group = "Uncertain predictions (extrapolation)") %>%
  addRasterImage(pred_current_rf, colors = pal_default, opacity = 1,
                 group = "Climatic suitability: current (1981-2010)") %>%
  addRasterImage(pred_future_rf, colors = pal_default, opacity = 1,
                 group = "Climatic suitability: future (2071-2100)") %>%
  addRasterImage(pred_diff_rf, colors = pal_diff, opacity = 1,
                 group = "Change in climatic suitability") %>%
  addCircleMarkers(
    color = "black", 
    data = presence_records,
    stroke = FALSE,
    fillOpacity = 0.7,
    radius = 3,
    group = "Occurrence records"
  ) %>%
  addLayersControl(
    overlayGroups = c("Uncertain predictions (extrapolation)",
                      "Climatic suitability: current (1981-2010)", 
                      "Climatic suitability: future (2071-2100)",
                      "Change in climatic suitability",
                      "Occurrence records"),
    options = layersControlOptions(collapsed = FALSE)
  )


### 9. EXPORT OUTPUTS

if(!dir.exists("outputs")){
  dir.create("outputs")
}
saveWidget(m, file="outputs/interactive_map.html")
mapshot(m, url = "outputs/m.html")

raster_outputs <- c(pred_current_rf, pred_future_rf, pred_diff_rf, as.numeric(extrapolation_nt1)) %>%
  setNames(c("current", "future", "difference", "extrapolation"))
writeRaster(raster_outputs, filename = "outputs/model_predictions.tif", overwrite = TRUE)
