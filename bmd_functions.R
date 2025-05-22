
sample_points_raster <- function(r, n, sample_factor = 10){

  cells_out <- vector(mode = "numeric", length = 0)
  cells_pool <- 1:ncell(r)

  while(length(cells_out) < n){

    sample_cells <- sample(cells_pool, n*sample_factor)
    probs <- extract(r, sample_cells)

    s <- tibble(cellnr = sample_cells,
                probs = probs[,1]) %>%
      filter(!is.na(probs)) %>%
      mutate(keep = rbernoulli(1, p = probs))


    batch <- s %>% filter(keep == TRUE) %>% pull(cellnr)

    cells_out <- unique(c(cells_out, batch))

  }


  cells_final <- sample(cells_out, n)

  cells_xy <- as.data.frame(xyFromCell(r, cells_final))

  out <- st_as_sf(cells_xy, crs = st_crs(r), coords = c("x", "y"))

  out
}


biasgrid_tg_extent <- function(x, r){

  counts <- terra::rasterizeGeom(vect(x), r, "count")
  biasgrid <- as.numeric(counts > 0)
  biasgrid <- terra::classify(biasgrid, matrix(c(0, NA), ncol = 2, byrow = TRUE))
  setNames(biasgrid, "biasgrid")
}


biasgrid_tg_count <- function(x, r, min_prob = 0, max_prob = 1, log_transform = TRUE){

  counts <- terra::rasterizeGeom(vect(x), r, "count")
  counts <- terra::mask(counts, r)
  if(log_transform){
    counts <- log(counts)
  }
  counts_reclass <- terra::classify(counts, matrix(c(-Inf, 0), ncol = 2, byrow = TRUE))
  biasgrid <- terra::stretch(counts_reclass, min_prob, max_prob)

  setNames(biasgrid, "biasgrid")
}


biasgrid_tg_kd <- function(x, r, bandwidth, agg_factor = 10,
                            n_cores = 1, min_prob = 0, max_prob = 1){

  if(agg_factor > 1){
    r_agg <- terra::aggregate(r, agg_factor)
  } else {
    r_agg <- r
  }

  # Estimate kernels
  xrange <- c(terra::xmin(r_agg), terra::xmax(r_agg))
  yrange <- c(terra::ymin(r_agg), terra::ymax(r_agg))
  rncol <- terra::ncol(r_agg)
  rnrow <- terra::nrow(r_agg)


  # Create Raster
  kde <- KernSmooth::bkde2D(st_coordinates(x), bandwidth = bandwidth,
                            range.x = list(xrange, yrange),
                            gridsize = c(rncol, rnrow))

  # Finish output
  out <- r_agg
  out[] <- scales::rescale(as.vector(kde$fhat[, ncol(kde$fhat):1]), to = c(min_prob, max_prob))

  out_p <- terra::project(out, r, method = "bilinear", threads = n_cores, use_gdal = TRUE)
  biasgrid <-  terra::mask(out_p, r)

  setNames(biasgrid, "biasgrid")
}


biasgrid_thickening <- function(x, r, min_prob = 0, max_prob = 1, quadsegs = 10){

    x_vect <- vect(x)
    buffer_size <- mean(terra::distance(x_vect))
    x_buffer <- terra::buffer(x_vect, buffer_size, quadsegs = quadsegs)
    buffer_count <- rast(fasterize::fasterize(sf::st_as_sf(x_buffer), raster::raster(r), fun = "sum"))
    r_zero <- terra::classify(r, matrix(c(0, Inf, 0), ncol = 3, byrow = TRUE))
    buffer_count_cover <- terra::cover(buffer_count, r_zero)
    biasgrid <- terra::stretch(buffer_count_cover, min_prob, max_prob)
    biasgrid <- terra::mask(biasgrid, r)

    setNames(biasgrid, "biasgrid")


}

thin_presence <- function(sf, thin_dist = 5000, runs = 10, ncores = 10){

  sample.vec <- function(x, ...) x[sample(length(x), ...)]

  sf_buffer <- sf::st_buffer(sf, thin_dist)
  buff_int <- sf::st_intersects(sf, sf_buffer)
  buff_int <- setNames(buff_int, 1:length(buff_int))

  n_int <- purrr::map_dbl(buff_int, length)

  furrr::plan(multisession, workers = ncores)

  seeds <- sample.int(n = runs)
  results_runs <- furrr::future_map(seeds, function(i){

    set.seed(i)
    while(max(n_int) > 1){
      max_neighbors <- names(which(n_int == max(n_int)))

      # remove point with max neighbors
      sampled_id <- sample.vec(max_neighbors, 1)

      purrr::pluck(buff_int, sampled_id) <- NULL
      buff_int <- purrr::map(buff_int, function(x) setdiff(x, as.numeric(sampled_id)))
      n_int <- purrr::map_dbl(buff_int, length)
    }

    unlist(buff_int) %>% unique()

  })

  lengths <- purrr::map_dbl(results_runs, length)

  selected_run <- results_runs[[sample.vec(which(lengths == max(lengths)), 1)]]

  out <- sf[selected_run,]

  out

}

remove_duplicates_per_cell <- function(sf, r){


  cellnrs <- terra::cells(r, vect(sf))[,"cell"]

  sf %>%
    mutate(cell = cellnrs) %>%
    group_by(cell) %>%
    sample_n(1) %>%
    dplyr::select(-cell)

}




fit_rf <- function(df, vars, ntree = 500, ...){

  prNum <- sum(df[["occ"]] == 1) # number of presence records
  spsize <- c("0" = prNum, "1" = prNum) # sample size for both classes

  f <- as.formula(paste0("occ~", paste(vars, collapse = "+")))
  m <- randomForest::randomForest(f,
                                  data = df,
                                  ntree = ntree,
                                  sampsize = spsize,
                                  replace = TRUE,
                                  ...)

  m

}


pred_rf <- function(model, newdata){

  require(randomForest, quietly = TRUE)
  predict(model, newdata, type = "prob")[,2]

}


fit_glm <- function(df, vars, ...){

  f <- as.formula(paste0("occ~", paste(vars, collapse = "+")))
  m <- glm(f,
          family = "binomial",
          data = df)

  m

}


pred_glm <- function(model, newdata){

  predict(model, newdata, type = "response")

}


fit_maxent <- function(df, vars, args =NULL,  features = c("hinge", "linear", "quadratic", "product"), beta = NA, clamp = TRUE){

  beta <- ifelse(is.na(beta), 1, beta)

  if(is.null(args)){

    h_par <- ifelse("hinge" %in% features, "true", "false")
    t_par <- ifelse("threshold" %in% features, "true", "false")
    l_par <- ifelse("linear" %in% features, "true", "false")
    q_par <- ifelse("quadratic" %in% features, "true", "false")
    p_par <- ifelse("product" %in% features, "true", "false")

    args <- c(
      paste0('linear=', l_par),
      paste0('quadratic=', q_par),
      paste0('product=', p_par),
      paste0('hinge=', h_par),
      paste0('threshold=', t_par),
      paste0('betamultiplier=', beta))

  } else {

    args <- args

  }

  if(!clamp){
    args <- c(args,  paste0('doclamp=false'))
  }


  predicts::maxent(df %>% dplyr::select(all_of(vars)),
                    df %>% pull(occ),
                    args = args)

}


pred_maxent <- function(model, newdata, args){

  out <- predict(model, as.data.frame(newdata), args = args)

  out

}

select_uncorrelated_features <- function(model_df, variable_names, univariate_metric = "p", cor_cut = 0.7){

  # step 1: univariate GLMs with performance metric
  univariate_performance <- map_dbl(variable_names, function(var_i){
    f <- as.formula(paste0("occ~", var_i))
    m <- glm(f, family = "binomial", data = model_df)

    if(univariate_metric == "p"){
      metric <- summary(m)$coefficients[2,4]
    } else if(univariate_metric == "aic"){
      metric <- AIC(m)
    }

    metric
  }) %>% setNames(variable_names) %>% sort()

  # step 2: remove highly correlated pairs based on univariate performance
  cor_mat <- cor(model_df %>% dplyr::select(all_of(variable_names)))

  cor_mat_long <-cor_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(-rowname) %>%
    setNames(c("var1", "var2", "r")) %>%
    filter(var1 != var2)

  highly_correlated <- cor_mat_long %>%
    filter(r > cor_cut)

  eliminated_vars <- map_chr(1:nrow(highly_correlated), function(i){

    var_pair <- unlist(highly_correlated[i,1:2])
    performances <- univariate_performance[var_pair]

    names(which.max(performances))

  }) %>% unique()

  setdiff(variable_names, eliminated_vars)

}



folds_spatial <- function(x, ...){

  folds <- blockCV::cv_spatial(x, column = "occ", progress = FALSE, report = FALSE, ...)

  x <- x %>%
    mutate(fold = folds$folds_ids)

  complete_folds <- x %>%
    group_by(fold) %>%
    summarize(l = length(unique(occ))) %>%
    filter(l == 2) %>%
    pull(fold)

  x %>%
    filter(fold %in% complete_folds)

}

folds_random <- function(x,
                          nfolds = 5){
  x %>%
    mutate(fold = sample(1:nfolds, nrow(.), TRUE))

}


tune_model <- function(fitfun, predictors, pred_fun, model_df, tune_grid, ncores = 1,  metrics = c("auc", "cor")){

  folds <- unique(model_df[["fold"]])

  ncores <- ifelse(ncores > nrow(tune_grid), nrow(tune_grid), ncores)

  future::plan(future::multisession, workers = ncores)

  tune_grid_split <- split(tune_grid, 1:nrow(tune_grid))
  tuning_results <- furrr::future_imap_dfr(tune_grid_split, function(pars_i, i){

    cv_results <- map_dfr(folds, function(fold_i){

      train_dat_i <- model_df %>% filter(fold != fold_i)
      test_dat_i <- model_df %>% filter(fold == fold_i)

      fit_args <- pars_i %>%
        as.list()

      fit_args$df <- train_dat_i
      fit_args$vars <- predictors

      m <- do.call(fitfun, fit_args)

      preds_i <- pred_fun(m, test_dat_i)

      eval <- predicts::pa_evaluate(preds_i[test_dat_i$occ == 1],
                            preds_i[test_dat_i$occ == 0])


      tibble(fold = fold_i) %>%
        bind_cols(eval@stats[metrics]) %>%
        bind_cols(pars_i) %>%
        mutate(iteration = as.numeric(i))

    })

    cv_results


  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))


  tuning_results %>%
    group_by_at(c("iteration", colnames(tune_grid))) %>%
    summarize_at(metrics, mean) %>%
    ungroup()


}


attribution_df <- function(m, feature_names, ref, newdata, pred_fun, nsim = 1, ncores = 1){


  ncores <- ifelse(ncores > length(feature_names), length(feature_names), ncores)
  future::plan(future::multisession, workers = ncores)

  exp <- furrr::future_map_dfc(feature_names, function(feature_i){

    exp <- fastshap::explain(
            object = m,
            feature_names = feature_i,
            X = as.data.frame(ref),
            nsim = nsim,
            pred_wrapper = pred_fun,
            newdata = as.data.frame(newdata),
            adjust = FALSE,
            shap_only = TRUE,
            parallel = FALSE)

  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))


  exp %>%
    as_tibble() %>%
    summarize_all(mean, na.rm = TRUE)


}

explain_pixel <- function(v, m, feature_names, pred_fun, nsim = 1){

  nvars <- length(v)/2


  if(!all(is.na(v))){

    x <- v[1:nvars]
    names(x) <- feature_names
    y <- v[(nvars+1):length(v)]
    names(y) <- feature_names

    exp <- fastshap::explain(
      object = m,
      X = t(as.matrix(x)),
      feature_names = feature_names,
      nsim = nsim,
      pred_wrapper = pred_fun,
      newdata = t(as.matrix(y)),
      shap_only = TRUE)

    out <- as.numeric(exp)

  } else {

    out <- rep(NA, nvars)

  }

  names(out) <- feature_names
  out


}


attribution_raster <- function(m, ref, newdata, feature_names = NULL, pred_fun, nsim = 1, ncores = 1){

  if(is.null(feature_names)){
    feature_names <- names(ref)
  }

  s <- c(ref, newdata)
  feature_names <- names(ref)

  if(ncores > 1){

    cl <- parallel::makeCluster(ncores)
    clusterEvalQ(cl, library(fastshap))
    parallel::clusterExport(cl, c("m", "pred_fun", "nsim", "feature_names"), envir=environment())
    out <- app(s, explain_pixel, m =m , feature_names = feature_names, pred_fun = pred_fun, nsim = nsim, cores=cl)
    stopCluster(cl)

  } else {

    out <- app(s, explain_pixel, m =m , feature_names = feature_names, pred_fun = pred_fun, nsim = nsim)

  }


  out

}


exdet_scores <- function (ref, p, mic = FALSE) {

  ud <- purrr::map2_dfr(ref, p, function(ref_i, p_i) {
    rng <- range(ref_i)
    intervals <- findInterval(p_i, rng)

    ifelse(intervals == 0, (p_i-rng[1])/diff(rng),
           ifelse(intervals==1, 0, (rng[2]-p_i)/diff(rng)))

  })
  nt1 <- d <- rowSums(ud) # type 1 novelty

  a <- colMeans(ref)
  b <- var(ref)
  d_ref <- stats::mahalanobis(as.matrix(ref), a, b, tol=1e-20)
  d_p <- stats::mahalanobis(as.matrix(p), a, b, tol=1e-20)

  d_max <- max(d_ref[is.finite(d_ref)])
  nt2 <- d_p/d_max # type 2 novelty
  e <- ifelse(nt1==0, nt2, nt1) # combine into exdet score

  out <-data.frame(score=e, nt1 = nt1, nt2 =nt2)

  if(mic) {
    out <- data.frame(score=e, nt1 = nt1, nt2 =nt2,MIC1=NA, MIC2=NA)

    # NT1
    if(any(nt1 != 0, na.rm = TRUE)){

      out$MIC1[which(nt1 != 0)] <- apply(ud[which(nt1 != 0), ], 1, which.min)


    } else {

      out$MIC1 <- NA

    }

    if(any(nt2 > 1, na.rm = TRUE)){

      # NT2
      d_v <- as.data.frame(lapply(seq_along(ref), function(i) {
        a <- colMeans(ref[, -i])
        b <- var(ref[, -i])
        d_ref_v <- stats::mahalanobis(x=as.matrix(ref)[, -i], a, b, tol = 1e-20) # reference distance with ith variable dropped
        d_p_v <-  stats::mahalanobis(x=as.matrix(p)[, -i], a, b, tol = 1e-20) # projection distance with ith variable dropped
        (d_p - d_p_v)/d_p
      }))
      names(d_v) <- names(ref)
      out$MIC2[which(nt2 > 1)] <- apply(d_v[which(nt2 > 1), ], 1, which.max)

    } else {

      out$MIC2 <- NA

    }


  }


  as.matrix(out)

}




extrapolation_exdet <- function(x, ref , ncores = 1, mic = FALSE){

  if(ncores > 1){

    if(ncores %% 2 == 1){
      ncores <- ncores-1
    }

    tiler <- rast(ncol=2,nrow=ncores/2,extent=ext(x))
    tiles <- makeTiles(x,tiler,overwrite=TRUE)

    future::plan(future::multisession, workers = ncores)
    results <- furrr::future_map(tiles, function(tile_i){

      p_i <- as.data.frame(values(rast(tile_i)))
      values_out_i <- exdet_scores(ref = ref, p = p_i, mic = mic)
      out_i <- rast(rast(tile_i), nlyrs = ncol(values_out_i), vals = FALSE)
      values(out_i) <- values_out_i
      wrap(out_i)


    }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))

    tiles_unwrap <- map(results, unwrap)
    out <- do.call(terra::merge,tiles_unwrap)
    unlink(tiles)

    out

  } else {

    p_i <- as.data.frame(values(x))
    values_out_i <- exdet_scores(ref = ref, p = p_i, mic = mic)
    out <- rast(x, nlyrs = ncol(values_out_i), vals = FALSE)
    values(out) <- values_out_i

  }

  names(out) <- c("score", "nt1", "nt2")
  out


}

download_occ_cube <- function(taxon_ids){


  taxon_ids_string <- paste0(taxon_ids, collapse =',')

  sql_string <- paste0("
  SELECT
    species,
    decimalLatitude,
    decimalLongitude
  FROM occurrence
  WHERE taxonKey IN (",taxon_ids_string,")
    AND \"year\" >= 2010
    AND hasCoordinate = TRUE
    AND (coordinateUncertaintyInMeters <= 1000 OR coordinateUncertaintyInMeters IS NULL)
    AND NOT ARRAY_CONTAINS(issue, 'ZERO_COORDINATE')
    AND NOT ARRAY_CONTAINS(issue, 'COORDINATE_OUT_OF_RANGE')
    AND NOT ARRAY_CONTAINS(issue, 'COORDINATE_INVALID')
    AND NOT ARRAY_CONTAINS(issue, 'COUNTRY_COORDINATE_MISMATCH')")

  dl <- occ_download_sql(sql_string)
  occ_download_wait(dl)
  dat <- occ_download_get(dl) %>% occ_download_import()

}


download_occ <- function(taxon_ids){

  dl <- occ_download(
    pred_in("taxonKey", taxon_ids),
    pred("occurrenceStatus","PRESENT"),
    pred_gte("year", 1980),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    pred_or(
      pred_lt("coordinateUncertaintyInMeters",1000),
      pred_isnull("coordinateUncertaintyInMeters")
    ),
    pred_in("basisOfRecord",c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN"))
  )

  occ_download_wait(dl)
  dat <- occ_download_get(dl) %>% occ_download_import()
}
