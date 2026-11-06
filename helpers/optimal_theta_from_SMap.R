optimal_theta_from_SMap <- function(ts_data, colname = "ts", drop_na = TRUE) {
  
  # --- extract a single series, stripping any names/dimnames ---
  if (is.data.frame(ts_data) || is.matrix(ts_data)) {
    if (ncol(ts_data) != 1) stop("Expected a single time series (one column).")
    x <- as.numeric(ts_data[, 1])
  } else {
    x <- as.numeric(ts_data)          # coercion drops names + dimnames
  }
  if (drop_na) x <- x[!is.na(x)]
  N <- length(x)
  
  # --- attach the name that 'columns'/'target' will reference ---
  df <- data.frame(x)
  names(df) <- colname
  
  # --- train/test split 70/30 ---
  mid  <- floor(N * 0.7)
  lib  <- c(1, mid)
  pred <- c(mid + 1, N)
  
  # --- 1) optimal embedding dimension ---
  e <- tryCatch({
    emb <- suppressWarnings(EmbedDimension(
      dataFrame = df, lib = lib, pred = pred,
      maxE = 15, Tp = 1,
      columns = colname, target = colname,
      noTime = TRUE, showPlot = FALSE
    ))
    emb$E[which.max(emb$rho)]
  }, error = function(e) NA_real_)
  
  if (is.na(e)) return(NA_real_)
  
  # --- 2) S-map rho over theta, take argmax ---
  thetas <- seq(0, 8, by = 0.5)
  rho_theta <- vapply(thetas, function(th) {
    sm <- tryCatch(
      suppressWarnings(SMap(
        dataFrame = df, lib = lib, pred = pred,
        columns = colname, target = colname,
        E = e, Tp = 1, theta = th, noTime = TRUE
      )),
      error = function(e) NULL
    )
    if (is.null(sm)) return(NA_real_)
    obs <- sm$predictions$Observations
    prd <- sm$predictions$Predictions
    ok  <- is.finite(obs) & is.finite(prd)
    if (sum(ok) < 3) return(NA_real_)
    cor(obs[ok], prd[ok])
  }, numeric(1))
  
  if (all(is.na(rho_theta))) return(NA_real_)
  
  thetas[which.max(rho_theta)]
}
