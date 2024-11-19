GenSamplingPoints <- function(center, radius, numpts, dist.metric = "Euclidean", cType){
  if (cType == "edge"){
    if (dist.metric == "Euclidean"){
      return(cbind(center[1] + cos(((1:numpts)/numpts)*2*pi) * radius,
                   center[2] + sin(((1:numpts)/numpts)*2*pi) * radius))
    }else if (dist.metric == "Geodesic"){
      return(destPoint(center, b = ((1:numpts)/numpts)*360, d = radius)) 
    }
  }else if (cType == "disk" | cType == "donut"){
    if (dist.metric == "Euclidean"){
      radius <- seq(0, radius, dist_unit)
      all_points <- lapply(radius, function(x) cbind(center[1] + cos(((1:numpts)/numpts)*2*pi) * x, 
                                                     center[2] + sin(((1:numpts)/numpts)*2*pi) * x))
      return(do.call("rbind", all_points))
    }else if (dist.metric == "Geodesic"){
      all_points <- lapply(radius, function(x) destPoint(center, b = ((1:numpts)/numpts)*360, d = x))
      return(do.call("rbind", all_points))
    }
  }
}

GenSamplingPointsBuffer <- function(center, ras, radius, cType){
  one_buffer <- buffer(center, width=radius)
  if (cType == "edge"){
    one_buffer_raster <- boundaries(rasterize(one_buffer, ras))
  }else if (cType == "disk" | cType == "donut"){
    one_buffer_raster <- rasterize(one_buffer, ras)
  }
  one_buffer_raster_values <- getValues(one_buffer_raster)
  coords <- coordinates(one_buffer_raster)[(!is.na(one_buffer_raster_values) & 
                                            !is.nan(one_buffer_raster_values) & 
                                            one_buffer_raster_values == 1), ]
  rm(one_buffer)
  rm(one_buffer_raster)
  rm(one_buffer_raster_values)
  return(coords)
}

RimAvg <- function(ras, Yobs = NULL, ras_Z = NULL, nz, Zup, Zdata, x_coord_Z, y_coord_Z, 
                   treatment, dUp, numpts = NULL, gridRes = NULL, evalpts = 10, 
                   only.unique = 0, dtype = "raster", dist.metric = dist.metric, 
                   cType = cType){

  Ybard <- Ybard_sum <- Ybard_len <- rep(NA, nz)
  if (dtype == "krig"){
    if (is.null(numpts)){
      numpts <- ceiling(ceiling(dUp/(gridRes/sqrt(2))+1)*pi)*evalpts
    }
    for (i in 1:nz) {
      coords <- GenSamplingPoints(center = c(Zdata[i, x_coord_Z], Zdata[i, y_coord_Z]), 
                                  radius = dUp, numpts = numpts, dist.metric = dist.metric,
                                  cType = cType)
      grids <- predict(ras, coords)
      Ybard[i] <- mean(grids)
      Ybard_sum[i] <- sum(grids, na.rm = 1)
      Ybard_len[i] <- length(grids[!is.na(grids)])
    }
  }else{
    if (is.null(numpts)){
      if (is.null(gridRes)){
        gridRes <- min(res(ras))
      }
      numpts <- ceiling(ceiling(dUp/(gridRes/sqrt(2))+1)*pi)*evalpts
    }
    for (i in 1:nz){
      if (!is.null(ras_Z)){
        coords <- GenSamplingPointsBuffer(center = as_Spatial(ras_Z[i, ]), ras = ras, radius = dUp,
                                          cType = cType)
      }else{
        coords <- GenSamplingPoints(center = c(Zdata[i, x_coord_Z], Zdata[i, y_coord_Z]), 
                                    radius = dUp, numpts = numpts, dist.metric = dist.metric,
                                    cType = cType)
      }
      if (only.unique == 1){
        grid_index <- c(unique(na.omit(cellFromXY(ras, coords))))
      }else{
        grid_index <- c(na.omit(cellFromXY(ras, coords)))
      }
      grids <- Yobs[grid_index]
      Ybard[i] <- mean(grids, na.rm = TRUE)
      Ybard_sum[i] <- sum(grids, na.rm = TRUE)
      Ybard_len[i] <- length(grids[!is.na(grids)])
      rm(grid_index)
      rm(coords)
    }
  }
  Rim.list <- list(Ybard, Ybard_sum, Ybard_len)
  return(Rim.list)
}

# these two functions are based on Hainmueller et al. 2018
LocalReg <- function(dVec, Sdata, bw, bw_debias = NULL, Zup = NULL, xevals = NULL, smooth.conley.se = 1, 
                     kernel = "uni", cutoff = 0, dist = NULL, dist.metric = "Euclidean", bias_correction = TRUE){
  
  if (is.null(xevals)){
    xevals <- dVec
  }
  # xevals_std <- (xevals - mean(xevals)) / sd(xevals)
  all_names <- names(Sdata)
  wls_formula <- as.formula("outcome ~ treatment + dVec_demeaned + interaction")
  wls_formula_debias <- as.formula("outcome ~ treatment + dVec_demeaned + interaction + dVec_demeaned_2 + interaction_2")
  num_coefs <- 4
  num_debias_coefs <- 6
  if (length(all_names) > 4){
    cov_names <- all_names[5:length(all_names)]
    cov_formula <- paste0(paste0("treatment * ", cov_names), collapse = "+")
    wls_formula <- as.formula(paste0("outcome ~ treatment + dVec_demeaned + interaction", "+", cov_formula))
    wls_formula_debias <- as.formula(paste0("outcome ~ treatment + dVec_demeaned + interaction + dVec_demeaned_2 + interaction_2", "+", cov_formula))
    num_coefs <- num_coefs + 2*length(cov_names)
    num_debias_coefs <- num_debias_coefs + 2*length(cov_names)
  }
  wls_results <- list()
  wls_coefs <- wls_ses <- rep(NA, length(xevals))
  wls_coefs_all <- matrix(NA, length(xevals), num_coefs)
  wls_coefs_debias_all <- matrix(NA, length(xevals), num_debias_coefs)
  
  if (smooth.conley.se){
    if (dist.metric == "Euclidean"){
      metric <- 1
    }else if (dist.metric == "Geodesic"){
      metric <- 2
    }
    x_coord <- Zdata[, x_coord_Z]
    y_coord <- Zdata[, y_coord_Z]
    if (is.null(dist)){
      dist <- DistanceCalculation(as.vector(x_coord), as.vector(y_coord), as.integer(metric))[["Dist_mat"]]
    }
    # dist_rep <- diag(length(xevals)) %x% dist
    dist_rep <- matrix(rep(t(dist), length(xevals)), ncol=ncol(dist), byrow=TRUE)
    dist_rep <- matrix(rep(dist_rep, length(xevals)), nrow=nrow(dist_rep), byrow=FALSE)
    xevals_adj_mat <- sapply(xevals, function(x) abs(x - xevals))
    unit_mat <- diag(2, nrow(Zdata))
    unit_mat_other <- matrix(-1, nrow = nrow(Zdata), ncol = nrow(Zdata))
    unit_mat <- unit_mat + unit_mat_other
    unit_mat_rep <- xevals_adj_mat %x% unit_mat
    dist_rep <- dist_rep + unit_mat_rep
  }
  
  for (i in 1:length(xevals)){
    xeval <- xevals[i]
    dat1 <- Sdata
    # dat1$dVec <- (dat1$dVec - mean(dat1$dVec)) / sd(dat1$dVec)
    dat1$dVec_demeaned <- dat1$dVec - xeval
    dat1$interaction <- dat1$treatment * dat1$dVec_demeaned
    if (bw == 0){
      dat1$w <- as.numeric(dat1$dVec_demeaned == 0) * dat1$pweight
    }else{
      # dat1$w <- dnorm(dat1$dVec_demeaned/bw)
      dat1$w <- (1-abs(dat1$dVec_demeaned/bw)) * (abs(dat1$dVec_demeaned/bw) <= 1) * dat1$pweight
      # dat1$w <- 0.75 * (1-abs(dat1$dVec_demeaned/bw)^2) * as.numeric(abs(dat1$dVec_demeaned/bw) <= 1)
      # dat1$w <- abs(dat1$dVec_demeaned/bw <= 1)
      w_index <- (dat1$w > 0)
    }
    if (bias_correction){
      # dat1$w_debias <- dnorm(dat1$dVec_demeaned/bw_debias)
      dat1$w_debias <- (1-abs(dat1$dVec_demeaned/bw_debias)) * (abs(dat1$dVec_demeaned/bw_debias) <= 1) * dat1$pweight
      # dat1$w_debias <- abs(dat1$dVec_demeaned/bw_debias <= 1)
      # dat1$w_debias <- 0.75 * (1-abs(dat1$dVec_demeaned/bw_debias)^2) * as.numeric(abs(dat1$dVec_demeaned/bw_debias) <= 1)
      w_index <- (dat1$w_debias > 0)
      dat1$dVec_demeaned_2 <- (dat1$dVec_demeaned)^2
      dat1$interaction_2 <- (dat1$interaction)^2
      wls_reg_debias <- lm(wls_formula_debias, data = dat1, weights = w_debias)
      wls_coefs_debias_all[i, ] <- coef(wls_reg_debias)
    }
    wls_reg <- lm(wls_formula, data = dat1, weights = w)
    wls_coef <- coef(wls_reg)
    wls_coef[which(is.na(wls_coef))] <- 0 
    wls_coefs[i] <- wls_coef[2]
    wls_coefs_all[i, ] <- wls_coef
    
    W_K <- diag(sqrt(dat1$w[w_index]))
    res <- W_K %*% wls_reg$residuals[w_index]
    # X_mat <- W_K %*% as.matrix(cbind(rep(1, length(res)), wls_reg$model[, c("treatment", "dVec_demeaned", "interaction")]))
    X_mat <- W_K %*% as.matrix(model.matrix(wls_reg))[w_index, ]
    W_meat <- X_mat
    XX_mat_inv <- solve(t(X_mat) %*% X_mat)
    
    if (bias_correction){
      W_K_debias <- diag(sqrt(dat1$w_debias[w_index]))
      res_debias <- W_K_debias %*% wls_reg_debias$residuals[w_index]
      # X_mat_debias_all <- as.matrix(cbind(rep(1, length(res)), wls_reg_debias$model[, c("treatment", "dVec_demeaned", "interaction", "dVec_demeaned_2", "interaction_2")]))
      X_mat_debias_all <- W_K_debias %*% as.matrix(model.matrix(wls_reg_debias))[w_index, ]
      X_mat_debias <- X_mat_debias_all[, c(5:6)]
      nu <- t(cbind(rep(0, num_debias_coefs), rep(0, num_debias_coefs)))
      nu[1, 5] <- nu[2, 6] <- 1
      W_meat <- t(X_mat) - 0.5*t(X_mat) %*% X_mat_debias %*% (nu %*% solve(t(X_mat_debias_all) %*% X_mat_debias_all) %*% t(X_mat_debias_all))
      W_meat <- t(W_meat)
      coef_debias <- 0.5*solve(t(X_mat) %*% X_mat) %*% (t(X_mat) %*% X_mat_debias) %*% wls_coefs_debias_all[i, 5:6]
      wls_coefs[i] <- wls_coefs[i] - coef_debias[2]
      res <- res_debias
    }
    
    if (smooth.conley.se){
      if (kernel != "uni" & kernel != "uniform" & kernel != "tri" & kernel != "triangular" & kernel != "epa" & kernel != "epanechnikov" & kernel != "" ){
        stop("kernel incorrectly specified")
      }else if (kernel == "uni" | kernel == "uniform"){
        k <- 1
      }else if (kernel == "tri" | kernel == "triangular"){
        k <- 2
      }else if (kernel == "epa" | kernel == "epanechnikov"){
        k <- 3
      }
      if (cutoff > 0){
        c <- cutoff + xeval
      }
      dist_rep_s <- dist_rep[w_index, w_index]
      
      Conley_result <- ConleySE(as.vector(res), as.matrix(W_meat), as.matrix(dist_rep_s), 
                                as.double(c), as.integer(k))
      VCE <- t(XX_mat_inv) %*% Conley_result[["VCE_meat"]] %*% XX_mat_inv
      dist_kernel <- Conley_result[["Dist_kernel"]]
      wls_ses[i] <- sqrt(VCE[2, 2])
    }
  }
  wls_coefs <- cbind(xevals, wls_coefs)
  wls_results[["coefs"]] <- wls_coefs
  wls_results[["coefs_all"]] <- wls_coefs_all
  wls_results[["ses"]] <- wls_ses
  if (bias_correction){
    wls_results[["coefs_debias"]] <- wls_coefs_debias_all
  }
  return(wls_results)   
} 

CrossValidation <- function(Sdata, outcome, treatment, dVec, grid = NULL, nfold = 5, block_cv = TRUE,
                            parallel = FALSE, metric = "MSPE", kernel = "uni", bias_correction = TRUE){
 
  ## calculate error for testing set
  getError <- function(train, test, bw, outcome, treatment, dVec, bias_correction){
    
    one_reg <- LocalReg(dVec, train, bw = bw, bw_debias = bw, xevals = NULL, smooth.conley.se = 0, 
                        kernel = kernel, cutoff = 0, dist, dist.metric = "Euclidean", bias_correction = bias_correction)
    
    if (bias_correction){
      coef <- one_reg[["coefs_debias"]][, c(1:2)]
      coef[is.na(coef)] <- 0
      esCoef <- function(x){
        X.eval <- dVec
        Xnew <- abs(dVec - x)
        d1 <- min(Xnew)     ## distance between x[i] and the first nearest x in training set
        label1 <- which.min(Xnew)
        Xnew[label1] <- Inf
        d2 <- min(Xnew)     ## distance between x[i] and the second nearest x in training set
        label2 <- which.min(Xnew)
        if(d1==0){
          func <- coef[label1, ] # X.eval (1), intercept (2), d (3), xx (4), d:xx (5), z
        }else if(d2==0){
          func <- coef[label2, ]  
        }else{ ## weighted average 
          func <- (coef[label1, ] * d2 + coef[label2, ] * d1)/(d1 + d2) 
        }
        return(func)
      } 
      Knn <- t(sapply(test$dVec, esCoef)) ## coefficients for test  class==matrix
      
      ## predicting
      test.Y <- test[, outcome]
      test.X <- as.matrix(cbind(rep(1, dim(test)[1]), test[, c("treatment")])) 
      error_debias <- test.Y - rowSums(test.X * Knn)
    }
    coef <- one_reg[["coefs_all"]][, c(1:2)]
    coef[is.na(coef)] <- 0
    esCoef <- function(x){
      X.eval <- dVec
      Xnew <- abs(dVec - x)
      d1 <- min(Xnew)     ## distance between x[i] and the first nearest x in training set
      label1 <- which.min(Xnew)
      Xnew[label1] <- Inf
      d2 <- min(Xnew)     ## distance between x[i] and the second nearest x in training set
      label2 <- which.min(Xnew)
      if(d1==0){
        func <- coef[label1, ] # X.eval (1), intercept (2), d (3), xx (4), d:xx (5), z
      }else if(d2==0){
        func <- coef[label2, ]  
      }else{ ## weighted average 
        func <- (coef[label1, ] * d2 + coef[label2, ] * d1)/(d1 + d2) 
      }
      return(func)
    } 
    Knn <- t(sapply(test$dVec, esCoef)) ## coefficients for test  class==matrix
    
    ## predicting
    test.Y <- test[, outcome]
    test.X <- as.matrix(cbind(rep(1, dim(test)[1]), test[, c("treatment")])) 
    error <- test.Y - rowSums(test.X * Knn)
    
    if (bias_correction){
      return(c(mean(abs(error)), mean(error^2), mean(abs(error_debias)), mean(error_debias^2)))
    }else{
      return(c(mean(abs(error)), mean(error^2), 0, 0))
    }
  } 
  
  ## grid search and 5 fold cross validation
  cv <- function(bw){
    mse <- matrix(NA, nfold, 4)
    for(j in 1:nfold){ # 5-fold CV
      testid <- which(fold == j)
      mse[j, ] <- getError(train = Sdata[-testid, ], test = Sdata[testid, ],
                           bw = bw, outcome = "outcome", treatment = "treatment", 
                           dVec = dVec, bias_correction = bias_correction)
    }
    return(c(bw, apply(mse, 2, mean)))
  }
  
  
  cat("Cross-validating bandwidth ... ")
  ## generate 5 random folds
  n <- dim(Sdata)[1]
  nz <- n/length(dVec)
  range_dVec <- max(dVec) - min(dVec)
  if (is.null(grid)){
    grid <- exp(seq(log(range_dVec/10), log(range_dVec), length.out = 40))
  }
  # fold <- rep(0, n)
  # fold <- c(0:(n-1))%%nfold + 1
  # fold <- sample(fold, n, replace = FALSE)
  # fold <- rep(0, nz)
  if (block_cv){
    fold <- rep(0, nz)
    x_coord_Z_vec <- c(Zdata[, x_coord_Z])
    y_coord_Z_vec <- c(Zdata[, y_coord_Z])
    nfold_sqr <- round(sqrt(nfold))
    x_coord_Z_q <- quantile(x_coord_Z_vec, seq(0, 1, 1/nfold_sqr))
    y_coord_Z_q <- quantile(y_coord_Z_vec, seq(0, 1, 1/nfold_sqr))
    for (pp in 2:length(x_coord_Z_q)){
      for (qq in 2:length(y_coord_Z_q)){
        block_ind <- (pp - 2) * (length(x_coord_Z_q) - 1) + qq - 1
        cond <- (x_coord_Z_vec >= x_coord_Z_q[pp-1] & x_coord_Z_vec <= x_coord_Z_q[pp] & 
                   y_coord_Z_vec >= y_coord_Z_q[qq-1] & y_coord_Z_vec <= y_coord_Z_q[qq])
        fold[cond] <- block_ind
      }
    }
    fold <- rep(sample(fold, nz, replace = FALSE), length(dVec))
    nfold <- nfold_sqr^2
  }else{
    fold <- c(0:(nz-1))%%nfold + 1
    fold <- rep(sample(fold, nz, replace = FALSE), length(dVec))
  }
  
  ## calculation
  if (parallel == TRUE) {
    Error <- suppressWarnings(foreach(bw = grid, .combine = rbind,
                                      .export = c("getError"),
                                      .inorder = FALSE) %dopar% {cv(bw)})
  } else {
    Error <- matrix(NA, length(grid), 5)
    for (i in 1:length(grid)) {
      Error[i, ] <- cv(grid[i])
      cat(".")
    } 
  } 
  colnames(Error) <- c("bandwidth", "MAPE", "MSPE", "MAPE_debias", "MSPE_debias")
  rownames(Error) <- NULL
  
  if (bias_correction){
    if (metric=="MAPE") {
      opt.bw <- grid[which.min(Error[, 2])]
      opt.bw.debias <- grid[which.min(Error[, 4])]
    } else {
      opt.bw <- grid[which.min(Error[, 3])]
      opt.bw.debias <- grid[which.min(Error[, 5])]
    } 
  }else{
    if (metric=="MAPE") {
      opt.bw <- grid[which.min(Error[, 2])]
      opt.bw.debias <- NULL
    } else {
      opt.bw <- grid[which.min(Error[, 3])]
      opt.bw.debias <- NULL
    } 
  }
  
  output <- list(CV.out = round(Error, 3),
                 opt.bw = opt.bw,
                 opt.bw.debias = opt.bw.debias)
  cat(paste("Bandwidth =", round(output$opt.bw, 3), "\n"))
  if (bias_correction){
    cat(paste("Bandwidth for bias correction =", round(output$opt.bw.debias, 3), "\n"))
  }
  return(output)
}