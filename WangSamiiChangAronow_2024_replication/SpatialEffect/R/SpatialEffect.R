SpatialEffect <- function(ras = NULL, Ydata = NULL, outcome = NULL, x_coord_Y = NULL, y_coord_Y = NULL, ras_Z = NULL, 
                          Zdata = NULL, x_coord_Z = NULL, y_coord_Z = NULL, treatment, covs = NULL, prob_treatment = NULL, 
                          dVec, dist.metric = "Euclidean", cType = "edge", numpts = NULL, evalpts = 1, only.unique = 0, 
                          smooth = 0, bw = NULL, bw_debias = NULL, bias_correction = TRUE, smooth.conley.se = 1, conf.band = 0, 
                          per.se = 1, blockvar = NULL, clustvar = NULL, conley.se = 1, kernel = "uni", cutoff = 0, alpha = 0.05, 
                          edf = FALSE, m = 2, lambda = 0.02, nPerms = 1000){
  
  # check the existence of treatment
  if ((is.null(Zdata) | is.null(treatment))){
    stop("No treatment specified")
  }else{
    nz <- dim(Zdata)[1]
    Zup <- Zdata[, treatment]
  }
  
  # check the type of intervention
  if (is.null(ras_Z)){
    cat("Point intervention", "\n")
    if (is.null(x_coord_Z) | is.null(x_coord_Z)){
      stop("Coordinates of intervention nodes are not detected")
    }
  }else if (!is.null(ras_Z)){
    cat("Polygon intervention", "\n")
    if (class(ras_Z)[1] == "RasterLayer"){
      # values in ras_Z indicate which polygon a cell belongs to
      ras_Z <- rasterToPolygons(ras_Z)
      ras_Z_coords <- gCentroid(ras_Z, byid = TRUE)@coords
      ras_Z <- st_as_sf(ras_Z)
    }else if (class(ras_Z)[1] == "SpatialPolygonsDataFrame"){
      ras_Z <- st_as_sf(ras_Z)
      ras_Z_coords <- st_coordinates(st_centroid(ras_Z))
    }else if (class(ras_Z)[1] == "sf"){
      ras_Z_coords <- st_coordinates(st_centroid(ras_Z))
    }else {
      stop("Unrecognized map format")
    }
    Zdata$x_coord_Z <- ras_Z_coords[, 1]
    Zdata$y_coord_Z <- ras_Z_coords[, 2]
    x_coord_Z <- "x_coord_Z"
    y_coord_Z <- "y_coord_Z"
  }
  
  if (dist.metric != "Geodesic" & dist.metric != "Euclidean"){
    stop("Unrecognized distance metric")
  }
  
  if (cType != "edge" & cType != "disk" & cType != "donut"){
    stop("Unrecognized circle type")
  }
  
  if (is.null(prob_treatment)){
    pweight <- rep(1, nz)
  }else{
    pweight <- Zup / Zdata[, "prob_treatment"] + (1 - Zup) / (1 - Zdata[, "prob_treatment"])
  }
  
  if (!is.null(covs)){
    n_covs <- dim(as.matrix(covs))[2]
    if (dim(as.matrix(covs))[1] == nz){
      covs <- matrix(rep(t(as.matrix(covs)), nz), ncol = ncol(as.matrix(covs)), byrow = TRUE)
    }
    if (!is.null(blockvar)){
      blocks <- model.matrix(~ blockvar - 1)
      blocks <- matrix(rep(t(blocks), nz), ncol = ncol(blocks), byrow = TRUE)
      covs <- cbind(covs, blocks)
      n_covs <- dim(as.matrix(covs))[2]
    }
    colnames(covs) <- paste0("X_", seq(1:n_covs))
  }
  
  if (is.null(dVec)){
    stop("No distance values")
  }else{
    dist_unit <- min(abs(dVec - c(NA, dVec[-length(dVec)])), na.rm = 1)
  }
  
  # check outcome and decide whether to impute with krig
  if (is.null(ras)){
    if (is.null(outcome)){
      stop("Outcome variable is not specified")
    }
    if (is.null(Ydata) & !is.null(Zdata[, outcome])){
      ras <- Krig(cbind(Zdata[, x_coord_Z], Zdata[, y_coord_Z]), Zdata[, outcome], m = m, lambda = lambda)
    }else if (!is.null(Ydata) & !is.null(Ydata[, outcome])){
      ras <- Krig(cbind(Ydata[, x_coord_Y], Ydata[, y_coord_Y]), Ydata[, outcome], m = m, lambda = lambda)
    }else{
      stop("Outcome variable is not specified")
    }
    Yobs <- NULL
    gridRes <- 1
    dtype <- "krig"
    ras_coords <- NULL
  }else{
    gridRes = res(ras)[1]
    Yobs <- getValues(ras)
    dtype <- "raster"
  }
  
  # Sdata is (mu_i(d), Z_i) for each d
  Sdata <- data.frame()
  Sdata.list <- list()
  cov.list <- list()
  AMR <- data.frame(d = dVec, taud = NA)
  AMR_est <- data.frame(d = dVec, taud_est = NA)
  Ybards <- Ybard_sums <- Ybard_lens <- matrix(NA, nz, length(dVec))
  cat("Calculate circle average for each distance value ")
  for (d in 1:length(dVec)){
    dUp <- dVec[d]
    Rim.list <- RimAvg(ras, Yobs, ras_Z, nz, Zup, Zdata, x_coord_Z, y_coord_Z, treatment, dUp, numpts, 
                       gridRes, evalpts, only.unique, dtype, dist.metric, cType)
    quiet(gc())    
    Ybard <- Rim.list[[1]]
    Ybard_sum <- Rim.list[[2]]
    Ybard_len <- Rim.list[[3]]
    if ((cType == "donut") & (d > 1)){
      Ybard_sums[, d] <- Ybard_sum
      Ybard_lens[, d] <- Ybard_len
      Ybards[, d] <- (Ybard_sums[, d] - Ybard_sums[, (d-1)]) / (Ybard_lens[, d] - Ybard_lens[, (d-1)])
    }else{
      Ybards[, d] <- Ybard
      Ybard_sums[, d] <- Ybard_sum
      Ybard_lens[, d] <- Ybard_len
    }
    AMR$taud[d] <- mean(Ybard[!is.nan(Ybard)])
    
    if (is.null(covs)){
      one_data <- cbind("outcome" = Ybards[, d], "treatment" = Zup, "dVec" = rep(dUp, nz), "w" = pweight)
      Sdata.list[[d]] <- one_data
      ols_formula <- as.formula("outcome ~ treatment")
      ols_fit <- lm(ols_formula, as.data.frame(one_data), weights = w)
      AMR_est$taud_est[d] <- coef(ols_fit)[2]
    }else{
      one_covs <- as.matrix(covs[((d-1)*nz + 1):(d*nz), ])
      for (c in 1:ncol(one_covs)){
        one_covs[, c] <- one_covs[, c] - mean(one_covs[, c])
        colnames(one_covs)[c] <- paste0("X_", c)
      }
      one_data <- cbind("outcome" = Ybards[, d], "treatment" = Zup, "dVec" = rep(dUp, nz), one_covs, "w" = pweight)
      Sdata.list[[d]] <- one_data
      ols_formula <- as.formula(paste0("outcome ~ ", paste0(paste0("treatment * ", colnames(one_covs)), collapse = "+")))
      ols_fit <- lm(ols_formula, as.data.frame(one_data), weights = w)
      AMR_est$taud_est[d] <- coef(ols_fit)[2]
      cov.list[[d]] <- as.matrix(model.matrix(ols_fit))
    }
    cat(d, " ")
  }
  Sdata <- data.frame(do.call("rbind", Sdata.list))
  rm(Sdata.list)
  names(Sdata)[1:4] <- c("outcome", "treatment", "dVec", "pweight")
  
  result.list <- list()

  result.list[["AMR_est"]] <- AMR_est
  c_n <- NULL
  if (per.se == 1){
    if(!require(ri)){
      install.packages("ri")
      library(ri)
    }
    permMat <- genperms(Zup, blockvar = blockvar, clustvar = clustvar, maxiter = nPerms)
    VCE.per <- matrix(nrow = nPerms, ncol = length(dVec))
    for(i in 1:ncol(permMat)){
      VCE.per[i, ] <- permMat[, i]%*%Ybards/sum(permMat[, i]) - (1-permMat[, i])%*%Ybards/sum(1-permMat[, i])
    }
    Per.CI <- apply(VCE.per, 2, function (x) quantile(x, c(alpha/2, 1-alpha/2), na.rm = TRUE))
    result.list[["Per.CI"]] <- t(Per.CI)
  }
  if (conley.se == 1){
    if (dist.metric == "Euclidean"){
      metric <- 1
    }else if (dist.metric == "Geodesic"){
      metric <- 2
    }
    x_coord <- Zdata[, x_coord_Z]
    y_coord <- Zdata[, y_coord_Z]
    dist <- DistanceCalculation(as.vector(x_coord), as.vector(y_coord), as.integer(metric))[["Dist_mat"]]
    if (kernel != "uni" & kernel != "uniform" & kernel != "tri" & kernel != "triangular" & kernel != "epa" & kernel != "epanechnikov" & kernel != "" ){
      stop("kernel incorrectly specified")
    }else if (kernel == "uni" | kernel == "uniform"){
      k <- 1
    }else if (kernel == "tri" | kernel == "triangular"){
      k <- 2
    }else if (kernel == "epa" | kernel == "epanechnikov"){
      k <- 3
    }
    Conley.SE <- rep(NA, length(dVec))
    Conley.CI <- matrix(NA, length(dVec), 2)
    if (dist.metric == "Euclidean"){
      metric <- 1
    }else if (dist.metric == "Geodesic"){
      metric <- 2
    }
    X_mat <- diag(sqrt(Sdata$pweight[Sdata$dVec == dVec[1]])) %*% cbind(rep(1, nz), Zup)
    W_meat <- X_mat
    XX_mat_inv <- solve(t(X_mat) %*% X_mat)
    z_vec <- c(0, 1) %*% XX_mat_inv %*% t(X_mat)
    M_mat <- diag(1, nz, nz) - X_mat %*% XX_mat_inv %*% t(X_mat)
    dofs <- rep(NA, length(dVec))
    c <- cutoff
    for (d in 1:length(dVec)){
      dUp <- dVec[d]
      mu_d <- diag(sqrt(Sdata$pweight[Sdata$dVec == dVec[1]])) %*% Ybards[, d]
      mu_d[is.nan(mu_d)] <- 0
      if (!is.null(covs)){
        X_mat <- cov.list[[d]]
        W_meat <- X_mat
        XX_mat_inv <- solve(t(X_mat) %*% X_mat)
        z_vec <- c(0, 1, rep(0, (ncol(X_mat)-2))) %*% XX_mat_inv %*% t(X_mat)
        M_mat <- diag(1, nz, nz) - X_mat %*% XX_mat_inv %*% t(X_mat)
        dofs <- rep(NA, length(dVec))
        c <- cutoff
      }
      beta_d <- solve(t(X_mat) %*% (X_mat)) %*% (t(X_mat) %*% mu_d)
      res_d <- mu_d - X_mat %*% beta_d
      if (cutoff > 0){
        c <- cutoff + d
      }
      Conley_result <- ConleySE(as.vector(res_d), as.matrix(W_meat), as.matrix(dist), 
                                as.double(c), as.integer(k))
      VCE_d <- t(XX_mat_inv) %*% Conley_result[["VCE_meat"]] %*% XX_mat_inv
      dist_kernel <- Conley_result[["Dist_kernel"]]
      if (edf){
        dof <- nz/((nz-2)*XX_mat_inv[2, 2])*sum(diag(M_mat %*% ((t(z_vec) %*% z_vec) * dist_kernel) %*% M_mat))
        temp= sum(diag (M_mat %*% (dist_kernel * t(z_vec) %*%  z_vec) %*% M_mat))
        Conley.SE[d] <- sqrt(VCE_d[2, 2]/dof)
        dofs[d] <- dof
      }else{
        Conley.SE[d] <- sqrt(VCE_d[2, 2])
      }
    }
    Conley.CI <- cbind(AMR_est[, 2] - qnorm(1-alpha/2) * Conley.SE, AMR_est[, 2] + qnorm(1-alpha/2) * Conley.SE)
    result.list[["Conley.SE"]] <- Conley.SE
    result.list[["Conley.CI"]] <- Conley.CI
  }
  rm(cov.list)
  if (smooth == 1){
    if (is.null(bw) | is.null(bw_debias)){
      bw.result <- CrossValidation(Sdata, outcome = "outcome", treatment = "treatment", dVec, 
                                   grid = NULL, nfold = 5, block_cv = TRUE, parallel = FALSE, 
                                   metric = "MSPE", kernel = kernel, bias_correction = bias_correction)
      bw <- bw.result[[2]]
      bw_debias <- bw.result[[3]]
    }
    wls_results <- LocalReg(dVec, Sdata, bw, bw_debias, Zup, xevals = NULL, smooth.conley.se, kernel, cutoff, dist, dist.metric, bias_correction)
    result.list[["AMR_est_smoothed"]] <- wls_results[["coefs"]]
    if (smooth.conley.se){
      if (edf){
        result.list[["smoothed.Conley.SE"]] <- wls_results[["ses"]]/sqrt(dofs)
      }else{
        result.list[["smoothed.Conley.SE"]] <- wls_results[["ses"]]
      }
      result.list[["smoothed.Conley.CI"]] <- cbind(wls_results[["coefs"]][, 2] - qnorm(1-alpha/2) * result.list[["smoothed.Conley.SE"]], 
                                                   wls_results[["coefs"]][, 2] + qnorm(1-alpha/2) * result.list[["smoothed.Conley.SE"]])
      if (conf.band){
        dVec_range <- max(dVec) - min(dVec)
        a_n <- 2*log(dVec_range/bw) + 2*log(0.5^(0.5)/(2*pi))
        c_n <- (a_n - 2*log(log((1-alpha)^(-0.5))))^(0.5)
        result.list[["smoothed.Conley.CB"]] <- cbind(wls_results[["coefs"]][, 2] - c_n * result.list[["smoothed.Conley.SE"]], 
                                                     wls_results[["coefs"]][, 2] + c_n * result.list[["smoothed.Conley.SE"]])
      }else{
        c_n <- NULL
      }
    }else{
      c_n <- NULL
    }
  }
  result.list[["Parameters"]] <- list("dVec" = dVec, "Zup" = Zup, "Ybards" = Ybards, "Sdata" = Sdata, 
                                      "blockvar" = blockvar, "clustvar" = clustvar, "nPerms" = nPerms, 
                                      "conley.se" = conley.se, "per.se" = per.se, "smooth" = smooth, 
                                      "smooth.conley.se" = smooth.conley.se, "bw" = bw, "bw_debias" = bw_debias, "c_n" = c_n)
  result.list$call <- match.call()
  class(result.list) <- "SpatialEffect"
  return(result.list)
}

SpatialEffectTest <- function(result.list, test.range, smooth = 0, alpha = 0.05){
  
  dVec <- result.list[["Parameters"]][[1]]
  Zup <- result.list[["Parameters"]][[2]]
  Ybards <- result.list[["Parameters"]][[3]]
  Sdata <- result.list[["Parameters"]][[4]]
  blockvar <- result.list[["Parameters"]][[5]]
  clustvar <- result.list[["Parameters"]][[6]]
  nPerms <- result.list[["Parameters"]][[7]]
  AMR_est <- result.list[["AMR_est"]]
  
  if (!require(ri)){
    install.packages("ri")
    library(ri)
  }
  if (is.null(test.range)){
    warning("test.range is not specified.")
    test.range <- c(min(dVec), max(dVec))
  }else if (length(test.range) != 2){
    warning("test.range can take only two values.")
    test.range <- c(min(dVec), max(dVec))
  }else if (test.range[1] > test.range[2]){
    warning("the smaller value should come first.")
    test.range <- test.range[2, 1]
  }
  if (test.range[1] < min(dVec)){
    test.range[1] <- min(dVec)
  }else if (test.range[2] > max(dVec)){
    test.range[2] <- max(dVec)
  }
  test.stat <- sum(AMR_est[AMR_est[, 1] >= test.range[1] & AMR_est[, 1] <= test.range[2], 2])
  
  if (smooth == 1){
    AMR_est_smoothed <- result.list[["AMR_est_smoothed"]]
    test.stat <- sum(AMR_est_smoothed[AMR_est_smoothed[, 1] >= test.range[1] & AMR_est_smoothed[, 1] <= test.range[2], 2])
    Sdata$zIndex <- rep(seq(1, nz, 1), length(dVec))
  }
  
  permMat <- genperms(Zup, blockvar = blockvar, clustvar = clustvar, maxiter = nPerms)
  test.per <- rep(NA, nPerms)
  for(i in 1:ncol(permMat)){
    if (smooth == 0){
      AMR.per <- permMat[, i]%*%Ybards/sum(permMat[, i]) - (1-permMat[, i])%*%Ybards/sum(1-permMat[, i])
      test.per[i] <- sum(AMR.per[AMR_est[, 1] >= test.range[1] & AMR_est[, 1] <= test.range[2]])
    }else if (smooth == 1){
      treatment.per <- data.frame(cbind(1:(length(Zup)), permMat[, i]))
      names(treatment.per) <- c("zIndex", "treatment.per")
      Sdata.per <- merge(Sdata, treatment.per, by = "zIndex")
      Sdata.per$treatment <- Sdata.per$treatment.per
      AMR.per <- LocalReg(dVec, Sdata.per, bw, bw_debias, Zup, xevals = NULL, smooth.conley.se = 0, kernel, cutoff, dist, dist.metric, bias_correction)[["coefs"]][, 2]
      test.per[i] <- sum(AMR.per[AMR_est_smoothed[, 1] >= test.range[1] & AMR_est_smoothed[, 1] <= test.range[2]])
    }
  }
  test.CI <- quantile(test.per, c(alpha/2, 1-alpha/2), na.rm = TRUE)
  test.list <- list(test.stat, test.CI)
  return(test.list)
}

summary.SpatialEffect <- function(object, dVec.range = NULL) {
  x    <- object
  #if (is.null(args[['level']])) { level <- 0.05 } else { level <- args[['level']] }
  
  cat("Call: SpatialEffect\n\n")
  output.data <- x[["AMR_est"]]
  cnames <- c("dVec", "AMR_est")
  if (x[["Parameters"]][[8]] == 1){
    output.data <- cbind(output.data, x[["Conley.CI"]])
    cnames <- c(cnames, "Conley.CI.l", "Conley.CI.u")
  }
  if (x[["Parameters"]][[9]] == 1){
    output.data <- cbind(output.data, x[["Per.CI"]])
    cnames <- c(cnames, "Per.CI.l", "Per.CI.u")
  }
  if (x[["Parameters"]][[10]] == 1){
    output.data <- cbind(output.data, x[["AMR_est_smoothed"]][, 2])
    cnames <- c(cnames, "AMR_est_smoothed")
  }
  if (x[["Parameters"]][[11]] == 1){
    output.data <- cbind(output.data, x[["smoothed.Conley.CI"]])
    cnames <- c(cnames, "smoothed.Conley.CI.l", "smoothed.Conley.CI.u")
  }
  names(output.data) <- cnames
  output.data <- apply(output.data, 2, function(x) round(x, 3))
  if (is.null(dVec.range)){
    print(output.data)
  }else {
    if (length(dVec.range) != 2){
      stop("Invalide range of distance values")
    }
    if (dVec.range[1] > dVec.range[2]){
      d <- dVec.range[2]
      dVec.range[2] <- dVec.range[1]
      dVec.range[1] <- d
    }
    print(output.data[output.data[, 1] >= dVec.range[1] & output.data[, 1] <= dVec.range[2], ])
  } 
}

print.SpatialEffect <- function(object, dVec.range = NULL) {
  x    <- object
  #if (is.null(args[['level']])) { level <- 0.05 } else { level <- args[['level']] }
  
  cat("Call: SpatialEffect\n\n")
  output.data <- x[["AMR_est"]]
  cnames <- c("dVec", "AMR_est")
  if (x[["Parameters"]][[8]] == 1){
    output.data <- cbind(output.data, x[["Conley.CI"]])
    cnames <- c(cnames, "Conley.CI.l", "Conley.CI.u")
  }
  if (x[["Parameters"]][[9]] == 1){
    output.data <- cbind(output.data, x[["Per.CI"]])
    cnames <- c(cnames, "Per.CI.l", "Per.CI.u")
  }
  if (x[["Parameters"]][[10]] == 1){
    output.data <- cbind(output.data, x[["AMR_est_smoothed"]][, 2])
    cnames <- c(cnames, "AMR_est_smoothed")
  }
  if (x[["Parameters"]][[11]] == 1){
    output.data <- cbind(output.data, x[["smoothed.Conley.CI"]])
    cnames <- c(cnames, "smoothed.Conley.CI.l", "smoothed.Conley.CI.u")
  }
  names(output.data) <- cnames
  output.data <- apply(output.data, 2, function(x) round(x, 3))
  if (is.null(dVec.range)){
    print(output.data)
  }else {
    if (length(dVec.range) != 2){
      stop("Invalide range of distance values")
    }
    if (dVec.range[1] > dVec.range[2]){
      d <- dVec.range[2]
      dVec.range[2] <- dVec.range[1]
      dVec.range[1] <- d
    }
    print(output.data[output.data[, 1] >= dVec.range[1] & output.data[, 1] <= dVec.range[2], ])
  } 
}

