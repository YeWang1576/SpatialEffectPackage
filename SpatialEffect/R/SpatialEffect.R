SpatialEffect <- function(ras = NULL, Ydata = NULL, outcome = NULL, x_coord_Y = NULL, y_coord_Y = NULL, 
                          Zdata, x_coord_Z, y_coord_Z, treatment, dVec, dist.metric = "Euclidean", 
                          numpts = NULL, evalpts = 1, only.unique = 0, per.se = 1, blockvar = NULL, 
                          clustvar = NULL, conley.se = 1, cutoff = 0, alpha = 0.05, m = 2, lambda = 0.02, 
                          nPerms = 1000){
  
  if (is.null(Zdata) | is.null(treatment)){
    stop("No intervention point")
  }else{
    nz <- dim(Zdata)[1]
    Zup <- Zdata[, treatment]
  }
  
  if (dist.metric != "Geodesic" & dist.metric != "Euclidean"){
    stop("Unrecognized distance metric")
  }else if(dist.metric == "Geodesic"){
    if(!require(geosphere)){
      install.packages("geosphere")
      library(geosphere)
    }
  }
  
  if (is.null(ras)){
    if (is.null(outcome)){
      stop("Outcome variable is not specified")
    }
    if(!require(fields)){
      install.packages("fields")
      library(fields)
    }
    if (is.null(Ydata)){
      ras <- Krig(cbind(Zdata[, x_coord_Z], Zdata[, y_coord_Z]), Zdata[, outcome], m = m, lambda = lambda)
    }else{
      ras <- Krig(cbind(Ydata[, x_coord_Y], Ydata[, y_coord_Y]), Ydata[, outcome], m = m, lambda = lambda)
    }
    Yobs <- NULL
    gridRes <- 1
    dtype <- "krig"
    ras_coords <- NULL
  }else{
    if(!require(raster)){
      install.packages("raster")
      library(raster)
    }
    gridRes = res(ras)[1]
    Yobs <- getValues(ras)
    dtype <- "raster"
    ras_coords <- coordinates(ras)
  }
  
  Sdata <- data.frame()
  dist_mat <- data.frame()
  Sdata.list <- list()
  dist_mat.list <- list()
  MIRmat <- data.frame(d = dVec, tauda2 = NA)
  MIRhat.mat <- data.frame(d = dVec, taudhat = NA)
  Ybards <- matrix(NA, nz, length(dVec))
  for (d in 1:length(dVec)){
    dUp <- dVec[d]
    Rim.list <- RimAvg(ras, ras_coords, Yobs, Zdata, x_coord_Z, y_coord_Z, dUp, numpts, gridRes, evalpts, only.unique, dtype, dist.metric)
    Ybard <- Rim.list[[1]]
    Ybards[, d] <- Ybard
    MIRmat$tauda2[d] <- mean(Ybard)
    Sdata.list[[d]] <- Rim.list[[2]]
    dist_mat.list[[d]] <- Rim.list[[3]]
    MIRhat.mat$taudhat[d] <- Zup%*%Ybard/sum(Zup) - (1-Zup)%*%Ybard/sum(1-Zup)
  }
  Sdata <- data.frame(do.call("rbind", Sdata.list))
  rm(Sdata.list)
  dist_mat <- data.frame(do.call("rbind", dist_mat.list))
  rm(dist_mat.list)
  names(Sdata) <- c("X1", "X2", "X3")
  names(dist_mat) <- paste0("d", c(1:nz))
  
  result.list <- list()
  result.list[["MIRhat.mat"]] <- MIRhat.mat
  
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
    #Per.SE <- apply(VCE.per, 2, function (x) sd(x))
    #result.list[["Per.SE"]] <- Per.SE
    result.list[["Per.CI"]] <- t(Per.CI)
  }
  if (conley.se == 1){
    reg_list <- RegRep(Sdata, dist_mat, dVec, Zup)
    Sdata <- reg_list[[1]]
    X_mat <- reg_list[[2]]
    coefs <- reg_list[[3]]
    res <- reg_list[[4]]
    if (dist.metric == "Euclidean"){
      metric <- 1
    }else if (dist.metric == "Geodesic"){
      metric <- 2
    }
    VCE <- ConleySE(as.vector(res), as.matrix(X_mat), as.vector(Sdata$X2), 
                    as.vector(Sdata$X3), as.double(cutoff), as.double(metric))
    Conley.SE <- sqrt(diag(VCE))
    Conley.CI <- cbind(MIRhat.mat[, 2] - qnorm(1-alpha/2) * Conley.SE, MIRhat.mat[, 2] + qnorm(1-alpha/2) * Conley.SE)
    result.list[["Conley.SE"]] <- Conley.SE
    result.list[["Conley.CI"]] <- Conley.CI
  }
  return(result.list)
}
