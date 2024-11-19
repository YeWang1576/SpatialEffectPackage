# simulation results in Wang, Samii, Chang, and Aronow (2024)
# generate results for the appendix
rm(list=ls())

setwd("~/Dropbox/CurrentProjects/SpatialReplication/simulation/graphs")

library(Rcpp)
library(RcppArmadillo)
library(plotrix)
library(foreign)
library(sp)
library(fields)
library(ri)
library(sf)
library(raster)

sourceCpp("~/Dropbox/CurrentProjects/SpatialReplication/SpatialEffect/src/DistanceCalculation.cpp")
sourceCpp("~/Dropbox/CurrentProjects/SpatialReplication/SpatialEffect/src/Conley.cpp")
source("~/Dropbox/CurrentProjects/SpatialReplication/simulation/scripts/functions.R")

set.seed(2024)

# simulation results (polygon intervention)
type <- "nonmono"
data_path <- "~/Dropbox/CurrentProjects/SpatialReplication/simulation/data/"

load(paste0(data_path, "simData_", type, "_polygon_varyingsizes.RData"))

sample_sizes <- c(20, 40, 60, 80, 100, 120)

len <- 10
dVec <- seq(from=.5, to=len, by=.25)

# initial bandwidth
bw <- 2
bw_debias <- 4

AME_list <- list()
AME_est_list <- AME_est_list_smooth <- list()
all_VCEs_list <- all_VCEs_list_smooth <- list()
Allbias <- MSEs <- SupNorms <- matrix(NA, length(sample_sizes), length(dVec))
Allbias_smooth <- MSEs_smooth <- SupNorms_smooth <- matrix(NA, length(sample_sizes), length(dVec))

for (ss in 1:length(sample_sizes)){
  sample_size <- sample_sizes[ss]

  result_list <- data_list[[ss]]
  Ydata <- result_list[[1]]
  Zdata <- result_list[[2]]
  ras0 <- result_list[[3]]
  ras_Z <- result_list[[4]]
  names(Zdata)[2:3] <- c("x", "y")
  ny <- dim(Ydata)[1]
  nz <- dim(Zdata)[1]
  YobsAll <- Ydata[, grep("Y1_", names(Ydata))]
  ZMat <- t(Zdata[, grep("curZ", names(Zdata))])
  
  AMEmat <- data.frame(d = dVec, tauda2 = NA, nd = NA)
  
  for(dIndex in 1:length(dVec)){
    dUp <- dVec[dIndex]
    AMEmat$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz, numpts = NULL, only.unique = 0)
  }
  
  AME_list[[ss]] <- AMEmat[, 1:2]
  
  AMEhat_mean <- matrix(NA, length(dVec), nrow(ZMat))
  AMEhat_mean_smooth <- matrix(NA, length(dVec), nrow(ZMat))
  all_VCEs <- all_VCEs_smooth <- matrix(NA, length(dVec), nrow(ZMat))
  bias <- MSE <- SupNorm <- matrix(NA, length(dVec), nrow(ZMat))
  bias_smooth <- MSE_smooth <- SupNorm_smooth <- matrix(NA, length(dVec), nrow(ZMat))
  
  # parameters for estimation
  x_coord_Z <- "x"
  y_coord_Z <- "y"
  dist.metric <- "Euclidean"
  numpts <- NULL
  evalpts <- 1
  smooth <- 1
  smooth.conley.se <- 1
  only.unique <- 0
  alpha <- 0.05
  per.se <- 1
  blockvar <- NULL
  clustvar <- NULL
  nPerms <- 1000
  conley.se <- 1
  kernel <- "tri"
  cutoff <- 6
  bias_correction <- TRUE
  k <- 1
  
  grid_list <- circleMean(ras0 = ras0, ras_Z = ras_Z, dVec = dVec, nz = ncol(ZMat))
  
  x_coord <- Zdata[, x_coord_Z]
  y_coord <- Zdata[, y_coord_Z]
  dist <- DistanceCalculation(as.vector(x_coord), as.vector(y_coord), 1)[["Dist_mat"]]
  
  pweight <- rep(1, nz)
  for(zIndex in 1:nrow(ZMat)){
    nz <- dim(Zdata)[1]
    treatment <- paste0("curZ", zIndex)
    Zup <- Zdata[, treatment]
    Yobs <- Ydata[, paste0("Y1_", zIndex)]
    Sdata <- data.frame()
    Sdata_list <- list()
    AME <- data.frame(d = dVec, taud = NA)
    AME_est <- data.frame(d = dVec, taud_est = NA)
    Ybards <- matrix(NA, nz, length(dVec))
    for (d in 1:length(dVec)){
      dUp <- dVec[d]
      Ybard <- rep(NA, nz)
      for (i in 1:nz){
        grid_index <- grid_list[[(i-1)*length(dVec)+d]]
        grids <- Yobs[grid_index]
        Ybard[i] <- mean(grids, na.rm = TRUE)
      }
      Ybards[, d] <- Ybard
      AME$taud[d] <- mean(Ybard)
      Sdata_list[[d]] <- cbind(Ybard, Zup, rep(dUp, nz), pweight)
      AME_est$taud_est[d] <- Zup%*%Ybard/sum(Zup) - (1-Zup)%*%Ybard/sum(1-Zup)
    }
    Sdata <- data.frame(do.call("rbind", Sdata_list))
    names(Sdata) <- c("outcome", "treatment", "dVec", "pweight")
    
    AMEhat_mean[, zIndex] <- AME_est[, 2]
    
    X_mat <- cbind(rep(1, nz), Zup)
    W_meat <- X_mat
    XX_mat_inv <- solve(t(X_mat) %*% X_mat)
    z_vec <- c(0, 1) %*% XX_mat_inv %*% t(X_mat)
    M_mat <- diag(1, nz, nz) - X_mat %*% XX_mat_inv %*% t(X_mat)
    for (d in 1:length(dVec)){
      if (ss > 2 & d < 12) {
        c <- cutoff + 2*dVec[d]
      } else{
        c <- cutoff + 2*dVec[11]
      }
      dUp <- dVec[d]
      mu_d <- Ybards[, d]
      beta_d <- solve(t(X_mat) %*% (X_mat)) %*% (t(X_mat) %*% mu_d)
      res_d <- mu_d - X_mat %*% beta_d
      Conley_result <- ConleySE(as.vector(res_d), as.matrix(W_meat), as.matrix(dist), 
                                as.double(c), as.integer(1))
      VCE_d <- t(XX_mat_inv) %*% Conley_result[["VCE_meat"]] %*% XX_mat_inv
      dist_kernel <- Conley_result[["Dist_kernel"]]
      Conley_SE <- sqrt(VCE_d[2, 2])
      all_VCEs[d, zIndex] <- Conley_SE
    }
    bias[, zIndex] <- round(AMEmat$tauda2 - AME_est[, 2], 5)
    MSE[, zIndex] <- round((AMEmat$tauda2 - AME_est[, 2])^2, 5)
    SupNorm[, zIndex] <- round(abs(AMEmat$tauda2 - AME_est[, 2]), 5)
    
    if (zIndex %% 100 == 0){
      cat(zIndex, "\n")
    }
  }
  
  AME_est_list[[ss]] <- AMEhat_mean
  all_VCEs_list[[ss]] <- all_VCEs
  bias <- bias[!is.nan(bias[, 1]), ]
  Allbias[ss, ] <- apply(bias, 1, mean)
  MSEs[ss, ] <- apply(MSE, 1, mean)
  SupNorms[ss, ] <- apply(SupNorm, 1, max)
  cat("Sample size is: ", sample_size, "\n")
}

save(AME_list, file = paste0(data_path, "AME_", type, "_polygon_bernoulli_smooth_varyingsizes.RData"))
save(AME_est_list, file = paste0(data_path, "AMEest_", type, "_polygon_bernoulli_varyingsizes.RData"))
save(all_VCEs_list, file = paste0(data_path, "simVCEs_", type, "_polygon_bernoulli_varyingsizes.RData"))
save(Allbias, MSEs, SupNorms, file = paste0(data_path, "simResult_", type, "_polygon_bernoulli_varyingsizes.RData"))


########################################## HT estimator #######################################
load(paste0(data_path, "simData_", type, "_polygon_varyingsizes.RData"))
result_list <- data_list[[4]]
Ydata <- result_list[[1]]
Zdata <- result_list[[2]]
ras0 <- result_list[[3]]
names(Zdata)[2:3] <- c("x", "y")
ny <- dim(Ydata)[1]
nz <- dim(Zdata)[1]
YobsAll <- Ydata[, grep("Y1_", names(Ydata))]
ZMat <- t(Zdata[, grep("curZ", names(Zdata))])

len <- 10
dVec <- seq(from=.5, to=len, by=.25)

AMEmat <- data.frame(d = dVec, tauda2 = NA, nd = NA)

for(dIndex in 1:length(dVec)){
  dUp <- dVec[dIndex]
  AMEmat$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz, numpts = NULL, only.unique = 0)
}

AMEhat_mean <- matrix(NA, length(dVec), nrow(ZMat))

# parameters for estimation
x_coord_Z <- "x"
y_coord_Z <- "y"
dist.metric <- "Euclidean"
ras_Z <- NULL
cutoff <- 6
k <- 1

grid_list <- circleMean(ras0 = ras0, ras_Z = ras_Z, dVec = dVec, nz = ncol(ZMat))

x_coord <- Zdata[, x_coord_Z]
y_coord <- Zdata[, y_coord_Z]
dist <- DistanceCalculation(as.vector(x_coord), as.vector(y_coord), 1)[["Dist_mat"]]

pweight <- rep(1, nz)
for(zIndex in 1:nrow(ZMat)){
  nz <- dim(Zdata)[1]
  treatment <- paste0("curZ", zIndex)
  Zup <- Zdata[, treatment]
  Yobs <- Ydata[, paste0("Y1_", zIndex)]
  Sdata <- data.frame()
  Sdata_list <- list()
  AME_est <- data.frame(d = dVec, taud_est = NA)
  Ybards <- matrix(NA, nz, length(dVec))
  for (d in 1:length(dVec)){
    dUp <- dVec[d]
    Ybard <- rep(NA, nz)
    for (i in 1:nz){
      grid_index <- grid_list[[(i-1)*length(dVec)+d]]
      grids <- Yobs[grid_index]
      Ybard[i] <- mean(grids, na.rm = TRUE)
    }
    Ybards[, d] <- Ybard
    Sdata_list[[d]] <- cbind(Ybard, Zup, rep(dUp, nz), pweight)
    AME_est$taud_est[d] <- Zup%*%Ybard/(nz * 0.5) - (1-Zup)%*%Ybard/(nz * 0.5)
  }
  Sdata <- data.frame(do.call("rbind", Sdata_list))
  names(Sdata) <- c("outcome", "treatment", "dVec", "pweight")
  
  AMEhat_mean[, zIndex] <- AME_est[, 2]
  if (zIndex %% 100 == 0){
    cat(zIndex, "\n")
  }
}

AME <- AMEmat[, 1:2]
AME_est <- AMEhat_mean
save(AME, AME_est, file = paste0(data_path, "simResult_", type, "_bernoulli_HT.RData"))


############################## complete randomization ####################################
load(paste0(data_path, "simData_", type, "_cr.RData"))

Ydata <- result_list[[1]]
Zdata <- result_list[[2]]
ras0 <- result_list[[3]]
ras_Z <- NULL
names(Zdata)[2:3] <- c("x", "y")
ny <- dim(Ydata)[1]
nz <- dim(Zdata)[1]
YobsAll <- Ydata[, grep("Y1_", names(Ydata))]
ZMat <- t(Zdata[, grep("curZ", names(Zdata))])

len <- 10
dVec <- seq(from=.5, to=len, by=.25)

AMEmat <- data.frame(d = dVec, tauda2 = NA, nd = NA)

for(dIndex in 1:length(dVec)){
  dUp <- dVec[dIndex]
  AMEmat$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz, numpts = NULL, only.unique = 0)
}

AMEhat_mean <- matrix(NA, length(dVec), nrow(ZMat))

# parameters for estimation
x_coord_Z <- "x"
y_coord_Z <- "y"
dist.metric <- "Euclidean"
ras_Z <- NULL
cutoff <- 6
k <- 1

grid_list <- circleMean(ras0 = ras0, ras_Z = ras_Z, dVec = dVec, nz = ncol(ZMat))

x_coord <- Zdata[, x_coord_Z]
y_coord <- Zdata[, y_coord_Z]
dist <- DistanceCalculation(as.vector(x_coord), as.vector(y_coord), 1)[["Dist_mat"]]

pweight <- rep(1, nz)
for(zIndex in 1:nrow(ZMat)){
  nz <- dim(Zdata)[1]
  treatment <- paste0("curZ", zIndex)
  Zup <- Zdata[, treatment]
  Yobs <- Ydata[, paste0("Y1_", zIndex)]
  Sdata <- data.frame()
  Sdata_list <- list()
  AME_est <- data.frame(d = dVec, taud_est = NA)
  Ybards <- matrix(NA, nz, length(dVec))
  for (d in 1:length(dVec)){
    dUp <- dVec[d]
    Ybard <- rep(NA, nz)
    for (i in 1:nz){
      grid_index <- grid_list[[(i-1)*length(dVec)+d]]
      grids <- Yobs[grid_index]
      Ybard[i] <- mean(grids, na.rm = TRUE)
    }
    Ybards[, d] <- Ybard
    Sdata_list[[d]] <- cbind(Ybard, Zup, rep(dUp, nz), pweight)
    AME_est$taud_est[d] <- Zup%*%Ybard/sum(Zup) - (1-Zup)%*%Ybard/sum(1-Zup)
  }
  Sdata <- data.frame(do.call("rbind", Sdata_list))
  names(Sdata) <- c("outcome", "treatment", "dVec", "pweight")
  
  AMEhat_mean[, zIndex] <- AME_est[, 2]
  if (zIndex %% 100 == 0){
    cat(zIndex, "\n")
  }
}

AME <- AMEmat[, 1:2]
AME_est <- AMEhat_mean
save(AME, AME_est, file = paste0(data_path, "simResult_", type, "_cr.RData"))



########################### null effect #########################################
type <- "null"
load(paste0(data_path, "simData_", type, "_null.RData"))

Ydata <- result_list[[1]]
Zdata <- result_list[[2]]
ras0 <- result_list[[3]]
ras_Z <- NULL
names(Zdata)[2:3] <- c("x", "y")
ny <- dim(Ydata)[1]
nz <- dim(Zdata)[1]
YobsAll <- Ydata[, grep("Y1_", names(Ydata))]
ZMat <- t(Zdata[, grep("curZ", names(Zdata))])

AMEmat <- data.frame(d = dVec, tauda2 = NA, nd = NA)

for(dIndex in 1:length(dVec)){
  dUp <- dVec[dIndex]
  AMEmat$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz, numpts = NULL, only.unique = 0)
}

AMEhat_mean <- matrix(NA, length(dVec), nrow(ZMat))

# parameters for estimation
x_coord_Z <- "x"
y_coord_Z <- "y"
dist.metric <- "Euclidean"
ras_Z <- NULL
cutoff <- 6
k <- 1

grid_list <- circleMean(ras0 = ras0, ras_Z = ras_Z, dVec = dVec, nz = ncol(ZMat))

x_coord <- Zdata[, x_coord_Z]
y_coord <- Zdata[, y_coord_Z]
dist <- DistanceCalculation(as.vector(x_coord), as.vector(y_coord), 1)[["Dist_mat"]]

all_VCEs <- matrix(NA, length(dVec), nrow(ZMat))
pweight <- rep(1, nz)
for(zIndex in 1:nrow(ZMat)){
  nz <- dim(Zdata)[1]
  treatment <- paste0("curZ", zIndex)
  Zup <- Zdata[, treatment]
  Yobs <- Ydata[, paste0("Y1_", zIndex)]
  Sdata <- data.frame()
  Sdata_list <- list()
  AME_est <- data.frame(d = dVec, taud_est = NA)
  Ybards <- matrix(NA, nz, length(dVec))
  for (d in 1:length(dVec)){
    dUp <- dVec[d]
    Ybard <- rep(NA, nz)
    for (i in 1:nz){
      grid_index <- grid_list[[(i-1)*length(dVec)+d]]
      grids <- Yobs[grid_index]
      Ybard[i] <- mean(grids, na.rm = TRUE)
    }
    Ybards[, d] <- Ybard
    Sdata_list[[d]] <- cbind(Ybard, Zup, rep(dUp, nz), pweight)
    AME_est$taud_est[d] <- Zup%*%Ybard/sum(Zup) - (1-Zup)%*%Ybard/sum(1-Zup)
  }
  Sdata <- data.frame(do.call("rbind", Sdata_list))
  names(Sdata) <- c("outcome", "treatment", "dVec", "pweight")
  
  AMEhat_mean[, zIndex] <- AME_est[, 2]
  
  X_mat <- cbind(rep(1, nz), Zup)
  W_meat <- X_mat
  XX_mat_inv <- solve(t(X_mat) %*% X_mat)
  z_vec <- c(0, 1) %*% XX_mat_inv %*% t(X_mat)
  M_mat <- diag(1, nz, nz) - X_mat %*% XX_mat_inv %*% t(X_mat)
  for (d in 1:length(dVec)){
    c <- cutoff + 2*dVec[d]
    dUp <- dVec[d]
    mu_d <- Ybards[, d]
    beta_d <- solve(t(X_mat) %*% (X_mat)) %*% (t(X_mat) %*% mu_d)
    res_d <- mu_d - X_mat %*% beta_d
    Conley_result <- ConleySE(as.vector(res_d), as.matrix(W_meat), as.matrix(dist), 
                              as.double(c), as.integer(k))
    VCE_d <- t(XX_mat_inv) %*% Conley_result[["VCE_meat"]] %*% XX_mat_inv
    dist_kernel <- Conley_result[["Dist_kernel"]]
    Conley_SE <- sqrt(VCE_d[2, 2])
    all_VCEs[d, zIndex] <- Conley_SE
  }
  if (zIndex %% 100 == 0){
    cat(zIndex, "\n")
  }
}

AME <- AMEmat[, 1:2]
AME_est <- AMEhat_mean
save(AME, AME_est, all_VCEs, file = paste0(data_path, "simResult_", type, "_bernoulli_HA.RData"))


################################# FRT ##########################################
len <- 10
dVec <- seq(from=.5, to=len, by=.25)
sample_sizes <- c(20, 40, 60, 80, 100, 120)

type <- "nonmono"
data_path <- "~/Dropbox/CurrentProjects/SpatialReplication/simulation/data/"
load(paste0(data_path, "simData_", type, "_varyingsizes.RData"))

AMR_list <- list()
AMR_est_list <- list()
FRT_coverage <- matrix(NA, length(sample_sizes), length(dVec))

x_coord_Z <- "x"
y_coord_Z <- "y"
dist.metric <- "Euclidean"
only.unique <- 0
alpha <- 0.05
per.se <- 1
blockvar <- NULL
clustvar <- NULL
nPerms <- 1000

for (ss in 1:length(sample_sizes)){
  sample_size <- sample_sizes[ss]
  result_list <- data_list[[ss]]
  Ydata <- result_list[[1]]
  Zdata <- result_list[[2]]
  ras0 <- result_list[[3]]
  ras_Z <- NULL
  ny <- dim(Ydata)[1]
  nz <- dim(Zdata)[1]
  YobsAll <- Ydata[,grep("Y1_", names(Ydata))]
  ZMat <- t(Zdata[,grep("curZ", names(Zdata))])
  
  AMRmat <- data.frame(d = dVec, tauda2 = NA, nd = NA)
  
  for(dIndex in 1:length(dVec)){
    dUp <- dVec[dIndex]
    AMRmat$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz, numpts = NULL, only.unique = 0)
  }
  
  AMR_list[[ss]] <- AMRmat[, 1:2]
  
  AMRhat.mean <- matrix(NA, length(dVec), nrow(ZMat))
  
  grid_list <- circleMean(ras0 = ras0, ras_Z = ras_Z, dVec = dVec, nz = ncol(ZMat))
  
  x_coord <- Zdata[, x_coord_Z]
  y_coord <- Zdata[, y_coord_Z]
  
  FRT_covered <- matrix(NA, nrow(ZMat), length(dVec))
  for(zIndex in 1:nrow(ZMat)){
    nz <- dim(Zdata)[1]
    treatment <- paste0("curZ", zIndex)
    Zup <- Zdata[, treatment]
    Yobs <- Ydata[, paste0("Y1_", zIndex)]
    Sdata <- data.frame()
    Sdata.list <- list()
    AMR <- data.frame(d = dVec, taud = NA)
    AMR_est <- data.frame(d = dVec, taud_est = NA)
    Ybards <- matrix(NA, nz, length(dVec))
    for (d in 1:length(dVec)){
      dUp <- dVec[d]
      Ybard <- rep(NA, nz)
      for (i in 1:nz){
        grid_index <- grid_list[[(i-1)*length(dVec)+d]]
        grids <- Yobs[grid_index]
        Ybard[i] <- mean(grids, na.rm = TRUE)
      }
      Ybards[, d] <- Ybard
      AMR$taud[d] <- mean(Ybard)
      Sdata.list[[d]] <- cbind(Ybard, Zup, rep(dUp, nz))
      AMR_est$taud_est[d] <- Zup%*%Ybard/sum(Zup) - (1-Zup)%*%Ybard/sum(1-Zup)
    }
    Sdata <- data.frame(do.call("rbind", Sdata.list))
    names(Sdata) <- c("outcome", "treatment", "dVec")
    
    AMRhat.mean[, zIndex] <- AMR_est[, 2]
    
    permMat <- genperms(Zup, blockvar = blockvar, clustvar = clustvar, maxiter = nPerms)
    VCE.per <- matrix(nrow = nPerms, ncol = length(dVec))
    for(i in 1:ncol(permMat)){
      VCE.per[i, ] <- permMat[, i]%*%Ybards/sum(permMat[, i]) - (1-permMat[, i])%*%Ybards/sum(1-permMat[, i])
    }
    Per.CI <- apply(VCE.per, 2, function (x) quantile(x, c(alpha/2, 1-alpha/2), na.rm = TRUE))
    FRT_covered[zIndex, ] <- Per.CI[1, ] <= AMR_est[, 2] & Per.CI[2, ] >= AMR_est[, 2]
    
    if (zIndex %% 100 == 0){
      cat(zIndex, "\n")
    }
  }
  
  AMR_est_list[[ss]] <- AMRhat.mean
  FRT_coverage[ss, ] <- 1 - apply(FRT_covered, 2, mean)
  
  cat("Sample size is: ", sample_size, "\n")
}

save(FRT_coverage, file = "~/Dropbox/CurrentProjects/SpatialReplication/simulation/data/FRT_nonmono_point_bernoulli_varyingsizes.RData")
