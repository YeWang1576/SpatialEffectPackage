# simulation results in Wang, Samii, Chang, and Aronow (2024)
rm(list=ls())

setwd("~/Dropbox/CurrentProjects/SpatialReplication/simulation/graphs/appendix")

library(plotrix)
library(foreign)
library(sp)
library(fields)
library(raster)
library(sf)

set.seed(2024)
source("~/Dropbox/CurrentProjects/SpatialReplication/simulation/scripts/functions.R")

type <- "nonmono"
label <- "additive"
data_path <- "~/Dropbox/CurrentProjects/SpatialReplication/simulation/data/"

# plotting
load(paste0(data_path, "simData_", type, "_polygon_varyingsizes.RData"))
load(paste0(data_path, "AME_", type, "_polygon_bernoulli_smooth_varyingsizes.RData"))
load(paste0(data_path, "AMEest_", type, "_polygon_bernoulli_varyingsizes.RData"))
load(paste0(data_path, "simVCEs_", type, "_polygon_bernoulli_varyingsizes.RData"))
load(paste0(data_path, "simResult_", type, "_polygon_bernoulli_varyingsizes.RData"))

sample_sizes <- c(20, 40, 60, 80, 100, 120)
len <- 10
dVec <- seq(from=.5, to=len, by=.25)
bias_colors <- colorRampPalette(c("blue", "red"))(length(sample_sizes))

# Figure 6  
AME <- AME_list[[4]]
all_VCEs <- all_VCEs_list[[4]]
AME_est <- AME_est_list[[4]]
result_list <- data_list[[4]]
Ydata <- result_list[[1]]
Zdata <- result_list[[2]]
ras0 <- result_list[[3]]
ras_Z <- result_list[[4]]
ny <- dim(Ydata)[1]
nz <- dim(Zdata)[1]
YobsAll <- Ydata[,grep("Y1_", names(Ydata))]
ZMat <- t(Zdata[,grep("curZ", names(Zdata))])
true_var <- apply(AME_est, 1, sd)
est_var_avg <- apply(all_VCEs, 1, function(x) mean(x, na.rm = TRUE))

pdf(file=paste0("AME-unbiased-", type, "-polygon-bernoulli-HA.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(AME[,1:2],
     ylab=paste0("Estimates (", label, ")"),
     xlab="Distance (coordinate degrees)",
     type="n",
     ylim=c(-2, 5),
     main="Distribution of the Estimates (Hajek)",
     cex.lab=2,
     cex.main=2,
     cex.axis=2)
for(zIndex in 1:nrow(ZMat)){
  points(AME_est[, zIndex]~dVec, col="gray", type="l", lwd = 3)
}
points(AME[, 1:2], type="l", col = "red", lwd = 3)
legend("topright", col = c("red", "gray"), lty=c(1,1), lwd = c(3, 3),
       legend = c("True AME", "AME estimates"), cex=2, bty = "n")
dev.off()


AME_est_avg <- apply(AME_est, 1, mean)
CI95_real <- apply(AME_est, 1, function(x) quantile(x, c(0.025, 0.975)))
SE_Conley <- apply(all_VCEs, 1, function(x) mean(x[!is.nan(x)], na.rm = 1))
CI_Conley <- rbind(AME_est_avg - 1.96 * SE_Conley, AME_est_avg + 1.96 * SE_Conley)
pdf(file=paste0("AME-CI-comparison-", type, "-polygon-bernoulli-HA.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(	AME[,1:2],
      ylab=paste0("Estimates (", label, ")"),
      xlab="Distance (coordinate degrees)",
      type="l",
      ylim=c(-2, 5),
      main="Comparison of CIs (Spatial HAC vs. Simulated)",
      col="red",
      cex.lab=2,
      cex.main=2,
      cex.axis=2,
      lwd = 3)
points(AME_est_avg ~ dVec, col="black", type="l", lwd = 3)
points(CI_Conley[1, ] ~ dVec, col="black", type="l", lty = 3, lwd = 2)
points(CI_Conley[2, ] ~ dVec, col="black", type="l", lty = 3, lwd = 2)
points(CI95_real[1, ] ~ dVec, col="red", type="l", lty = 3, lwd = 2)
points(CI95_real[2, ] ~ dVec, col="red", type="l", lty = 3, lwd = 2)
legend("topright", col = c("red", "black", "red", "black"), lty=c(1,1,3,3), lwd = c(3,3,2,2),
       legend = c("True AME", "AME estimates", "Simulated 95% CIs", "Spatial HAC 95% CIs"), cex=2, bty = "n")
dev.off()

coverage <- matrix(NA, length(sample_sizes), length(dVec))
true_vars <- est_vars_avg <- matrix(NA, length(sample_sizes), length(dVec))
for (ss in 1:length(sample_sizes)){
  all_VCEs <- all_VCEs_list[[ss]]
  all_VCEs[is.nan(all_VCEs)] <- NA
  all_VCEs[all_VCEs == Inf] <- NA
  AME_est <- AME_est_list[[ss]]
  result_list <- data_list[[ss]]
  Ydata <- result_list[[1]]
  Zdata <- result_list[[2]]
  ras0 <- result_list[[3]]
  ny <- dim(Ydata)[1]
  nz <- dim(Zdata)[1]
  YobsAll <- Ydata[,grep("Y1_", names(Ydata))]
  ZMat <- t(Zdata[,grep("curZ", names(Zdata))])
  true_vars[ss, ] <- apply(AME_est, 1, sd)
  est_vars_avg[ss, ] <- apply(all_VCEs, 1, function(x) mean(x, na.rm = TRUE))

  all_CIs_upper <- all_CIs_lower <- matrix(NA, length(dVec), 1000)
  covered <- matrix(0, length(dVec), 1000)
  AME <- AME_list[[ss]]
  for (i in 1:1000){
    all_CIs_upper[, i] <- AME_est[, i] + 1.96*all_VCEs[, i]
    all_CIs_lower[, i] <- AME_est[, i] - 1.96*all_VCEs[, i]
    covered[, i] <- (AME$tauda2 >= all_CIs_lower[, i] & AME$tauda2 <= all_CIs_upper[, i])
  }
  coverage[ss, ] <- apply(covered, 1, function(x) mean(x, na.rm = TRUE))
}

# Figure 15
pdf(file=paste0("AME-CI-coverage-", type, "polygon-bernoulli-HA.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(0,
     ylab=paste0("Coverage rate (", label, ")"),
     xlab="Distance (coordinate degrees)",
     type="n",
     xlim=c(0, 10),
     ylim=c(0, 1),
     main="Coverage of the 95% Spatial HAC CIs",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2,
     cex.axis=2)
for (i in 1:nrow(coverage)){
  points(coverage[i, ]~dVec, type = "l", col = bias_colors[i], lwd = 2)
}
abline(h = 0.95, lty = 2)
legend("bottomleft", col = bias_colors, lty=c(1,1), lwd = 2,
       legend = paste0("N=", sample_sizes^2/100), cex=2, bty = "n")
dev.off()

# bias, MSE, SupNorm
pdf(file=paste0("AME-bias-", type, "-polygon-bernoulli-HA.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(4,3,4,2))
par(pty="s")
plot(0,
     ylab=paste0("Bias (", label, ")"),
     xlab="Distance (coordinate degrees)",
     type="n",
     xlim=c(0, 10),
     ylim=c(-.3, 1.5),
     main="Bias of the HA Estimator",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2)
for (i in 1:nrow(Allbias)){
  points(Allbias[i, ]~dVec, type = "l", col = bias_colors[i], lwd = 2)
}
legend("topright", col = bias_colors, lty=c(1,1), lwd = 2,
       legend = paste0("N=", sample_sizes^2/100), cex=2, bty = "n")
dev.off()

pdf(file=paste0("AME-MSE-", type, "-polygon-bernoulli-HA.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(0,
     ylab=paste0("MSE (", label, ")"),
     xlab="Distance (coordinate degrees)",
     type="n",
     xlim=c(0, 10),
     ylim=c(-.3, 1.5),
     main="MSE of the Hajek Estimator",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2,
     cex.axis=2)
for (i in 1:nrow(Allbias)){
  points(MSEs[i, ]~dVec, type = "l", col = bias_colors[i], lwd = 2)
}
legend("topright", col = bias_colors, lty=c(1,1), lwd = 2,
       legend = paste0("N=", sample_sizes^2/100), cex=2, bty = "n")
dev.off()

pdf(file=paste0("AME-SupNorm-", type, "-polygon-bernoulli-HA.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(4,3,4,2))
par(pty="s")
plot(0,
     ylab=paste0("Supreme norm of the bias (", label, ")"),
     xlab="Distance (coordinate degrees)",
     type="n",
     xlim=c(0, 10),
     ylim=c(-.3, 3),
     main="Supreme Norm of the HA Estimator's Bias",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2)
for (i in 1:nrow(Allbias)){
  points(SupNorms[i, ]~dVec, type = "l", col = bias_colors[i], lwd = 2)
}
legend("topright", col = bias_colors, lty=c(1,1), lwd = 2,
       legend = paste0("N=", sample_sizes^2/100), cex=2, bty = "n")
dev.off()


# Figure 13 top-left 
len <- 10
dVec <- seq(from=.5, to=len, by=.25)

load(paste0(data_path, "simResult_", type, "_bernoulli_HT.RData"))

pdf(file=paste0("AME-unbiased-", type, "-bernoulli-HT.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(AME,
     ylab=paste0("Estimates (additive)"),
     xlab="Distance (coordinate degrees)",
     type="n",
     ylim=c(-2, 5),
     main="Distribution of the Estimates (Horvitz-Thompson)",
     cex.lab=2,
     cex.main=2,
     cex.axis=2)
for(zIndex in 1:nrow(ZMat)){
  points(AME_est[, zIndex]~dVec, col="gray", type="l", lwd = 3)
}
points(AME, type="l", col = "red", lwd = 3)
legend("topright", col = c("red", "gray"), lty=c(1,1), lwd = c(3, 3),
       legend = c("True AME", "AME estimates"), cex=2, bty = "n")
dev.off()


# Figure 13 top-right 
len <- 10
dVec <- seq(from=.5, to=len, by=.25)

load(paste0(data_path, "simResult_", type, "_cr.RData"))

pdf(file=paste0("AME-unbiased-", type, "-cr.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(AME,
     ylab=paste0("Estimates (additive)"),
     xlab="Distance (coordinate degrees)",
     type="n",
     ylim=c(-2, 5),
     main="Distribution of the Estimates (Hajek)",
     cex.lab=2,
     cex.main=2,
     cex.axis=2)
for(zIndex in 1:nrow(ZMat)){
  points(AME_est[, zIndex]~dVec, col="gray", type="l", lwd = 3)
}
points(AME, type="l", col = "red", lwd = 3)
legend("topright", col = c("red", "gray"), lty=c(1,1), lwd = c(3, 3),
       legend = c("True AME", "AME estimates"), cex=2, bty = "n")
dev.off()


# Figure 13 bottom
len <- 10
dVec <- seq(from=.5, to=len, by=.25)

load(paste0(data_path, "simResult_", type, "_bernoulli_HA.RData"))

pdf(file=paste0("AME-unbiased-", type, "-bernoulli-HA.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(AME,
     ylab=paste0("Estimates (null)"),
     xlab="Distance (coordinate degrees)",
     type="n",
     ylim=c(-2, 5),
     main="Distribution of the Estimates (Hajek)",
     cex.lab=2,
     cex.main=2,
     cex.axis=2)
for(zIndex in 1:nrow(ZMat)){
  points(AME_est[, zIndex]~dVec, col="gray", type="l", lwd = 3)
}
points(AME, type="l", col = "red", lwd = 3)
legend("topright", col = c("red", "gray"), lty=c(1,1), lwd = c(3, 3),
       legend = c("True AME", "AME estimates"), cex=2, bty = "n")
dev.off()

AME_est_avg <- apply(AME_est, 1, mean)
CI95_real <- apply(AME_est, 1, function(x) quantile(x, c(0.025, 0.975)))
SE_Conley <- apply(all_VCEs, 1, function(x) mean(x[!is.nan(x)], na.rm = 1))
CI_Conley <- rbind(AME_est_avg - 1.96 * SE_Conley, AME_est_avg + 1.96 * SE_Conley)
pdf(file=paste0("AME-CI-comparison-", type, "-bernoulli-HA.pdf"), height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(	AME[,1:2],
      ylab=paste0("Estimates (", type, ")"),
      xlab="Distance (coordinate degrees)",
      type="l",
      ylim=c(-2, 5),
      main="Comparison of CIs (Spatial HAC vs. Simulated)",
      col="red",
      cex.lab=2,
      cex.main=2,
      cex.axis=2,
      lwd = 3)
points(AME_est_avg ~ dVec, col="black", type="l", lwd = 3)
points(CI_Conley[1, ] ~ dVec, col="black", type="l", lty = 3, lwd = 2)
points(CI_Conley[2, ] ~ dVec, col="black", type="l", lty = 3, lwd = 2)
points(CI95_real[1, ] ~ dVec, col="red", type="l", lty = 3, lwd = 2)
points(CI95_real[2, ] ~ dVec, col="red", type="l", lty = 3, lwd = 2)
legend("topright", col = c("red", "black", "red", "black"), lty=c(1,1,3,3), lwd = c(3,3,2,2),
       legend = c("True AME", "AME estimates", "Simulated 95% CIs", "Spatial HAC 95% CIs"), cex=2, bty = "n")
dev.off()


# Figure 16
load("~/Dropbox/CurrentProjects/SpatialReplication/simulation/data/FRT_nonmono_point_bernoulli_varyingsizes.RData")
bias_colors <- colorRampPalette(c("blue", "red"))(length(sample_sizes))
pdf(file="~/Dropbox/CurrentProjects/SpatialReplication/simulation/graphs/appendix/AMR-perCI-coverage-nonmono-bernoulli-HA.pdf", height=8, width=8)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(0,
     ylab="Probability of rejection",
     xlab="Distance (coordinate degrees)",
     type="n",
     xlim=c(0, 10),
     ylim=c(0, 1),
     main="Probability of rejecting the sharp null",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2)
for (i in 1:nrow(FRT_coverage)){
  points(FRT_coverage[i, ]~dVec, type = "b", pch = 17, col = bias_colors[i], cex=.8, lwd = 1.5)
}
abline(h = 0.95, lty = 2)
legend("topright", col = bias_colors, pch = 17,
       legend = paste0("N=", sample_sizes^2/100), cex=1.5, bty = "n")
dev.off()
