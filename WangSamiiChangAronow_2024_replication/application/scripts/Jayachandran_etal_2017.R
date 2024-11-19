rm(list=ls())

#setwd("~/Dropbox/CurrentProjects/SpatialReplication/application/graphs")

library(Rcpp)
library(RcppArmadillo)
library(plotrix)
library(foreign)
library(sp)
library(fields)
library(ri)
library(sf)
library(raster)
library(terra)
library(ddpcr)
library(geosphere)
library(pbapply)
set.seed(2024)


###############################################################################
###################Set Path####################################################
###############################################################################

root="C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/SpatialReplication/"
setwd(root)
data_path <- "application/data/"

###############################################################################
###################Import Function#############################################
###############################################################################
source('simulation/scripts/6_functions.R')


###############################################################################
###################Import Data#################################################
###############################################################################

#get outcome data
ras <- raster(paste0(data_path, "/rct_uganda_gfc_updated.tif"), na.rm = TRUE, band = 3)
#ras_2012 <- raster(paste0(data_path, "/rct_uganda_gfc_updated.tif"), na.rm = TRUE, band = 1)
#ras_2013 <- raster(paste0(data_path, "/rct_uganda_gfc_updated.tif"), na.rm = TRUE, band = 2)
ras = aggregate(ras,20,fun=mean) #coarsen the raster for simplicity
outcome <- getValues(ras)
Yxy=coordinates(ras)
Yxy=st_as_sf(SpatialPoints(Yxy,proj4string = crs(ras)))
Yindex <- 1:nrow(Yxy)
Ydata <- data.frame(Yindex,cbind(outcome))

#get treatment data
map <- st_read(dsn = paste0(data_path, "/digitized_rct"), layer = "rct_studyarea_vectorized")
Zxy = st_as_sf(map)
centr <- st_coordinates(st_centroid(map))
Zdata <- cbind(st_set_geometry(map, NULL)[, c(1:3)], centr)
Zdata <- Zdata[, c(4, 5, 2)]
names(Zdata) <- c("x_coord", "y_coord", "treatment")

#create raster
ras_Z <- st_as_sf(as(map, "Spatial"))
#names(ras_Z@data)[2] <- "treatment"

##############################################################################
###############Calculate Distances ###########################################
##############################################################################

nz <- dim(Zdata)[1]

#distance matrix from evaluation points to intervention nodes

# start.time=proc.time()
# dist_matrix_YZ=matrix(0,nrow(Yxy),nrow(Zxy))
# 
# #this one took about two hours to run
# for (i in 1:nrow(Yxy)){
#   print(i)
#   dist_matrix_YZ[i,] = st_distance(Yxy[i,],Zxy)
# }
# print(proc.time()-start.time)
#save(list=c('dist_matrix_YZ'),file='dist_matrix_YZ.Rdata')

load('application/dist_matrix_YZ.Rdata')

#dist_matrix_YZ = pbsapply(1:nrow(Yxy),function(i) return(st_distance(Yxy[i,],Zxy)))
Ydata=cbind(Ydata,dist_matrix_YZ)
names(Ydata)[3:ncol(Ydata)] <- sapply(1:nz,function(i) return(paste0("dZ",i)))

#distance matrix from intervention nodes to intervention nodes
dist_matrix_ZZ <-pbsapply(1:nrow(Zxy),function(i) return(st_distance(Zxy,Zxy[i,])))
Zdata=cbind(Zdata,dist_matrix_ZZ)
names(Zdata)[4:ncol(Zdata)] <- sapply(1:nz,function(i) return(paste0("dZ",i)))


##############################################################################
###############Calculate Circle Averages #####################################
##############################################################################

len <- 30
dVec <- seq(from=0, to=len, by=1)*1000 #in meters

grid_list <- circleMean(Zdata=Zdata,ras0 = ras, ras_Z = ras_Z, dVec = dVec, nz = nrow(Zdata)) 



ras <- ras
Zdata <- Zdata
x_coord_Z <- "x_coord"
y_coord_Z <- "y_coord"
treatment <- "treatment"
numpts <- 1000
evalpts <- 1
only.unique <- 0
per.se <- 1
blockvar <- NULL
clustvar <- NULL
conley.se <- 1
alpha <- 0.1
cutoff <- 6
nPerms <- 500
dist.metric <- "Geodesic"
smooth <- 1
bw <- 1
kernel <- "uni"
edf <- 0
trim <- 0
smooth.conley.se <- 1
cType <- "donut"
covs <- NULL
prob_treatment <- NULL

len <- 30
dVec <- seq(from=0, to=len, by=1)

# maps
pdf(file="J2017-forest-map-polygon-circles-treated.pdf", height=8, width=8)
par(mfrow=c(1,1))
plot(ras, bty="n", box=FALSE, col = c("white", "darkgoldenrod4"), legend=FALSE, 
     cex.lab=2,
     cex.main=2,
     cex.axis = 2)
lines(ras_Z[Zdata$treatment == 0, ], col="darkturquoise", lwd = 2)
lines(ras_Z[Zdata$treatment == 1, ], col="goldenrod1", lwd = 2)
for (i in which(Zdata$treatment == 1)){
  buffer_treated <- buffer(ras_Z[i, ], width=0.015 * 111000)
  lines(buffer_treated, col="goldenrod1", lwd = 2, lty = 2)
}
legend("bottomright", col = c("darkturquoise", "goldenrod1", "darkgoldenrod4"),
       legend = c("Control", "Treated", "Deforestation"), 
       lty=c(1, 1, NA), pch=c(NA, NA, 19), cex=1.5, pt.cex = c(.8, .8, .4), bty = "n")
dev.off()

pdf(file="J2017-forest-map-polygon-circles-control.pdf", height=8, width=8)
par(mfrow=c(1,1))
plot(ras, bty="n", box=FALSE, col = c("white", "darkgoldenrod4"), legend=FALSE,
     cex.lab=2,
     cex.main=2,
     cex.axis = 2)
lines(map[Zdata$treatment == 0, ], col="darkturquoise", lwd = 2)
lines(map[Zdata$treatment == 1, ], col="goldenrod1", lwd = 2)
for (i in which(Zdata$treatment == 0)){
  buffer_control <- buffer(ras_Z[i, ], width=0.015 * 111000)
  lines(buffer_control, col="darkturquoise", lwd = 2, lty = 2)
}
legend("bottomright", col = c("darkturquoise", "goldenrod1", "darkgoldenrod4"),
       legend = c("Control", "Treated", "Deforestation"), 
       lty=c(1, 1, NA), pch=c(NA, NA, 19), cex=1.5, pt.cex = c(.8, .8, .4), bty = "n")
dev.off()

sr <- "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
ras <- projectRaster(ras, crs = sr)
ras_Z <- spTransform(ras_Z, crs(ras))

centr <- st_coordinates(st_centroid(map))
Zdata <- cbind(st_set_geometry(map, NULL)[, c(1:3)], centr)
Zdata <- Zdata[, c(4, 5, 2)]
names(Zdata) <- c("x_coord", "y_coord", "treatment")

# analysis
result.list <- SpatialEffect(ras = ras, Zdata = Zdata, x_coord_Z = x_coord_Z, y_coord_Z = y_coord_Z, ras_Z = ras_Z, 
                             treatment = treatment, dVec = dVec, numpts = numpts, only.unique = only.unique, smooth = smooth, 
                             per.se = per.se, conley.se = conley.se, cutoff = cutoff, alpha = alpha, edf = edf, trim = trim, bw = 10.996, bw_debias = 25.13,
                             nPerms = nPerms, smooth.conley.se = smooth.conley.se, dist.metric = dist.metric, cType = "donut", kernel = kernel)

save(result.list, file = "../data/forestResult_reproj.RData")

load("../data/forestResult_reproj.RData")
AMR_est <- result.list[["AMR_est"]]
Per.CI <- result.list[["Per.CI"]]
Conley.SE <- result.list[["Conley.SE"]]
Conley.CI <- result.list[["Conley.CI"]]
AMR_est_smoothed <- result.list[["AMR_est_smoothed"]]
smoothed.Conley.SE <- result.list[["smoothed.Conley.SE"]]
smoothed.Conley.CI <- result.list[["smoothed.Conley.CI"]]
smoothed.Conley.CB <- result.list[["smoothed.Conley.CB"]]

test.result <- SpatialEffectTest(result.list, test.range = c(1, 5), smooth = 0, alpha = 0.1)

dVec_real <- dVec[-1]


pdf(file="J2017_AME_polygon_unsmoothed_reproj_90CI.pdf", height=8, width=8)
par(mar = c(5, 5, 5, 5))
plot(	AMR_est[1:17, 2] ~ dVec_real[1:17],
      ylab="Estimates",
      xlab="Distance (kilometers)",
      type="l",
      ylim=c(-0.04, 0.04),
      main = "Hajek Estimator with Its CIs",
      cex.lab=2,
      cex.main=2,
      cex.axis = 2,
      lwd = 3)
points(Conley.CI[1:17, 1]~dVec_real[1:17], type = "l", col = "black", lty = "dotted", lwd = 3)
points(Conley.CI[1:17, 2]~dVec_real[1:17], type = "l", col = "black", lty = "dotted", lwd = 3)
points(Per.CI[1:17, 1]~dVec_real[1:17], type = "l", col = "blue", lty = "dotted", lwd = 3)
points(Per.CI[1:17, 2]~dVec_real[1:17], type = "l", col = "blue", lty = "dotted", lwd = 3)
abline(h = 0, lty = 2, lwd = 3, col = "gray")
legend("topright", col = c("black", "black", "blue"), lty=c(1,3,3), lwd = c(3, 3, 3),
       legend = c("AME estimates", "Spatial HAC 95% CI", "Quantiles under sharp null"), cex=2, bty = "n")
dev.off()

pdf(file="J2017_AME_polygon_smoothed_reproj_90CI.pdf", height=8, width=8)
par(mar = c(5, 5, 5, 5))
plot(	AMR_est_smoothed[1:17, 2] ~ dVec_real[1:17],
      ylab="Estimates",
      xlab="Distance (kilometers)",
      type="l",
      ylim=c(-0.04, 0.04),
      main = "Kernel Regression Estimator with Its CIs",
      cex.lab=2,
      cex.main=2,
      cex.axis = 2,
      lwd = 3)
points(smoothed.Conley.CI[1:17, 1]~dVec_real[1:17], type = "l", col = "black", lty = "dotted", lwd = 3)
points(smoothed.Conley.CI[1:17, 2]~dVec_real[1:17], type = "l", col = "black", lty = "dotted", lwd = 3)
abline(h = 0, lty = 2, lwd = 3, col = "gray")
legend("topright", col = c("black", "black"), lty=c(1,3), lwd = c(3, 3),
       legend = c("AME estimates", "Spatial HAC 95% CI"), cex=2, bty = "n")
dev.off()


ha_ests <- AMR_est[1:17, 2]
ha_CI95_upper <- Conley.CI[1:17, 2]
ha_CI95_lower <- Conley.CI[1:17, 1]
ha_CI90_upper <- AMR_est[1:17, 2] + qnorm(0.95) * Conley.SE[1:17]
ha_CI90_lower <- AMR_est[1:17, 2] - qnorm(0.95) * Conley.SE[1:17]

names <- dVec[1:17]

ha.data <- data.frame("estimates" = ha_ests * 100, "CI95_u" = ha_CI95_upper * 100, 
                     "CI95_l" = ha_CI95_lower * 100, "CI90_u" = ha_CI90_upper * 100, 
                     "CI90_l" = ha_CI90_lower * 100, "names" = names)
ha.data$names <- factor(ha.data$names, levels = rev(rev(ha.data$names)))

colors <- rep("red", dim(ha.data)[1])

require(ggplot2)
coefs_p <- ggplot(ha.data) + 
  geom_pointrange(aes(x=names, y=estimates, ymin=CI95_l, ymax=CI95_u), size = 0.3, colour=colors, lty = "dashed") + 
  geom_pointrange(aes(x=names, y=estimates, ymin=CI90_l, ymax=CI90_u), size = 0.6, colour=colors, lty = "solid") + 
  theme_bw() + geom_hline(aes(yintercept = 0), colour = "black", lty=2) + xlab('Distance (kilometers)') +
  scale_linetype_manual(labels = c("90% CI", "95% CI"), values=c(1, 2)) + ylim(-3, 1) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "bottom", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + ylab("Estimates")
ggsave(paste0("J2017_graphs/ha_coefs_reproj_90CI.pdf"), coefs_p, device = "pdf", height = 8, width = 8) 


smooth_ests <- AMR_est_smoothed[1:17, 2]
smooth_CI95_upper <- smoothed.Conley.CI[1:17, 2]
smooth_CI95_lower <- smoothed.Conley.CI[1:17, 1]
smooth_CI90_upper <- AMR_est_smoothed[1:17, 2] + qnorm(0.95) * smoothed.Conley.SE[1:17]
smooth_CI90_lower <- AMR_est_smoothed[1:17, 2] - qnorm(0.95) * smoothed.Conley.SE[1:17]
colors <- rep("red", length(smooth_ests))

names <- dVec[1:17]

smooth.data <- data.frame("estimates" = smooth_ests * 100, "CI95_u" = smooth_CI95_upper * 100, 
                      "CI95_l" = smooth_CI95_lower * 100, "CI90_u" = smooth_CI90_upper * 100, 
                      "CI90_l" = smooth_CI90_lower * 100, "names" = names)
smooth.data$names <- factor(smooth.data$names, levels = rev(rev(smooth.data$names)))

# varying cutoff values
Sdata <- result.list$Parameters$Sdata
nz <- nrow(Zdata)
Zup <- Zdata$treatment

Zmat <- Zdata[, c(1:2)]

dist <- DistanceCalculation(as.vector(Zmat[, 1]), as.vector(Zmat[, 2]), as.integer(2))[["Dist_mat"]]
dist <- dist / 1000

cutoffs <- seq(2, 10, 2)

smoothed.Conley.CIs <- Conley.CIs <- list()
cc <- 1
bw <- 18.707
smoothed.est <- rep(NA, length(dVec))
for (cutoff in cutoffs){
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

  smoothed.Conley.SE <- rep(NA, length(dVec))

  if (dist.metric == "Euclidean"){
    metric <- 1
  }else if (dist.metric == "Geodesic"){
    metric <- 2
  }
  X_mat <- cbind(rep(1, nz), Zup)
  W_meat <- X_mat
  XX_mat_inv <- solve(t(X_mat) %*% X_mat)
  z_vec <- c(0, 1) %*% XX_mat_inv %*% t(X_mat)
  M_mat <- diag(1, nz, nz) - X_mat %*% XX_mat_inv %*% t(X_mat)
  dofs <- rep(NA, length(dVec))
  cvs <- rep(NA, length(dVec))
  for (d in 1:length(dVec)){
    dUp <- dVec[d]
    mu_d <- Sdata[Sdata$dVec == dVec[d], "outcome"]
    beta_d <- solve(t(X_mat) %*% (X_mat)) %*% (t(X_mat) %*% mu_d)
    res_d <- mu_d - X_mat %*% beta_d
    c <- 2*cutoff + 2*dUp
    Conley_result <- ConleySE(as.vector(res_d), as.matrix(W_meat), as.matrix(dist),
                              as.double(c), as.integer(k), as.integer(1))
    VCE_d <- t(XX_mat_inv) %*% Conley_result[["VCE_meat"]] %*% XX_mat_inv
    # if (VCE_d[2, 2] < 0){
    #   Conley_result <- ConleySE(as.vector(res_d), as.matrix(W_meat), as.matrix(dist),
    #                             as.double(c), as.integer(k), as.integer(1))
    #   VCE_d <- t(XX_mat_inv) %*% Conley_result[["VCE_meat"]] %*% XX_mat_inv
    # }
    dist_kernel <- Conley_result[["Dist_kernel"]]
    # Conley.SE[d] <- sqrt(VCE_d[2, 2])
    dof_mat <- M_mat %*% ((t(z_vec) %*% z_vec) * dist_kernel) %*% M_mat
    dofs[d] <- nz/((nz-2)*XX_mat_inv[2, 2])*sum(diag(dof_mat))
    nu <- (nz/((nz-2)*XX_mat_inv[2, 2]))^2*sum(diag(dof_mat %*% dof_mat))
    cvs[d] <- (2*dofs[d]^2) / nu
    Conley.SE[d] <- sqrt(VCE_d[2, 2]/dofs[d])
    
    # smoothed outcome
    # mu_all <- matrix(Sdata[, "outcome"], nz, length(dVec)) 
    # d_weight <- (1-abs((dVec - dUp)/bw)) * (abs((dVec - dUp)/bw) <= 1)
    # mu_d_smoothed <- mu_all %*% d_weight
    # beta_d_smoothed <- solve(t(X_mat) %*% (X_mat)) %*% (t(X_mat) %*% mu_d_smoothed)
    # smoothed.est[d] <- beta_d_smoothed[2]
    # res_d_smoothed <- mu_d_smoothed - X_mat %*% beta_d_smoothed
    # Conley_result_smoothed <- ConleySE(as.vector(res_d_smoothed), as.matrix(W_meat), as.matrix(dist),
    #                                    as.double(c), as.integer(k), as.integer(1))
    # VCE_d_smoothed <- t(XX_mat_inv) %*% Conley_result_smoothed[["VCE_meat"]] %*% XX_mat_inv
    # smoothed.Conley.SE[d] <- sqrt(VCE_d_smoothed[2, 2]/dofs[d])
  }
  # cv <- qnorm(1-alpha/2)
  cv <- qt(1-alpha/2, cvs)
  Conley.CI <- cbind(AMR_est[, 2] - cv * Conley.SE, AMR_est[, 2] + cv * Conley.SE)
  Conley.CIs[[cc]] <- Conley.CI
  
  # smoothed.Conley.CI <- cbind(smoothed.est - cv * smoothed.Conley.SE, smoothed.est + cv * smoothed.Conley.SE)
  # smoothed.Conley.CIs[[cc]] <- smoothed.Conley.CI
  wls_results <- LocalReg(dVec, Sdata, bw = 18.707, bw_debias = 25.13, dist = dist, Zup, 
                          xevals = NULL, smooth.conley.se = 1, kernel, as.numeric(cutoff), 
                          as.integer(1), edf = TRUE, dist.metric, bias_correction = 0)
  smoothed.est <- wls_results[["coefs"]][, 2]
  smoothed.Conley.CIs[[cc]] <- cbind(wls_results[["coefs"]][, 2] - qt(1-alpha/2, wls_results[["edof"]]) * wls_results[["ses"]],
                                     wls_results[["coefs"]][, 2] + qt(1-alpha/2, wls_results[["edof"]]) * wls_results[["ses"]])

  
  cc <- cc + 1
  cat(cc, "\n")
}

save(smoothed.Conley.CIs, Conley.CIs, smoothed.est, file = "../data/forestResult_reproj_varyingbandwidth.RData")

dVec_real <- dVec
CI_colors <- colorRampPalette(c("blue", "red"))(length(cutoffs))
pdf(file="J2017_AME_polygon_HA_varyingbandwidth_proj_90CI.pdf", height=8, width=8)
par(mar = c(5, 5, 5, 5))
plot(	AMR_est[1:17, 2] ~ dVec_real[1:17],
      ylab="Estimates",
      xlab="Distance (kilometers)",
      type="l",
      ylim=c(-0.04, 0.05),
      main = "Hajek Estimator with Its CIs",
      cex.lab=2,
      cex.main=2,
      cex.axis = 2,
      lwd = 3)
for (cc in 1:5){
  points(Conley.CIs[[cc]][1:17, 1]~dVec_real[1:17], type = "l", col = CI_colors[cc], lty = "dotted", lwd = 3)
  points(Conley.CIs[[cc]][1:17, 2]~dVec_real[1:17], type = "l", col = CI_colors[cc], lty = "dotted", lwd = 3)
}
legend("topright", col = c("black", CI_colors[1:5]), lty=c(1, rep(3, 5)), lwd = rep(3, 6),
       legend = c("AME estimates", paste0("90% CI with cutoff at ", cutoffs[1:5], " km")), cex=2, bty = "n")
dev.off()


pdf(file="J2017_AME_polygon_smoothed_varyingbandwidth_proj_90CI.pdf", height=8, width=8)
par(mar = c(5, 5, 5, 5))
plot(	smoothed.est[1:17] ~ dVec[1:17],
      ylab="Estimates",
      xlab="Distance (kilometers)",
      type="l",
      ylim=c(-0.04, 0.05),
      main = "Kernel Regression Estimator with Its CIs",
      cex.lab=2,
      cex.main=2,
      cex.axis = 2,
      lwd = 3)
for (cc in 1:5){
  points(smoothed.Conley.CIs[[cc]][1:17, 1]~dVec[1:17], type = "l", col = CI_colors[cc], lty = "dotted", lwd = 3)
  points(smoothed.Conley.CIs[[cc]][1:17, 2]~dVec[1:17], type = "l", col = CI_colors[cc], lty = "dotted", lwd = 3)
}
legend("topright", col = c("black", CI_colors[1:5]), lty=c(1, rep(3, 5)), lwd = rep(3, 6),
       legend = c("AME estimates", paste0("90% CI with cutoff at ", cutoffs[1:5], " km")), cex=2, bty = "n")
dev.off()
