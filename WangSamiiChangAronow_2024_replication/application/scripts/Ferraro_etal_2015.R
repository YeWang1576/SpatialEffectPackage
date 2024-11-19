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

set.seed(2024)

data_path <- "~/Dropbox/CurrentProjects/SpatialReplication/application/data/"

protected_areas <- read_sf("/Users/yewang/Documents/GitHub/spatial/data/application/Costa Rica/PAs_Reproj/PAs_before80_1.shp")
forest <- st_read("/Users/yewang/Documents/GitHub/spatial/data/application/Costa Rica/Deforestation and Carbon/Spatial/final_pixels_w_covs.shp")
forest <- st_transform(forest, st_crs(protected_areas)) # transform forest to the crs of protected areas

invalid_protected_areas <- st_is_valid(protected_areas, NA_on_exception = TRUE)

sum(!invalid_protected_areas, na.rm = TRUE)  # Report the number of invalid geometries for protected areas

protected_areas <- st_make_valid(protected_areas) # If there are invalid geometries, attempt to repair them
plot(protected_areas[, 2])

crs_protected_areas <- st_crs(protected_areas) # Extract the crs of protected areas
proj_string <- crs_protected_areas$proj4string

buffer_distance <- 5000 # Create buffers around the protected areas with the radius of 5km
protected_areas_buffers <- st_buffer(protected_areas, dist = buffer_distance)

full_extent_sf <- st_union(st_union(protected_areas, protected_areas_buffers)) # Convert to Spatial for raster compatibility
full_extent_sp <- as(full_extent_sf, "Spatial")
res <- 5000 # Create a blank raster with desired resolution; resolution in units of the CRS
blank_raster <- raster(extent(full_extent_sp), res=res)
crs(blank_raster) <- CRS(proj_string)
st_crs(protected_areas) <- crs(blank_raster)
st_crs(protected_areas_buffers) <- crs(blank_raster)

protected_areas_raster <- rasterize(as(protected_areas, 'Spatial'), blank_raster, field=1, fun='mean', background=0)
buffers_raster <- rasterize(as(protected_areas_buffers, 'Spatial'), blank_raster, field=1, fun='mean', background=0)
protected_areas_raster[protected_areas_raster == 0] <- NA
buffers_raster[buffers_raster == 0] <- NA
buffers_raster[buffers_raster == 1] <- 0

combined_raster <- cover(protected_areas_raster, buffers_raster, cover=0)
table(getValues(combined_raster))
plot(combined_raster)

boundary_function <- function(x) {
  if (length(unique(x)) > 1) {
    return(1)  # Returns 1 for boundary cells
  } else {
    return(NA)  # Returns NA for non-boundary cells
  }
}

# Apply the focal operation with a 3x3 window to check all 8 neighbors around each cell
window_mat <- matrix(1, 3, 3)

boundary_tiles <- focal(combined_raster, w = window_mat, fun = boundary_function, pad = TRUE, padValue = 0)
boundary_tiles[is.na(combined_raster)] <- NA
boundary_tiles[combined_raster == 0] <- 0

plot(boundary_tiles)
table(getValues(boundary_tiles))


mean(forest$defor97[forest$prot_b80 == 1]) - mean(forest$defor97[forest$prot_b80 == 0])

res <- 1000
rast_template <- raster(extent(forest), res=res, crs=st_crs(forest))

ras <- rasterize(forest, rast_template, field="defor97", fun=function(x,...) mean(x, na.rm = TRUE))
ras[is.na(ras)] <- 0
plot(ras)

# centr <- st_coordinates(st_centroid(districts))
forest$luc_high <- as.numeric(forest$luc1 == 1 | forest$luc2 == 1 | forest$luc3 == 1)
ras_cov1 <- rasterize(forest, rast_template, field="luc_high", fun=function(x,...) mean(x, na.rm = TRUE))
ras_cov1[is.na(ras_cov1)] <- 0
ras_cov2 <- rasterize(forest, rast_template, field="dist_rd_69", fun=function(x,...) mean(x, na.rm = TRUE))
ras_cov2[is.na(ras_cov2)] <- 0
ras_cov3 <- rasterize(forest, rast_template, field="dist_mcity", fun=function(x,...) mean(x, na.rm = TRUE))
ras_cov3[is.na(ras_cov3)] <- 0
values_cov1 <- getValues(ras_cov1)
values_cov2 <- getValues(ras_cov2)
values_cov3 <- getValues(ras_cov3)

ras_Z <- rasterToPolygons(boundary_tiles)
com_index <- which(!is.na(getValues(boundary_tiles)) == TRUE) # select tiles with values
cov1 <- cov2 <- cov3 <- c()
num_tiles <- c()
for (i in 1:length(ras_Z)){
  center <- ras_Z[i, ]
  one_buffer_raster <- rasterize(center, ras)
  getValues(one_buffer_raster)
  one_buffer_raster_values <- getValues(one_buffer_raster)
  coords <- coordinates(one_buffer_raster)[(!is.na(one_buffer_raster_values) & 
                                              !is.nan(one_buffer_raster_values) & 
                                              one_buffer_raster_values != 0), ]
  if (is.na(coords[1])){
    boundary_tiles[com_index[i]] <- NA
  }else{
    num_tiles <- c(num_tiles, nrow(coords))
    grid_index1 <- c(na.omit(cellFromXY(ras_cov1, coords)))
    grid_index2 <- c(na.omit(cellFromXY(ras_cov2, coords)))
    grid_index3 <- c(na.omit(cellFromXY(ras_cov3, coords)))
    cov1 <- c(cov1, mean(values_cov1[grid_index1]))
    cov2 <- c(cov2, mean(values_cov2[grid_index2]))
    cov3 <- c(cov3, mean(values_cov3[grid_index3]))
  }
}

outcome <- getValues(ras)
table(outcome)
Yxy <- coordinates(ras)
Yindex <- 1:nrow(Yxy)
Ydata <- data.frame(Yindex,cbind(outcome, Yxy))
rm(Yxy)


Zdata <- data.frame("treatment" = getValues(boundary_tiles))
Zdata <- as.data.frame(Zdata[!is.na(Zdata$treatment), ])
names(Zdata)[1] <- "treatment"
Zdata$cov1 <- cov1
Zdata$cov2 <- cov2
Zdata$cov3 <- cov3

pscore_formula <- "treatment ~ cov1 + cov2 + cov3"
pscore_fit <- glm(as.formula(pscore_formula), Zdata, family = "binomial")
pscore_est <- predict(pscore_fit, type = "response")
Zdata$prob_treatment <- pscore_est

pdf(file="F2015-forest-pscore.pdf", height=8, width=8)
par(mfrow=c(1,1))
par(mar = c(5, 5, 4, 3))
plot(density(Zdata$prob_treatment),
     cex.lab=2,
     cex.main=2,
     cex.axis = 2,
     xlab = "Estimated propensity score",
     main = "Density of propensity score estimates")
dev.off()

ras_Z <- rasterToPolygons(boundary_tiles)
ras_Z <- st_as_sf(ras_Z)

# Zdata <- as.data.frame(Zdata[, c(10, 11, 9)])
# names(Zdata) <- c("x_coord", "y_coord", "treatment")

ras <- ras
Zdata <- Zdata
x_coord_Z <- "x_coord"
y_coord_Z <- "y_coord"
treatment <- "treatment"
prob_treatment <- "prob_treatment"
numpts <- 1000
evalpts <- 1
only.unique <- 0
per.se <- 1
blockvar <- NULL
clustvar <- NULL
conley.se <- 1
alpha <- 0.05
cutoff <- 15000
nPerms <- 500
# dist.metric <- "Geodesic"
dist.metric <- "Euclidean"
smooth <- 1
bw <- NULL
bw_debias <- NULL
bias_correction <- 1
kernel <- "tri"
edf <- 0
smooth.conley.se <- 1
cType <- "donut"
covs <- NULL

# len <- 500/111
# dVec <- seq(from=0, to=len, by=20/111)
len <- 20000
dVec <- seq(from=0, to=len, by=500)

# maps
tr_colors <- rep("goldenrod1", nrow(Zdata))
tr_colors[Zdata$treatment == 0] <- "darkturquoise"
pdf(file="F2015-forest-map-polygon.pdf", height=8, width=8)
par(mfrow=c(1,1))
plot(ras, bty="n", box=FALSE, col = c("white", "darkgoldenrod4"), legend=FALSE, 
     cex.lab=2,
     cex.main=2,
     cex.axis = 2)
lines(ras_Z, col=tr_colors, lwd = 2)
legend("bottomleft", col = c("darkturquoise", "goldenrod1", "darkgoldenrod4"),
       legend = c("Control", "Treated", "Deforestation"), 
       lty=c(1, 1, NA), pch=c(NA, NA, 19), cex=1.5, pt.cex = c(.8, .8, .4), bty = "n")
dev.off()

pdf(file="F2015-forest-map-polygon-circles-treated.pdf", height=8, width=8)
par(mfrow=c(1,1))
plot(ras, bty="n", box=FALSE, col = c("white", "darkgoldenrod4"), legend=FALSE, 
     cex.lab=2,
     cex.main=2,
     cex.axis = 2)
lines(ras_Z, col=tr_colors, lwd = 2)
for (i in which(Zdata$treatment == 1)){
  buffer_treated <- buffer(as_Spatial(ras_Z[i, ]), width=10000)
  lines(buffer_treated, col="goldenrod1", lwd = 2, lty = 2)
}
legend("bottomleft", col = c("darkturquoise", "goldenrod1", "darkgoldenrod4"),
       legend = c("Control", "Treated", "Deforestation"), 
       lty=c(1, 1, NA), pch=c(NA, NA, 19), cex=1.5, pt.cex = c(.8, .8, .4), bty = "n")
dev.off()

pdf(file="F2015-forest-map-polygon-circles-control.pdf", height=8, width=8)
par(mfrow=c(1,1))
plot(ras, bty="n", box=FALSE, col = c("white", "darkgoldenrod4"), legend=FALSE,
     cex.lab=2,
     cex.main=2,
     cex.axis = 2)
lines(ras_Z, col=tr_colors, lwd = 2)
for (i in which(Zdata$treatment == 0)){
  buffer_control <- buffer(as_Spatial(ras_Z[i, ]), width=10000)
  lines(buffer_control, col="darkturquoise", lwd = 2, lty = 2)
}
legend("bottomleft", col = c("darkturquoise", "goldenrod1", "darkgoldenrod4"),
       legend = c("Control", "Treated", "Deforestation"), 
       lty=c(1, 1, NA), pch=c(NA, NA, 19), cex=1.5, pt.cex = c(.8, .8, .4), bty = "n")
dev.off()


# analysis
result.list <- SpatialEffect(ras = ras, Zdata = Zdata, x_coord_Z = x_coord_Z, y_coord_Z = y_coord_Z, ras_Z = NULL, 
                             treatment = treatment, dVec = dVec[-1], numpts = 10000, only.unique = only.unique, smooth = smooth, 
                             per.se = per.se, conley.se = conley.se, cutoff = cutoff, alpha = alpha, edf = edf,
                             nPerms = nPerms, smooth.conley.se = smooth.conley.se, dist.metric = dist.metric, cType = "edge", kernel = kernel)

result.list <- SpatialEffect(ras = ras, Zdata = Zdata, x_coord_Z = x_coord_Z, y_coord_Z = y_coord_Z, ras_Z = ras_Z, 
                             treatment = treatment, dVec = dVec, numpts = numpts, only.unique = only.unique, smooth = smooth, 
                             per.se = per.se, conley.se = conley.se, cutoff = cutoff, alpha = alpha, edf = edf, bw = NULL, bw_debias = NULL,
                             nPerms = nPerms, smooth.conley.se = smooth.conley.se, dist.metric = dist.metric, cType = "donut", kernel = kernel)

result.list <- SpatialEffect(ras = ras, Zdata = Zdata, ras_Z = ras_Z, alpha = alpha,
                             treatment = treatment, dVec = dVec, numpts = numpts, cutoff = cutoff, dist.metric = dist.metric, smooth = 1, kernel = kernel)


summary(result.list, dVec.range = c(1000, 5000))

load("/Users/yewang/Documents/GitHub/spatial/data/forest/Ferraro_etal_2015.RData")
AMR_est <- result.list[["AMR_est"]]
Per.CI <- result.list[["Per.CI"]]
Conley.SE <- result.list[["Conley.SE"]]
Conley.CI <- result.list[["Conley.CI"]]
AMR_est_smoothed <- result.list[["AMR_est_smoothed"]]
smoothed.Conley.SE <- result.list[["smoothed.Conley.SE"]]
smoothed.Conley.CI <- result.list[["smoothed.Conley.CI"]]
smoothed.Conley.CB <- result.list[["smoothed.Conley.CB"]]

save(result.list, file = "/Users/yewang/Documents/GitHub/spatial/data/forest/Ferraro_etal_2015.RData")

test.result <- SpatialEffectTest(result.list, test.range = c(1000, 5000), smooth = 0)
test.result <- SpatialEffectTest(result.list, test.range = c(1, 5), smooth = 0, alpha = 0.1)

dVec_real <- dVec / 1000


pdf(file="Ferraro_etal_AMR_polygon_unsmoothed.pdf", height=8, width=8)
par(mar = c(5, 5, 5, 5))
plot(	AMR_est[, 2] ~ dVec_real,
      ylab="Estimates",
      xlab="Distance (kilometers)",
      type="l",
      ylim=c(-0.04, 0.04),
      main = "Hajek Estimator with Its CIs",
      cex.lab=2,
      cex.main=2,
      cex.axis = 2,
      lwd = 3)
points(Conley.CI[, 1]~dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
points(Conley.CI[, 2]~dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
points(Per.CI[, 1]~dVec_real, type = "l", col = "blue", lty = "dotted", lwd = 3)
points(Per.CI[, 2]~dVec_real, type = "l", col = "blue", lty = "dotted", lwd = 3)
abline(h = 0, lty = 2, lwd = 3, col = "gray")
legend("topright", col = c("black", "black", "blue"), lty=c(1,3,3), lwd = c(3, 3, 3),
       legend = c("AME estimates", "Spatial HAC 95% CI", "Quantiles under sharp null"), cex=2, bty = "n")
dev.off()

pdf(file="Ferraro_etal_AMR_polygon_smoothed.pdf", height=8, width=8)
par(mar = c(5, 5, 5, 5))
plot(	AMR_est_smoothed[, 2] ~ dVec_real,
      ylab="Estimates",
      xlab="Distance (kilometers)",
      type="l",
      ylim=c(-0.04, 0.04),
      main = "Kernel Regression Estimator with Its CIs",
      cex.lab=2,
      cex.main=2,
      cex.axis = 2,
      lwd = 3)
points(smoothed.Conley.CI[, 1]~dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
points(smoothed.Conley.CI[, 2]~dVec_real, type = "l", col = "black", lty = "dotted", lwd = 3)
abline(h = 0, lty = 2, lwd = 3, col = "gray")
legend("topright", col = c("black", "black"), lty=c(1,3), lwd = c(3, 3),
       legend = c("AME estimates", "Spatial HAC 95% CI"), cex=2, bty = "n")
dev.off()


ha_ests <- AMR_est[, 2]
ha_CI95_upper <- Conley.CI[, 2]
ha_CI95_lower <- Conley.CI[, 1]
ha_CI90_upper <- AMR_est[, 2] + qnorm(0.95) * Conley.SE
ha_CI90_lower <- AMR_est[, 2] - qnorm(0.95) * Conley.SE

names <- dVec_real

ha.data <- data.frame("estimates" = ha_ests * 100, "CI95_u" = ha_CI95_upper * 100, 
                     "CI95_l" = ha_CI95_lower * 100, "CI90_u" = ha_CI90_upper * 100, 
                     "CI90_l" = ha_CI90_lower * 100, "names" = names)
ha.data$names <- factor(ha.data$names, levels = rev(rev(ha.data$names)))

colors <- "red"

require(ggplot2)
coefs_p <- ggplot(ha.data[seq(1, nrow(ha.data) - 1, 2), ]) + 
  geom_pointrange(aes(x=names, y=estimates, ymin=CI95_l, ymax=CI95_u), size = 0.3, colour=colors, lty = "dashed") + 
  geom_pointrange(aes(x=names, y=estimates, ymin=CI90_l, ymax=CI90_u), size = 0.6, colour=colors, lty = "solid") + 
  theme_bw() + geom_hline(aes(yintercept = 0), colour = "black", lty=2) + xlab('Distance (kilometers)') +
  scale_linetype_manual(labels = c("90% CI", "95% CI"), values=c(1, 2)) + ylim(-5, 2) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = "bottom", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + ylab("Estimates")
ggsave(paste0("F2015_graphs/ha_coefs_reproj_90CI.pdf"), coefs_p, device = "pdf", height = 8, width = 8) 


smooth_ests <- AMR_est_smoothed[, 2]
smooth_CI95_upper <- smoothed.Conley.CI[, 2]
smooth_CI95_lower <- smoothed.Conley.CI[, 1]
smooth_CI90_upper <- AMR_est_smoothed[, 2] + qnorm(0.95) * smoothed.Conley.SE
smooth_CI90_lower <- AMR_est_smoothed[, 2] - qnorm(0.95) * smoothed.Conley.SE
colors <- "red"

names <- dVec_real

smooth.data <- data.frame("estimates" = smooth_ests * 100, "CI95_u" = smooth_CI95_upper * 100, 
                      "CI95_l" = smooth_CI95_lower * 100, "CI90_u" = smooth_CI90_upper * 100, 
                      "CI90_l" = smooth_CI90_lower * 100, "names" = names)
smooth.data$names <- factor(smooth.data$names, levels = rev(rev(smooth.data$names)))


require(ggplot2)
coefs_p <- ggplot(smooth.data[seq(1, nrow(smooth.data) - 1, 2), ]) + 
  geom_pointrange(aes(x=names, y=estimates, ymin=CI95_l, ymax=CI95_u), size = 0.3, colour=colors, lty = "dashed") + 
  geom_pointrange(aes(x=names, y=estimates, ymin=CI90_l, ymax=CI90_u), size = 0.6, colour=colors, lty = "solid") + 
  theme_bw() + geom_hline(aes(yintercept = 0), colour = "black", lty=2) + xlab('Distance (kilometers)') +
  scale_linetype_manual(labels = c("90% CI", "95% CI"), values=c(1, 2)) + ylim(-5, 2) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position = c(18, -0.015), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + ylab("Estimates")
ggsave(paste0("F2015_graphs/smooth_coefs_reproj_90CI.pdf"), coefs_p, device = "pdf", height = 8, width = 8) 



