\name{SpatialEffect}
\alias{SpatialEffect}
\title{SpatialEffect}
\description{The function estimate the average marginalized response (AMR) 
curve defined in Aronow, Samii and Wang (2019) and its standard error.}
\usage{
  SpatialEffect(ras, Ydata, outcomme = "outcome", x_coord_Y = "x_coord_Y", 
  y_coord_Y = "y_coord_Y", Zdata, x_coord_Z = "x_coord_Z", y_coord_Z = "y_coord_Z", 
  treatment = "treatment", dVec, dist.metric = "Euclidean", numpts = NULL, 
  evalpts = 1, only.unique = 0, per.se = 1, blockvar = NULL, clustvar = NULL, 
  conley.se = 1, alpha = 0.05, cutoff = 0, m = 2, lambda = 0.02, nPerms = 1000)
}
\arguments{
\item{ras}{a raster object from which outcome values are derived.}
\item{Ydata}{a dataset containing values of outcome units. When it is set to be NULL, kriging                 will be used to extrapolate the value of each outcome unit.}
\item{outcome}{a character string indicating the name of the outcome variable in Ydata/Zdata.}
\item{x_coord_Y}{a character string indicating the name of the variable that represents the first                  coordinate in Ydata.}
\item{y_coord_Y}{a character string indicating the name of the variable that represents the second                 coordinate in Ydata.}
\item{Zdata}{a dataset containing the coordinates of intervention units.}
\item{x_coord_Z}{a character string indicating the name of the variable that represents the first                  coordinate in Zdata.}
\item{y_coord_Z}{a character string indicating the name of the variable that represents the second                 coordinate in Zdata.}
\item{treatment}{a character string indicating the name of the treatment variable.}
\item{dVec}{a vector of distances at which to estimate the MIR effects.}
\item{dist.metric}{the metric used to calculate distance, can be either "Euclidean" (default) or "Geodesic".} 
\item{numpts}{the number of sampling points around each circle we take to calculate the average                 outcome.Its value is set to be proportional to the circle's radius by default}
\item{evalpts}{the density of sampling points around each circle we take to calculate the average                 outcome.It is not used unless numpts = NULL.}
\item{only.unique}{when a cell is taken as a sampling point for multiple times, whether we                           calculate its value only once.}
\item{smooth}{whether the estimates at each distance value should be smoothed via local polynomical regression.}
\item{bw}{the bandwidth for the local polynomical regression.}
\item{per.se}{a boolean variable indicating whether to calculate the permutation confidence                     interval.}
\item{blockvar}{a positive integer vector, with unique values indicating different blocks}
\item{clustervar}{a positive integer vector, with unique values indicating different clusters}
\item{conley.se}{a boolean variable indicating whether to calculate the Conley standard error and                  confidence interval.}
\item{alpha}{significant level, with default 0.05}
\item{cutoff}{the bandwidth distance (in the coordinate space of the underlying data) for the                   Conley spatial HAC standard error.}
\item{m}{an integer, the degree of polynomial used in kriging estimation.}
\item{lambda}{a real number, the degree of overfitting used in kriging estimation.}
\item{nPerms}{the number of permutation used when calculating the permuation standard error.}
}

\value{
  a list of AMR estimates and (if required) its standard error
  \item{AMR_est}{the vector of AMR estimates at each distance.}
  \item{Per.CI}{confidence interval estimates from permuation}
  \item{Conley.SE}{Conley standard error estimates}
  \item{Conley.CI}{confidence interval estimates based on the Conley standard error}
  \item{AMR_est_smoothed}{the vector of AMR estimates at each distance from kernel regression.}
}
\examples{
data(sim_Zdata)
data(sim_ras)
len <- 10
dVec <- seq(from=.5, to=len, by=.25)

result.list <- SpatialEffect(ras = ras, Zdata = sim_Zdata, x_coord_Z = "x", y_coord_Z = 
"y", treatment = "treatment", dVec = dVec, numpts = NULL, evalpts = 1, only.unique = 0, 
per.se = 1, conley.se = 1, alpha = 0.05, cutoff = 3, m = 2, lambda = 0.02, nPerms= 1000)

AMR_est <- result.list[[1]]
Per.CI <- result.list[[2]]
Conley.SE <- result.list[[3]]
Conley.CI <- result.list[[4]]
plot(AMR_est,
     ylab="AMR estimates",
     xlab="Distance",
     type="l",
     ylim=c(-1,1.5))
points(Conley.CI[, 1]~dVec, type = "l", col = "red")
points(Conley.CI[, 2]~dVec, type = "l", col = "red")
points(Per.CI[, 1]~dVec, type = "l", col = "blue", lty = "dotted")
points(Per.CI[, 2]~dVec, type = "l", col = "blue", lty = "dotted")

data(worms)
dVec <- seq(.01,.4,.01)

result.list <- SpatialEffect(outcome = "any_ics99", Zdata = wormsprep, 
x_coord_Z = "east", y_coord_Z = "north", treatment = "wgrp1", dVec = dVec, 
numpts = 10, evalpts = 1, only.unique = 0, cutoff = 0.1, per.se = 1, 
conley.se = 1, alpha = 0.05, m = 2, lambda = 0.02, nPerms = 1000)

AMR_est <- result.list[[1]]
Per.CI <- result.list[[2]]
Conley.SE <- result.list[[3]]
Conley.CI <- result.list[[4]]

plot(AMR_est,
     ylab="AMR estimates",
     xlab="Distance",
     type="l",
     ylim=c(-.3, 0.2))
points(Conley.CI[, 1]~dVec, type = "l", col = "red")
points(Conley.CI[, 2]~dVec, type = "l", col = "red")
points(Per.CI[, 1]~dVec, type = "l", col = "blue", lty = "dotted")
points(Per.CI[, 2]~dVec, type = "l", col = "blue", lty = "dotted")
}