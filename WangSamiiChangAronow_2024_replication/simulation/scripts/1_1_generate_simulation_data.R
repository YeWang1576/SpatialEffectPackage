# simulation results in Wang, Samii, Chang, and Aronow (2024)
rm(list=ls())

##############################################################################
###################Import Libraries###########################################
##############################################################################
library(foreign)
library(sp)
library(fields)
library(raster)
library(sf)
library(Rcpp)
library(RcppArmadillo)
library(pbapply)



#set.seed(2024)
set.seed(2025)

##############################################################################
###################Set Root Directory#########################################
##############################################################################
root_dir = "C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/SpatialReplication/"
setwd(root_dir)
data_path="simulation/data_new/"


##############################################################################
#################Number of Simulations #######################################
##############################################################################
nsim=2000


##############################################################################
################Import Functions##############################################
##############################################################################
sourceCpp("./SpatialEffect/src/DistanceCalculation2.cpp",verbose=0)
sourceCpp("./SpatialEffect/src/calculate_tau.cpp",verbose=0)
sourceCpp("./SpatialEffect/src/calculate_tau_separate.cpp",verbose=0)
sourceCpp("./SpatialEffect/src/effect_nonmo.cpp",verbose=0)
sourceCpp("./SpatialEffect/src/effect_interactive.cpp",verbose=0)
source("./simulation/scripts/5_functions_generate_simulation_data.R")
source("./simulation/scripts/6_functions.R")


##############################################################################
#################Today's Date#################################################
##############################################################################
today_date=format(Sys.Date(), format="%b_%d")


##############################################################################
#################Simulation: Polygon Interventions (paired) ##################
##############################################################################

start.time=proc.time()
# nonmono effect + YZratio = 10 + paired polygon 
# On my laptop it took 2397.17 seconds
TYPE <- "nonmono"
sample_sizes <- c(80, 100, 120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=1)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
print(proc.time()-start.time)


start.time=proc.time()
# interactive effect + YZratio = 10 + paired polygon 
# On my laptop it took 3125.52 seconds
TYPE <- "interactive"
sample_sizes <- c(80, 100, 120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=1)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
 print(proc.time()-start.time)

start.time=proc.time()
# null effect + YZratio = 10 + paired polygon 
# On my laptop it took 3107.11 seconds
TYPE <- "null"
sample_sizes <- c(80, 100, 120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=1)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
print(proc.time()-start.time)


##############################################################################
#################Simulation: Point Interventions (paired) ####################
##############################################################################

start.time=proc.time()
# nonmono effect + YZratio = 10 + paired points
# On my laptop it took 150.68 seconds
TYPE <- "nonmono"
sample_sizes <- c(80, 100, 120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=1)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
print(proc.time()-start.time)

start.time=proc.time()
# interactive effect + YZratio = 10 + paired points
# On my laptop it took 142.52 seconds
TYPE <- "interactive"
sample_sizes <- c(80, 100, 120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=1)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
print(proc.time()-start.time)

start.time=proc.time()
# null effect + YZratio = 10 + paired points
# On my laptop it took 169.28 seconds
TYPE <- "null"
sample_sizes <- c(80, 100, 120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=1)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
print(proc.time()-start.time)


##############################################################################
#################Simulation: Point Interventions (for ploting)################
##############################################################################

start.time=proc.time()
# nonmono effect + YZratio = 10 + paired points
# On my laptop it took 30 seconds
TYPE <- "nonmono"
sample_sizes <- c(120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator_for_ploting, type=TYPE,nsim=5000,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=3)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_point_ploting.RData"))
print(proc.time()-start.time)

start.time=proc.time()
# interactive effect + YZratio = 10 + paired points
# On my laptop it took 30 seconds
TYPE <- "interactive"
sample_sizes <- c(120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator_for_ploting, type=TYPE,nsim=1000,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=3)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_point_ploting.RData"))
print(proc.time()-start.time)

start.time=proc.time()
# nonmono effect + YZratio = 10 + paired points
# On my laptop it took 2000 seconds
TYPE <- "nonmono"
sample_sizes <- c(120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator_for_ploting, type=TYPE,nsim=1000,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=3)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon_ploting.RData"))
print(proc.time()-start.time)

start.time=proc.time()
# interactive effect + YZratio = 10 + paired points
# On my laptop it took 2197.68 seconds
TYPE <- "interactive"
sample_sizes <- c(120)
YZratio=10
data_list=lapply(sample_sizes,dataGenerator_for_ploting, type=TYPE,nsim=1000,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=TRUE,polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=3)
save(data_list, file = paste0(data_path,"Paired_simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon_ploting.RData"))
print(proc.time()-start.time)

# ##############################################################################
# #################Simulation: Point Interventions #############################
# ##############################################################################
# 
# # start.time=proc.time()
# # generate data with non-monotonic effect function and YZratio=3
# # This should take at least around 8 mins
# TYPE <- "nonmono"
# sample_sizes <- c(80)
# YZratio=3
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=FALSE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# # print(proc.time()-start.time)
# 
# 
# # start.time=proc.time()
# # generate data with interactive effect function and YZratio=3
# # This should take at least around 8 mins
# TYPE <- "interactive"
# sample_sizes <- c(80)
# YZratio=3
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=FALSE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# # print(proc.time()-start.time)
# 
# # start.time=proc.time()
# # generate data with non-monotonic effect function and YZratio=3
# # This should take at least around 8 mins
# TYPE <- "nonmono"
# sample_sizes <- c(80)
# YZratio=5
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=FALSE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# # print(proc.time()-start.time)
# 
# 
# # start.time=proc.time()
# # generate data with interactive effect function and YZratio=3
# # This should take at least around 8 mins
# TYPE <- "interactive"
# sample_sizes <- c(80)
# YZratio=5
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=FALSE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# # print(proc.time()-start.time)
# 
# # start.time=proc.time()
# # generate data with non-monotonic effect function and YZratio=10
# # This should take at least around 30 seconds
# 
# TYPE <- "nonmono"
# sample_sizes <- c(40,60,80, 100, 120)
# YZratio=10
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=FALSE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# # print(proc.time()-start.time)
# 
# 
# # start.time=proc.time()
# # generate data with interactive effect function and YZratio=10
# # This should take at least around 30 seconds
# 
# TYPE <- "interactive"
# sample_sizes <- c(40,60,80, 100, 120)
# YZratio=10
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=FALSE,polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# # print(proc.time()-start.time)
# 
# 
# 
# 
# 
# ##############################################################################
# #################Simulation: Polygon Interventions ###########################
# ##############################################################################
# 
# start.time=proc.time()
# # generate data with non-monotonic effect function and YZratio=10 (polygon intervention)
# # This should take at least around 30 seconds
# TYPE <- "nonmono"
# sample_sizes <- c(40,60,80, 100, 120)
# YZratio=10
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE,pair=FALSE, polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# print(proc.time()-start.time)
# 
# start.time=proc.time()
# # generate data with non-monotonic effect function and YZratio=10 (polygon intervention)
# # This should take at least around 30 seconds
# TYPE <- "interactive"
# sample_sizes <- c(40,60,80, 100, 120)
# YZratio=10
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE,pair=FALSE, polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# print(proc.time()-start.time)
# 
# # start.time=proc.time()
# # generate data with non-monotonic effect function and YZratio=3
# # This should take at least around 8 mins
# TYPE <- "nonmono"
# sample_sizes <- c(80)
# YZratio=5
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=FALSE,polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# # print(proc.time()-start.time)
# 
# 
# # start.time=proc.time()
# # generate data with interactive effect function and YZratio=3
# # This should take at least around 8 mins
# TYPE <- "interactive"
# sample_sizes <- c(80)
# YZratio=5
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=FALSE,polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# # print(proc.time()-start.time)
# 
# # start.time=proc.time()
# # generate data with non-monotonic effect function and YZratio=3
# # This should take at least around 8 mins
# TYPE <- "nonmono"
# sample_sizes <- c(80)
# YZratio=3
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, pair=FALSE,polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# # print(proc.time()-start.time)
# 
# 
# # start.time=proc.time()
# # generate data with interactive effect function and YZratio=3
# # This should take at least around 8 mins
# TYPE <- "interactive"
# sample_sizes <- c(80)
# YZratio=3
# data_list=lapply(sample_sizes,dataGenerator, type=TYPE,nsim=nsim,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE,pair=FALSE, polygon=1, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
# save(data_list, file = paste0(data_path,"simData_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# # print(proc.time()-start.time)
# 
