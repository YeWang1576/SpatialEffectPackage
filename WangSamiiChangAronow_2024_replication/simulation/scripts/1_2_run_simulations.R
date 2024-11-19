# simulation results in Wang, Samii, Chang, and Aronow (2024)
rm(list=ls())

library(Rcpp)
library(Matrix)
library(RcppArmadillo)
library(plotrix)
library(foreign)
library(sp)
library(fields) #for kriging
library(ri)
library(sf)
library(raster) #used to create rasters for simulations
library(units) #used to transform some data objects to standard matriices
library(ggplot2) #ggplots
library(data.table) #for resahpe data for ploting
##############################################################################
##################If ri package is not installed##############################
##############################################################################
# #ri package no longer available from repository
# install.packages('./simulation/ri_0.9.tar.gz',type='source')


##############################################################################
##################Set Root Directories########################################
##############################################################################
root_dir = "C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/SpatialReplication/"
setwd(root_dir)
data_path="simulation/data_new/"


##############################################################################
##################Import Functions############################################
##############################################################################
sourceCpp("./SpatialEffect/src/DistanceCalculation.cpp")
sourceCpp("./SpatialEffect/src/ConleySE_kernel_matrix.cpp")
sourceCpp("./SpatialEffect/src/Conley2.cpp")
source("./simulation/scripts/5_functions_run_simulations.R")
source("./simulation/scripts/6_functions.R")
source("./simulation/scripts/6_functions2.R")
source("./simulation/scripts/6_functions3.R")


##############################################################################
#################Today's Date#################################################
##############################################################################
today_date=format(Sys.Date(), format="%b_%d")


##############################################################################
##################Simulation##################################################
##############################################################################



##############################################################################
##################Point Intervention (Pair)###################################
##############################################################################


# simulation with non-monotonic effect function and YZratio=10
start.time=proc.time()
sim1_raw_data_file = "simulation/data_new/Paired_simData_nonmono_10_Jun_06_points.RData"
TYPE='nonmono'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
cutoff = 2*6
simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=0)
save(simulation_result, file = paste0(data_path,"Paired_simulation_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
print(proc.time()-start.time)


# simulation with non-monotonic effect function and YZratio=10
start.time=proc.time()
sim1_raw_data_file = "simulation/data_new/Paired_simData_interactive_10_Jun_06_points.RData"
TYPE='interactive'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
cutoff = NULL
simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=0)
save(simulation_result, file = paste0(data_path,"Paired_simulation_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
print(proc.time()-start.time)

##############################################################################
##################Polygon Intervention (Paired) ##############################
##############################################################################
# simulation with non-monotonic effect function and YZratio=10
start.time=proc.time()
sim1_raw_data_file = "simulation/data_new/Paired_simData_nonmono_10_Jun_09_polygon.RData"
TYPE='nonmono'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
cutoff = 2*6
simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=1)
save(simulation_result, file = paste0(data_path,"Paired_simulation_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
print(proc.time()-start.time)


# simulation with interactive effect function and YZratio=10
start.time=proc.time()
sim1_raw_data_file = "simulation/data_new/Paired_simData_interactive_10_Jun_09_polygon.RData"
TYPE='interactive'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
cutoff = NULL
simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=1)
save(simulation_result, file = paste0(data_path,"Paired_simulation_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
print(proc.time()-start.time)




# ##############################################################################
# ##################Point Intervention##########################################
# ##############################################################################
# 
# 
# # simulation with non-monotonic effect function and YZratio=10
# start.time=proc.time()
# sim1_raw_data_file = "simulation/data_new/simData_nonmono_10_May_23_points.RData"
# TYPE='nonmono'
# YZratio=10
# sample_sizes <-c(40,60,80, 100, 120)
# dVec <- seq(from=.5, to=10, by=.25)
# bw <- 1
# cutoff = 2*6
# simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=0)
# save(simulation_result, file = paste0(data_path,"simulation_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# print(proc.time()-start.time)
# 
# # simulation with interactive effect function and YZratio=10
# start.time=proc.time()
# sim1_raw_data_file = "simulation/data_new/simData_interactive_10_May_23_points.RData"
# TYPE='interactive'
# YZratio=10
# sample_sizes <-c(40,60,80, 100, 120)
# dVec <- seq(from=.5, to=10, by=.25)
# bw <- 1
# cutoff = NULL
# simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=0)
# save(simulation_result, file = paste0(data_path,"simulation_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# print(proc.time()-start.time)
# 
# # simulation with non-monotonic effect function and YZratio=5
# start.time=proc.time()
# sim1_raw_data_file = "simulation/data_new/simData_nonmono_5_May_23_points.RData"
# TYPE='nonmono'
# YZratio=5
# sample_sizes <-c(80)
# dVec <- seq(from=.5, to=10, by=.25)
# bw <- 1
# cutoff = 2*6
# simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=0)
# save(simulation_result, file = paste0(data_path,"simulation_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# print(proc.time()-start.time)
# 
# # simulation with interactive effect function and YZratio=5
# start.time=proc.time()
# sim1_raw_data_file = "simulation/data_new/simData_interactive_5_May_23_points.RData"
# TYPE='interactive'
# YZratio=5
# sample_sizes <-c(80)
# dVec <- seq(from=.5, to=10, by=.25)
# bw <- 1
# cutoff = NULL
# simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=0)
# save(simulation_result, file = paste0(data_path,"simulation_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
# print(proc.time()-start.time)
# 
# ##############################################################################
# ##################Polygon Intervention########################################
# ##############################################################################
# # simulation with non-monotonic effect function and YZratio=10
# start.time=proc.time()
# sim1_raw_data_file = "simulation/data_new/simData_nonmono_10_May_23_polygon.RData"
# TYPE='nonmono'
# YZratio=10
# sample_sizes <-c(40,60,80, 100, 120)
# dVec <- seq(from=.5, to=10, by=.25)
# bw <- 1
# cutoff = 2*6
# simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=1)
# save(simulation_result, file = paste0(data_path,"simulation_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# print(proc.time()-start.time)
# 
# 
# # simulation with interactive effect function and YZratio=10
# start.time=proc.time()
# sim1_raw_data_file = "simulation/data_new/simData_interactive_10_May_23_polygon.RData"
# TYPE='interactive'
# YZratio=10
# sample_sizes <-c(40,60,80, 100, 120)
# dVec <- seq(from=.5, to=10, by=.25)
# bw <- 1
# cutoff = NULL
# simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=1)
# save(simulation_result, file = paste0(data_path,"simulation_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# print(proc.time()-start.time)
# 
# 
# start.time=proc.time()
# sim1_raw_data_file = "simulation/data_new/simData_nonmono_5_May_23_polygon.RData"
# TYPE='nonmono'
# YZratio=5
# sample_sizes <-c(80)
# dVec <- seq(from=.5, to=10, by=.25)
# bw <- 1
# cutoff = 2*6
# simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=1)
# save(simulation_result, file = paste0(data_path,"simulation_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# print(proc.time()-start.time)
# 
# # simulation with interactive effect function and YZratio=10
# start.time=proc.time()
# sim1_raw_data_file = "simulation/data_new/simData_interactive_5_May_23_polygon.RData"
# TYPE='interactive'
# YZratio=5
# sample_sizes <-c(80)
# dVec <- seq(from=.5, to=10, by=.25)
# bw <- 1
# cutoff = NULL
# simulation_result=run_simulation(sim1_raw_data_file,cutoff,bw,sample_sizes,polygon=1)
# save(simulation_result, file = paste0(data_path,"simulation_", TYPE, "_" ,YZratio,'_',today_date,"_polygon.RData"))
# print(proc.time()-start.time)
# 
