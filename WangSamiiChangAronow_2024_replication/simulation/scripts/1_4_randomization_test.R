# simulation results in Wang, Samii, Chang, and Aronow (2024)
rm(list=ls())

##############################################################################
##################If ri package is not installed##############################
##############################################################################
# #ri package no longer available from repository
# install.packages('./simulation/ri_0.9.tar.gz',type='source')


##############################################################################
##################Set Root Directories########################################
##############################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data_path="../data_new/"

library(foreign)
library(sp)
library(fields)
library(raster)
library(sf)
library(Rcpp)
library(RcppArmadillo)
library(pbapply)

##############################################################################
##################Import Functions############################################
##############################################################################
source("5_functions_randomization_test.R")
source("6_functions.R")



##############################################################################
#################Today's Date#################################################
##############################################################################
today_date=format(Sys.Date(), format="%b_%d")

# simulation with non-monotonic effect function and YZratio=10
start.time=proc.time()
sim1_raw_data_file = "../data_new/Paired_simData_null_10_Jun_06_points.RData"
TYPE='null'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
simulation_result=run_simulation_randomization_test(sim1_raw_data_file,bw,sample_sizes,polygon=0,nperm=1000)
save(simulation_result, file = paste0(data_path,"Paired_randomization_test_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
print(proc.time()-start.time)

# simulation with non-monotonic effect function and YZratio=10
start.time=proc.time()
sim1_raw_data_file = "../data_new/Paired_simData_nonmono_10_Jun_06_points.RData"
TYPE='nonmono'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
simulation_result=run_simulation_randomization_test(sim1_raw_data_file,bw,sample_sizes,polygon=0,nperm=1000)
save(simulation_result, file = paste0(data_path,"Paired_randomization_test_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
print(proc.time()-start.time)


# simulation with non-monotonic effect function and YZratio=10
#This took 177722 seconds on my laptop
start.time=proc.time()
sim1_raw_data_file = "../data_new/Paired_simData_interactive_10_Jun_06_points.RData"
TYPE='interactive'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
simulation_result=run_simulation_randomization_test(sim1_raw_data_file,bw,sample_sizes,polygon=0,nperm=1000)
save(simulation_result, file = paste0(data_path,"Paired_randomization_test_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
print(proc.time()-start.time)




##############################################################################
#################Randomization Test with sup norm#############################
##############################################################################
today_date=format(Sys.Date(), format="%b_%d")

# simulation with non-monotonic effect function and YZratio=10
#This one took 2735 seconds on my laptop
start.time=proc.time()
sim1_raw_data_file = "../data_new/Paired_simData_null_10_Jun_06_points.RData"
TYPE='null'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
simulation_result=run_simulation_randomization_test_sup(sim1_raw_data_file,bw,sample_sizes,polygon=0,nperm=1000)
save(simulation_result, file = paste0(data_path,"Paired_randomization_test_sup_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
time_rand_sup_null=proc.time()-start.time
print(time_rand_sup_null)

# simulation with non-monotonic effect function and YZratio=10
#This one took 2800 seconds on my laptop
start.time=proc.time()
sim1_raw_data_file = "../data_new/Paired_simData_nonmono_10_Jun_06_points.RData"
TYPE='nonmono'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
simulation_result=run_simulation_randomization_test_sup(sim1_raw_data_file,bw,sample_sizes,polygon=0,nperm=1000)
save(simulation_result, file = paste0(data_path,"Paired_randomization_test_sup_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
time_rand_sup_nonmono=proc.time()-start.time
print(time_rand_sup_nonmono)


# simulation with non-monotonic effect function and YZratio=10
#This one took 2828 seconds on my laptop
start.time=proc.time()
sim1_raw_data_file = "../data_new/Paired_simData_interactive_10_Jun_06_points.RData"
TYPE='interactive'
YZratio=10
sample_sizes <-c(80, 100, 120)
dVec <- seq(from=.5, to=10, by=.25)
bw <- 1
simulation_result=run_simulation_randomization_test_sup(sim1_raw_data_file,bw,sample_sizes,polygon=0,nperm=1000)
save(simulation_result, file = paste0(data_path,"Paired_randomization_test_sup_", TYPE, "_" ,YZratio,'_',today_date,"_points.RData"))
time_rand_interactive_nonmono=proc.time()-start.time
print(time_rand_interactive_nonmono)
