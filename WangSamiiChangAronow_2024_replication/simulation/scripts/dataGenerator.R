# simulation results in Wang, Samii, Chang, and Aronow (2024)
rm(list=ls())

library(foreign)
library(sp)
library(fields)
library(raster)
library(sf)
library(Rcpp)
library(pbapply)



set.seed(2024)

##############################################################################
###################Set Root Directory#########################################
##############################################################################
root_dir = "C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/SpatialReplication/"
setwd(root_dir)
data_path <- "/simulation/data_new/"

##############################################################################
################Import Functions##############################################
##############################################################################
sourceCpp("./SpatialEffect/src/DistanceCalculation2.cpp")
sourceCpp("./SpatialEffect/src/calculate_tau.cpp")
sourceCpp("./SpatialEffect/src/effect_nonmo.cpp")
sourceCpp("./SpatialEffect/src/effect_interactive.cpp")
source("./simulation/scripts/functions_generatingdata.R")


##############################################################################
#################Today's Date#################################################
##############################################################################
today_date=format(Sys.Date(), format="%b_%d")



# generate data with non-monotonic effect function and YZratio=3
TYPE <- "nonmono"
sample_sizes <- c(80, 100, 120)
YZratio=3
data_list1=pblapply(sample_sizes,dataGenerator, type=TYPE,nsim=1000,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
save(data_list1, file = paste0("simData_", TYPE, "_" ,YZratio,'_',today_date,"_varyingsizes.RData"))


# generate data with interactive effect function and YZratio=3
TYPE <- "interactive"
sample_sizes <- c(80, 100, 120)
YZratio=3
data_list2=pblapply(sample_sizes,dataGenerator, type=TYPE,nsim=1000,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
save(data_list2, file = paste0("simData_", TYPE, "_" ,YZratio,'_',today_date,"_varyingsizes.RData"))

# generate data with non-monotonic effect function and YZratio=10
TYPE <- "nonmono"
sample_sizes <- c(40,60,80, 100, 120)
YZratio=10
data_list3=pblapply(sample_sizes,dataGenerator, type=TYPE,nsim=1000,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
save(data_list3, file = paste0("simData_", TYPE, "_" ,YZratio,'_',today_date,"_varyingsizes.RData"))

# generate data with interactive effect function and YZratio=10
TYPE <- "interactive"
sample_sizes <- c(40,60,80, 100, 120)
YZratio=10
data_list4=pblapply(sample_sizes,dataGenerator, type=TYPE,nsim=1000,tr_prob=0.5,YZratio=YZratio, unit_size=1,  cr=FALSE, polygon=0, trim_ras=1, effect=3, Y_sd=1, Zjitter=0.1)
save(data_list4, file = paste0("simData_", TYPE, "_" ,YZratio,'_',today_date,"_varyingsizes.RData"))
