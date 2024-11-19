# simulation results in Wang, Samii, Chang, and Aronow (2024)
rm(list=ls())


root_dir = "C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/SpatialReplication/"
setwd(root_dir)
#setwd("~/Dropbox/CurrentProjects/SpatialReplication/simulation/graphs")

#ri package no longer available from repository
install.packages('./simulation/ri_0.9.tar.gz',type='source')

library(Rcpp)
library(Matrix)
library(RcppArmadillo)
library(plotrix)
library(geosphere)
library(foreign)
library(sp)
library(fields)
library(ri)
library(sf)
library(raster)

#############################################################################
#################Import Source Functions#####################################
#############################################################################
sourceCpp("./SpatialEffect/src/DistanceCalculation.cpp")
sourceCpp("./SpatialEffect/src/ConleySE_kernel_matrix.cpp")
sourceCpp("./SpatialEffect/src/Conley2.cpp")
source("./simulation/scripts/functions.R")
source("./simulation/scripts/functions2.R")
source("./simulation/scripts/functions3.R")



set.seed(2024)

# simulation results
type <- "nonmono"
data_path <- "./simulation/data/"
data_save_path <- "./simulation/data_new/"
load(paste0(data_save_path, "simData_nonmono_May_09_varyingsizes.RData"))
#load(paste0(data_path, "simData_nonmono_varyingsizes.RData"))

#sample_sizes <- c(20, 40, 60, 80, 100, 120)
sample_sizes <- c(80, 100, 120)

today_date=format(Sys.Date(), format="%b_%d")

len <- 10
dVec <- seq(from=.5, to=len, by=.25)

# for smoothed estimator
bw <- 1


#lists for store results
output_AME_list <- output_AME_sm_list <-list()
output_AME_est_list <- output_AME_est_sm_list <- list()
output_all_VCEs_list <- output_all_VCEs_sm_list <- list()
output_all_VCEs_pd_list <- output_all_VCEs_pd_sm_list<- list()
output_all_SAH_list <- output_all_SAH_sm_list <-list()
output_all_edof_list <- output_all_edof_sm_list <-list()


for (ss in 1){
  
  sample_size <- sample_sizes[ss]
  
  result_list <- data_list[[ss]]
  Ydata <- result_list[[1]] #outcome data
  Zdata <- result_list[[2]] #intervention node level data
  ras0 <- result_list[[3]]  #raster data
  ras_Z <- NULL 
  ny <- dim(Ydata)[1] #number of outcome points
  nz <- dim(Zdata)[1] #number of intervention nodes
  YobsAll <- Ydata[, grep("Y1_", names(Ydata))] #observed outcome
  ZMat <- t(Zdata[, grep("curZ", names(Zdata))]) #observed assignment

  cat("Sample size is: ", nz, "\n")
  
  #to hold the true parameter
  AMEmat <- data.frame(d = dVec, tauda2 = NA, nd = NA)
  AMEmat_sm <- data.frame(d = dVec, tauda2 = NA, nd = NA)
  
  #Average Effect Function for each z
  AMEmat_perz = matrix(NA,nrow=nz,ncol=length(dVec))
  
  for(dIndex in 1:length(dVec)){
    #calculate true AME over d indices
    dUp <- dVec[dIndex]
    AMEmat$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz, numpts = NULL, only.unique = 0)
    #AMEmat_perz[,dIndex] <- circleMean_oracle_perz(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz, numpts = NULL, only.unique = 0)
    
  }
  
  for (d in 1:length(dVec)){
    
    #difference in distance
    diff = dVec[d]-dVec
    
    #triangular kernel
    triangular_kernel = sapply(diff, function(x) as.numeric((abs(x/bw) <= 1)*(1-(abs(x)/bw))))
    
    #renormalized_the_weights
    triangular_kernel = triangular_kernel/sum(triangular_kernel)
    
    AMEmat_sm$tauda2[d]  = AMEmat$tauda2 %*% triangular_kernel 
    
    
  }
  
  
  output_AME_list[[ss]] <- AMEmat[, 1:2]
  output_AME_sm_list[[ss]] <- AMEmat_sm[, 1:2]
  
  #############################################################################
  ##################Create Empty Matrices to Store Data########################
  #############################################################################
  
  output_AMEhat_mean <- matrix(NA, length(dVec), nrow(ZMat))
  output_AMEhat_mean_sm <- matrix(NA, length(dVec), nrow(ZMat))
  output_all_VCEs <- matrix(NA, length(dVec), nrow(ZMat)) #HAC
  output_all_VCEs_sm <- matrix(NA, length(dVec), nrow(ZMat)) #HAC for sAME
  output_all_VCEs_pd <- matrix(NA, length(dVec), nrow(ZMat)) #HAC_PD
  output_all_VCEs_pd_sm <- matrix(NA, length(dVec), nrow(ZMat)) #HAC_PD for sAME
  output_all_SAH<- matrix(NA, length(dVec), nrow(ZMat)) #SAH for AME
  output_all_SAH_sm<- matrix(NA, length(dVec), nrow(ZMat)) #SAH for sAME
  output_all_edof<- matrix(NA, length(dVec), nrow(ZMat)) #edof for AME
  output_all_edof_sm<- matrix(NA, length(dVec), nrow(ZMat)) #edof for sAME
  
  
  #############################################################################
  ##################Set Simulation Parameters##################################
  #############################################################################
  
  
  # parameters for estimation
  x_coord_Z <- "x" #column names that has the x-coordinate information
  y_coord_Z <- "y" #column names that has the y-coordinate information
  dist.metric <- "Euclidean" #metric used to calculate distances
  numpts <- NULL #number of points to include
  only.unique <- 0 # 1 if each raster is included only once when calculating the circle averages
  cutoff <- 12 #this is the 2*bar{d} in the paper
  k <- 1 #use a uniform kernel
  tr_prob <- 0.5 #treatment probability is 0.5 in the simulation
  
  #This is an important object: This list indicates which raster crosses the d-circle of an intervention node
  #Rasters that cross the dVec[d]-circle of the ith intervention node are in the list (i-1)*length(dVec)+ d 
  grid_list <- circleMean(ras0 = ras0, ras_Z = ras_Z, dVec = dVec, nz = ncol(ZMat)) 
  
  
  #####Create Distance Matrix########################
  x_coord <- Zdata[, x_coord_Z]
  y_coord <- Zdata[, y_coord_Z]
  dist <- DistanceCalculation(as.vector(x_coord), as.vector(y_coord), 1)[["Dist_mat"]]
  
  
  ###################################################
  ####Create Kernel Matrix for HAC Estimation
  ####################################################
  HAC_kernel_list=list()
  HAC_pd_kernel_list=list()
  HAC_kernel_sm_list=list()
  HAC_pd_kernel_sm_list=list()  
    for (d in 1:length(dVec)){
     
     print(d)
     c = cutoff + 2*dVec[d]
     HAC_kernel_list[[d]]=ConleySE_kernel(dist,nz,c,as.integer(k),0)[['Dist_kernel']]
     HAC_pd_kernel_list[[d]]=ConleySE_kernel(dist,nz,c,as.integer(k),1)[['Dist_kernel']]
     
     c_sm =c + 2*bw
     HAC_kernel_sm_list[[d]]=ConleySE_kernel(dist,nz,c_sm,as.integer(k),0)[['Dist_kernel']]
     HAC_pd_kernel_sm_list[[d]]=ConleySE_kernel(dist,nz,c_sm,as.integer(k),1)[['Dist_kernel']]
     
  }
  
  for(dIndex in 1:length(dVec)){
    #calculate true AME over d indices
    dUp <- dVec[dIndex]
    AMEmat$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz, numpts = NULL, only.unique = 0)
    #AMEmat_perz[,dIndex] <- circleMean_oracle_perz(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz, numpts = NULL, only.unique = 0)
    
  }

  #simulation
  for(zIndex in 1:nrow(ZMat)){
    #for(zIndex in 1:100){
    
    print(zIndex)
    
    nz <- dim(Zdata)[1] #number of intervention nodes
    treatment <- paste0("curZ", zIndex) #extract the treatment information for this simulation
    Zup <- Zdata[, treatment]  #extract the treatment information for this simulation
    
    
    Yobs <- Ydata[, paste0("Y1_", zIndex)] #extract observed outcomes at evaluation points
    Ybards <- matrix(NA, nz, length(dVec)) #observed d circle mean
    Ybard2s <- matrix(NA,nz, length(dVec)) #observed smoothed d circle mean
    
    ##########################################################################
    ##################Create and Estimate Circle Averages ####################
    ##########################################################################

    #Circle Averages 
    for (d in 1:length(dVec)){
      
      Ybard <- rep(NA, nz)
      
      #calculate outcome
      for (i in 1:nz){
        #grid for calculate unit i's circle outcome
        grid_for_i = ((i-1)*length(dVec)+1):(i*length(dVec))
        
        #calculate circle average
        grid_index <- grid_list[[(i-1)*length(dVec)+d]]
        grids <- Yobs[grid_index]
        Ybard[i] <- mean(grids, na.rm = TRUE)
        
      }
      
      
      Ybards[, d] <- Ybard
      
      output_AMEhat_mean[d, zIndex] <- Zup%*%Ybard/sum(Zup) - (1-Zup)%*%Ybard/sum(1-Zup)
      
      
    }

    #Circle Averages for the smoothed AME

    for (d in 1:length(dVec)){
      
      #difference in distance
      diff = dVec[d]-dVec
      
      #triangular kernel
      triangular_kernel = sapply(diff, function(x) as.numeric((abs(x/bw) <= 1)*(1-abs(x/bw))))
      
      #renormalized_the_weights
      triangular_kernel = triangular_kernel/sum(triangular_kernel)
      
      Ybard2s[,d] = Ybards %*% triangular_kernel 
      
      output_AMEhat_mean_sm[d, zIndex] <- Zup%*%Ybard2s[,d]/sum(Zup) - (1-Zup)%*%Ybard2s[,d]/sum(1-Zup)
      
    }
    ##########################################################################
    ##################Some Preliminary Data Creation #########################
    ##########################################################################
    
    X_mat <- cbind(rep(1, nz), Zup)
    W_meat <- X_mat
    XX_mat_inv <- solve(t(X_mat) %*% X_mat)
    
    ##########################################################################
    ########################Regression Estimator##############################
    ########################Conley Standard Error#############################
    ##########################################################################

    
    Conley_standard_error=sapply(1:length(dVec),Regression_AME,Y=Ybards,X_mat=X_mat,XX_mat_inv=XX_mat_inv,W_meat=W_meat,cutoff=cutoff,dVec=dVec,dist=dist,k=k,Zdata=Zdata,prob=tr_prob,D=ZMat[zIndex,],N=nz,if_edof=1,kernel=HAC_kernel_list,kernel_pd=HAC_pd_kernel_list)
    
    output_all_VCEs[,zIndex] = Conley_standard_error[1,] #HAC
    output_all_VCEs_pd[,zIndex] = Conley_standard_error[2,] #HAC_PD
    output_all_SAH[,zIndex] = Conley_standard_error[3,] #SAH
    output_all_edof[,zIndex] = Conley_standard_error[4,] #edof
    

    ##########################################################################
    ########################Regression Estimator Smoothed#####################
    ########################Conley Standard Error#############################
    ##########################################################################
    

    c=cutoff+2*bw
    Conley_standard_error_sm=sapply(1:length(dVec),Regression_AME,Y=Ybard2s,X_mat=X_mat,XX_mat_inv=XX_mat_inv,W_meat=W_meat,cutoff=c,dVec=dVec,dist=dist,k=k,Zdata=Zdata,prob=tr_prob,D=ZMat[zIndex,],N=nz,if_edof=1,kernel=HAC_kernel_sm_list,kernel_pd=HAC_pd_kernel_sm_list)

    output_all_VCEs_sm[,zIndex] = Conley_standard_error_sm[1,] #HAC
    output_all_VCEs_pd_sm[,zIndex] = Conley_standard_error_sm[2,] #HAC_PD
    output_all_SAH_sm[,zIndex] = Conley_standard_error_sm[3,] #SAH 
    output_all_edof_sm[,zIndex] = Conley_standard_error_sm[4,] #edof
    

    ###################################################################################
    #####Print Progress Index##########################################################
    ###################################################################################
    # if (zIndex %% 100 == 0){
    #   cat(zIndex, "\n")
    #   
    # }
  }
  
  
  #Attch outputs to a list
  output_AME_est_list[[ss]] <- output_AMEhat_mean
  output_AME_est_sm_list[[ss]] <- output_AMEhat_mean_sm
  
  output_all_VCEs_list[[ss]] <- output_all_VCEs
  output_all_VCEs_sm_list[[ss]] <- output_all_VCEs_sm
  
  output_all_VCEs_pd_list[[ss]] <- output_all_VCEs_pd
  output_all_VCEs_pd_sm_list[[ss]] <- output_all_VCEs_pd_sm
  
  output_all_SAH_list[[ss]] <- output_all_SAH
  output_all_SAH_sm_list[[ss]] <- output_all_SAH_sm
  
  output_all_edof_list[[ss]] <- output_all_edof
  output_all_edof_sm_list[[ss]] <- output_all_edof_sm
  
  
  
}

#check homophly condition



save.image(file=paste0(data_save_path,"all_" ,today_date,".RData")) 
