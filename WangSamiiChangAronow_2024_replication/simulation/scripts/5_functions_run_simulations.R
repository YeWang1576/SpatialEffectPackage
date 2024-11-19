run_simulation=function(sim_data_path,cutoff,bw,sample_sizes,polygon){
  
  ############################################
  #########Load Data##########################
  ############################################
  load(sim_data_path)
  
  
  ############################################
  #######Initiate Empty Lists for Outputs#####
  ############################################
  output_AME_list <- output_AME_sm_list <-list()
  output_AME_est_list <- output_AME_est_sm_list <- list()
  output_all_VCEs_list <- output_all_VCEs_sm_list <- list()
  output_all_VCEs_pd_list <- output_all_VCEs_pd_sm_list<- list()
  output_all_SAH_list <- output_all_SAH_sm_list <-list()
  output_all_edof_list <- output_all_edof_sm_list <-list()
  output_all_edof_sd_list <- output_all_edof_sd_sm_list <-list()
  
  #print('remember to change it back later')
  for (ss in 1:length(sample_sizes)){
  #for (ss in 5){
  
    sample_size <- sample_sizes[ss]
    
    result_list <- data_list[[ss]]
    Ydata <- result_list[[1]] #outcome data
    Zdata <- result_list[[2]] #intervention node level data
    ras0 <- result_list[[3]]  #raster data
    ras_Z <- NULL
    
    
    #########################################################
    #######Calculate cutoff for the interacted case##########
    #########################################################
    
    if (polygon==1){
      ras_Z <- result_list[[4]]
      
      #for the polygon interventions with interactive case use the following 
      #note for the interative case this list has length 7, for the nonmonotonic case this list
      #has length 6
      if (length(result_list)==7){
        
        print('Polygon Intervention with interactive effects')
        cutoff <- 2*max(result_list[[5]])+ 2 * 6
        print(paste0('Cutoff is set to: ',cutoff))
        
      }
    }
    
    if (polygon==0){
      #for the point interventions with interactive case use the following 
      #note for the interative case this list has length 7, for the nonmonotonic case this list
      #has length 6
      if (length(result_list)==7){
        
        print('Polygon Intervention with interactive effects')
        cutoff <- 2*max(result_list[[5]])+ 2*6
        print(paste0('Cutoff is set to: ',cutoff))
        
      }
    }
    
    
    
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
    
    #####Create Distance Matrix########################
    x_coord_Z <- "x" #column names that has the x-coordinate information
    y_coord_Z <- "y" #column names that has the y-coordinate information
    

    if (polygon==0){
      x_coord <- Zdata[, x_coord_Z]
      y_coord <- Zdata[, y_coord_Z]
      dist <- DistanceCalculation(as.vector(x_coord), as.vector(y_coord), 1)[["Dist_mat"]]
    }
    
    if (polygon==1){
      print('Using the correct distance matrix')
      dist = drop_units(result_list[['dist_matrix_ZZ']])
      diag(dist)=0
      print(paste0('Maximum Distnace:',apply(dist,1,max)))
    }

    print('Calculate True AME')
    for(dIndex in 1:length(dVec)){
      #calculate true AME over d indices
      print(dIndex)
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
    output_all_edof_sd<- matrix(NA, length(dVec), nrow(ZMat)) #edof for AME
    output_all_edof_sd_sm<- matrix(NA, length(dVec), nrow(ZMat)) #edof for sAME
    
    
    #############################################################################
    ##################Set Simulation Parameters##################################
    #############################################################################
    
    
    # parameters for estimation
    dist.metric <- "Euclidean" #metric used to calculate distances
    numpts <- NULL #number of points to include
    only.unique <- 0 # 1 if each raster is included only once when calculating the circle averages
    k <- 1 #use a uniform kernel
    tr_prob <- 0.5 #treatment probability is 0.5 in the simulation
    
    #This is an important object: This list indicates which raster crosses the d-circle of an intervention node
    #Rasters that cross the dVec[d]-circle of the ith intervention node are in the list (i-1)*length(dVec)+ d 
    print('Calculate Circle Means')
    grid_list <- circleMean(Zdata=Zdata,ras0 = ras0, ras_Z = ras_Z, dVec = dVec, nz = ncol(ZMat)) 
    
    


    
    ###################################################
    ####Create Kernel Matrix for HAC Estimation
    ####################################################
    HAC_kernel_list=list()
    HAC_pd_kernel_list=list()
    HAC_kernel_sm_list=list()
    HAC_pd_kernel_sm_list=list()  
    for (d in 1:length(dVec)){
      
      c = cutoff + 2*dVec[d]
      HAC_kernel_list[[d]]=ConleySE_kernel(dist,nz,c,as.integer(k),0)[['Dist_kernel']]
      HAC_pd_kernel_list[[d]]=ConleySE_kernel(dist,nz,c,as.integer(k),1)[['Dist_kernel']]
      
      c_sm =c + 2*bw
      HAC_kernel_sm_list[[d]]=ConleySE_kernel(dist,nz,c_sm,as.integer(k),0)[['Dist_kernel']]
      HAC_pd_kernel_sm_list[[d]]=ConleySE_kernel(dist,nz,c_sm,as.integer(k),1)[['Dist_kernel']]
      
    }
    
    #simulation
#    print('remember to change it back later')
    
    for(zIndex in 1:nrow(ZMat)){
 #   for(zIndex in 1:10){
      

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
      output_all_edof_sd[,zIndex]=Conley_standard_error[5,] #edof
      print(paste0('edof:',      Conley_standard_error[4,]))
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
      output_all_edof_sd_sm[,zIndex]=Conley_standard_error_sm[5,] #edof
      
      
      ###################################################################################
      #####Print Progress Index##########################################################
      ###################################################################################
       if (zIndex %% 10 == 0){
         cat(zIndex, "\n")
         
       }
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
    
    output_all_edof_sd_list[[ss]] <- output_all_edof_sd
    output_all_edof_sd_sm_list[[ss]] <- output_all_edof_sd_sm   
    
  }
  
  #####################################################################
  #####################Output Simulation Files#########################
  #####################################################################
  
  return(  list(output_AME_list=output_AME_list , output_AME_sm_list=output_AME_sm_list ,
                output_AME_est_list=output_AME_est_list , output_AME_est_sm_list=output_AME_est_sm_list,
                output_all_VCEs_list=output_all_VCEs_list , output_all_VCEs_sm_list=output_all_VCEs_sm_list ,
                output_all_VCEs_pd_list =output_all_VCEs_pd_list , output_all_VCEs_pd_sm_list=output_all_VCEs_pd_sm_list,
                output_all_SAH_list=output_all_SAH_list ,output_all_SAH_sm_list=output_all_SAH_sm_list ,
                output_all_edof_list=output_all_edof_list , output_all_edof_sm_list=output_all_edof_sm_list,
                output_all_edof_sd_list=output_all_edof_sd_list,output_all_edof_sd_sm_list=output_all_edof_sd_sm_list ))
}
