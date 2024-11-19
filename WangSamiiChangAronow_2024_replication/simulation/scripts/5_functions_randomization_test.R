run_simulation_randomization_test=function(sim_data_path,bw,sample_sizes,polygon,nperm){
  
  ############################################
  #########Load Data##########################
  ############################################
  load(sim_data_path)
  
  
  ############################################
  #######Initiate Empty Lists for Outputs#####
  ############################################
  output_rt_list <- list()
  output_rt_sm_list <- list()

  
  for (ss in 1:length(sample_sizes)){

    sample_size <- sample_sizes[ss]
    
    result_list <- data_list[[ss]]
    Ydata <- result_list[[1]] #outcome data
    Zdata <- result_list[[2]] #intervention node level data
    ras0 <- result_list[[3]]  #raster data
    ras_Z <- NULL
    
    if (polygon==1){
      ras_Z <- result_list[[4]]
    }
    #########################################################
    #######Calculate cutoff for the interacted case##########
    #########################################################
    
    
    
    ny <- dim(Ydata)[1] #number of outcome points
    nz <- dim(Zdata)[1] #number of intervention nodes
    YobsAll <- Ydata[, grep("Y1_", names(Ydata))] #observed outcome
    ZMat <- t(Zdata[, grep("curZ", names(Zdata))]) #observed assignment
    
    cat("Sample size is: ", nz, "\n")
    
    
    #############################################################################
    ##################Create Empty Matrices to Store Data########################
    #############################################################################
    
    output_all_rt <- matrix(NA, length(dVec), nrow(ZMat)) #test results 
    output_all_rt_sm <- matrix(NA, length(dVec), nrow(ZMat)) #testing results

    
    #############################################################################
    ##################Set Simulation Parameters##################################
    #############################################################################
    
    
    # parameters for estimation
    tr_prob <- 0.5 #treatment probability is 0.5 in the simulation
    alpha=0.05 
    
    #This is an important object: This list indicates which raster crosses the d-circle of an intervention node
    #Rasters that cross the dVec[d]-circle of the ith intervention node are in the list (i-1)*length(dVec)+ d 
    print('Calculate Circle Means')
    grid_list <- circleMean(Zdata=Zdata,ras0 = ras0, ras_Z = ras_Z, dVec = dVec, nz = ncol(ZMat)) 
    
    
    
    ###########################################################################
    ##########Permutation Tests ###############################################
    ###########################################################################
    
    start.time=proc.time()
    for(zIndex in 1:nrow(ZMat)){

      
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
        
        #estimates
        est <- abs(Zup%*%Ybard/sum(Zup) - (1-Zup)%*%Ybard/sum(1-Zup))
        
        #quantiles of the randomization distribution
        alpha_quantile = randomization_test(Ybard,nperm,tr_prob,alpha)
        
        #whether the sharp null hypothesis is rejected
        output_all_rt[d,zIndex] = est > alpha_quantile

        
      }
      
      #Circle Averages for the smoothed AME
      for (d in 1:length(dVec)){
        
        #difference in distance
        diff = dVec[d]-dVec
        
        #triangular kernel
        triangular_kernel = sapply(diff, function(x) as.numeric((abs(x/bw) <= 1)*(1-abs(x/bw))))
        
        #renormalized_the_weights
        triangular_kernel = triangular_kernel/sum(triangular_kernel)
        
        Ybard2 = Ybards %*% triangular_kernel 
        
        Ybard2s[,d] = Ybard2
        
        
        #estimates
        est_sm <- abs(Zup%*%Ybard2/sum(Zup) - (1-Zup)%*%Ybard2/sum(1-Zup))
        
        #quantiles of the randomization distribution
        alpha_quantile_sm = randomization_test(Ybard2,nperm,tr_prob,alpha)
        
        #whether the sharp null hypothesis is rejected
        output_all_rt_sm[d,zIndex] = est_sm > alpha_quantile_sm       
        
      }

      
      ###################################################################################
      #####Print Progress Index##########################################################
      ###################################################################################
      if (zIndex %% 10 == 0){
        cat(zIndex, "\n")
        
        end.time=proc.time()
        # print(end.time-start.time)
        # 
        # start.time=proc.time()
      }
    }
    
    
    #Attch outputs to a list
    output_rt_list[[ss]] <- output_all_rt
    output_rt_sm_list[[ss]] <- output_all_rt_sm
    

    
  }
  
  #####################################################################
  #####################Output Simulation Files#########################
  #####################################################################
  
  return(  list(output_rt_list=output_rt_list , output_rt_sm_list=output_rt_sm_list))
}


randomization_test=function(Y,nperm,tr_prob,alpha){
  
  #the number of permutations
  n_node=length(Y)
  output=rep(NA,nperm)
  
  for (i in 1:nperm){
    
    #randomly draw treatments
    Zperm = rbinom(n_node,1,tr_prob)
    
    #Hajek estimators for this random draw
    est_perm=abs(mean(Y[which(Zperm==1)])-mean(Y[which(Zperm==0)]))
    
    #return estimates
    output[i]=est_perm
  }
  #return the 1-alpha quantile of the permutatoin estimates
  return(quantile(output,1-alpha) )
  
}

run_simulation_randomization_test_sup=function(sim_data_path,bw,sample_sizes,polygon,nperm){
  
  ############################################
  #########Load Data##########################
  ############################################
  load(sim_data_path)
  
  
  ############################################
  #######Initiate Empty Lists for Outputs#####
  ############################################
  output_rt_list <- list()
  output_rt_sm_list <- list()
  
  
  for (ss in 1:length(sample_sizes)){
    
    sample_size <- sample_sizes[ss]
    
    result_list <- data_list[[ss]]
    Ydata <- result_list[[1]] #outcome data
    Zdata <- result_list[[2]] #intervention node level data
    ras0 <- result_list[[3]]  #raster data
    ras_Z <- NULL
    
    if (polygon==1){
      ras_Z <- result_list[[4]]
    }
    #########################################################
    #######Calculate cutoff for the interacted case##########
    #########################################################
    
    
    
    ny <- dim(Ydata)[1] #number of outcome points
    nz <- dim(Zdata)[1] #number of intervention nodes
    YobsAll <- Ydata[, grep("Y1_", names(Ydata))] #observed outcome
    ZMat <- t(Zdata[, grep("curZ", names(Zdata))]) #observed assignment
    
    cat("Sample size is: ", nz, "\n")
    
    
    #############################################################################
    ##################Create Empty Matrices to Store Data########################
    #############################################################################
    
    output_all_rt_sup <- rep(NA,nrow(ZMat)) #test results 
    output_all_rt_sup_sm <- rep(NA,  nrow(ZMat)) #testing results
    
    
    #############################################################################
    ##################Set Simulation Parameters##################################
    #############################################################################
    
    
    # parameters for estimation
    tr_prob <- 0.5 #treatment probability is 0.5 in the simulation
    alpha=0.05 
    
    #This is an important object: This list indicates which raster crosses the d-circle of an intervention node
    #Rasters that cross the dVec[d]-circle of the ith intervention node are in the list (i-1)*length(dVec)+ d 
    print('Calculate Circle Means')
    grid_list <- circleMean(Zdata=Zdata,ras0 = ras0, ras_Z = ras_Z, dVec = dVec, nz = ncol(ZMat)) 
    
    
    
    ###########################################################################
    ##########Permutation Tests ###############################################
    ###########################################################################
    
    start.time=proc.time()
    for(zIndex in 1:nrow(ZMat)){
      
      
      nz <- dim(Zdata)[1] #number of intervention nodes
      treatment <- paste0("curZ", zIndex) #extract the treatment information for this simulation
      Zup <- Zdata[, treatment]  #extract the treatment information for this simulation
      
      
      Yobs <- Ydata[, paste0("Y1_", zIndex)] #extract observed outcomes at evaluation points
      Ybards <- matrix(NA, nz, length(dVec)) #observed d circle mean
      Ybard2s <- matrix(NA,nz, length(dVec)) #observed smoothed d circle mean
      
      #Hajek estimates at different d
      est = rep(NA,length(dVec))
      est_sm = rep(NA,length(dVec))
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
        
        #estimates
        est[d] <- abs(Zup%*%Ybard/sum(Zup) - (1-Zup)%*%Ybard/sum(1-Zup))
        
  
        
      }
      
      #Circle Averages for the smoothed AME
      for (d in 1:length(dVec)){
        
        #difference in distance
        diff = dVec[d]-dVec
        
        #triangular kernel
        triangular_kernel = sapply(diff, function(x) as.numeric((abs(x/bw) <= 1)*(1-abs(x/bw))))
        
        #renormalized_the_weights
        triangular_kernel = triangular_kernel/sum(triangular_kernel)
        
        Ybard2 = Ybards %*% triangular_kernel 
        
        Ybard2s[,d] = Ybard2
        
        
        #estimates
        est_sm[d] <- abs(Zup%*%Ybard2/sum(Zup) - (1-Zup)%*%Ybard2/sum(1-Zup))
        
 

      }
      #######################################################################
      #######################Randomization Tests#############################
      #######################################################################
      
      #quantiles of the randomization distribution
      alpha_quantile = randomization_test_sup(Ybards,nperm,tr_prob,alpha,dVec)
      
      
      #quantiles of the randomization distribution
      alpha_quantile_sm = randomization_test_sup(Ybard2s,nperm,tr_prob,alpha,dVec)
      
      
      #whether the sharp null hypothesis is rejected
      
      output_all_rt_sup[zIndex] = max(abs(est)) > alpha_quantile
      output_all_rt_sup_sm[zIndex] = max(abs(est_sm)) > alpha_quantile_sm       
      ###################################################################################
      #####Print Progress Index##########################################################
      ###################################################################################
      if (zIndex %% 10 == 0){
        cat(zIndex, "\n")
        
        end.time=proc.time()
        # print(end.time-start.time)
        # 
        # start.time=proc.time()
      }
    }
    
    
    #Attch outputs to a list
    output_rt_list[[ss]] <- output_all_rt_sup
    output_rt_sm_list[[ss]] <- output_all_rt_sup_sm
    
    
    
  }
  
  #####################################################################
  #####################Output Simulation Files#########################
  #####################################################################
  
  return(  list(output_rt_list=output_rt_list , output_rt_sm_list=output_rt_sm_list))
}


randomization_test=function(Y,nperm,tr_prob,alpha){
  
  #the number of permutations
  n_node=length(Y)
  output=rep(NA,nperm)
  
  for (i in 1:nperm){
    
    #randomly draw treatments
    Zperm = rbinom(n_node,1,tr_prob)
    
    #Hajek estimators for this random draw
    est_perm=abs(mean(Y[which(Zperm==1)])-mean(Y[which(Zperm==0)]))
    
    #return estimates
    output[i]=est_perm
  }
  #return the 1-alpha quantile of the permutatoin estimates
  return(quantile(output,1-alpha) )
  
}

randomization_test_sup=function(Y,nperm,tr_prob,alpha,dVec){
  
  #the number of permutations
  n_node=nrow(Y)
  output=rep(NA,nperm)
  
  #est
  est=rep(0,length(dVec))
  
  for (i in 1:nperm){
    
    #randomly draw treatments
    Zperm = rbinom(n_node,1,tr_prob)
    
    for (d in 1:length(dVec)){
      
      
      #extract information at distance value d
      Yd = Y[,d]
      
      #estimates
      est[d] <- abs(Zperm%*%Yd/sum(Zperm) - (1-Zperm)%*%Yd/sum(1-Zperm))
      
      
    }
    
    output[i] =  max(abs(est))
    
  }
  
  #return the 1-alpha quantile of the permutation estimates
  return(quantile(output,1-alpha) )
  
}
