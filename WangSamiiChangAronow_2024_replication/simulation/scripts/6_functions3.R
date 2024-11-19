

#this function calculate outcomes for circles centered at d (smoothed or not smoothed)
construct_circle_outcome_inner = function(Yobs,indices,d,d_index,grid_list,dVec,smoothed=TRUE,bw=1){
  
  
  #Yobs: observed outcome data
  #indices: indices associated with intervention node i
  #d: the d value for which the circle average is constructed
  #d_index: the index of the d
  #grid_list: a list indicating which raster belongs to the d-cricel
  #dVec: a list of distances
  
  #average circle outcome by d
  outcome_by_d=sapply(indices, function(x) return(     mean(Yobs[grid_list[[x]]],na.rm=TRUE)))
  
  if (smoothed==TRUE){
    
    #difference in distance
    diff = d-dVec
    
    #triangular kernel
    triangular_kernel = sapply(diff, function(x) as.numeric((abs(x/bw) <= 1)*(1-(x/bw))))
    
    #renormalized_the_weights
    triangular_kernel = triangular_kernel/sum(triangular_kernel)
    
    #weighted averages
    outcome = sum(triangular_kernel * outcome_by_d)
    
  }else{
    
    #if not return pointwise estimate
    outcome = outcome_by_d[d_index]
    
  }
  
  
  
  return(outcome)
  
  
}



#this function calculate outcomes for circles centered at d (smoothed or not smoothed) for all units
construct_circle_outcome = function(i,Yobs,d,d_index,grid_list,dVec,smoothed=TRUE,bw=1){
  
  
  #Yobs: observed outcome data
  #indices: indices associated with intervention node i
  #d: the d value for which the circle average is constructed
  #d_index: the index of the d
  #grid_list: a list indicating which raster belongs to the d-cricel
  #dVec: a list of distances
  
  #grid for i
  grid_for_i = ((i-1)*length(dVec)+1):(i*length(dVec)) 
  
  #average circle outcome by d
  outcome_by_d=sapply(grid_for_i, function(x) return(     mean(Yobs[grid_list[[x]]],na.rm=TRUE)))
  
  if (smoothed==TRUE){
    
    #difference in distance
    diff = d-dVec
    
    #triangular kernel
    triangular_kernel = sapply(diff, function(x) as.numeric((abs(x/bw) <= 1)*(1-(x/bw))))
    
    #renormalized_the_weights
    triangular_kernel = triangular_kernel/sum(triangular_kernel)
    # print(triangular_kernel)
    #weighted averages
    outcome = sum(triangular_kernel * outcome_by_d)
    
  }else{
    
    #if not return pointwise estimate
    outcome = outcome_by_d[d_index]
    
  }
  
  
  
  return(outcome)
  
  
}


#Regression
Regression_AME= function(d,Y,X_mat,XX_mat_inv,W_meat,cutoff,dVec,dist,k,Zdata,prob,D,N,if_edof=1,kernel,kernel_pd){
  
  # #cid vector
  # 
  # if (ss > 2 & d < 12) {
  #   c <- cutoff + 2*dVec[d]
  # } else{
  #   c <- cutoff + 2*dVec[11]
  # }
  # 
  # 
  # 
  #kernel matrix
  kernel_matrix=kernel[[d]]
  kernel_pd_matrix=kernel_pd[[d]] 
  
  c = cutoff + 2*dVec[d]
  #cid vector

  nz=dim(Zdata)[1]
  cid = cbind(1:nz,apply(dist,1, function(x) sum(x<=c)))

  

  
  mu_d <- Y[, d]
  beta <- solve(t(X_mat) %*% (X_mat)) %*% (t(X_mat) %*% mu_d)
  res <- mu_d - X_mat %*% beta

  #Conley_result <- ConleySE(as.vector(res), as.matrix(W_meat), as.matrix(dist), as.matrix(XX_mat_inv),
  #                            as.double(c), as.integer(k),as.integer(0),as.integer(if_edof))

  Conley_result2 <- ConleySE2(as.vector(res), as.matrix(W_meat), as.matrix(kernel_matrix), as.matrix(XX_mat_inv),
                             as.integer(k),as.integer(0),as.integer(if_edof))

  
  #meat = Conley_result2[["VCE_meat"]]
  #meat2= Conley_result2[['VCE_meat2']]
  # print('Old Implemtation')
  # print(meat)
  # print('New Implemtation')
  # print(meat2)
  # 
   # VCE <- t(XX_mat_inv) %*% Conley_result[["VCE_meat"]] %*% XX_mat_inv
  VCE2 <- t(XX_mat_inv) %*% Conley_result2[["VCE_meat2"]] %*% XX_mat_inv
  
   # Conley_var <- VCE[2, 2]
  Conley_var2 <- VCE2[2, 2]  

  # 
  #for Young's correction
  # mu=Conley_result[['mu']]
  # v=Conley_result[['v']]

  mu2=Conley_result2[['mu']]
  v2=Conley_result2[['v']]
  
  # print('Old Implemtation')
  # print(mu)
  # print(v)
  # 
  # print('New Implemtation')
  # print(mu2)
  # print(v2)
  #edof
  # z_vec <- c(0, 1, rep(0, (ncol(X_mat)-2))) %*% XX_mat_inv %*% t(X_mat)
  # M_mat <- diag(1, nz, nz) - X_mat %*% XX_mat_inv %*% t(X_mat)
  # dist_kernel=Conley_result[['Dist_kernel']]
  # dof <- 1/(XX_mat_inv[2, 2])*sum(diag(M_mat %*% ((t(z_vec) %*% z_vec) * dist_kernel) %*% M_mat))
  # print('r implemtation')
  # print(dof)
  # print('cpp implementation')
  # print(Conley_result[['mu']])
  # print('distance')
  # print(dVec[d])
  # print("edof")
  # print(2*mu^2/v)
  # 
  # if (abs(dof-Conley_result[['mu']])>0.01){
  #   browser()
  # }
  #positive semidefinite variances
  #Conley_result_pd <- ConleySE(as.vector(res), as.matrix(W_meat), as.matrix(dist),as.matrix(XX_mat_inv),
  #                               as.double(c), as.integer(k),as.integer(1),as.integer(0))

  Conley_result2_pd <- ConleySE2(as.vector(res), as.matrix(W_meat), as.matrix(kernel_pd_matrix), as.matrix(XX_mat_inv),
                               as.integer(k),as.integer(0),as.integer(if_edof))
   # VCE_pd <- t(XX_mat_inv) %*% Conley_result_pd[["VCE_meat"]] %*% XX_mat_inv
  VCE_pd2 <- t(XX_mat_inv) %*% Conley_result2_pd[["VCE_meat2"]] %*% XX_mat_inv
  
  mu2_pd=Conley_result2_pd[['mu']]
  v2_pd=Conley_result2_pd[['v']]
   # Conley_var_pd <- VCE_pd[2, 2]
  Conley_var2_pd <- VCE_pd2[2,2]
  
  # print('Old Implemtation')
  # print(Conley_var_pd)
  # print('New Implemtation')
  # print(Conley_var2_pd)

  SAHvar=(sum(D* res^2 / prob * cid[,2]) + sum( (1-D)* (res^2 / (1-prob) * cid[,2])))/N^2
  
  return(c(Conley_var2/mu2,Conley_var2_pd/mu2_pd,SAHvar,edof=2*mu2^2/v2,edof_pd=mu2_pd^2/v2_pd))
}

Regression_AME_obs= function(d,Y,X_mat,XX_mat_inv,W_meat,cutoff,dVec,dist,k,Zdata,D,N,if_edof=1,kernel,kernel_pd){
  
  # #cid vector
  # 
  # if (ss > 2 & d < 12) {
  #   c <- cutoff + 2*dVec[d]
  # } else{
  #   c <- cutoff + 2*dVec[11]
  # }
  # 
  # 
  # 
  #kernel matrix
  kernel_matrix=kernel[[d]]
  kernel_pd_matrix=kernel_pd[[d]] 
  
  c = cutoff + 2*dVec[d]
  #cid vector
  
  nz=dim(Zdata)[1]
  cid = cbind(1:nz,apply(dist,1, function(x) sum(x<=c)))

  #calculate variance: preparation
  Z_tr = Zdata$treatment
  mu_d <- Y[, d]
  O=cbind(rep(1,nz),Zdata$cov1,Zdata$cov2,Zdata$cov3)
  
  #propensity score
  ps = Zdata$prob_treatment
  ps_w = ps*(1-ps)
  
  #coefficient
  OOw = t(O) %*% diag(ps_w) %*% O * (1/nz)
  Oy1w = t(O) %*% diag(ps_w) %*% (Ybards[,d] * Zdata$treatment * (1/ps)) * (1/nz)
  Oy0w = t(O) %*% diag(ps_w) %*% (Ybards[,d] * (1-Zdata$treatment) * (1/(1-ps))) * (1/nz)
  
  beta1=solve(OOw,Oy1w)
  beta0=solve(OOw,Oy0w)
  
  res = Z_tr * (mu_d - O %*% beta1) * 1/ps  * (sum(Z_tr))/nz + (1-Z_tr) * (mu_d - O %*% beta0) * (1/(1-ps)) * (nz-sum(Z_tr))/nz
  #Conley_result <- ConleySE(as.vector(res), as.matrix(W_meat), as.matrix(dist), as.matrix(XX_mat_inv),
  #                            as.double(c), as.integer(k),as.integer(0),as.integer(if_edof))
  
  Conley_result2 <- ConleySE2(as.vector(res), as.matrix(W_meat), as.matrix(kernel_matrix), as.matrix(XX_mat_inv),
                              as.integer(k),as.integer(0),as.integer(if_edof))
  
  
  #meat = Conley_result2[["VCE_meat"]]
  #meat2= Conley_result2[['VCE_meat2']]
  # print('Old Implemtation')
  # print(meat)
  # print('New Implemtation')
  # print(meat2)
  # 
  # VCE <- t(XX_mat_inv) %*% Conley_result[["VCE_meat"]] %*% XX_mat_inv
  VCE2 <- t(XX_mat_inv) %*% Conley_result2[["VCE_meat2"]] %*% XX_mat_inv
  
  # Conley_var <- VCE[2, 2]
  Conley_var2 <- VCE2[2, 2]  
  
  # 
  #for Young's correction
  # mu=Conley_result[['mu']]
  # v=Conley_result[['v']]
  
  mu2=Conley_result2[['mu']]
  v2=Conley_result2[['v']]
  
  # print('Old Implemtation')
  # print(mu)
  # print(v)
  # 
  # print('New Implemtation')
  # print(mu2)
  # print(v2)
  #edof
  # z_vec <- c(0, 1, rep(0, (ncol(X_mat)-2))) %*% XX_mat_inv %*% t(X_mat)
  # M_mat <- diag(1, nz, nz) - X_mat %*% XX_mat_inv %*% t(X_mat)
  # dist_kernel=Conley_result[['Dist_kernel']]
  # dof <- 1/(XX_mat_inv[2, 2])*sum(diag(M_mat %*% ((t(z_vec) %*% z_vec) * dist_kernel) %*% M_mat))
  # print('r implemtation')
  # print(dof)
  # print('cpp implementation')
  # print(Conley_result[['mu']])
  # print('distance')
  # print(dVec[d])
  # print("edof")
  # print(2*mu^2/v)
  # 
  # if (abs(dof-Conley_result[['mu']])>0.01){
  #   browser()
  # }
  #positive semidefinite variances
  #Conley_result_pd <- ConleySE(as.vector(res), as.matrix(W_meat), as.matrix(dist),as.matrix(XX_mat_inv),
  #                               as.double(c), as.integer(k),as.integer(1),as.integer(0))
  
  Conley_result2_pd <- ConleySE2(as.vector(res), as.matrix(W_meat), as.matrix(kernel_pd_matrix), as.matrix(XX_mat_inv),
                                 as.integer(k),as.integer(0),as.integer(if_edof))
  # VCE_pd <- t(XX_mat_inv) %*% Conley_result_pd[["VCE_meat"]] %*% XX_mat_inv
  VCE_pd2 <- t(XX_mat_inv) %*% Conley_result2_pd[["VCE_meat2"]] %*% XX_mat_inv
  
  mu2_pd=Conley_result2_pd[['mu']]
  v2_pd=Conley_result2_pd[['v']]
  # Conley_var_pd <- VCE_pd[2, 2]
  Conley_var2_pd <- VCE_pd2[2,2]
  
  # print('Old Implemtation')
  # print(Conley_var_pd)
  # print('New Implemtation')
  # print(Conley_var2_pd)
  
  SAHvar=(sum(D* res^2 / ps * cid[,2]) + sum( (1-D)* (res^2 / (1-ps) * cid[,2])))/nz^2
  
  return(c(Conley_var2/mu2,Conley_var2_pd/mu2_pd,SAHvar,edof=2*mu2^2/v2,edof_pd=mu2_pd^2/v2_pd))
}

#this function calculate outcomes for circles centered at d (smoothed or not smoothed)
construct_circle_outcome = function(Yobs,indices,d,d_index,grid_list,dVec,smoothed=TRUE,bw=1){
  

  #Yobs: observed outcome data
  #indices: indices associated with intervention node i
  #d: the d value for which the circle average is constructed
  #d_index: the index of the d
  #grid_list: a list indicating which raster belongs to the d-cricel
  #dVec: a list of distances
  
  #average circle outcome by d
  outcome_by_d=sapply(indices, function(x) return(     mean(Yobs[grid_list[[x]]],na.rm=TRUE)))

  if (smoothed==TRUE){
      
      #difference in distance
      diff = d-dVec
    
      #triangular kernel
      triangular_kernel = sapply(diff, function(x) as.numeric((abs(x/bw) <= 1)*0.75*(1-(x/bw)^2)))
      
      #renormalized_the_weights
      triangular_kernel = triangular_kernel/sum(triangular_kernel)
      
      #weighted averages
      outcome = sum(triangular_kernel * outcome_by_d)
      
  }else{
      
      #if not return pointwise estimate
      outcome = outcome_by_d[d_index]
    
  }

  

  return(outcome)
  
  
}

