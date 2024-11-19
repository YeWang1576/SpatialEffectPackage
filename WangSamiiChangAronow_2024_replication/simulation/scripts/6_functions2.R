################################################################################
##########This script contains functions for variance estimation################
################################################################################

calculate_cid_eachd=function(Zdata,cutoff,dtemp){
  
  #function calculate cid as defined in Section 4 of the paper for each distance d.
  
  cid = matrix(NA,nrow=nrow(Zdata),ncol=2)
  
  cid[,1]=Zdata$index
  
  #coordinate information
  x_coord=Zdata[,'x']
  y_coord=Zdata[,'y']
  
  
  c = cutoff + 2*dtemp
  
  cid[,2]=sapply(1:nrow(Zdata),function(i) return(calculate_cid_inner(Zdata[i,'x'],Zdata[i,'y'],c,x_coord,y_coord)))  
  
  
  
  return(cid)
  
}


weight_matrix=function(Zdata,cutoff,dtemp,pos_semidef=TRUE){
  
  
  #function calculate the weight matrix for HAC variance estimation
  
  num_node=nrow(Zdata)
  
  weights=matrix(0,num_node,num_node)
  
  c=cutoff+2*dtemp 
  
  for (i in 1:num_node){
    
    for (j in 1:num_node){
      
      distij = sqrt((Zdata[i,'x']-Zdata[j,'x'])^2 + (Zdata[i,'y']-Zdata[j,'y'])^2) 
      
      weights[i,j]=  (distij <= c)
      
    }
    
  }
  
  
  if (pos_semidef==TRUE){
    
    #eigen decomposition
    eigen_decomposition=eigen(weights)
    
    #eigen vectors
    eig_vector = eigen_decomposition$vectors
    #eigen values: first threshold out the negative values then create a diaognal matrix
    eig_value = eigen_decomposition$values
    eig_value[which(eig_value<0)] = 0
    eig_value=diag(eig_value)
    
    weights= eig_vector %*% eig_value %*% t(eig_vector)
    
  }
  
  return(weights)
}


calculate_var = function(lm_obj,cid,HAC_weight,HAC_weight_pos,prob,N){
  
  #Extract design matrix
  X=model.matrix(lm_obj)
  D=X[,'D']
  #inverse design matrix
  XXinv = solve(t(X)%*% X)
  
  #residuals
  res=lm_obj$residuals
  
  error_vcov = res %*% t(res)
  
  #sandwidch regular HAC
  meat_reg = error_vcov * HAC_weight
  HACvar = XXinv %*% ( t(X) %*% meat_reg %*% X) %*% XXinv
  
  #sandwidch posdef HAC
  meat_reg = error_vcov * HAC_weight_pos
  HACvar_pos = XXinv %*% ( t(X) %*% meat_reg %*% X) %*% XXinv
  
  #SAH var
  SAHvar=(sum(D* res^2 / prob * cid[,2]) + sum( (1-D)* (res^2 / (1-prob) * cid[,2])))/N^2
  
  return(list(SAH_var=SAHvar,HAC_var=HACvar[2,2],HAC_var_pos=HACvar_pos[2,2]))
  
}

calculate_var = function(lm_obj,cid,HAC_weight,HAC_weight_pos,prob,N){
  
  #Extract design matrix
  X=model.matrix(lm_obj)
  D=X[,'D']
  #inverse design matrix
  XXinv = solve(t(X)%*% X)
  
  #residuals
  res=lm_obj$residuals
  
  error_vcov = res %*% t(res)
  
  #sandwidch regular HAC
  #meat_reg = error_vcov * HAC_weight
  #HACvar = XXinv %*% ( t(X) %*% meat_reg %*% X) %*% XXinv
  
  #sandwidch posdef HAC
  #meat_reg = error_vcov * HAC_weight_pos
  #HACvar_pos = XXinv %*% ( t(X) %*% meat_reg %*% X) %*% XXinv
  
  #SAH var
  SAHvar=(sum(D* res^2 / prob * cid[,2]) + sum( (1-D)* (res^2 / (1-prob) * cid[,2])))/N^2
  
  return(SAHvar)
  
}




calculate_cid_alld=function(Zdata,cutoff,dVec){
  
  #function calculate cid as defined in Section 4 of the paper for all distances d in the paper
  
  cid = matrix(NA,nrow=nrow(Zdata),ncol=length(dVec)+1)
  
  cid[,1]=Zdata$index
  
  #coordinate information
  x_coord=Zdata[,'x']
  y_coord=Zdata[,'y']
  
  for (i in 1:length(dVec)){
    
    dtemp=dVec[i]
    
    
    c = cutoff + 2*dtemp
    
    cid[,i+1]=sapply(1:nrow(Zdata),function(i) return(calculate_cid_inner(Zdata[i,'x'],Zdata[i,'y'],c,x_coord,y_coord)))  
    
  }
  
  return(cid)
  
}


calculate_cid_inner=function(x,y,dist,x_coord_temp,y_coord_temp){
  
  #this function counts how many points are of distance "dist" away from x_coord_temp, y_coord_temp
  
  dist_to_xy= sqrt( (x_coord_temp  - x)^2 + (y_coord_temp-y)^2)
  
  cid= sum (dist_to_xy<=dist)
  
  return(cid)  
  
}