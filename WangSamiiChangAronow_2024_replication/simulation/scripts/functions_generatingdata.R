##############################################################################
#The following function generates data we use for simulation################## 
##############################################################################


dataGenerator <- function(nsim,scale, tr_prob, YZratio, unit_size,  cr = FALSE, polygon = 0, trim_ras = 0, effect = 3, Y_sd = 0.1, Zjitter = 1, type = "non-monotonic"){
  
  
  #nsim: number of simulations
  #tr_prob: treatment probability
  #YZratio: number of unit distance (per side) for one intervention node
  #
  
  s <- scale
  s2 <- scale^2
  p <- round(s2*tr_prob)
  
  #generate a raster
  ras0 <- raster(	nrows=s, ncols=s, 
                  xmn=0, xmx=s*unit_size, 
                  ymn=0, ymx=s*unit_size, 
                  vals=rnorm(s2, mean=0, sd=Y_sd))
  
  # Dataset with Y0 values and then centroids of the raster cells
  Y0 <- getValues(ras0)
  Yxy <- coordinates(ras0)
  Yindex <- 1:nrow(Yxy)
  Ydata <- data.frame(Yindex,cbind(Y0,Yxy))
  ny=length(Y0)
  
  #############################################################################
  ###########Generate Control Outcomes with Kriging############################
  #############################################################################
  FExy <- coordinates(aggregate(ras0, fact=(s/4), fun=mean)) #intervention nodes
  FExy[, 1] <- jitter(FExy[, 1], amount=1)
  FExy[, 2] <- jitter(FExy[, 2], amount=1)
  index <- 1:nrow(FExy)
  FE0 <- abs(rnorm(nrow(FExy), mean=0, sd=1))
  FEdata <- as.data.frame(cbind(index, FExy, FE0))
  krigfit <- Krig(FEdata[, c(2:3)], FEdata[, 4], m = 2, lambda = 1)
  alpha <- predict(krigfit, Yxy)
  
  #############################################################################
  ###########For polygon interventions#########################################
  #############################################################################
  
  if (polygon == 1){
    

    Y0_matrix <- matrix(Y0, s, s)
    crs(ras0) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
    if ((trim_ras == 1) & (s >= 40)){
      sub_ras0 <- raster(nrows=(s-20),
                         ncols=(s-20),
                         xmn=11*unit_size, xmx=(s-10)*unit_size,
                         ymn=11*unit_size, ymx=(s-10)*unit_size,
                         vals=c(Y0_matrix[c(11:(s-10)), c(11:(s-10))]))
    }else{
      sub_ras0 <- ras0
    }
    crs(sub_ras0) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
    ras_polygon <- as(sub_ras0, 'SpatialPolygonsDataFrame')
    sf_raster <- st_as_sf(ras_polygon)
    polygon_points <- st_sample(sf_raster, round(s2/(YZratio^2/10), 0)) %>%
      st_sf() %>%
      st_transform(crs = crs(sub_ras0))
    
    all_polygons <- polygon_points %>%  # consider the master points
      st_geometry() %>% # ... as geometry only (= throw away the data items)
      st_union() %>% # unite them ...
      st_voronoi() %>% # ... and perform the voronoi tessellation
      st_collection_extract(type = "POLYGON") %>% # select the polygons
      st_sf(crs = crs(sub_ras0)) 
    
    all_polygons <- all_polygons[st_is_valid(all_polygons), ]
    all_polygons <- st_intersection(all_polygons, st_union(sf_raster))
    
    sample_points <- sample(round(s2/(YZratio^2/10), 0), (s2/YZratio^2))
    ras_Z <- all_polygons[sample_points, ]
    Zxy <- st_coordinates(st_centroid(ras_Z))
    index <- 1:nrow(Zxy)
    Zdata <- as.data.frame(cbind(index, Zxy))
    z2 <- dim(Zdata)[1]
    
    Ydata_coords <- st_as_sf(Ydata, crs = st_crs(ras_Z), coords = c("x", "y"))
    YZ_dist <- st_distance(Ydata_coords, ras_Z, which = "Euclidean")
    for (i in 1:dim(Zdata)[1]){
      Ydata[, i+4] <- YZ_dist[, i]
      names(Ydata)[i+4] <- paste0("dZ", i)
    }
  }
  #############################################################################
  ###########For point interventions###########################################
  ############################################################################# 
  
  if (polygon==0){
   

    #generate intervention notes
    ras_Z <- aggregate(ras0, fact=YZratio, fun=mean)
    Zxy <- coordinates(ras_Z)
    Zxy[, 1] <- jitter(Zxy[, 1], amount=Zjitter)
    Zxy[, 2] <- jitter(Zxy[, 2], amount=Zjitter)
    index <- 1:nrow(Zxy)
    Zdata <- as.data.frame(cbind(index, Zxy))
    nz <- dim(Zdata)[1]
    

    #Calculate distances from intervention nodes to evaluation points and intervention nodes
    node_x = Zdata[,2]
    node_y = Zdata[,3]
    outcome_x= Ydata[,3]
    outcome_y= Ydata[,4]
    
    #distance matrix from evaluation points to intervention nodes
    dist_matrix_YZ=DistanceCalculation2(outcome_x,outcome_y,node_x,node_y,as.integer(1))[['Dist_mat']]
    Ydata=cbind(Ydata,dist_matrix_YZ)
    names(Ydata)[5:ncol(Ydata)] <- sapply(1:nz,function(i) return(paste0("dZ",i)))
    
    #distance matrix from intervention nodes to intervention nodes
    dist_matrix_ZZ=DistanceCalculation2(node_x,node_y,node_x,node_y,as.integer(1))[['Dist_mat']]
    Zdata=cbind(Zdata,dist_matrix_ZZ)
    names(Zdata)[4:ncol(Zdata)] <- sapply(1:nz,function(i) return(paste0("dZ",i)))
    

  }
  
  #############################################################################
  ###########Generate Simualtion Points########################################
  ############################################################################# 
  
  curZ <- matrix(NA, nz, nsim)
  for (c in 1:nsim){
    if (cr == FALSE){
      #redraw if there is no treated or control nodes
      tempZ <- rbinom(nz, 1, tr_prob)
      while (sum(tempZ) == 0 | sum(tempZ) == nz){
        tempZ <- rbinom(nz, 1, tr_prob)
      }
    }else if (cr == TRUE){
      
      tempZ <- rep(0, nz) 
      tempZ[sample(nz, nz*mean(tr_prob))] <- 1
    }
    curZ[, c] <- tempZ
  }
  colnames(curZ) <- paste0("curZ",1:ncol(curZ))
  


  ############################################################################
  ###################Calculate Effects########################################
  ############################################################################
  cat("Calculate effects", "\n")
  
  
  ############################################################################
  ###################Nonmonotonic Effect######################################
  ############################################################################  
  if (type=="nonmono"){
    
    ################
    ##parameters####
    ################      
    sh=1
    sc=0.5
    a=1
    
    ########################
    ##calculate outcomes####
    ########################       
    Y1s = effect_nonmo1(Y0,dist_matrix_YZ,curZ,alpha, effect,sh,  sc,  1)[['Sim_data']]
    colnames(Y1s)=paste0("Y1_",1:ncol(curZ))
    

  }
  
  
  ############################################################################
  ###################Interactive Effect######################################
  ############################################################################
  

  if (type=="interactive"){
    
    #create a vector of nearest neighbor
    diag(dist_matrix_ZZ)=Inf
    Zneighbor = apply(dist_matrix_ZZ,1,which.min)
    Zneighbot_distance = apply(dist_matrix_ZZ,1,min)
    ################
    ##parameters####
    ################      
    sh=1
    sc=0.5
    a=1
    
    ########################
    ##calculate outcomes####
    ########################       
    Y1s = effect_interactive(Y0,dist_matrix_YZ,curZ,Zneighbor,alpha, effect,sh,  sc, a)[['Sim_data']]
    colnames(Y1s)=paste0("Y1_",1:ncol(curZ))
    

  }
  
  # for(k in 1:ncol(curZ)){
  #   print(k)
  #  if (type == "monotonic"){
  #     Y1s[, k] <- Ydata$Y0 + effect*alpha*apply(Ydata[,grep("dZ", names(Ydata))],
  #                                               1,
  #                                               function(x){sum(monoEffect(as.numeric(x), effect)*curZ[,k])})
  #   }else if (type == "null"){
  #     Y1s[, k] <- Ydata$Y0
  #   }
  # 
  # 
  # }
  # 
  ############################################################################
  ###################This block below generate (5) in the paper################
  ############################################################################

  tau=calculate_tau(Y1s,curZ)[['tau']]
  colnames(tau)=sapply(1:nz,function(i) return(paste("tau",i,sep="")))
  Ydata = cbind(Ydata,tau)
  
  ############################################################################
  ###################Bind Data################################################
  ############################################################################  
  
  Ydata <- cbind(Ydata, Y1s)
  Zdata <- cbind(Zdata, curZ)
  
  if (type=='interactive'){
    result_list <- list(Ydata, Zdata, ras0, ras_Z,Zneighbot_distance)
    
  }else{
    result_list <- list(Ydata, Zdata, ras0, ras_Z)
    
  }
  return(result_list)
}



