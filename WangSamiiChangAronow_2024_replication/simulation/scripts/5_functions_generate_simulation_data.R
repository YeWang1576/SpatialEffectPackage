##############################################################################
#The following function generates data we use for simulation################## 
##############################################################################


dataGenerator <- function(nsim,scale, tr_prob, YZratio, unit_size,  cr = FALSE, pair=FALSE,polygon = 0, trim_ras = 0, effect = 3, Y_sd = 0.1, Zjitter = 1, type = "non-monotonic"){
  
  
  ##########################################################################
  s <- scale #length of one side
  s2 <- scale^2 #length of the square
  p <- round(s2*tr_prob)
  
  ##########################################################################
  ############Generate a raster#############################################
  ##########################################################################
  
  #generate a raster
  ras0 <- raster(	nrows=s, ncols=s, 
                  xmn=0, xmx=s*unit_size, 
                  ymn=0, ymx=s*unit_size, 
                  vals=rnorm(s2, mean=0, sd=Y_sd))
  
  # Dataset with Y0 values and then centroids of the raster cells
  Y0 <- getValues(ras0) #control outcomes
  Yxy <- coordinates(ras0) 
  Yindex <- 1:nrow(Yxy)
  Ydata <- data.frame(Yindex,cbind(Y0,Yxy))
  ny=length(Y0)
  
  #############################################################################
  ###########Generate alpha_x with Kriging#####################################
  #############################################################################
  FExy <- coordinates(aggregate(ras0, fact=(s/4), fun=mean)) #coarsened rasters
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
    
    #control outcome matrix
    Y0_matrix <- matrix(Y0, s, s)
    
    #add a coordinate reference system
    crs(ras0) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
    
    if ((trim_ras == 1) & (s >= 40)){
    
      #trim the rasters on the sides: this way we make sure that the intervention nodes and their sphere of influence are all within the rasters  
      sub_ras0 <- raster(nrows=(s-20),
                         ncols=(s-20),
                         xmn=11*unit_size, xmx=(s-10)*unit_size,
                         ymn=11*unit_size, ymx=(s-10)*unit_size,
                         vals=c(Y0_matrix[c(11:(s-10)), c(11:(s-10))]))
      
      
    }else{
      
      sub_ras0 <- ras0
    
    }
    
    #inherit the coordinate reference system
    crs(sub_ras0) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
    
    ras_polygon <- as(sub_ras0, 'SpatialPolygonsDataFrame')
    sf_raster <- st_as_sf(ras_polygon)
    
    ############################################################################
    ##################Subsample Points for Creating Polygon#####################
    ############################################################################    
    polygon_points <- st_sample(sf_raster, round(s2/10, 0)) %>%
      st_sf() %>%
      st_transform(crs = crs(sub_ras0))
    
    ############################################################################
    ##################Creating Polygons#########################################
    ############################################################################    
    all_polygons <- polygon_points %>%  # consider the master points
      st_geometry() %>% # ... as geometry only (= throw away the data items)
      st_union() %>% # unite them ...
      st_voronoi() %>% # ... and perform the voronoi tessellation
      st_collection_extract(type = "POLYGON") %>% # select the polygons
      st_sf(crs = crs(sub_ras0)) 
    
    all_polygons <- all_polygons[st_is_valid(all_polygons), ]
    all_polygons <- st_intersection(all_polygons, st_union(sf_raster))
    
    ############################################################################
    ##################Sample Intervention Polygons##############################
    ############################################################################  

    if (pair){
      
      #randomly sample some polygons (half the size of the number of intervention nodes)
      
      sample_points <- sample(round(s2/10, 0), (s2/YZratio^2)/2)
      
      index_polygons=1:nrow(all_polygons)
      index_to_include=c()
      
      for (i in sample_points){
        


        #distance between polytopes
        dist_temp=sapply( 1:nrow(all_polygons), function(j) return(st_distance(all_polygons[j,],all_polygons[i,])))
    
        #set distance to itself as infinity    
        dist_temp[i]=Inf
        
        #find the adjacent neighbors and randomly select one
        neigbor_temp=(1:nrow(all_polygons))[which(dist_temp<1e-8)]
        
        #randomly sample a pair
        i_pair_temp=sample(neigbor_temp,1)
        
        #append the pair of polygons
        index_to_include=c(index_to_include,i,i_pair_temp)
        
      }  
      
      #create rasters
      ras_Z = all_polygons[index_to_include,]
      
    }else{
      
      #randomly sample some polygons
      sample_points <- sample(round(s2/(YZratio^2/10), 0), (s2/YZratio^2))
      ras_Z <- all_polygons[sample_points,]
      
    }

    
    Zxy <- st_coordinates(st_centroid(ras_Z))
    index <- 1:nrow(Zxy)
    Zdata <- as.data.frame(cbind(index, Zxy))
    names(Zdata)[2:3] <- c("x", "y")
    
    
    nz <- dim(Zdata)[1]
    Ydata_coords <- st_as_sf(Ydata, crs = st_crs(ras_Z), coords = c("x", "y"))
  
    #distance matrix from evaluation points to intervention nodes
    dist_matrix_YZ <- st_distance(Ydata_coords, ras_Z, which = "Euclidean")
    Ydata=cbind(Ydata,dist_matrix_YZ)
    names(Ydata)[5:ncol(Ydata)] <- sapply(1:nz,function(i) return(paste0("dZ",i)))

    #distance matrix from intervention nodes to intervention nodes
    dist_matrix_ZZ <- st_distance(ras_Z, ras_Z, which = "Euclidean")
    Zdata=cbind(Zdata,dist_matrix_ZZ)
    names(Zdata)[4:ncol(Zdata)] <- sapply(1:nz,function(i) return(paste0("dZ",i)))
    
    #diamaeter of polygons
    
    diameter = rep(0,nrow(ras_Z))
    for (i in 1:nrow(ras_Z)){
     diameter[i]=ras_Z[i,] %>% 
      st_cast('MULTIPOINT') %>% 
      st_cast('POINT') %>% 
      st_distance %>% 
      max()
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
    
    if (pair){
      Zxy1 <- coordinates(ras_Z)
      Zxy1[, 1] <- jitter(Zxy1[, 1], amount=Zjitter)
      Zxy1[, 2] <- jitter(Zxy1[, 2], amount=Zjitter)
      sample_ind <- sort(sample(index, length(index)/2))
      Zxy <- t(matrix(c(t(cbind(Zxy[sample_ind, ], Zxy1[sample_ind, ]))), nrow = 2))
    }

    
    #create node data
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
    
    names(Zdata)[2:3] <- c("x", "y")
    

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
    Y1s = effect_nonmo1(Y0,dist_matrix_YZ,curZ,alpha, effect,sh,  sc,  a)[['Sim_data']]
    colnames(Y1s)=paste0("Y1_",1:ncol(curZ))
    

  }
  
  
  ############################################################################
  ###################Interactive Effect#######################################
  ############################################################################
  

  if (type=="interactive"){
    
    
    #create a vector of nearest neighbor
    diag(dist_matrix_ZZ)=Inf
    
    #In the case there are multiple adjacent polygons, which.min returns the first polygon by index
    
    Zneighbor = apply(dist_matrix_ZZ,1,which.min)
    if (polygon==0){
      Zneighbor_distance = apply(dist_matrix_ZZ,1,min) 
      
    }

    if (polygon==1){

      #for each polygon: it is influenced by its nearest neighbors
      #the spehre of influnece is upper bounded by the distance between polygons 
      #plus the diameter of the *affected* polygon
      Zneighbor_distance = apply(dist_matrix_ZZ,1,min) + diameter
    }
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
  

  ############################################################################
  ###################Null Effect##############################################
  ############################################################################
  
  
  if (type=='null'){
    
    #initialize values
    Y1s=matrix(0,nrow(Ydata),ncol(curZ))
    
    for(k in 1:ncol(curZ)){

        Y1s[, k] <- Ydata$Y0
      }
    colnames(Y1s)=paste0("Y1_",1:ncol(curZ))
    

    }
  
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
    result_list <- list(Ydata=Ydata, Zdata=Zdata, ras0=ras0, ras_Z=ras_Z, Zneighbor_distance=Zneighbor_distance,
                        dist_matrix_ZZ=dist_matrix_ZZ,dist_matrix_YZ=dist_matrix_YZ)
    
  }else{
    result_list <- list(Ydata=Ydata, Zdata=Zdata, ras0=ras0, ras_Z=ras_Z,
                        dist_matrix_ZZ=dist_matrix_ZZ,dist_matrix_YZ=dist_matrix_YZ)
    
  }
  return(result_list)
}



dataGenerator_for_ploting <- function(nsim,scale, tr_prob, YZratio, unit_size,  cr = FALSE, pair=FALSE,polygon = 0, trim_ras = 0, effect = 3, Y_sd = 0.1, Zjitter = 1, type = "non-monotonic"){
  
  
  ##########################################################################
  s <- scale #length of one side
  s2 <- scale^2 #length of the square
  p <- round(s2*tr_prob)
  
  ##########################################################################
  ############Generate a raster#############################################
  ##########################################################################
  
  #generate a raster
  ras0 <- raster(	nrows=s, ncols=s, 
                  xmn=0, xmx=s*unit_size, 
                  ymn=0, ymx=s*unit_size, 
                  vals=rnorm(s2, mean=0, sd=Y_sd))
  
  # Dataset with Y0 values and then centroids of the raster cells
  Y0 <- getValues(ras0) #control outcomes
  Yxy <- coordinates(ras0) 
  Yindex <- 1:nrow(Yxy)
  Ydata <- data.frame(Yindex,cbind(Y0,Yxy))
  ny=length(Y0)
  
  #############################################################################
  ###########Generate alpha_x with Kriging#####################################
  #############################################################################
  FExy <- coordinates(aggregate(ras0, fact=(s/4), fun=mean)) #coarsened rasters
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
    
    #control outcome matrix
    Y0_matrix <- matrix(Y0, s, s)
    
    #add a coordinate reference system
    crs(ras0) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
    
    if ((trim_ras == 1) & (s >= 40)){
      
      #trim the rasters on the sides: this way we make sure that the intervention nodes and their sphere of influence are all within the rasters  
      sub_ras0 <- raster(nrows=(s-20),
                         ncols=(s-20),
                         xmn=11*unit_size, xmx=(s-10)*unit_size,
                         ymn=11*unit_size, ymx=(s-10)*unit_size,
                         vals=c(Y0_matrix[c(11:(s-10)), c(11:(s-10))]))
      
      
    }else{
      
      sub_ras0 <- ras0
      
    }
    
    #inherit the coordinate reference system
    crs(sub_ras0) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
    
    ras_polygon <- as(sub_ras0, 'SpatialPolygonsDataFrame')
    sf_raster <- st_as_sf(ras_polygon)
    
    ############################################################################
    ##################Subsample Points for Creating Polygon#####################
    ############################################################################    
    polygon_points <- st_sample(sf_raster, round(s2/10, 0)) %>%
      st_sf() %>%
      st_transform(crs = crs(sub_ras0))
    
    ############################################################################
    ##################Creating Polygons#########################################
    ############################################################################    
    all_polygons <- polygon_points %>%  # consider the master points
      st_geometry() %>% # ... as geometry only (= throw away the data items)
      st_union() %>% # unite them ...
      st_voronoi() %>% # ... and perform the voronoi tessellation
      st_collection_extract(type = "POLYGON") %>% # select the polygons
      st_sf(crs = crs(sub_ras0)) 
    
    all_polygons <- all_polygons[st_is_valid(all_polygons), ]
    all_polygons <- st_intersection(all_polygons, st_union(sf_raster))
    
    ############################################################################
    ##################Sample Intervention Polygons##############################
    ############################################################################  
    
    if (pair){
      
      #randomly sample some polygons (half the size of the number of intervention nodes)
      

      
      sample_points <- sample(round(s2/10, 0), (s2/YZratio^2)/2)
      
      index_polygons=1:nrow(all_polygons)
      index_to_include=c()
      
      for (i in sample_points){
        
        
        
        #distance between polytopes
        dist_temp=sapply( 1:nrow(all_polygons), function(j) return(st_distance(all_polygons[j,],all_polygons[i,])))
        
        #set distance to itself as infinity    
        dist_temp[i]=Inf
        
        #find the adjacent neighbors and randomly select one
        neigbor_temp=(1:nrow(all_polygons))[which(dist_temp<1e-8)]
        
        #randomly sample a pair
        i_pair_temp=sample(neigbor_temp,1)
        
        #append the pair of polygons
        index_to_include=c(index_to_include,i,i_pair_temp)
        
      }  
      
      #create rasters
      ras_Z = all_polygons[index_to_include,]
      
    }else{
      
      #randomly sample some polygons
      sample_points <- sample(round(s2/(YZratio^2/10), 0), (s2/YZratio^2))
      ras_Z <- all_polygons[sample_points,]
      
    }
    
    
    Zxy <- st_coordinates(st_centroid(ras_Z))
    index <- 1:nrow(Zxy)
    Zdata <- as.data.frame(cbind(index, Zxy))
    names(Zdata)[2:3] <- c("x", "y")
    
    
    nz <- dim(Zdata)[1]
    Ydata_coords <- st_as_sf(Ydata, crs = st_crs(ras_Z), coords = c("x", "y"))
    
    #distance matrix from evaluation points to intervention nodes
    dist_matrix_YZ <- st_distance(Ydata_coords, ras_Z, which = "Euclidean")
    Ydata=cbind(Ydata,dist_matrix_YZ)
    names(Ydata)[5:ncol(Ydata)] <- sapply(1:nz,function(i) return(paste0("dZ",i)))
    
    #distance matrix from intervention nodes to intervention nodes
    dist_matrix_ZZ <- st_distance(ras_Z, ras_Z, which = "Euclidean")
    Zdata=cbind(Zdata,dist_matrix_ZZ)
    names(Zdata)[4:ncol(Zdata)] <- sapply(1:nz,function(i) return(paste0("dZ",i)))
    
    #diamaeter of polygons
    
    diameter = rep(0,nrow(ras_Z))
    for (i in 1:nrow(ras_Z)){
      diameter[i]=ras_Z[i,] %>% 
        st_cast('MULTIPOINT') %>% 
        st_cast('POINT') %>% 
        st_distance %>% 
        max()
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
    
    if (pair){
      Zxy1 <- coordinates(ras_Z)
      Zxy1[, 1] <- jitter(Zxy1[, 1], amount=Zjitter)
      Zxy1[, 2] <- jitter(Zxy1[, 2], amount=Zjitter)
      sample_ind <- sort(sample(index, length(index)/2))
      Zxy <- t(matrix(c(t(cbind(Zxy[sample_ind, ], Zxy1[sample_ind, ]))), nrow = 2))
    }
    
    
    #create node data
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
    
    names(Zdata)[2:3] <- c("x", "y")
    
    ras_Z=NULL
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
    Y1s = effect_nonmo1(Y0,dist_matrix_YZ,curZ,alpha, effect,sh,  sc,  a)[['Sim_data']]
    colnames(Y1s)=paste0("Y1_",1:ncol(curZ))
    
    tau = calculate_tau(Y1s,curZ)
    Ydata =  cbind(Ydata,tau)
    
    #########################
    ####Calculate AME########
    #########################
    dVec=seq(0.5,10,by=0.25)
    AMEmat <- data.frame(d = dVec, tauda2 = NA, nd = NA)
    
    for(dIndex in 1:length(dVec)){
      #calculate true AME over d indices
      print(dIndex)
      dUp <- dVec[dIndex]
      
      AMEmat$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz=nrow(Zdata), numpts = NULL, only.unique = 0)
       
    }
    
    return(list(tau=AMEmat,ras0=ras0,Zdata=Zdata,ras_Z=ras_Z))    
  }
  
  
  ############################################################################
  ###################Interactive Effect#######################################
  ############################################################################
  
  
  if (type=="interactive"){
    
    
    #create a vector of nearest neighbor
    diag(dist_matrix_ZZ)=Inf
    
    #In the case there are multiple adjacent polygons, which.min returns the first polygon by index
    
    Zneighbor = apply(dist_matrix_ZZ,1,which.min)
    if (polygon==0){
      Zneighbor_distance = apply(dist_matrix_ZZ,1,min) 
      
    }
    
    if (polygon==1){
      
      #for each polygon: it is influenced by its nearest neighbors
      #the spehre of influnece is upper bounded by the distance between polygons 
      #plus the diameter of the *affected* polygon
      Zneighbor_distance = apply(dist_matrix_ZZ,1,min) + diameter
    }
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
    
    ###############################
    ##calculate effect function####
    ###############################
    
    tau_temp=calculate_tau_separate(Y1s,curZ,Zneighbor)
    tau_t=tau_temp[['tau_t']]
    tau_c=tau_temp[['tau_c']]
    
    tau = calculate_tau(Y1s,curZ)
    
    Ydata_t=  cbind(Ydata,tau_t)
    Ydata_c=  cbind(Ydata,tau_c)
    Ydata =  cbind(Ydata,tau)
    
    #########################
    ####Calculate AME########
    #########################
    dVec=seq(0.5,10,by=0.25)
    
    
    AMEmat_t <- data.frame(d = dVec, tauda2 = NA, nd = NA)
    AMEmat_c <- data.frame(d = dVec, tauda2 = NA, nd = NA)
    AMEmat <- data.frame(d = dVec, tauda2 = NA, nd = NA)
    
    for(dIndex in 1:length(dVec)){
      #calculate true AME over d indices
      print(dIndex)
      dUp <- dVec[dIndex]
      
      AMEmat$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata, Zdata, ras_Z = ras_Z, dUp, nz=nrow(Zdata), numpts = NULL, only.unique = 0)
      AMEmat_t$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata_t, Zdata, ras_Z = ras_Z, dUp, nz=nrow(Zdata), numpts = NULL, only.unique = 0)
      AMEmat_c$tauda2[dIndex] <- circleMean_oracle(ras0, Ydata_c, Zdata, ras_Z = ras_Z, dUp, nz=nrow(Zdata), numpts = NULL, only.unique = 0)
      
    }
    
    return(list(tau=AMEmat,tau_t=AMEmat_t,tau_c=AMEmat_c,ras0=ras0,Zdata=Zdata,ras_Z=ras_Z))
  }
  
  

  


}
