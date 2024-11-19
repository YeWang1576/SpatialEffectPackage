monoEffect <- function(r, a){exp(-(r/a)^2)}
nonmonoEffect <- function(d, sh, sc, a){(dgamma(abs(d), shape=sh, scale=2*sc) - a*dgamma(abs(d), shape=5*sh, scale=sc))}

# functions to calculate the circle average
circleMean_oracle <- function(ras0, Ydata, Zdata, ras_Z = NULL, dUp, nz, numpts = NULL, gridRes = res(ras0)[1], evalpts = 10, only.unique = 0){
  
  dUpMean <- 0
  if (is.null(numpts)){
    numpts <- ceiling(ceiling(dUp/(gridRes/sqrt(2))+1)*pi)*evalpts #discretize evaluation points on a circle
    
     }
  for (i in 1:nz){
    if (is.null(ras_Z)){
      
      #construct coordinates for each i
      coords <- cbind(Zdata$x[i] + cos(((1:numpts)/numpts)*2*pi) * dUp, 
                      Zdata$y[i] + sin(((1:numpts)/numpts)*2*pi) * dUp)
    }else{
      one_buffer <- buffer(as_Spatial(ras_Z[i, ]), width=dUp)
      one_buffer_raster <- boundaries(rasterize(one_buffer, ras0))
      one_buffer_raster_values <- getValues(one_buffer_raster)
      coords <- coordinates(one_buffer_raster)[!is.na(one_buffer_raster_values) & !is.nan(one_buffer_raster_values) & one_buffer_raster_values == 1, ]
    }
    
    if (only.unique == 1){
      grid_index <- c(unique(na.omit(cellFromXY(ras0, coords))))
    }else{
      grid_index <- c(na.omit(cellFromXY(ras0, coords)))
    }
    grids <- Ydata[grid_index, (4+nz+i)] #these Ydata are on the circle
    dUpMean <- dUpMean + mean(grids)
  }
  dUpMean <- dUpMean / nz #this is equation (7) of the paper
  return(dUpMean)
}


circleMean_oracle_perz <- function(ras0, Ydata, Zdata, ras_Z = NULL, dUp, nz, numpts = NULL, gridRes = res(ras0)[1], evalpts = 10, only.unique = 0){
  
  
  dUpMean <- rep(NA,nz)
  if (is.null(numpts)){
    numpts <- ceiling(ceiling(dUp/(gridRes/sqrt(2))+1)*pi)*evalpts #discretize evaluation points on a circle
  }
  for (i in 1:nz){
    
    assignment_temp=Zdata[which(Zdata$index==i),4:ncol(Zdata)] #this is the assignment variable associated with the ith node
    
    if (is.null(ras_Z)){
      
      #construct coordinates for each i
      coords <- cbind(Zdata$x[i] + cos(((1:numpts)/numpts)*2*pi) * dUp, 
                      Zdata$y[i] + sin(((1:numpts)/numpts)*2*pi) * dUp)
    }else{
      one_buffer <- buffer(as_Spatial(ras_Z[i, ]), width=dUp)
      one_buffer_raster <- boundaries(rasterize(one_buffer, ras0))
      one_buffer_raster_values <- getValues(one_buffer_raster)
      coords <- coordinates(one_buffer_raster)[!is.na(one_buffer_raster_values) & !is.nan(one_buffer_raster_values) & one_buffer_raster_values == 1, ]
    }
    
    if (only.unique == 1){
      grid_index <- c(unique(na.omit(cellFromXY(ras0, coords))))
    }else{
      grid_index <- c(na.omit(cellFromXY(ras0, coords)))
    }
    grids <- Ydata[grid_index, (4+nz+i)] #these Ydata are on the circle
    
    
    dUpMean[i]=mean(grids) #This is equation (6) of the paper
  }
  return(dUpMean)
}

# it calculates the coordinates of points located on each circle
# we only need to calculate these coordinates once
circleMean <- function(Zdata,ras0, ras_Z = NULL, dVec, nz, numpts = NULL, gridRes = res(ras0)[1], evalpts = 10, only.unique = 0){
  
  #create a grid list for each node at different distances
  grid_list <- list()
  for (i in 1:nz){
    print(i)
    for (dIndex in 1:length(dVec)){
      dUp <- dVec[dIndex]
      print(dUp)
      if (is.null(numpts)){
        nts <- ceiling(ceiling(dUp/(gridRes/sqrt(2))+1)*pi)*evalpts
      } else{
        nts <- numpts
      }
      if (is.null(ras_Z)){
        coords <- cbind(Zdata$x[i] + cos(((1:nts)/nts)*2*pi) * dUp, 
                        Zdata$y[i] + sin(((1:nts)/nts)*2*pi) * dUp)
      }else{

        one_buffer <- buffer(as_Spatial(ras_Z[i, ]), width=dUp)
        one_buffer_raster <- boundaries(rasterize(one_buffer, ras0))
        one_buffer_raster_values <- getValues(one_buffer_raster)
        coords <- coordinates(one_buffer_raster)[(!is.na(one_buffer_raster_values)) & (one_buffer_raster_values == 1), ]
      }
      if (only.unique == 1){
        grid_index <- c(unique(na.omit(cellFromXY(ras0, coords))))
      }else{
        grid_index <- c(na.omit(cellFromXY(ras0, coords)))
      }
      grid_list[[(i-1)*length(dVec)+dIndex]] <- grid_index
      
    }
  }
  return(grid_list)
}
