GenSamplingPoints <- function(center, radius, numpts, dist.metric = "Euclidean"){
  if (dist.metric == "Euclidean"){
    return(cbind(center[1] + cos(((1:numpts)/numpts)*2*pi) * radius,
                 center[2] + sin(((1:numpts)/numpts)*2*pi) * radius))
  }else if (dist.metric == "Geodesic"){
    return(destPoint(center, b = ((1:numpts)/numpts)*360, d = radius)) 
  }
}

RimAvg <- function(ras, ras_coords = NULL, Yobs = NULL, Zdata, x_coord_Z, y_coord_Z, dUp, 
                   numpts = NULL, gridRes = NULL, evalpts = 10, only.unique = 0, 
                   dtype = "raster", dist.metric = dist.metric){

  nz <- nrow(Zdata)
  Ybard <- rep(NA, nz)
  if (dtype == "krig"){
    if (is.null(numpts)){
      numpts <- ceiling(ceiling(dUp/(gridRes/sqrt(2))+1)*pi)*evalpts
    }
    S.rims <- matrix(NA, nz*numpts, 3)
    dist.rims <- matrix(0, nz*numpts, nz)
    for (i in 1:nz) {
      coords <- GenSamplingPoints(center = c(Zdata[i, x_coord_Z], Zdata[i, y_coord_Z]), 
                                  radius = dUp, numpts = numpts, dist.metric = dist.metric)
      # coords <- cbind(Zdata[i, x_coord_Z] + cos(((1:numpts)/numpts)*2*pi) * dUp,
      #                 Zdata[i, y_coord_Z] + sin(((1:numpts)/numpts)*2*pi) * dUp)
      grids <- predict(ras, coords)
      S.rim <- cbind(grids, coords)
      index.rim <- c(((i - 1) * numpts + 1) : (i * numpts))
      S.rims[index.rim, ] <- S.rim
      dist.rim <- matrix(0, dim(S.rim)[1], nz)
      dist.rim[, i] <- dUp
      dist.rims[index.rim, ] <- dist.rim
      Ybard[i] <- mean(grids)
    }
  }else{
    if (is.null(numpts)){
      if (is.null(gridRes)){
        gridRes <- res(ras)
      }
      numpts <- ceiling(ceiling(dUp/(gridRes/sqrt(2))+1)*pi)*evalpts
    }
    S.rims <- data.frame()
    dist.rims <- data.frame()
    S.rims.list <- list()
    dist.rims.list <- list()
    for (i in 1:nz){
      coords <- GenSamplingPoints(center = c(Zdata[i, x_coord_Z], Zdata[i, y_coord_Z]), 
                                  radius = dUp, numpts = numpts, dist.metric = dist.metric)
      # 
      # coords <- cbind(Zdata[i, x_coord_Z] + cos(((1:numpts)/numpts)*2*pi) * dUp,
      #                 Zdata[i, y_coord_Z] + sin(((1:numpts)/numpts)*2*pi) * dUp)
      if (only.unique == 1){
        grid_index <- c(unique(na.omit(cellFromXY(ras, coords))))
      }else{
        grid_index <- c(na.omit(cellFromXY(ras, coords)))
      }
      grids <- Yobs[grid_index]
      if (only.unique == 1){
        S.rim <- cbind(grids, ras_coords[unique(na.omit(cellFromXY(ras, coords))), ])
      }else{
        S.rim <- cbind(grids, ras_coords[na.omit(cellFromXY(ras, coords)), ])
      }
      S.rims.list[[i]] <- S.rim[!is.na(S.rim[, 1]), ]
      dist.rim <- matrix(0, dim(S.rim)[1], nz)
      dist.rim[, i] <- dUp
      dist.rims.list[[i]] <- dist.rim[!is.na(S.rim[, 1]), ]
      Ybard[i] <- mean(grids, na.rm = TRUE)
    }
    S.rims <- do.call("rbind", S.rims.list)
    dist.rims <- do.call("rbind", dist.rims.list)
    rm(S.rims.list)
    rm(dist.rims.list)
  }
  Rim.list <- list(Ybard, S.rims, dist.rims)
  return(Rim.list)
}

# this function generates regression representation for any given dataset and assignment
RegRep <- function(Sdata, dist_mat, dVec, Zup){

  nz <- length(Zup)
  # SD <- cbind(Sdata, dist_mat, rep(1, nrow(Sdata)))
  # names(SD)[ncol(SD)] <- "cons"
  # SDcol <- ncol(SD)
  # SD.agg.formula <- paste0("cons ~ ", paste0(names(SD)[-ncol(SD)], collapse = "+"))
  # SD.agg <- aggregate(as.formula(SD.agg.formula), SD, length)
  # SD <- SD.agg #SD.agg[, c(1:(nz+2))]
  # names(SD)[ncol(SD)] <- "weight"
  # Sdata.new <- SD.agg[, 1:3]
  # dist_mat.new <- SD[, c(4:(nz+3))]
  # weights <- diag(unlist(SD[ncol(SD)]))
  # weights_2 <- weights^2
  # weighted_Y <- weights %*% Sdata.new[, 1]

  dist_mat[dist_mat > dVec[length(dVec)]] <- 0
  p <- mean(Zup)
  Zup_normal <- (Zup - p) / (p * (1-p))
  W_mat <- matrix(0, nrow(dist_mat), length(dVec))
  tau <- rep(NA, length(dVec))

  for (l in 1:length(dVec)){
    # if (mean(diag(weights) %in% 1) == 1){
    #   dd_temp <- d_temp <- (dist_mat == dVec[l])
    # }else{
    dd_temp <- d_temp <- (dist_mat == dVec[l])
    # }
    for (k in 1:ncol(d_temp)){
      if (sum(d_temp[, k]) == 0){
        dd_temp[, k] <- 0
      }else{
        dd_temp[, k] <- dd_temp[, k] / sum(dd_temp[, k])
      }
    }
    X_temp <- dd_temp %*% Zup_normal
    tau[l] <- t(X_temp) %*% Sdata[, 1]
    W_mat[, l] <- X_temp
  }

  #W_mat <- W_mat / nz
  X_mat <- nz * W_mat %*% solve(t(W_mat) %*% W_mat)
  coefs <- solve(t(X_mat) %*% X_mat) %*% (t(X_mat) %*% Sdata[, 1])
  res <- Sdata[, 1] - X_mat %*% coefs
  reg_list <- list(Sdata, X_mat, coefs, res)
  return(reg_list)
}

