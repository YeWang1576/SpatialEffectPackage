#This script requires the input of mainText.R, cincluding the variables
# 1. Zdata
# 2. AMEmat
# 3. AMEmat_perz
# 4. dVec


dbar=6
homophily=rep(0,length(dVec))

for (d in 1:length(dVec)){
  print(d)
  print(mean(AMEmat_perz[,d]))
  print(AMEmat[d,2])
  
  if (abs(AMEmat[d,2]-mean(AMEmat_perz[,d]))>1e-6){
    print('Significant Differences')
  }

}


calculate_homophily = function(d,dbar,Zdata,AMEmat_perz,AMEmat,dVec){
  
  
  #d: distance of the effect of intereset
  #dbar: maximum distance of influence
  #Zdata: data on intervention notes
  #AMEmat_perz: AME on the node level
  #AMEmat: overall AME
  nz=nrow(Zdata)
  
  #AME at distance d
  AMEd= AMEmat[d,2]
  sum=0
  
  #threshold
  hd= 2 * (dbar + dVec[d])

  for (i in 1:nz){
    
    x_coord=Zdata[i,'x'] #x-coordinate of ith node
    y_coord=Zdata[i,'y'] #y-coordinate of ith node
    
    dist_temp= sqrt((Zdata[,'x']-x_coord)^2 + (Zdata[,'y']-y_coord)^2)
    
    AME_neigh = AMEmat_perz[ which(dist_temp<=hd), d]
    
    sum = sum + (AMEmat_perz[i,d]-AMEd) * sum( AME_neigh-AMEd)
  }
  
  return(sum/nz^2)
}


for (d in 1:length(dVec)){
  
  homophily[d]=calculate_homophily(d,dbar,Zdata,AMEmat_perz,AMEmat,dVec)
  print(paste0('Homophily value at ', dVec[d], ' is:',calculate_homophily(d,dbar,Zdata,AMEmat_perz,AMEmat,dVec) ))
  
}
