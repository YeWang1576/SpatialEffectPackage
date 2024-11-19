# simulation results in Wang, Samii, Chang, and Aronow (2024)
rm(list=ls())


root_dir = "C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/SpatialReplication/"
setwd(root_dir)

#Old data 
#root_dir = "C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/"

#data directory
data_address = "./simulation/data_new/"
data_file="Paired_simulation_nonmono_10_May_28_points"

n_ss=5
sample_sizes <- c(16,36,64,100,144)

# n_ss=1
# sample_sizes <- c(729)

#load data
load(paste0(data_address,data_file,'.Rdata'))


#title=paste0("AME-CI-coverage-", type, "-bernoulli-smooth.pdf")
dVec = seq(0.5,10,0.25)

true_parameter_list = simulation_result[['output_AME_list']]
true_parameter_est_list=simulation_result[['output_AME_est_list']]
HAC_est_list = simulation_result[['output_all_VCEs_list']]
HAC_PD_est_list = simulation_result[['output_all_VCEs_pd_list']]
SAH_est_list = simulation_result[['output_all_SAH_list']]
edof_list = simulation_result[['output_all_edof_list']]
edof_sd_list = simulation_result[['output_all_edof_sd_list']]



# # 
# true_parameter_list = simulation_result[['output_AME_sm_list']]
# true_parameter_est_list=simulation_result[['output_AME_est_sm_list']]
# HAC_est_list = simulation_result[['output_all_VCEs_sm_list']]
# HAC_PD_est_list = simulation_result[['output_all_VCEs_pd_sm_list']]
# SAH_est_list = simulation_result[['output_all_SAH_sm_list']]
# edof_sm_list = simulation_result[['output_all_edof_sm']]

####################################################################
#########Destination for Plots######################################
####################################################################
####################################################################
title_MSE= paste0('./simulation/graph_new/',data_file,'_MSE.pdf')
title_coverage= paste0('./simulation/graph_new/',data_file,'_coverage_N=',sample_sizes[n_ss],'.pdf')
title_length=paste0('./simulation/graph_new/',data_file,'_length_N=',sample_sizes[n_ss],'.pdf')
title_coverage_all= paste0('./simulation/graph_new/',data_file,'_coverage_all.pdf')
title_length_all= paste0('./simulation/graph_new/',data_file,'_length_all.pdf')

####################################################################
#########Extract Data###############################################
####################################################################
####################################################################
result_list=list()
num_neg_var_list=list()
for (ss in 1:length(sample_sizes)){
  
  #initilize a matrix to hold data
  result_matrix=matrix(0,nrow=length(dVec),ncol=12)
  result_matrix[,1]=dVec
  
  #num of negative standard errors
  num_neg_var_matrix=matrix(0,nrow=length(dVec),ncol=2)
  num_neg_var_matrix[,1]=dVec
  
  for (d in 1:length(dVec)){
    
    
    tau=true_parameter_list[[ss]][d,2]
    
    tau_est=true_parameter_est_list[[ss]][d,] 
    
    std1=sqrt(HAC_est_list[[ss]][d,])
    std2=sqrt(HAC_PD_est_list[[ss]][d,])
    std3=sqrt(SAH_est_list[[ss]][d,])
    edof=edof_list[[ss]][d,]
    edof_pd=edof_sd_list[[ss]][d,]
    
    #for those with edof <=1, make it 1
    edof[which(edof<=1)]=1
    edof_pd[which(edof_pd<=1)]=1
    
    
    num_neg_var_matrix[d,2]=sum(is.na(std1))
    
    #coverages
    result_matrix[d,2] =  mean ( abs(tau_est-tau) <= std1 * qt(0.975,edof),na.rm=T) #HAC
    #result_matrix[d,2] =  mean ( abs(tau_est-tau) <= std1 * 1.96,na.rm=T) #HAC
    #print('Do not forget to change back')
    result_matrix[d,3] =  mean ( abs(tau_est-tau) <= std2 * qt(0.975,edof_pd),na.rm=T) #HAC_PD
    result_matrix[d,4] =  mean ( abs(tau_est-tau) <= std3 * 1.96,na.rm=T) #SAH
    
    #median half length
    result_matrix[d,5] =  median(std1*qt(0.975,edof),na.rm=T) #HAC
    #result_matrix[d,5] =  mean(std1*1.96,na.rm=T) #HAC
    
    result_matrix[d,6] =  median(std2*qt(0.975,edof_pd),na.rm=T) #HAC_PD
    result_matrix[d,7] =  median(std3*1.96,na.rm=T) #SAH
    
    #MSE
    result_matrix[d,8] = mean( (tau_est-tau)^2)
    result_matrix[d,9] = var(tau_est)
    result_matrix[d,10]= mean(HAC_est_list[[ss]][d,])
    result_matrix[d,11]= mean(HAC_PD_est_list[[ss]][d,])
    result_matrix[d,12]= mean(SAH_est_list[[ss]][d,])
  }
  result_list[[ss]] = result_matrix
  num_neg_var_list[[ss]]=num_neg_var_matrix
}



#############################################################################
#####################MSE################################################
#############################################################################
n_color <- colorRampPalette(c("blue", "red"))(length(sample_sizes))
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(0,
     ylab=paste0("Mean Square Errors"),
     xlab="Distance (d)",
     type="n",
     xlim=c(0, 10),
     ylim=c(0, 0.11),
     main="Mean Square Errors",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2)
for (i in 1:length(sample_sizes)){
  points(result_list[[i]][,8]~dVec, type = "l", col = n_color[i], lwd = 2)
}
abline(h = 0.95, lty = 2)
legend("topright", col = n_color, lty=c(1,1), lwd = 2,
       legend = sample_sizes, cex=2, bty = "n")


#############################################################################
#####################Coverage################################################
#############################################################################
var_type_color <- colorRampPalette(c("blue", "red"))(3)
var_type = c('HAC','HAC_PD','SAH')
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(0,
     ylab=paste0("Coverage rate"),
     xlab="Distance (d)",
     type="n",
     xlim=c(0, 10),
     ylim=c(0, 1),
     main="Coverage of the 95% CIs",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2)
for (i in 1:3){
  points(result_list[[n_ss]][,i+1]~dVec, type = "l", col = var_type_color[i], lwd = 2)
}
abline(h = 0.95, lty = 2)
legend("bottomleft", col = var_type_color, lty=c(1,1), lwd = 2,
       legend = var_type, cex=2, bty = "n")



#############################################################################
#####################Half Length#############################################
#############################################################################
var_type_color <- colorRampPalette(c("blue", "red"))(3)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(0,
     ylab=paste0("Length"),
     xlab="Distance (d)",
     type="n",
     xlim=c(0, 10),
     ylim=c(0, 1),
     main="Median Half Length of the 95% CIs",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2)
for (i in 1:3){
  points(result_list[[n_ss]][,i+4]~dVec, type = "l", col = var_type_color[i], lwd = 2)
}
abline(h = 0.95, lty = 2)
legend("bottomleft", col = var_type_color, lty=c(1,1), lwd = 2,
       legend = var_type, cex=2, bty = "n")


#############################################################################
#####################All HAC Coverage########################################
#############################################################################
n_color <- colorRampPalette(c("blue", "red"))(length(sample_sizes))
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(0,
     ylab=paste0("Coverage rate"),
     xlab="Distance (d)",
     type="n",
     xlim=c(0, 10),
     ylim=c(0, 1),
     main="Coverage of the 95% CIs",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2)
for (i in 1:length(sample_sizes)){
  points(result_list[[i]][,2]~dVec, type = "l", col = n_color[i], lwd = 2)
}
abline(h = 0.95, lty = 2)
legend("topright", col = n_color, lty=c(1,1), lwd = 2,
       legend = sample_sizes, cex=2, bty = "n")



#############################################################################
#####################ALL HAC Half Length#####################################
#############################################################################
n_color <- colorRampPalette(c("blue", "red"))(length(sample_sizes))
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
par(pty="s")
plot(0,
     ylab=paste0("Length"),
     xlab="Distance (d)",
     type="n",
     xlim=c(0, 10),
     ylim=c(0, 1),
     main="Median Half Length of the 95% CIs ",
     lty = 2, col = "black",
     cex.lab=2,
     cex.main=2)
for (i in 1:length(sample_sizes)){
  points(result_list[[i]][,5]~dVec, type = "l", col = n_color[i], lwd = 2)
}
abline(h = 0.95, lty = 2)
legend("topright", col = n_color, lty=c(1,1), lwd = 2,
       legend = sample_sizes, cex=2, bty = "n")
