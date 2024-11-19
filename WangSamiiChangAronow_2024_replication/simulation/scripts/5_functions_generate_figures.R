generate_figures=function(data_files,sample_sizes,n_ss,data_address,est='regular',case){


  
load(paste0(data_address,data_file,'.Rdata'))
  
dVec = seq(0.5,10,0.25)

if (est=='regular'){
true_parameter_list = simulation_result[['output_AME_list']]
true_parameter_est_list=simulation_result[['output_AME_est_list']]
HAC_est_list = simulation_result[['output_all_VCEs_list']]
HAC_PD_est_list = simulation_result[['output_all_VCEs_pd_list']]
SAH_est_list = simulation_result[['output_all_SAH_list']]
edof_list = simulation_result[['output_all_edof_list']]
edof_sd_list = simulation_result[['output_all_edof_sd_list']]

}

if (est=='sm'){
  
#
true_parameter_list = simulation_result[['output_AME_sm_list']]
true_parameter_est_list=simulation_result[['output_AME_est_sm_list']]
HAC_est_list = simulation_result[['output_all_VCEs_sm_list']]
HAC_PD_est_list = simulation_result[['output_all_VCEs_pd_sm_list']]
SAH_est_list = simulation_result[['output_all_SAH_sm_list']]
edof_list = simulation_result[['output_all_edof_sm_list']]
edof_sd_list = simulation_result[['output_all_edof_sd_sm_list']]

}

####################################################################
#########Destination for Plots######################################
####################################################################
####################################################################
title_effect = paste0('./simulation/graph_new/',case,'_',est,'_',data_file,'_effect.pdf')
title_MSE= paste0('./simulation/graph_new/',case,'_',est,'_',data_file,'_MSE.pdf')
title_coverage= paste0('./simulation/graph_new/',case,'_',est,'_',data_file,'_coverage.pdf')
title_length=paste0('./simulation/graph_new/',case,'_',est,'_',data_file,'_length.pdf')
# title_coverage= paste0('./simulation/graph_new/',case,'_',est,'_',data_file,'_coverage_N=',sample_sizes[n_ss],'.pdf')
# title_length=paste0('./simulation/graph_new/',case,'_',est,'_',data_file,'_length_N=',sample_sizes[n_ss],'.pdf')
# title_coverage_all= paste0('./simulation/graph_new/',case,'_',est,'_',data_file,'_coverage_all.pdf')
# title_length_all= paste0('./simulation/graph_new/',case,'_',est,'_',data_file,'_length_all.pdf')


####################################################################
#########Extract Data###############################################
####################################################################
####################################################################

result_list=list()
num_neg_var_list=list()
for (ss in 1:length(sample_sizes)){
  
  #initilize a matrix to hold data
  result_matrix=matrix(0,nrow=length(dVec),ncol=16)
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
    result_matrix[d,3] =  mean ( abs(tau_est-tau) <= std1 * 1.96,na.rm=T) #HAC
    result_matrix[d,4] =  mean ( abs(tau_est-tau) <= std2 * qt(0.975,edof_pd),na.rm=T) #HAC_PD
    result_matrix[d,5] =  mean ( abs(tau_est-tau) <= std2 * 1.96,na.rm=T) #HAC_PD
    
    result_matrix[d,6] =  mean ( abs(tau_est-tau) <= std3 * 1.96,na.rm=T) #SAH
    
    #median half length
    result_matrix[d,7] =  median(std1*qt(0.975,edof),na.rm=T) #HAC
    result_matrix[d,8] =  median(std1*1.96,na.rm=T) #HAC
    

    result_matrix[d,9] =  median(std2*qt(0.975,edof_pd),na.rm=T) #HAC_PD
    result_matrix[d,10] =  median(std2*1.96,na.rm=T)
    result_matrix[d,11] =  median(std3*1.96,na.rm=T) #SAH
    
    #MSE
    result_matrix[d,12] = mean( (tau_est-tau)^2)
    result_matrix[d,13] = var(tau_est)
    result_matrix[d,14]= mean(HAC_est_list[[ss]][d,])
    result_matrix[d,15]= mean(HAC_PD_est_list[[ss]][d,])
    result_matrix[d,16]= mean(SAH_est_list[[ss]][d,])
  }
  result_list[[ss]] = result_matrix
  num_neg_var_list[[ss]]=num_neg_var_matrix
}

#############################################################################
#####################MSE##################################################
#############################################################################


#data frame for ploting effects
MSE_plot=as.data.frame(cbind(dVec,result_list[[n_ss]][,12]))
colnames(MSE_plot)=c('dVec','MSE')
MSE_plot = ggplot(MSE_plot,aes(dVec,MSE))+
  geom_smooth(se=FALSE) +
  xlab('Distance (d)')+
  ylab('MSE')+
  ggtitle('Mean Squared Error by Distance')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio=1)
ggsave(MSE_plot,file=title_MSE,height=4,width=7)
#development


#############################################################################
#####################Effect################################################
#############################################################################
effect_plot=as.data.frame(cbind(dVec,true_parameter_list[[n_ss]][,2]))
colnames(effect_plot)=c('dVec','Effect')
effect_plot = ggplot(effect_plot,aes(dVec,Effect))+
  geom_smooth(se=FALSE) +
  xlab('Distance (d)')+
  ylab('AME')+
  ggtitle('True AME by Distance')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(aspect.ratio=1)
ggsave(effect_plot,file=title_effect,height=4,width=7)


#############################################################################
#####################Coverage################################################
#############################################################################
cov_data=as.data.frame(cbind(dVec,result_list[[n_ss]][,2:6]))
names(cov_data)=c('dVec','HAC (edof)', 'HAC', 'HAC_PD (edof)','HAC_PD','SAH')

#reshape from wide to long
cov_data=melt(setDT(cov_data),id.vars='dVec',variable.name='var_type')

cov_plot=ggplot(cov_data,aes(dVec,value,group=var_type)) +
  geom_line(aes(color=var_type),linewidth=1)+
  geom_point(aes(shape=var_type),size=2)+
  scale_shape_manual(values=c(1,2,3,4,5),name='Variance Estimator')+
  ylim(0.7,1.05)+
  scale_y_continuous(breaks = sort(c(0.7,0.8,0.9,0.95,1.0)))+
  geom_hline(yintercept=0.95)+
  ylab('Coverage Rate')+
  xlab('Distance (d)')+
  ggtitle('Coverage Rate by Distance')+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(color="Variance Estimator")+
  theme(aspect.ratio=1)
ggsave(cov_plot,file=title_coverage,height=4,width=7)
  


#############################################################################
#####################Half Length#############################################
#############################################################################
length_data=as.data.frame(cbind(dVec,result_list[[n_ss]][,7:11]))
names(length_data)=c('dVec','HAC (edof)', 'HAC', 'HAC_PD (edof)','HAC_PD','SAH')
  
  #reshape from wide to long
length_data=melt(setDT(length_data),id.vars='dVec',variable.name='var_type')
  
length_plot=ggplot(length_data,aes(dVec,value,group=var_type)) +
    geom_line(aes(color=var_type),linewidth=1)+
  geom_point(aes(shape=var_type),size=2)+
  scale_shape_manual(values=c(1,2,3,4,5),name='Variance Estimator')+
  ylim(0,1.5)+
    ylab('Half Length')+
    xlab('Distance (d)')+
    ggtitle('Half Length by Distance')+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(color="Variance Estimator")+
   theme(aspect.ratio=1)
ggsave(length_plot,file=title_length,height=4,width=7)
   
# 
# #############################################################################
# #####################All HAC Coverage########################################
# #############################################################################
# n_color <- colorRampPalette(c("blue", "red"))(length(sample_sizes))
# pdf(file=title_coverage_all, height=8, width=8)
# par(mfrow=c(1,1))
# par(mar=c(5,5,5,5))
# par(pty="s")
# plot(0,
#      ylab=paste0("Coverage rate"),
#      xlab="Distance (d)",
#      type="n",
#      xlim=c(0, 10),
#      ylim=c(0, 1),
#      main="Coverage of the 95% CIs",
#      lty = 2, col = "black",
#      cex.lab=2,
#      cex.main=2)
# for (i in 1:length(sample_sizes)){
#   points(result_list[[i]][,2]~dVec, type = "l", col = n_color[i], lwd = 2)
# }
# abline(h = 0.95, lty = 2)
# legend("topright", col = n_color, lty=c(1,1), lwd = 2,
#        legend = sample_sizes, cex=2, bty = "n")
# dev.off()


# 
# #############################################################################
# #####################ALL HAC Half Length#####################################
# #############################################################################
# n_color <- colorRampPalette(c("blue", "red"))(length(sample_sizes))
# pdf(file=title_length_all, height=8, width=8)
# par(mfrow=c(1,1))
# par(mar=c(5,5,5,5))
# par(pty="s")
# plot(0,
#      ylab=paste0("Length"),
#      xlab="Distance (d)",
#      type="n",
#      xlim=c(0, 10),
#      ylim=c(0, 1),
#      main="Median Half Length of the 95% CIs ",
#      lty = 2, col = "black",
#      cex.lab=2,
#      cex.main=2)
# for (i in 1:length(sample_sizes)){
#   points(result_list[[i]][,5]~dVec, type = "l", col = n_color[i], lwd = 2)
# }
# abline(h = 0.95, lty = 2)
# legend("topright", col = n_color, lty=c(1,1), lwd = 2,
#        legend = sample_sizes, cex=2, bty = "n")
# dev.off()


return(list(result_list,num_neg_var_list))
}


output_MSE=function(data_files,sample_sizes,n_ss,data_address,est='regular',case){
  
  
  
  load(paste0(data_address,data_file,'.Rdata'))
  
  dVec = seq(0.5,10,0.25)
  
  if (est=='regular'){
    true_parameter_list = simulation_result[['output_AME_list']]
    true_parameter_est_list=simulation_result[['output_AME_est_list']]
    HAC_est_list = simulation_result[['output_all_VCEs_list']]
    HAC_PD_est_list = simulation_result[['output_all_VCEs_pd_list']]
    SAH_est_list = simulation_result[['output_all_SAH_list']]
    edof_list = simulation_result[['output_all_edof_list']]
    edof_sd_list = simulation_result[['output_all_edof_sd_list']]
    
  }
  
  if (est=='sm'){
    
    #
    true_parameter_list = simulation_result[['output_AME_sm_list']]
    true_parameter_est_list=simulation_result[['output_AME_est_sm_list']]
    HAC_est_list = simulation_result[['output_all_VCEs_sm_list']]
    HAC_PD_est_list = simulation_result[['output_all_VCEs_pd_sm_list']]
    SAH_est_list = simulation_result[['output_all_SAH_sm_list']]
    edof_list = simulation_result[['output_all_edof_sm_list']]
    edof_sd_list = simulation_result[['output_all_edof_sd_sm_list']]
    
  }
  

  
  ####################################################################
  #########Extract Data###############################################
  ####################################################################
  ####################################################################

  #initilize a matrix to hold data
  result_matrix=c()
  
  for (ss in 1:length(sample_sizes)){
    

    for (d in 1:length(dVec)){
      
      
      tau=true_parameter_list[[ss]][d,2]
      
      tau_est=true_parameter_est_list[[ss]][d,] 
      
      entry_temp = c(dVec[d], mean( (tau_est-tau)^2),sample_sizes[ss])
      #MSE
      result_matrix = rbind(result_matrix,entry_temp)
    }
    
  }
  
  #############################################################################
  #####################MSE##################################################
  #############################################################################
  

  
   #data frame for ploting MSE
   MSE_plot=as.data.frame(result_matrix)
   colnames(MSE_plot)=c('dVec','MSE','Sample_Size')
   MSE_plot$Sample_Size=as.factor(MSE_plot$Sample_Size)
   
   #plotting
   title_MSE= paste0('./simulation/graph_new/',case,'_',est,'_',data_file,'_MSE_all_cases.pdf')
   
   MSE_plot=ggplot(MSE_plot,aes(dVec,MSE,group=Sample_Size))+
    geom_line(aes(color=Sample_Size,linetype=Sample_Size)) +
    xlab('Distance (d)')+
    ylab('MSE')+
    ggtitle('Mean Squared Error by Distance')+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(aspect.ratio=1)
    ggsave(MSE_plot,file=title_MSE,height=4,width=7)

  return(result_matrix)
}
