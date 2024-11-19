rm(list=ls())

library(ggplot2)
library(reshape2) #for melt function
###############################################################################
###################Plot AME Curve #############################################
###############################################################################

root_dir = 'C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/SpatialReplication'
setwd(root_dir)

dVec=seq(0.5,10,by=0.25)
################################################################################
###############Randomization Test###############################################
################################################################################

######################
#Null Case ###########
######################

figure_title_null =paste0('simulation/graph_new/Randomization_test_null_', today_date,'.png') #figure name

load('simulation/data_new/Paired_randomization_test_null_10_Jun_12_points.RData')
temp=simulation_result[['output_rt_list']]
cov_64=apply(temp[[1]],1,mean)
cov_100=apply(temp[[2]],1,mean)
cov_144=apply(temp[[3]],1,mean)

plot_data = data.frame(cbind(dVec,cov_64,cov_100,cov_144))
colnames(plot_data)=c('d','64','100','144')
plot_data=melt(plot_data,id.vars='d')
colnames(plot_data)=c('d','sample_size','cov')

plot_null=ggplot(plot_data,aes(x=d,y=cov,group=sample_size))+
  geom_line(aes(linetype=sample_size,color=sample_size))+
  ylim(0,1)+
  ylab('Rejection Probability')+
  ggtitle('Rejection Probability (Null Effect)')+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))



ggsave(plot_null,file=figure_title_null,height=4,width=7)

######################
#Interactive Case ####
######################

figure_title_interactive =paste0('simulation/graph_new/Randomization_test_interactive_', today_date,'.png') #figure name

load('simulation/data_new/Paired_randomization_test_interactive_10_Jun_12_points.RData')
temp=simulation_result[['output_rt_list']]
cov_64=apply(temp[[1]],1,mean)
cov_100=apply(temp[[2]],1,mean)
cov_144=apply(temp[[3]],1,mean)

plot_data = data.frame(cbind(dVec,cov_64,cov_100,cov_144))
colnames(plot_data)=c('d','64','100','144')
plot_data=melt(plot_data,id.vars='d')
colnames(plot_data)=c('d','sample_size','cov')

plot_interactive=ggplot(plot_data,aes(x=d,y=cov,group=sample_size))+
  geom_line(aes(linetype=sample_size,color=sample_size))+
  ylim(0,1)+
  ylab('Rejection Probability')+
  ggtitle('Rejection Probability (Interactive Effect)')+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))


ggsave(plot_interactive,file=figure_title_interactive,height=4,width=7)


######################
#Interactive Case ####
######################


figure_title_nonmono =paste0('simulation/graph_new/Randomization_test_nonmono_', today_date,'.png') #figure name

temp=simulation_result[['output_rt_list']]
cov_64=apply(temp[[1]],1,mean)
cov_100=apply(temp[[2]],1,mean)
cov_144=apply(temp[[3]],1,mean)

plot_data = data.frame(cbind(dVec,cov_64,cov_100,cov_144))
colnames(plot_data)=c('d','64','100','144')
plot_data=melt(plot_data,id.vars='d')
colnames(plot_data)=c('d','sample_size','cov')

plot_nonmono=ggplot(plot_data,aes(x=d,y=cov,group=sample_size))+
  geom_line(aes(linetype=sample_size,color=sample_size))+
  ylim(0,1)+
  ylab('Rejection Probability')+
  ggtitle('Rejection Probability (Additive Effect)')+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))

ggsave(plot_nonmono,file=figure_title_nonmono,height=4,width=7)

################################################################################
###############Randomization Test: sup norm#####################################
################################################################################

rej_probability=matrix(NA,3,3)
colnames(rej_probability)=c('Null','Nonmono','Interactive')
rownames(rej_probability)=c(64,100,144)

######################
#Null Case ###########
######################

load('simulation/data_new/Paired_randomization_test_sup_null_10_Jun_14_points.RData')
temp=simulation_result[['output_rt_list']]
rej_probability[1,1]=mean(temp[[1]])
rej_probability[2,1]=mean(temp[[2]])
rej_probability[3,1]=mean(temp[[3]])


load('simulation/data_new/Paired_randomization_test_sup_nonmono_10_Jun_14_points.RData')
temp=simulation_result[['output_rt_list']]
rej_probability[1,2]=mean(temp[[1]])
rej_probability[2,2]=mean(temp[[2]])
rej_probability[3,2]=mean(temp[[3]])

load('simulation/data_new/Paired_randomization_test_sup_interactive_10_Jun_14_points.RData')
temp=simulation_result[['output_rt_list']]
rej_probability[1,3]=mean(temp[[1]])
rej_probability[2,3]=mean(temp[[2]])
rej_probability[3,3]=mean(temp[[3]])

