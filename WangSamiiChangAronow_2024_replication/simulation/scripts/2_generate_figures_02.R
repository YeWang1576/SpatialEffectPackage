rm(list=ls())

library(ggplot2)
library(reshape2) #for melting function
###############################################################################
###################Plot AME Curve #############################################
###############################################################################

root_dir = 'C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/SpatialReplication'
setwd(root_dir)
###############################
#######Interactive Case########
###############################
rm(list=ls())

#load data
dVec=seq(0.5,10,by=0.25)
load('simulation/data_new/Paired_simData_interactive_10_Jun_06_point_ploting.RData')
tau_t=data_list[[1]]["tau_t"]$tau_t$tauda2
tau_c=data_list[[1]]['tau_c']$tau_c$tauda2
tau=data_list[[1]]['tau']$tau$tauda2

#make data into the long format
data_for_plot=data.frame(cbind(dVec,tau,tau_t,tau_c))
names(data_for_plot)=c('d','tau','tau_t','tau_c')
data_for_plot_long <- melt(data_for_plot, id = c("d")) 


#ploting
today_date = format(Sys.Date(), format="%b_%d")
title_interactive = paste0('simulation/graph_new/AME_interactive_', today_date,'.png') #figure name

plot_AME_interactive=ggplot(data_for_plot_long, aes(x=d,y=value)) + 
  geom_line(aes(linetype = variable)) +
  xlab('Distance (d)')+
  ylab('Effect Size')+
  labs(title='AME Effect Curve (Interactive)')+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))+
  scale_linetype_manual(values=c("solid","longdash", "dotted"), label=c('AME','Neighbor Treated','Neighbor Control'))

ggsave(plot_AME_interactive,file=title_interactive,height=4,width=7)

###############################
#######Additive Case###########
###############################
rm(list=ls())

#load data
dVec=seq(0.5,10,by=0.25)
load('simulation/data_new/Paired_simData_nonmono_10_Jun_06_point_ploting.RData')
tau=data_list[[1]]['tau']$tau$tauda2

#make data into data frame
data_for_plot=data.frame(cbind(dVec,tau,rep(1,length(dVec))))
data_for_plot[,3]=as.factor(data_for_plot[,3])
names(data_for_plot)=c('d','tau','value')


today_date = format(Sys.Date(), format="%b_%d")
title_additive = paste0('simulation/graph_new/AME_nonmo_', today_date,'.png') #figure name

plot_AME_nonmono=ggplot(data_for_plot, aes(x=d,y=tau,group=value)) + 
  geom_line(aes(linetype=value))+
  xlab('Distance (d)')+
  ylab('Effect Size')+
  labs(title='AME Effect Curve (Additive)',color="AME")+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))+
  scale_linetype_manual(values=c("solid"), label=c('AME'))

ggsave(plot_AME_nonmono,file=title_additive,height=4,width=6)


rm(list=ls())

library(ggplot2)
library(reshape2)
###############################################################################
###################Plot Intervention Nodes ####################################
###############################################################################


###############################
#######Polygon Case###########
###############################
rm(list=ls())

today_date = format(Sys.Date(), format="%b_%d")
title_polygon = paste0('simulation/graph_new/polygon_', today_date,'.png') #figure name

#load data
load('simulation/data_new/Paired_simData_nonmono_10_Jun_06_polygon_ploting.RData')
ras0=data_list[[1]]['ras0']$ras0
ras_Z=data_list[[1]]['ras_Z']$ras_Z

png(filename=title_polygon)
plot(ras0)
plot(ras_Z,add=TRUE)
dev.off()

###############################
#######Point Case#############
###############################
rm(list=ls())

today_date = format(Sys.Date(), format="%b_%d")
title_polygon = paste0('simulation/graph_new/point_', today_date,'.png') #figure name

#load data
load('simulation/data_new/Paired_simData_nonmono_10_Jun_06_point_ploting.RData')
ras0=data_list[[1]]['ras0']$ras0
point_data=data_list[[1]]['Zdata']$Zdata[,c(2,3)]

png(filename=title_polygon)
plot(ras0)
points(point_data)
dev.off()


################################################################################
###############Randomization Test###############################################
################################################################################

dVec <- seq(from=.5, to=10, by=.25)


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

load('simulation/data_new/Paired_randomization_test_nonmono_10_Jun_12_points.RData')
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

