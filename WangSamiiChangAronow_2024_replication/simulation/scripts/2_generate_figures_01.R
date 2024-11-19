# simulation results in Wang, Samii, Chang, and Aronow (2024)
rm(list=ls())



###############################################################################
##########################Root Directories#####################################
###############################################################################
root_dir = "C:/Users/haogechang/OneDrive - Microsoft/Desktop/SpatialReplication/SpatialReplication/"
setwd(root_dir)
data_address = "./simulation/data_new/"
###############################################################################
##########################Source Functions#####################################
###############################################################################
source('simulation/scripts/5_functions_generate_figures.R')



# ###############################################################################
# ##########################Generate Figures: Point Interventions################
# ###############################################################################
# 
# 
# #data directory
# case=1
# data_file="simulation_nonmono_10_May_23_points"
# n_ss=5
# sample_sizes <- c(16,36,64,100,144)
# pt_10_nonmono=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)
# 
# 
# #data directory
# case=2
# data_address = "./simulation/data_new/"
# data_file="simulation_interactive_10_May_23_points"
# n_ss=5
# sample_sizes <- c(16,36,64,100,144)
# pt_10_interactive=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)
# 
# 
# #data directory
# case=3
# data_address = "./simulation/data_new/"
# data_file="simulation_nonmono_5_May_23_points"
# n_ss=1
# sample_sizes <- c(729)
# pt_5_nonmono=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)
# 
# #data directory
# case=4
# data_address = "./simulation/data_new/"
# data_file="simulation_interactive_5_May_23_points"
# n_ss=1
# sample_sizes <- c(729)
# pt_5_interactive=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)
# 
# #data directory
# case=1
# data_file="simulation_nonmono_10_May_23_points"
# n_ss=5
# sample_sizes <- c(16,36,64,100,144)
# pt_10_nonmono_sm=generate_figures(data_files,sample_sizes,n_ss,data_address,est='sm',case)
# 
# 
# #data directory
# case=2
# data_address = "./simulation/data_new/"
# data_file="simulation_interactive_10_May_23_points"
# n_ss=5
# sample_sizes <- c(16,36,64,100,144)
# pt_10_interactive_sm=generate_figures(data_files,sample_sizes,n_ss,data_address,est='sm',case)
# 
# 
# #data directory
# case=3
# data_address = "./simulation/data_new/"
# data_file="simulation_nonmono_5_May_23_points"
# n_ss=1
# sample_sizes <- c(729)
# pt_5_nonmono_sm=generate_figures(data_files,sample_sizes,n_ss,data_address,est='sm',case)
# 
# #data directory
# case=4
# data_address = "./simulation/data_new/"
# data_file="simulation_interactive_5_May_23_points"
# n_ss=1
# sample_sizes <- c(729)
# pt_5_interactive_sm=generate_figures(data_files,sample_sizes,n_ss,data_address,est='sm',case)
# 
# 
# ###############################################################################
# ##########################Generate Figures: Polygon Interventions##############
# ###############################################################################
# 
# 
# #data directory
# case=5
# data_file="simulation_nonmono_10_May_23_polygon"
# n_ss=5
# sample_sizes <- c(16,36,64,100,144)
# pl_10_nonmono=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)
# 
# 
# #data directory
# case=6
# data_address = "./simulation/data_new/"
# data_file="simulation_interactive_10_May_23_polygon"
# n_ss=5
# sample_sizes <- c(16,36,64,100,144)
# pl_10_interactive=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)
# 
# 
# #data directory
# case=7
# data_address = "./simulation/data_new/"
# data_file="simulation_nonmono_5_May_23_polygon"
# n_ss=1
# sample_sizes <- c(729)
# pl_5_nonmono=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)
# 
# #data directory
# case=8
# data_address = "./simulation/data_new/"
# data_file="simulation_interactive_5_May_23_polygon"
# n_ss=1
# sample_sizes <- c(729)
# pl_5_interactive=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)

###############################################################################
##############Generate Figures: Point Interventions (Pair) Coverage and Length#
###############################################################################

#data directory
case=9
data_file="Paired_simulation_nonmono_10_Jun_07_points"
n_ss=3
sample_sizes <- c(64,100,144)
pr_pt_10_nonmono=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)

#data directory
case=10
data_file="Paired_simulation_interactive_10_Jun_07_points"
n_ss=3
sample_sizes <- c(64,100,144)
pr_pt_interactive=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)

#data directory
case=9
data_file="Paired_simulation_nonmono_10_Jun_07_points"
n_ss=3
sample_sizes <- c(64,100,144)
pr_pt_10_nonmono_sm=generate_figures(data_files,sample_sizes,n_ss,data_address,est='sm',case)

#data directory
case=10
data_file="Paired_simulation_interactive_10_Jun_07_points"
n_ss=3
sample_sizes <- c(64,100,144)
pr_pt_interactive_sm=generate_figures(data_files,sample_sizes,n_ss,data_address,est='sm',case)


#################################################################################
#########Generate Figures: Polygon Interventions (Pair) Coverage and Length######
#################################################################################

#data directory
case=11
data_file="Paired_simulation_nonmono_10_Jun_13_polygon"
n_ss=3
sample_sizes <- c(64,100,144)
pr_pl_10_nonmono=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)

#data directory
case=12
data_file="Paired_simulation_interactive_10_Jun_13_polygon"
n_ss=3
sample_sizes <- c(64,100,144)
pr_pl_interactive=generate_figures(data_files,sample_sizes,n_ss,data_address,est='regular',case)


#data directory
case=11
data_file="Paired_simulation_nonmono_10_Jun_13_polygon"
n_ss=3
sample_sizes <- c(64,100,144)
pr_pl_10_nonmono_sm=generate_figures(data_files,sample_sizes,n_ss,data_address,est='sm',case)

#data directory
case=12
data_file="Paired_simulation_interactive_10_Jun_13_polygon"
n_ss=3
sample_sizes <- c(64,100,144)
pr_pl_interactive_sm=generate_figures(data_files,sample_sizes,n_ss,data_address,est='sm',case)


###############################################################################
##############Generate Figures: Point Interventions (Pair) MSE ################
###############################################################################


#data directory
case=9
data_file="Paired_simulation_nonmono_10_Jun_07_points"
sample_sizes <- c(64,100,144)
pr_pt_10_nonmono=output_MSE(data_files,sample_sizes,n_ss,data_address,est='regular',case)

#data directory
case=10
data_file="Paired_simulation_interactive_10_Jun_07_points"
sample_sizes <- c(64,100,144)
pr_pt_interactive=output_MSE(data_files,sample_sizes,n_ss,data_address,est='regular',case)

#data directory
case=9
data_file="Paired_simulation_nonmono_10_Jun_07_points"
sample_sizes <- c(64,100,144)
pr_pt_10_nonmono_sm=output_MSE(data_files,sample_sizes,n_ss,data_address,est='sm',case)

#data directory
case=10
data_file="Paired_simulation_interactive_10_Jun_07_points"
sample_sizes <- c(64,100,144)
pr_pt_interactive_sm=output_MSE(data_files,sample_sizes,n_ss,data_address,est='sm',case)


###############################################################################
##############Generate Figures: Polygon Interventions (Pair) MSE###############
###############################################################################

#data directory
case=11
data_file="Paired_simulation_nonmono_10_Jun_13_polygon"
sample_sizes <- c(64,100,144)
pr_pl_10_nonmono=output_MSE(data_files,sample_sizes,n_ss,data_address,est='regular',case)

#data directory
case=12
data_file="Paired_simulation_interactive_10_Jun_13_polygon"
sample_sizes <- c(64,100,144)
pr_pl_interactive=output_MSE(data_files,sample_sizes,n_ss,data_address,est='regular',case)


#data directory
case=11
data_file="Paired_simulation_nonmono_10_Jun_13_polygon"
sample_sizes <- c(64,100,144)
pr_pl_10_nonmono_sm=output_MSE(data_files,sample_sizes,n_ss,data_address,est='sm',case)

#data directory
case=12
data_file="Paired_simulation_interactive_10_Jun_13_polygon"
sample_sizes <- c(64,100,144)
pr_pl_interactive_sm=output_MSE(data_files,sample_sizes,n_ss,data_address,est='sm',case)

