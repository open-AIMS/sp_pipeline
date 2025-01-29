# For HPC
## Pipeline sp_models V2
rm(list=ls())

# Indicate the working directory
setwd(paste0(getwd(),"/SP_models"))

# Loading packages and functions
source("R/packages.R")
source("R/functions.R")

# Creating folder structure associated with the name of the simulation

title_of_run <- "Low_location"
make_folders(title_of_run)

# Type of monitoring surveys
surveys <- "fixed" # or "random" (to update codes for random survey)

##################################
############ Create synthetic data
##################################

## Config of the spatio-temporal model
config_sp <- list(
  seed = 1,
  crs = 4326,
  model = "Exp",
  psill = 1,
  range = 15,
  nugget = 0,
  alpha = 2,
  kappa = 1,
  variance = 1,
  patch_threshold = 1.75,
  reef_width = 0.01,
  years = 1:15,
  dhw_weight = 0.8,
  cyc_weight = 0.4,
  other_weight = 0.1,
  hcc_growth = 0.3,
  sc_growth =  0.3
)

## Config of sampling design for large scale details
config_lrge <- list(n_locs = 8, n_sites = 2, seed = 123)

## Config for sampling detaisl for fine scale detais 
config_fine <- list(
  years =  1:15,
  Number_of_transects_per_site = 5,
  Depths = 1,
  Number_of_frames_per_transect = 100,
  Points_per_frame = 5,
  ## Note, the following are on the link scale
  hcc_site_sigma = 0.5, # variability in Sites within Locations
  hcc_transect_sigma = 0.2, # variability in Transects within Sites
  hcc_sigma = 0.1, # random noise

  sc_site_sigma = 0.05, # variability in Sites within Locations
  sc_transect_sigma = 0.02, # variability in Transects within Sites
  sc_sigma = 0.01, # random noise

  ma_site_sigma = 0.5, # variability in Sites within Locations
  ma_transect_sigma = 0.2, # variability in Transects within Sites
  ma_sigma = 0.1 # random noise
)

## Generate point-based data 
config_pt <- list(
  Depths = 2,
  Depth_effect_multiplier = 2,
  Number_of_transects_per_site = 5,
  Number_of_frames_per_transect = 100,
  Points_per_frame = 5
)

config_list <- list(config_sp, config_lrge, config_fine, config_pt)
save(config_list, file = paste0(title_of_run,"/lists_of_parameters.RData"))


source(paste0("scripts/make_data_",surveys,"_synthos_custom.R"))

##################################
############ Create grid
##################################

# Choose the grid size in degree
grid_size <- .1

source("scripts/make_grid.R")

##################################
############ Create predictive layer
##################################

source("scripts/make_predictive_layer.R")

##################################
############ INLA model
##################################

# Without covariates 
#source("scripts/INLA_model_without.R")

# With covariates 
source("scripts/INLA_model.R")

##################################
############ INLA VARITIONAL BAYES model
##################################

# Without covariates
#source("scripts/INLA_VB_model_without.R")

# With covariates
#source("scripts/INLA_VB_model.R")

##################################
############ FRK model
##################################

# Without covariates
#source("scripts/FRK_model_without.R")

# With covariates
source("scripts/FRK_model.R")

##################################
############ BRMS model
##################################

#niter <- 10000
#nwarm <- niter / 2

# Without covariates
#source("scripts/brms_model_without.R")

# With covariates
#source("scripts/brms_model.R")

##################################
############ ML model(s)  
##################################

# Config ----

options(future.globals.maxSize = 100000*1024^2) # 100 Gb
plan(strategy = multisession, workers = 2)

# Bootstrap runs number ----
#n_bootstrap <- 2 #(for dev)
n_bootstrap <- 1200

# Folding for tuning 
#vfold <- 2 #(for dev) 
vfold <- 5

source("scripts/ML_model_no.R")
#source("scripts/ML_model_gp.R")
#source("scripts/ML_model_bs.R")

##################################
############ model predictions (incl. regional level)
##################################

source("scripts/model_predictions.R")

##################################
############ model predictive checks 
##################################

source("scripts/model_checks.R")

################################################## OLD 

##################################
############ leave-out data analysis 
##################################

#source("scripts/leave_out_data_run.R")

##################################
############ model diagnostics 
##################################

#source("scripts/model_performances.R")

##################################
############ rendering the report 
##################################

#source("scripts/report.qmd")
