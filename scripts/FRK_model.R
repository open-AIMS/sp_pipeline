# For HPC
## FRK model 

## Import data
dat <- read.csv(paste0(title_of_run,"/data/reef_data_NAs_with_cov_",surveys,".csv")) %>%
  filter(!is.na(COUNT)) %>%
  filter(fGROUP == "HCC") %>%
  arrange(REEF_NAME, SITE_NO, TRANSECT_NO, fYEAR)

## Import predictive layer
hexpred <- st_read(paste0(title_of_run,"/data/hexpred_cov.shp")) 

## Prepare objects for the model 
obj_frk <- frk_prep(dat) 

## Fit FRK model

start_time <- Sys.time()
fit <- FRK(f = COUNT ~ 1 + DHW + CYC + OT + (1|Reef), 
         data = list(obj_frk$STObj), 
         BAUs = obj_frk$ST_BAUs, 
         basis = obj_frk$basis, 
         response = "binomial", 
         link = "logit", 
         K_type = "precision", 
         method = "TMB", 
         est_error = FALSE)

end_time <- Sys.time()
computing_time <- end_time - start_time

## Save FRK model 

mod <- list(fit = fit, computing_time = computing_time, obj_frk = obj_frk)
save(mod, file = paste0(title_of_run,"/model_outputs/full_run/FRK.Rdata"))


