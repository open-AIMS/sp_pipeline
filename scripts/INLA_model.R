## INLA model for HPC
formula <- "cov"

# Import data
dat <- read.csv(paste0(title_of_run,"/data/reef_data_NAs_with_cov_",surveys,".csv")) %>%
 # filter(!is.na(COUNT)) %>%
  filter(fGROUP == "HCC") %>% 
  mutate(across(REEF_NAME, as.factor)) %>%
  rename(Reef = REEF_NAME) %>%
  mutate(tier_cat = case_when(is.na(COUNT) ~ "tier_test", TRUE ~ "tier_train")) 

# Import predictive layer
hexpred <- st_read(paste0(title_of_run,"/data/hexpred_cov.shp"))

# Prepare objects for the model 
obj_inla <- inla_prep(dat, hexpred) 

start_time <- Sys.time()
#fit0 <-  inla(obj_inla$form,data = inla.stack.data(obj_inla$stack.est),
#                  family= 'binomial',
#                  Ntrials=Total,
#                  control.predictor = list(compute = TRUE,
#                                           link = 1,
#                                           A = inla.stack.A(obj_inla$stack.est)
#                  ),
#                  control.compute = list(return.marginals.predictor=TRUE,config = TRUE, dic= TRUE,
#                                         return.marginals = TRUE), 
#                  verbose = FALSE)

fit <-  inla(obj_inla$form,data = inla.stack.data(obj_inla$stack.full),
             family= 'binomial',
             Ntrials=Total,
             control.predictor = list(compute = TRUE,
                                      link = 1,
                                      A = inla.stack.A(obj_inla$stack.full)
             ),
             control.compute = list(return.marginals.predictor=TRUE,config = TRUE, dic= TRUE,
                                    return.marginals = TRUE), 
             verbose = FALSE)

#index <- inla.stack.index(stack = obj_inla$stack.full, tag = "pred")$data

end_time <- Sys.time()
computing_time <- end_time - start_time

#mod <- list(fit = fit0, computing_time = computing_time, obj_inla = obj_inla)
mod <- list(fit = fit, computing_time = computing_time, obj_inla = obj_inla) #, index = index
save(mod, file = paste0(title_of_run,"/model_outputs/full_run/INLA.Rdata"))
