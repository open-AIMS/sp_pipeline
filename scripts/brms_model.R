# For HPC

## Import data
dat <- read.csv(paste0(title_of_run,"/data/reef_data_NAs_with_cov_",surveys,".csv")) %>%
  filter(!is.na(COUNT)) %>%
  filter(fGROUP == "HCC")%>% 
  mutate(across(REEF_NAME, as.factor)) %>%
  rename(Reef = REEF_NAME)

dat$fYEAR <- as.numeric(as.character(dat$fYEAR))

## Fit brms model

start_time <- Sys.time()
fit <- brm(COUNT | trials(TOTAL) ~ 1 + DHW + CYC + OT + gp(LONGITUDE, LATITUDE, by = fYEAR, gr = TRUE) + # not working when removing "by = fYEAR"
          (1 | Reef),
           data = dat, family = binomial, iter = niter, 
           warmup = nwarm, thin = 5, chains = 3,
           verbose = F,
           seed = 100,
           refresh = max(niter / 100, 1))

end_time <- Sys.time()
computing_time <- end_time - start_time

## Save brms model 

mod <- list(fit = fit, computing_time = computing_time)
save(mod, file = paste0(title_of_run,"/model_outputs/full_run/brms.Rdata"))