# Year : factor
# Space : Reef

## Created by: Hina
## Modified from: Jeremy
## Date: November 2024 
## Contact: hinagluza@icloud.com

# Model name ----
model_name = "ML_original"

# Load data at known locations used for model tuning ----
dat_raw <- read.csv(paste0(title_of_run,"/data/reef_data_NAs_with_cov_fixed.csv")) |> 
    filter(fGROUP == "HCC") |> filter(!is.na(COUNT)) |>
    dplyr::select(- c(LONGITUDE,LATITUDE)) # use lng/lat of tier centroid instead

# Add coordinates associated with tier-level
hexpred <- st_read(paste0(title_of_run,"/data/hexpred_cov.shp"))

hexpred_unique <- hexpred %>%
  group_by(tier) %>% 
  filter(row_number() == 1) %>%  
  ungroup() %>%                
  st_centroid() %>%            
  mutate(
    LONGITUDE = st_coordinates(geometry)[, 1],  
    LATITUDE = st_coordinates(geometry)[, 2]   
  ) %>%
  dplyr::select(tier, Reef, LONGITUDE, LATITUDE) |>
  st_drop_geometry()

dat_cov <- dat_raw |>
  left_join(hexpred_unique) |> 
  dplyr::select(-c(P_CODE, fDEPTH, TRUE_COUNT, TOTAL, SITE_NO, TRANSECT_NO, tier, fGROUP, REEF_NAME
                   )) |> 
  mutate(across(c(Reef, fYEAR), as_factor)) |> #LONGITUDE, LATITUDE, 
  mutate(across(c(CYC, DHW, OT), ~ c(scale(.)))) 

# Load the data at unknown locations used for model bootstrapping using predictive layer  
hexpred_dat <- hexpred %>%
    st_centroid() %>%
    mutate(LONGITUDE = st_coordinates(.)[,1],
           LATITUDE = st_coordinates(.)[,2]) %>%
    dplyr::select(fYEAR, LONGITUDE, LATITUDE, DHW, CYC, OT, Reef, tier) %>%
    st_drop_geometry()%>% 
    mutate(TOTAL = round(mean(dat_raw$TOTAL),0)) %>%
  dplyr::select(-c(tier, TOTAL
                   )) |> 
  mutate(across(c(Reef, fYEAR), as_factor)) |> # LONGITUDE, LATITUDE, 
  mutate(across(c(CYC, DHW, OT), ~ c(scale(.))))

# Model tunning 
model_tuning <- ml_tuning(dat_cov, vfold)

# Save tuning parameters 
saveRDS(model_tuning, file = paste0(title_of_run,"/model_outputs/tuning/",model_name,"_tunning.Rdata"))

# Bootstrapping 
model_hyperparams <- model_tuning$model_fit$model_hyperparams
  
# bootstrap model
start_time <- Sys.time()
bootstrap_model <- map(1:n_bootstrap,~model_bootstrap_run(bootstrap_i = .,
                                             pdp = FALSE)) %>% 
  map_df(., ~ as.data.frame(map(.x, ~ unname(nest(.))))) %>% 
  map(., bind_rows)

names(bootstrap_model) <- c("model_performance", "result_vip", "result_trends")

end_time <- Sys.time()
computing_time <- end_time - start_time
  
mod <- list(fit = bootstrap_model, computing_time = computing_time)
save(mod, file = paste0(title_of_run,"/model_outputs/full_run/", model_name,".Rdata"))







