# For HPC
## Model predictions extract and compute model performances

## Import data
dat <- read.csv(paste0(title_of_run,"/data/reef_data_NAs_with_cov_",surveys,".csv")) %>%
  filter(!is.na(COUNT)) %>%
  filter(fGROUP == "HCC") %>%
  arrange(REEF_NAME, SITE_NO, TRANSECT_NO, fYEAR)

# Import data for INLA (to fix later)
dat_inla <- read.csv(paste0(title_of_run,"/data/reef_data_NAs_with_cov_",surveys,".csv")) %>%
 # filter(!is.na(COUNT)) %>%
  filter(fGROUP == "HCC") %>% 
  mutate(across(REEF_NAME, as.factor)) %>%
  rename(Reef = REEF_NAME) %>%
  mutate(tier_cat = case_when(is.na(COUNT) ~ "tier_test", TRUE ~ "tier_train")) 

# Import predictive layer

hexpred <- st_read(paste0(title_of_run,"/data/hexpred_cov.shp"))

hexpred_unique <- hexpred %>%
                group_by(tier) %>%
                filter(row_number() == 1) %>%
                dplyr::select(tier) 
                
# Read model_outputs 

model_list <- list.files(paste0(title_of_run,"/model_outputs/full_run"), recursive = TRUE, pattern = "Rdata")

##################################
############ INLA model(s) 
##################################

model_inla <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "INLA"))

for (k in 1:nrow(model_inla)){
  if(nrow(model_inla)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_inla[k,]))
  
# Making spatial prediction and saving data table 
#obj_inla <- mod$obj_inla
#pred_sum_sf <- predictions_INLA(model.out = mod$fit, n.sim = 1000) 

pred_sum_sf <- predictions_INLA_stacked(model.out = mod$fit, dat = dat_inla, hexpred = hexpred) 
head(pred_sum_sf)

# Making plots and saving them 

# Mean
p_pred <- plot_predictions(pred_sum_sf)

ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred",str_remove(model_inla[k,], ".Rdata"),".png"),
       width=6, height=6)

# Uncertainty 
p_unc <- plot_predictions_unc(pred_sum_sf)

ggsave(p_unc, filename = paste0(title_of_run,"/report/extra/pred_unc",str_remove(model_inla[k,], ".Rdata"),".png"),
       width=6, height=6)

# Save predictions 
save(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_inla[k,], ".Rdata"),".Rdata"),
          row.names = F)
}

##################################
############ FRK model(s)  
##################################

model_frk <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "FRK"))

for (l in 1:nrow(model_frk)){
  if(nrow(model_frk)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_frk[l,]))

# Making spatial prediction and saving data table 
obj_frk <- mod$obj_frk
pred_sum_sf <- predictions_FRK(model.out = mod$fit) 
head(pred_sum_sf)

# Making plots and saving them 

# Mean
p_pred <- plot_predictions(pred_sum_sf)

ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred",str_remove(model_frk[l,], ".Rdata"),".png"),
       width=6, height=6)

# Uncertainty 
p_unc <- plot_predictions_unc(pred_sum_sf)

ggsave(p_unc, filename = paste0(title_of_run,"/report/extra/pred_unc",str_remove(model_frk[l,], ".Rdata"),".png"),
       width=6, height=6)

# Save predictions
save(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_frk[l,], ".Rdata"),".Rdata"),
          row.names = F)

}

##################################
############ BRMS model 
##################################

model_brms <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "brms"))

for (i in 1:nrow(model_brms)){
  if(nrow(model_brms)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_brms[i,]))
  
# Making spatial prediction and saving data table 
pred_sum_sf <-  predictions_brms(model.out = mod$fit, n.sim = 500) 
head(pred_sum_sf)
  
# Making plots and saving them 
  
# Mean
p_pred <- plot_predictions(pred_sum_sf)
  
ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred",str_remove(model_brms[i,], ".Rdata"),".png"),
         width=6, height=6)
  
  # Uncertainty 
p_unc <- plot_predictions_unc(pred_sum_sf)
  
ggsave(p_unc, filename = paste0(title_of_run,"/report/extra/pred_unc",str_remove(model_brms[i,], ".Rdata"),".png"),
         width=6, height=6)

# Save predictions
save(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_brms[i,], ".Rdata"),".Rdata"),
          row.names = F)

}

##################################
############ ML model(s) 
##################################

model_ML <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "ML"))

for (l in 1:nrow(model_ML)){
  if(nrow(model_ML)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_ML[l,]))

# Making spatial prediction and saving data table 
#obj_frk <- mod$obj_frk
pred_sum_sf <- predictions_ML(model.out =  mod$fit$result_trends, dat = dat) 
head(pred_sum_sf)

# Making plots and saving them 

# Mean
p_pred <- plot_predictions(pred_sum_sf)

ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred",str_remove(model_ML[l,], ".Rdata"),".png"),
       width=6, height=6)

# Uncertainty 
p_unc <- plot_predictions_unc(pred_sum_sf)

ggsave(p_unc, filename = paste0(title_of_run,"/report/extra/pred_unc",str_remove(model_ML[l,], ".Rdata"),".png"),
       width=6, height=6)

# Save predictions
save(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_ML[l,], ".Rdata"),".Rdata"),
          row.names = F)

}

###########################################
########################################### BROAD SPATIAL SCALE
###########################################

# Import GRMF layer
GRMF <- st_read(paste0(title_of_run,"/data/true_sp_field.shp")) 

GRMF_all <- st_join(hexpred_unique, GRMF) %>%
    data.frame %>%
    filter(!is.na(Year)) %>%
    ungroup() %>%
    group_by(Year) %>% 
    median_hdci(HCC) %>%
    mutate(fYEAR = Year + min(as.numeric(dat$fYEAR)) - 1) %>%
    mutate(true = plogis(HCC)) #%>%
    #st_as_sf() 
     
##################################
############ INLA model(s) 
##################################

model_inla <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "INLA"))

for (k in 1:nrow(model_inla)){
  if(nrow(model_inla)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_inla[k,]))
  
# Making spatial prediction and saving data table 
obj_inla <- mod$obj_inla
#pred_sum_sf <- predictions_INLA_broad(model.out = mod$fit, n.sim = 500) 

pred_sum_sf <- predictions_INLA_broad_stacked(model.out = mod$fit, dat = dat_inla, hexpred = hexpred)  #index = mod$index, 
head(pred_sum_sf)

# Making plots and saving them 
pred_sum_sf$fYEAR <- as.factor(pred_sum_sf$fYEAR)

# Mean
p_pred <- plot_traj_broad(pred_sum_sf, GRMF_all)
ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred_broad",str_remove(model_inla[k,], ".Rdata"),".png"),
       width=6, height=4)

# Save predictions 
write.csv(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_inla[k,], ".Rdata"),"_broad.csv"),
          row.names = F)
}

##################################
############ FRK model(s)  
##################################

model_frk <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "FRK"))

for (l in 1:nrow(model_frk)){
  if(nrow(model_frk)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_frk[l,]))

# Making spatial prediction and saving data table 
obj_frk <- mod$obj_frk
pred_sum_sf <- predictions_FRK_broad(model.out = mod$fit) %>% data.frame()
head(pred_sum_sf)

# Making plots and saving them 
pred_sum_sf$fYEAR <- as.factor(pred_sum_sf$fYEAR)

p_pred <- plot_traj_broad(pred_sum_sf, GRMF_all)
ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred_broad",str_remove(model_frk[l,], ".Rdata"),".png"),
       width=6, height=4)

# Save predictions 
write.csv(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_frk[l,], ".Rdata"),"_broad.csv"),
          row.names = F)

}

##################################
############ BRMS model(s) 
##################################

model_brms <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "brms"))

for (i in 1:nrow(model_brms)){
  if(nrow(model_brms)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_brms[i,]))
  
# Making spatial prediction and saving data table 
pred_sum_sf <-  predictions_brms_broad(model.out = mod$fit, n.sim = 500) 
head(pred_sum_sf)
  
# Making plots and saving them 

p_pred <- plot_traj_broad(pred_sum_sf, GRMF_all)
ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred_broad",str_remove(model_brms[k,], ".Rdata"),".png"),
       width=6, height=4)

# Save predictions 
write.csv(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_brms[k,], ".Rdata"),"_broad.csv"),
          row.names = F)

}

##################################
############ ML model(s) 
##################################

model_ML <- model_list %>% 
  data.frame() %>%
  filter(str_detect(., "ML"))

for (w in 1:nrow(model_ML)){
  if(nrow(model_ML)==0) next

load(paste0(title_of_run,"/model_outputs/full_run/", model_ML[w,]))

# Making spatial prediction and saving data table 
pred_sum_sf <-  predictions_ML_broad(model.out = mod$fit$result_trends, dat = dat) 
head(pred_sum_sf)
  
# Making plots and saving them 

p_pred <- plot_traj_broad(pred_sum_sf, GRMF_all)
ggsave(p_pred, filename = paste0(title_of_run,"/report/extra/pred_broad",str_remove(model_ML[w,], ".Rdata"),".png"),
       width=6, height=4)

# Save predictions 
write.csv(pred_sum_sf, file = paste0(title_of_run,"/model_outputs/predictions/",str_remove(model_ML[w,], ".Rdata"),"_broad.csv"),
          row.names = F)

}

### Effect of disturbances

# Full table 
# coef_table_all <- coef_uncertainty(M, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = TRUE)%>%
#   data.frame() %>%
#   tibble::rownames_to_column() %>%
#   tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value")%>%
#   mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
#   tidyr::pivot_wider(names_from = rowname, values_from = value)
# 
# # Fixed effects only 
# 
# coef_table_fixed <- coef_table_all %>%
#   filter(Type == "fixed") 
# 
# p_coef <- ggplot(coef_table_fixed[-1,], aes(y=param, x=`50%`))+ geom_point() +
#   geom_errorbar(aes(y = param, xmin = `2.5%`, xmax = `97.5%`), width=.1) + 
#   geom_vline(xintercept = 0, linetype = "dashed") +theme_bw() +
#   xlab("Effect size") + ylab("")
# 
# ggsave(plot =  p_coef,width=6, height=4, file = "../figures/fixed_effects_model.png")
# 

##################################
############ BRMS model 
##################################
# brms models to be read in the same way - TO DO 

#model_brms <- model_list %>% 
#  data.frame() %>%
#  filter(str_detect(., "brms"))

#for ( i in 1:nrow(model_brms)){
 # mod <- load(paste0("../model_outputs/", model_brms[i,]))
  
  
# sdmTMB

# 
# model_sdmTMB <- model_list %>% 
#   data.frame() %>%
#   filter(str_detect(., "sdmTMB"))
# 
# mod <- load(paste0("../model_outputs/", model_sdmTMB[1,]))
# 
# 
# ################## Fit on new data
# 
# hexpred_dat <- hexpred %>%
#   st_centroid() %>%
#   mutate(Longitude = st_coordinates(.)[,1],
#          Latitude = st_coordinates(.)[,2]) %>%
#   dplyr::select(fYEAR, Longitude, Latitude, Reef) %>%
#   st_drop_geometry()#%>%
# #  mutate(Reef = factor(NA))
# 
# colnames(hexpred_dat) <- c("fYEAR", "LONGITUDE", "LATITUDE","Reef")
# hexpred_dat$fYEAR <- as.numeric(as.character(hexpred_dat$fYEAR))
# 
# M_new_sf <- predict(fit_mle, newdata = hexpred_dat, type = 'link', re_form_iid = ~ 0) %>%
#   mutate(pred = plogis(est))
# 
# # Prediction
# pal_pred <- lacroix_palette("Pamplemousse", n = 10, type = "continuous")
# pal_unc<- wes_palette("Zissou1", 10, type = "continuous")
# 
# fit_plot <- M_new_sf %>%
#   ggplot() +
#   geom_point(aes(y=LATITUDE, x=LONGITUDE, color=mean), size=1) +
#   facet_wrap(~fYEAR)+
#   scale_color_gradientn(colours = pal_pred, limits = c(0,1)) +
#   coord_fixed() +
#   theme_bw() +
#   theme(axis.title = element_blank()) +
#   theme(legend.position="right", legend.justification=c(1,1))
# 
# #ggsave(filename = paste0(model_name,"/Model_prediction.pdf"), plot = fit_plot, width = 14, height = 10)
# 
# M_new_sf %>%
#   ggplot() +
#   geom_point(aes(y=LATITUDE, x=LONGITUDE, color=epsilon_st), size=1) +
#   facet_wrap(~fYEAR)+
#   scale_color_gradientn(colours = pal_pred, limits = c(0,1)) +
#   coord_fixed() +
#   theme_bw() +
#   theme(axis.title = element_blank()) +
#   theme(legend.position="right", legend.justification=c(1,1))+
#   ggtitle("Spatiotemporal random effects only") 
# 
