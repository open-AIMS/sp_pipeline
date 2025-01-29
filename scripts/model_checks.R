# Model predictions extract and compare with true values 

# Extract computing time for each model 
model_list_run <- list.files(paste0(title_of_run,"/model_outputs/full_run"), recursive = TRUE, pattern = "Rdata") %>% 
  data.frame() 

computing_time_list <- list()

for (k in 1:nrow(model_list_run)){
load(paste0(title_of_run,"/model_outputs/full_run/", model_list_run[k,]))
computing_time_list[[k]] <- cbind(round(as.numeric(mod$computing_time, units = "mins"),3), str_remove(model_list_run[k,], ".Rdata"))
}

computing_time <- data.frame(do.call(rbind,computing_time_list))
colnames(computing_time) <- c("time", "Model")

##################################
############ Leave_out tier analysis using true count  
################################## SPATIAL 

#Import data
dat <- read.csv(paste0(title_of_run,"/data/reef_data_NAs_with_cov_",surveys,".csv")) %>%
  mutate(tier_cat = case_when(is.na(COUNT) ~ "tier_test", TRUE ~ "tier_train")) 

# Group_by tier and year 
dat_grp <- dat %>%
  filter(!is.na(LONGITUDE)) %>% # needed because no values of TRUE_COUNT in every tier 
  #filter(fGROUP == "HCC") %>% 
  group_by(tier, fYEAR) %>%
  mutate(TRUE_COUNT_tier = mean(TRUE_COUNT),
          TOTAL_tier = mean(TOTAL)) %>%
  filter(row_number() == 1) %>%
  dplyr::select(tier, fYEAR, TRUE_COUNT_tier, TOTAL_tier, tier_cat) %>%
  mutate(COVER_tier = TRUE_COUNT_tier / TOTAL_tier) %>%
  data.frame() 

dat_grp$fYEAR <- as.character(dat_grp$fYEAR) 
dat_grp$tier <- as.character(dat_grp$tier) 

# Read model_outputs from model fit 

model_list_pred <- list.files(paste0(title_of_run,"/model_outputs/predictions"), recursive = TRUE, pattern = "Rdata") %>% 
  data.frame() 

# List to save outputs 
pred_dat_true.list <- list()
list.indicators_tier <- list()

for (j in 1:nrow(model_list_pred)){ 
  
load(paste0(title_of_run,"/model_outputs/predictions/", model_list_pred[j,]))

pred_dat <- pred_sum_sf %>%
  st_drop_geometry() %>%
  data.frame() 

pred_dat$fYEAR <- as.character(pred_dat$fYEAR) 
pred_dat$tier <- as.character(pred_dat$tier) 

 # Join with true values 

pred_dat_true.list[[j]] <- pred_dat %>%
 inner_join(dat_grp) %>%
 mutate(Model = str_remove(model_list_pred[j,], ".Rdata")) %>%
 mutate(Diff = COVER_tier - pred)

 # Measures

test_pred_sum <- pred_dat_true.list[[j]] %>%
 mutate(true = TRUE_COUNT_tier / TOTAL_tier) %>%
 filter(tier_cat == "tier_test")

list.indicators_tier[[j]] <- c(round(coverage95(test_pred_sum$true, test_pred_sum$.lower, test_pred_sum$.upper),2),
                              round(IS95(test_pred_sum$true, test_pred_sum$.lower, test_pred_sum$.upper),2),
                              round(RMSPE(test_pred_sum$true, test_pred_sum$pred),2),
                              round(crps(test_pred_sum$true, data.frame(test_pred_sum$pred, test_pred_sum$Unc))$CRPS,2), 
                              str_remove(model_list_pred[j,], ".Rdata"), title_of_run)


}

pred_dat_true <- data.frame(do.call(rbind,pred_dat_true.list)) %>%
                 mutate(Title_of_run = title_of_run)

# Vizualisation 
pal_model <- lacroix_palette("Pamplemousse", n = length(unique(pred_dat_true$Model)), type = "discrete")

p_check <- ggplot(pred_dat_true %>% filter(tier_cat == "tier_test"), aes(x = COVER_tier, y = pred)) + 
 geom_point(alpha=0.4) + 
 facet_wrap(~Model, nrow = 1) +
 geom_abline(linetype = "dashed") +
 labs(title = unique(pred_dat_true$Title_of_run)) +
 xlab("True values") +
 ylab("Predicted values") +
 theme_bw() +
 scale_color_manual(values = pal_model)

ggsave(p_check, filename = paste0(title_of_run,"/report/extra/check_true.png"),
       width=8, height=6)

# Map of differences 

hexpred_unique <- st_read(paste0(title_of_run,"/data/hexpred_cov.shp")) %>%
group_by(tier) %>%
filter(row_number() == 1) %>%
dplyr::select(tier) 

hexpred_unique$tier <- as.character(hexpred_unique$tier) 

pred_dat_true_sf <- pred_dat_true %>%
left_join(hexpred_unique, by = "tier") %>%
st_as_sf() 

# Vizualisation 

#pred_dat_true_sf %>%
#group_split(Model) %>%
#map(~ { plot_diff(.x) -> p_plot 
#        ggsave(p_plot, filename = paste0(title_of_run,"/report/extra/pred_diff",first(.x$Model),".png"),
#       width=6, height=6)
#       })

# Get performances measures 
indicators_table_tier <- data.frame(do.call(rbind,list.indicators_tier)) 
colnames(indicators_table_tier) <- c("Cvg", "IS", "RMSPE", "CRPS",
                                     "Model", "Title_of_run")

indicators_table_tier  <- indicators_table_tier %>%  left_join(computing_time) %>%
                          mutate(across(c(Cvg, IS, RMSPE, CRPS, time), as.numeric))

write.csv(indicators_table_tier, file = paste0(title_of_run,"/model_outputs/leave_out/table_performances_tier_true.csv"), row.names = F)

# Viz model performances
scale <- 4.5 # scale is an adjustment for a 4k screen

p_perf <- ggspider(indicators_table_tier %>% dplyr::select(Model, RMSPE, `Cvg`, `IS`, `CRPS`, time), axis_name_offset = 0.15,
               background_color = "gray98", fill_opacity = 0.15, polygon = FALSE) +
  labs(col = "Model", title = unique(indicators_table_tier$Title_of_run)) +
  theme(legend.position = "bottom", 
        plot.title = element_text(size = 5*scale, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 4*scale, hjust = 0.5),
        plot.caption = element_text(size = 3*scale),
        legend.text = element_text(size = 2*scale),
        legend.title = element_text(size = 2*scale, face = "bold"))

ggsave(p_perf, filename = paste0(title_of_run,"/report/extra/viz_perf.png"),
       width=5.5, height=6)


##################################
############ Leave_out tier analysis using true count  
################################## TEMPORAL TRENDS 

# Coral cover trajectories of tier with data 
dat_plot <- dat %>%
filter(fGROUP == "HCC") %>%
filter(!is.na(COUNT)) 

dat_plot$fYEAR <- as.character(dat_plot$fYEAR)

# Pick 9 tiers randomly 
tier.observed <- dat_plot  %>%
  dplyr::select(tier) %>%
  distinct() %>%
  sample_n(size = 9) %>%
  pull(tier)

# Vizualisation - not working
#pred_dat_true_sf %>%
#filter(tier_cat == "tier_train") %>%
#filter(tier %in% tier.observed) %>%
#group_split(Model) %>%
#map(~ { plot_traj(dat_plot, .x) -> p_plot 
#       ggsave(p_plot, filename = paste0(title_of_run,"/report/extra/pred_traj_tier_data",first(.x$Model),".png"),
#       width=6 height=8)
#       })

dat_plot <- dat_plot %>%
  filter(tier %in% tier.observed) 

dat_pred <- pred_dat_true_sf %>%
# filter(tier_cat == "tier_train") %>%
  filter(tier %in% tier.observed) %>%
  group_split(Model)

for ( i in 1:length(dat_pred)){
 p_traj <- ggplot() + 
  geom_line(data = dat_plot, aes(x = fYEAR, y = (COUNT/TOTAL), group = interaction(as.factor(TRANSECT_NO), REEF_NAME)), 
            show.legend = FALSE, linewidth=.1, col="grey30") + 
  geom_ribbon(data = dat_pred[[i]], aes(x=fYEAR,ymin=.lower, ymax=.upper, group=1),alpha=.2, fill ="#0072B2FF") +
  geom_line(data = dat_pred[[i]], aes(x=fYEAR, y=pred, group=1),linewidth=.4) +
  facet_wrap(~tier, ncol=3) +
  ylab("Cover") + xlab("Year")+theme_bw()+
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
  scale_x_discrete(breaks= nth_element(unique(dat_pred[[i]]$fYEAR),1,4))

ggsave(p_traj, filename = paste0(title_of_run,"/report/extra/pred_traj_tier_data", unique(dat_pred[[i]]$Model),".png"),
       width=6, height=8)
}

# Coral cover trajectories of tier without data 
# Pick 9 tiers randomly 
dat_plot_no <- dat %>%
filter(fGROUP == "HCC") %>%
filter(is.na(COUNT)) 

dat_plot_no$fYEAR <- as.character(dat_plot_no$fYEAR)

tier.without <- dat_plot_no %>%
  dplyr::select(tier) %>%
  distinct() %>%
  sample_n(size = 9) %>%
  pull(tier)

dat_plot_no <- dat_plot_no %>%
  filter(tier %in% tier.without) 

dat_pred_without <- pred_dat_true_sf %>%
# filter(tier_cat == "tier_test") %>%
  filter(tier %in% tier.without) %>%
  group_split(Model)

for ( i in 1:length(dat_pred_without)){
 p_traj <- ggplot() + 
   geom_line(data = dat_plot_no, aes(x = fYEAR, y = (TRUE_COUNT/TOTAL)*100, group = interaction(as.factor(TRANSECT_NO), REEF_NAME)), 
            show.legend = FALSE, linewidth=.1, col="red", linetype = "dashed") +
  geom_ribbon(data = dat_pred_without[[i]], aes(x=fYEAR,ymin=.lower*100, ymax=.upper*100, group=1),alpha=.2, fill ="#D68A8A") +
  geom_line(data = dat_pred_without[[i]], aes(x=fYEAR, y=pred*100, group=1), linewidth=.4) +
  facet_wrap(~tier, ncol=3) +
  ylab("Cover") + xlab("Year")+theme_bw()+
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
  scale_x_discrete(breaks= nth_element(unique(dat_pred_without[[i]]$fYEAR),1,4))

ggsave(p_traj, filename = paste0(title_of_run,"/report/extra/pred_traj_tier_no_data", unique(dat_pred_without[[i]]$Model),".png"),
       width=6, height=8)
}

##################################
############ Leave_out tier analysis using Gaussian Random Markov Field (GRMF) values  
################################## SPATIAL 

# # Import GRMF layer
# GRMF <- st_read(paste0(title_of_run,"/data/true_sp_field.shp"))

# # Plot 
# p_g <- ggplot() + geom_sf(data = GRMF, aes(fill = plogis(HCC))) +
#         scale_color_gradientn(colors = terrain.colors(10))

# ggsave(p_g, filename = paste0(title_of_run,"/report/extra/true_sp_field.png"),
#        width=6, height=6)

# # Combine 

# GRMF_tier <- st_join(hexpred_unique, GRMF) %>%
#     group_by(Year,tier) %>% 
#     median_hdci(HCC) %>%
#     mutate(fYEAR = Year + 2009) %>%
#     st_as_sf()
     
# # Plot 
# p_g_tier <- ggplot() + geom_sf(data = GRMF_tier, aes(fill = plogis(HCC))) +
#         scale_color_gradientn(colors = terrain.colors(10)) +
#         facet_wrap(~fYEAR)

# ggsave(p_g_tier, filename = paste0(title_of_run,"/report/extra/true_sp_field_tier.png"),
#        width=6, height=6)

# # Plot together 
# #missing_tier <- GRMF_tier %>% group_by(tier) %>% tally() %>% filter(n<12)

# #p_g_both <- ggplot() + geom_sf(data = GRMF_tier %>% filter(tier %in% missing_tier$tier), fill = "transparent") +
# #        geom_sf(data = GRMF, aes(fill = plogis(HCC))) +
# #        scale_color_gradientn(colors = terrain.colors(10)) +
# #        facet_wrap(~fYEAR)

# #ggsave(p_g_both, filename = paste0(title_of_run,"/report/extra/true_sp_field_missing_tier.png"),
# #       width=6, height=6)

# GRMF_tier <- GRMF_tier %>%
#     dplyr::select(fYEAR, tier, HCC) %>%
#     rename(true = HCC) 

# GRMF_tier$fYEAR <- as.character(GRMF_tier$fYEAR) 
# GRMF_tier$tier <- as.character(GRMF_tier$tier) 

# # add tier information 
# dat_grp_tier <- dat_grp %>%
# group_by(tier) %>%
# filter(row_number() == 1) %>%
# dplyr::select(tier,tier_cat) 

# GRMF_tier <- GRMF_tier %>%
#     left_join(dat_grp_tier) %>%
#     mutate(tier_cat = case_when(is.na(tier_cat) ~ "tier_test", TRUE ~ "tier_train")) %>%
#     data.frame()

# # List to save outputs 
# pred_GRMF_true.list <- list()
# list.indicators_GRMF <- list()

# for (j in 1:nrow(model_list_pred)){ 
  
# load(paste0(title_of_run,"/model_outputs/predictions/", model_list_pred[j,]))

# pred_dat <- pred_sum_sf %>%
#   st_drop_geometry() %>%
#   data.frame() 

# pred_dat$fYEAR <- as.character(pred_dat$fYEAR) 
# pred_dat$tier <- as.character(pred_dat$tier) 

#  # Join with true values 

# pred_GRMF_true.list[[j]] <- pred_dat %>%
#  left_join(GRMF_tier) %>%
#  mutate(true = plogis(true)) %>%
#  mutate(Model = str_remove(model_list_pred[j,], ".Rdata")) %>%
#  filter(!is.na(true)) 

# # Measures

# test_pred_sum <- pred_GRMF_true.list[[j]] %>%
# mutate(Unc = .upper - .lower) %>%
#  filter(tier_cat == "tier_test")

# list.indicators_GRMF[[j]]<- c(round(coverage95(test_pred_sum$true, test_pred_sum$.lower, test_pred_sum$.upper),2),
#                               round(IS95(test_pred_sum$true, test_pred_sum$.lower, test_pred_sum$.upper),2),
#                               round(RMSPE(test_pred_sum$true, test_pred_sum$pred),2),
#                               round(crps(test_pred_sum$true, data.frame(test_pred_sum$pred, test_pred_sum$Unc))$CRPS,2), 
#                               str_remove(model_list_pred[j,], ".Rdata"))

# }

# pred_GRMF_true <- data.frame(do.call(rbind,pred_GRMF_true.list)) %>%
#  mutate(Diff = true - pred)

# # Vizualisation 
# pal_model2 <- lacroix_palette("Pamplemousse", n = length(unique(pred_GRMF_true$Model)), type = "discrete")

# p_check2 <- ggplot(pred_GRMF_true, aes(x = true, y = pred, col = Model)) + 
#  geom_point() + 
#  geom_abline(linetype = "dashed") +
#  xlab("GRMF values") +
#  ylab("Predictions") +
#  scale_color_manual(values = pal_model2) +
#  facet_wrap(Model~tier_cat, nrow = 2)

# ggsave(p_check2, filename = paste0(title_of_run,"/report/extra/check_GRMF.png"),
#        width=6, height=6)

# # Vizualisaiton 
# pred_GRMF_true_sf <- pred_GRMF_true %>%
# left_join(hexpred_unique, by = "tier") %>%
# st_as_sf() 

# pred_GRMF_true_sf %>%
# group_split(Model) %>%
# map(~ { plot_diff(.x) -> p_plot 
#         ggsave(p_plot, filename = paste0(title_of_run,"/report/extra/pred_diff_GRMF",first(.x$Model),".png"),
#        width=6, height=6)
#        })

# # Get performances measures 

# indicators_table_GRMF <- data.frame(do.call(rbind,list.indicators_GRMF))
# colnames(indicators_table_GRMF) <- c("95% coverage", "95% interval score", "Root-mean-squared prediction error", "Continuous ranked probability score",
#                                      "Model")

# write.csv(indicators_table_GRMF %>%  left_join(computing_time), file = paste0(title_of_run,"/model_outputs/leave_out/table_performances_tier_GRMF.csv"), row.names = F)

# ##################################
# ############ Leave_out tier analysis using Gaussian Random Markov Field (GRMF) values  
# ################################## TEMPORAL TRENDS 

# dat_pred_GRMF <- pred_GRMF_true_sf %>%
#  filter(tier_cat == "tier_train") %>%
#  filter(tier %in% tier.observed) %>%
#  group_split(Model)

# for ( i in 1:length(dat_pred_GRMF)){
#  p_traj <- ggplot() + 
#   geom_line(data = dat_plot, aes(x = fYEAR, y = (COUNT/TOTAL)*100, group = interaction(as.factor(TRANSECT_NO), REEF_NAME)), 
#             show.legend = FALSE, linewidth=.1, col="grey30") + 
#   geom_ribbon(data = dat_pred_GRMF[[i]], aes(x=fYEAR,ymin=.lower*100, ymax=.upper*100, group=1),alpha=.2, fill ="#0072B2FF") +
#   geom_line(data = dat_pred_GRMF[[i]], aes(x=fYEAR, y=pred*100, group=1),size=.4) +
#     geom_line(data = dat_pred_GRMF[[i]], aes(x=fYEAR, y=true*100, group=1),size=.4, col = "blue", linetype = "dashed") +
#   facet_wrap(~tier, ncol=3) +
#   ylab("Cover") + xlab("Year")+theme_bw()+
#   theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
#         axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
#         axis.title.x=element_text(size=11),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill = 'white')) +
#   scale_x_discrete(breaks= nth_element(unique(dat_pred_GRMF[[i]]$fYEAR),1,4))

# ggsave(p_traj, filename = paste0(title_of_run,"/report/extra/pred_traj_tier_data_GRMF", unique(dat_pred[[i]]$Model),".png"),
#        width=6, height=8)
# }

# # Coral cover trajectories of tier without data 
# # Pick 9 tiers randomly 

# dat_pred_GRMF_without <- pred_GRMF_true_sf  %>%
#  filter(tier_cat == "tier_test") 
 
#  tier.without.GRMF <- dat_pred_GRMF_without %>%
#   dplyr::select(tier) %>%
#   distinct() %>%
#   sample_n(size = 9) %>%
#   pull(tier)


#  dat_pred_GRMF_without <- dat_pred_GRMF_without %>%
#   filter(tier %in% tier.without.GRMF) %>%
#   group_split(Model)

# for ( i in 1:length(dat_pred_GRMF_without)){
#  p_traj <- ggplot() + 
#   # geom_line(data = dat_plot_no, aes(x = fYEAR, y = (TRUE_COUNT/TOTAL)*100, group = interaction(as.factor(TRANSECT_NO), REEF_NAME)), 
#   #          show.legend = FALSE, linewidth=.1, col="red", linetype = "dashed") +
#   geom_ribbon(data = dat_pred_GRMF_without[[i]], aes(x=fYEAR,ymin=.lower*100, ymax=.upper*100, group=1),alpha=.2, fill ="#D68A8A") +
#   geom_line(data = dat_pred_GRMF_without[[i]], aes(x=fYEAR, y=pred*100, group=1),size=.4) +
#   geom_line(data = dat_pred_GRMF_without[[i]], aes(x=fYEAR, y=true*100, group=1),size=.4, col = "blue", linetype = "dashed") +
#   facet_wrap(~tier, ncol=3) +
#   ylab("Cover") + xlab("Year")+theme_bw()+
#   theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
#         axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
#         axis.title.x=element_text(size=11),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill = 'white')) +
#   scale_x_discrete(breaks= nth_element(unique(dat_pred_GRMF_without[[i]]$fYEAR),1,4))

# ggsave(p_traj, filename = paste0(title_of_run,"/report/extra/pred_traj_tier_no_data_GRMF", unique(dat_pred_without[[i]]$Model),".png"),
#        width=6, height=8)
# }




