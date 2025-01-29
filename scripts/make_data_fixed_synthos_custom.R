# For HPC
## Generate the spatial and temporal domains ==================================
## Create the full spatial domain
spatial.domain <- st_geometry(
  st_multipoint(
    x = rbind(
      c(0, -11),
      c(3,-11),
      c(6,-14),
      c(1,-15),
      c(2,-12),
      c(0,-11)
    )
  )
) %>%
  st_set_crs(4326) %>%
  st_cast('POLYGON')

## ----end

# Using Synthos functions 
## ---- SpatialPoints
set.seed(config_sp$seed)
spatial.grid <- spatial.domain %>%
  st_set_crs(NA) %>%
  st_sample(size = 10000, type = 'regular') %>%
  st_set_crs(4236)

## Compile the spatial data frame - note it is important for RFsimulate
## that the data be sorted by Longitude then Latitude
spatial.grid.pts.df <- spatial.grid %>%
  st_coordinates() %>%
  as.data.frame() %>%
  dplyr::rename(Longitude = X, Latitude = Y) %>%
  arrange(Longitude, Latitude)

simulated_field_sf  <- generate_field(spatial.grid, config_sp)
simulated_patches_sf <- generate_patches(simulated_field_sf, config_sp)
reefs.sf <- generate_reefs(simulated_patches_sf, config_sp)

st_write(reefs.sf$simulated_reefs_poly_sf, paste0(title_of_run,"/data/reef_location.shp"), delete_dsn = TRUE)

## ----end

## Generate synthetic (simulated) data ==================================
## The response data will be effected by the following:
## - base coral cover - the global average coral cover (pooled over space and time)
## - spatial pattern in this base cover which reflects the spatial pattern at T0
## - annual growth (e.g. 5-10% annual increase) 
## - influence of covariates (spatio-temporal effects)
## - random noise

## ---- SyntheticData_Spatial.mesh
mesh <- create_spde_mesh(spatial.grid,config_sp)
spde <- create_spde(spatial.grid, config_sp)

######### DHW 
# pts effects 
dhw.pts.effects.df <- disturbance_dhw(spatial.grid,spde,config_sp)$dhw_pts_effects_df

g_spat_dhw <- ggplot(dhw.pts.effects.df, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  facet_wrap(~Year, ncol = 4) +
  scale_fill_gradientn(colors = rev(heat.colors(10))) + 
  coord_sf(crs = 4236) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank())

ggsave(filename = paste0(title_of_run,"/report/extra/dhw.png"),
       plot = g_spat_dhw, width=14, height=6)

write.csv(dhw.pts.effects.df, file = paste0(title_of_run,"/data/dhw.pts.df.csv"), row.names = F)

# effects 
dhw.effects <- disturbance_dhw(spatial.grid,spde,config_sp)$dhw_effects

######### CYC
cyc.pts.effects <- disturbance_cyc(spatial.grid,spde,config_sp)$cyc_pts_effects

g_spat_cyc <- ggplot(cyc.pts.effects, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  facet_wrap(~Year, ncol = 4) +
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  coord_sf(crs = 4236) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank())

ggsave(filename = paste0(title_of_run,"/report/extra/cyclone.png"),
       plot = g_spat_cyc, width=14, height=6)

write.csv(cyc.pts.effects, file = paste0(title_of_run,"/data/cyc.pts.df.csv"), row.names = F)

# effects
cyc.effects <- disturbance_cyc(spatial.grid,spde,config_sp)$cyc_effects

########## OT
# pts effects
other.pts.effects <- disturbance_other(spatial.grid,spde,config_sp)$other_pts_effects

g_spat_ot <- ggplot(other.pts.effects, aes(y = Latitude, x = Longitude)) +
  geom_tile(aes(fill = Value)) +
  facet_wrap(~Year, ncol = 4) +
  scale_fill_gradientn(colors = terrain.colors(10)) + 
  coord_sf(crs = 4236) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank())

ggsave(filename = paste0(title_of_run,"/report/extra/other_dist.png"),
       plot = g_spat_ot, width=14, height=6)

write.csv(other.pts.effects, file = paste0(title_of_run,"/data/other.pts.df.csv"), row.names = F)

# effects
other.effects <- disturbance_other(spatial.grid,spde,config_sp)$other_effects


## Compile all disturbances------------------------------------------------
disturb.effects <-
  (config_sp$dhw_weight*dhw.effects) + (config_sp$cyc_weight*cyc.effects) + (config_sp$other_weight*other.effects) %>%
  as.data.frame()
all.effects.df <- mesh$loc[,1:2] %>%
  as.data.frame() %>%
  dplyr::rename(Longitude = V1, Latitude = V2) %>%
  cbind(disturb.effects) %>%
  pivot_longer(cols = c(-Longitude, -Latitude),
               names_to = 'Year',
               names_pattern = 'sample:(.*)',
               values_to = 'Y') %>%
  mutate(Year = factor(Year, levels=sort(unique(as.numeric(as.character(Year)))))) %>%
  group_by(Longitude, Latitude) %>%
  mutate(
    Growth.HCC =  config_sp$hcc_growth,                ## Add growth onto this
    Growth.SC =  config_sp$sc_growth,
    Y.HCC = cumsum(-Y + Growth.HCC), ## cumsum on link scale will accumulate effects
    Y.SC = cumsum(-Y + Growth.SC)
  )
## ----end

## Synthetic Baselines ----------------------------------------------------
## ---- SyntheticData_BaselineSpatial.HCC
baseline.sample.hcc <- mesh$loc[,1:2] %>%
  as.data.frame() %>%
  dplyr::select(Longitude = V1, Latitude = V2) %>%
  mutate(clong = as.vector(scale(Longitude, scale = FALSE)),
         clat = as.vector(scale(Latitude, scale = FALSE)),
         Y = clong + sin(clat) + #rnorm(1,0,1) +
           1.5*clong + clat) %>%
  mutate(Y = scales::rescale(Y, to = c(-2, 0.8)))

baseline.effects.hcc <- baseline.sample.hcc %>%
  dplyr::select(Y) %>%
  as.matrix
baseline.pts.sample.hcc <- inla.mesh.project(mesh,
                                             loc=as.matrix(spatial.grid.pts.df [,1:2]),
                                             baseline.effects.hcc)
baseline.pts.effects.hcc = baseline.pts.sample.hcc %>% 
  cbind() %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  cbind(spatial.grid.pts.df ) %>% 
  pivot_longer(cols = c(-Longitude,-Latitude),
               names_to = c('Year'),
               names_pattern = 'sample:(.*)',
               values_to = 'Value') %>%
  mutate(Year = as.numeric(Year))
## ----end

## ---- SyntheticData_BaselineSpatial.SC
baseline.sample.sc <- mesh$loc[,1:2] %>%
  as.data.frame() %>%
  dplyr::select(Longitude = V1, Latitude = V2) %>%
  mutate(clong = as.vector(scale(Longitude, scale = FALSE)),
         clat = as.vector(scale(Latitude, scale = FALSE)),
         Y = clong + sin(clat) + #rnorm(1,0,1) +
           1.5*clong + -1.5*clat) %>%
  mutate(Y = scales::rescale(Y, to = c(-4, -2)))

baseline.effects.sc <- baseline.sample.sc %>%
  dplyr::select(Y) %>%
  as.matrix
baseline.pts.sample.sc <- inla.mesh.project(mesh,
                                            loc=as.matrix(spatial.grid.pts.df [,1:2]),
                                            baseline.effects.sc)
baseline.pts.effects.sc = baseline.pts.sample.sc %>% 
  cbind() %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  cbind(spatial.grid.pts.df ) %>% 
  pivot_longer(cols = c(-Longitude,-Latitude),
               names_to = c('Year'),
               names_pattern = 'sample:(.*)',
               values_to = 'Value') %>%
  mutate(Year = as.numeric(Year))
## ----end

## Synthetic data ---------------------------------------------------------
## Hard coral

all.effects.hcc <- all.effects.df %>%
  full_join(baseline.sample.hcc %>% dplyr::select(Longitude, Latitude, BASE.HCC=Y)) %>%
  group_by(Longitude, Latitude) %>%
  mutate(HCC = BASE.HCC + Y.HCC) %>%
  ungroup() %>%
  dplyr::select(Longitude, Latitude, Year,HCC) %>%
  tidyr::pivot_wider(
    names_prefix = 'sample:',
    names_from = Year, 
    values_from = HCC) %>%
  dplyr::select(-Longitude, -Latitude) %>%
  as.matrix

all.pts.sample.hcc <- inla.mesh.project(mesh,
                                        loc = as.matrix(spatial.grid.pts.df [,1:2]),
                                        all.effects.hcc)
all.pts.effects.hcc = all.pts.sample.hcc %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  cbind(spatial.grid.pts.df ) %>% 
  pivot_longer(cols = c(-Longitude, -Latitude),
               names_to = c('Year'),
               names_pattern = 'sample:(.*)',
               values_to = 'Value') %>%
  mutate(Year = as.numeric(Year),
         Value = Value)
## ----end

## Soft coral
## ---- SyntheticData_CompileSyntheticData.SC
## Do all this on the link scale so that can use cumsum
all.effects.sc <- all.effects.df %>%
  full_join(baseline.sample.sc %>% dplyr::select(Longitude, Latitude, BASE.SC=Y)) %>%
  group_by(Longitude, Latitude) %>%
  mutate(SC = BASE.SC + Y.SC) %>%
  ungroup() %>%
  dplyr::select(Longitude, Latitude, Year,SC) %>%
  pivot_wider(id_cols = c(Longitude, Latitude), 
              names_prefix = 'sample:',
              names_from = Year, 
              values_from = SC) %>%
  dplyr::select(-Longitude, -Latitude) %>%
  as.matrix()

all.pts.sample.sc <- inla.mesh.project(mesh,
                                       loc = as.matrix(spatial.grid.pts.df [,1:2]),
                                       all.effects.sc)
all.pts.effects.sc = all.pts.sample.sc %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  cbind(spatial.grid.pts.df ) %>% 
  pivot_longer(cols = c(-Longitude, -Latitude),
               names_to = c('Year'),
               names_pattern = 'sample:(.*)',
               values_to = 'Value') %>%
  mutate(Year = as.numeric(Year),
         Value = Value)
## ----end

## Macroalgae
## ---- SyntheticData_CompileSyntheticData.MA
## Do all this on the link scale so that can use cumsum
all.pts.effects.ma <- all.pts.effects.hcc %>%
  dplyr::rename(HCC=Value) %>% 
  bind_cols(all.pts.effects.sc %>%
              dplyr::select(SC=Value)) %>%
  mutate(Total.Avail = 0.8 - plogis(HCC) + plogis(SC),
         MA = Total.Avail,
         Value = qlogis(MA)) %>%
  dplyr::select(-HCC, -SC, -Total.Avail, -MA)
## ----end

## Broad scale reef patterns ==================================================
## INLA ----------------------------------------------------------
## - rasterize the reefs frame
## - convert to points (centroids of raster cells)
## - filter to the values of 1
## - extract coordinates
## - convert to data frame

## ---- SyntheticData_PointsInReefs
data.reefs.sf <- reefs.sf$simulated_reefs_sf %>%
  st_as_stars(dx = 0.01) %>%  # rasterize
  st_as_sf(as_points = TRUE) %>%
  filter(values == 1L)

data.reefs.df = data.reefs.sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  dplyr::rename(Longitude = X, Latitude = Y)
## ----end

## ---- SyntheticData_ProjectOntoReefs.HCC
data.reefs.sample.hcc <- inla.mesh.project(mesh,
                                           loc = as.matrix(data.reefs.df[,1:2]),
                                           all.effects.hcc)
data.reefs.hcc = data.reefs.sample.hcc %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  cbind(data.reefs.df) %>% 
  pivot_longer(cols = c(-Longitude,-Latitude),
               names_to = c('Year'),
               names_pattern = 'sample:(.*)',
               values_to = 'Value') %>%
  mutate(Year = as.numeric(Year),
         Value = Value)

data.reefs.pts.hcc.sf = data.reefs.hcc %>%
  st_as_sf(coords = c('Longitude', 'Latitude')) %>%
  st_set_crs(st_crs(data.reefs.sf))
sf_use_s2(FALSE)

data.reefs.pts.hcc.sf <- data.reefs.pts.hcc.sf %>%
  st_intersection(reefs.sf$simulated_reefs_poly_sf)
sf_use_s2(TRUE)
## ----en

## ---- SyntheticData_ProjectOntoReefs.SC
data.reefs.sample.sc <- inla.mesh.project(mesh,
                                          loc = as.matrix(data.reefs.df[,1:2]),
                                          all.effects.sc)

data.reefs.sc = data.reefs.sample.sc %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  cbind(data.reefs.df) %>% 
  pivot_longer(cols = c(-Longitude,-Latitude),
               names_to = c('Year'),
               names_pattern = 'sample:(.*)',
               values_to = 'Value') %>%
  mutate(Year = as.numeric(Year),
         Value = Value)

data.reefs.pts.sc.sf = data.reefs.sc %>%
  st_as_sf(coords = c('Longitude', 'Latitude')) %>%
  st_set_crs(st_crs(data.reefs.sf))
sf_use_s2(FALSE)

data.reefs.pts.sc.sf <- data.reefs.pts.sc.sf %>%
  st_intersection(reefs.sf$simulated_reefs_poly_sf)
sf_use_s2(TRUE)
## ----end

## ---- SyntheticData_ProjectOntoReefs.MA
data.reefs.ma <- data.reefs.hcc %>%
  rename(HCC = Value) %>%
  full_join(data.reefs.sc %>% rename(SC = Value)) %>%
  mutate(Total.Avail = 0.8 - plogis(HCC) + plogis(SC),
         MA = Total.Avail,
         Value = qlogis(MA)) %>%
  dplyr::select(-HCC, -SC, -Total.Avail, -MA)

data.reefs.pts.ma.sf = data.reefs.ma %>%
  st_as_sf(coords = c('Longitude', 'Latitude')) %>%
  st_set_crs(st_crs(data.reefs.sf))
sf_use_s2(FALSE)

data.reefs.pts.ma.sf <- data.reefs.pts.ma.sf %>%
  st_intersection(reefs.sf$simulated_reefs_poly_sf)
sf_use_s2(TRUE)
## ----end

## Sampling designs (large scale components) ==================================

## ---- ProjectOntoReefs

## ---- SyntheticData_ProjectOntoReefs
data.reefs.pts.sf <- data.reefs.pts.hcc.sf %>%
  rename(HCC = Value) %>%
  bind_cols(data.reefs.pts.sc.sf %>%
              dplyr::select(SC = Value) %>%
              st_drop_geometry()
  ) %>% 
  bind_cols(data.reefs.pts.ma.sf %>%
              dplyr::select(MA = Value) %>%
              st_drop_geometry()
  ) 

st_write(data.reefs.pts.sf, paste0(title_of_run,"/data/true_sp_field.shp"), delete_dsn = TRUE)

## ----end

## Fixed locs ....................................................
## Selecting n_sites for all Reefs

data.fixed_locs.sf <- data.reefs.pts.sf %>%
  dplyr::select(Reef, geometry) %>%
  distinct(.keep_all=TRUE) %>% 
  group_by(Reef) %>%
  sample_n(config_lrge$n_sites) %>%
  mutate(Site = paste0('S',1:n())) %>%
  ungroup %>% 
  st_join(data.reefs.pts.sf %>% 
            dplyr::select(-Reef))

## ----end

## Finer sampling design components ===============================
## Fixed locations ....................................................
## ---- SyntheticData_fixedLocsObs
set.seed(config_lrge$seed)
data.fixed_locs.obs <- data.fixed_locs.sf %>%
  bind_cols(data.fixed_locs.sf %>%
              st_coordinates() %>%
              as.data.frame %>%
              dplyr::rename(Longitude = X, Latitude = Y)) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  group_by(Longitude, Latitude, Reef) %>%
  crossing(Transect = paste0('T',1:config_fine$Number_of_transects_per_site)) %>%
  group_by(Site, .add = TRUE) %>%
  mutate(
    SiteEffects.HCC = rnorm(1, 0, config_fine$hcc_site_sigma),
    SiteEffects.SC = rnorm(1, 0, config_fine$sc_site_sigma),
    SiteEffects.MA = rnorm(1, 0, config_fine$ma_site_sigma)
  ) %>%
  group_by(Transect, .add = TRUE) %>%
  mutate(
    TransectEffects.HCC = rnorm(1, 0, config_fine$hcc_transect_sigma),
    TransectEffects.SC = rnorm(1, 0, config_fine$sc_transect_sigma),
    TransectEffects.MA = rnorm(1, 0, config_fine$ma_transect_sigma)
  ) %>%
  ungroup() %>%
  mutate(
    HCC1 = HCC + SiteEffects.HCC + TransectEffects.HCC + rnorm(n(), 0, config_fine$hcc_sigma),
    HCC2 = 100*plogis(HCC1),
    SC1 = SC + SiteEffects.SC + TransectEffects.SC + rnorm(n(), 0, config_fine$sc_sigma),
    SC2 = 100*plogis(SC1),
    MA1 = MA + SiteEffects.MA + TransectEffects.MA + rnorm(n(), 0, config_fine$ma_sigma),
    MA2 = 100*plogis(MA1)
  ) %>%
  arrange(Reef, Site, Transect, Year) %>%
  dplyr::select(Reef, Longitude, Latitude, Site, Transect, Year, HCC = HCC2, SC = SC2, MA = MA2) %>%
  mutate(Year = as.numeric(format(Sys.Date(), "%Y")) - max(config_fine$years) + Year,
         Date = as.POSIXct(paste0(Year, '-01-01 14:00:00')))
## ----end
## ---- SyntheticData_fixedLocsObsDepths
## The following are on a fold scale.
## Hence a value of 0.8, indicates that 
Depth_effect.multiplier <- 2

#Extract information about the depth 
d_info <- depth_info(config_fine$Depths)

data.fixed_locs.obs <- data.fixed_locs.obs %>%
  tidyr::crossing(Depth= d_info) %>% 
  pivot_longer(cols = c('HCC', 'SC', 'MA'), names_to = 'Group', values_to = 'Value') %>%
  group_by(Reef, Site, Transect, Year, Date) %>%
  mutate(Value = Value + rev(sort(Depth_effect.multiplier * scale(rnorm(config_fine$Depths), center = TRUE, scale = FALSE)))) %>%
  ungroup()
## ----end

## ---- SyntheticData_fixedLocsObsFortifyData
## Need to split the percentage cover into point and frames
data.fixed_locs.obs <- data.fixed_locs.obs %>%
  group_by(Reef,Site,Transect,Year,Depth,Date) %>%
  mutate(Points = round(config_pt$Number_of_frames_per_transect * config_pt$Points_per_frame * (Value/sum(Value)),0),
         Points = ifelse(Points<0, 0, Points)) %>%
  tidyr::uncount(Points) %>%
  sample_n(n(), replace=FALSE) %>%
  mutate(POINT_NO = rep_len(1:config_pt$Points_per_frame, length = n()),
         FRAME = rep(1:config_pt$Number_of_frames_per_transect, each=config_pt$Points_per_frame, length = n())) %>%
  ungroup() 

reef_data.synthetic_fixed <- data.fixed_locs.obs %>% 
  mutate(P_CODE = paste0("SYNTHETIC-", surveys),
         ID = 1:n(),
         REEF_NAME = Reef,
         LATITUDE = Latitude,
         LONGITUDE = Longitude,
         SITE_NO = Site,
         TRANSECT_NO = Transect,
         SITE_DEPTH = Depth,
         REEF_ZONE = '-',
         REPORT_YEAR = Year,
         SURVEY_DATE = Date,
         FRAME = paste0(P_CODE, '/', REEF_NAME, '/', REEF_ZONE, '/', SITE_NO, '/', SITE_DEPTH, '/', TRANSECT_NO, '/', REPORT_YEAR, '/', FRAME),
         POINT_NO = POINT_NO,
         GROUP_DESC = Group,
         BENTHIC_CATEGORY = paste0(Group,'_alt')
  ) %>%
  dplyr::select(P_CODE, ID, REEF_NAME,
                LATITUDE, LONGITUDE, SITE_NO, TRANSECT_NO, SITE_DEPTH,
                REEF_ZONE, REPORT_YEAR, SURVEY_DATE, FRAME, POINT_NO,
                GROUP_DESC, BENTHIC_CATEGORY)

# Prepare format for the model - create response variable TRUE_COUNT
reef_data.synthetic_fixed_ready_all <- reef_data.synthetic_fixed  %>%
  mutate(
    P_CODE = factor(P_CODE),
    ID = factor(ID),
    fYEAR = factor(REPORT_YEAR),
    SITE_DEPTH = ifelse(is.na(SITE_DEPTH),'_', SITE_DEPTH),  # replace missing depths with a non NA value
    fDEPTH = factor(SITE_DEPTH),
    across(c(SITE_NO, TRANSECT_NO, fYEAR, fDEPTH), function(x) factor(as.character(x))),
    DATE = as.Date(SURVEY_DATE, format = '%Y-%m-%d %h:%m:%s'),
    fGROUP = factor(GROUP_DESC)) %>%
  group_by(P_CODE, REEF_NAME, SITE_NO, TRANSECT_NO,
           LATITUDE, LONGITUDE,
           REPORT_YEAR, DATE, fYEAR, fDEPTH, REEF_ZONE,
           fGROUP, GROUP_DESC) %>%
  summarise(TRUE_COUNT = n()) %>%
  ungroup(fGROUP, GROUP_DESC) %>%
  mutate(TOTAL=sum(TRUE_COUNT)) %>%
  ungroup %>%
  filter(!is.na(REPORT_YEAR)) %>% droplevels() 

# Select nReef to be treated as observed locations and NAs to the others 
## Selecting n Reefs
Reefs.observed <- reef_data.synthetic_fixed_ready_all  %>%
  dplyr::select(REEF_NAME) %>%
  distinct() %>%
  sample_n(size = config_lrge$n_locs) %>%
  pull(REEF_NAME)

reef_data.synthetic_fixed_ready <- reef_data.synthetic_fixed_ready_all %>%
 mutate(COUNT = case_when(.$REEF_NAME %in% Reefs.observed ~ TRUE_COUNT, TRUE ~ NA_real_))  %>%
  dplyr::select(P_CODE:LONGITUDE,fYEAR,fDEPTH, fGROUP, TRUE_COUNT,TOTAL, COUNT)

map(reef_data.synthetic_fixed_ready, ~sum(is.na(.)))

write_csv(reef_data.synthetic_fixed_ready, file=paste0(title_of_run,"/data/reef_data_fixed_with_NAs.csv"))

## ----end

## Data viz 

p_vis_data_fixed <- ggplot(reef_data.synthetic_fixed_ready %>% filter(!is.na(COUNT)) %>% filter(fGROUP == "HCC")) + 
  geom_line(aes(x = fYEAR, y = (COUNT/TOTAL)*100, group = interaction(as.factor(TRANSECT_NO), as.factor(SITE_NO), REEF_NAME),
  col = as.factor(SITE_NO)), 
   show.legend = FALSE) + 
  facet_wrap(~REEF_NAME, ncol=4) + theme_bw() +
  labs(x = "Year", y = "Coral cover") +
  ylab("Coral cover") + xlab("Year")+theme_bw()+
  theme(axis.text.x = element_text(size=10, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=10),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10, margin = margin())) + 
  ggtitle("Fixed design")

ggsave(filename = paste0(title_of_run,"/report/extra/trend_data_",surveys,".png"),
       plot = p_vis_data_fixed, width=13, height=12)  
