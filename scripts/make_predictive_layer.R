## Make predictive layer 

# Import data
dat <- read.csv(paste0(title_of_run,"/data/reef_data_",surveys,"_with_NAs.csv")) %>%
  filter(fGROUP == "HCC") 
  
# Import reef locations 
 reef_loc <- st_read(paste0(title_of_run,"/data/reef_location.shp"))
 
# Import grid 
 
pred_layer <- st_read(paste0(title_of_run,"/data/grid_",grid_size*100,"K.shp"))%>%
   rename(tier = FID) 

# Import covariates 

dhw_pts <- read.csv(paste0(title_of_run,"/data/dhw.pts.df.csv")) # heatwave
cyc_pts <- read.csv(paste0(title_of_run,"/data/cyc.pts.df.csv")) # cyclone
other_pts <- read.csv(paste0(title_of_run,"/data/other.pts.df.csv")) # other 

# Manipulate initial grid 

# Rename tier 
pred_layer$tier <- pred_layer$tier + 1000

# Format repeated n year times  
pred.list <- list()

for (i in 1:length(unique(dat$fYEAR))){
  pred.list[[i]] <- pred_layer
}

pred_layer_all <- do.call(rbind, pred.list) %>%
  data.frame() %>%
  mutate(fYEAR = rep(unique(sort(dat$fYEAR)), each = nrow(pred_layer))) 

# Add reefid 

pred_reefid_prep <-  pred_layer %>%
  st_join(reef_loc) %>%
  dplyr::select(tier, Reef, geometry)

# Add other years
pred_reefid <- inner_join(pred_layer_all %>% as.data.frame(), 
                          pred_reefid_prep %>% as.data.frame(), by = "tier")%>%
  dplyr::select(.,-geometry.y)

pred_reefid <-  pred_reefid %>% 
  st_sf(sf_column_name = 'geometry.x')

rename(pred_reefid, geometry = geometry.x)

pred_reefid  <- pred_reefid %>% group_by(tier, fYEAR) %>%
  mutate(reefid_merged = as.factor(paste0(Reef, collapse = "_"))) %>%
  dplyr::select(tier, reefid_merged) %>%
  distinct() %>%
  st_drop_geometry()

HexPred_reefid2 <-inner_join(pred_layer_all %>% data.frame() , pred_reefid) %>%
  group_by(tier, fYEAR) %>% 
  filter(row_number()==1) %>%
  replace(is.na(.), 0) %>%
  st_as_sf(sf_column_name = "geometry") %>%
  rename(Reef = reefid_merged) %>%
  arrange(tier)

# Create covariates 

# DHW
dat_dhw <- extract_cov(pred_layer, dhw_pts) %>%
  mutate(cov = "DHW")

# Cyclone
dat_cyc <- extract_cov(pred_layer, cyc_pts) %>%
  mutate(cov = "CYC")

# DHW
dat_ot <- extract_cov(pred_layer, other_pts) %>%
  mutate(cov = "OT")

# All cov together 

cov_table <- rbind(dat_dhw, dat_cyc, dat_ot) %>%
  st_drop_geometry() %>%
        pivot_wider(
            names_from = cov, 
            values_from = mean_value) 


# Merge with predictive layer 

HexPred <- HexPred_reefid2 %>%
  filter(tier %in% cov_table$tier) %>% # some tiers are missing in the covariates layers due to the edging effect - need to be fixed. 
  left_join(cov_table)

st_write(HexPred, paste0(title_of_run,"/data/hexpred_cov.shp"), delete_dsn = T)

# Extract covariates at data locations 

# Step 1 - add tier
dat_sf <- dat %>%
dplyr::select(P_CODE:LONGITUDE, fYEAR, fDEPTH, fGROUP, COUNT, TOTAL, TRUE_COUNT) %>%
st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = st_crs(4326)) %>%
mutate_at("fYEAR", as.numeric)  %>%
st_join(pred_layer)


dat_sf2 <- dat %>%
dplyr::select(P_CODE:LONGITUDE, fYEAR, fDEPTH, fGROUP, COUNT, TOTAL, TRUE_COUNT) %>%
st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = st_crs(4326)) %>%
mutate_at("fYEAR", as.numeric)  %>%
st_join(HexPred)

# Step 2 - add covariates associated to tier level 

hex <- HexPred %>% 
dplyr::select(., - Reef) %>%
  st_drop_geometry() %>%
  data.frame()

dat_join <- dat_sf %>%
mutate(LONGITUDE = st_coordinates(.)[,1],
       LATITUDE = st_coordinates(.)[,2]) %>% 
       st_drop_geometry() %>%
left_join(hex) %>%
  arrange(REEF_NAME, SITE_NO, TRANSECT_NO, fYEAR) %>%
filter(! is.na(DHW)) # some data may be outside the spatial domain of the predlayer (at the edge)

map_dbl(dat_join, ~ sum(is.na(.)))

write_csv(dat_join, file=paste0(title_of_run,"/data/reef_data_NAs_with_cov_", surveys, ".csv"))
