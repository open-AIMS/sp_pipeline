## Make synthetic grid

## Read data and transform in sf 
X <- read.csv(paste0(title_of_run,"/data/reef_data_",surveys,"_with_NAs.csv")) %>%
  filter(!is.na(COUNT))
X_sf <- st_as_sf(X, coords = c("LONGITUDE", "LATITUDE"), crs = st_crs(4326))

## Import spatial domain and reefs 

reef_loc <- st_read(paste0(title_of_run,"/data/reef_location.shp"))

## Make grid - 10k

pred_layer <- st_make_grid(reef_loc, cellsize = grid_size, square = FALSE)[reef_loc]

p_reef <- ggplot() + 
  geom_sf(data = reef_loc, fill = "gray66") + 
  geom_sf(data = X_sf, col = "red", size = 0.5) +
  theme_bw() + 
  xlab("longitude") + ylab("latitude")

p_grid <- ggplot() + 
  geom_sf(data = pred_layer, fill = NA) + 
  geom_sf(data = X_sf, col = "red", size = 0.5) +
  theme_bw() + 
  xlab("longitude") + ylab("latitude")

p_loc <- p_reef + p_grid + plot_annotation(tag_levels = 'A')


ggsave(filename = paste0(title_of_run,"/figures/grid_reef_loc_",grid_size*100,"K.png"),
       plot = p_loc,  width=12, height=8)

st_write(pred_layer, paste0(title_of_run,"/data/grid_",grid_size*100,"K.shp"), delete_dsn = TRUE)

