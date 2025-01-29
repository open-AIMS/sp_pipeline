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

spatial.domain %>% 
  ggplot() +
  geom_sf() +
  theme_bw()
## ----end

## ---- SpatialPoints
set.seed(seed)
grid <- spatial.domain %>%
  st_set_crs(NA) %>%
  st_sample(size = 10000, type = 'regular') %>%
  st_set_crs(4236)

## Compile the spatial data frame - note it is important for RFsimulate
## that the data be sorted by Longitude then Latitude
spatial.grid.pts.df <- grid %>%
  st_coordinates() %>%
  as.data.frame() %>%
  dplyr::rename(Longitude = X, Latitude = Y) %>%
  arrange(Longitude, Latitude)

## Generate reefs ------------------------------------------------------
## Create a random field................................................
## ---- reefsRF
RFoptions(seed = 1)
threshold = 1.75
model <- RMexp(var = 1, scale = 0.1)
sim <- RFsimulate(model,
                  x = as.vector(scale(spatial.grid.pts.df$Longitude, 
                                      scale = FALSE)),
                  y = as.vector(scale(spatial.grid.pts.df$Latitude, 
                                      scale = FALSE)))
## combine with spatial data
reefs <- spatial.grid.pts.df  %>%
  mutate(Y = as.vector(sim))
reefs %>%
  ggplot() +
  geom_tile(aes(y = Latitude, x = Longitude, fill = Y)) +
  coord_sf(crs = 4326) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = c(0.95,0.95),
        legend.justification = c(1,1))

reefs %>%
  mutate(Y = Y > threshold) %>%
  ggplot() +
  geom_tile(aes(y = Latitude, x = Longitude, fill = Y)) +
  coord_sf(crs = 4326) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = c(0.95,0.95),
        legend.justification = c(1,1))
## ----end

## ---- reefsRFpolygons
reefs.sf <- reefs %>% 
  st_as_sf(coords = c('Longitude','Latitude')) %>%
  filter(Y > threshold) %>%
  st_buffer(0.05, endCapStyle = 'SQUARE') %>%
  st_cast('POLYGON') %>%
  st_union() %>%
  st_set_crs(4326)
## ----end

## ---- reefsRFpolygonsHollow
sf_use_s2(FALSE)  # negative buffers dont work if this is true
reefs.sf <- reefs.sf %>%
  st_buffer(0.01) %>%
  st_difference(reefs.sf %>% st_buffer(-0.01))
sf_use_s2(TRUE)
reefs.sf %>% ggplot() +
  geom_sf() +
  coord_sf(xlim = c(2.4,2.9), ylim = c(-12.05, -11.25)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1))

reefs.poly.sf <- reefs.sf %>% 
  st_cast('POLYGON') %>% 
  st_as_sf() %>%
  mutate(Reef = paste0('Reef', 1:n()))


g_domain <-  ggplot() +
  geom_sf(data = spatial.domain, fill = NA) + 
  geom_sf(data = reefs.poly.sf) +
  theme_bw() +
  theme(axis.title = element_text(size=12),
        legend.position = c(0.95,0.95),
        legend.justification = c(1,1)) + 
  xlab("longitude") + ylab("latitude")

st_write(reefs.poly.sf, paste0(title_of_run,"/data/reef_location.shp"), delete_dsn = TRUE)

## ----end

## Generate synthetic (simulated) data ==================================
## The response data will be effected by the following:
## - base coral cover - the global average coral cover (pooled over space and time)
## - spatial pattern in this base cover which reflects the spatial pattern at T0
## - annual growth (e.g. 5-10% annual increase) 
## - influence of covariates (spatio-temporal effects)
## - random noise

## ---- SyntheticData_Spatial.mesh
variance <- 1
kappa <- 1
alpha <- 2
mesh.pars <- c(1, 0.5, 0.1, 1, 0.5)*sqrt(alpha-ncol(spatial.grid.pts.df )/2)/kappa 
s = inla.mesh.segment(spatial.grid.pts.df [chull(spatial.grid.pts.df ), ])
mesh = inla.mesh.2d(
  spatial.grid.pts.df [chull(spatial.grid.pts.df ), ], 
  max.edge = mesh.pars[1:2], 
  cutoff = mesh.pars[3], 
  offset = mesh.pars[4:5],
  boundary = s)
## ----end

## ---- SyntheticData_Spatial.spde2
spde <- inla.spde2.matern(mesh, alpha = alpha)
## ----end

## ---- SyntheticData_Spatial.precision.matrix
theta <- c(-0.5*log(4*pi*variance*kappa^2), log(kappa))
Q <- inla.spde2.precision(spde, theta = theta)
## ----end

## Calculate a lattice projection to and from the mesh
## ---- SyntheticData_Spatial.A
A <- inla.spde.make.A(mesh = mesh, 
                      loc = as.matrix(spatial.grid.pts.df ))
## ----end

## Synthetic DHW --------------------------------------------------------
## ---- SyntheticData_CovariatesDHW.temporal.trend
set.seed(seed)
dhw.temporal <- data.frame(Year = years) %>%
  mutate(cYear = Year-1, 
         Y = 0.2*cYear + sin(cYear),
         Y = Y*rbeta(length(years), Y, 1),
         Y = scales::rescale(Y-min(Y), to = c(0,5)))
g_temp_dhw <- dhw.temporal %>%
  ggplot(aes(y = Y, x = Year)) +
  geom_line() +
  theme_bw(base_size = 12) + 
  ylab("DHW")
## ----end

## ---- SyntheticData_CovariatesDHW.effect
set.seed(seed)
dhw.sample <- inla.qsample(length(years), 
                           Q,
                           seed = seed, 
                           constr = spde$f$extraconstr)

rho <- rep(0.7,length(years))
rho <- rbeta(length(years), 0.2, 1)
x <- dhw.sample
for (j in 2:length(years)) {
  x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * dhw.sample[, j]
}
x <- sweep(x, 2, dhw.temporal$Y, FUN = '+')
dhw.effects <- scales::rescale(x, to = c(0,1))
dhw.pts.sample <- inla.mesh.project(mesh,
                                    loc = as.matrix(spatial.grid.pts.df [,1:2]),
                                    dhw.effects)

dhw.pts.effects.df = dhw.pts.sample %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  cbind(spatial.grid.pts.df) %>% 
  pivot_longer(cols = c(-Longitude,-Latitude),
               names_to = c('Year'),
               names_pattern = 'sample:(.*)',
               values_to = 'Value') %>%
  mutate(Year = as.numeric(Year))

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

## Synthetic CYC --------------------------------------------------------
## ---- CovariatesCYC.temporal.trend
set.seed(seed)
cyc <- vector('list', length(years))

for (yr in years) {
  #cat(paste("Year:", yr, '\n'))
  cyc.occur <- rbinom(1, 1, prob = min(0.05*yr^2, 0.6))
  #cat(paste("Cyclone Occurance:", cyc.occur, '\n'))
  cyc.intensity <- rbeta(1, 2, 1) %>% round(2)
  #cat(paste("Cyclone intensity:", cyc.intensity, '\n'))
  lat_offset <- runif(1,0,5)
  cyc.spatial <- mesh$loc[,1:2] %>%
    as.data.frame() %>%
    dplyr::select(Longitude = V1, Latitude = V2) %>%
    mutate(clong = as.vector(scale(Longitude, scale = FALSE)),
           clat = as.vector(scale(Latitude, scale = FALSE)),
           Y = lat_offset + runif(1,-1,1)*clong + runif(1,-1,1)*clat + sin(clat),
           Y = abs(Y),
           Y = ifelse(Y > cyc.intensity, cyc.intensity, Y),
           Y = cyc.intensity-Y,
           Value = Y*cyc.occur
    )
  cyc[[yr]] <- cyc.spatial %>% mutate(Year = yr)
}
cyc = do.call('rbind', cyc)
cyc.effects.df <- cyc %>%
  mutate(Value = scales::rescale(Value, to = c(0, 1)))

g <- cyc.effects.df %>%
  group_by(Year) %>%
  summarise(Mean = mean(Value),
            Median = median(Value)) %>%
  ggplot(aes(x = Year)) +
  geom_line(aes(y = Mean), color = 'blue') +
  theme_bw()
## ----end

## ---- SyntheticData_CovariatesCYC.effects
cyc.effects <- cyc.effects.df %>%
  dplyr::select(-clong, -clat, -Y) %>%
  pivot_wider(id_cols = c(Longitude, Latitude), 
              names_prefix = 'sample:',
              names_from = Year, 
              values_from = Value) %>%
  dplyr::select(-Longitude, -Latitude) %>%
  as.matrix

cyc.pts.sample <- inla.mesh.project(mesh,
                                    loc = as.matrix(spatial.grid.pts.df [,1:2]),
                                    cyc.effects)

cyc.pts.effects <- cyc.pts.sample %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  cbind(spatial.grid.pts.df) %>% 
  pivot_longer(cols = c(-Longitude,-Latitude),
               names_to = c('Year'),
               names_pattern = 'sample:(.*)',
               values_to = 'Value') %>%
  mutate(Year = as.numeric(Year))

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
## ----end

## Synthetic Other --------------------------------------------------------
## ---- SyntheticData_CovariatesOther.effects
set.seed(seed+1)
other.sample <- inla.qsample(length(years), 
                             Q,
                             seed = seed+1, 
                             constr = spde$f$extraconstr)

rho <- rep(0.7, length(years))
rho <- rbeta(length(years), 0.2, 1)
x <- other.sample
for (j in 2:length(years)) {
  x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * other.sample[, j]
}
other.effects <- scales::rescale(x, to = c(0,1))
other.pts.sample <- inla.mesh.project(mesh,
                                      loc = as.matrix(spatial.grid.pts.df [,1:2]),
                                      other.effects)

other.pts.effects = other.pts.sample %>% 
  as.matrix() %>% 
  as.data.frame() %>%
  cbind(spatial.grid.pts.df ) %>% 
  pivot_longer(cols = c(-Longitude, -Latitude),
               names_to = c('Year'),
               names_pattern = 'sample:(.*)',
               values_to = 'Value') %>%
  mutate(Year = as.numeric(Year))

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
## ----end

## Compile all disturbances------------------------------------------------
disturb.effects <-
  (dhw_eff*dhw.effects) + (cyc_eff*cyc.effects) + (ot_eff*other.effects) %>%
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
    Growth.HCC = 0.3,                ## Add growth onto this
    Growth.SC = 0.3,
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
data.reefs.sf <- reefs.sf %>%
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
  st_intersection(reefs.poly.sf)
sf_use_s2(TRUE)
## ----end

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
  st_intersection(reefs.poly.sf)
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
  st_intersection(reefs.poly.sf)
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

## Random locs ....................................................
# ## ---- randomLocs
 set.seed(seed)

## Select nSites within each of the reefs per year
data.random_locs.sf <- data.reefs.pts.sf %>%
  group_by(Reef, Year) %>%
  sample_n(nSites) %>%
  mutate(Site = paste0("S", 1:n())) %>%
  arrange(Year, Reef) %>%
  ungroup()

#p_grid <- ggplot() + 
#  geom_sf(data = data.reefs.pts.sf, fill = NA) + 
#  facet_wrap(~Year) +
#  geom_sf(data = data.random_locs.sf, col = "red", size = 0.5) +
#  theme_bw() + 
#  xlab("longitude") + ylab("latitude")


#ggsave(filename = paste0(title_of_run,"/figures/random_loc.png"),
#       plot = p_grid,  width=12, height=8)

# Finer sampling design components ===============================
## Random locations ....................................................

set.seed(seed)
 
data_random_locs_obs <- data.random_locs.sf %>%
  bind_cols(data.random_locs.sf %>%
              st_coordinates() %>%
              as.data.frame() %>%
              dplyr::rename(Longitude = X, Latitude = Y)) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  group_by(Longitude, Latitude, Reef) %>%
  crossing(Transect = paste0("T",
    1:Number_of_transects_per_site)) %>%
  group_by(Site, .add = TRUE) %>%
  mutate(
    SiteEffects_HCC = rnorm(1, 0, hcc_site_sigma),
    SiteEffects_SC = rnorm(1, 0, sc_site_sigma),
    SiteEffects_MA = rnorm(1, 0, ma_site_sigma)
  ) %>%
  group_by(Transect, .add = TRUE) %>%
  mutate(
    TransectEffects_HCC = rnorm(1, 0, hcc_transect_sigma),
    TransectEffects_SC = rnorm(1, 0, sc_transect_sigma),
    TransectEffects_MA = rnorm(1, 0, ma_transect_sigma)
  ) %>%
  ungroup() %>%
  mutate(
    HCC1 = HCC + SiteEffects_HCC + TransectEffects_HCC +
      rnorm(n(), 0, hcc_sigma),
    HCC2 = 100*plogis(HCC1),
    SC1 = SC + SiteEffects_SC + TransectEffects_SC +
      rnorm(n(), 0, sc_sigma),
    SC2 = 100*plogis(SC1),
    MA1 = MA + SiteEffects_MA + TransectEffects_MA +
      rnorm(n(), 0, ma_sigma),
    MA2 = 100*plogis(MA1)
  ) %>%
  arrange(Reef, Site, Transect, Year) %>%
  dplyr::select(Reef, Longitude, Latitude, Site, Transect,
    Year, HCC = HCC2, SC = SC2, MA = MA2) %>%
  mutate(Year = 2021 - max(years) + Year,
    Date = as.POSIXct(paste0(Year, "-01-01 14:00:00")))


#
## The following are on a fold scale.
## Hence a value of 0.8, indicates that 
Depth_effect_multiplier <- 2

#Extract information about the depth 
d_info <- depth_info(Depths)

data_random_locs_obs <-
  data_random_locs_obs %>%
  tidyr::crossing(Depth= d_info) %>% 
  pivot_longer(
    cols = c(HCC, SC, MA),
    names_to = "Group",
    values_to = "Value"
  ) %>%
  group_by(Reef, Site, Transect, Year, Date) %>%
    mutate(Value = Value + rev(sort(Depth_effect_multiplier * scale(rnorm(Depths), center = TRUE, scale = FALSE)))) %>%
  ungroup()

data_random_locs_obs %>%
  head()

## ----end

## ---- SyntheticData_randomLocsObsFortifyData
## Need to split the percentage cover into point and frames
data_random_locs_obs <- data_random_locs_obs %>%
  group_by(Reef,Site,Transect,Year,Depth,Date) %>%
  mutate(Points = round(Number_of_frames_per_transect * Points_per_frame * (Value/sum(Value)),0),
         Points = ifelse(Points<0, 0, Points)) %>%
  tidyr::uncount(Points) %>%
  sample_n(n(), replace=FALSE) %>%
  mutate(POINT_NO = rep_len(1:Points_per_frame, length = n()),
         FRAME = rep(1:Number_of_frames_per_transect, each=Points_per_frame, length = n())) %>%
  ungroup() 

reef_data.synthetic_random <- data_random_locs_obs %>% 
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

# write_csv(reef_data.synthetic_fixed, file=("data/reef_data_random_raw.csv"))

# Prepare format for the model - create response variable TRUE_COUNT
reef_data.synthetic_random_ready_all <- reef_data.synthetic_random  %>%
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

Reefs_random <- reef_data.synthetic_random_ready_all  %>%
  dplyr::select(REEF_NAME, fYEAR) %>%
  distinct() %>%
  group_by(fYEAR) %>%
  sample_n(nLocs) %>%
  ungroup() 
  
Reefs_rand <- Reefs_random %>%
  mutate(REEF_YEAR = paste0(REEF_NAME, fYEAR)) %>%
  pull(REEF_YEAR)

reef_data.synthetic_random_ready <- reef_data.synthetic_random_ready_all %>%
mutate(REEF_YEAR = paste0(REEF_NAME, fYEAR)) %>%
mutate(COUNT = case_when(.$REEF_YEAR %in% Reefs_rand ~ TRUE_COUNT, TRUE ~ NA_real_))  %>%
  dplyr::select(P_CODE:LONGITUDE,fYEAR,fDEPTH, fGROUP, TRUE_COUNT,TOTAL, COUNT)

map(reef_data.synthetic_random_ready, ~sum(is.na(.)))

write_csv(reef_data.synthetic_random_ready, file=paste0(title_of_run,"/data/reef_data_random_with_NAs.csv"))

## ----end

## Data viz 

p_vis_data_random <- ggplot(reef_data.synthetic_random_ready %>% filter(!is.na(COUNT)) %>% filter(fGROUP == "HCC")) + 
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
  ggtitle("Random design")

ggsave(filename = paste0(title_of_run,"/report/extra/trend_data_",surveys,".png"),
       plot = p_vis_data_random, width=13, height=12)  
