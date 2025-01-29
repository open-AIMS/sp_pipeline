
# Functions 

##################################
############ Create folders structure of the pipeline 
##################################

make_folders <- function(title_of_run){
  
  wd <- getwd() 
  
  # lEVEL 1 - name of run 
  pathway = paste0(wd,"/", title_of_run, "/")
  dir.create(pathway)
  ifelse(dir.exists(title_of_run)!=TRUE,print("directory already exists - is it a new simulation?"), FALSE)
  
  # LEVEL 2 - subdirectories 
  
  create_subdir = function (x) {
    for (i in x){
      if(!dir.exists(paste0(pathway, i))){dir.create(paste0(pathway, i))}
    }}
  
  
  create_subdir(c("data","figures", "model_outputs", "report"))
  
  # LEVEL 3 - subsubdirectories 
  
  # within model_outputs 
  pathway_2 = paste0(wd,"/", title_of_run, "/model_outputs/")
  
  create_subsubdir = function (x) {
    for (i in x){
      if(!dir.exists(paste0(pathway_2, i))){dir.create(paste0(pathway_2, i))}
    }}
  
  create_subsubdir(c("full_run","leave_out", "predictions", "tuning"))
  
  # within report 
  pathway_3 = paste0(wd,"/", title_of_run, "/report/")

  create_subsubdir2 = function (x) {
  for (i in x){
    if(!dir.exists(paste0(pathway_3, i))){dir.create(paste0(pathway_3, i))}
  }}

  create_subsubdir2(c("extra","resources", "data"))
  
  # copying file in new folder 
  
  file.copy("toc_logo.html", paste0(pathway_3,"/resources/"))

}

##################################
########### Functions from the synthos package
#https://github.com/open-AIMS/synthos/tree/main
####################################

generate_field <- function(spatial_grid, config) {
  testthat::expect(
    inherits(spatial_grid, c("sf", "sfc")),
    "spatial_grid must be an sf object"
  )
  testthat::expect_in(
    sort(c("seed", "psill", "model", "range", "nugget")),
    sort(names(config))
  )

  set.seed(config$seed)
  ## create a variogram model object that can be used to simulate a
  ## random field
  vgm_model <- gstat::vgm(
    psill = config$psill,
    model = config$model,
    range = config$range,
    nugget = config$nugget
  )
  ## Simulate a random field using gstat with the sf object
  sim <- gstat::gstat(formula = z ~ 1, locations = spatial_grid,
    model = vgm_model,
    ## beta: for an intercept only model, it represents the intercept
    beta = 0,
    ## nmax: the number of nearest observations that should be used for a
    ## kriging prediction or simulation, where nearest is defined in
    ## terms of the space of the spatial locations
    nmax = 20,
    ## dummy: consider these data as a dummy variable
    dummy = TRUE
  )
  ## predict the random field
  simulated_field_sf <- predict(sim, newdata = spatial_grid, nsim = 1)
  return(simulated_field_sf)
}

generate_patches <- function(simulated_field_sf, config) {
  testthat::expect(
    inherits(simulated_field_sf, c("sf")),
    "simulated_field_sf must be an sf object"
  )
  testthat::expect_match(
    names(simulated_field_sf),
    "sim[0-9]*|geometry"
  )
  testthat::expect_in(
    sort(c("patch_threshold")),
    sort(names(config))
  )
  ## create a thresholded field
  sf_use_s2(FALSE)
  simulated_patches_sf <- simulated_field_sf |>
    dplyr::rename(Y = sim1) |> 
    sf::st_as_sf(coords = c("Longitude", "Latitude")) |>
    dplyr::filter(Y > config$patch_threshold) |>
    sf::st_buffer(0.05, endCapStyle = "SQUARE") |>
    sf::st_cast("POLYGON") |>
    sf::st_union() |>
    ## st_set_crs(4326)
    suppressMessages() |>
    suppressWarnings()
  sf_use_s2(TRUE)
  return(simulated_patches_sf)
}

generate_reefs <- function(simulated_patch_sf, config) {
  testthat::expect(
    inherits(simulated_patch_sf, c("sfc")),
    "simulated_patch_sf must be an sfc object"
  )
  testthat::expect_in(
    sort(c("reef_width")),
    sort(names(config))
  )
  ## buffer the outline to create a ribbon
  sf_use_s2(FALSE)
  simulated_reefs_sf <- simulated_patch_sf |>
    sf::st_buffer(config$reef_width) |>
    sf::st_difference(simulated_patch_sf |>
                        sf::st_buffer(-config$reef_width)) |> 
    sf::st_transform(crs = 4326) |>
    suppressMessages() |>
    suppressWarnings()
  simulated_reefs_poly_sf <- simulated_reefs_sf |>
    sf::st_cast("POLYGON") |>
    sf::st_sf(geometry = _) |>
    dplyr::mutate(Reef = paste0("Reef", 1:dplyr::n())) |>
    sf::st_transform(crs = 4326) |>
    suppressMessages() |>
    suppressWarnings()
  sf_use_s2(TRUE)
  return(list(
    simulated_reefs_sf = simulated_reefs_sf,
    simulated_reefs_poly_sf = simulated_reefs_poly_sf
  ))
}

#pointify_polygons <- function(reefs_sf) {
 # testthat::expect(
#    inherits(reefs_sf, c("sfc")),
#    "reefs_sf must be a sfc object"
#  )
#  data_reefs_sf <- reefs_sf |>
#    stars::st_as_stars(dx = 0.01) |>  # rasterize
#    sf::st_as_sf(as_points = TRUE) |>
#    dplyr::filter(values == 1L)

#  data_reefs_df <- data_reefs_sf |>
#    sf::st_coordinates() |>
 #   as.data.frame() |>
 #   dplyr::rename(Longitude = X, Latitude = Y)
 # list(data_reefs_sf = data_reefs_sf, data_reefs_df = data_reefs_df)
#}

spatial_grid_sfc_to_df <- function(spatial_grid) {
  spatial_grid |>
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y) |>
    dplyr::arrange(Longitude, Latitude)
}

create_spde_mesh <- function(spatial_grid, config) {
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)
  mesh_pars <- c(1, 0.5, 0.1, 1, 0.5) *
    sqrt(config$alpha - ncol(spatial_grid_pts_df) / 2) / config$kappa
  s <- INLA::inla.mesh.segment(
    spatial_grid_pts_df[chull(spatial_grid_pts_df), ]
  )
  mesh <- INLA::inla.mesh.2d(
    spatial_grid_pts_df[chull(spatial_grid_pts_df), ],
    max.edge = mesh_pars[1:2],
    cutoff = mesh_pars[3],
    offset = mesh_pars[4:5],
    boundary = s
  )
  return(mesh)
}

create_spde_matern <- function(spatial_grid, mesh, config) {
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  spde <- INLA::inla.spde2.matern(mesh, alpha = config$alpha)
  ## calculate precision matrix from the parameter values (theta)
  theta <- c(-0.5 * log(4 * pi * config$variance * config$kappa^2), log(config$kappa))
  Q <- INLA::inla.spde2.precision(spde, theta = theta)
  ## calculate a lattic projection to and from mesh
  A <- INLA::inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(spatial_grid_pts_df)
  )
  list(spde = spde, Q = Q, A = A)
}


create_spde <- function(spatial_grid, config) {
  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("alpha", "kappa", "variance")),
    sort(names(config))
  )
  mesh <- create_spde_mesh(spatial_grid, config)
  spde <- create_spde_matern(spatial_grid, mesh, config)
  list(mesh = mesh, spde = spde$spde, Q = spde$Q, A = spde$A)
}


disturbance_dhw <- function(spatial_grid, spde, config) {
  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("mesh", "spde", "Q", "A")),
    sort(names(spde))
  )
  testthat::expect(
    inherits(spde$spde, c("inla.spde")),
    "spde$spde must be a inla.spde object"
  )
  testthat::expect_in(
    sort(c("years", "seed")),
    sort(names(config))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  ## Overall temporal trend in DHW
  set.seed(config$seed)
  dhw_temporal <- data.frame(Year = config$years) |>
    dplyr::mutate(
      cYear = Year - 1, # as.vector(scale(Year, scale=FALSE)),
      Y = 0.2 * cYear + sin(cYear),
      Y = Y * rbeta(length(config$years), Y, 1),
      Y = scales::rescale(Y - min(Y), to = c(0, 5))
    )
  ## Propagate this temporal trend across a random rield with a time
  ## varying autocorrelation coefficient drawn from a beta distribution
  ## with shape parameters of 0.2 and 1
  set.seed(config$seed)
  dhw_sample <- INLA::inla.qsample(length(config$years),
    spde$Q,
    seed = config$seed,
    constr = spde$spde$f$extraconstr
  ) |>
    suppressMessages() |>
    suppressWarnings()

  rho <- rep(0.7, length(config$years))
  rho <- rbeta(length(config$years), 0.2, 1)
  x <- dhw_sample
  for (j in 2:length(config$years)) {
    x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * dhw_sample[, j]
  }
  x <- sweep(x, 2, dhw_temporal$Y, FUN = "+")
  dhw_effects <- scales::rescale(x, to = c(0, 1))
  dhw_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    dhw_effects
  )

  dhw_pts_effects_df <- dhw_pts_sample %>%
    as.matrix() %>%
    as.data.frame() %>%
    cbind(spatial_grid_pts_df) %>%
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) %>%
    dplyr::mutate(Year = as.numeric(Year))
  ## Value=scales::rescale(Value, to=c(0,1)))
  list(
    dhw_temporal = dhw_temporal,
    dhw_effects = dhw_effects,
    dhw_pts_sample = dhw_pts_sample,
    dhw_pts_effects_df = dhw_pts_effects_df
  )
}

disturbance_cyc <- function(spatial_grid, spde, config) {

  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("mesh", "spde", "Q", "A")),
    sort(names(spde))
  )
  testthat::expect(
    inherits(spde$spde, c("inla.spde")),
    "spde$spde must be a inla.spde object"
  )
  testthat::expect_in(
    sort(c("years", "seed")),
    sort(names(config))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  set.seed(config$seed)
  cyc <- vector("list", length(config$years))

  for (yr in config$years) {
    ## cat(paste("Year:", yr, "\n"))
    cyc_occur <- rbinom(1, 1, prob = min(0.05 * yr^2, 0.6))
    ## cat(paste("Cyclone Occurance:", cyc_occur, "\n"))
    cyc_intensity <- rbeta(1, 2, 1) |> round(2)
    ## cat(paste("Cyclone intensity:", cyc_intensity, "\n"))
    ## cyc_spatial <- spatial_grid_pts_df  |>
    lat_offset <- runif(1, 0, 5)
    cyc_spatial <- spde$mesh$loc[, 1:2] |>
      as.data.frame() |>
      dplyr::select(Longitude = V1, Latitude = V2) |>
      dplyr::mutate(
        clong = as.vector(scale(Longitude, scale = FALSE)),
        clat = as.vector(scale(Latitude, scale = FALSE)),
        Y = lat_offset + runif(1, -1, 1) * clong + runif(1, -1, 1) *
          clat + sin(clat),
        # Y= Y - runif(1,-10,10),
        Y = abs(Y),
        Y = ifelse(Y > cyc_intensity, cyc_intensity, Y),
        Y = cyc_intensity - Y,
        Value = Y * cyc_occur
      )
    cyc[[yr]] <- cyc_spatial |>
      dplyr::mutate(Year = yr)
  }
  cyc <- do.call("rbind", cyc)
  cyc_effects_df <- cyc |>
    dplyr::mutate(Value = scales::rescale(Value, to = c(0, 1)))

  cyc_effects <- cyc_effects_df |>
    dplyr::select(-clong, -clat, -Y) |>
    tidyr::pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = Value
    ) |>
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

  cyc_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    cyc_effects
  )

  cyc_pts_effects <- cyc_pts_sample |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(Year = as.numeric(Year))
  list(
    cyc_effects = cyc_effects,
    cyc_effects_df = cyc_effects_df,
    cyc_pts_sample = cyc_pts_sample,
    cyc_pts_effects = cyc_pts_effects
  )
}

disturbance_other <- function(spatial_grid, spde, config) {
  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("mesh", "spde", "Q", "A")),
    sort(names(spde))
  )
  testthat::expect(
    inherits(spde$spde, c("inla.spde")),
    "spde$spde must be a inla.spde object"
  )
  testthat::expect_in(
    sort(c("years", "seed")),
    sort(names(config))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)
  set.seed(config$seed + 1)
  other_sample <- INLA::inla.qsample(length(config$years),
    spde$Q,
    seed = config$seed + 1,
    constr = spde$spde$f$extraconstr
  ) |>
    suppressMessages() |>
    suppressWarnings()

  rho <- rep(0.7, length(config$years))
  rho <- rbeta(length(config$years), 0.2, 1)
  x <- other_sample
  for (j in 2:length(config$years)) {
    x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * other_sample[, j]
  }
  other_effects <- scales::rescale(x, to = c(0, 1))
  other_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    other_effects
  )

  other_pts_effects <- other_pts_sample |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(Year = as.numeric(Year)) # ,
  ## Value=scales::rescale(Value, to=c(0,1)))

  list(
    other_effects = other_effects,
    other_pts_sample = other_pts_sample,
    other_pts_effects = other_pts_effects
  )
}

#################################
##################################
############ Information about the depth(s) 
##################################
depth_info <- function(x) {
if (x > 1){
  Depth_info = seq(3, 10, length=x)
}else{
  Depth_info = 10
}
return(Depth_info)
}

##################################
############ Select nth element from a vector 
##################################

nth_element <- function(vector, starting_position, n) { 
  vector[seq(starting_position, length(vector), n)] 
  }

##################################
############ Extract covariates 
##################################
extract_cov <- function(predictive_layer, cov_name) {

intersectlist <- st_intersects(predictive_layer, cov_name %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(4326)))
  
datlist <- list()
for (i in 1:nrow(pred_layer)){
  intersect <-  intersectlist[[i]]
  datlist[[i]] <-  cov_name[intersect,] %>%
       group_by(Year) %>%
       summarize(mean_value = mean(Value)) %>%
       mutate(tier = i + 999) %>%
       mutate(fYEAR = Year + min(dat$fYEAR) - 1)
  }
  
dat_all <- do.call(rbind, datlist) %>%
    merge(pred_layer) %>%
    st_as_sf() %>%
    dplyr::select(tier, fYEAR, mean_value,geometry)

return(dat_all)
}

##################################
############ GAM model 
##################################

# Picks number of spatial and temporal knots from mgcv models 

pick_knots_mgcv <- function(dat) {

# Degree of freedom = number of unique location x number of years 
 df <- dat %>%
    group_by(LONGITUDE, LATITUDE) %>%
    summarize()
  
  max_kspat <- ceiling(nrow(df) * .75)
  max_ktemp <- length(unique(dat$fYEAR))
  
  kspat <- seq(5, max_kspat, by = 10) # minium of 5 knots on the spatial dimension
  ktemp <- seq(4, max_ktemp, by = 5)  # minimum of 2 knots of the temporal dimension
  
  knbre<- expand.grid(kspat,ktemp)
  
  mod_list <- list()
  
  for ( i in 1 : nrow(knbre)){
    mod0 <- mgcv::gam(COUNT/TOTAL ~ te(LONGITUDE,LATITUDE, fYEAR, # inputs over which to smooth
                                       bs = c("tp", "cr"), # types of bases
                                       k=c(knbre[i,1],knbre[i,2]), # knot count in each dimension
                                       d=c(2,1)), # (s,t) basis dimension
                      data = dat,
                      control =  gam.control(scalePenalty = FALSE),
                      method = "GCV.Cp", family = binomial("logit"),
                      weights = TOTAL)
    
    mod_list[[i]] <- cbind(as.numeric(summary(mod0)$r.sq), as.numeric(summary(mod0)$s.table[[1]]), as.numeric(summary(mod0)$sp.criterion),
                           as.numeric(AIC(mod0)), knbre[i,1],knbre[i,2])
  }
  
  table_result <- do.call(rbind, mod_list) 
  colnames(table_result) <- c("Rsquared", "edf", "GCV", "AIC", "kspat", "ktemp")
  
  ## Criteria 
  # edf cannot be greater than degree of freedom 
  ## lowest GCV
  ## highest r2
  ## lowest AIC
  
  table_result <- table_result %>%
    data.frame() %>%
    arrange(desc(Rsquared), GCV, desc(AIC))
  
  return(table_result)
}

##################################
############ INLA model 
##################################

inla_prep <- function(dat, hexpred){
  
   ## ---- meshINLADataGrid
full.grid <- hexpred %>%
              st_centroid() %>%
              mutate(Longitude = st_coordinates(.)[,1],
                     Latitude = st_coordinates(.)[,2]) %>%
              st_drop_geometry()
  
full.coords <- full.grid %>%
              dplyr::select(Longitude,Latitude) %>%
              distinct()

max.edge =diff(range(full.coords[,1]))/15
bound.outer = diff(range(full.coords[,1]))/15
pol =   st_bbox(hexpred) %>%  st_as_sfc() 

mesh = inla.mesh.2d(loc.domain = st_coordinates(pol)[,1:2],
                     n=4,
                     max.edge = c(1,2)*max.edge,
                     offset=c(max.edge, bound.outer),
                     cutoff = max.edge/5)

 class(mesh) <- "inla.mesh"

spde <- inla.spde2.matern(mesh, alpha = 2) #, prior.range = c(.02, .01), prior.sigma = c(3, 0.01)
  i.spatial <- inla.spde.make.index('spatial.field',
                                    n.spde = spde$n.spde,
                                    n.group = length(unique(full.grid$fYEAR)))
# Data grid
hexpred_unique <- hexpred %>%
                group_by(tier) %>%
                filter(row_number() == 1) %>%
                dplyr::select(tier) 

data.grid_sf <- dat %>%
               filter(tier_cat == "tier_train") %>%
               left_join(hexpred_unique) %>% 
              mutate(
                P_CODE = factor(P_CODE),
                SITE_NO = factor(interaction(Reef, SITE_NO)),
                TRANSECT_NO = factor(interaction(SITE_NO, TRANSECT_NO)),
                Reef = factor(interaction(tier, Reef))) %>%
              dplyr::select(P_CODE, tier, fYEAR, fDEPTH, Reef, SITE_NO, TRANSECT_NO, DHW, CYC, OT, COUNT, TOTAL, geometry) %>%
              distinct() %>%
              st_as_sf()

# Filter NAs in covariate because of issues with edging (to be fixed later) 

data.grid_sf <- data.grid_sf %>% 
filter(is.na(st_dimension(.)) == FALSE )

coords <- data.grid_sf %>%
             st_centroid() %>%
             st_coordinates() %>%
              `[`(,c('X','Y'))

data.grid_sf$fYEAR <- as.factor(data.grid_sf$fYEAR)
  
data.grid <- data.grid_sf %>%
              mutate(fYEAR = factor(fYEAR, levels = rev(sort(levels(fYEAR))))) %>%
              st_drop_geometry() %>%
              data.frame()

 A.est <- inla.spde.make.A(mesh = mesh,
                            loc = coords,
                            group = as.numeric(data.grid$fYEAR))
  
  covariate <- data.grid %>%
    dplyr::select(fYEAR, Reef, DHW, CYC, OT, SITE_NO, TRANSECT_NO) %>%
    st_drop_geometry()
  
  covariate$Reef <- as.numeric(as.factor(covariate$Reef))
  covariate$DHW <- scale(covariate$DHW, scale = FALSE)
  covariate$CYC <- scale(covariate$CYC, scale = FALSE)
  covariate$OT <- scale(covariate$OT, scale = FALSE)

# Testing objects

data.grid_p_sf <- dat %>% 
              group_by(tier, fYEAR) %>%
              mutate(TRUE_COUNT_tier = mean(TRUE_COUNT),
              TOTAL_tier = mean(TOTAL)) %>%
              filter(row_number() == 1) %>%
              dplyr::select(tier, fYEAR, TRUE_COUNT_tier, TOTAL_tier, tier_cat, Reef, DHW, CYC, OT) %>%
              filter(tier_cat == "tier_test") %>%
              left_join(hexpred_unique) %>%
              ungroup() %>%
              st_as_sf()

# Filter NAs in covariate because of issues with edging (to be fixed later) 

data.grid_p_sf <- data.grid_p_sf %>% 
filter(is.na(st_dimension(.)) == FALSE )

coords_p <- data.grid_p_sf %>%
             st_centroid() %>%
             st_coordinates() %>%
              `[`(,c('X','Y'))

data.grid_p_sf$fYEAR <- as.factor(data.grid_p_sf$fYEAR)
  
data.grid_p <- data.grid_p_sf %>%
              mutate(fYEAR = factor(fYEAR, levels = rev(sort(levels(fYEAR))))) %>%
              mutate(TOTAL_tier = 500) %>%
              st_drop_geometry() %>%
              data.frame()

A.pred <- inla.spde.make.A(mesh = mesh,
                            loc = coords_p,
                            group = as.numeric(data.grid_p$fYEAR))

covariate_p <- data.grid_p %>%
  dplyr::select(fYEAR, Reef, DHW, CYC, OT) %>%
  st_drop_geometry()

covariate_p$Reef <- as.numeric(as.factor(covariate_p$Reef))
covariate_p$DHW <- scale(covariate_p$DHW, scale = FALSE)
covariate_p$CYC <- scale(covariate_p$CYC, scale = FALSE)
covariate_p$OT <- scale(covariate_p$OT, scale = FALSE)

# Stacks 

rprior<- list(theta = list(prior = "pccor1", param = c(0,0.9)))

stack.est <- inla.stack(data=list(y = data.grid$COUNT,
                                    Total=data.grid$TOTAL),
                          A=list(A.est, 1, 1, 1, 1, 1),
                          effects=list(c(i.spatial, list(b0=1)),
                          list(fYEAR = covariate$fYEAR),
                          list(Reef = covariate$Reef),
                          list(DHW = covariate$DHW),
                          list(CYC = covariate$CYC),
                          list(OT = covariate$OT)),
                          tag = 'est')


stack.pred <- inla.stack(data=list(y = NA,
                                  Total= data.grid_p$TOTAL_tier),
                        A=list(A.pred, 1, 1, 1, 1, 1),
                        effects=list(c(i.spatial, list(b0=1)),
                                     list(fYEAR = covariate_p$fYEAR),
                                     list(Reef = covariate_p$Reef),
                                     list(DHW = covariate_p$DHW),
                                     list(CYC = covariate_p$CYC),
                                     list(OT = covariate_p$OT)),
                        tag = 'pred')

stack.full <- inla.stack(stack.est, stack.pred)
                      
  if(formula == "cov"){
    form <- y ~ -1 + b0 + fYEAR + DHW + CYC + OT +
    f(spatial.field, model = spde, group = spatial.field.group, 
    control.group = list(model = "ar1", hyper=rprior)) + 
    f(Reef, model = "iid") 
  }else{
    form <- y ~ -1 + b0 + fYEAR +
    f(spatial.field, model = spde, group = spatial.field.group, 
    control.group = list(model = "ar1", hyper=rprior)) + 
    f(Reef, model = "iid") 
  }

  obj_inla <- list("stack.est" = stack.est, "stack.pred" = stack.pred, "spde" = spde, "form" = form, "Total" = dat$TOTAL, "mesh" = mesh,
                   "i.spatial" = i.spatial, "formula" = formula, "stack.full" = stack.full)
  
  return(obj_inla)
}

##################################
############ FRK model 
##################################

frk_prep <- function(dat){
  
  ## Construct STIDF object from data
  dat$Year <- as.Date(paste0(as.character(dat$fYEAR),"-01-01")) 
  dat$k_Z <- dat$TOTAL                                         
  lon_idx <- which(names(dat) == "LONGITUDE")                  
  lat_idx <- which(names(dat) == "LATITUDE")
  STObj <- stConstruct(x = dat,                               
                       space = c(lon_idx, lat_idx), 
                       time = "Year",                      
                       interval = TRUE)     
  
  ## Predictive layer
  HexPred_sp <- as_Spatial(hexpred)                                   
  nHEX <- nrow(subset(HexPred_sp, fYEAR == unique(dat$fYEAR)[1]))       
  nYEAR <- length(unique(HexPred_sp@data$fYEAR))            
  
  HexPred_sp@data$n_spat <- rep(1:nHEX, each = nYEAR)   
  BAUs_spat <- subset(HexPred_sp, fYEAR == unique(dat$fYEAR)[1])        
  coordnames(BAUs_spat) <- c("LONGITUDE", "LATITUDE")
  
  nrow(BAUs_spat@data)
  
  ## Construct spatio-temporal BAUs (will not contain covariate information for now)
  ST_BAUs <- auto_BAUs(manifold = STplane(),
                       data = STObj,
                       spatial_BAUs = BAUs_spat,
                       tunit = "years")
  
  ST_BAUs <- ST_BAUs[, 1:nYEAR, 1:2]                 
  ST_BAUs$fYEAR <- as.character(ST_BAUs$t + unique(dat$fYEAR)[1]-1)    
  ST_BAUs$n_spat <- rep(1:nHEX, nYEAR)              
  
  HexPred_sp@data$fYEAR <- as.character(HexPred_sp@data$fYEAR) 
  HexPred_sp@data$tier <- as.factor(HexPred_sp@data$tier) 
  
  ST_BAUs@data <- left_join(ST_BAUs@data, HexPred_sp@data , by = c("fYEAR","n_spat")) 
  
  ST_BAUs$fs <- 1                  
  ST_BAUs@sp@proj4string  <- CRS()  
  
  head(ST_BAUs@data)
  head(HexPred_sp@data)
  
  ## Covariates must only be in BAUs, so remove covariates associated with data
  overlapping_fields <- intersect(names(STObj@data), names(ST_BAUs@data))
  STObj@data[,overlapping_fields] <- NULL
  
  ## Create basis functions
  basis <- auto_basis(STplane(),
                      STObj,
                      tunit = "years",
                      #nres = 2L, # for development
                      nres = 3L, # for final run
                      regular = TRUE)

  
  obj_frk <- list("ST_BAUs" = ST_BAUs, "STObj" = STObj, "basis" = basis)
  return(obj_frk)
  
}

##################################
############ ML model 
##################################

ml_tuning <- function(dat_cov, vfold){

start_time <- Sys.time()
# Split dataset 
## Split into training and testing data

data_split <- initial_split(dat_cov, prop = .999)
data_train <- training(data_split)
data_test <- testing(data_split)

# 2. Hyperparameters tuning

## 2.1 Define the recipe

boosted_recipe <- recipe(COUNT ~ ., data = data_train) %>% 
  step_dummy(all_nominal_predictors())

## 2.2 Define the model

boosted_model <- boost_tree(learn_rate = tune(),
                            trees = tune(), 
                            min_n = tune(), 
                            tree_depth = tune()) %>% # Model type
  set_engine("xgboost") %>% # Model engine
  set_mode("regression") # Model mode

## 2.3 Define the workflow

boosted_workflow <- workflow() %>%
  add_recipe(boosted_recipe) %>% 
  add_model(boosted_model)

## 2.4 Create the grid - plus long quand size parametre est grand

tune_grid <- grid_space_filling(learn_rate(),
                                trees(),
                                tree_depth(),
                                min_n(),
                                size = 30) # 5 for dev

## 2.5 Run the hyperparameters tuning - plus long en presence des covariates

tuned_results <- tune_grid(boosted_workflow,
                           resamples = vfold_cv(data_train, v = vfold), 
                           grid = tune_grid)

## 2.6 Get best set of parameters

model_hyperparams <- select_best(tuned_results, metric = "rmse") %>% 
  select(-".config") %>% 
  as_tibble(.) %>%
  mutate(nb_training = nrow(data_train),
         nb_testing = nrow(data_test)) 

# 3. Predicted vs observed

## 3.1 Redefine the model (with hyperparameters)

boosted_model <- boost_tree(learn_rate = model_hyperparams$learn_rate,
                            trees = model_hyperparams$trees, 
                            min_n = model_hyperparams$min_n, 
                            tree_depth = model_hyperparams$tree_depth) %>% # Model type
  set_engine("xgboost") %>% # Model engine
  set_mode("regression") # Model mode

## 3.2 Redefine the workflow

boosted_workflow <- workflow() %>%
  add_recipe(boosted_recipe) %>% 
  add_model(boosted_model)

## 3.3 Fit the final model

final_model <- boosted_workflow %>%
  last_fit(data_split)

final_fitted <- final_model$.workflow[[1]]

## 3.4 Model performance

model_performance <- collect_metrics(final_model) %>% 
  select(-".estimator", -".config") %>% 
  pivot_wider(names_from = ".metric", values_from = ".estimate") 

## 3.4 Predicted vs Observed

result_pred_obs <- data_test %>% 
  mutate(yhat = predict(final_fitted, data_test)$.pred) %>% 
  rename(y = COUNT) %>% 
  select( y, yhat)

p <- ggplot(result_pred_obs, aes(x = y, y = yhat)) +
  geom_point(alpha = .4) + 
  geom_abline(col = "red", linetype = "dashed") + 
  ggtitle(model_name) +
  theme_bw() +
  ylim(0, max(result_pred_obs$y+5)) +
  xlim(0, max(result_pred_obs$y+5)) +
  xlab("Observed values") +
  ylab("Predicted values") 

figure_name <- paste(model_name, "fit.png", sep = "_")
figure_path <-paste0(getwd(),"/",title_of_run,"/figures/")

# 3.5 Save the plot
ggsave(p, file = paste(figure_path, figure_name, sep="/"))

# 4. Results
tuning_result <- lst(model_hyperparams,
                     model_performance,
                     result_pred_obs)

end_time <- Sys.time()
computing_time <- end_time - start_time

model_tuning <- list(name = model_name, model_fit = tuning_result, time = computing_time)
return(model_tuning)
}

##################################
############ model predictions 
##################################

predictions_INLA <- function(model.out, n.sim){
  
draws <- inla.posterior.sample(n.sim, result=model.out, seed=123) 
  
hexpred_dat <- hexpred %>%
    st_centroid() %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2]) %>%
    dplyr::select(fYEAR, Longitude, Latitude, Reef, tier, DHW, CYC, OT) %>%
    st_drop_geometry()%>% 
    mutate(across(fYEAR, as.factor)) 
  
  full.coords <- hexpred_dat %>%
    dplyr::select(Longitude,Latitude) %>%
    distinct()

  proj.grid <- inla.mesh.projector(obj_inla$mesh, loc=as.matrix(full.coords))
  cellmeans = sapply(draws, function(x) x[['latent']])
  
  if(obj_inla$formula == "cov"){
  i.mod <- sapply(c('APredictor','^Predictor','spatial.field','Reef','fYEAR[0-9]*[:][0-9]*$', 'DHW', 'CYC', 'OT', 'b0'),
                  function(x) grep(x, draws[[1]]$latent %>% rownames))
  }else{
  i.mod <- sapply(c('APredictor','^Predictor','spatial.field','Reef','fYEAR[0-9]*[:][0-9]*$','b0'),
                  function(x) grep(x, draws[[1]]$latent %>% rownames))
  }

  # retrieve the spatial.fields posteriors
  
  cellmeans.full <- cellmeans[i.mod[[3]],] %>%         
    as.data.frame %>%                                 
    mutate(fYEAR = rep(as.numeric(levels(obj_inla$stack.est$effects$data$fYEAR)),
                       each = which(obj_inla$i.spatial$spatial.field.group == 1) %>% length)) %>%
    group_by(fYEAR) %>%
    tidyr::nest() %>% 
    mutate(Spatial = map(.x = data,
                         .f = function(x)
                           as.matrix(inla.mesh.project(proj.grid, x)))) %>%
    mutate(geometry = list(hexpred_dat %>%
                             dplyr::select(tier, Longitude, Latitude) %>%
                             distinct())) 
  
  cellmeans.full$fYEAR <- as.factor(cellmeans.full$fYEAR)
  
  # retrieve the fixed effects 
  
  if(obj_inla$formula == "cov"){
  Xmat <- cbind(1, model.matrix(~ -1 + fYEAR + DHW + CYC + OT, data=hexpred_dat))
  wch <- grep(paste0(c("b0","fYEAR", "DHW", "CYC", "OT"), collapse="|"), names(i.mod))
  ii = unlist(i.mod[wch])
  }else{
  Xmat <- cbind(1, model.matrix(~ -1 + fYEAR, data=hexpred_dat)) 
  wch <- grep(paste0(c("b0","fYEAR"), collapse="|"), names(i.mod))
  ii = unlist(i.mod[wch])}

 cellmeans.full.1 <- t(cellmeans[ii,]) %*% t(Xmat)
    
  cellmeans.fixed <- hexpred_dat %>%
    dplyr::select(Longitude, Latitude, fYEAR) %>%
    cbind(V = t(cellmeans.full.1)) %>%
    as.data.frame %>%                           
    dplyr::select(starts_with("V"))%>%
    mutate(fYEAR = rep(unique(hexpred_dat$fYEAR),each=nrow(full.coords)))%>%
    group_by(fYEAR) %>%
    tidyr::nest() 

  # Add the posteriors together
  cellmeans.full.c <-
    cellmeans.full %>%
    full_join(cellmeans.fixed %>%
                rename(data1 = data)) %>%
    mutate(value = map2(.x = Spatial, .y = data1,
                        .f = function(.x, .y) as.data.frame(.x + .y))) %>%
    dplyr::select(fYEAR, geometry, value) %>%
    tidyr::unnest(cols = c(geometry, value)) %>%
    tidyr::pivot_longer(c = starts_with("V"), names_to = "Rep") %>%
    mutate(Rep = gsub('\\.','',Rep)) %>%
    ungroup()
  
  cellmeans.full.cc <- cellmeans.full.c %>%
    mutate(pred = plogis(value))
  
  # ###################################### Sum across tiers 
  
  pred_sum_sf <- cellmeans.full.cc %>%
    group_by(fYEAR,tier) %>% 
    median_hdci(pred) %>%
    inner_join(hexpred %>% group_by(tier) %>% slice(1) %>% dplyr::select(geometry,tier)) %>% 
    st_as_sf(sf_column_name = "geometry") %>%
    mutate(Unc = .upper - .lower) %>%
    mutate(tier_fYEAR = paste0(tier,fYEAR)) %>%
    dplyr::select(fYEAR, tier, pred, .lower, .upper, Unc, tier_fYEAR)

  return(pred_sum_sf)
 
}

predictions_INLA_stacked <- function(model.out, hexpred, dat, index){

hexpred_unique <- hexpred %>%
                group_by(tier) %>%
                filter(row_number() == 1) %>%
                dplyr::select(tier) 

data.grid_p_sf <- dat %>% 
              group_by(tier, fYEAR) %>%
              filter(row_number() == 1) %>%
              filter(tier_cat == "tier_test")  %>%
              left_join(hexpred_unique) %>%
              ungroup() %>%
              st_as_sf()

# Filter NAs in covariate because of issues with edging (to be fixed later) 

data.grid_p <- data.grid_p_sf %>% 
filter(is.na(st_dimension(.)) == FALSE )

index_pred <- inla.stack.index(stack = obj_inla$stack.full, tag = "pred")$data

data.grid_p$pred <- model.out$summary.fitted.values[index_pred,"mean"]
data.grid_p$.lower <- model.out$summary.fitted.values[index_pred,"0.025quant"]
data.grid_p$.upper <- model.out$summary.fitted.values[index_pred,"0.975quant"]

pred_sum_sf <- data.grid_p %>% 
    st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = st_crs(4326)) %>%
    mutate(Unc = .upper - .lower) %>%
    mutate(tier_fYEAR = paste0(tier,fYEAR)) %>%
    dplyr::select(fYEAR, tier, pred, .lower, .upper, Unc, tier_fYEAR) 

  return(pred_sum_sf)
  
}

predictions_FRK <- function(model.out){
  
  pred <- predict(model.out, type = c("mean"))
  
  # Extracting posterior distributions of predictive locations 
  
  post_dist_df <- as.data.frame(pred$MC$mu_samples) %>% 
    mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(tier = obj_frk$ST_BAUs@data$tier) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR,tier,id_loc),
                        names_to = "draw", 
                        values_to = "pred"
    )
  
  # Summary predictions at tier level
  hexpred$tier <- as.factor(hexpred$tier)
  
  pred_sum_sf <- post_dist_df %>% group_by(fYEAR,tier) %>% 
    median_hdci(pred)%>%
    inner_join(hexpred %>% group_by(tier) %>% 
    summarize() %>% dplyr::select(geometry,tier)) %>% 
    st_as_sf(sf_column_name = "geometry") %>%
    mutate(Unc = .upper - .lower) %>%
    mutate(tier_fYEAR = paste0(tier,fYEAR)) %>%
    dplyr::select(fYEAR, tier, pred, .lower, .upper, Unc, tier_fYEAR)

  return(pred_sum_sf)
}

predictions_brms <- function(model.out, n.sim){
  
  hexpred_dat <- hexpred %>%
    st_centroid() %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2]) %>%
    dplyr::select(fYEAR, Longitude, Latitude, DHW, CYC, OT, Reef, tier) %>%
    st_drop_geometry()%>% 
    mutate(TOTAL = round(mean(dat$TOTAL),0)) %>%
    rename(LONGITUDE = Longitude) %>%
    rename(LATITUDE = Latitude)
  
  
  pred <- predict(model.out, hexpred_dat, allow_new_levels = TRUE, ndraws = n.sim,
                  incl_autocor = TRUE, nug =  1e-12) %>%
    data.frame() %>%
    mutate(across(everything()), . / unique(hexpred_dat$TOTAL)) 
 
  # Summary predictions at tier level
  hexpred$tier <- as.factor(hexpred$tier)
  
  pred_sum_sf <- pred %>%
    data.frame() %>%
    cbind(hexpred) %>% 
    st_as_sf(sf_column_name = "geometry") %>%
    mutate(Unc = Q97.5 - Q2.5) %>%
    mutate(tier_fYEAR = paste0(tier,fYEAR)) %>%
    rename(pred = Estimate) %>%
    rename(.upper = Q97.5) %>%
    rename(.lower = Q2.5) %>%
    mutate(tier_fYEAR = paste0(tier,fYEAR)) %>%
    dplyr::select(fYEAR, tier, pred, .lower, .upper, Unc, tier_fYEAR)

  return(pred_sum_sf)
  
}

##################################
############ model predictions - broad spatial scale
##################################

predictions_INLA_broad <- function(model.out, n.sim){
  
draws <- inla.posterior.sample(n.sim, result=model.out, seed=123) 
  
hexpred_dat <- hexpred %>%
    st_centroid() %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2]) %>%
    dplyr::select(fYEAR, Longitude, Latitude, Reef, tier) %>%
    st_drop_geometry()%>% 
    mutate(across(fYEAR, as.factor)) 
  
  full.coords <- hexpred_dat %>%
    dplyr::select(Longitude,Latitude) %>%
    distinct()
   
  proj.grid <- inla.mesh.projector(obj_inla$mesh, loc=as.matrix(full.coords))
  cellmeans = sapply(draws, function(x) x[['latent']])
  
  i.mod <- sapply(c('APredictor','^Predictor','spatial.field','Reef','fYEAR[0-9]*[:][0-9]*$'),
                  function(x) grep(x, draws[[1]]$latent %>% rownames))
  
  # retrieve the spatial.fields posteriors
  
  cellmeans.full <- cellmeans[i.mod[[3]],] %>%         
    as.data.frame %>%                                 
    mutate(fYEAR = rep(as.numeric(levels(obj_inla$stack.est$effects$data$fYEAR)),
                       each = which(obj_inla$i.spatial$spatial.field.group == 1) %>% length)) %>%
    group_by(fYEAR) %>%
    tidyr::nest() %>% 
    mutate(Spatial = map(.x = data,
                         .f = function(x)
                           as.matrix(inla.mesh.project(proj.grid, x)))) %>%
    mutate(geometry = list(hexpred_dat %>%
                             dplyr::select(tier, Longitude, Latitude) %>%
                             distinct())) 
  
  cellmeans.full$fYEAR <- as.factor(cellmeans.full$fYEAR)
  
  # retrieve the fixed effects 
  
  Xmat <- model.matrix(reformulate(c("0 + fYEAR")), data=hexpred_dat)
  
  #wch <- c(6,7,8)
  #ii = unlist(i.mod[wch])
  
  ii = unlist(i.mod[5])
  
  cellmeans.full.1 <- t(cellmeans[ii,]) %*% t(Xmat)
  
  cellmeans.fixed <- hexpred_dat %>%
    dplyr::select(Longitude, Latitude, fYEAR) %>%
    cbind(V = t(cellmeans.full.1)) %>%
    as.data.frame %>%                           
    dplyr::select(starts_with("V"))%>%
    mutate(fYEAR = rep(unique(hexpred_dat$fYEAR),each=nrow(full.coords)))%>%
    group_by(fYEAR) %>%
    tidyr::nest() 
  
  # Add the posteriors together
  cellmeans.full.c <-
    cellmeans.full %>%
    full_join(cellmeans.fixed %>%
                rename(data1 = data)) %>%
    mutate(value = map2(.x = Spatial, .y = data1,
                        .f = function(.x, .y) as.data.frame(.x + .y))) %>%
    dplyr::select(fYEAR, geometry, value) %>%
    tidyr::unnest(cols = c(geometry, value)) %>%
    tidyr::pivot_longer(c = starts_with("V"), names_to = "Rep") %>%
    mutate(Rep = gsub('\\.','',Rep)) %>%
    ungroup()
  
  cellmeans.full.cc <- cellmeans.full.c %>%
    mutate(pred = plogis(value))
  
  # ###################################### Sum across tiers 
  
  pred_sum_sf <- cellmeans.full.cc %>%
    group_by(fYEAR) %>% 
    median_hdci(pred) 

  return(pred_sum_sf)
 
}

predictions_INLA_broad_stacked <- function(model.out, dat, hexpred){

hexpred_unique <- hexpred %>%
                group_by(tier) %>%
                filter(row_number() == 1) %>%
                dplyr::select(tier) 

data.grid_p <- dat %>% 
              group_by(tier, fYEAR) %>%
              filter(row_number() == 1) %>%
              filter(tier_cat == "tier_test")  %>%
              left_join(hexpred_unique) %>%
              ungroup() %>%
              st_as_sf()

index_pred <- inla.stack.index(stack = obj_inla$stack.full, tag = "pred")$data

data.grid_p$pred <- model.out$summary.fitted.values[index_pred,"mean"]
data.grid_p$.lower <- model.out$summary.fitted.values[index_pred,"0.025quant"]
data.grid_p$.upper <- model.out$summary.fitted.values[index_pred,"0.975quant"]

# Filter NAs in covariate because of issues with edging (to be fixed later) 

data.grid_p <- data.grid_p %>% 
filter(is.na(st_dimension(.)) == FALSE )

pred_sum_sf <- data.grid_p %>%
    mutate(Unc = .upper - .lower) %>%
    mutate(tier_fYEAR = paste0(tier,fYEAR)) %>%
    dplyr::select(fYEAR, tier, pred, .lower, .upper, Unc, tier_fYEAR) %>%
    data.frame() %>%
    group_by(fYEAR) %>% 
    median_hdci(pred) 

  return(pred_sum_sf)

}

predictions_FRK_broad <- function(model.out){
  
  pred <- predict(model.out, type = c("mean"))
  
  # Extracting posterior distributions of predictive locations 
  
  post_dist_df <- as.data.frame(pred$MC$mu_samples) %>% 
    mutate(fYEAR = obj_frk$ST_BAUs@data$fYEAR) %>%
    mutate(tier = obj_frk$ST_BAUs@data$tier) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR,tier,id_loc),
                        names_to = "draw", 
                        values_to = "pred"
    )
  
  # Summary predictions at tier level
  hexpred$tier <- as.factor(hexpred$tier)
  
  pred_sum_sf <- post_dist_df %>% group_by(fYEAR) %>% 
    median_hdci(pred)

  return(pred_sum_sf)
}

predictions_brms_broad <- function(model.out, n.sim){
  
  hexpred_dat <- hexpred %>%
    st_centroid() %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2]) %>%
    dplyr::select(fYEAR, Longitude, Latitude, DHW, CYC, OT, Reef, tier) %>%
    st_drop_geometry()%>% 
    mutate(TOTAL = round(mean(dat$TOTAL),0)) %>%
    rename(LONGITUDE = Longitude) %>%
    rename(LATITUDE = Latitude)
  
  
  pred <- predict(model.out, hexpred_dat, allow_new_levels = TRUE, ndraws = n.sim,
                  incl_autocor = TRUE) %>%
    data.frame() %>%
    mutate(across(everything()), . / unique(hexpred_dat$TOTAL)) 
 
  # Summary predictions at tier level
  hexpred$tier <- as.factor(hexpred$tier)
  
  pred_sum_sf <- pred %>%
    data.frame() %>%
    cbind(hexpred) %>% 
    st_as_sf(sf_column_name = "geometry") %>%
    mutate(Unc = Q97.5 - Q2.5) %>%
    mutate(tier_fYEAR = paste0(tier,fYEAR)) %>%
    rename(pred = Estimate) %>%
    rename(.upper = Q97.5) %>%
    rename(.lower = Q2.5) %>% group_by(fYEAR) %>% 
    median_hdci(pred)

  return(pred_sum_sf)
  
}

model_bootstrap_run <- function(bootstrap_i, pdp){
  
  # Data preparation
  
    ## Sample with replacement (TO CHECK)
  
  data_split <- slice_sample(dat_cov, n = nrow(dat_cov), replace = TRUE)
  data_split <- initial_split(data_split, prop = .99)
  data_train <- training(data_split)
  data_test <- testing(data_split)
  
  # Fit the model
  
  ## Define the recipe
  
  boosted_recipe <- recipe(COUNT ~ ., data = data_train) |> 
    step_dummy(all_nominal_predictors())
  
  ## Define the model
  
  boosted_model <- boost_tree(learn_rate = model_hyperparams$learn_rate,
                              trees = model_hyperparams$trees, 
                              min_n = model_hyperparams$min_n, 
                              tree_depth = model_hyperparams$tree_depth) |> # Model type
    set_engine("xgboost") |> # Model engine
    set_mode("regression") # Model mode
  
  ## Define the workflow
  
  boosted_workflow <- workflow() |>
    add_recipe(boosted_recipe) |> 
    add_model(boosted_model)
  
  ## Fit the final model
  
  final_model <- boosted_workflow |>
    last_fit(data_split)
  
  # Model outputs
  
  ## Model performance
  
  model_performance <- collect_metrics(final_model) |> 
    select(-".estimator", -".config") |> 
    pivot_wider(names_from = ".metric", values_from = ".estimate") |>
    mutate(bootstrap = bootstrap_i)
  
  ## Variable importance
  
  result_vip <- final_model |> 
    extract_fit_parsnip() |> 
    vip(num_features = 100) %>% 
    .$data %>% 
    rename(predictor = 1, importance = 2) |> 
    mutate(bootstrap = bootstrap_i)
  
  ## Partial Dependence Plots
  
  final_fitted <- final_model$.workflow[[1]]
  pdp = FALSE
  
  # if(pdp == TRUE){
  #   
  #   model_explain <- explain_tidymodels(model = final_fitted, 
  #                                       data = dplyr::select(data_train), 
  #                                       y = data_train$COUNT)
  #   
  #   result_pdp <- model_profile(explainer = model_explain,
  #                               N = NULL, 
  #                               center = FALSE,
  #                               type = "partial",
  #                               variables = c("year", "decimalLongitude", "decimalLatitude", "nb_cyclones",
  #                                             "wind_speed_y5", "nb_cyclones_y5", "pred_elevation", "pred_reefextent",
  #                                             "pred_land", "pred_gravity", "pred_enso", "pred_sst_sd", "pred_sst_skewness",
  #                                             "pred_sst_max", "pred_sst_max_y1", "pred_sst_mean",
  #                                             "pred_sst_mean_y1", "pred_dhw_max", "pred_dhw_max_y1",
  #                                             "pred_chla_mean", "pred_chla_sd", "pred_population", "verbatimDepth"),
  #                               variable_splits_type = "uniform") |> 
  #     .$agr_profiles |> 
  #     as_tibble(.) |> 
  #     select(-"_label_", -"_ids_") |> 
  #     rename(x = "_x_", y_pred = "_yhat_", predictor = "_vname_") #|> 
  #     #mutate(fGROUP = fGROUP_i, bootstrap = bootstrap_i)
  #   
  # }
  # 
  # 4. Predictions
  
  ## 4.1 Predict values for new locations
  
  result_trends <-  hexpred_dat |> 
  #  dplyr::select(-c(P_CODE, fDEPTH, TRUE_COUNT, TOTAL, SITE_NO, TRANSECT_NO, tier, fGROUP)) |>
    mutate(pred = predict(final_fitted, new_data = hexpred_dat)$.pred) |>
    mutate(bootstrap = bootstrap_i) |>
    mutate(Name = model_name)
  
  # 5. Return the results

  return(list(model_performance = model_performance, result_vip = result_vip, result_trends = result_trends) )
  
}

predictions_ML <- function(model.out, dat){

# Extracting posterior distributions of predictive locations 
hexpred_join <- hexpred %>% 
    st_centroid() %>%
     mutate(LONGITUDE = st_coordinates(.)[,1],
           LATITUDE = st_coordinates(.)[,2]) %>%
  mutate(across(c(Reef, fYEAR), as_factor)) %>%
  mutate(TOTAL = round(mean(dat$TOTAL),0))

model.out <- as.data.frame(model.out) #%>%
#mutate(
#    across(
#      c(LONGITUDE, LATITUDE),
#      ~ case_when(
 #       is.factor(.) ~ as.numeric(levels(.))[.], 
#        TRUE ~ as.numeric(.)                        
#      )))

tt <- hexpred_join %>% 
dplyr::select(LONGITUDE, LATITUDE, tier) %>%
distinct() %>%
st_drop_geometry()

post_dist_df <- model.out %>% 
left_join(hexpred_join %>% dplyr::select(tier, LONGITUDE, LATITUDE), by = c("LONGITUDE", "LATITUDE"),
                                         relationship = "many-to-many")  

pred_sum_sf <- post_dist_df %>% group_by(fYEAR,tier) %>% 
    median_hdci(pred) %>%
    mutate(across(c(pred, .lower, .upper), ~ . / unique(hexpred_join$TOTAL))) %>%
   left_join(hexpred %>% dplyr::select(tier, geometry) %>% distinct())  %>% 
   st_as_sf(sf_column_name = "geometry", crs = 4326) %>%
    mutate(Unc = .upper - .lower) %>%
    mutate(tier_fYEAR = paste0(tier,fYEAR)) %>%
    dplyr::select(fYEAR, tier, pred, .lower, .upper, Unc, tier_fYEAR)

  return(pred_sum_sf)
}

predictions_ML_broad <- function(model.out, dat){

# Extracting posterior distributions of predictive locations 
hexpred_join <- hexpred %>% 
    st_centroid() %>%
     mutate(LONGITUDE = st_coordinates(.)[,1],
           LATITUDE = st_coordinates(.)[,2]) %>%
  mutate(across(c(Reef, fYEAR), as_factor)) %>%
  mutate(TOTAL = round(mean(dat$TOTAL),0))

model.out <- as.data.frame(model.out) #%>%
#mutate(
#    across(
 #     c(LONGITUDE, LATITUDE),
#      ~ case_when(
#        is.factor(.) ~ as.numeric(as.character(.)), 
#        TRUE ~ as.numeric(.)                        
#      )))

post_dist_df <- model.out %>% 
left_join(hexpred_join %>% dplyr::select(tier, LONGITUDE, LATITUDE), by = c("LONGITUDE", "LATITUDE"),
                                         relationship = "many-to-many")  

pred_sum_sf <- post_dist_df %>% group_by(fYEAR) %>% 
    median_hdci(pred) %>%
    mutate(across(c(pred, .lower, .upper), ~ . / unique(hexpred_join$TOTAL))) %>%
    mutate(Unc = .upper - .lower) %>%
    dplyr::select(fYEAR, pred, .lower, .upper, Unc)

  return(pred_sum_sf)
}

##################################
############ plotting 
##################################

plot_predictions <- function(dat){

pal_pred <- lacroix_palette("Pamplemousse", n = 10, type = "continuous")
  
  p_pred <- ggplot() + 
    geom_sf(data = dat, aes(fill = pred)) +
    facet_wrap(~fYEAR) +  scale_fill_gradientn(colours = pal_pred) + 
    theme_bw() + 
    xlab("longitude") +
    ylab("latitude")

return(p_pred)
}

plot_predictions_unc <- function(dat){
  
pal_unc<- wes_palette("Zissou1", 10, type = "continuous")

  p_unc <- ggplot() + 
    geom_sf(data = dat, aes(fill = Unc)) +
    facet_wrap(~fYEAR) +  scale_fill_gradientn(colours = pal_unc) + 
    theme_bw() + 
    xlab("longitude") +
    ylab("latitude")

return(p_unc)
}

plot_diff <- function(dat){
  
  p_diff <- ggplot() + 
    geom_sf(data = dat, aes(fill = Diff)) +
    geom_sf(data = dat %>% filter(tier_cat == "tier_train"), col = "black") +
    facet_wrap(~fYEAR) +  
    theme_bw() + 
    xlab("longitude") +
    ylab("latitude") + 
    scale_fill_continuous(
    type = "viridis",
    guide = "colorbar",
    na.value = "transparent",
    limits = c(-0.7, 0.7)
  )

return(p_diff)
}

plot_traj <- function(dat, dat_pred){
p_traj <- ggplot() + 
  geom_line(data = dat, aes(x = fYEAR, y = (COUNT/TOTAL)*100, group = interaction(as.factor(TRANSECT_NO), REEF_NAME)), 
            show.legend = FALSE, linewidth=.1, col="grey30") + 
  geom_ribbon(data = dat_pred, aes(x=fYEAR,ymin=.lower*100, ymax=.upper*100, group=1),alpha=.2, fill ="#0072B2FF") +
  geom_line(data = dat_pred, aes(x=fYEAR, y=pred*100, group=1),size=.4) +
  facet_wrap(~tier, ncol=3) +
  ylab("Cover") + xlab("Year")+theme_bw()+
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white')) +
  scale_x_discrete(breaks= nth_element(unique(dat$fYEAR),1,4))

return(p_traj)
}

plot_traj_broad <- function(dat, GRMF_all){
plot_traj_broad <- ggplot() +
  geom_ribbon(data = pred_sum_sf %>% data.frame(), aes(x = fYEAR, ymin=.lower, ymax=.upper, group=1), alpha=.2, fill="#83BFA9")+
  geom_line(data = pred_sum_sf %>% data.frame(), aes(x=fYEAR, y=pred,group=1), col="black", linewidth=1.1)+
  geom_point(data = pred_sum_sf %>% data.frame(), aes(x=fYEAR, y=pred), col="black", size=2.1)+
  geom_line(data = GRMF_all, aes(x=as.factor(fYEAR), y=true, group=1),linewidth = .4, col = "blue", linetype = "dashed") +
  xlab("Year") +ylab("Coral cover")+theme_bw()+
  ylim(0,1) + 
  theme(axis.text.x = element_text(size=13),legend.position = "none",
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15))+
  scale_x_discrete(breaks= nth_element(unique(dat$fYEAR),1,4))

return(plot_traj_broad)
}
##################################
############ predictive indicators  
##################################

## 95% coverage
coverage95 <- function(z, lower, upper) {
  abs(0.95 - (sum((z < upper) & (z > lower)) / length(z)))
}

## 95% interval score
IS95 <- function(true, lower, upper) {
  alpha = 0.05
  pred95l <- lower 
  pred95u <- upper
  ISs <- (pred95u - pred95l) + 2/alpha * (pred95l - true) * (true < pred95l) +
    2/alpha * (true - pred95u) * (true > pred95u)
  mean(ISs)
}

## Root-mean-squared prediction error
RMSPE <- function(z,pred) {
  Y <- (z - pred)^2
  sqrt(mean(Y))
}

crps <- function(obs, pred, ...)
  ## Tilmann Gneiting's crps code, assumes pred is either a vector of length
  ## 2 (mu, sig) or a matrix of mu and sig if each forcast is different
{
  if(is.null( dim(pred)) & length(pred)==2){mu <- pred[1];
  sigma <- pred[2]} else {
    mu<- as.numeric( pred[,1] ); sigma <- as.numeric( pred[,2]) }
  
  z <- (obs-mu)/sigma ## center and scale
  
  crps<- sigma * (z*(2*pnorm(z,0,1)-1) + 2*dnorm(z,0,1) - 1/sqrt(pi))
  ign <-  0.5*log(2*pi*sigma^2) + (obs - mu)^2/(2*sigma^2)
  pit <- pnorm(obs, mu,sigma )
  
  return(list(crps = crps, CRPS = mean(crps), ign = ign, IGN = mean(ign), pit = pit) )
  
}

#################################
############# automated report
#################################

## Grab png files associated with name_plot and add them together 

print_attr_plots_level1 <- function(list_folders, name_plot) {


  # List files matching the pattern in the specified path
  file_mod <- list.files(path, pattern = name_plot, full.names = TRUE)
  
  # Generate plots for each file
  p_mod <- lapply(file_mod, function(file_mod) {
    ggplot2::ggplot() +
      ggplot2::annotation_custom(grid::rasterGrob(png::readPNG(file_mod), interpolate = TRUE)) +
      ggplot2::theme_void() 
  })
  
  # Combine the plots using patchwork
  combined_mod <- patchwork::wrap_plots(p_mod, ncol = 2)
  return(combined_mod)
}

print_attr_plots_level2 <- function(list_folders, name_plot) {

  # List files matching the pattern in the specified path
  file <- list.files(paste0(path, list_folders, FIGURES_PATH),
    pattern = name_plot,
    full.names = TRUE
  )

  # Generate plots for each file
  p <- lapply(file, function(file) {
    ggplot2::ggplot() +
      ggplot2::annotation_custom(grid::rasterGrob(png::readPNG(file), interpolate = TRUE)) +
      ggplot2::theme_void()
  })
  
  # Combine the plots using patchwork
  combined_mod <- patchwork::wrap_plots(p, ncol = 2)
  return(combined_mod)
}


print_attr_plots_level3 <- function(list_folders, name_plot, subtitle = FALSE, plot_labels = FALSE, label = "", n_cols) {

  # Debug: print the list_folders argument
  cat("####", list_folders, "\n\n")
  
  # Define file paths
  file <- list.files(paste0(path, list_folders, FIGURES_PATH),
    pattern = name_plot,
    full.names = TRUE
  )
  
  # Generate plots for each file
  p <- lapply(file, function(file) {
    plot <- ggplot2::ggplot() +
      ggplot2::annotation_custom(grid::rasterGrob(png::readPNG(file), interpolate = TRUE)) +
      ggplot2::theme_void() + 
      ggplot2::theme(plot.subtitle = element_text(size = 18, hjust = 0.5) )
    
    # Add subtitle if subtitle = TRUE
    if (subtitle) {
      plot <- plot + ggplot2::labs(
        subtitle = sub(
          pattern = paste0("^.*",str_replace(name_plot, "\\[.*", ""), "(.*)\\.png$"), # Use name_plot dynamically
          replacement = "\\1",
          x = file
        )
      )
    }
   return(plot)
  })
  
  # Combine the plots using patchwork
  combined_mod <- patchwork::wrap_plots(p, ncol = n_cols)
  
  # Print plots with labels if plot_labels = TRUE
  if (plot_labels) {
    div_label <- paste0("#fig-pred-data", list_folders)
    cat("\n::: {", div_label, "}\n", sep = "")
    
    print(combined_mod) # Print combined plot
    cat(paste0("\n\n", label, " ", list_folders, "\n"), sep = "")
    cat("\n:::\n")
  }
  
  # Return combined plot
  return(combined_mod)
}

# Function to load an RData file and extract a list
load_and_extract <- function(file) {
  env <- new.env()
  load(file, envir = env)

# Function to find "years" and replace it with its length
  config_list <- map( env$config_list, ~ {
  if (!is.null(.x$years)) {
    .x$years <- length(.x$years)  # Replace "years" with its length
  }
  .x 
})

data_tables <- map(seq_along(config_list[1:4]), ~ {
  config_list[[.x]] %>%
    as_tibble() %>%
    mutate(simulation = map_chr(str_split(file, "/"), 2)) 
})

names(data_tables) <- c("config_sp", "config_large", "config_fine", "config_pt")
return(data_tables)
}


# # MATERN SPDE MODEL IN MGCV 
# # These functions define the Matern SPDE model as a basis-penalty smoother
# # Setup the SPDE Matern smooth in mgcv
# # See ?smooth.construct in mgcv for details of the input and output 
# # Special note: the xt argument in s(..., bs = "spde", xt = ...) can be used
# # to specify a mesh, if NULL, a mesh with regular knots is constructed. 
# smooth.construct.spde.smooth.spec <- function(object, data, knots){
#   # observation locations
#   dim <- length(object$term) 
#   if (dim > 2 | dim < 1) stop("SPDE Matern can only be fit in 1D or 2D.")
#   if (dim == 1) {
#     x <- data[[object$term]]
#   } else {
#     x <- matrix(0, nr = length(data[[1]]), nc = 2) 
#     x[,1] <- data[[object$term[1]]]
#     x[,2] <- data[[object$term[2]]]  
#   }
#   # setup mesh or accept user mesh
#   if (is.null(object$xt)) {
#     if (dim == 1) {
#       t <- seq(min(x), max(x), len=object$bs.dim)
#       mesh <- inla.mesh.1d(loc=t, degree=2, boundary="free") 
#     } else {
#       stop("For 2D, mesh must be supplied as argument xt$mesh in s(...,xt = )")
#     }
#   } else {
#     if (class(object$xt$mesh) != "inla.mesh") stop("xt must be NULL or an inla.mesh object")
#     mesh <- object$xt$mesh 
#   }
#   # model matrix: projects parameters to observation locations on mesh 
#   object$X <- as.matrix(inla.spde.make.A(mesh, x))
#   # compute finite element matrices used as smoothing penalty matrices 
#   inlamats <- inla.mesh.fem(mesh)
#   object$S <- list()
#   object$S[[1]] <- as.matrix(inlamats$c1)
#   object$S[[2]] <- 2 * as.matrix(inlamats$g1)
#   object$S[[3]] <- as.matrix(inlamats$g2)
#   # L is a matrix with a column for each smoothing parameter (tau, kappa) 
#   # and a row for each smoothing matrix (c1, g1, g2). 
#   # The (i,j)^th entry of L contains the power that smoothing parameter i 
#   # is computed to before being multiplied by smoothing matrix j. 
#   # E.g. If (1, 2) has value 4, then smoothing parameter 2 (kappa) is taken
#   # to the power 4 before being multiplied by smoothing matrix 1 (c1): i.e. kappa^4*c1
#   # All of these computations for each element of L are then summed to create a single
#   # smoothing matrix. 
#   object$L <- matrix(c(2,2,2,4,2,0), ncol = 2)
#   # Rank is the basis dimension, it is repeated for each smoothing matrix 
#   object$rank <- rep(object$bs.dim,3)
#   # As kappa > 0, the null space of the Matern SPDE is empty 
#   object$null.space.dim <- 0 
#   # Save the mesh
#   object$mesh <- mesh
#   object$df <- ncol(object$X)     # maximum DoF (if unconstrained)
#   # Give object a class
#   class(object) <- "spde.smooth" 
#   return(object)
# }
# 
# # Prediction function for the `spde' smooth class
# # See ?smooth.construct in mgcv for details on input and output 
# Predict.matrix.spde.smooth <- function(object, data){
#   dim <- length(object$term) 
#   if (dim > 2 | dim < 1) stop("SPDE Matern can only be fit in 1D or 2D.")
#   if (dim == 1) {
#     x <- data[[object$term]]
#   } else {
#     x <- matrix(0, nr = length(data[[1]]), nc = 2) 
#     x[,1] <- data[[object$term[1]]]
#     x[,2] <- data[[object$term[2]]]  
#   }
#   Xp <- inla.spde.make.A(object$mesh, x)
#   return(as.matrix(Xp))
# }
