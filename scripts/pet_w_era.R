
# CALCULATION OF PET (PENMAN MONTEITH) WITH ERA5 DATA



# SETUP ------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(stars)
library(furrr)
library(units)

options(future.fork.enable = T)
plan(multicore)


# era monthly data (1-day means)
dir_era <- "/mnt/bucket_mine/era/monthly"

# temporary directory
dir_tmp <- "/mnt/pers_disk/tmp"
dir.create(dir_tmp)

# results directory
dir_res <- "/mnt/bucket_mine/results/global_drought_ww/pet/era5"

# variables
vars <- c("windspeed", 
          "dewpoint-temp", 
          "tasmax", 
          "tasmin", 
          "solar-radiation", 
          "thermal-radiation", 
          "pressure")




# DOWNLOAD + PRE-PROCESS VARIABLES ---------------------------------------------

walk(vars, function(var){
  
  print(str_glue("Downloading + pre-processing {var}"))
  
  # Download files to temp directory
  
  list.dirs(dir_era) %>% 
    str_subset(var) %>% 
    list.files(full.names = T) %>% 
    str_subset(seq(1970,2020) %>% str_flatten("|")) %>%
    suppressWarnings() %>%
    
    future_walk(function(f){
      
      f <- f %>% str_replace("/mnt/bucket_mine/", "gs://clim_data_reg_useast1/")
      
      "gsutil cp {f} {dir_tmp}" %>% 
        str_glue() %>% 
        system(ignore.stdout = T, ignore.stderr = T)
      
    })
  
  
  # Concatenate years and shift central meridian
  
  ff <- 
    dir_tmp %>% 
    list.files(full.names = T) %>% 
    str_subset("cat", negate = T)
  
  "cdo -b F32 sellonlatbox,-180,180,-90.25,90 -cat {str_flatten(ff, ' ')} {dir_tmp}/cat_{var}.nc" %>% 
    str_glue() %>% 
    system(ignore.stdout = T, ignore.stderr = T)
  
  ff %>% 
    file.remove()
  
})



# PET CALCULATION --------------------------------------------------------------

# extract time dimension
t <- 
  dir_tmp %>% 
  list.files(full.names = T) %>% 
  str_subset("windspeed") %>%  
  read_ncdf(proxy = T) %>%
  suppressMessages() %>% 
  st_get_dimension_values("time") %>% 
  as_date()

# extract years
yrs <- 
  t %>% 
  year() %>% 
  unique()


# loop across years
walk(yrs, function(yr){
  
  print(yr)
  
  # index years
  i <- which(year(t) == yr)
  
  
  # LOAD VARS -----------------------------------------------------------------
  s_vars <- 
    map(vars, function(var){
      
      s <- 
        dir_tmp %>% 
        list.files(full.names = T) %>% 
        str_subset(var) %>% 
        
        # subset 1 year
        read_ncdf(ncsub = cbind(start = c(1, 1, first(i)),
                                count = c(NA,NA,length(i)))) %>% 
        
        suppressMessages() %>%
        setNames("v") %>% 
        st_set_dimensions("time", 
                          t[i], 
                          point = NA)
      
      
      # convert units
      if(var %in% c("dewpoint-temp", "tasmax", "tasmin")){
        
        s <- 
          s %>% 
          mutate(v = set_units(v, degC),
                 v = set_units(v, NULL))
        
      } else if(str_detect(var, "radiation")){
        
        s <- 
          s %>%  
          mutate(v = set_units(v, MJ/m^2),
                 v = set_units(v, NULL))
        
      } else if(var == "pressure"){
        
        s <- 
          s %>%  
          mutate(v = set_units(v, kPa),
                 v = set_units(v, NULL))
        
      } else {
        
        s <- 
          s %>%  
          mutate(v = set_units(v, NULL))
        
      }
      
      s %>% 
        setNames(var)
      
    }) %>%
    
    do.call(c, .)
  
  s_vars <- 
    s_vars %>% 
    setNames(vars %>% str_replace("-", "_"))
  
  
  # CALCULATE -----------------------------------------------------------------
  
  # prepare equation elements
  s_vars_f <- 
    s_vars %>% 
    
    mutate(
      windspeed = windspeed * 4.87 / log(67.8 * 10 - 5.42),
      
      tasmean = (tasmax + tasmin) / 2,
      
      e_tmax = 0.6108 * exp((17.27 * tasmax) / (tasmax + 237.3)),
      e_tmin = 0.6108 * exp((17.27 * tasmin) / (tasmin + 237.3)),
      es = (e_tmax + e_tmin) / 2,
      
      ea = 0.6108 * exp((17.27 * dewpoint_temp) / (dewpoint_temp + 237.3)),
      
      delta = 4098 * es / (tasmean + 237.3)^2,
      
      gamma = 1.013e-3 * pressure / (0.622 * 2.45),
      
      Rn = solar_radiation + thermal_radiation,
      
      G = 0.3 * Rn) %>% 
    
    select(delta, Rn, tasmean, gamma, windspeed, es, ea, G)
  
  
  # apply equation
  s_ET <- 
    s_vars_f %>% 
    mutate(
      
      numerator = 0.408 * delta * (Rn - G) + gamma * (900 / (tasmean + 273)) * windspeed * (es - ea),
      
      denominator = delta + gamma * (1 + 0.34 * windspeed),
      
      ET = numerator/denominator,
      
      ET = ifelse(ET < 0, 0, ET)
      
    ) %>% 
    select(ET)
  
  
  # SAVE ----------------------------------------------------------------------
  
  # define dimensions
  dim_lon <- ncdf4::ncdim_def(name = "longitude", 
                              units = "degrees_east", 
                              vals = s_ET %>% 
                                st_get_dimension_values(1))
  
  dim_lat <- ncdf4::ncdim_def(name = "latitude", 
                              units = "degrees_north", 
                              vals = s_ET %>% 
                                st_get_dimension_values(2))
  
  dim_time <- ncdf4::ncdim_def(name = "time", 
                               units = "days since 1970-01-01", 
                               vals = s_ET %>% 
                                 st_get_dimension_values(3) %>% 
                                 as_date() %>%
                                 as.integer())
  
  # define variable
  varis <- ncdf4::ncvar_def(name = "pet",
                            units = "mm",
                            dim = list(dim_lon, dim_lat, dim_time), 
                            missval = -99999)
  
  
  # create file
  f_name <- 
    str_glue("{dir_res}/era5_mon_mean-pet_{yr}.nc")
  
  if(file.exists(f_name)){
    file.remove(f_name)
    print("file replaced")
  }
  
  ncnew <- ncdf4::nc_create(filename = f_name, 
                            vars = varis,
                            force_v4 = TRUE)
  
  # write data
  ncdf4::ncvar_put(nc = ncnew, 
                   varid = varis, 
                   vals = s_ET %>% pull(1))
  
  ncdf4::nc_close(ncnew)
  
  
})

# delete temporary dir
unlink(dir_tmp, recursive = T)











