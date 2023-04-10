
# CALCULATION OF PET (PENMAN MONTEITH) WITH NEX CMIP6 DATA

model <- c("GFDL-ESM4",
           "MPI-ESM1-2-HR",
           "MRI-ESM2-0",
           "UKESM1-0-LL",
           "IPSL-CM6A-LR")[5]


# SETUP ------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(stars)
library(furrr)
library(units)

options(future.fork.enable = T)
plan(multicore)


# nex daily data
dir_nex_d <- str_glue("/mnt/bucket_mine/cmip6/nex/daily/{model}")

# nex monthly data (1-day means)
dir_nex_m <- str_glue("/mnt/bucket_mine/cmip6/nex/monthly/{model}")

# temporary directory
dir_tmp <- "/mnt/pers_disk/tmp"
dir.create(dir_tmp)

# results directory
dir_res <- str_glue("/mnt/bucket_mine/results/global_drought_ww/pet/nex/{model}")

if(!dir.exists(dir_res)){
  dir.create(dir_res)
}

#variables
vars <- c("maximum_temperature", 
          "minimum_temperature", 
          "relative_humidity", 
          "wind_speed", 
          "surf_solar_radiation_down")




{

# DOWNLOAD + PRE-PROCESS VARIABLES ---------------------------------------------

walk(vars, function(v){
  
  print(str_glue("Downloading + pre-processing {v}"))
  
  # Download files to temporary directory
  
  list.dirs(dir_nex_m) %>% 
    str_subset(v) %>% 
    list.files(full.names = T) %>% 
    # str_subset(seq(1970,2020) %>% str_flatten("|")) %>%
    
    future_walk(function(f){
      
      f <- 
        f %>% 
        str_replace("/mnt/bucket_mine/", "gs://clim_data_reg_useast1/")
      
      "gsutil cp {f} {dir_tmp}" %>% 
        str_glue() %>% 
        system(ignore.stdout = T, ignore.stderr = T)
      
    })
  
  
  # Concatenate years and shift central meridian
  
  ff <- 
    dir_tmp %>% 
    list.files(full.names = T) %>% 
    str_subset("cat", negate = T)
  
  "cdo sellonlatbox,-180,180,-60,90 -cat {str_flatten(ff, ' ')} {dir_tmp}/cat_{v}.nc" %>% 
    str_glue() %>% 
    system(ignore.stdout = T, ignore.stderr = T)
    
  ff %>% file.remove()
  
})


# verify time dimension
tb <- 
  imap_dfr(vars %>% set_names, function(v, i){
    
    t <- 
      dir_tmp %>% 
      list.files(full.names = T) %>% 
      str_subset(v) %>% 
      
      read_ncdf(proxy = T) %>% 
      suppressMessages() %>% 
      st_get_dimension_values(3) %>% 
      as_date()
    
    tibble(var = i,
           timesteps = length(t),
           t_i = first(t),
           t_f = last(t))
  })

print(tb)



# PREPARE CONSTANTS -----------------------------------------------------------

s_proxy <- 
  dir_tmp %>% 
  list.files(full.names = T) %>% 
  .[1] %>% 
  read_ncdf(ncsub = cbind(start = c(1, 1, 1),
                          count = c(NA,NA,1))) %>% 
  suppressMessages() %>%
  setNames("v") %>% 
  adrop() %>% 
  drop_units()


dir_constants <- "/mnt/pers_disk/constants"

# dir.create(dir_constants)


# Elevation 

# z_file <- "/mnt/bucket_mine/era/era5_geopotential.nc" 
#   
# "cdo sellonlatbox,-180,180,-90.25,90 {z_file} {dir_constants}/z.nc" %>% 
#   str_glue() %>% 
#   system(ignore.stdout = T, ignore.stderr = T)

z <- 
  str_glue("{dir_constants}/z.nc") %>% 
  read_ncdf() %>% 
  suppressMessages() %>% 
  adrop() %>% 
  drop_units() %>% 
  st_warp(s_proxy)

gamma <- 
  z %>% 
  mutate(z = ifelse(z > 10000, 10000, z),
         P = 101.3 * ((293 - 0.0065 * z) / 293)^5.26,
         gamma = 1.013e-3 * P / (0.622 * 2.45)) %>% 
  select(z,
         gamma)

gamma <- 
  rep(list(gamma), 12) %>% 
  {do.call(c, c(., along = "time"))}
  


# Extraterrestrial radiation

# # Download
# 
# "/mnt/bucket_mine/era/monthly/mean-toa-incident-solar-radiation" %>%
#   list.files(full.names = T) %>%
#   str_subset(seq(1970,2000) %>% str_flatten("|")) %>%
# 
#   future_walk(function(f){
# 
#     f <- f %>% str_replace("/mnt/bucket_mine/", "gs://clim_data_reg_useast1/")
# 
#     "gsutil cp {f} {dir_constants}" %>%
#       str_glue() %>%
#       system(ignore.stdout = T, ignore.stderr = T)
# 
#   })
# 
# # Concatenate years, calculate ymon means, and shift central meridian
# 
# ff <-
#   dir_constants %>%
#   list.files(full.names = T) %>%
#   str_subset("radiation")

# "cdo sellonlatbox,-180,180,-90.25,90 -ymonmean -cat {str_flatten(ff, ' ')} {dir_constants}/tisr_ymonmean.nc" %>%
#   str_glue() %>%
#   system(ignore.stdout = T, ignore.stderr = T)
# 
# ff %>% file.remove()

tisr <- 
  str_glue("{dir_constants}/tisr_ymonmean.nc") %>% 
  read_ncdf() %>%
  suppressMessages() %>% 
  mutate(tisr = set_units(tisr, MJ/m^2),
         tisr = set_units(tisr, NULL)) %>% 
  st_warp(s_proxy)




# CALCULATE PET ---------------------------------------------------------------

# time vector
t <- 
  dir_tmp %>% 
  list.files(full.names = T) %>% 
  .[1] %>% 
  read_ncdf(proxy = T) %>%
  suppressMessages() %>% 
  st_get_dimension_values("time") %>% 
  as_date()

# extract years
yrs <- 
  t %>% 
  year() %>% 
  unique()


# Loop through years

walk(yrs, function(yr){
  
  print(str_glue("Processing {yr}"))
  
  # index years
  i <- which(year(t) == yr)
  
  
  # LOAD VARS -----------------------------------------------------------------
  
  s_vars <- 
    map(vars, function(v){
      
      s <- 
        dir_tmp %>% 
        list.files(full.names = T) %>% 
        str_subset(v) %>% 
        read_ncdf(ncsub = cbind(start = c(1, 1, first(i)),
                                count = c(NA,NA,length(i)))) %>% 
        suppressMessages() %>%
        setNames("v") 
      
      # corroborate
      yr_ref <- st_get_dimension_values(s, 3) %>% year() %>% unique()
      if(!identical(yr, yr_ref)){
        print(str_glue("   wrong years!"))
      }
      
      s <- 
        s %>% 
        st_set_dimensions("time", 
                          t[i], 
                          point = NA)
      
      # conversions
      if(str_detect(v, "temperature")){
        
        s <- 
          s %>% 
          mutate(v = set_units(v, degC),
                 v = set_units(v, NULL))
        
      } else if(str_detect(v, "radiation")){
        
        s <- 
          s %>%  
          mutate(v = set_units(v, MJ/d/m^2),
                 v = set_units(v, NULL))
        
      } else {
        
        s <- 
          s %>%  
          mutate(v = set_units(v, NULL))
        
      }
      
      s %>% 
        setNames(v)
      
    }) %>%
    
    do.call(c, .)
  
  
  s_gamma <- 
    gamma %>% 
    st_set_dimensions(3, values = st_get_dimension_values(s_vars, 3))
  
  s_tisr <- 
    tisr %>% 
    st_set_dimensions(3, values = st_get_dimension_values(s_vars, 3))
  
  
  s_vars_f <- 
    
    s_vars %>%
    
    c(s_gamma) %>% 
    c(s_tisr) %>% 
    
    mutate(
      tasmean = (maximum_temperature + minimum_temperature) / 2,
      
      e_tmax = 0.6108 * exp((17.27 * maximum_temperature) / (maximum_temperature + 237.3)),
      e_tmin = 0.6108 * exp((17.27 * minimum_temperature) / (minimum_temperature + 237.3)),
      es = (e_tmax + e_tmin) / 2,
      
      relative_humidity = ifelse(relative_humidity < 0, 0, relative_humidity),
      ea = es * relative_humidity / 100,
      
      delta = 4098 * es / (tasmean + 237.3)^2,
      
      Rns = (1 - 0.23) * surf_solar_radiation_down,
      
      Rso = (0.75 + 2e-5 * z) * tisr,
      Rso = ifelse(Rso == 0, 1e-10, Rso),
      Rnl = 4.903e-9 * ((maximum_temperature + 273.16)^4 + (minimum_temperature + 273.16)^4) / 2 *
        (0.34 - 0.14 * sqrt(ea)) *
        (1.35 * surf_solar_radiation_down / Rso - 0.35),
      
      Rn = Rns - Rnl,
      
      G = 0.3 * Rn) %>% 
    
    select(delta, Rn, tasmean, gamma, wind_speed, es, ea, G)
  
  
  s_ET <- 
    s_vars_f %>% 
    mutate(
      
      numerator = 0.408 * delta * (Rn - G) + gamma * (900 / (tasmean + 273)) * wind_speed * (es - ea),
      
      denominator = delta + gamma * (1 + 0.34 * wind_speed),
      
      ET = numerator/denominator,
      
      ET = ifelse(ET < 0, 0, ET)
      
    ) %>% 
    select(ET)
  
  
  # SAVE ----------------------------------------------------------------------
  
  # define dimensions
  dim_lon <- ncdf4::ncdim_def(name = "lon", 
                              units = "degrees_east", 
                              vals = s_ET %>% 
                                st_get_dimension_values(1))
  
  dim_lat <- ncdf4::ncdim_def(name = "lat", 
                              units = "degrees_north", 
                              vals = s_ET %>% 
                                st_get_dimension_values(2))
  
  dim_time <- ncdf4::ncdim_def(name = "time", 
                               units = "days since 1970-01-01", 
                               vals = s_ET %>% 
                                 st_get_dimension_values(3) %>% 
                                 as_date() %>%
                                 as.integer())
  
  # define variables
  varis <- ncdf4::ncvar_def(name = "pet",
                            units = "mm",
                            dim = list(dim_lon, dim_lat, dim_time), 
                            missval = -99999)
  
  
  # create file
  f_name <- 
    str_glue("{dir_res}/nex-{model}_mon_mean-pet_{yr}.nc")
  
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


unlink(dir_tmp, recursive = T)

}




