
# CALCULATION OF scPDSI WITH ERA5 DATA



# SETUP ------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(stars)
library(furrr)
library(units)

options(future.fork.enable = T)
plan(multicore)


# era monthly precip data (1-day means)
dir_pr <- "/mnt/bucket_mine/era/monthly/mean-daily-precip"

# era monthly pet data
dir_pet <- "/mnt/bucket_mine/results/global_drought_ww/pet/era5"

# temporary directory
dir_tmp <- "/mnt/pers_disk/tmp"
dir.create(dir_tmp)

# tiles directory
dir_tiles <- "/mnt/pers_disk/tiles"
dir.create(dir_tiles)

# results directory
dir_res <- "/mnt/bucket_mine/results/global_drought_ww/pdsi/era5"




# DOWNLOAD + PRE-PROCESS VARIABLES ---------------------------------------------


vars <- c("pr", "pet")

walk(vars, function(var){
  
  print(str_glue("Downloading + pre-processing {var}"))
  
  # Download files to temporary dir
  
  if(var == "pet"){
    dir_data <- dir_pet
    
  } else {
    dir_data <- dir_pr
    
  } 
  
  list.dirs(dir_data) %>% 
    list.files(full.names = T) %>%
    str_subset(seq(1970,2020) %>% str_flatten("|")) %>%
    
    future_walk(function(f){
      
      f <- f %>% str_replace("/mnt/bucket_mine/", "gs://clim_data_reg_useast1/")
      
      "gsutil cp {f} {dir_tmp}" %>% 
        str_glue() %>% 
        system(ignore.stdout = T, ignore.stderr = T)
      
    })
  
  
  # Concatenate years
  
  ff <- 
    dir_tmp %>% 
    list.files(full.names = T) %>% 
    str_subset("cat", negate = T)
  
  if(var == "pr"){
    
    # shift central meridian in case of precip
    "cdo sellonlatbox,-180,180,-90.25,90 -cat {str_flatten(ff, ' ')} {dir_tmp}/cat_{var}.nc" %>% 
      str_glue() %>% 
      system(ignore.stdout = T, ignore.stderr = T)
  
  } else {
    
    "cdo cat {str_flatten(ff, ' ')} {dir_tmp}/cat_{var}.nc" %>% 
      str_glue() %>% 
      system(ignore.stdout = T, ignore.stderr = T)
    
  }
  
  ff %>% file.remove()
  
})




# TILE ------------------------------------------------------------------------

s_proxy <- 
  dir_tmp %>% 
  list.files(full.names = T) %>% 
  first() %>% 
  read_ncdf(ncsub = cbind(start = c(1, 1, 1),
                          count = c(NA,NA,1))) %>% 
  suppressMessages() %>% 
  adrop()


# Create land layer

rast_reference_0.05 <- 
  s_proxy %>%
  st_bbox() %>% 
  st_as_stars(dx = 0.05, dy = 0.05, values = -9999) 

land <- 
  "/mnt/bucket_mine/misc_data/ne_50m_land/ne_50m_land.shp" %>%
  st_read(quiet = T) %>%
  mutate(a = 1) %>%
  select(a) %>%
  st_rasterize(rast_reference_0.05)

land <- 
  land %>%
  st_warp(s_proxy, use_gdal = T, method = "max") %>%
  suppressWarnings() %>% 
  setNames("a") %>%
  mutate(a = ifelse(a == -9999, NA, a))


sz <- 25 # chunk size

# split into chunks
lims <- 
  map(c(1,2), function(dim_id){
    
    d <- 
      dim(s_proxy)[dim_id] %>%
      seq_len()
    
    n <- 
      round(dim(s_proxy)[dim_id]/sz)
    
    l <- split(d, 
               ceiling(d/(length(d)/n))) %>% 
      map(~c(first(.x), last(.x)))
    
    return(l)
    
  })

# identify tiles with land
tiles <-
  
  future_imap(lims[[1]], function(lon_ch, lon_i){
    imap(lims[[2]], function(lat_ch, lat_i){
      
      s_proxy_tile <- 
        s_proxy %>%
        slice(longitude, lon_ch[1]:lon_ch[2]) %>%
        slice(latitude, lat_ch[1]:lat_ch[2])
      
      land_tile <- 
        st_warp(land,
                s_proxy_tile)
      
      pol_tile <- 
        s_proxy_tile %>%
        mutate(a = 1) %>% 
        select(a) %>% 
        st_as_sf(as_points = F, merge = T) %>%
        mutate(lon_ch = lon_i,
               lat_ch = lat_i)
      
      pol_tile <- 
        pol_tile %>% 
        mutate(cover = ifelse(all(is.na(pull(land_tile))) | all(is.na(pull(s_proxy_tile))), F, T)) %>% 
        select(-a)
      
      return(pol_tile)
      
    }) %>% 
      bind_rows()
    
  }) %>% 
  bind_rows()



# PDSI CALCULATION ------------------------------------------------------------

# whole time vector
d <- 
  dir_tmp %>% 
  list.files(full.names = T) %>% 
  first() %>% 
  read_ncdf(proxy = T) %>% 
  suppressMessages() %>% 
  st_get_dimension_values(3) %>% 
  as_date()

# index final year of calibration period
cal_f <- 
  d %>% 
  year %>% 
  unique() %>% 
  {. == 2000} %>% 
  which()



# Process tiles with land

future_pwalk(st_drop_geometry(tiles) %>% filter(cover == T), function(lon_ch, lat_ch, ...){
  
  
  # LOAD VARS -----------------------------------------------------------------
  s_vars <- 
    map(vars, function(var){
      
      # subset matrix
      nc <- 
        cbind(start = c(lims[[1]][[lon_ch]][1],
                        lims[[2]][[lat_ch]][1],
                        1),
              count = c(lims[[1]][[lon_ch]][2] - lims[[1]][[lon_ch]][1] + 1,
                        lims[[2]][[lat_ch]][2] - lims[[2]][[lat_ch]][1] + 1,
                        NA))
      
      
      s <- 
        dir_tmp %>% 
        list.files(full.names = T) %>% 
        str_subset(var) %>% 
        
        read_ncdf(ncsub = nc) %>% 
        
        suppressMessages() %>% 
        setNames("v") %>% 
        st_set_dimensions(3, 
                          st_get_dimension_values(., 3) %>% as_date(), 
                          point = NA)
      
      
      # convert units
      if(var == "pr"){
        
        s <- 
          s %>% 
          mutate(v = set_units(v, mm),
                 v = set_units(v, NULL))
        
      } else {
        
        s <- 
          s %>%  
          mutate(v = set_units(v, NULL))
        
        
      }
      
      s %>% 
        setNames(var)
      
    }) %>%
    
    {do.call(c, c(., along = "v"))} %>% 
    setNames("vv")
  
  
  # CALCULATE ------------------------------------------------------------------
  
  s_pdsi <-     
    s_vars %>% 
    st_apply(c(1,2), function(x){
      
      scPDSI::pdsi(P = x[,1],
                   PE = x[,2],
                   sc = T,
                   cal_start = 1,
                   cal_end = cal_f)$X %>% 
        as.vector()
      
    },
    .fname = "time") %>% 
    st_set_dimensions("time", values = d) %>% 
    aperm(c(2,3,1)) %>% 
    setNames("pdsi")
  
  
  lon_ch_ <- str_pad(lon_ch, 2, side = "left", pad = "0")
  lat_ch_ <- str_pad(lat_ch, 2, side = "left", pad = "0")
  
  saveRDS(s_pdsi, str_glue("{dir_tiles}/s_pdsi_{lon_ch_}_{lat_ch_}.rds"))
  
})


# Process tiles with no land (empty)

future_pwalk(st_drop_geometry(tiles) %>% filter(cover == F), function(lon_ch, lat_ch, ...){
  
  nc <- 
    cbind(start = c(lims[[1]][[lon_ch]][1],
                    lims[[2]][[lat_ch]][1],
                    1),
          count = c(lims[[1]][[lon_ch]][2] - lims[[1]][[lon_ch]][1] + 1,
                    lims[[2]][[lat_ch]][2] - lims[[2]][[lat_ch]][1] + 1,
                    NA))
  
  s_pdsi <- 
    
    dir_tmp %>% 
    list.files(full.names = T) %>% 
    first() %>% 
    
    read_ncdf(ncsub = nc) %>% 
    
    suppressMessages() %>% 
    st_set_dimensions(3, 
                      st_get_dimension_values(., 3) %>% as_date(), 
                      point = NA) %>% 
    
    setNames("pdsi") %>% 
    drop_units() %>% 
    mutate(pdsi = NA_real_)
  
  lon_ch_ <- str_pad(lon_ch, 2, side = "left", pad = "0")
  lat_ch_ <- str_pad(lat_ch, 2, side = "left", pad = "0")
  
  saveRDS(s_pdsi, str_glue("{dir_tiles}/s_pdsi_{lon_ch_}_{lat_ch_}.rds"))
  
})



# MOSAICK ----------------------------------------------------------------------

mos <- 
  
  future_map(seq_along(lims[[1]]) %>% str_pad(2, "left", "0"), function(col_){
    
    dir_tiles %>% 
      list.files(full.names = T) %>% 
      str_subset(str_glue("_{col_}_")) %>%
      map(readRDS) %>% 
      {do.call(c, c(., along = 2))}
    
  }) %>%
  
  {do.call(c, c(., along = 1))}


mos <- 
  mos %>% 
  st_set_dimensions(1, 
                    values = st_get_dimension_values(s_proxy, 1, center = F)) %>% 
  st_set_dimensions(2, 
                    values = st_get_dimension_values(s_proxy, 2, center = F)) %>% 
  st_set_crs(4326)



# SAVE ----------------------------------------------------------------------

walk(d %>% year() %>% unique(), function(yr){
  
  print(yr)
  
  s <- 
    mos %>% 
    filter(year(time) == yr)
  
  # define dimensions
  dim_lon <- ncdf4::ncdim_def(name = "longitude", 
                              units = "degrees_east", 
                              vals = s %>% 
                                st_get_dimension_values(1))
  
  dim_lat <- ncdf4::ncdim_def(name = "latitude", 
                              units = "degrees_north", 
                              vals = s %>% 
                                st_get_dimension_values(2))
  
  dim_time <- ncdf4::ncdim_def(name = "time", 
                               units = "days since 1970-01-01", 
                               vals = s %>% 
                                 st_get_dimension_values(3) %>% 
                                 as_date() %>%
                                 as.integer())
  
  # define variables
  varis <- ncdf4::ncvar_def(name = "pdsi",
                            units = "",
                            dim = list(dim_lon, dim_lat, dim_time), 
                            missval = -99999)
  
  
  # create file
  f_name <- 
    str_glue("{dir_res}/era5_mon_pdsi_{yr}.nc")
  
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
                   vals = s %>% pull(1))
  
  ncdf4::nc_close(ncnew)
  
  
})

unlink(dir_tiles, recursive = T)
unlink(dir_tmp, recursive = T)





