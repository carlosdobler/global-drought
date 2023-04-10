
library(tidyverse)
library(furrr)
options(future.fork.enable = T)
plan(multicore, workers = 8)

vars_long <- c("precipitation", "wind_speed", "maximumum_temperature", "minimum_temperature", "relative_humidity", "surf_solar_radiation_down")
vars <- c("pr", "sfcWind", "tasmax", "tasmin", "hurs", "rsds")


model <- "MRI-ESM2-0"
g <- "gn"
r <- "r1i1p1f1"

walk2(vars_long, vars, function(var_long, var){
  
  # var_long = vars_long[6]
  # var = vars[6]
  
  print(" ")
  print(var_long)
  
  dd <- str_glue("/mnt/bucket_mine/cmip6/nex/daily/{model}")
  if(!file.exists(dd)){
    dir.create(dd)
  }
  
  
  d <- str_glue("/mnt/bucket_mine/cmip6/nex/daily/{model}/{var_long}")
  if(!file.exists(d)){
    dir.create(d)
  }
  
  
  
  future_walk(seq(1970,2014), function(yr){
    
    tictoc::tic(yr)
    
    root_sc <- 
      str_glue("https://ds.nccs.nasa.gov/thredds/fileServer/AMES/NEX/GDDP-CMIP6/{model}/historical/{r}/{var}")
    
    f <-
      str_glue("{var}_day_{model}_historical_{r}_{g}_{yr}.nc")
    
    download.file(str_glue("{root_sc}/{f}"),
                  destfile = str_glue("{d}/{f}"),
                  method = "wget",
                  quiet = T)
    
    tictoc::toc()
  })
  
  
  
  future_walk(seq(2015,2099), function(yr){
    # future_walk(seq(2029,2099), function(yr){
    
    tictoc::tic(yr)
    
    root_sc <- 
      str_glue("https://ds.nccs.nasa.gov/thredds/fileServer/AMES/NEX/GDDP-CMIP6/{model}/ssp585/{r}/{var}")
    
    f <-
      str_glue("{var}_day_{model}_ssp585_{r}_{g}_{yr}.nc")
    
    
    download.file(str_glue("{root_sc}/{f}"),
                  destfile = str_glue("{d}/{f}"),
                  method = "wget",
                  quiet = T)
    
    tictoc::toc()
  })
  
  
  
})

