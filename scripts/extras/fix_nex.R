
library(tidyverse)
library(furrr)
library(stars)
options(future.fork.enable = T)
plan(multicore, workers = 8)

vars_long <- c("precipitation", "wind_speed", "maximum_temperature", "minimum_temperature", "relative_humidity", "surf_solar_radiation_down")
vars <- c("pr", "sfcWind", "tasmax", "tasmin", "hurs", "rsds")


model <- "GFDL-ESM4"
g <- "gr1"
r <- "r1i1p1f1"

# nex monthly data (1-day means)
dir_nex_m <- str_glue("/mnt/bucket_mine/cmip6/nex/monthly/{model}")
dir.create(dir_nex_m)

dir_nex_d <- str_glue("/mnt/bucket_mine/cmip6/nex/daily/{model}")
dir.create(dir_nex_d)



for(i in 1:6){
  
  var_long = vars_long[i]
  var = vars[i]
  
  print(str_glue(" "))
  print(str_glue("{var_long}"))
  
  
  # ************
  
  dir_nex_m_v <- str_glue("{dir_nex_m}/{var_long}")
  
  if(!dir.exists(dir_nex_m_v)){
    
    print(str_glue("Aggregating"))
    
    dir.create(dir_nex_m_v)
    
    
    ff <- character()
    i <- 0  
    
    while(length(ff) < 1){
      
      ff <- 
        "{dir_nex_d}/{var_long}" %>% 
        str_glue() %>% 
        list.files()
      
      if(i > 0){
        print(str_glue("...retrying reading files..."))
        Sys.sleep(5)
      }
      
      i <- i+1
      
    }
    
    future_walk(ff, function(f){
      
      final_name <- str_glue("{dir_nex_m}/{var_long}/{str_replace(f, '_day_', '_mon_')}")
      
      "cdo -setcalendar,standard -settunits,days -setday,1 -monmean {dir_nex_d}/{var_long}/{f} {final_name}" %>% 
        str_glue() %>% 
        system(ignore.stdout = T, ignore.stderr = T)
      
    })
    
    Sys.sleep(8)
    print(str_glue("...done"))
    
  }
  
  
  # ************
  
  
  walk(seq(1970,2014), function(yr){
    
    f <- 
      str_glue("{dir_nex_m}/{var_long}/{var}_mon_{model}_historical_{r}_{g}_{yr}.nc")
    
    if(!file.exists(f)){
      print(yr)
    }
    
    len_time <- 
      f %>% 
      read_ncdf(proxy = T) %>% 
      st_get_dimension_values(3) %>%
      suppressMessages() %>% 
      length()
    
    while(len_time < 12){
      
      print(str_glue("downloading {yr} daily data"))
      
      root_sc <- 
        str_glue("https://ds.nccs.nasa.gov/thredds/fileServer/AMES/NEX/GDDP-CMIP6/{model}/historical/{r}/{var}")
      
      f_day <-
        str_glue("{var}_day_{model}_historical_{r}_{g}_{yr}.nc")
      
      download.file(str_glue("{root_sc}/{f_day}"),
                    destfile = str_glue("{dir_nex_d}/{var_long}/{f_day}"),
                    method = "wget",
                    quiet = T)
      
      
      "cdo -setcalendar,standard -settunits,days -setday,1 -monmean {dir_nex_d}/{var_long}/{f_day} {f}" %>% 
        str_glue() %>% 
        system(ignore.stdout = T, ignore.stderr = T)
      
      len_time <- 
        f %>% 
        read_ncdf(proxy = T) %>% 
        st_get_dimension_values(3) %>%
        suppressMessages() %>% 
        length()
      
      
    }
    
  })
  
  
  # *********
  
  walk(seq(2015,2099), function(yr){
    
    # print(yr)
    
    f <- 
      str_glue("{dir_nex_m}/{var_long}/{var}_mon_{model}_ssp585_{r}_{g}_{yr}.nc")
    
    if(!file.exists(f)){
      print(yr)
    }
    
    # if(!file.exists(f)){
    #   
    #   print(str_glue("aggregating {yr}"))
    # 
    #   d_name <- f %>% str_replace("_mon_", "_day_") %>% str_replace("/monthly/", "/daily/")
    # 
    #   "cdo -setcalendar,standard -settunits,days -setday,1 -monmean {d_name} {f}" %>%
    #     str_glue() %>%
    #     system(ignore.stdout = T, ignore.stderr = T)
    # 
    # }
    
    len_time <- 
      f %>% 
      read_ncdf(proxy = T) %>% 
      st_get_dimension_values(3) %>%
      suppressMessages() %>% 
      length()
    
    
    while(len_time < 12){
      
      print(str_glue("downloading {yr} daily data"))
      
      root_sc <- 
        str_glue("https://ds.nccs.nasa.gov/thredds/fileServer/AMES/NEX/GDDP-CMIP6/{model}/ssp585/{r}/{var}")
      
      f_day <-
        str_glue("{var}_day_{model}_ssp585_{r}_{g}_{yr}.nc")
      
      download.file(str_glue("{root_sc}/{f_day}"),
                    destfile = str_glue("{dir_nex_d}/{var_long}/{f_day}"),
                    method = "wget",
                    quiet = T)
      
      
      "cdo -setcalendar,standard -settunits,days -setday,1 -monmean {dir_nex_d}/{var_long}/{f_day} {f}" %>% 
        str_glue() %>% 
        system(ignore.stdout = T, ignore.stderr = T)
      
      
      len_time <- 
        f %>% 
        read_ncdf(proxy = T) %>% 
        st_get_dimension_values(3) %>%
        suppressMessages() %>% 
        length()
      
      
    }
    
  })
  
}

