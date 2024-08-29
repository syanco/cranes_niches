#---- Input Parameters ----#
if(interactive()) {
  # library(here)
  
  .wd <- '~/projects/dynamic_niches'
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .table <- "niche_prep"
  .out <- file.path(.wd, "out/wi_ind_ax_yrs.csv")
  
} else {
  library(docopt)
  # library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  .table <- ag$table
  .out <- makePath(ag$out)
}

#---- Initialize Environment ----#

# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'src/init/startup.r'))

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    # library(MVNH)
    # library(INLA)
    library(ggfortify)
    library(zoo)
    library(patchwork)
    library(colorspace)
    library(tidyverse)
    library(sf)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)


#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

# create event table
evt0 <- tbl(db, .table) %>% 
  collect()

phen_summ <- read_csv("out/phenology_summary.csv")

# # extract phenology info from segmentaitons
# species_list <- c("Anthropoides_virgo", "Grus_grus", "Grus_nigricollis", "Grus_vipio", "Balearica_pavonina")

segP <- "~/projects/dynamic_niches/data/shiny_track_segmentation/"

segFs <- list.files(segP, full.names = T)

#---- Get Seasonal Centroids ----#

# make vector of inds over which to loop
inds <- evt0 %>% 
  pull(individual_id) %>% 
  unique()

cent_out <- list()
for(i in 1:length(inds)){
  
  ind <- inds[i]
  
  # Extract individual events
  evt_foc <- evt0 %>% 
    filter(individual_id == ind)
  
  # extract species
  sp <- evt_foc %>% 
    summarize(sp = species[1]) %>% 
    pull(sp)
  
  #-- Get Individual Segmentation --#
  
  # get seg fike ind
  file_idx <- grep(ind, segFs)
  
  # make sure segmentaito file exists
  if(identical(file_idx, integer(0))){
    message(glue("No segmentation matching individual {ind}. Better luck next time brah..."))
    next
  } #if
  
  # load in segmentaion
  segs <- read_csv(segFs[file_idx]) %>% 
    mutate(yr = year(Date))
  
  # Loop over years in the data
  yrs <- evt_foc %>% pull(yr) %>% unique()
  
  cents_tmp <- list()
  for(j in 1:length(yrs)){
    
    #-- Summer centroid --#
    
    s_start <- segs %>% 
      filter(yr == yrs[j],
             Status == "End Spring")
    s_stop <- segs %>% 
      filter(yr == yrs[j],
             Status == "Start Fall")
    
    if(nrow(s_start) > 0 & nrow(s_stop) > 0){
      
      # extract summer events
      evt_summ <- evt_foc %>% 
        filter(timestamp > s_start$Date & timestamp < s_stop$Date)
      
      
    if(nrow(evt_summ) > 0){
      # get centroid
      summ_cent <- evt_summ %>% 
        st_as_sf(coords = c("lat", "lon"), crs = 4326) %>% 
        summarize(geometry = st_union(geometry)) %>% 
        st_centroid() %>% 
        mutate(lon = unlist(map(.$geometry,1)),
               lat = unlist(map(.$geometry,2))) %>% 
        as.data.frame() %>% 
        select(-geometry)  
      
      
      summ_df <- data.frame(
        species = sp,
        individual_id = ind,
        season = "summer",
        yr = yrs[j]
      ) %>% 
        cbind(summ_cent)
     } else{# if nrow(evt_summ)
        message("No summer data for individual {ind} in {yrs[j]}...")
      }
    } else{ # if segs exist
      message(glue("Unable to segment summer period for individual {ind}..."))
    } #else

  
  #-- Winter Centroid --#
  
  w_start <- segs %>% 
    filter(yr == yrs[j],
           Status == "End Fall")
  w_stop <- segs %>% 
    filter(yr == yrs[j],
           Status == "Start Spring")
  
  if(nrow(w_start) > 0 & nrow(w_stop) > 0){
    
    # extract summer events
    evt_wint <- evt_foc %>% 
      filter(timestamp > w_start$Date & timestamp < w_stop$Date)
    
    if(nrow(evt_wint) > 0){
    # get centroid
    wint_cent <- evt_wint %>% 
      st_as_sf(coords = c("lat", "lon"), crs = 4326) %>% 
      summarize(geometry = st_union(geometry)) %>% 
      st_centroid() %>% 
      mutate(lon = unlist(map(.$geometry,1)),
             lat = unlist(map(.$geometry,2))) %>% 
      as.data.frame() %>% 
      select(-geometry)  
    
    wint_df <- data.frame(
      species = sp,
      individual_id = ind,
      season = "winter",
      yr = yrs[j]
    ) %>% 
      cbind(wint_cent)
    
    } else{ # if nrow(evt_wint)
      message("No winter data for individual {ind} in {yrs[j]}...")
    } # if 
    
  }else{ # if segs exist
    message(glue("Unable to segment winter period for individual {ind}..."))
  } #else
  
  if(exists("summ_df") & exists("wint_df")){
  cents_tmp[[j]] <- rbind(summ_df, wint_df)
  }else{
    if(exists("summ_df")){
      cents_tmp[[j]] <- summ_df
    }else{
      if(exists("wint_df")){
        cents_tmp[[j]] <- wint_df
      }
    }
  }
  } # j
  
  cent_out[[i]] <- do.call("rbind", cents_tmp)
} #i

out_df <- do.call("rbind", cent_out)

write_csv(out_df, file = "out/individual_seasonal_centroids.csv")
