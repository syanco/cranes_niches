#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script imputes missing niche data

'
Prepare and impute missng niche data
Usage:
prep_data.r <db> <table> 
clean_movement.r (-h | --help)

Control files:

Conda Environment: 
    niche_mix

Options:
-h --help     Show this screen.
-v --version  Show version.
-c --sp_cut   Minimum # of individuals to retain a species in analysis (defaults to 0)  
-d --dur_cut  Minimum indivdiual track duration in days (deafults to 365)
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  # library(here)
  
  .wd <- getwd()
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .table <- "event_clean_new"
  # .out <- file.path(.wd, "out/wi_ind_ax_yrs.csv")
  
} else {
  library(docopt)
  # library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  .table <- ag$table
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
    library(MVNH)
    library(zoo)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

# # declare vars of interest
# vars <- c("dist2water", "ghm", "prop_crop", "evi", "altitude")

#---- Initialize database ----#
message("Initializing database...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF)
invisible(assert_that(length(dbListTables(db))>0))


#---- Prep and Impute Data ----#

#-- Create Event Table
message("Extracting events from database...")
evt0 <- tbl(db, .table) %>% 
  collect()


#-- Rename vars
message("renaming variables...")

evt1 <- evt0 %>% 
  rename(
    #dist2water = dist2water,
    # ghm = value_modification,
    evi = EVI,
    # lst = value_lst_day_1km
    #pfor_h = pfor_h
    # altitude = value_altitude
  ) 

#-- Calculate Proportion Crops
message("Calculating proportion crops...")

# extract habitat cols
habs <- evt1 %>% 
  select(tail(names(.), 31)) %>% # the last 31 cols should be the ESACCI data
  select(!`QA`) # drop the QA column

# calc total habitat
habtots <- rowSums(habs, na.rm = T) %>% 
  as.data.frame()  
colnames(habtots) <- "tot_hab"

# calc totla crops
crops <- evt1 %>% 
  select(`10`, `20`, `30`) 
croptots <- rowSums(crops, na.rm = T) %>% 
  as.data.frame()  
colnames(croptots) <- "crops"

# add crops to DF
evt2 <- evt1 %>% 
  bind_cols(habtots, croptots) %>%
  mutate(prop_crop = crops/tot_hab) %>% 
  filter(!is.na(prop_crop)) #remove missing hab datat (crops and total == 0...)


#-- Filter data insufficiencies
message("Filtering out data insufficiencies...")

evt3 <- evt2 %>%

  # filter out shoulder seasons
  # filter(!is.na(season)) %>%

  # check for group level sufficiency
  group_by(individual_id) %>%

  #remove groups with only a single obs per individual
  filter(n() > 1) %>%

  # remove groups with at least one variable entirely NA or NaN
  # filter(if_any(all_of(vars), ~all(!is.na(.)))) %>%

  ungroup()



# #-- Impute Missing Data
# message("Imputng missing data...")
# 
# # define maximum gap over which to impute
# maxgap <- 7 # don't interpolate over any gaps lareger than a week
# 
# evt4 <- evt3 %>%
#   group_by(individual_id) %>% # impute w/i individuals
#   arrange(timestamp) %>% # sort by time
#   # impute all niche vars
#   mutate(dist2water = na.approx(dist2water, maxgap = maxgap, rule = 2),
#          lst = na.approx(lst, maxgap = maxgap, rule = 2),
#          prop_crop = na.approx(prop_crop, maxgap = maxgap, rule = 2),
#          evi = na.approx(evi, maxgap = maxgap, rule = 2)) %>%
#   ungroup()


#-- Scale and center vars
message("Scalign and centering variables...")

evt5 <- evt3 %>%
  group_by(individual_id) %>% 
  arrange(timestamp) %>% # sort by time
  ungroup() %>%
  # group_by(species) %>% 
  mutate(prop_crop_scale = scale(prop_crop),
         dist2water_scale = scale(dist2water),
         evi_scale = scale(evi),
         lst_scale = scale(lst)
         # altitude_scale = scale(altitude),
         # pfor_h_scale = scale(pfor_h),
         # v_scale = scale(v),
         # rad = bearing*(pi/180),
         # rad_scale = sclae(rad)
  )


#---- Write back to DB ----#
message("writing output back to database...")

dbWriteTable(db, "niche_prep", evt5, append = F, overwrite = T)


#---- Finalize Script ----#
message("Closing database connection...")

dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
