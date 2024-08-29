#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script preps and cleans data for the Dynamic Niches Project
# See project documentation for details about anticipated directory structure.

'
Prep and clean data.
Expects db to be of format mosey_db, with a table named "anno_join_2022-11-16" 
which contains the already-annotated events.

Usage:
clean_movement.r <db> <table> <out_table> [<sp_cut> <dur_cut>]
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
  # .table <- "anno_join_2024-07-01"
  .table <- "niche_prep"
  .out_table <- "event_clean"
  .sp_cut <- 5
  .min_dur <- 50

} else {
  library(docopt)
  # library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()

  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  .table <- ag$table
  .out_table <- ag$out_table
  .sp_cut <- ag$sp_cut
  .min_dur <- ag$dur_cut
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
    library(sf)
    library(geosphere)
    library(lubridate)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

# Defualts for optional agrd
message("Minimum no. of individuals per species:")
(.sp_cut <- ifelse(is.null(.sp_cut), 5, .sp_cut))

message("Minimum track duration:")
(.min_dur <- ifelse(is.null(.min_dur), 50, .min_dur))

#---- Load control files ----#



#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

# create event table
evt0 <- tbl(db, .table)

# get list of inidviduals to remove from study base don bad coords
indtb <- tbl(db, "individual") %>%  collect()

# get key to link inds to spp and studies
spp_key <- indtb %>% 
  select(taxon_canonical_name, individual_id, study_id)

# dbBegin(db)

#---- Perform analysis ----#


#-- Clean Outliers

message("Beginning outlier removal...")

cnt <- evt0 %>% 
  group_by(individual_id) %>% 
  summarise(n = n())

# vec of individuals with bad spp info
rm_inds <- indtb %>% 
  filter(is.na(taxon_canonical_name)) %>% 
  pull(individual_id)

# get step lengths and turn angles across dataset
evt1 <- evt0 %>% 
  collect() %>% 
  filter(individual_id %notin% rm_inds) %>%                                 # rm inds w/ bad species info
  left_join(spp_key) %>%                                                    # add species names to event table
  group_by(individual_id) %>% 
  arrange(timestamp) %>% 
  mutate(rn = row_number(),                                                 # create a row number
         grp = 1,                                                           # create a constant group # for annotations
         lag_lon = dplyr::lag(lon, 1),                                      # get lat and lon from last step
         lag_lat = dplyr::lag(lat, 1),                                      # ...
         sl = distGeo(cbind(lon,lat), cbind(lag_lon, lag_lat)),             # step dist
         dt = as.numeric(difftime(timestamp, dplyr::lag(timestamp, 1)), units='secs'), # time diff (secs)
         v = sl/dt, # velocity (m/s)
         bearing = bearing(cbind(lon,lat), cbind(lag_lon, lag_lat)),        # step bearing
         ta = 180-abs(180 - abs(bearing - dplyr::lag(bearing, 1)) %% 360),  # turn angle
         doy = yday(timestamp),                                             # day of year
         yr = year(timestamp),                                              # year
         month = month(timestamp),                                          # month
         season = case_when(
           month == 1 ~ "Winter",
           month == 2 ~ NA_character_,
           month == 3 ~ "Spring",
           month == 4 ~ "Spring",
           month == 5 ~ NA_character_,
           month == 6 ~ "Summer",
           month == 7 ~ "Summer",
           month == 8 ~ NA_character_,
           month == 9 ~ "Fall",
           month == 10 ~ "Fall",
           month == 11 ~ NA_character_, 
           month == 12 ~ "Winter"
         ),
         # Make Feb the 1st month of the year so winter doesn't split across years
         animal_year = case_when(month < 2 ~ (yr - 1),
                                 TRUE ~ yr),
         species = case_when(                                               # make new species var
           taxon_canonical_name == "Gruidae" ~ "Grus nigricollis",            # change "Gruidae" entries to correct sp
           TRUE ~ taxon_canonical_name)                                       # keep the rest as is
  )

# calculate qunatile-based cutoffs
cuts <- evt1 %>% 
  # as.data.frame() %>% 
  group_by(individual_id) %>% 
  summarize(
    qta = quantile(ta, probs = 0.95, na.rm = T, names = F),
    qsl = quantile(sl, probs = 0.95, na.rm = T, names = F),
    qv  = quantile(v,  probs = 0.95, na.rm = T, names = F)
  ) 

# filter outliers
evt2 <- evt1 %>% 
  # just join the cutpoints back to the dataset
  left_join(cuts) %>% 
  # conservative outlier thresh, must be past 95% quant for either sl and TA
  filter(v < qv & ta < qta) %>%
  filter(gps_hdop < 5 | is.na(gps_hdop)) %>% 
  filter(gps_dop < 5 | is.na(gps_dop)) %>% 
  filter(horizontal_accuracy < 25 | is.na(horizontal_accuracy)) %>% 
  ungroup()

message("Data summary after outlier removal:")
(ind_count1 <- evt2 %>% 
  group_by(species) %>% 
  summarise(n_ind = n_distinct(individual_id),
            n_obs = n())
)
nobs1 <- nrow(evt2)
message(glue("{nobs1} unique events."))


#-- Check minimum inds per sp --#

message("Beginning species-level data sufficiency check...")

rm_sp <- ind_count1 %>% 
  filter(n_ind < .sp_cut) %>% 
  pull(species)

evt3 <- evt2 %>% 
  filter(species %notin% rm_sp)

message("Data summary after removing data poor species:")
(ind_count2 <- evt3 %>% 
    group_by(species) %>% 
    summarise(n_ind = n_distinct(individual_id),
              n_obs = n())
)
nobs2 <- nrow(evt3)
message(glue("{nobs2} unique events."))


#-- Check minimum duration per ind --#

message("Beginning individual-level data sufficiency check...")

dur_summ <- evt3 %>% 
  group_by(individual_id) %>% 
  summarize(dur = difftime(max(timestamp), min(timestamp), units = "days"))

rm_ind2 <- dur_summ %>% 
  filter(dur < .min_dur) %>% 
  pull(individual_id)

out <- evt3 %>% 
  filter(individual_id %notin% rm_ind2)

message("Data summary after removing data poor individuals:")
(ind_count3 <- out %>% 
    group_by(species) %>% 
    summarise(n_ind = n_distinct(individual_id),
              n_obs = n())
)
nobs3 <- nrow(out)
message(glue("{nobs3} unique events."))


#---- Finalize script ----#

message("Writing results back to database...")

# write table back to db
dbWriteTable(conn = db, name = .out_table, value = out, append = FALSE, overwrite = T)

# disconnect from db
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
