
# library(here)
if(interactive()) {
  .wd <- getwd()
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .table <- "niche_prep"
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

ind0 <- tbl(db, "individual") %>% 
  collect()

study0 <- tbl(db, "study") %>% 
  collect()

deploy0 <- tbl(db, "deployment") %>% 
  collect()

tag0 <- tbl(db, "tag") %>% 
  collect()

# Process to shareable format

evt_out <- evt0 %>% 
  filter(taxon_canonical_name != "Balearica pavonina") %>% 
  mutate(prox2water = case_when(dist2water == 0 ~ 1, # transform prox to water var
                                TRUE ~ 1/dist2water),
         prox2water_scale = as.numeric(scale(prox2water))) %>% 
  select(event_id, individual_id, lon, lat, timestamp, tag_id, sensor_type_id,
         ground_speed, gps_speed_accuracy_estimate, gps_dop, gps_hdop, gps_vdop,
         gps_satellite_count, horizontal_accuracy, time_to_fix, fix_type, 
         taxon_canonical_name, study_id, prop_crop_scale, prox2water_scale, 
         evi_scale, lst_scale) #reduce columns

# Make vectors of inds and studies to keep
inds_keep <- unique(evt_out$individual_id)
stds_keep <- unique(evt_out$study_id)

# Reduce ind table to included animals
ind_out <- ind0 %>% 
  filter(individual_id %in% inds_keep)

# Reduce study table to included studies
std_out <- study0 %>% 
  filter(study_id %in% stds_keep)

# Reduce deployments
deploy_out <- deploy0 %>% 
  filter(individual_id %in% inds_keep)

# Get used tag_id
tags_keep <- unique(deploy_out$tag_id)

tag_out <- tag0 %>% 
  filter(tag_id %in% tags_keep)

# Write Out

# dir.create("./data/data_to_share")
write_csv(evt_out, "./data/data_to_share/cranes_event.csv")
write_csv(std_out, "./data/data_to_share/cranes_study.csv")
write_csv(ind_out, "./data/data_to_share/cranes_individual.csv")
write_csv(deploy_out, "./data/data_to_share/cranes_deployment.csv")
write_csv(tag_out, "./data/data_to_share/cranes_tag.csv")
