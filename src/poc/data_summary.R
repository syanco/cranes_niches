# Interactive script to explore and summarize data and to inform requisite cleaning steps.

#---- Input Parameters ----#
if(interactive()) {
  # library(here)
  
  .wd <- '~/projects/dynamic_niches'
  .dbPF <- file.path(.wd,'data/anno_move.db')
  
} else {
  library(docopt)
  # library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  
}

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
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

`%notin%` <- Negate(`%in%`)

.anno <- "anno_join_2022-11-16"

db <- dbConnect(RSQLite::SQLite(), .dbPF)

dbListTables(db)

evttb <- tbl(db, .anno) %>%  collect()
indtb <- tbl(db, "individual") %>%  collect()

spp_key <- indtb %>% 
  select(taxon_canonical_name, individual_id)

# Individuals per species
indtb %>% 
  group_by(taxon_canonical_name) %>% 
  summarise(n_distinct(individual_id))


# What are the NA species??
nas <- indtb %>% 
  filter(is.na(taxon_canonical_name)) %>% 
  collect()

Grudiae <- indtb %>% 
  filter(taxon_canonical_name == "Gruidae") %>% 
  collect()

durations <- evttb %>% 
  mutate(doy = yday(timestamp),
         yr = year(timestamp),
         month = month(timestamp)) %>% 
  group_by(individual_id) %>% 
  summarize(tot_dur = difftime(max(timestamp), min(timestamp), units = "days")) %>% 
  left_join(spp_key)

durations %>% 
  group_by(taxon_canonical_name) %>% 
  summarize(n = n(),
            mean(tot_dur),
            sum(tot_dur > 365))
            