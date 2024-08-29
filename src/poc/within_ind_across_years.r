#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script compares season-matched niches within an individual across years.
# See project documentation for details about anticipated directory structure.

'
Compare season-matched niches within an individual across years.
Anticipates a cleaned event table in a mosey-formatted db.

Usage:
clean_movement.r <db> <table> <out>
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
    library(MVNH)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

# vector of seasons we'll need later
seasons <- c("Winter", "Spring", "Summer", "Fall")

#---- Load control files ----#

hab_vars <- c("dist2water_scale", 
              # "pfor_h_scale",
              # "evi_scale", 
              "ghm_scale", 
              # "altitude_scale", 
              "prop_crop_scale"
)

grin_vars <- c(
  # "dist2water_scale", 
  # "pfor_h_scale",
  "evi_scale",
  # "ghm_scale", 
  "altitude_scale"
  # "prop_crop_scale"
)


#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

# create event table
evt0 <- tbl(db, .table) %>% 
  collect()

#---- Perform Analyses ----#

inds <- evt0 %>% 
  distinct(individual_id) %>% 
  pull(individual_id)

# test with i = 1

counter <- 0

# Collect Habitat Niche Metrics
hab_out <- list()
for(i in 1:length(inds)){
  
  message(glue("Starting individual {inds[i]}..."))
  
  dat <- evt0 %>% 
    filter(individual_id == inds[i])
  
  years <- dat %>% 
    distinct(animal_year) %>%
    pull(animal_year) 
  
  if(length(years) < 2) {
    message("Only a single year in data, skipping...")
    next
  } else {
  combos <- years%>% 
    combn(2)
  
  message("Potential comparisons: ")
  print(combos)             
  
  combos_out <- list()
  for(y in 1:ncol(combos)){
    
    # d1 <- dat %>% 
    #   filter(animal_year == combos[1,y]) %>% 
    #   select(lst, evi) %>% 
    #   filter(complete.cases(.))
    # 
    # d2 <- dat %>% 
    #   filter(animal_year == combos[2,y]) %>% 
    #   select(lst, evi) %>% 
    #   filter(complete.cases(.))
    # 
    # if(nrow(d1) > 0 & nrow(d2) > 0){
    #   annual <- MVNH_dissimilarity(d1, d2) %>% 
    #     unlist() %>% 
    #     t() %>% 
    #     data.frame() %>% 
    #     mutate(individual_id = inds[i],
    #            yr1 = combos[1,y],
    #            yr2 = combos[2,y],
    #            species = unique(dat$taxon_canonical_name),
    #            season = "annual")
    # } else{
    #   message(glue("Not enough data to compare years {combos[1,y]} and {combos[2,y]} for individual {inds[i]}, moving on..."))
    # } # else  
    
    # Do comparisons matched within seasons
    season_out <- list()
    
    for(s in 1:4){
      # tes s <- 2
      
      d1_hab <- dat %>% 
        filter(animal_year == combos[1,y],
               season == seasons[s]) %>% 
        select(!!hab_vars) %>% 
        filter(complete.cases(.))
      
      d2_hab <- dat %>% 
        filter(animal_year == combos[2,y],
               season == seasons[s]) %>% 
        select(!!hab_vars) %>% 
        filter(complete.cases(.))

      if(nrow(d1_hab) > 0 & nrow(d2_hab) > 0){
        tryCatch({
          season_out[[s]] <- MVNH_dissimilarity(d1_hab, d2_hab) %>% 
            unlist() %>% 
            t() %>% 
            data.frame() %>% 
            mutate(individual_id = inds[i],
                   yr1 = combos[1,y],
                   yr2 = combos[2,y],
                   species = unique(dat$species),
                   season = seasons[s])
          
          counter <- counter + 1
        }, error = function(e){cat(glue("ERROR CAUGHT: unspecified error in `MVNH_dissimilarity` for the comparison: {seasons[s]} {combos[1,y]} and {seasons[s]} {combos[2,y]} for individual {inds[i]}, moving on...", 
                                        "\n"))})
      } else{
        message(glue("Not enough data to compare {seasons[s]} {combos[1,y]} and {seasons[s]} {combos[2,y]} for individual {inds[i]}, moving on..."))
      } # else
      
    } #s
    
    if(length(season_out) > 0){
      combos_out[[y]] <- do.call("rbind", season_out) 
    } #fi
    
  } # y
  
  if(length(combos_out) > 0){
    hab_out[[i]] <- do.call("rbind", combos_out) 
  } # fi
  } 
} #i

hab_out_df <- do.call("rbind", hab_out) %>% 
  mutate(niche_type = "habitat")

message(glue("{counter} Habitat comparisons made..."))


# Collect Grinellian Niche Metrics
grin_out <- list()
counter <- 0
for(i in 1:length(inds)){
  
  message(glue("Starting individual {inds[i]}..."))
  
  dat <- evt0 %>% 
    filter(individual_id == inds[i])
  
  years <- dat %>% 
    distinct(animal_year) %>%
    pull(animal_year) 
  
  if(length(years) < 2) {
    message("Only a single year in data, skipping...")
    next
  } else {
    combos <- years%>% 
      combn(2)
    
  
  message("Potential comparisons: ")
  print(combos)             
  
  combos_out <- list()
  for(y in 1:ncol(combos)){
    
    # Do comparisons matched within seasons
    season_out <- list()
    
    for(s in 1:4){
      # tes s <- 2
      
      d1_grin <- dat %>% 
        filter(animal_year == combos[1,y],
               season == seasons[s]) %>% 
        select(!!grin_vars) %>% 
        filter(complete.cases(.))
      
      d2_grin <- dat %>% 
        filter(animal_year == combos[2,y],
               season == seasons[s]) %>% 
        select(!!grin_vars) %>% 
        filter(complete.cases(.))
      
      if(nrow(d1_grin) > 0 & nrow(d2_grin) > 0){
        tryCatch({
          season_out[[s]] <- MVNH_dissimilarity(d1_grin, d2_grin) %>% 
            unlist() %>% 
            t() %>% 
            data.frame() %>% 
            mutate(individual_id = inds[i],
                   yr1 = combos[1,y],
                   yr2 = combos[2,y],
                   species = unique(dat$species),
                   season = seasons[s])
          
          counter <- counter + 1
        }, error = function(e){cat(glue("ERROR CAUGHT: unspecified error in `MVNH_dissimilarity` for the comparison: {seasons[s]} {combos[1,y]} and {seasons[s]} {combos[2,y]} for individual {inds[i]}, moving on...", 
                                        "\n"))})
      } else{
        message(glue("Not enough data to compare {seasons[s]} {combos[1,y]} and {seasons[s]} {combos[2,y]} for individual {inds[i]}, moving on..."))
      } # else
      
    } #s
    
    if(length(season_out) > 0){
      combos_out[[y]] <- do.call("rbind", season_out) 
    } #fi
    
  } # y
  
  if(length(combos_out) > 0){
    grin_out[[i]] <- do.call("rbind", combos_out) 
  } # fi
  }
} #i

grin_out_df <- do.call("rbind", grin_out) %>% 
  mutate(niche_type = "grinellian")


message(glue("{counter} Grinellian comparisons made..."))

out <- bind_rows(hab_out_df, grin_out_df)

#---- Finalize script ----#

message(glue("Writing results to {.out}..."))

# write table back to db
write_csv(out, .out, append = FALSE)

# disconnect from db
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
