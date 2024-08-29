#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script compares seasonasl niches within an individual
# See project documentation for details about anticipated directory structure.

'
Compare seasonal niches within an individual.
Anticipates a cleaned event table in a mosey-formatted db.

Usage:
within_ind_across_seasons.r <db> <table> <out>
within_ind_across_seasons.r (-h | --help)

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
  .out <- file.path(.wd, "out/wi_ind_ax_seas.csv")
  
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
list.files(file.path(.wd,'src/funs/auto'), full.names=TRUE) %>%
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

# loop over habitat comparisons
counter <- 0
hab_out <- list()
for(i in 1:length(inds)){ # test with i = 4
  
  message(glue("Starting individual {inds[i]}..."))
  
  dat <- evt0 %>% 
    filter(individual_id == inds[i])
  
  yrs <- dat %>% 
    distinct(animal_year) %>%
    pull(animal_year)
  
  combos_out <- list()
  for(y in 1:length(yrs)){
    
    season_out <- list()
    
    seasons_yr <- dat %>% 
      pull(season) %>% 
      unique() %>% 
      na.omit()
    
    if(length(seasons_yr) > 1){
      combos <- combn(seasons, 2)
      
      
      for(s in 1:ncol(combos)){
        
        if(combos[1,s] %in% seasons_yr & combos[2,s] %in% seasons_yr){
          
          d1 <- dat %>% 
            filter(animal_year == yrs[y],
                   season == combos[1,s]) %>% 
            select(!!hab_vars) %>% 
            filter(complete.cases(.))
          
          d2 <- dat %>% 
            filter(animal_year == yrs[y],
                   season == combos[2,s]) %>% 
            select(!!hab_vars) %>% 
            filter(complete.cases(.))
          
          if(nrow(d1) > 0 & nrow(d2) > 0){
            season_out[[s]] <- tryCatch({
              MVNH_dissimilarity(d1, d2) %>% 
                unlist() %>% 
                t() %>% 
                data.frame() %>% 
                mutate(individual_id = inds[i],
                       season1 = combos[1,s],
                       season2 = combos[2,s],
                       species = unique(dat$species),
                       year = yrs[y])
            }, error=function(e){})
            
            # counter <-ifelse(is.null(season_out[[s]]), counter + 1, counter)
            
          } else{ #fi (d1 and d2 have data)
            message(glue("Not enough data to compare {combos[1,s]} and {combos[2,s]} for individual {inds[i]} in {yrs[y]}, moving on..."))
          } # else
        } else {NULL} # fi (combos in season_yr)
      } #s
      
      if(length(season_out) > 0){
        combos_out[[y]] <- do.call("rbind", season_out) 
      } #fi
    }
  } # y
  
  if(length(combos_out) > 0){
    hab_out[[i]] <- do.call("rbind", combos_out) 
  } # fi
  
} #i

hab_out_df <- do.call("rbind", hab_out) %>% 
  mutate(niche_type = "habitat")

message(glue("{nrow(hab_out_df)} habitat comparisons made..."))

# loop over habitat comparisons
counter <- 0
grin_out <- list()
for(i in 1:length(inds)){ # test with i = 4
  
  message(glue("Starting individual {inds[i]}..."))
  
  dat <- evt0 %>% 
    filter(individual_id == inds[i])
  
  yrs <- dat %>% 
    distinct(animal_year) %>%
    pull(animal_year)
  
  combos_out <- list()
  for(y in 1:length(yrs)){
    
    season_out <- list()
    
    seasons_yr <- dat %>% 
      pull(season) %>% 
      unique() %>% 
      na.omit()
    if(length(seasons_yr) > 1){
      combos <- combn(seasons, 2)
      
      
      for(s in 1:ncol(combos)){
        
        if(combos[1,s] %in% seasons_yr & combos[2,s] %in% seasons_yr){
          
          d1 <- dat %>% 
            filter(animal_year == yrs[y],
                   season == combos[1,s]) %>% 
            select(!!grin_vars) %>% 
            filter(complete.cases(.))
          
          d2 <- dat %>% 
            filter(animal_year == yrs[y],
                   season == combos[2,s]) %>% 
            select(!!grin_vars) %>% 
            filter(complete.cases(.))
          
          if(nrow(d1) > 0 & nrow(d2) > 0){
            season_out[[s]] <- tryCatch({
              MVNH_dissimilarity(d1, d2) %>% 
                unlist() %>% 
                t() %>% 
                data.frame() %>% 
                mutate(individual_id = inds[i],
                       season1 = combos[1,s],
                       season2 = combos[2,s],
                       species = unique(dat$species),
                       year = yrs[y])
            }, error=function(e){})
            
            # counter <-ifelse(is.null(season_out[[s]]), counter + 1, counter)
            
          } else{ #fi (d1 and d2 have data)
            message(glue("Not enough data to compare {combos[1,s]} and {combos[2,s]} for individual {inds[i]} in {yrs[y]}, moving on..."))
          } # else
        } else {NULL} # fi (combos in season_yr)
      } #s
      
      if(length(season_out) > 0){
        combos_out[[y]] <- do.call("rbind", season_out) 
      } #fi
    }
  } # y
  
  if(length(combos_out) > 0){
    grin_out[[i]] <- do.call("rbind", combos_out) 
  } # fi
  
} #i

grin_out_df <- do.call("rbind", grin_out) %>% 
  mutate(niche_type = "grinellian")

message(glue("{nrow(grin_out_df)} Grinellian comparisons made..."))

out <- bind_rows(hab_out_df, grin_out_df)


#---- Finalize script ----#

message(glue("Writing results to {.out}..."))

# write table back to db
write_csv(out, .out, append = FALSE)

# disconnect from db
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
