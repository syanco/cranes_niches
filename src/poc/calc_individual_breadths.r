#!/usr/bin/env Rscript 

# Script to loop through species and individuals to calculate niche dissimilarity 
# between seasonal ranges.  Relies on segmentations having been pre-recorded and 
# stored in individual-specific control files. Returns a single data frame with 
# individual-specific dissimilarities.
# 
# TODO:  Update docopt to accommodate inputs used. (after breezy header is final)


'
calc_individual_dissims

Usage:
script_template <dat> <out> [--db=<db>] [--seed=<seed>] [-b] [-t]
script_template (-h | --help)

Control files:
  ctfs/individual.csv

Parameters:
  dat: path to input csv file. 
  out: path to output directory.

Options:
-h --help     Show this screen.
-v --version     Show version.
-d --db=<db> Path to movement database. Defaults to <wd>/data/move.db
-s --seed=<seed>  Random seed. Defaults to 5326 if not passed
-b --rollback   Rollback transaction if set to true.
-t --test         Indicates script is a test run, will not save output parameters or commit to git
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  
  .wd <- '~/projects/dynamic_niches'
  
  # .datPF <- file.path(.wd,'data/dat.csv')
  .outP <- file.path(.wd,'out')
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .segsP <- file.path(.wd, 'ctfs/segmentations')
  .sppPF <- file.path(.wd, 'ctfs/species.csv')
  .envsPF <- file.path(.wd, 'ctfs/anno_vars.csv')
  .nc <- 2
  
} else {
  library(docopt)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  # .datPF <- makePath(ag$dat)
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .segsP <- makePath(ag$segs)
  .sppPF <- makePath(ag$spp)
  .nc <- ag$nc
  
}

#---- Initialize Environment ----#
message("Initializing environment...")
t0 <- Sys.time()

source(file.path(.wd, 'src/init/startup.r'))
source(file.path(.wd, 'src/funs/seg2anno.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(MVNH)
    library(foreach)
    library(doMC)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd, 'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

# define vec of environmental vars
envs <- c("value_evi", "value_lst_day_1km", "value_modification")

registerDoMC(.nc)

#---- Load control files ----#
message("Loading control files...")

# load species ctf
spp <- read.csv(.sppPF) %>% filter(run == 1)
segs <- list.files(.segsP, pattern = "*.csv")

# load environments ctf
# envs <- read.csv(.envsPF) %>% filter(run == 1)

#---- Initialize database ----#
message("Connecting to database...")
invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF)
invisible(assert_that(length(dbListTables(db))>0))

#---- Load data ----#

# Create df of individuals from species in ctf 
inds <- tbl(db, "individual") %>% 
  collect() %>% 
  filter(taxon_canonical_name %in% spp$species_name)

#---- Perform analysis ----#
# create df of annotations
anno0 <- tbl(db, "anno_join_2022-11-16") %>% 
  filter(individual_id %in% !!inds$individual_id) %>% # filter to tagrget ind     
  filter(is.na(`ground_speed`) | `ground_speed`<10) %>% # rm speedy obs...
  filter(is.na(gps_hdop) | gps_hdop < 5) %>% # rm high hdop
  mutate(value_evi = `value_derived:evi`) %>% 
  collect() # collect into mem

# add season annotations based on the manual segemntation cutpoints
anno_seasons <- cuts2anno(df = anno0, segs = segs, inds = inds$individual_id) %>% 
  #only retain good segmentations
  filter(checksum == 1)

# filter out migratory periods
anno_res <- anno_seasons %>% 
  filter(winter == 1 | summer == 1)

ind_out <- list()

# Toggle `%do%` to `%dopar%` for HPC, %do% for local
foreach(i = 1:nrow(inds), .errorhandling = "pass", .inorder = F) %do% {
  # foreach(i = 1:nrow(inds), .errorhandling = "pass", .inorder = F) %dopar% {
  
  # ind_out_breadth <- data.frame()
  

    summ <- anno_res %>% 
      filter(individual_id==inds$individual_id[i] & summer == 1) %>% 
      mutate(evi_sc = scale(`value_derived:evi`),
             lst_sc = scale(value_lst_day_1km),
             # slope_sc = scale(value_slope),
             # altitude_sc = scale(value_altitude),
             humod_sc = scale(value_modification)) %>% 
      select(evi_sc, 
             lst_sc, 
             # slope_sc, 
             # altitude_sc,
             humod_sc) %>% 
      na.omit()
    
    wint <- anno_res %>% 
      filter(individual_id==inds$individual_id[i] & winter == 1)%>% 
      mutate(evi_sc = scale(`value_derived:evi`),
             lst_sc = scale(value_lst_day_1km),
             # slope_sc = scale(value_slope),
             # altitude_sc = scale(value_altitude),
             humod_sc = scale(value_modification)) %>% 
      select(evi_sc, 
             lst_sc, 
             # slope_sc, 
             # altitude_sc,
             humod_sc ) %>%
      na.omit()
    
    if(nrow(summ) < 1 | nrow(wint) < 1){
      next
    } else {
      
      out_summ <- data.frame(ind = inds$individual_id[i],
                             species = inds$taxon_canonical_name[i],
                             season = "summer",
                             total = as.numeric(MVNH_det(summ, log = T)['total']),
                             cor = as.numeric(MVNH_det(summ)['cor']),
                             evi = as.numeric(MVNH_det(summ)['evi_sc']),
                             lst = as.numeric(MVNH_det(summ)['lst_sc']),
                             # slope = as.numeric(MVNH_det(summ)['slope_sc']),
                             # altitude = as.numeric(MVNH_det(summ)['altitude_sc']),
                            humod = as.numeric(MVNH_det(summ)['humod_sc'])
      )
      out_wint <- data.frame(ind = inds$individual_id[i],
                             species = inds$taxon_canonical_name[i],
                             season = "winter",
                             total = as.numeric(MVNH_det(wint, log = T)['total']),
                             cor = as.numeric(MVNH_det(wint)['cor']),
                             evi = as.numeric(MVNH_det(wint)['evi_sc']),
                             lst = as.numeric(MVNH_det(wint)['lst_sc']),
                             # slope = as.numeric(MVNH_det(wint)['slope_sc']),
                             # altitude = as.numeric(MVNH_det(wint)['altitude_sc']),
                             humod = as.numeric(MVNH_det(summ)['humod_sc'])
      )
      
      ind_out[[i]] <- rbind(out_summ, out_wint)

    }

}

out <- do.call("rbind", ind_out)


#---- Save output ---#
message(glue('Saving to {.outP}'))

write_csv(out, path = file.path(.outP, "niche_breadth.csv"))

#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
