#!/usr/bin/env Rscript 

# Script to unpack the JSON habitat proportion var

'
unpack_habitat.r

Usage:
unpack_habitat.r <db> <var>
unpack_habitat.r (-h | --help)

Parameters:
  db: path to mosey database. 
  var: target variable name.

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
  
  .wd <- getwd()
  
  .dbPF <- file.path(.wd,'data/anno_move.db')

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

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(future.apply)
    library(future)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd, 'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

# Source expandESACCI function
source(file = "src/funs/expandESACCI.r")

#---- Initialize database ----#
message("Connecting to database...")
invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF)
invisible(assert_that(length(dbListTables(db))>0))

#---- Load data ----#

lc <- tbl(db, "cranes_esacci_300_1000") %>% 
  filter(s_buff == 300) %>% 
  collect()


# # Create df of individuals from species in ctf 
# evt0 <- tbl(db, "event_clean") %>% 
#   collect()


#---- Perform analysis ----#
# lc_test <- lc %>% filter(event_id == 1154917924 | event_id ==1154742303)

# Set up future plan
plan(multisession, workers = 8)

# expand ESA CCI habitat
lc_exp <- expandESACCI_par(lc) 
# %>% 
  # rowwise() %>%
  # dplyr::mutate(total = rowSums(., na.rm = T)) # calc total pixels

# put df back together
dat_new <- cbind(lc, lc_exp)

#---- Save output ---#
message(glue('Saving to table back to database...'))

dbWriteTable(db, "cranes_esacci_300_1000_expanded", dat_new, overwrite = T)
#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
