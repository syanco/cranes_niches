#!/usr/bin/env Rscript 
# chmod 744 add_unpacked.r #Use to make executable

# This script can be used to join annotations back to a mosey_db-style database
# and add then joined table back to the db.  it anticipates a control file that 
# selects which elements of the database should be included in the join actions.

# ==== Setup ====

'
Join the unpacked ESA CCI annotations to the other movement and environmental data data.

Usage:
add_unpacked.r <db> 
add_unpacked.r (-h | --help)

Parameters:
  dat: path to input csv file. 

Options:
-h --help     Show this screen.
-v --version     Show version.
-b --rollback   Rollback transaction if set to true.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- getwd()
  .rollback <- TRUE
  rd <- here::here
  
  .outPF <- file.path(.wd,'analysis/cranes_anno.csv')
  .dbPF <- file.path(.wd,'data/anno_move.db')
  .ctfPF <- file.path(.wd, "ctfs/anno_vars.csv")
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  .script <-  thisfile()
  .seed <- ag$seed
  .rollback <- as.logical(ag$rollback)
  rd <- is_rstudio_project$make_fix_file(.script)
  
  source(rd('src/funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .ctfPF <- makePath(ag$ctf)
}

#---- Initialize Environment ----#
t0 <- Sys.time()

source('src/init/startup.r')

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
  }))

# #Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Control Files ----#

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#

#get cleaned event table
evt <- tbl(db, "event_clean") %>% 
  mutate(event_id = as.numeric(event_id)) %>% 
  collect()

esacci <- tbl(db, "cranes_esacci_300_1000_expanded") %>% 
  collect() 

out <- evt %>% 
  # select(!`...1`) %>% 
  left_join(esacci)

#---- Save output ---#
message(glue('Saving to table back to database...'))

dbWriteTable(db, "event_clean_new", out, overwrite = T)
#---- Finalize script ----#
dbDisconnect(db)
message(glue('Script complete in {diffmin(t0)} minutes'))
