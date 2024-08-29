
#---- Input Parameters ----#
if(interactive()) {
  
  .wd <- '~/projects/dynamic_niches'
  
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
    library(rjson)
    library(tictoc)
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

X <- lc[1:10000,]

tic()
expandESACCI(X)
toc()

tic()
expandESACCI2(X)
toc()

tic()
expandESACCIgpt(X)
toc()

tic()
# set up parallel backend
plan(multisession, workers = 8)
expandESACCI_parallel(X)
toc()