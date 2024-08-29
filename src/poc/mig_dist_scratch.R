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
    library(sf)
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

inds <- evt0 %>%  
  pull(individual_id) %>% 
  unique()

# get individual centroid-centroid diffs
dist_out <- tibble()
x <- evt0 %>% 
  filter(individual_id==inds[i])

for (i in 1:length(unique(evt0$individual_id))){
  summ <- evt0 %>% 
    filter(individual_id==inds[i] & season == "Summer") 
  
  if(nrow(summ) > 0) {
    summ_cent <- st_as_sf(summ, coords = c("lon", "lat")) %>% 
      st_union() %>% 
      st_centroid()
  }
  
  wint <- evt0 %>% 
    filter(individual_id==inds[i] & season == "Winter")
  
  if(nrow(wint) > 0) {
    wint_cent <- st_as_sf(wint, coords = c("lon", "lat")) %>% 
      st_union() %>% 
      st_centroid()
  }
  
  if(exists("wint_cent") & exists("summ_cent")){
    dist <- st_distance(summ_cent, wint_cent)
    out <- data.frame(ind = inds[i],
                      dist = dist,
                      wint_cent = wint_cent,
                      summ_cent = summ_cent
    )
    
    dist_out <- rbind(dist_out, out)
  } else {
    message(glue("Data for {i} do not contain position information for both seasons, distance not calculated."))
  }
  
  #clear env so that control flow doesn't respond to last loop iteration
  if(exists("wint_cent")){rm(list = "wint_cent")}
  if(exists("summ_cent")){rm(list = "summ_cent")}
  if(exists("wint")){rm(list = "wint")}
  if(exists("summ")){rm(list = "summ")}
  
  # rm(list = "wint_cent", "wint_cent", "wint", "summ")
  
}


seas <- read_csv("out/wi_ind_ax_seas.csv") %>% 
  filter(
    # Bhattacharyya_distance.total < 100,
    species != "Anthropoides paradiseus") %>% 
  unite(comp, season1:season2, remove = F)

seas %>% 
  group_by(species) %>% 
  summarize(n = n(),
            inds = n_distinct(individual_id))
evt0 %>% 
  group_by(species) %>% 
  summarize(n = n(),
            inds = n_distinct(individual_id))

comb <- seas %>% 
  inner_join(dist_out, by = c("individual_id" = "ind"))
comb %>% 
  group_by(species) %>% 
  summarize(n = n(),
            inds = n_distinct(individual_id))
comb %>% 
  filter(comp == "Winter_Summer") %>% 
ggplot() +
  geom_point(aes(x=dist, y = Bhattacharyya_distance.total, color = niche_type)) +
  facet_wrap(~species)


comb %>% 
  filter(comp == "Winter_Summer") %>% 
  pull(Bhattacharyya_distance.total) %>% 
  summary()
