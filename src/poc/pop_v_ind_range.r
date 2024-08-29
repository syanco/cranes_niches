###########################
#                         #
# Pop vs Ind Niche Dissim #
#                         #
###########################

# Scratch script to test out ideas of comparing population niche dissim with individual.

# LIBRARIES

library(MVNH)
library(jsonlite)
library(httr)
library(RSQLite)
library(DBI)
library(tidyverse)
library(lubridate)
library(glue)
library(patchwork)
library(ggplot2)

# FUNCTIONS

get_scenarios <- function () {
  resp <- httr::GET(paste0(BASE_URL, '/list/scenarios'))
  jsonlite::flatten(jsonlite::fromJSON(httr::content(resp, "text", encoding = enc), simplifyVector = T))
}

get_species_scenario <- function (sciname, product, variable, s_buff, t_buff, limit=50000) {
  resp <- httr::POST(paste0(BASE_URL, '/species/metrics/scatter'), body = list(
    mode = 'temporal',
    scientificname = sciname,
    limit= limit,
    yaxis = list(
      product = product,
      variable = variable,
      temporal = t_buff, 
      spatial = s_buff
    )
  ), encode = "json",)
  jsonlite::flatten(jsonlite::fromJSON(httr::content(resp, "text", encoding = enc))$rows)
}

`%notin%` <- Negate(`%in%`)

source("src/funs/seg2anno.r")

# INITS

BASE_URL <- 'https://stoat-otf.azurewebsites.net'
enc = "UTF-8"
.dbPF <- file.path('data/anno_move.db')
.anno <- "anno_join_2022-11-16"
.ctfs <- file.path("ctfs/segmentations")


# Pull STOAT Pre-Annotations
#
# View possible annotations
# get_scenarios()
envs <- read_csv("ctfs/anno_vars.csv") %>% 
  filter(run == 1)

pop_env <- list()
for(i in 1:nrow(envs)){
  pop_env[[i]] <- get_species_scenario('Anthropoides virgo', envs$product[i], 
                                       envs$variable_code[i], envs$s_buff[i], 
                                       envs$t_buff[i]) %>% 
    mutate(variable = envs$variable[i])
}

pop_total <- do.call("rbind", pop_env) %>% 
  distinct(eventdate, latitude, longitude, variable, .keep_all = T) %>% 
  pivot_wider(names_from = variable, values_from = value) %>% 
  mutate(doy = yday(eventdate))


# Read in cranes

#get segmentation filenames
segs <- list.files(.ctfs, pattern = "*.csv")

db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)

# Get target inds from one sp
inds <- tbl(db, "individual") %>%
  filter(taxon_canonical_name == "Anthropoides virgo") %>% 
  pull(individual_id)

# create df of annotations
anno0 <- tbl(db, .anno) %>% 
  filter(individual_id %in% inds) %>%     
  filter(is.na(`ground_speed`) | `ground_speed`<10) %>% 
  filter(is.na(gps_hdop) | gps_hdop < 5) %>% 
  collect()

# add season annotations based on the manual segemntation cutpoints
anno_seasons <- cuts2anno(df = anno0, segs = segs, inds = inds) %>% 
  #only retain good segmentations
  filter(checksum == 1)

anno_res <- anno_seasons %>% 
  filter(winter == 1 | summer == 1)

# get individual dissimilarities

# init empty df
ind_out <- data.frame()

# get reduced list of inds (only those who made it through above)
ind2 <- unique(anno_res$individual_id)

for (i in 1:length(ind2)){
  summ <- anno_res %>% 
    filter(individual_id==ind2[i] & summer == 1) %>% 
    select(`value_derived:evi`, value_slope, value_modification, value_altitude) %>% 
    na.omit()
  wint <- anno_res %>% 
    filter(individual_id==ind2[i] & winter == 1)%>% 
    select(`value_derived:evi`, value_slope, value_modification, value_altitude) %>% 
    na.omit()
  
  if(nrow(summ) < 1 | nrow(wint) < 1){
    next
  } else {
    
    out <- data.frame(ind = ind2[i],
                      ind_dis = as.numeric(MVNH_dissimilarity(summ, wint)$Bhattacharyya_distance['total']),
                      md = as.numeric(MVNH_dissimilarity(summ, wint)$Mahalanobis_distance['total']),
                      dr = as.numeric(MVNH_dissimilarity(summ, wint)$Determinant_ratio['total']),
                      mean_summ = mean(summ$`value_derived:evi`),
                      mean_wint= mean(wint$`value_derived:evi`)
    )
    
    ind_out <- rbind(ind_out, out)
    
  }
}

# Get niche breadths

ind_out_breadth <- data.frame()

for (i in 1:length(ind2)){
  summ <- anno_res %>% 
    filter(individual_id==ind2[i] & summer == 1) %>% 
    mutate(evi_sc = scale(`value_derived:evi`),
           # lst_sc = scale(value_lst_day_1km),
           slope_sc = scale(value_slope),
           altitude_sc = scale(value_altitude),
           mod_sc = scale(value_modification)) %>% 
    select(evi_sc, slope_sc, altitude_sc, mod_sc) %>% 
    na.omit()
  wint <- anno_res %>% 
    filter(individual_id==ind2[i] & winter == 1)%>% 
    mutate(evi_sc = scale(`value_derived:evi`),
           # lst_sc = scale(value_lst_day_1km),
           slope_sc = scale(value_slope),
           altitude_sc = scale(value_altitude),
           mod_sc = scale(value_modification)) %>% 
    select(evi_sc, slope_sc, altitude_sc, mod_sc) %>% 
    na.omit()
  
  if(nrow(summ) < 1 | nrow(wint) < 1){
    next
  } else {
    
    out_summ <- data.frame(ind = ind2[i],
                           season = "summer",
                           total = as.numeric(MVNH_det(summ, log = T)['total']),
                           cor = as.numeric(MVNH_det(summ)['cor']),
                           evi = as.numeric(MVNH_det(summ)['evi_sc']),
                           # lst = as.numeric(MVNH_det(summ)['lst_sc']),
                           slope = as.numeric(MVNH_det(summ)['slope_sc']),
                           altitude = as.numeric(MVNH_det(summ)['altitude_sc']),
                           mod = as.numeric(MVNH_det(summ)['mod_sc'])
    )
    out_wint <- data.frame(ind = ind2[i],
                           season = "winter",
                           total = as.numeric(MVNH_det(wint, log = T)['total']),
                           cor = as.numeric(MVNH_det(wint)['cor']),
                           evi = as.numeric(MVNH_det(wint)['evi_sc']),
                           # lst = as.numeric(MVNH_det(wint)['lst_sc']),
                           slope = as.numeric(MVNH_det(wint)['slope_sc']),
                           altitude = as.numeric(MVNH_det(wint)['altitude_sc']),
                           mod = as.numeric(MVNH_det(wint)['mod_sc'])
    )
    
    out <- rbind(out_summ, out_wint)
    ind_out_breadth <- rbind(ind_out_breadth, out)
    
  }
}


# Get phenology
phen <- data.frame()
for(j in 1:length(ind2)){
  # grab the segmentation file by matching filename to ind
  if(ind2[j] %in% str_extract(segs, "[0-9]+")){ #if the file exists...
    seg_temp <- read_csv(file.path( # ...read it in...
      .ctfs,segs[which(str_extract(segs, "[0-9]+") == ind2[j])[1]] # ...matching individual_id
    )) %>%
      arrange(Date) %>% # sort by date
      #remove stopovers
      filter(Status %in% c("Start Fall", "End Fall", "Start Spring", "End Spring")) %>% 
      mutate(yr = year(Date)) # add year col
    
  } else {
    # otherwise skip to next individual
    message(glue("No segmentation file found for individual {inds[j]}, moving on!"))
    next
  }
  phen <- rbind(phen, seg_temp)
}

# Add day of year
phen <- phen %>% 
  mutate(doy = yday(Date))

# calc median transition days
mssp <- median(phen %>% filter(Status == "Start Spring") %>% pull(doy))
mesp <- median(phen %>% filter(Status == "End Spring") %>% pull(doy))
msfa <- median(phen %>% filter(Status == "Start Fall") %>% pull(doy))
mefa <- median(phen %>% filter(Status == "End Fall") %>% pull(doy))




# get pop-scale seasonal data
gbif_summ <- pop_total %>% 
  filter(doy > mesp & doy < msfa)
gbif_wint <- pop_total %>% 
  filter(doy < mssp | doy > mefa)

#TODO:  make this work for more than 1 var...
# # calc pop means
# pop_stats <- data.frame(season = c("summer", "winter"),
#                         mean = c(mean(gbif_summ$value), mean(gbif_wint$value)),
#                         variance = c(var(gbif_summ$value), var(gbif_wint$value)),
#                         stdv = c(sd(gbif_summ$value), sd(gbif_wint$value)))
# pop_mean_summ <- mean(gbif_summ$value)
# pop_mean_wint <- mean(gbif_wint$value)

# Calc pop dissimilarity
pop_diss <- MVNH_dissimilarity(gbif_summ %>% select(modification, `derived:evi`, altitude, slope) %>% na.omit(),
                               gbif_wint %>% select(modification, `derived:evi`, altitude, slope) %>% na.omit())$Bhattacharyya_distance['total']
pop_breadth_summ <- MVNH_det(gbif_summ %>% select(modification, `derived:evi`, altitude, slope) %>% na.omit(), log = T)
pop_breadth_wint <- MVNH_det(gbif_wint %>% select(modification, `derived:evi`, altitude, slope) %>% na.omit(), log = T)

# Plots
# \
# DISSIM
ggplot()+
  geom_density(data = ind_out, aes(x = ind_dis))+
  geom_vline(aes(xintercept = pop_diss), color = "darkred")+
  theme_minimal()

ind_out_long <- ind_out %>% 
  pivot_longer(c(mean_summ, mean_wint), names_to = "case", values_to = "mean") %>% 
  mutate(season = case_when(case == "mean_summ" ~ "summer",
                            case == "mean_wint" ~ "winter"),
         ind = as.factor(ind))

ggplot(pop_stats)+
  geom_point(aes(x=season, y=mean)) +
  geom_errorbar(aes(x= season, ymin = mean-stdv, ymax = mean+stdv), width = 0.1)+
  theme_minimal()+
  ylim(0,0.5)+
  theme(legend.position = "none")

ggplot(ind_out_long)+
  # geom_point(aes(x=season, y=mean, color = ind)) +
  geom_path(aes(x= season, y = mean, group = ind, color = ind))+
  ylim(0,0.5)+
  theme_minimal()+
  theme(legend.position = "none")

#  BREADTH
ind_out_breadth %>%
  pivot_wider(id_cols = ind, names_from = season, values_from = total) %>% 
  mutate(diff = summer-winter,
         ind = as.factor(ind)) %>% 
  ggplot(aes(y = reorder(ind, diff)))+
  geom_point(aes(x = summer), color = "red")+
  geom_point(aes(x = winter), color = "blue") +
  geom_segment(aes(x= summer, xend = winter, yend = reorder(ind, diff)))+
  theme_minimal()+
  ylab("Individual")+
  xlab("log(Niche Breadth)")+
  ggtitle("Seasonal niche breadth for Anthropoides virgo")+
  theme(axis.text.y = element_blank())

ind_out_breadth %>% 
  filter(season == "summer") %>% 
  ggplot() +
  geom_density(aes(x=total)) +
  geom_vline(aes(xintercept = pop_breadth_summ["total"]), color = "darkred")+
  theme_minimal()

ind_out_breadth %>% 
  filter(season == "winter") %>% 
  ggplot() +
  geom_density(aes(x=total)) +
  geom_vline(aes(xintercept = pop_breadth_wint["total"]), color = "darkred")+
  theme_minimal()
