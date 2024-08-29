#---- Input Parameters ----#
if(interactive()) {
  # library(here)
  
  .wd <- getwd()
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
    library(ggfortify)
    library(zoo)
    library(patchwork)
    library(colorspace)
    library(tidyverse)
    library(rnaturalearth)
    library(sf)
    library(raster)
    library(spatial)
    library(ggnewscale)
    library(ggspatial)
    library(scales)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)


#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

# create event table
evt0 <- tbl(db, .table) %>% 
  collect() %>% 
  mutate(species_common = case_when(species == "Anthropoides virgo" ~ "Demoiselle Crane",
                                    species == "Grus grus" ~ "Common Crane",
                                    species == "Grus nigricollis" ~ "Black-necked Crane",
                                    species == "Grus vipio" ~ "White-naped Crane"))

# Species level phenology
phen_summ <- read_csv("out/phenology_summary.csv")

# Individual centroids data
# ind_cents <- read_csv("out/individual_seasonal_centroids.csv")
pal <- c("evi" = "#416000", "prox2water" = "#3B75A9", "lst" = "#FF5733", 
         "propcrop" = "#FFD700")
out_plots <- list()
sp_ls <- unique(evt0$species_common)[1:4]

for(i in 1:length(sp_ls)){
  sp_name <- sp_ls[i]
  
  #-- Env Plot --#
  
  # Extract species data
  #TODO: below is a frankenstein - probably some lines I could remove
  evt_sp <- evt0 %>% 
    filter(species_common == sp_name) %>% 
    mutate(grp = glue("{individual_id}_{yr}"),
           ind_f = as.factor(individual_id),
           ind_num = as.numeric(ind_f),
           wk = week(timestamp),
           prox2water = case_when(dist2water == 0 ~ 1,
                                  TRUE ~ 1/dist2water),
           prox2water_scale = scale(prox2water)) %>% 
    group_by(grp, wk) %>% 
    summarize(dist2water = mean(dist2water, na.rm = T), 
              prox2water = mean(prox2water, na.rm = T),
              evi = mean(evi, na.rm = T), 
              propcrop = mean(prop_crop, na.rm = T), 
              lst = mean(lst, na.rm = T),
              # v_scale = mean(v_scale, na.rm = T),
              # rad_scale = mean(rad_scale, na.rm = T),
              individual_id = individual_id[1],
              yr = yr[1]) %>% 
    group_by(wk) %>% # average individuals (for now...)
    summarize(
      # dist2water_scale = median(dist2water_scale, na.rm = T), 
      prox2water = median(prox2water, na.rm = T),
      evi = median(evi, na.rm = T), 
      propcrop = median(propcrop, na.rm = T),
      lst = median(lst, na.rm = T)
      # v_scale = mean(v_scale, na.rm = T),
      # rad_scale = mean(rad_scale, na.rm = T)) 
    ) %>% ungroup()

  
  # create empirical dist functions
  evi_f <- ecdf(evt_sp$evi)
  lst_f <- ecdf(evt_sp$lst)
  prop_crop_f <- ecdf(evt_sp$propcrop)
  prox2water_f <- ecdf(evt_sp$prox2water)
  
  
  evt_sp2 <- evt_sp %>% 
    mutate(evi = evi_f(evi), # scale w/i the individual
           lst = lst_f(lst),
           prop_crop = prop_crop_f(propcrop),
           prox2water = prox2water_f(prox2water)) %>% 
    pivot_longer(cols = c(evi, lst, propcrop, prox2water), names_to = "variable") %>% 
    filter(complete.cases(.)) 
  
  df_bounds <- evt_sp2 %>%
    summarize(minv = min(value),
              maxv = max(value))
  
  # phen_sp <- phen_summ %>%
  #   # phen_sp <- ind_phen %>%
  #   filter(Species == sp_name) %>%
  #   mutate(Season = case_when(Status == "End Fall" ~ "Fall",
  #                             Status == "Start Fall" ~ "Fall",
  #                             Status == "End Spring" ~ "Spring",
  #                             Status == "Start Spring" ~ "Spring"),
  #          Type = case_when(Status == "End Fall" ~ "Stop",
  #                           Status == "Start Fall" ~ "Start",
  #                           Status == "End Spring" ~ "Stop",
  #                           Status == "Start Spring" ~ "Start")) %>% 
  #   pivot_wider(id_cols = c(Species, Season), names_from = Type, values_from = med) %>% 
  #   bind_cols(df_bounds)
  
  (out_plots[[i]] <- ggplot() +
      
      geom_smooth(data = evt_sp2, aes(x = wk, y = value, group = variable, color = variable, fill = variable),
                  inherit.aes = F, method = "gam")+
      geom_point(data = evt_sp2, aes(x = wk, y = value, color = variable)) +
      scale_y_continuous(limits = c(0,1), oob=rescale_none)+
      scale_color_manual(values = pal, labels = c("EVI", "Temperature", 
                                                  "Crop Proportion",  "Water Proximity"),
                         name = "Niche Component")+
      scale_fill_manual(values = pal, labels = c("EVI", "Temperature", 
                                                 "Crop Proportion",  "Water Proximity"),
                        name = "Niche Component")+
      ggtitle(sp_name)+
      xlab("Week") +
      ylab("Scaled Value")+
      ylim(c(-0.1,1.15))+
      theme_classic() +
      theme(
        # legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
  )
}

plot(out_plots[[1]])
plot(out_plots[[2]])
plot(out_plots[[3]])
plot(out_plots[[4]])


(vars_plot_fig <- out_plots[[1]] + out_plots[[2]] + out_plots[[3]] + out_plots[[4]] + 
    plot_layout(ncol = 2, nrow = 2, guides = "collect", widths = c(4,4), heights = c(3,3))) &
  theme(legend.position = "bottom")

