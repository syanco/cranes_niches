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
    # library(MVNH)
    # library(INLA)
    library(ggfortify)
    library(zoo)
    library(patchwork)
    library(colorspace)
    library(tidyverse)
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
  collect()

phen_summ <- read_csv("out/phenology_summary.csv")

# sample sizes
(samples <- evt0 %>% 
    group_by(species) %>% 
    summarise(n = n(),
              ind = n_distinct(individual_id)))

species_list <- samples %>% 
  filter(species != "Balearica pavonina") %>%
  pull(species)

#---- Perform Analysis ----#

# init list to store output
sp_out <- list()

# init vec of vars
vars <- c("dist2water", "evi", "propcrop", "lst")
labels <- c("Distance to Water", "EVI", "Proportion Crops", "Temp.")

# Extract a single species
# j <- 2
for(j in 1:length(species_list)){
  
  sp_name <- species_list[j]
  
  
  # extract species
  evt_sp <- evt0 %>% 
    mutate(grp = glue("{individual_id}_{yr}"),
           ind_f = as.factor(individual_id),
           ind_num = as.numeric(ind_f)) %>% 
    filter(species == sp_name) %>% 
    mutate(dt = as.numeric(difftime(timestamp, dplyr::lag(timestamp, 1)), units='secs'), # time diff (secs)
           v = sl/dt, # velocity (m/s)
           v = case_when(is.infinite(v) ~ NA,
                         !is.infinite(v) ~ v)
    ) %>% 
    ungroup() %>% 
    mutate(v_scale = scale(v),
           rad = bearing*(pi/180),
           rad_scale = scale(rad))
  
  # Make UNSCALED weekly vars
  evt_indwk <- evt_sp %>% 
    mutate(wk = week(timestamp)) %>% 
    group_by(grp, wk) %>% 
    # get weekly individual means
    summarize(dist2water = mean(dist2water, na.rm = T), 
              evi = mean(evi, na.rm = T), 
              propcrop = mean(prop_crop, na.rm = T), 
              lst = mean(lst_scale, na.rm = T),
              v = mean(v, na.rm = T),
              rad = mean(rad, na.rm = T),
              individual_id = individual_id[1],
              yr = yr[1]) 
  
  evt_wk <- evt_indwk %>% 
    group_by(wk) %>% # average individuals (for now...)
    # get group-wide means
    summarize(dist2water = mean(dist2water, na.rm = T), 
              evi = mean(evi, na.rm = T), 
              propcrop = mean(propcrop, na.rm = T), 
              lst = mean(lst, na.rm = T),
              v = mean(v, na.rm = T),
              rad = mean(rad, na.rm = T)) 
  
  evt_max <- evt_indwk %>% 
    ungroup() %>%
    summarize(dist2water_max = max(dist2water, na.rm = T),
              evi_max = max(evi, na.rm = T),
              propcrop_max = max(propcrop, na.rm = T),
              lst_max = max(lst, na.rm = T),
              dist2water_min = min(dist2water, na.rm = T),
              evi_min = min(evi, na.rm = T),
              propcrop_min = min(propcrop, na.rm = T),
              lst_min = min(lst, na.rm = T)) %>% 
    mutate(Species = sp_name)
  
  phen_sp <- phen_summ %>% 
    filter(Species == sp_name) %>%
    mutate(Season = case_when(Status == "End Fall" ~ "Fall",
                              Status == "Start Fall" ~ "Fall",
                              Status == "End Spring" ~ "Spring",
                              Status == "Start Spring" ~ "Spring"),
           Type = case_when(Status == "End Fall" ~ "Stop",
                            Status == "Start Fall" ~ "Start",
                            Status == "End Spring" ~ "Stop",
                            Status == "Start Spring" ~ "Start")) %>% 
    pivot_wider(id_cols = c(Species, Season), names_from = Type, values_from = med) %>% 
    left_join(evt_max)
  
  
  
  # Color Pallette
  n_ind <- evt_sp %>% summarize(n=n_distinct(grp)) %>% pull(n)
  pal <- sequential_hcl(n_ind, palette = "ag_Sunset")
  
  
  # init list to store var plots
  out_vars <- list()
  for(i in 1:length(vars)){
    
    p <- ggplot() +
      geom_rect(data = phen_sp, aes(xmin = Start, xmax = Stop, 
                                    # ymin = min(!!sym(vars[i])), 
                                    # ymax = max(!!sym(vars[i]))
                                    ymin = !!sym(glue("{vars[i]}_min")),
                                    ymax = !!sym(glue("{vars[i]}_max")),
                                    fill = Season), alpha = 0.5) +
      scale_fill_manual(values = c("orange", "lightblue")) +
      geom_line(data = evt_indwk, aes(x = wk, y = !!sym(vars[i]), color = grp), 
                # color = "gray", 
                alpha = 0.5) +
      scale_color_manual(values = pal)+
      geom_line(data = evt_wk, aes(x = wk, y = !!sym(vars[i])), size = 3) +
      ylab(labels[i]) +
      xlab("Week") +
      theme_minimal() +
      theme(legend.position = "none")
    
    out_vars[[i]] <- p 
  } # i
  
  sp_out[[j]] <- out_vars
  # [[1]] / out_vars[[2]] / out_vars[[3]]
  
  # names(sp_out[[j]]) <- sp_name
  
} # j


ß#-- Assemble plots --#

layout <- "AAABBB
           AAABBB
           CCCDDD
           CCCDDD"

all_plots <- unlist(sp_out, recursive = F)

rawplot <- wrap_elements(all_plots[[1]] / all_plots[[2]] / all_plots[[3]] / all_plots[[4]] + plot_annotation(title = species_list[1])) +
  wrap_elements(all_plots[[5]] / all_plots[[6]] / all_plots[[7]] / all_plots[[8]] + plot_annotation(title = species_list[2])) +
  wrap_elements(all_plots[[9]] / all_plots[[10]] / all_plots[[11]] / all_plots[[12]] + plot_annotation(title = species_list[3]))+
  wrap_elements(all_plots[[13]] / all_plots[[14]] / all_plots[[15]] / all_plots[[16]] + plot_annotation(title = species_list[4]))+
  plot_layout(design = layout)

rawplot

# Save out
ggsave(plot = rawplot, file = "out/rawplots.png", width = 14, height = 12, units = "in")
