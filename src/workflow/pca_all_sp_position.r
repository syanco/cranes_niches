#########################################
####----    PCA Among Species    ----####
#########################################

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

##################################
#---- Initialize Environment ----#
##################################

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
  }))

# init vec of vars & labels
vars <- c("prox2water", "evi", "propcrop", "lst")
labels <- c("Prox. to Water", "EVI", "Prop. Crops", "Temp.")

# Set color palettes
# pal2 <- diverging_hcl(2, palette = "Blue-Red2", rev = T)
pal <- c("#E69F00", "#009E73", "#8B0000", "#CC79A7")
pal2 <- c("N" = "#D33F6A", "P" = "#4A6FE3")
#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)


###########################
#---- Initialize Data ----#
###########################

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

# species_list <- samples %>% pull(species)


##########################
#---- Run & Plot PCA ----#
##########################

# Prep Data
evt1 <- evt0 %>% 
  filter(species != "Balearica pavonina") %>% # Remove B pavonina
  mutate(grp = glue("{individual_id}_{yr}"),
         ind_f = as.factor(individual_id),
         ind_num = as.numeric(ind_f),
         dt = as.numeric(difftime(timestamp, dplyr::lag(timestamp, 1)), units='secs'), # time diff (secs)
         v = sl/dt, # velocity (m/s)
         v = case_when(is.infinite(v) ~ NA,
                       !is.infinite(v) ~ v),
         v_scale = scale(v),
         rad = bearing*(pi/180),
         rad_scale = scale(rad),
         wk = week(timestamp),
         prox2water = case_when(dist2water == 0 ~ 1,
                                TRUE ~ 1/dist2water),
         prox2water_scale = scale(prox2water)) %>% 
  group_by(species, grp, wk) %>% 
  summarize(dist2water_scale = mean(dist2water_scale, na.rm = T), 
            prox2water_scale = mean(prox2water_scale, na.rm = T),
            evi_scale = mean(evi_scale, na.rm = T), 
            propcrop_scale = mean(prop_crop_scale, na.rm = T), 
            lst_scale = mean(lst_scale, na.rm = T),
            individual_id = individual_id[1],
            yr = yr[1]) %>% 
  group_by(species, wk) %>% # average individuals (for now...)
  summarize(dist2water_scale = median(dist2water_scale, na.rm = T), 
            prox2water_scale = median(prox2water_scale, na.rm = T),
            evi_scale = median(evi_scale, na.rm = T), 
            propcrop_scale = median(propcrop_scale, na.rm = T),
            lst_scale = median(lst_scale, na.rm = T)) %>% 
  filter_at(vars(
    prox2water_scale,
    evi_scale, 
    propcrop_scale,
    lst_scale), all_vars(!is.na(.))) %>% 
  rename("Prox. to Water" = "prox2water_scale",
         "EVI" = "evi_scale",
         "Prop. Crops" = "propcrop_scale",
         "Temp." = "lst_scale") %>% 
  ungroup() %>% 
  mutate("ID" = glue("{species}-{wk}"))


evt_pca <- evt1 %>%   
  column_to_rownames("ID") %>% 
  select(c(`Prox. to Water`, EVI, `Prop. Crops`, Temp.))


# Run PCA

pca_fit <- prcomp(evt_pca, retx=TRUE, scale=F)


plot_dat <- evt_pca %>%
  rownames_to_column("ID") %>%
  left_join(evt1)

coord_plot <- autoplot(pca_fit, data = plot_dat,
                       variance_percentage = F,
                       frame = T,
                       colour = 'species',
                       loadings.size = 3,
                       loadings = F, 
                       loadings.label = F,
                       loadings.label.colour = "black",
                       loadings.label.size = 4,
                       loadings.label.repel = T,
                       loadings.colour = "black",
                       loadings.label.vjust = 1.8) +
  scale_colour_manual(values = pal, name = "Species")+
  scale_fill_manual(values=pal, name = "Species")+
  scale_x_continuous(breaks = seq(-0.15, 0.1, by = 0.1))+
  scale_y_continuous(breaks = seq(-0.1, 0.3, by = 0.1))+
  # geom_path(aes(color = wk))+
  # scale_colour_viridis_c(limits = c(0, 53), name = "Week")+
  # ggtitle(sp_name) +
  theme_linedraw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title =  element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
coord_plot$layers[[2]]$aes_params$size <- 1

coord_plot
ggsave(plot = coord_plot, file = "out/all_sp_pca.png", width = 2.5, height = 2.5)
#this combines with the output of the map_fg.r script too....  sloppy i know...

comb <- (all_map + coord_plot+plot_layout(guides = "auto", widths = c(10,4))) / guide_area() + plot_layout(guides = "collect", heights = c(5,1)) 

ggsave(plot = comb, filename = "out/draft_figs/map_pca.png")


