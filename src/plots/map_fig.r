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
  collect()

# Species level phenology
phen_summ <- read_csv("out/phenology_summary.csv")

# Individual centroids data
ind_cents <- read_csv("out/individual_seasonal_centroids.csv")


#---- All indInds Map ----#

#get world map of country boundaries 

world <- ne_countries(scale = "small", returnclass = "sf")

# hypso <- raster("data/HYP_50M_SR_W/HYP_50M_SR_W.tif")
# hypso_mask <- mask(hypso, world)
# plot(hypso_mask)
# plot(world)
# hypso_spdf <- as(hypso, "SpatialPixelsDataFrame")
# convert hypsos to data.frames (for ggplot2)
# hypso_df <- as.data.frame(hypso, xy=T)

evt_sf <- evt0 %>% 
  filter(species != "Balearica pavonina") %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) 

evt_thin <- evt_sf %>% 
  group_by(individual_id, doy) %>% 
  # summarize(geometry = geometry[1])
  filter(row_number(geometry) == 1) %>% 
  mutate(species_common = case_when(species == "Anthropoides virgo" ~ "Demoiselle Crane",
                                                  species == "Grus grus" ~ "Common Crane",
                                                  species == "Grus nigricollis" ~ "Black-necked Crane",
                                                  species == "Grus vipio" ~ "White-naped Crane"))

evt_line <- evt_thin %>% 
  group_by(individual_id) %>%
  summarise(species = species_common[1],
            do_union = F) %>%
  st_cast("LINESTRING")

# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pal <- c("#E69F00", "#009E73", "#8B0000", "#CC79A7")

(all_map <- ggplot()+
    geom_sf(data = world, fill = "#ECECEC")+
    geom_sf(data = evt_line, aes(color = species), alpha = 0.6)+
    scale_color_manual(values = pal, name = "Species") +
    coord_sf(xlim = c(-15, 150), ylim = c(5, 65)) +
    theme_minimal() +
    theme(plot.margin = margin(-10,0,-10,0),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "bottom")+
    guides(color = guide_legend(override.aes = list(linewidth = 2))))

leg <- cowplot::get_legend(all_map)
plot(leg)
ggsave(leg, file = "out/map_lg.png",  dpi = 600)
ggsave(leg, file = "out/map_lg.pdf")

(all_map <- all_map +
    theme(legend.position = "none"))


ggsave(all_map, file = "out/all_map.png", width = 5.5, height = 2.5, dpi = 600)
ggsave(all_map, file = "out/all_map.pdf", width = 5.5, height = 2.5)

#---- Individual Inset ----#

evt_line %>%
  filter(species == "Grus nigricollis") %>%
  # filter(species == "Grus grus") %>%
  pull(individual_id) %>%
  unique()
# unique(evt_ind1_pts$yr)
# 10722504 
# 28714952 2015
# focal_ind <-   146932155
focal_ind <-   55754621

focal_yr <- 2015

ind_phen <- read_csv(glue("data/shiny_track_segmentation/Grus_nigricollis_55754621_2022-11-17_20_24_29_Lia.csv")) %>%
  mutate(Status = case_when(Status == "Presumed Start Fall" ~"Start Fall",
                            Status == "Presumed Start Spring" ~"Start Spring",
                            Status == "Presumed End Fall" ~"End Fall",
                            Status == "Presumed End Spring" ~"End Spring",
                            T ~ Status
  ),
  yr = year(Date),
  wk = week(Date)) %>% 
  filter(yr == focal_yr)

evt_ind1_pts <- evt_sf %>% 
  filter(individual_id == focal_ind) %>%
  filter(yr == focal_yr)

evt_ind1_ln <- evt_ind1_pts %>% 
  group_by(individual_id) %>%
  summarise(do_union = F) %>%
  st_cast("LINESTRING")


# evt_thin <- evt_ind1_pts %>% 
#   group_by(doy) %>% 
#   summarise(geometry = geometry[1],
#             individual_id = individual_id[1])
# evt_thin_ln <- evt_thin %>% 
#   group_by(individual_id) %>%
#   summarise(do_union = F) %>%
#   st_cast("LINESTRING")
# 
# mapview::mapview(evt_thin, zcol = "doy")+
# mapview::mapview(evt_thin_ln)

(indmap1 <- ggplot()+
    geom_sf(data = world, fill = "#ECECEC")+
    # geom_sf(data = evt_ind1_pts,  aes(color = doy), alpha = 0.4, size = 0.5)+
    # scale_color_viridis_c()+
    geom_sf(data = evt_ind1_pts,  color = "#E69F00", alpha = 0.4, size = 0.5)+
    geom_sf(data = evt_ind1_ln, color = "#E69F00", alpha = 1, linewidth = 0.8)+
    # scale_color_manual(values = "#009E73", name = "Species") +
    scale_x_continuous(breaks = seq(89, 92, by = 1)) +
    # coord_sf(xlim = c(-10, 30), ylim = c(35, 65), expand = F) +
    coord_sf(xlim = c(89, 92), ylim = c(27, 29.5), expand = F) +
    # theme_bw() +
    annotation_scale(location = "bl", width_hint = 0.4) +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering,
                           height = unit(1, "cm"),
                           width = unit(1, "cm"),
                           text_size = 5)+
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    NULL)

ggsave(indmap1, file = "out/ind_map1.png", width = 2, height = 2, units = "in")
ggsave(indmap1, file = "out/ind_map1.pdf", width = 2, height = 2, units = "in")


#-- EVI Plot --3

df_ind1 <- evt0 %>% 
  filter(individual_id == focal_ind) %>% 
  filter(yr == focal_yr) %>% 
  mutate(wk = week(timestamp),
         prox2water = 1/dist2water) %>% 
  group_by(wk) %>% 
  summarize(evi = mean(evi, na.rm = T),
            lst = mean(lst, na.rm = T),
            prop_crop = mean(prop_crop, na.rm = T),
            dist2water = mean(dist2water, na.rm = T),
            prox2water = case_when(dist2water == 0 ~ 1, 
                                   TRUE ~ 1/dist2water)) 

# create empirical dist functions
evi_f <- ecdf(df_ind1$evi)
lst_f <- ecdf(df_ind1$lst)
prop_crop_f <- ecdf(df_ind1$prop_crop)
dist2water_f <- ecdf(df_ind1$dist2water)
prox2water_f <- ecdf(df_ind1$prox2water)



df_ind2 <- df_ind1 %>% 
  mutate(evi = evi_f(evi), # scale w/i the individual
         lst = lst_f(lst),
         prop_crop = prop_crop_f(prop_crop),
         dist2water = dist2water_f(dist2water),
         prox2water = prox2water_f(prox2water)) %>% 
  pivot_longer(cols = c(evi, lst, prop_crop, dist2water, prox2water), names_to = "variable") %>% 
  filter(complete.cases(.)) %>% 
  # filter(variable == "evi" | variable == "prop_crop")
  filter(variable == "evi" | variable == "prox2water")

df_bounds <- df_ind2 %>%
  summarize(minv = min(value),
            maxv = max(value))

phen_sp <- phen_summ %>%
  # phen_sp <- ind_phen %>%
  filter(Species == "Grus nigricollis") %>%
  mutate(Season = case_when(Status == "End Fall" ~ "Fall",
                            Status == "Start Fall" ~ "Fall",
                            Status == "End Spring" ~ "Spring",
                            Status == "Start Spring" ~ "Spring"),
         Type = case_when(Status == "End Fall" ~ "Stop",
                          Status == "Start Fall" ~ "Start",
                          Status == "End Spring" ~ "Stop",
                          Status == "Start Spring" ~ "Start")) %>% 
  pivot_wider(id_cols = c(Species, Season), names_from = Type, values_from = med) %>% 
  bind_cols(df_bounds)

(evi_plot <- ggplot() +
    # geom_rect(data = phen_sp, aes(xmin = Start, xmax = Stop, 
    #                               ymin = 1,
    #                               ymax = 0,
    #                               fill = Season), alpha = 0.5, inherit.aes = F) +
    # scale_fill_manual(values = c("#F3872F", "#5B9E48")) +
    # new_scale_fill()+
    geom_smooth(data = df_ind2, aes(x = wk, y = value, group = variable, color = variable, fill = variable),
                inherit.aes = F, method = "gam")+
    geom_point(data = df_ind2, aes(x = wk, y = value, color = variable)) +
    scale_y_continuous(limits = c(0,1), oob=rescale_none)+
    scale_color_manual(values = c("#416000", "#3B75A9"), labels = c("EVI", "Dist. to Water"), 
                       name = "Niche Component")+
    scale_fill_manual(values = c("#416000", "#3B75A9"), labels = c("EVI", "Dist. to Water"), 
                      name = "Niche Component")+
    # geom_line(data = df_ind2, aes(x = wk, y = value, color = variable)) +
    xlab("Week") +
    ylab("Scaled Value")+
    # ylim(c(-0.1,1.15))+
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12))
)
  
  ggsave(evi_plot, file = "out/evi_plot.png", width = 6, height = 2, units = "in")
  ggsave(evi_plot, file = "out/evi_plot.pdf", width = 6, height = 2, units = "in")
  
  # #---- Again fro GRus grus ----#
  # # focal_ind <-   146932155
  # # focal_yr <- 2016
  # 
  # focal_ind <-   28714952
  # focal_yr <- 2015
  # 
  # 
  # ind_phen <- read_csv(glue("data/shiny_track_segmentation/Grus_grus_28714952_2022-11-09_16_45_38_Lia.csv
  #                         ")) %>%
  #   mutate(Status = case_when(Status == "Presumed Start Fall" ~"Start Fall",
  #                             Status == "Presumed Start Spring" ~"Start Spring",
  #                             Status == "Presumed End Fall" ~"End Fall",
  #                             Status == "Presumed End Spring" ~"End Spring",
  #                             T ~ Status
  #   ),
  #   yr = year(Date),
  #   wk = week(Date)) %>% 
  #   filter(yr == focal_yr)
  # 
  # evt_ind1_pts <- evt_sf %>% 
  #   filter(individual_id == focal_ind) %>%
  #   filter(yr == focal_yr)
  # 
  # unique(evt_ind1_pts$yr)
  # evt_ind1_ln <- evt_ind1_pts %>% 
  #   
  #   group_by(individual_id) %>%
  #   summarise(do_union = F) %>%
  #   st_cast("LINESTRING")
  # 
  # 
  # (indmap1 <- ggplot()+
  #     geom_sf(data = world, fill = "#ECECEC")+
  #     geom_sf(data = evt_ind1_pts,  color = "#009E73", alpha = 0.4, size = 0.5)+
  #     geom_sf(data = evt_ind1_ln, color = "#009E73", alpha = 1, linewidth = 0.8)+
  #     # scale_color_manual(values = "#009E73", name = "Species") +
  #     coord_sf(xlim = c(-10, 30), ylim = c(36, 57), expand = F) +
  #     # coord_sf(xlim = c(89, 92), ylim = c(27.4, 30), expand = F) +
  #     theme_minimal() +
  #     theme(legend.position = "none"))
  # 
  # ggsave(indmap1, file = "out/ind_map2.png")
  # 
  # 
  # #-- EVI Plot --3
  # 
  # df_ind1 <- evt0 %>% 
  #   filter(individual_id == focal_ind) %>% 
  #   filter(yr == focal_yr) %>% 
  #   mutate(wk = week(timestamp)) %>% 
  #   group_by(wk) %>% 
  #   summarize(evi = mean(evi, na.rm = T),
  #             lst = mean(lst, na.rm = T),
  #             prop_crop = mean(prop_crop, na.rm = T),
  #             dist2water = mean(dist2water, na.rm = T)) 
  # 
  # # create empirical dist functions
  # evi_f <- ecdf(df_ind1$evi)
  # lst_f <- ecdf(df_ind1$lst)
  # prop_crop_f <- ecdf(df_ind1$prop_crop)
  # dist2water_f <- ecdf(df_ind1$dist2water)
  # 
  # 
  # 
  # df_ind2 <- df_ind1 %>% 
  #   mutate(evi = evi_f(evi), # scale w/i the individual
  #          lst = lst_f(lst),
  #          prop_crop = prop_crop_f(prop_crop),
  #          dist2water = rev(dist2water_f(dist2water))) %>% 
  #   pivot_longer(cols = c(evi, lst, prop_crop, dist2water), names_to = "variable") %>% 
  #   filter(complete.cases(.)) %>% 
  #   # filter(variable == "evi" | variable == "prop_crop")
  #   filter(variable == "dist2water" | variable == "lst")
  # 
  # df_bounds <- df_ind2 %>%
  #   summarize(minv = min(value),
  #             maxv = max(value))
  # 
  # phen_sp <- phen_summ %>%
  #   # phen_sp <- ind_phen %>%
  #   filter(Species == "Grus nigricollis") %>%
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
  # 
  # (evi_plot <- ggplot() +
  #     # geom_rect(data = phen_sp, aes(xmin = Start, xmax = Stop, 
  #     #                               ymin = 1,
  #     #                               ymax = 0,
  #     #                               fill = Season), alpha = 0.5, inherit.aes = F) +
  #     # scale_fill_manual(values = c("#F3872F", "#5B9E48")) +
  #     # new_scale_fill()+
  #     geom_smooth(data = df_ind2, aes(x = wk, y = value, group = variable, color = variable, fill = variable),
  #                 inherit.aes = F, method = "gam")+
  #     geom_point(data = df_ind2, aes(x = wk, y = value, color = variable)) +
  #     scale_color_manual(values = c("#990000", "#A38A00"), labels = c("Temperature", "Prop. Crops"), 
  #                        name = "Niche Component")+
  #     scale_fill_manual(values = c("#990000", "#A38A00"), labels = c("Temperature", "Prop. Crops"), 
  #                       name = "Niche Component")+
  #     # geom_line(data = df_ind2, aes(x = wk, y = value, color = variable)) +
  #     xlab("Week") +
  #     ylab("Scaled Value") +
  #     ylim(c(-0.1,1.25))+
  #     theme_classic() +
  #     # theme(legend.position = "none") +
  #     NULL)
  # 
  # ggsave(evi_plot, file = "out/lst_plot.png", width = 3, height = 1.5, units = "in")
  
  #---- Combine Plots ----#
  
  # Currently done by hand in slides
  
  