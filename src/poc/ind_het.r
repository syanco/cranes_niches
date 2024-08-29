



#---- Input Parameters ----#
if(interactive()) {
  # library(here)
  
  .wd <- '~/projects/dynamic_niches'
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
    library(MVNH)
    library(INLA)
    library(ggfortify)
    library(zoo)
    library(patchwork)
    library(colorspace)
  }))

# init vec of vars & labels
vars <- c("dist2water", "evi", "propcrop", "lst")
labels <- c("Dist. to H2O", "EVI", "Prop. Crops", "Temp.")

# Set color palettes
# pal2 <- diverging_hcl(2, palette = "Blue-Red2", rev = T)
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

#load PCAs
load("out/fitted_pcas.rdata")

species_list <- samples %>% 
  # filter(species != "Balearica pavonina") %>% 
  pull(species)


out_ls <- list()
for(j in 1:length(species_list)){
  
  sp <- species_list[j]
  
  # Get PCA
  pca <- pca_out[[j]]
  
  #check the species match
  if(is.null(pca$species)){
    message("Species missing from PCA data, probably wasn't fitted...")
    next
  } else{
    if(pca$species != sp){
      message(glue("The species did not match up!!!  The results are wrong!!!"))
      next
    }
  }
  # species level data
  # Extract species data
  evt_sp <- evt0 %>% 
    mutate(grp = glue("{individual_id}_{yr}"),
           ind_f = as.factor(individual_id),
           ind_num = as.numeric(ind_f)) %>% 
    filter(species == sp) %>% 
    mutate(dt = as.numeric(difftime(timestamp, dplyr::lag(timestamp, 1)), units='secs'), # time diff (secs)
           v = sl/dt, # velocity (m/s)
           v = case_when(is.infinite(v) ~ NA,
                         !is.infinite(v) ~ v)
    ) %>% 
    ungroup() %>% 
    mutate(v_scale = scale(v),
           rad = bearing*(pi/180),
           rad_scale = scale(rad))
  
  #TODO - this part doesn't work yet
  # Make weekly vars
  evt_wk <- evt_sp %>% 
    mutate(wk = week(timestamp)) %>% 
    group_by(grp, wk) %>% 
    summarize(dist2water_scale = mean(dist2water_scale, na.rm = T), 
              evi_scale = mean(evi_scale, na.rm = T), 
              propcrop_scale = mean(prop_crop_scale, na.rm = T), 
              lst_scale = mean(lst_scale, na.rm = T),
              # v_scale = mean(v_scale, na.rm = T),
              # rad_scale = mean(rad_scale, na.rm = T),
              individual_id = individual_id[1],
              yr = yr[1]) %>% 
    # group_by(wk) %>% # average individuals (for now...)
    # column_to_rownames("wk") %>%
    # mutate(wk = as.factor(wk)) %>% 
    filter_at(vars(
      dist2water_scale,
      evi_scale, 
      propcrop_scale,
      lst_scale
      # v_scale, 
      # rad_scale
    ),
    all_vars(!is.na(.))) %>% 
    rename("Dist. to H2O" = "dist2water_scale",
           "EVI" = "evi_scale",
           "Prop. Crops" = "propcrop_scale",
           "Temp." = "lst_scale")
  
  #get vec of inds
  inds <- evt_sp %>% 
    pull(individual_id) %>% 
    unique()
  
  # Loop over inds to make predictions
  inds_out <- list()
  for(i in 1:length(inds)){
    
    #ind data
    ind_dat <- evt_wk %>% 
      filter(individual_id == inds[i]) %>% 
      arrange(wk) 
    
    yrs <- ind_dat %>% 
      pull(yr) %>% 
      unique()
    
    yr_out <- list()
    for(y in 1:length(yrs)){
      
      yr_dat <- ind_dat %>% 
        filter(yr == yrs[y])
      
      yr_out[[y]] <- cbind(ind_dat, as.data.frame(predict(pca$pca, newdata = ind_dat)))
    } # y
    
    # gather years output for individual
    inds_out[[i]] <- do.call("rbind", yr_out)
    
  } # 1
  
  out_ls[[j]] <- do.call("rbind", inds_out)
  
} #j

# out_ls <- out_ls[-2] 

#---- Plots ----#


plot_out <- list()
for(j in 1:length(species_list)){
  
  sp_name <- species_list[j]
  
  preds <- out_ls[[j]] 
  
  pred_var <- preds %>% 
    mutate(wk_f = as.factor(wk)) %>%
    group_by(wk_f) %>% 
    summarize(pc_var1 = var(PC1),
              pc_var2 = var(PC2),
              wk = wk[1])
  
  # get min and max var values
  var_max1 <- pred_var %>% 
    summarize(var_max = max(pc_var1, na.rm = T),
              var_min = min(pc_var1, na.rm = T)) %>% 
    mutate(Species = sp_name)
  
  var_max2 <- pred_var %>% 
    summarize(var_max = max(pc_var2, na.rm = T),
              var_min = min(pc_var2, na.rm = T)) %>% 
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
    pivot_wider(id_cols = c(Species, Season), names_from = Type, values_from = med) 
  
  phen_sp1 <- phen_sp %>% 
    left_join(var_max1)
  
  phen_sp2 <- phen_sp %>% 
    left_join(var_max2)
  
  # # Line Plots
  # 
  # (  pc1_ln <- ggplot(preds)+
  #     geom_line(aes(x = wk, y = PC1, group = grp, color = grp))+
  #     scale_color_viridis_d(option = "magma")+
  #     theme_linedraw()+
  #     theme(legend.position = "none")
  #   
  # )    
  # 
  # (  pc2_ln <- ggplot(preds)+
  #     geom_line(aes(x = wk, y = PC2, group = grp, color = grp))+
  #     scale_color_viridis_d(option = "magma")+
  #     theme_linedraw()+
  #     theme(legend.position = "none")
  #   
  # )    
  # 
  # # Box Plots
  # 
  # (  pc1_bx <- preds %>% 
  #     mutate(wk_f = as.factor(wk)) %>% 
  #     ggplot()+
  #     geom_boxplot(aes(x = wk_f, y = PC1))+
  #     scale_color_viridis_d(option = "magma")+
  #     theme_linedraw()+
  #     theme(legend.position = "none")
  #   
  # ) 
  # 
  # (  pc2_bx <- preds %>% 
  #     mutate(wk_f = as.factor(wk)) %>% 
  #     ggplot()+
  #     geom_boxplot(aes(x = wk_f, y = PC2))+
  #     scale_color_viridis_d(option = "magma")+
  #     theme_linedraw()+
  #     theme(legend.position = "none")
  #   
  # ) 
  # 
  
  # straigt variance yo
  pc1_v <- ggplot()+
    geom_point(data = pred_var, aes(x = wk, y = pc_var1),
               alpha = 1, size = 2, pch = 16) +
    geom_line(data = pred_var, aes(x = wk, y = pc_var1))+
    geom_rect(data = phen_sp1, aes(xmin = Start, xmax = Stop, 
                                  ymin = !!sym("var_min")-0.5,
                                  ymax = !!sym("var_max")+0.5,
                                  fill = Season), alpha = 0.5,
              inherit.aes = F) +
    scale_fill_manual(values = c("orange", "lightblue")) +      
    scale_x_continuous(name = "Time step (week)") +
    scale_y_continuous(name = "Variance") +
    ggtitle(glue("PC1")) +
    theme_classic()+
    theme(legend.position = "none")
  
  pc2_v <- ggplot()+
    geom_point(data = pred_var, aes(x = wk, y = pc_var2),
               alpha = 1, size = 2, pch = 16) +
    geom_line(data = pred_var, aes(x = wk, y = pc_var2))+
    geom_rect(data = phen_sp2, aes(xmin = Start, xmax = Stop, 
                                  ymin = !!sym("var_min")-0.5,
                                  ymax = !!sym("var_max")+0.5,
                                  fill = Season), alpha = 0.5,
              inherit.aes = F) +
    scale_fill_manual(values = c("orange", "lightblue")) +      
    scale_x_continuous(name = "Time step (week)") +
    scale_y_continuous(name = "Variance") +
    ggtitle(glue("PC2")) +
    theme_classic()+
    theme(legend.position = "none")

  tmp_out <- list("PC1" = pc1_v,
                  "PC2" = pc2_v)
  plot_out[[j]] <- tmp_out
} # j


var_plots_ul <- unlist(plot_out, recursive = F)

a <- wrap_plots(var_plots_ul[1:2], ncol = 1) + plot_annotation(title = species_list[1])
b <- wrap_plots(var_plots_ul[3:4], ncol = 1) + plot_annotation(title = species_list[2])
c <- wrap_plots(var_plots_ul[5:6], ncol = 1) + plot_annotation(title = species_list[3])
d <- wrap_plots(var_plots_ul[7:8], ncol = 1) + plot_annotation(title = species_list[4])

var_plot <- (wrap_elements(a)/wrap_elements(b))|(wrap_elements(c)/wrap_elements(d))

ggsave(var_plot, file = "out/pc_var_plot.png", width = 13, height = 8)
