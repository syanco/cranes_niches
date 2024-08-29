



#---- Input Parameters ----#
if(interactive()) {
  # library(here)
  
  .wd <- '.'
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
pal2 <- diverging_hcl(2, palette = "Blue-Red2", rev = T)

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

species_list <- samples %>% pull(species)


##########################
#---- Run & Plot PCA ----#
##########################

# Init lists to hold outputs
sp_plot <- list()
coord_plots <- list()
scree_plots <- list()
time_plots <- list()
scores_plots_sp <- list()
loadings_plots_sp <- list()

#-- Loop over species --#
for(j in 1:length(species_list)){
  
  # Focal sp name
  sp_name <- species_list[j]
  
  # Extract species data
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
  
  # Make weekly vars
  evt_wk <- evt_sp %>% 
    mutate(wk = week(timestamp)) %>% 
    group_by(grp, wk) %>% 
    summarize(dist2water_scale = var(dist2water_scale, na.rm = T), 
              evi_scale = var(evi_scale, na.rm = T), 
              propcrop_scale = var(prop_crop_scale, na.rm = T), 
              lst_scale = var(lst_scale, na.rm = T),
              # v_scale = mean(v_scale, na.rm = T),
              # rad_scale = mean(rad_scale, na.rm = T),
              individual_id = individual_id[1],
              yr = yr[1]) %>% 
    group_by(wk) %>% # average individuals (for now...)
    summarize(dist2water_scale = mean(dist2water_scale, na.rm = T), 
              evi_scale = mean(evi_scale, na.rm = T), 
              propcrop_scale = mean(propcrop_scale, na.rm = T),
              lst_scale = mean(lst_scale, na.rm = T)
              # v_scale = mean(v_scale, na.rm = T),
              # rad_scale = mean(rad_scale, na.rm = T)) 
    )%>% 
    column_to_rownames("wk") %>%
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
  
  # If no complete cases, then break
  if(nrow(evt_wk) == 0){
    message(glue("Not enough data for {sp_name}, so sad, too bad..."))
    next}
  # # show data coverage
  # 
  # time_cov <- evt_wk %>% 
  #   # mutate(ind_f = as.factor(individual_id)) %>% 
  #   group_by(grp) %>% 
  #   summarize(start = min(wk),
  #             stop = max(wk),
  #             dur = stop-start,
  #             individual_id = individual_id[1])
  # ggplot(time_cov) +
  #   geom_segment(aes(y = grp, yend = grp, 
  #                    x = start, xend = stop))
  # 
  # # find duration cutoff
  # # TODO: need to find a better way to do this - should I average over individuals to make a population scale metric?
  # #     or impute within individuals?  Data missingness is a an issue here...
  # nrow(time_cov %>% filter(dur > 45) )
  # inds <- time_cov %>% 
  #   filter(dur >= 40) %>% 
  #   pull(grp)
  # 
  # # Format for PCA
  # evt_wide <- evt_wk %>%
  #     # filter(grp %in% inds) %>% 
  #     pivot_wider(id_cols = c("wk"), names_from = c("grp", ),
  #               values_from = c("dist2water_scale", "evi_scale", "propcrop_scale", "ghm_scale")) %>% 
  #   arrange(wk) %>% 
  #   column_to_rownames("wk")
  
  # #   ### Find number of NAs
  # na_col_count <- apply(is.na(evt_wide),2,sum)
  # na_col_test <- na_col_count > 0
  # data_mat <- evt_wide[,!na_col_test]
  # # data_mat <- evt_wide
  # 
  # data_mat <- as.matrix(data_mat)
  
  
  #--  Perform PCA --#
  pca_fit <- prcomp(evt_wk, retx=TRUE, scale=F)
  
  
  #-- Create Scree plot --#
  #TODO: remove this section of code?
  
  # get PCA Variance stats
  pc_importance <- summary (pca_fit)
  
  # # Format data
  # Scree <- data.frame(cbind(Component=seq(1,dim(pc_importance$importance)[2]),t(pc_importance$importance)))
  # Scree$EigenVal <- Scree$Standard.deviation^2
  # Scree_portion <- Scree[1:3,]
  # 
  # # Plot (code based on Jarzyna et al)
  # sp <- ggplot(Scree_portion, aes(x = Component, y  = Cumulative.Proportion*100)) %>%
  #   + geom_line(size=2) %>%
  #   + geom_point(alpha=1,size=4, pch=16) %>%
  #   + scale_x_continuous(name = "Principal Component", limits=c(1,12),breaks=c(seq(1,12,1))) %>%
  #   + scale_y_continuous(name = "Cummulative Variance (%)", limits=c(0,100),breaks=c(0,20,40,60,80,100)) %>%
  #   #TODO: change these theme settings
  #   + theme(panel.background = element_rect(fill = 'white', colour = 'grey85'),
  #           panel.border = element_rect(fill=NA, colour = "white", size=1),
  #           axis.line = element_line(color = 'black', size=1.5),
  #           plot.title = element_text(size=15, vjust=2, family="sans"),
  #           axis.text.x = element_text(colour='black',size=22),
  #           axis.text.y = element_text(colour='black',size=22),
  #           axis.title.x = element_text(colour='black',size=27),
  #           axis.title.y = element_text(colour='black',size=27),
  #           axis.ticks = element_line(color = 'black', size=1.5),
  #           axis.ticks.length=unit(0.3,"cm"),
  #           legend.position="none",
  #           legend.text=element_text(size=20),
  #           legend.title=element_blank(),
  #           panel.grid.major = element_blank(), 
  #           panel.grid.minor = element_blank())
  # scree_plots[[j]] <- sp
  
  
  #-- PC Scores Over Time --#
  
  # Init list
  scores_plots <- list()
  
  # No. of axes to consider
  axes <- 2
  
  # Loop over PC Axes
  for(i in 1:axes){
    
    # Extract scores
    plot_df <- data.frame(week = seq(1,nrow(evt_wk)), score = pca_fit$x[,i]) # this last subset grabs the PC
    #plot_df$score <- plot_df$score*(-1)   ##do this only for PC2
    
    
    #- Get info from migration boxes
    
    # get min and max var values
    score_max <- plot_df %>% 
      summarize(score_max = max(score, na.rm = T),
                score_min = min(score, na.rm = T)) %>% 
      mutate(Species = sp_name)
    
    # mod phenology df and join with min/max
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
      left_join(score_max)
    
    ### Create plot
    scores_plots[[i]] <- ggplot(plot_df, aes(x=week, y=score)) +
      geom_line(linewidth=1) +
      geom_point(alpha=1,size=2, pch=16) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_rect(data = phen_sp, aes(xmin = Start, xmax = Stop, 
                                    ymin = !!sym("score_min")-0.5,
                                    ymax = !!sym("score_max")+0.5,
                                    fill = Season), alpha = 0.5,
                inherit.aes = F) +
      scale_fill_manual(values = c("orange", "lightblue")) +      
      # geom_rect(aes(xmin = spring_start, xmax = spring_stop, ymin = min(score-0.1), ymax = max(score+0.1)), alpha = 0.01, fill = "lightblue") +
      # geom_rect(aes(xmin = fall_start, xmax = fall_stop, ymin = min(score-0.1), ymax = max(score+0.1)), alpha = 0.01, fill = "orange") +
      scale_x_continuous(name = "Time step (week)") +
      scale_y_continuous(name = "Score") +
      ggtitle(glue("PC{i}")) +
      theme_classic()+
      theme(legend.position = "none")
    
    
  }
  
  # # make coord plot
  plot_dat <- evt_wk %>%
    rownames_to_column("wk") %>%
    mutate(wk = as.numeric(wk))
  
  tmpplt <- autoplot(pca_fit, data = plot_dat,
                     colour = 'wk',
                     loadings.size = 3,
                     loadings = T, 
                     loadings.label = T,
                     loadings.label.colour = "black",
                     loadings.label.repel = T,
                     loadings.colour = "black",
                     loadings.label.vjust = 1.8) +
    geom_path(aes(color = wk))+
    scale_colour_viridis_c(limits = c(0, 53), name = "Week")+
    ggtitle(sp_name) +
    theme_linedraw()
  
  tmpplt$layers[[2]]$aes_params$size <- 1
  
  coord_plots[[j]] <- tmpplt
  
  loadings_plots <- list()
  for(i in 1:axes){
    # laodings across ind-vars
    ### Create temporary loading dataframe
    loading_temp <- data.frame(var = rownames(pca_fit$rotation), 
                               load = pca_fit$rotation[,i]) %>% 
      mutate(dir = case_when(load <= 0 ~ "N",
                             load > 0 ~ "P"))
    
    
    loadings_plots[[i]] <- ggplot(loading_temp) +
      geom_col(aes(x = load, y = var, fill = dir), 
               width = 0.5) + 
      scale_fill_manual(limits = levels(loading_temp$dir),values = pal2) +
      # scale_y_discrete(label = c("Distance to Water", "EVI", "Proportion Crops"))  +
      geom_vline(aes(xintercept = 0)) +
      xlab("") +
      ylab("") +
      theme_classic() +
      theme(legend.position = "none")
    
  } # i
  
  scores_plots_sp[[j]] <- scores_plots
  loadings_plots_sp[[j]] <- loadings_plots
  
} # j


# Save out PCAs
save(pca_out, file = "out/fitted_var_pcas.rdata")


###--- Plots ---###
(coord_plot_fig <- coord_plots[[1]] + coord_plots[[3]] + coord_plots[[4]] + coord_plots[[5]] + 
    plot_layout(ncol = 2, nrow = 2, guides = "collect") & 
    theme(legend.position = 'bottom'))

ggsave(coord_plot_fig, filename = "out/coord_plot_fig.png", width = 7, height = 7, units = "in")


# Time plot fig

# Unpack plots
scores_plots_sp  <- unlist(scores_plots_sp, recursive = F)

loadings_plots_sp <- unlist(loadings_plots_sp, recursive = F)

a <- (wrap_plots(scores_plots_sp[1:2], ncol = 1) | wrap_plots(loadings_plots_sp[1:2], ncol = 1)) +
  plot_layout(widths = c(4, 1)) + plot_annotation(title = species_list[1])
b <- (wrap_plots(scores_plots_sp[3:4], ncol = 1) | wrap_plots(loadings_plots_sp[3:4], ncol = 1)) +
  plot_annotation(title = species_list[3]) + plot_layout(widths = c(4, 1))
c <- (wrap_plots(scores_plots_sp[5:6], ncol = 1) | wrap_plots(loadings_plots_sp[5:6], ncol = 1)) +
  plot_annotation(title = species_list[4]) + plot_layout(widths = c(4, 1))
d <- (wrap_plots(scores_plots_sp[7:8], ncol = 1) | wrap_plots(loadings_plots_sp[7:8], ncol = 1)) +
  plot_annotation(title = species_list[5]) + plot_layout(widths = c(4, 1))

(scores_fig <- (wrap_elements(a)/wrap_elements(b))|(wrap_elements(c)/wrap_elements(d)))

ggsave(scores_fig, file = "out/scores-fig.png", width = 15, height = 8, units = "in")