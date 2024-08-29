# scratch plotting code

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggthemes)
    library(patchwork)
    library(dplyr)
    library(glue)
    library(tidyr)
  }))

#-- Interannual Plot --#

dat_yrs0 <- read.csv("out/wi_ind_ax_yrs.csv")

# simple outlier filtering for now...
dat_yrs <- dat_yrs0 %>% 
  filter(
    # Bhattacharyya_distance.total < 100,
         species != "Anthropoides paradiseus")

dat_yrs %>% 
  group_by(species) %>% 
  summarise(n = n())

(pl_years <- ggplot(dat_yrs)+
    geom_boxplot(aes(y=Bhattacharyya_distance.total, x = season, color = niche_type), 
                 outlier.shape = NA)+
    facet_wrap(~species, nrow = 1) +
    ylab("Niche Dissimilarity")+
    ylim(0,100)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

# ggsave("out/draft_figs/wi_ind_across_yrs.png", pl_years)

#-- Interseasonal Plot --#

dat_seas0 <- read.csv("out/wi_ind_ax_seas.csv")

# simple outlier filtering for now...
dat_seas <- dat_seas0 %>% 
  filter(
    # Bhattacharyya_distance.total < 100,
    species != "Anthropoides paradiseus") %>% 
  unite(comp, season1:season2, remove = F)

dat_seas %>% 
  group_by(species) %>% 
  summarise(n = n()) 

(pl_seasons <- dat_seas %>% 
    filter(comp == "Spring_Fall" | comp == "Winter_Summer") %>% 
    ggplot()+
    geom_boxplot(aes(y=Bhattacharyya_distance.total, x = comp, color = niche_type), 
                 outlier.shape = NA)+
    facet_wrap(~species, nrow = 1) +
    ylab("Niche Dissimilarity")+
    xlab("")+
    ylim(0,100)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))


#-- Combine Plots --#

(comb <- pl_years + pl_seasons + plot_layout(guides = 'collect'))

ggsave(comb, filename = "out/draft_figs/comb_plot.png")

#-- Data Summary --#

dat_seas %>% 
  filter(comp == "Spring_Fall" | comp == "Winter_Summer") %>%
  group_by(species, comp) %>% 
  summarize(n = n(),
            inds = n_distinct(individual_id),
            yrs = n_distinct(year))
