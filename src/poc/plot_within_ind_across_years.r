# scratch plotting code

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggthemes)
    library(patchwork)
    library(dplyr)
  }))

dat0 <- read.csv("out/wi_ind_ax_yrs.csv")

# simple outlier filtering for now...
dat1 <- dat0 %>% 
  filter(Bhattacharyya_distance.total < 100,
         species != "Anthropoides paradiseus")

dat1 %>% 
  group_by(species) %>% 
  summarise(n = n())

(pl_years <- ggplot(dat1)+
    geom_boxplot(aes(y=Bhattacharyya_distance.total, x = season, color = niche_type),
                 outlier.shape = NA)+
    facet_wrap(~species, nrow = 1) +
    ylab("Niche Dissimilarity")+
    ylim(0, 60)+
    theme_classic())

ggsave("out/draft_figs/wi_ind_across_yrs.png", pl_years)
