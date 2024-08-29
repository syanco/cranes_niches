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

dat0 <- read.csv("out/wi_ind_ax_seas.csv")

# simple outlier filtering for now...
dat1 <- dat0 %>% 
  filter(
    Bhattacharyya_distance.total < 100,
         species != "Anthropoides paradiseus") %>% 
  unite(comp, season1:season2, remove = F)

dat1 %>% 
  group_by(species) %>% 
  summarise(n = n()) 

(pl_seasons <- ggplot(dat1)+
  geom_boxplot(aes(y=Bhattacharyya_distance.total, x = comp, color = niche_type))+
  facet_wrap(~species, nrow = 1) +
  ylab("Niche Dissimilarity")+
    xlab("")+
    ylim(0,100)+
  theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

ggsave("out/draft_figs/wi_ind_across_seas.png", pl_seasons)


# comb plots
#TODO: need to merge the plotting scripts
library(patchwork)

comb <- pl_years + pl_seasons + plot_layout(guides = 'collect')

ggsave(comb, filename = "out/draft_figs/comb_plot.png")
