library(brms)
library(tidyverse)
library(lme4)

#-- Seasons --#
dat_seas0 <- read.csv("out/wi_ind_ax_seas.csv")

# simple outlier filtering for now...
dat_seas <- dat_seas0 %>% 
  filter(
    # Bhattacharyya_distance.total < 100,
    species != "Anthropoides paradiseus",
    !is.infinite(Bhattacharyya_distance.total)) %>% 
  unite(comp, season1:season2, remove = F) %>% 
  filter(comp == "Spring_Fall" | comp == "Winter_Summer") %>%  
  mutate(individual_id = as.factor(individual_id))

# set model formulat
form <- bf(Bhattacharyya_distance.total ~ 0 + comp*niche_type + (1|individual_id))

#get spp list
sp_ls <- unique(dat_seas$species)

#-- Grus grus --#
mod_out <- list()
for(i in 1:length(sp_ls)){
  dat <- dat_seas %>% 
    filter(species == !!sp_ls[i])


# get_prior(form, data = dat_seas)

fit <- brm(form, data = dat, 
              chains = 3, iter = 10000, cores = getOption("mc.cores", 3),
              control = list(adapt_delta = 0.95))
fit

tmp <- list("species" = sp_ls[i],
            "data" = dat,
            "model" = fit,
            "summary" = summary(fit))
mod_out[[i]] <- tmp
}


# Conditional Effects Plot for interaction
ce <- conditional_effects(x=fit_seas)
ce


effects = "sg_norm:ghm_scale",
int_conditions = list(ghm_scale = ghmq),
re_formula = NA)


summary(lmer(Bhattacharyya_distance.total ~ 0 + comp*niche_type + (1|individual_id), data = dat_seas))
summary(dat_seas$Bhattacharyya_distance.total)
