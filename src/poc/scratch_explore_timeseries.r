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
    library(MVNH)
    library(INLA)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

# vector of seasons we'll need later
seasons <- c("Winter", "Spring", "Summer", "Fall")

#---- Load control files ----#

hab_vars <- c("dist2water_scale", 
              # "pfor_h_scale",
              # "evi_scale", 
              "ghm_scale", 
              # "altitude_scale", 
              "prop_crop_scale"
)

grin_vars <- c(
  # "dist2water_scale", 
  # "pfor_h_scale",
  "evi_scale",
  # "ghm_scale", 
  "altitude_scale"
  # "prop_crop_scale"
)


#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

# create event table
evt0 <- tbl(db, .table) %>% 
  collect()


evt1 <- evt0 %>% 
  mutate(grp = glue("{individual_id}-{yr}"),
         season = fct_relevel(season, "Spring", "Summer", "Fall")) %>% 
  pivot_longer(cols = c("dist2water_scale", "altitude_scale", "prop_crop_scale", "evi_scale", "ghm_scale"), values_to = "niche_val")



# formula.fixed <-  0 + used ~  var1 + 
#   f(stratum, model="iid", hyper = list(theta = list(initial = log(1e-6),
#                                                     fixed=T))) 

evt_av <- evt0 %>% 
  mutate(grp = glue("{individual_id}-{yr}"),
         ind_f = as.factor(individual_id),
         ind_num = as.numeric(ind_f)) %>% 
  filter(species == "Anthropoides virgo")

# INLA model
# form <- evi_scale ~ season + f(individual_id, model="iid") 

# mod_inla <- inla(evi_scale ~ 0 + season + f(individual_id, model="iid") + f(doy, replicate = individual_id, model="ar1"), 
                 # data=evt_gg)
mod_evi_av <- inla(evi_scale ~ 0 + season + f(ind_num, model="iid") + f(doy, replicate = ind_num, model="ar1"),
                   control.family = list(prec = list(param = c(1, 0.01))), 
                   formula_tau = 0+season,
                 data=evt_av)
summary(mod_evi_av)

evt0 %>% 
  group_by(species) %>% 
  summarise(n = n())

evt_av <- evt0 %>% 
  mutate(grp = glue("{individual_id}-{yr}"),
         ind_f = as.factor(individual_id),
         ind_num = as.numeric(ind_f),
         t_num = as.numeric(timestamp)) %>% 
  filter(species == "Anthropoides virgo")

# INLA model
# form <- evi_scale ~ season + f(individual_id, model="iid") 

# mod_inla <- inla(evi_scale ~ 0 + season + f(individual_id, model="iid") + f(doy, replicate = individual_id, model="ar1"), 
# data=evt_gg)
mod_evi_gg <- inla(evi_scale ~ 0 + season + f(ind_num, model="iid") + f(doy, replicate = ind_num, model="ar1"), 
                   data=evt_gg)
summary(mod_evi_gg)


# BRMS Model
library(brms)
library(tictoc)
conflicts_prefer(purrr::accumulate)
conflicts_prefer(brms::ar)
conflicts_prefer(tidyr::expand)
conflicts_prefer(tibble::has_name)
conflicts_prefer(tidyr::pack)
conflicts_prefer(tidyr::unpack)
conflicts_prefer(purrr::when)

tic()
fit_evi <- brm(bf(evi_scale ~ 0 + season + (1|ind_f),
              # + ar(time = doy, gr = ind_f, cov = F), 
                   sigma ~ 0 + season, family = gaussian()), data = evt_av, 
                chains = 3, iter = 1000, cores = getOption("mc.cores", 3))
toc()

tic()
fit_d2w <- brm(bf(dist2water_scale ~ 0 + season + (1|ind_f),
                  # + ar(time = doy, gr = ind_f, cov = F), 
                  sigma ~ 0 + season, family = gaussian()), data = evt_av, 
               chains = 3, iter = 1000, cores = getOption("mc.cores", 3))
toc()




fit
conditional_effects(fit)


(raw <- evt1 %>% 
  # filter(species == "Grus grus") %>%
  mutate(doy = yday(timestamp)) %>% 
  ggplot(aes(x = doy, y = niche_val, color = season)) +
  geom_point()+
  facet_grid(
    species
    ~ name,
    scales = "free") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank())
)

ggsave(filename = "out/draft_figs/raw_plots.png", raw)


# test model