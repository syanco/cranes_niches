# scratch dissim plots
library(ggplot2)

# FUNCTIONS

get_scenarios <- function () {
  resp <- httr::GET(paste0(BASE_URL, '/list/scenarios'))
  jsonlite::flatten(jsonlite::fromJSON(httr::content(resp, "text", encoding = enc), simplifyVector = T))
}

get_species_scenario <- function (sciname, product, variable, s_buff, t_buff, limit=50000) {
  resp <- httr::POST(paste0(BASE_URL, '/species/metrics/scatter'), body = list(
    mode = 'temporal',
    scientificname = sciname,
    limit= limit,
    yaxis = list(
      product = product,
      variable = variable,
      temporal = t_buff, 
      spatial = s_buff
    )
  ), encode = "json",)
  jsonlite::flatten(jsonlite::fromJSON(httr::content(resp, "text", encoding = enc))$rows)
}

`%notin%` <- Negate(`%in%`)
dissim <- read_csv("out/niche_dissim.csv")
breadth <- read_csv("out/niche_breadth.csv")

BASE_URL <- 'https://stoat-otf.azurewebsites.net'
enc = "UTF-8"
.dbPF <- file.path('data/anno_move.db')
.anno <- "anno_join_2022-11-16"
.ctfs <- file.path("ctfs/segmentations")


# Pull STOAT Pre-Annotations
#
# View possible annotations
# get_scenarios()
envs <- read_csv("ctfs/anno_vars.csv") %>% 
  filter(run == 1)

spp <- unique(dissim$species)

sp_out <- list()
for(i in 1:length(spp)){
  spp_env <- list()
  
  for(j in 1:nrow(envs)){
    spp_env[[j]] <- get_species_scenario(spp[i], envs$product[j], 
                                         envs$variable_code[j], envs$s_buff[j], 
                                         envs$t_buff[j]) %>% 
      mutate(variable = envs$variable[j])
  }
  
  sp_out[i] <- do.call("rbind", pop_env)
}

out <- do.call("rbind", sp_out)

ggplot(dissim) +
  geom_density(aes(x=ind_dis))+
  theme_minimal()+
  facet_wrap(~species, ncol = 1)


ggplot(breadth) +
  geom_density(aes(x=total, color = season))+
  theme_minimal()+
  xlab("Niche Breadth") +
  facet_wrap(~species, ncol = 1)
