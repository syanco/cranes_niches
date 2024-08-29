library(tidyverse)


# extract phenology info from segmentaitons
species_list <- c("Anthropoides_virgo", "Grus_grus", "Grus_nigricollis", "Grus_vipio", "Balearica_pavonina")

segP <- "~/projects/dynamic_niches/data/shiny_track_segmentation/"

segFs <- list.files(segP, full.names = T)
# i <- 1
out_tmp <- list()
for(i in 1:length(species_list)){
  
  sp <- species_list[i]
  file_idx <- grep(sp, segFs)
  
  start_fall_ls <- list()
  stop_fall_ls <- list()
  
  start_spring_ls <- list()
  stop_spring_ls <- list()
  
  # sp_out <- list()
  for(j in 1:length(file_idx)){
  seg_temp <- read.csv(segFs[file_idx[j]])
  
  start_fall_ls[[j]] <- seg_temp %>% filter(Status == "Start Fall")
  stop_fall_ls[[j]] <- seg_temp %>% filter(Status == "End Fall")
  start_spring_ls[[j]] <- seg_temp %>% filter(Status == "Start Spring")
  stop_spring_ls[[j]] <- seg_temp %>% filter(Status == "End Spring")
  }

  out_tmp[[i]] <- rbind(do.call("rbind", start_fall_ls),
        do.call("rbind", stop_fall_ls),
        do.call("rbind", start_spring_ls),
        do.call("rbind", stop_spring_ls)
  )
}

phen <- do.call("rbind", out_tmp)


phen_summ <- phen %>% 
  mutate(wk= week(Date)) %>% 
  group_by(Species, Status) %>% 
summarize(med = median(wk),
          mn = mean(wk, na.rm = T))

write_csv(phen_summ, file = "out/phenology_summary.csv")
