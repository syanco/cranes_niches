load("./out/fitted_pcas_old.rdata")
old_pca <- pca_out

.wd <- getwd()
.dbPF <- file.path(.wd,'data/anno_move.db')
.dbPF <- file.path(.wd,'data/anno_move_copy2.db')

db <- dbConnect(RSQLite::SQLite(), .dbPF)
dbListTables(db)


event <- tbl(db, "event")
colnames(event)


event_new <- tbl(db, "event_new")
colnames(event_new)

event_clean <- tbl(db, "event_clean")
colnames(event_clean)

anno_mar28 <- tbl(db, "anno_join_2023-03-27")
colnames(anno_mar28)

event_clean_new <- tbl(db, "event_clean_new")
colnames(event_clean_new)


x <- evt0 %>% 
  mutate(ts = ymd_hms(timestamp),
         wk = week(ts)) %>% 
  group_by(individual_id) %>% 
  summarize(n = n(),
            species = species[1],
            min_ts = min(ts),
            max_ts = max(ts)) %>% 
  mutate(dur = difftime(max_ts, min_ts, units = 'weeks'))


plot(x = evt_wk$`Prop. Crops`, y = old_pca[[j]]$data$`Prop. Crops`)
abline(a=0, b=1, col="red", lwd=2, lty=2)
plot(x = evt_wk$`Prox. to Water`, y = old_pca[[j]]$data$`Prox. to Water`)
abline(a=0, b=1, col="red", lwd=2, lty=2)
plot(x = evt_wk$EVI, y = old_pca[[j]]$data$EVI)
abline(a=0, b=1, col="red", lwd=2, lty=2)
plot(x = evt_wk$Temp., y = old_pca[[j]]$data$Temp.)
abline(a=0, b=1, col="red", lwd=2, lty=2)
  
  
d <- tibble(x = c(1,1,1,2,2,2),
            y = c(0.5, 0.5, 0.5, 3, 3, 3))
str(d)

d %>% 
  group_by(x) %>% 
  summarise(mu = mean(y))


dur_summ <- evt0 %>% 
  group_by(individual_id) %>% 
  summarize(dur = difftime(max(timestamp), min(timestamp), units = "days"))

rm_ind2 <- dur_summ %>% 
  filter(dur < 50) %>% 
  pull(individual_id)

out <- evt0 %>% 
  filter(individual_id %notin% rm_ind2)
