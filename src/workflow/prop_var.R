# Summarize PCAs

library(glue)
library(ggfortify)

load(file = "out/fitted_pcas.rdata")

pca_out <- pca_out[!sapply(pca_out,is.null)]

for(i in 1:length(pca_out)){
  print(pca_out[[i]]$species)
  print(glue("PC 1: {summary(pca_out[[i]]$pca)$importance['Proportion of Variance', 1]}"))
  print(glue("PC 2: {summary(pca_out[[i]]$pca)$importance['Proportion of Variance', 2]}"))
  print(glue("PC 1+2: {summary(pca_out[[i]]$pca)$importance['Cumulative Proportion', 2]} \n
             \n"))
  
  
}

load(file = "out/fitted_pcas_var.rdata")

pca_out <- pca_out[!sapply(pca_out,is.null)]

for(i in 1:length(pca_out)){
  print(pca_out[[i]]$species)
  print(glue("PC 1: {summary(pca_out[[i]]$pca)$importance['Proportion of Variance', 1]}"))
  print(glue("PC 2: {summary(pca_out[[i]]$pca)$importance['Proportion of Variance', 2]}"))
  print(glue("PC 1+2: {summary(pca_out[[i]]$pca)$importance['Cumulative Proportion', 2]} \n
             \n"))
  
  
}

