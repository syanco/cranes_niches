expandESACCI <- function(dat){
  x <- dat$notes %>% 
    str_split(pattern = "\\|", simplify = F)
  
  
  out <- list()
  for(i in 1:length(x)){
    if(x[[i]][1] == "Values masked using quality control flag."){
      out[[i]] <- data.frame(QA = "masked")
    }else{
      
      out[[i]] <- x[i] %>% 
        lapply(., str_split, pattern = ":", simplify = T) %>% 
        as.data.frame() %>% 
        dplyr::rename("cat" = X1, "pixels" = X2) %>% 
        dplyr::mutate(pixels = as.numeric(pixels)) %>% 
        pivot_wider(names_from = cat, values_from = pixels)
    } #else
  } #i
  tot <- do.call("bind_rows", out)
  
  return(tot)
} #fxn

# Parallel version

expandESACCI_par <- function(dat){
  x <- dat$notes 
  
  out <- future_lapply(x, function(x_i){
    if(x_i[1] == "Values masked using quality control flag."){
      return(data.frame(QA = "masked"))
    } else {
      x_i %>% 
        str_split(pattern = "\\|", simplify = F) %>% 
        lapply(., str_split, pattern = ":", simplify = T) %>% 
        as.data.frame() %>% 
        dplyr::rename("cat" = X1, "pixels" = X2) %>% 
        dplyr::mutate(pixels = as.numeric(pixels)) %>% 
        pivot_wider(names_from = cat, values_from = pixels)
    }
  })
  
  tot <- do.call("bind_rows", out)
  
  return(tot)
}



# # test with example data
# dat <- data.frame(notes = rep(c("Values masked using quality control flag.", "A:1|B:2|C:3"), 1000))
# system.time(expandESACCI_parallel(dat))
