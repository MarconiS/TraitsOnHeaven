stack_spectra = function(f_pt =  "./AOP_from_coords/outputs/outputs/Spectra/"){
  library(readr)
  library(tidyverse)
  files_to_stack <- list.files(f_pt, pattern = "csv")
  f0 <- read_csv(paste(f_pt, files_to_stack[1], sep="/"))
  colnames(f0) <- c("band_chm", paste("band", seq(1,369), sep="_"))
  data_id <- gsub('.{4}$', '', files_to_stack[1]) 
    
  ee <- lapply(strsplit(data_id, split = "_"), function(x) return(c(x[4], x[2], x[3], x[1]))) %>%
    data.frame(stringsAsFactors = F) %>%
    t 
  f0$individualID <- ee[,1]
  f0$band_site <- ee[,2]
  f0$band_species <- ee[,3]
  f0$flightpath <- ee[,4]
  for(f in files_to_stack[-1]){
    foo <-  read_csv(paste(f_pt, f, sep="/"))
    if(dim(foo)[1]!=0){
    colnames(foo) <- c("band_chm", paste("band", seq(1,369), sep="_"))
    data_id <- gsub('.{4}$', '', f)
    ee <- lapply(strsplit(data_id, split = "_"), function(x) return(c(x[4], x[2], x[3], x[1]))) %>%
      data.frame(stringsAsFactors = F) %>%
      t 
    foo$individualID <- ee[,1]
    foo$band_site <- ee[,2]
    foo$band_species <- ee[,3]
    foo$flightpath <- ee[,4]
    
    f0 <- rbind(f0, foo)
    }else{
      warning(paste(f, "is an empty file"))
    }
  }
  write_csv(f0, "./AOP_from_coords/outputs/neon_spectra.csv")
  write_csv(f0, "/ufrc/ewhite/s.marconi/Chapter1/jMIR/data/spectra/buffered_features.csv")
  return(f0)
}


