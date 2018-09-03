stack_spectra = function(f_pt =  "./AOP_from_coords/outputs/outputs/Spectra/"){
  files_to_stack <- list.files(f_pt, pattern = "csv")
  f0 <- read_csv(paste(f_pt, files_to_stack[1], sep="/"))
  colnames(f0) <- c("band_chm", paste("band", seq(1,369), sep="_"))
  f0$individualID <- gsub('.{4}$', '', files_to_stack[1])
  for(f in files_to_stack[-1]){
    foo <-  read_csv(paste(f_pt, f, sep="/"))
    if(dim(foo)[1]!=0){
    colnames(foo) <- c("band_chm", paste("band", seq(1,369), sep="_"))
    foo$individualID <- gsub('.{4}$', '', f)
    
    f0 <- rbind(f0, foo)
    }else{
      warning(paste(f, "is an empty file"))
    }
  }
  write_csv(f0, "./AOP_from_coords/outputs/neon_spectra.csv")
  return(f0)
}


