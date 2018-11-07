pls_transform <- function(traindata, nComp){
  
  #use the pls framework to transform the data
  Xt = traindata[grepl("band", names(traindata))]
  Yt = traindata[!grepl("band", names(traindata))] %>%
    select(-one_of(c("species_ID", "site_ID")))
  
  library(mice)
  imp <- mice(Yt[-1], m = 1, print = FALSE) %>%
    complete(1)
  
  mintPLS <- plsdepot::plsreg2(Xt, imp, comps = nComp)
  
  return(mintPLS)
}