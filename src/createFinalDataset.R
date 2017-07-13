# # geto to take the best 20 models, use it to define the pixel to use for the 5 methods to try ensamble (migrate it ot python?)
# rm(list=ls(all=TRUE))   # clear workspace
# rounds = 10
# NeonSite <- "TALL"
# setwd('~/Documents/Projects/TraitsOnHeaven/')
# out.dir = paste(getwd(), NeonSite,"outputs/", sep="/")
# in.dir = paste(getwd(),NeonSite,  "inputs/", sep="/")
# names <- c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct")
# path = paste(getwd(), NeonSite, sep="/")
# setwd(path)
# nentries = 40

storeFinalSet <- function(in.dir, out.dir, rounds = 10, names = c("LMA_g.m2", "d13C","d15N","C_pct","N_pct", "P_pct"),nentries = 40){
  bestPixes.r2 <- read_csv(paste(out.dir, 'P_pct_Boost_', rounds, '_R2.csv', sep = ""))
  for(j in c(seq(3,13, 2))){
    foo <- cbind(bestPixes.r2[,1], bestPixes.r2[,j])
    foo <- foo[order(foo[,2], decreasing = T),]  
    foo <- foo[1:nentries,]
    files.ls <-paste('onePix1Crown_', names[(j-1)/2], foo$X1, sep = '')
    file.create(paste(in.dir, 'FinalSet/',"finalAugmentedMatrix_", names[(j-1)/2], '.csv', sep = ""))
    for(f in files.ls){
      file.copy(list.files(paste(in.dir, 'Bootstrap_',rounds, sep = ""), f), paste(in.dir, 'FinalSet', sep=''))
      #now what you want is to gathere the right pixels in one single file per trait
      file.append(paste(in.dir, 'FinalSet/',"finalAugmentedMatrix_", names[(j-1)/2], '.csv', sep = ""), 
                  paste(in.dir, 'Bootstrap_',rounds,'/', f, '.csv',sep=""))
      
    }
  }
}
