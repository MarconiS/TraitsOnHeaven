# script to read .asd files
# Sarah Graves, November 15, 2016

# files stored in separate folder that is not part of the repo
# in total there are 20,005 files at 700 MB

library(dplyr)
library(prospectr)


#### FUNCTION "average_files" ####
# function arguments
# file_list = string vector of the file paths to average over
# output is a vector of the mean reflectance for all files in the file_list
average_files <- function(file_list){
  
  # open 1 file to check dimensions and get wavelengths
  s <- readASD(file_list[1],out_format = "list")[[1]]
  wl_um <- s$wavelength
  
  # create matrix to store spectra
  store <- matrix(nrow=length(file_list),ncol=length(wl_um))  
  
  # loop over each file in the list of files
  # extract spectra, store in matrix
  for (i in 1:length(file_list)){
    
    spec <- as.numeric(readASD(file_list[i]))
    store[i,] <- spec
    
  } # end loop of file names
  
  # average across the rows of the matrix - want an average of all rows for each column
  # other functions can be used in this step
  ave_spec <- apply(store,2,FUN=mean)
  
  names(ave_spec) <- wl_um
  
  return(ave_spec)
  
} # end function

#### LOAD DATA ####

# set reflectance file names
asd_files <- list.files("../../data/dimensions_field_data/2015/OSBS/spectra/sample_refl/")


# set up data frame with tree and leaf ids
df <- data.frame(tree = substr(asd_files,11,13),
                 leaf = substr(asd_files,15,15))

# check out data
summary(df)

unique(df$tree)
unique(df$leaf)

# create unique data frame of tree and leaf ids
# this is used to subset list of files
df_unique <- unique(df)


#### GENERATE LIST OF FILES TO AVERAGE ####

# what is the file name - generate pattern based on unique_df
pat <- paste("osbs_refl_",df_unique$tree[1],"_",df_unique$leaf[1],"_",sep="")

# generate file list for 1 row of unique_df
f <- list.files("../../data/dimensions_field_data/2015/OSBS/spectra/sample_refl/",
                pattern=pat,
                full.names = T)

#### APPLY FUNCTION TO EXTRACT LEAF REFLECTANCE ####

# apply extract and average function
# output is a named vector
t <- average_files(f)

plot(as.numeric(names(t)),t,
     main=paste0("tree-",df_unique$tree[1]," leaf-",df_unique$leaf[1]),
     ylim=c(0,.8))

#### EXTRACT LEAF REFL FOR ALL LEAVES ####

# loop over each row of the df_unique
# store average leaf reflectance in matrix

# set up matrix
leaf_refl <- matrix(NA,nrow=nrow(df_unique),ncol=length(t))

# loop
for (i in 1:nrow(df_unique)){

  # what is the file name - generate pattern based on unique_df
  pat <- paste("osbs_refl_",df_unique$tree[i],"_",df_unique$leaf[i],"_*",sep="")
  
  print(pat)
  
  # generate file list for 1 row of unique_df
  f <- list.files("../../data/dimensions_field_data/2015/OSBS/spectra/sample_refl/",
                  pattern=pat,
                  full.names = T)
  
  # apply extract and average function
  # output is a named vector
  leaf_refl[i,] <- average_files(f)
  
}

# add column names as wavelength
wl_names <- paste("um",names(average_files(f)),sep="_")
colnames(leaf_refl) <- wl_names

# convert decimal of reflectance to integer (scale factor of 10000)
# smaller data storage size
leaf_refl <- round(leaf_refl*10000)

# add matrix to data frame
df_leaf_refl <- data.frame(df_unique,leaf_refl)

# write to file
write.csv(df_leaf_refl,"Outputs/leaf_refl.csv")


