stratified_data_split <- function(in_dir = "./Model_builder/inputs/raw", 
                                  out_dir = "./Model_builder/inputs/spectra", seed = 1){
  # divide data in training, validation, and test
  library(readr)
  library(dplyr)
  set.seed(seed)
  features = read_csv(paste(in_dir, "CrownPix.csv", sep="/"))
  variables = read_csv(paste(in_dir, "CrownTraits.csv", sep="/"))
  allData <- inner_join(variables, features, by="pixel_crownID") %>%
    unique
  
  test_labels <- allData %>%
    group_by(pixel_crownID) %>%
    summarise_at(.vars = c("pixel_crownID", "name", "site"),.funs = c(unique="unique"))
  test_labels <- test_labels %>%
    group_by(name_unique, site_unique) %>%
    sample_frac(0.1)
  
  test_variables <- variables[variables$pixel_crownID %in% test_labels$pixel_crownID,]
  test_features <- features[features$pixel_crownID %in% test_variables$pixel_crownID,]
  train_variables <- variables[!(variables$pixel_crownID %in% test_labels$pixel_crownID), ]
  train_features <- features[features$pixel_crownID %in% train_variables$pixel_crownID, ]
  
    
  write_csv(test_features, paste(out_dir, '/CrownPix_outBag.csv', sep=""))
  write_csv(test_variables, paste(out_dir, '/CrownTraits_outBag.csv', sep=""))
  write_csv(train_features, paste(out_dir, '/CrownPix_norm.csv', sep=""))
  write_csv(train_variables, paste(out_dir, '/CrownTraits.csv', sep=""))
  
}
