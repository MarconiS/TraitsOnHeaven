# impute where NAs
#import datasets
data_pre_processing <- function(){
  features_l3 <- readr::read_csv("AOP_from_coords/outputs/buffer_l3_spectra.csv")[-371]
  features_l1 <- readr::read_csv("/Users/sergiomarconi/Documents/GitHub/jMIR/outputs/buffer_neon_spectra.csv")[-371]
  
  mean.no.na <- function(x){mean(x, na.rm = T)}
  traits<- readr::read_csv("TOS_retriever/out/utm_dataset_nooutliers.csv") %>% #siteID.x, taxonID.x,
    dplyr::select(individualID, leafMassPerArea, ligninPercent, 
                  cellulosePercent, foliarPhosphorusConc, foliarCalciumConc, 
                  foliarMagnesiumConc, foliarSulfurConc, foliarManganeseConc, 
                  foliarIronConc, foliarPotassiumConc, extractCarotConc,
                  foliarCopperConc, foliarBoronConc, foliarZincConc, 
                  extractChlAConc, extractChlBConc, nitrogenPercent,carbonPercent,
                  CNratio) %>%
    group_by(individualID) %>%
    summarise_all(mean) %>%
    unique
  
  traits[-1] <- log(traits[-1])
  # #OR
  # tmp_log <- traits %>% 
  #   dplyr::select(-one_of(c("individualID","leafMassPerArea","foliarZincConc", 
  #                           "foliarPhosphorusConc","foliarManganeseConc"))) %>%
  #   mutate_all(log) 
  # traits[!colnames(traits) %in% c("individualID","leafMassPerArea","foliarZincConc", 
  #                                 "foliarPhosphorusConc","foliarManganeseConc")] <- tmp_log
  # 
  # Scale features to be centered in 0 and standardized
  mean_traits <-  apply(traits[-1],2, mean.no.na)
  sd_traits <-  apply(traits[-1],2, sd, na.rm = T)
  traits[-1] <- apply(traits[-1],2, scale)
  
  source("Model_builder/src/utility.R")
  library(dplyr)
  
  
  
  L3_set <- features_l3 %>%
    dplyr::select(-one_of(c("band_chm", "band_1","band_2","band_3", 
                            "band_4","band_5","band_6",
                            "band_7", "band_8", "band_362",
                            "band_363","band_364","band_365",
                            "band_366","band_367","band_368","band_369"))) %>%
    no.val.rm %>%
    na.rm %>%
    optical.filter %>%
    na.rm %>%
    na.omit %>%
    #spectral_transform %>%
    no.val.rm 
  
  L1_set <- features_l1 %>%
    dplyr::select(-one_of(c("band_chm", "band_1","band_2","band_3", 
                            "band_4","band_5","band_6",
                            "band_7", "band_8", "band_362",
                            "band_363","band_364","band_365",
                            "band_366","band_367","band_368","band_369"))) %>%
    no.val.rm %>%
    na.rm %>%
    optical.filter %>%
    na.rm %>%
    na.omit %>%
    #spectral_transform %>%
    no.val.rm 
  L1_set <- L1_set[!L1_set$individualID %in% L3_set$individualID,]
  L1_set <- L1_set[L1_set$individualID %in% traits$individualID,]
  
  final_features_set <- rbind(L1_set, L3_set)
  #scale features to be centered in zero and standardized
  col_use = c("individualID","flightpath", "band_species", "band_site")
  # final_features_set[!(colnames(final_features_set) %in% col_use)] <- 
  #   final_features_set[!(colnames(final_features_set) %in% col_use)] %>%
  #   #group_by(band_site) %>% 
  #   mutate_all(scale)
  # 
  final_features_set[final_features_set$band_site=="STEI","band_site"] <- "CHEQ"
  augmented_matrix <- inner_join(traits, final_features_set, by = "individualID")
  
  set.seed(1)
  test_labels <- augmented_matrix[which(colnames(augmented_matrix) %in% colnames(final_features_set))] %>%
    group_by(individualID) %>%
    summarise_at(.vars = c("individualID", "band_species", "band_site"),.funs = c(unique="unique"))
  test_labels <- test_labels %>%
    group_by(band_species_unique, band_site_unique) %>%
    sample_frac(0.2)
  
  augmented_matrix$band_species <- factor(augmented_matrix$band_species)
  augmented_matrix$band_site <- factor(augmented_matrix$band_site)
  test_data <- augmented_matrix[augmented_matrix$individualID %in% test_labels$individualID,]
  train_data <- augmented_matrix[!(augmented_matrix$individualID %in% test_data$individualID), ]
  
  readr::write_csv(test_data, 'Model_builder/dat/Crown_outBag.csv')
  readr::write_csv(train_data, 'Model_builder/dat/Crown_norm.csv')
}


