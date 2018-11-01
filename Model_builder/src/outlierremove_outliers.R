outlierKD <- function(dt, var) {
  var_name <- dt %>% dplyr::select(var) %>% unlist
  tot <- sum(!is.na(var_name))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  message("Outliers identified: ", na2 - na1, " from ", tot, " observations")
  message("Proportion (%) of outliers: ", (na2 - na1) / tot*100)
  #message("Mean of the outliers: ", mo)
  m2 <- mean(var_name, na.rm = T)
  message("Mean without removing outliers: ", m1)
  message("Mean if we remove outliers: ", m2)
  response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if(response == "y" | response == "yes"){
    #dt[as.character(substitute(var))] <- invisible(var_name)
    #assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    #message("Outliers successfully removed", "\n")
    return(var_name)
  } else{
    message("Nothing changed", "\n")
    return(invisible(var_name))
  }
}

dat <- readr::read_csv("TOS_retriever/out/utm_dataset.csv") %>%
  unique

traits <- readr::read_csv("TOS_retriever/out/utm_dataset.csv") %>% #siteID.x, taxonID.x,
  dplyr::select(individualID, leafMassPerArea, ligninPercent, 
                cellulosePercent, foliarPhosphorusConc, foliarPotassiumConc, foliarCalciumConc, 
                foliarMagnesiumConc, foliarSulfurConc, foliarManganeseConc, 
                foliarIronConc, foliarCopperConc, foliarBoronConc, foliarZincConc, 
                extractChlAConc, extractChlBConc, extractCarotConc, nitrogenPercent,carbonPercent,
                CNratio) 
#traits[-1] <- log(traits[-1])
out_filtered <- traits
for(ii in 1:19){
  out_filtered[[ii+1]] <- outlierKD(dt = traits[ii+1], var = colnames(traits)[ii+1])
}

readr::write_csv(out_filtered, "TOS_retriever/out/utm_dataset_nooutliers.csv")

mean.no.na <- function(x){mean(x, na.rm = T)}
traits_noout <- readr::read_csv("TOS_retriever/out/utm_dataset_nooutliers.csv") %>% #siteID.x, taxonID.x,
  group_by(individualID) %>%
  summarise_all(mean.no.na) %>%
  unique

summary(traits_noout)
