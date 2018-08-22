get_crown_area <- function(tos_data){
  data_no_crown <- tos_data[is.na(tos_data$maxCrownDiameter),]
  data_to_allometry <- tos_data[!is.na(tos_data$maxCrownDiameter),]
  maxCrown <- glm(maxCrownDiameter ~ stemDiameter, data = data_to_allometry )
  minCrown <- glm(ninetyCrownDiameter ~  stemDiameter, data = data_to_allometry )
  #pred_df <- data.frame(cbind(data_no_crown$siteID, data_no_crown$stemDiameter))
  #colnames(pred_df) <- c("siteID", "stemDiameter")
  new <- data.frame(stemDiameter = data_to_allometry$stemDiameter)
  data_to_allometry$maxCrownDiameter <- predict.glm(maxCrown, newdata = new)
  data_to_allometry$ninetyCrownDiameter <- predict.glm(minCrown, newdata = new)
  tos_data[!is.na(tos_data$maxCrownDiameter),] <- data_to_allometry
  return(tos_data)
}
