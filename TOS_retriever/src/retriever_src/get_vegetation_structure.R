get_vegetation_structure <- function(){
  file_tos_coordinates = read_csv("./TOS_Retriever/tmp/filesToStack10098/stackedFiles/vst_perplotperyear.csv") %>%
    dplyr::select(c("plotID","plotType", "utmZone", "easting", "northing", "coordinateUncertainty", "nlcdClass", "elevation"))
  
  file_mapping = read_csv("./TOS_Retriever/tmp/filesToStack10098/stackedFiles/vst_mappingandtagging.csv") %>%
    dplyr::select(c("uid", "eventID", "domainID","siteID","plotID","subplotID",
             "nestedSubplotID","pointID","stemDistance","stemAzimuth",
             "cfcOnlyTag","individualID","supportingStemIndividualID","previouslyTaggedAs",
             "taxonID","scientificName"))
  pointID <- readr::read_csv("./TOS_retriever/dat/pointID_scheme.csv")
  dat = inner_join(file_mapping,file_tos_coordinates,  by = "plotID") %>%
    inner_join(pointID,  by = "pointID") %>%
    drop_na(stemAzimuth) %>%
    unique
  dat$easting = dat$easting + dat$easting_offset
  dat$northing = dat$northing + dat$northing_offset
  
  # get tree coordinates
  dat_apply <- dat %>%
    dplyr::select(c(stemDistance, stemAzimuth, easting, northing)) 
  coords <- apply(dat_apply,1,function(params)from_dist_to_utm(params[1],params[2], params[3], params[4])) %>%
    t %>%
    data.frame
  colnames(coords) <- c('UTM_E', 'UTM_N')
  field_tag <- cbind(dat, coords)
  
  apparent = read_csv("./TOS_Retriever/tmp/filesToStack10098/stackedFiles/vst_apparentindividual.csv") %>%
    dplyr::select("individualID", "stemDiameter", "height", "baseCrownHeight", 
                  "maxCrownDiameter",  "ninetyCrownDiameter") %>%
    group_by(individualID) %>%
    summarise_all(funs(max))
  
  crown_attributes = left_join(field_tag, apparent, by="individualID") %>%
    unique
  
  write_csv(field_tag, './TOS_Retriever/out/field_data_no_attribute.csv')
  write_csv(crown_attributes, './TOS_Retriever/out/field_data.csv')
}

