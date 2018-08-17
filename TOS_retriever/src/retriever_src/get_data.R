#prd = "10026"
get_data <- function(prd=NULL){
  req <- GET(paste("http://data.neonscience.org/api/v0/products/DP1.",prd,".001", sep=""))
  # make this JSON readable -> "text"
  unlink(paste("./TOS_Retriever/tmp/filesToStack", prd, "/*",sep=""))
  req.text <- content(req, as="text")
  # Flatten data frame to see available data. 
  avail <- fromJSON(req.text, simplifyDataFrame=T, flatten=T)
  sitesID <- unlist(avail$data$siteCodes$siteCode)
  
  for(id in sitesID){
    tryCatch({
      zipsByProduct(dpID=paste("DP1.", prd, ".001", sep=""), site=id, package="basic", savepath="./TOS_Retriever/tmp/", check.size=F)
    },error=function(e){})
  }
  stackByTable(filepath=paste("./TOS_Retriever/tmp/filesToStack", prd, "/",sep=""), folder=T)
  
}
