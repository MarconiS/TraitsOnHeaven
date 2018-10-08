get_epsg_from_utm <- function(utm){
  dictionary <- cbind(32616, 32615, 32617, 32617, 32616, 32616, 32612, 32613, 32617, 32617) 
  colnames(dictionary) <- c("STEI", "CHEQ", "SCBI", "GRSM", "ORNL", "TALL", "MOAB", "JORN", "OSBS", "MLBS")
  return(dictionary[colnames(dictionary)==utm])
}
