get_joint_dataset <- function(){
  library(dplyr)
  chemical <- read_csv("./TOS_Retriever/out/chemical_data.csv")
  isotope <- read_csv("./TOS_Retriever/out/isotopes_data.csv")
  structure <- read_csv("./TOS_Retriever/out/field_data.csv")
  
  #join the products all available traits data first
  dat = inner_join(chemical, isotope,  by = "sampleID") %>%
    unique %>%
    write_csv('./TOS_Retriever/out/field_traits_dataset.csv')
  
  # just the geolocalized data
  dat %>%
    inner_join(structure,  by = "individualID") %>%
    unique %>%
    write_csv('./TOS_Retriever/out/utm_dataset.csv')
  chemical %>%
    inner_join(structure,  by = "individualID") %>%
    unique %>%
    write_csv('./TOS_Retriever/out/utm_data_chemical.csv')
}
