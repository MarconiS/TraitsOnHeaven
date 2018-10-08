get_joint_dataset <- function(){
  library(dplyr)
  chemical <- read_csv("./TOS_Retriever/out/chemical_data.csv") %>%
    dplyr::select("individualID", colnames(chemical)[sapply(chemical, class)=="numeric"]) %>%
    dplyr::select(-one_of("elevation"))
  #isotope <- read_csv("./TOS_Retriever/out/isotopes_data.csv")
  structure <- read_csv("./TOS_Retriever/out/field_data.csv")
  
  #join the products all available traits data first
  # dat = inner_join(chemical, isotope) %>%
  #   unique %>%
  #   write_csv('./TOS_Retriever/out/field_traits_dataset.csv')
  # 
  # just the geolocalized data
  dat <- inner_join(chemical, structure) %>% 
    unique 
  write_csv(dat, './TOS_Retriever/out/utm_dataset.csv')
}
