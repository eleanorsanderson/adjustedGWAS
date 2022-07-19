library(dplyr)
library(readr)
fullresults_con_pos <- list.files(path="./Positive_Confounder/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

write.csv(fullresults_con_pos, file = "fullresults_con_pos.csv")