library(dplyr)
library(readr)
fullresults_med_pos <- list.files(path="./Positive_Mediator/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

write.csv(fullresults_med_pos, file = "fullresults_med_pos.csv")