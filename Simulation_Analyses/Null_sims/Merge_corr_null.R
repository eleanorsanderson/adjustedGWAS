library(dplyr)
library(readr)
fullresults_corr_pos <- list.files(path="./Null_Correlation/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

write.csv(fullresults_corr_pos, file = "fullresults_corr_null.csv")