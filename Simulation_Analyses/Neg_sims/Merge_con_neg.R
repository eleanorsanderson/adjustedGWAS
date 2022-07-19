library(dplyr)
library(readr)
fullresults_con_neg <- list.files(path="./Negative_Confounder/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

write.csv(fullresults_con_neg, file = "fullresults_con_neg.csv")