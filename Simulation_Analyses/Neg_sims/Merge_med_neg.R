library(dplyr)
library(readr)
fullresults_med_neg <- list.files(path="./Negative_Mediator/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

write.csv(fullresults_med_neg, file = "fullresults_med_neg.csv")