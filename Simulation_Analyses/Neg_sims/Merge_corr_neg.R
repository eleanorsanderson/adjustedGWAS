library(dplyr)
library(readr)
fullresults_corr_neg <- list.files(path="./Negative_Correlation/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

write.csv(fullresults_corr_neg, file = "fullresults_corr_neg.csv")