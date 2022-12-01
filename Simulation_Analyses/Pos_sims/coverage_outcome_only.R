setwd("/mnt/storage/scratch/bu19525/BristolPhD/MVMRproject/output/pos/outcome/")
library(ccostr)
library(ggplot2)
library(knitr)
library(parallel)
library(msm)
results_med1<-read.csv("/mnt/storage/scratch/bu19525/BristolPhD/MVMRproject/output/pos/outcome/results_med_pos_onlyout_1.csv")
results_con1<-read.csv("/mnt/storage/scratch/bu19525/BristolPhD/MVMRproject/output/pos/outcome/results_con_pos_onlyout_1.csv")
results_corr1<-read.csv("/mnt/storage/scratch/bu19525/BristolPhD/MVMRproject/output/pos/outcome/results_corr_pos_onlyout_1.csv")

results_true<-data.frame("bmi_effect" = 0.8, "sbp_effect" = 0.4, "sbpaj_effect" = 0.4, "mvmr.bmi_effect"=0.8, "mvmr.sbp_effect"=0.4,"mvmraj.bmi_effect"=0.8,"mvmraj.sbpaj_effect"=0.4)
row.names(results_true) <- "true_results"
#results_bias<-data.frame("bmi_effect" = (colMeans(results_med1$bmi_effect)), "sbp_effect" = (colMeans(results_med1$sbp_effect)), "sbpaj_effect" = (colMeans(results_med1$sbpaj_effect)), "mvmr.bmi_effect"=(colMeans(results_med1$mvmr.bmi_effect)), "mvmr.sbp_effect"=(colMeans(results_med1$mvmr.sbp_effect)),"mvmraj.bmi_effect"=(colMeans(results_med1$mvmraj.bmi_effect)),"mvmraj.sbpaj_effect"=(colMeans(results_med1$mvmraj.sbpaj_effect)))

results_bias_med<-data.frame(colMeans(results_med1[,c(2,4,6,8,10,13,15)]))
results_bias_con<-data.frame(colMeans(results_con1[,c(2,4,6,8,10,13,15)]))
results_bias_corr<-data.frame(colMeans(results_corr1[,c(2,4,6,8,10,13,15)]))
names(results_bias_med)[1] <- "med_results"
names(results_bias_corr)[1] <- "corr_results"
names(results_bias_con)[1] <- "con_results"
results_bias_med<- as.data.frame(t(results_bias_med))
results_bias_con<- as.data.frame(t(results_bias_con))
results_bias_corr<- as.data.frame(t(results_bias_corr))
results <- rbind(results_true, results_bias_med)
results <- rbind(results, results_bias_con)
results <- rbind(results, results_bias_corr)
#Bias
results_bias <- cbind(round(results[,c(1,4,6)] - 0.8,2), 
                      round(results[,c(2,3,5,7)] - 0.4,2))
#coverage
#BMI
results_bias_med<-data.frame(results_med1)
results_bias_con<-data.frame(results_con1)
results_bias_corr<-data.frame(results_corr1)

results_bias_med$bmi_effect_lower<-results_bias_med$bmi_effect - (1.96*results_bias_med$bmi_se)
results_bias_med$bmi_effect_upper<-results_bias_med$bmi_effect + (1.96*results_bias_med$bmi_se)
results_bias_med$mvmr.bmi_effect_lower<-results_bias_med$mvmr.bmi_effect - (1.96*results_bias_med$mvmr.bmi_se)
results_bias_med$mvmr.bmi_effect_upper<-results_bias_med$mvmr.bmi_effect + (1.96*results_bias_med$mvmr.bmi_se)
results_bias_med$mvmraj.bmi_effect_lower<-results_bias_med$mvmraj.bmi_effect - (1.96*results_bias_med$mvmraj.bmi_se)
results_bias_med$mvmraj.bmi_effect_upper<-results_bias_med$mvmraj.bmi_effect + (1.96*results_bias_med$mvmraj.bmi_se)
#sbp
results_bias_med$sbp_effect_lower<-results_bias_med$sbp_effect - (1.96*results_bias_med$sbp_se)
results_bias_med$sbp_effect_upper<-results_bias_med$sbp_effect + (1.96*results_bias_med$sbp_se)
results_bias_med$sbpaj_effect_lower<-results_bias_med$sbpaj_effect - (1.96*results_bias_med$sbpaj_se)
results_bias_med$sbpaj_effect_upper<-results_bias_med$sbpaj_effect + (1.96*results_bias_med$sbpaj_se)
results_bias_med$mvmr.sbp_effect_lower<-results_bias_med$mvmr.sbp_effect - (1.96*results_bias_med$mvmr.sbp_se)
results_bias_med$mvmr.sbp_effect_upper<-results_bias_med$mvmr.sbp_effect + (1.96*results_bias_med$mvmr.sbp_se)
results_bias_med$mvmraj.sbpaj_effect_lower<-results_bias_med$mvmraj.sbpaj_effect - (1.96*results_bias_med$mvmraj.sbpaj_se)
results_bias_med$mvmraj.sbpaj_effect_upper<-results_bias_med$mvmraj.sbpaj_effect + (1.96*results_bias_med$mvmraj.sbpaj_se)
#CON

results_bias_con$bmi_effect_lower<-results_bias_con$bmi_effect - (1.96*results_bias_con$bmi_se)
results_bias_con$bmi_effect_upper<-results_bias_con$bmi_effect + (1.96*results_bias_con$bmi_se)
results_bias_con$mvmr.bmi_effect_lower<-results_bias_con$mvmr.bmi_effect - (1.96*results_bias_con$mvmr.bmi_se)
results_bias_con$mvmr.bmi_effect_upper<-results_bias_con$mvmr.bmi_effect + (1.96*results_bias_con$mvmr.bmi_se)
results_bias_con$mvmraj.bmi_effect_lower<-results_bias_con$mvmraj.bmi_effect - (1.96*results_bias_con$mvmraj.bmi_se)
results_bias_con$mvmraj.bmi_effect_upper<-results_bias_con$mvmraj.bmi_effect + (1.96*results_bias_con$mvmraj.bmi_se)
#sbp
results_bias_con$sbp_effect_lower<-results_bias_con$sbp_effect - (1.96*results_bias_con$sbp_se)
results_bias_con$sbp_effect_upper<-results_bias_con$sbp_effect + (1.96*results_bias_con$sbp_se)
results_bias_con$sbpaj_effect_lower<-results_bias_con$sbpaj_effect - (1.96*results_bias_con$sbpaj_se)
results_bias_con$sbpaj_effect_upper<-results_bias_con$sbpaj_effect + (1.96*results_bias_con$sbpaj_se)
results_bias_con$mvmr.sbp_effect_lower<-results_bias_con$mvmr.sbp_effect - (1.96*results_bias_con$mvmr.sbp_se)
results_bias_con$mvmr.sbp_effect_upper<-results_bias_con$mvmr.sbp_effect + (1.96*results_bias_con$mvmr.sbp_se)
results_bias_con$mvmraj.sbpaj_effect_lower<-results_bias_con$mvmraj.sbpaj_effect - (1.96*results_bias_con$mvmraj.sbpaj_se)
results_bias_con$mvmraj.sbpaj_effect_upper<-results_bias_con$mvmraj.sbpaj_effect + (1.96*results_bias_con$mvmraj.sbpaj_se)
#Corr
results_bias_corr$bmi_effect_lower<-results_bias_corr$bmi_effect - (1.96*results_bias_corr$bmi_se)
results_bias_corr$bmi_effect_upper<-results_bias_corr$bmi_effect + (1.96*results_bias_corr$bmi_se)
results_bias_corr$mvmr.bmi_effect_lower<-results_bias_corr$mvmr.bmi_effect - (1.96*results_bias_corr$mvmr.bmi_se)
results_bias_corr$mvmr.bmi_effect_upper<-results_bias_corr$mvmr.bmi_effect + (1.96*results_bias_corr$mvmr.bmi_se)
results_bias_corr$mvmraj.bmi_effect_lower<-results_bias_corr$mvmraj.bmi_effect - (1.96*results_bias_corr$mvmraj.bmi_se)
results_bias_corr$mvmraj.bmi_effect_upper<-results_bias_corr$mvmraj.bmi_effect + (1.96*results_bias_corr$mvmraj.bmi_se)
#sbp
results_bias_corr$sbp_effect_lower<-results_bias_corr$sbp_effect - (1.96*results_bias_corr$sbp_se)
results_bias_corr$sbp_effect_upper<-results_bias_corr$sbp_effect + (1.96*results_bias_corr$sbp_se)
results_bias_corr$sbpaj_effect_lower<-results_bias_corr$sbpaj_effect - (1.96*results_bias_corr$sbpaj_se)
results_bias_corr$sbpaj_effect_upper<-results_bias_corr$sbpaj_effect + (1.96*results_bias_corr$sbpaj_se)
results_bias_corr$mvmr.sbp_effect_lower<-results_bias_corr$mvmr.sbp_effect - (1.96*results_bias_corr$mvmr.sbp_se)
results_bias_corr$mvmr.sbp_effect_upper<-results_bias_corr$mvmr.sbp_effect + (1.96*results_bias_corr$mvmr.sbp_se)
results_bias_corr$mvmraj.sbpaj_effect_lower<-results_bias_corr$mvmraj.sbpaj_effect - (1.96*results_bias_corr$mvmraj.sbpaj_se)
results_bias_corr$mvmraj.sbpaj_effect_upper<-results_bias_corr$mvmraj.sbpaj_effect + (1.96*results_bias_corr$mvmraj.sbpaj_se)

#COVERAGE
#BMI
results_bias_med$bmi_coverage<-NULL
for(i in 1:nrow(results_bias_med)){
if (results_bias_med$bmi_effect_lower[i]<=0.8 & results_bias_med$bmi_effect_upper[i]>=0.8){
results_bias_med$bmi_coverage[i]<-1} else {results_bias_med$bmi_coverage[i]<-0}}
results_bias_med$sbp_coverage<-NULL
for(i in 1:nrow(results_bias_med)){
if (results_bias_med$sbp_effect_lower[i]<=0.4 & results_bias_med$sbp_effect_upper[i]>=0.4){
results_bias_med$sbp_coverage[i]<-1} else {results_bias_med$sbp_coverage[i]<-0}}
results_bias_med$sbpaj_coverage<-NULL
for(i in 1:nrow(results_bias_med)){
if (results_bias_med$sbpaj_effect_lower[i]<=0.4 & results_bias_med$sbpaj_effect_upper[i]>=0.4){
results_bias_med$sbpaj_coverage[i]<-1} else {results_bias_med$sbpaj_coverage[i]<-0}}
results_bias_med$mvmr.bmi_coverage<-NULL
for(i in 1:nrow(results_bias_med)){
if (results_bias_med$mvmr.bmi_effect_lower[i]<=0.8 & results_bias_med$mvmr.bmi_effect_upper[i]>=0.8){
results_bias_med$mvmr.bmi_coverage[i]<-1} else {results_bias_med$mvmr.bmi_coverage[i]<-0}}
results_bias_med$mvmraj.bmi_coverage<-NULL
for(i in 1:nrow(results_bias_med)){
if (results_bias_med$mvmraj.bmi_effect_lower[i]<=0.8 & results_bias_med$mvmraj.bmi_effect_upper[i]>=0.8){
results_bias_med$mvmraj.bmi_coverage[i]<-1} else {results_bias_med$mvmraj.bmi_coverage[i]<-0}}
results_bias_med$mvmraj.sbpaj_coverage<-NULL
for(i in 1:nrow(results_bias_med)){
if (results_bias_med$mvmraj.sbpaj_effect_lower[i]<=0.4 & results_bias_med$mvmraj.sbpaj_effect_upper[i]>=0.4){
results_bias_med$mvmraj.sbpaj_coverage[i]<-1} else {results_bias_med$mvmraj.sbpaj_coverage[i]<-0}}
results_bias_med$mvmr.sbp_coverage<-NULL
for(i in 1:nrow(results_bias_med)){
if (results_bias_med$mvmr.sbp_effect_lower[i]<=0.4 & results_bias_med$mvmr.sbp_effect_upper[i]>=0.4){
results_bias_med$mvmr.sbp_coverage[i]<-1} else {results_bias_med$mvmr.sbp_coverage[i]<-0}}
results_coverage_med <- data.frame("sbp_effect" = (mean(results_bias_med$sbp_coverage, na.rm = T)),"sbpaj_effect" = (mean(results_bias_med$sbpaj_coverage, na.rm = T)),
                               "bmi_effect" = (mean(results_bias_med$bmi_coverage, na.rm = T)),
                              "mvmr.sbp_effect"  = (mean(results_bias_med$mvmr.sbp_coverage,  na.rm = T)),
                               "mvmr.bmi_effect"  = (mean(results_bias_med$mvmr.bmi_coverage,  na.rm = T)),
                               "mvmraj.sbpaj_effect"  = (mean(results_bias_med$mvmraj.sbpaj_coverage,  na.rm = T)),
                               "mvmraj.bmi_effect"  = (mean(results_bias_med$mvmraj.bmi_coverage,  na.rm = T)))
#Con

results_bias_con$bmi_coverage<-NULL
for(i in 1:1000){
if (results_bias_con$bmi_effect_lower[i]<=0.8 & results_bias_con$bmi_effect_upper[i]>=0.8){
results_bias_con$bmi_coverage[i]<-1} else {results_bias_con$bmi_coverage[i]<-0}}
results_bias_con$sbp_coverage<-NULL
for(i in 1:1000){
if (results_bias_con$sbp_effect_lower[i]<=0.4 & results_bias_con$sbp_effect_upper[i]>=0.4){
results_bias_con$sbp_coverage[i]<-1} else {results_bias_con$sbp_coverage[i]<-0}}
results_bias_con$sbpaj_coverage<-NULL
for(i in 1:nrow(results_bias_con)){
if (results_bias_con$sbpaj_effect_lower[i]<=0.4 & results_bias_con$sbpaj_effect_upper[i]>=0.4){
results_bias_con$sbpaj_coverage[i]<-1} else {results_bias_con$sbpaj_coverage[i]<-0}}
results_bias_con$mvmr.bmi_coverage<-NULL
for(i in 1:1000){
if (results_bias_con$mvmr.bmi_effect_lower[i]<=0.8 & results_bias_con$mvmr.bmi_effect_upper[i]>=0.8){
results_bias_con$mvmr.bmi_coverage[i]<-1} else {results_bias_con$mvmr.bmi_coverage[i]<-0}}
results_bias_con$mvmraj.bmi_coverage<-NULL
for(i in 1:1000){
if (results_bias_con$mvmraj.bmi_effect_lower[i]<=0.8 & results_bias_con$mvmraj.bmi_effect_upper[i]>=0.8){
results_bias_con$mvmraj.bmi_coverage[i]<-1} else {results_bias_con$mvmraj.bmi_coverage[i]<-0}}
results_bias_con$mvmraj.sbpaj_coverage<-NULL
for(i in 1:1000){
if (results_bias_con$mvmraj.sbpaj_effect_lower[i]<=0.4 & results_bias_con$mvmraj.sbpaj_effect_upper[i]>=0.4){
results_bias_con$mvmraj.sbpaj_coverage[i]<-1} else {results_bias_con$mvmraj.sbpaj_coverage[i]<-0}}
results_bias_con$mvmr.sbp_coverage<-NULL
for(i in 1:1000){
if (results_bias_con$mvmr.sbp_effect_lower[i]<=0.4 & results_bias_con$mvmr.sbp_effect_upper[i]>=0.4){
results_bias_con$mvmr.sbp_coverage[i]<-1} else {results_bias_con$mvmr.sbp_coverage[i]<-0}}
results_coverage_con <- data.frame("sbp_effect" = (mean(results_bias_con$sbp_coverage, na.rm = T)),"sbpaj_effect" = (mean(results_bias_con$sbpaj_coverage, na.rm = T)),
                               "bmi_effect" = (mean(results_bias_con$bmi_coverage, na.rm = T)),
                               "mvmr.sbp_effect"  = (mean(results_bias_con$mvmr.sbp_coverage,  na.rm = T)),
                               "mvmr.bmi_effect"  = (mean(results_bias_con$mvmr.bmi_coverage,  na.rm = T)),
                               "mvmraj.sbpaj_effect"  = (mean(results_bias_con$mvmraj.sbpaj_coverage,  na.rm = T)),
                               "mvmraj.bmi_effect"  = (mean(results_bias_con$mvmraj.bmi_coverage,  na.rm = T)))

#COVERAGE
#Corr
results_bias_corr$bmi_coverage<-NULL
for(i in 1:1000){
if (results_bias_corr$bmi_effect_lower[i]<=0.8 & results_bias_corr$bmi_effect_upper[i]>=0.8){
results_bias_corr$bmi_coverage[i]<-1} else {results_bias_corr$bmi_coverage[i]<-0}}
results_bias_corr$sbp_coverage<-NULL
for(i in 1:1000){
if (results_bias_corr$sbp_effect_lower[i]<=0.4 & results_bias_corr$sbp_effect_upper[i]>=0.4){
results_bias_corr$sbp_coverage[i]<-1} else {results_bias_corr$sbp_coverage[i]<-0}}
results_bias_corr$sbpaj_coverage<-NULL
for(i in 1:nrow(results_bias_corr)){
if (results_bias_corr$sbpaj_effect_lower[i]<=0.4 & results_bias_corr$sbpaj_effect_upper[i]>=0.4){
results_bias_corr$sbpaj_coverage[i]<-1} else {results_bias_corr$sbpaj_coverage[i]<-0}}
results_bias_corr$mvmr.bmi_coverage<-NULL
for(i in 1:1000){
if (results_bias_corr$mvmr.bmi_effect_lower[i]<=0.8 & results_bias_corr$mvmr.bmi_effect_upper[i]>=0.8){
results_bias_corr$mvmr.bmi_coverage[i]<-1} else {results_bias_corr$mvmr.bmi_coverage[i]<-0}}
results_bias_corr$mvmraj.bmi_coverage<-NULL
for(i in 1:1000){
if (results_bias_corr$mvmraj.bmi_effect_lower[i]<=0.8 & results_bias_corr$mvmraj.bmi_effect_upper[i]>=0.8){
results_bias_corr$mvmraj.bmi_coverage[i]<-1} else {results_bias_corr$mvmraj.bmi_coverage[i]<-0}}
results_bias_corr$mvmraj.sbpaj_coverage<-NULL
for(i in 1:1000){
if (results_bias_corr$mvmraj.sbpaj_effect_lower[i]<=0.4 & results_bias_corr$mvmraj.sbpaj_effect_upper[i]>=0.4){
results_bias_corr$mvmraj.sbpaj_coverage[i]<-1} else {results_bias_corr$mvmraj.sbpaj_coverage[i]<-0}}
results_bias_corr$mvmr.sbp_coverage<-NULL
for(i in 1:1000){
if (results_bias_corr$mvmr.sbp_effect_lower[i]<=0.4 & results_bias_corr$mvmr.sbp_effect_upper[i]>=0.4){
results_bias_corr$mvmr.sbp_coverage[i]<-1} else {results_bias_corr$mvmr.sbp_coverage[i]<-0}}
results_coverage_corr <- data.frame("sbp_effect" = (mean(results_bias_corr$sbp_coverage, na.rm = T)),
"sbpaj_effect" = (mean(results_bias_corr$sbpaj_coverage, na.rm = T)),
                               "bmi_effect" = (mean(results_bias_corr$bmi_coverage, na.rm = T)),
                               "mvmr.sbp_effect"  = (mean(results_bias_corr$mvmr.sbp_coverage,  na.rm = T)),
                               "mvmr.bmi_effect"  = (mean(results_bias_corr$mvmr.bmi_coverage,  na.rm = T)),
                               "mvmraj.sbpaj_effect"  = (mean(results_bias_corr$mvmraj.sbpaj_coverage,  na.rm = T)),
                               "mvmraj.bmi_effect"  = (mean(results_bias_corr$mvmraj.bmi_coverage,  na.rm = T)))
row.names(results_coverage_med) <- "med_results"
row.names(results_coverage_corr) <- "corr_results"
row.names(results_coverage_con) <- "con_results"

results_bias_onlyout_corr<-data.frame("sbp_effect" = (mean(results_bias_corr$sbp_effect, na.rm = T)),
"sbpaj_effect" = (mean(results_bias_corr$sbpaj_effect, na.rm = T)),
                               "bmi_effect" = (mean(results_bias_corr$bmi_effect, na.rm = T)),
                               "mvmr.sbp_effect"  = (mean(results_bias_corr$mvmr.sbp_effect,  na.rm = T)),
                               "mvmr.bmi_effect"  = (mean(results_bias_corr$mvmr.bmi_effect,  na.rm = T)),
                               "mvmraj.sbpaj_effect"  = (mean(results_bias_corr$mvmraj.sbpaj_effect,  na.rm = T)),
                               "mvmraj.bmi_effect"  = (mean(results_bias_corr$mvmraj.bmi_effect,  na.rm = T)))

results_bias_onlyout_med<-data.frame("sbp_effect" = (mean(results_bias_med$sbp_effect, na.rm = T)),
"sbpaj_effect" = (mean(results_bias_med$sbpaj_effect, na.rm = T)),
                               "bmi_effect" = (mean(results_bias_med$bmi_effect, na.rm = T)),
                               "mvmr.sbp_effect"  = (mean(results_bias_med$mvmr.sbp_effect,  na.rm = T)),
                               "mvmr.bmi_effect"  = (mean(results_bias_med$mvmr.bmi_effect,  na.rm = T)),
                               "mvmraj.sbpaj_effect"  = (mean(results_bias_med$mvmraj.sbpaj_effect,  na.rm = T)),
                               "mvmraj.bmi_effect"  = (mean(results_bias_med$mvmraj.bmi_effect,  na.rm = T)))

results_bias_onlyout_con<-data.frame("sbp_effect" = (mean(results_bias_con$sbp_effect, na.rm = T)),
"sbpaj_effect" = (mean(results_bias_con$sbpaj_effect, na.rm = T)),
                               "bmi_effect" = (mean(results_bias_con$bmi_effect, na.rm = T)),
                               "mvmr.sbp_effect"  = (mean(results_bias_con$mvmr.sbp_effect,  na.rm = T)),
                               "mvmr.bmi_effect"  = (mean(results_bias_con$mvmr.bmi_effect,  na.rm = T)),
                               "mvmraj.sbpaj_effect"  = (mean(results_bias_con$mvmraj.sbpaj_effect,  na.rm = T)),
                               "mvmraj.bmi_effect"  = (mean(results_bias_con$mvmraj.bmi_effect,  na.rm = T)))


row.names(results_bias_onlyout_med) <- "med_results"
row.names(results_bias_onlyout_corr) <- "corr_results"
row.names(results_bias_onlyout_con) <- "con_results"
results_effect<-rbind(results_bias_onlyout_med,results_bias_onlyout_corr)
results_effect<-rbind(results_effect,results_bias_onlyout_con)
write.csv(results_effect, "results_effect_only_out.csv")



write.csv(results_coverage_con, "Coverage_results_con_onlyout.csv")
write.csv(results_coverage_corr, "Coverage_results_corr_onlyout.csv")
write.csv(results_coverage_med, "Coverage_results_med_onlyout.csv")

results_coverage<-rbind(results_coverage_med,results_coverage_corr)
results_coverage<-rbind(results_coverage,results_coverage_con)
write.csv(results_coverage, "results_coverage_only_out.csv")

