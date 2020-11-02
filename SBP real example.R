setwd("D:/Documents/Bristol PhD/Miniproject 2/Univariate MR analyses")
library(metafor)
library(plyr) 
library(meta) 
library(rmeta) 
library(TwoSampleMR)
library(MRInstruments)
library(MRPracticals)
vignette("MRBase")

data(gwas_catalog)
ao<-available_outcomes()

#unadjusted SBP MR



exposure_data<-extract_instruments("ukb-b-20175")
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ebi-a-GCST005413")
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median"))

mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))


#BMI estimate 



exposure_data<-extract_instruments("ukb-b-2303")
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ebi-a-GCST005413")
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median"))

mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

#adjusted SBP MR

exposure_data<-extract_instruments("ieu-b-38")
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ebi-a-GCST005413")
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median"))

mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))



#MVMR unadjusted
exposures<-mv_extract_exposures(c("ukb-b-20175","ukb-b-2303") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = exposures$SNP, outcomes = "ebi-a-GCST005413")
H_data<-mv_harmonise_data(
  exposure_dat = exposures, 
  outcome_dat = outcome_data
)
betaexp<-data.frame(H_data[[c(1)]]) 
hdata<-data.frame(H_data)
library(data.table)
hdata<-setDT(hdata, keep.rownames = TRUE)[]
MVMRdata<-format_mvmr(BXGs = hdata[,c(2,3)],hdata[,8],seBXGs = hdata[,c(6,7)],hdata[,10],RSID=hdata[,1])
sres<-strength_mvmr(r_input=MVMRdata,gencov=0)
sres<-strength_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov=0)
res<-ivw_mvmr(r_input=MVMRdata)
res<-data.frame(res)
res
exp(cbind(coef(res[1,]), confint(res[1,])))

#mvmradjusted

exposures<-mv_extract_exposures(c("ieu-b-38","ukb-b-2303") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = exposures$SNP, outcomes = "ebi-a-GCST005413")
H_data<-mv_harmonise_data(
  exposure_dat = exposures, 
  outcome_dat = outcome_data
)
betaexp<-data.frame(H_data[[c(1)]]) 
hdata<-data.frame(H_data)
library(data.table)
hdata<-setDT(hdata, keep.rownames = TRUE)[]
MVMRdata<-format_mvmr(BXGs = hdata[,c(2,3)],hdata[,8],seBXGs = hdata[,c(6,7)],hdata[,10],RSID=hdata[,1])
sres<-strength_mvmr(r_input=MVMRdata,gencov=0)
sres<-strength_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov=0)
res<-ivw_mvmr(r_input=MVMRdata)
res<-data.frame(res)
res
