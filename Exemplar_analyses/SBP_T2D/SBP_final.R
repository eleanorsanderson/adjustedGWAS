setwd("D:/Documents/Bristol PhD/Miniproject 2/Univariate MR analyses")
library(metafor)
library(plyr) 
library(meta) 
library(rmeta) 
library(TwoSampleMR)
library(MRInstruments)
library(MRPracticals)

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
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))
mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
#Extract SNPs for later analyses
unadjustedSNP<-data.frame(exposure_data$SNP)

#BMI estimate 



exposure_data<-extract_instruments("ukb-b-2303")
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ebi-a-GCST005413")
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))

mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

#adjusted SBP MR
#SD extracted from paper
SD<-20.7
exposure_data<-extract_instruments("ieu-b-38")
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ebi-a-GCST005413")
exposure_data$beta.exposure<-exposure_data$beta.exposure/SD
exposure_data$se.exposure<-exposure_data$se.exposure/SD
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))

mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
#Snps extracted for later analyses

#adjusted outcome SBP MR
exposure_data<-extract_instruments("ukb-b-20175")
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ebi-a-GCST005413")
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))
mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
#SNPs extracted for later analyses
unadjustedSNP<-data.frame(exposure_data$SNP)


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
SBPconfintupper<-(res$Estimate[1]+1.96*res$Std..Error[1])
SBPconfintlower<-(res$Estimate[1]-1.96*res$Std..Error[1])
BMIconfintupper<-(res$Estimate[2]+1.96*res$Std..Error[2])
BMIconfintlower<-(res$Estimate[2]-1.96*res$Std..Error[2])


#mvmradjusted


exposures<-mv_extract_exposures(c("ieu-b-38","ukb-b-2303") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
exposures$beta.exposure    <- ifelse(exposures$id.exposure == "ieu-b-38", exposures$beta.exposure/SD, exposures$beta.exposure)
exposures$se.exposure    <- ifelse(exposures$id.exposure == "ieu-b-38", exposures$se.exposure/SD, exposures$se.exposure)

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
SBPconfintupper<-(res$Estimate[1]+1.96*res$Std..Error[1])
SBPconfintlower<-(res$Estimate[1]-1.96*res$Std..Error[1])
BMIconfintupper<-(res$Estimate[2]+1.96*res$Std..Error[2])
BMIconfintlower<-(res$Estimate[2]-1.96*res$Std..Error[2])



#Reduced SNP MVMR analyses



#MVMR unadjusted
exposures<-mv_extract_exposures(c("ukb-b-20175","ukb-b-2303") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = unadjustedSNP$exposure_data.SNP, outcomes = "ebi-a-GCST005413")
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
SBPconfintupper<-(res$Estimate[1]+1.96*res$Std..Error[1])
SBPconfintlower<-(res$Estimate[1]-1.96*res$Std..Error[1])
BMIconfintupper<-(res$Estimate[2]+1.96*res$Std..Error[2])
BMIconfintlower<-(res$Estimate[2]-1.96*res$Std..Error[2])
#mvmradjusted


  exposures<-mv_extract_exposures(c("ieu-b-38","ukb-b-2303") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
exposures$beta.exposure    <- ifelse(exposures$id.exposure == "ieu-b-38", exposures$beta.exposure/SD, exposures$beta.exposure)
exposures$se.exposure    <- ifelse(exposures$id.exposure == "ieu-b-38", exposures$se.exposure/SD, exposures$se.exposure)

outcome_data <- extract_outcome_data(
  snps = adjustedSNP$exposure_data.SNP, outcomes = "ebi-a-GCST005413")
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
SBPconfintupper<-(res$Estimate[1]+1.96*res$Std..Error[1])
SBPconfintlower<-(res$Estimate[1]-1.96*res$Std..Error[1])
BMIconfintupper<-(res$Estimate[2]+1.96*res$Std..Error[2])
BMIconfintlower<-(res$Estimate[2]-1.96*res$Std..Error[2])

#BP

#unadj BMI 
exposure_data<-extract_instruments("ieu-b-40")
outcome_data <- extract_outcome_data(
  snps = unadjustedSNP$exposure_data.SNP, outcomes = "ukb-b-20175")
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))

mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

#adj BMI 
exposure_data<-extract_instruments("ieu-b-40")
outcome_data <- extract_outcome_data(
  snps = adjustedSNP$exposure_data.SNP, outcomes = "ukb-b-20175")
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))

mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))


