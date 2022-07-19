setwd("D:/Documents/Bristol PhD/Miniproject 2/Univariate MR analyses")
library(metafor)
library(plyr) 
library(meta) 
library(rmeta) 
library(TwoSampleMR)
library(MRInstruments)
library(MRPracticals)
library(MVMR)

data(gwas_catalog)
ao<-available_outcomes()

#unadjusted WC MR
exposure_data<-extract_instruments("ieu-a-61")
#SD
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ukb-b-20175")
#SD
unadjustedSNP<-data.frame(exposure_data$SNP)
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))

mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
#BMI estimate 
exposure_data<-extract_instruments("ieu-b-40")
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ukb-b-20175")
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))
mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

#adjusted WC MR
exposure_data<-extract_instruments("ieu-a-66")
#SD
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ukb-b-20175")
#Extract SNPs used for later analyses
adjustedSNP<-data.frame(exposure_data$SNP)
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))

mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
#adjusted outcome WC
#SD extracted from paper
SD<-20.7
exposure_data<-extract_instruments("ieu-a-61")
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ieu-b-38")
outcome_data$beta.outcome<-outcome_data$beta.outcome/SD
outcome_data$se.outcome<-outcome_data$se.outcome/SD
unadjustedSNP<-data.frame(exposure_data$SNP)
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))
mr_results
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))


#MVMR unadjusted
exposures<-mv_extract_exposures(c("ieu-a-61","ieu-b-40") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = exposures$SNP, outcomes = "ukb-b-20175")
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

WCconfintupper<-(res$Estimate[1]+1.96*res$Std..Error[1])
WCconfintlower<-(res$Estimate[1]-1.96*res$Std..Error[1])
BMIconfintupper<-(res$Estimate[2]+1.96*res$Std..Error[2])
BMIconfintlower<-(res$Estimate[2]-1.96*res$Std..Error[2])

#mvmradjusted

exposures<-mv_extract_exposures(c("ieu-a-66","ieu-b-40") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = exposures$SNP, outcomes = "ukb-b-20175")
H_data<-mv_harmonise_data(
  exposure_dat = exposures, 
  outcome_dat = outcome_data
)
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
WCconfintupper<-(res$Estimate[1]+1.96*res$Std..Error[1])
WCconfintlower<-(res$Estimate[1]-1.96*res$Std..Error[1])
BMIconfintupper<-(res$Estimate[2]+1.96*res$Std..Error[2])
BMIconfintlower<-(res$Estimate[2]-1.96*res$Std..Error[2])


#mvmradjustedoutcome

exposures<-mv_extract_exposures(c("ieu-a-61","ieu-b-40") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = exposures$SNP, outcomes = "ieu-b-38")
outcome_data$beta.outcome<-outcome_data$beta.outcome/SD
outcome_data$se.outcome<-outcome_data$se.outcome/SD
outcome_data<-outcome_data[-nrow(outcome_data),]
H_data<-mv_harmonise_data(
  exposure_dat = exposures, 
  outcome_dat = outcome_data
)
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
WCconfintupper<-(res$Estimate[1]+1.96*res$Std..Error[1])
WCconfintlower<-(res$Estimate[1]-1.96*res$Std..Error[1])
BMIconfintupper<-(res$Estimate[2]+1.96*res$Std..Error[2])
BMIconfintlower<-(res$Estimate[2]-1.96*res$Std..Error[2])


#reduced SNP analyses

#MVMR unadjusted
exposures<-mv_extract_exposures(c("ieu-a-61","ieu-b-40") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = unadjustedSNP$exposure_data.SNP, outcomes = "ukb-b-20175")
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
WCconfintupper<-(res$Estimate[1]+1.96*res$Std..Error[1])
WCconfintlower<-(res$Estimate[1]-1.96*res$Std..Error[1])
BMIconfintupper<-(res$Estimate[2]+1.96*res$Std..Error[2])
BMIconfintlower<-(res$Estimate[2]-1.96*res$Std..Error[2])
#mvmradjusted

exposures<-mv_extract_exposures(c("ieu-a-66","ieu-b-40") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = adjustedSNP$exposure_data.SNP, outcomes = "ukb-b-20175")
H_data<-mv_harmonise_data(
  exposure_dat = exposures, 
  outcome_dat = outcome_data
)
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

WCconfintupper<-(res$Estimate[1]+1.96*res$Std..Error[1])
WCconfintlower<-(res$Estimate[1]-1.96*res$Std..Error[1])
BMIconfintupper<-(res$Estimate[2]+1.96*res$Std..Error[2])
BMIconfintlower<-(res$Estimate[2]-1.96*res$Std..Error[2])





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
