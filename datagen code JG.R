#function to generate data under different scenarios
#input (confounder, mediator, correlated) - relates to the position of the adjustment variable X2 relative to X1 and the outcome. 
#requires functions simulateGP and dplyr
setwd("D:/Documents/Bristol PhD/Miniproject 2/Univariate MR analyses")
install.packages("ellipsis")
devtools::install_github("explodecomputer/simulateGP")
library(dplyr)
library(simulateGP)
library(TwoSampleMR)
library(MVMR)
  
  
datagen <- function(model,nsnps,beta1,beta2){

n=100000
effs_x1 <- abs(rnorm(snps/2,0,0.1))
effs_x2 <- abs(rnorm(snps/2,0,0.1))

G <- make_geno(n,snps,0.4)
v_x1 <- rnorm(n,0,1)
v_x2 <- rnorm(n,0,1)
C <- rnorm(n,0,1)

#second sample
G2 <- make_geno(n,snps,0.4)
v_x12 <- rnorm(n,0,1)
v_x22 <- rnorm(n,0,1)
C2 <- rnorm(n,0,1)

ind_exp <- data.frame(G)
ind_exp <- ind_exp %>% rename_at(vars(starts_with("X")), 
                                       funs(str_replace(., "X", "G.X")))
ind_out <- data.frame(G2)
ind_out <- ind_out %>% rename_at(vars(starts_with("X")), 
                                 funs(str_replace(., "X", "G.Y")))

if(model == "confounder"){
  
ind_exp[,"X2"] <- G[,1:(snps/2)]%*%(effs_x2) + C + v_x2
ind_exp[,"X1"]  <- G[,((snps/2)+1):snps]%*%effs_x1 + 0.25*ind_exp[,"X2"] + 0.8*C  + v_x1

x22 <- G2[,1:(snps/2)]%*%(effs_x2) + C2 + v_x22
x12 <- G2[,((snps/2)+1):snps]%*%effs_x1 + 0.25*x22 + 0.8*C2  + v_x12

}

if(model == "correlated"){
  
  ind_exp[,"X2"] <- G[,1:(snps/2)]%*%(effs_x2) + C + v_x2
  ind_exp[,"X1"]  <- G[,((snps/2)+1):snps]%*%effs_x1 + 0.8*C  + v_x1
  
  x22 <- G2[,1:(snps/2)]%*%(effs_x2) + C2 + v_x22
  x12 <- G2[,((snps/2)+1):snps]%*%effs_x1 + 0.8*C2  + v_x12
  
}


if(model == "mediator"){
 
  ind_exp[,"X1"]  <- G[,((snps/2)+1):snps]%*%effs_x1 +  0.8*C  + v_x1
  ind_exp[,"X2"] <- G[,1:(snps/2)]%*%(effs_x2) + C + 0.25*ind_exp[,"X1"] + v_x2
  
  x12 <- G2[,((snps/2)+1):snps]%*%effs_x1 + + 0.8*C2  + v_x12
  x22 <- G2[,1:(snps/2)]%*%(effs_x2) + C2 + 0.25*x12 + v_x22
  
}


ind_out[,"Y"] <- beta1*x12 + beta2*x22 + 0.5*C2 + rnorm(n,0,1)  

data <- cbind.data.frame(ind_exp,ind_out)
return(data)
}


#test
#confounder
Condata<-datagen("confounder", 250, 0.6, -0.5)
Corrdata<- datagen("correlated", 250, 0.6, -0.5)
Meddata<-datagen("mediator", 250, 0.6, -0.5)

Condata<-tibble::rowid_to_column(Condata, "ID")
Corrdata<-tibble::rowid_to_column(Corrdata, "ID")
Meddata<-tibble::rowid_to_column(Meddata, "ID")

#EXTRACT SNPS
#COnfounding data structure
X1SNP<- get_effs(Condata$X1,Condata$Y,Condata[,c(1:250)],xname="X1",yname="Y")
X2SNP<- get_effs(Condata$X2,Condata$Y,Condata[,c(1:250)],xname="X2",yname="Y")
mr(subset(X1SNP, pval.exposure < 5e-8))
mr(subset(X2SNP, pval.exposure < 5e-8))

mvdat <- make_mvdat(list(Condata$X1, Condata$X2), Condata$Y,Condata[,c(2:253)])
# Perform MV MR
mv_multiple(mvdat,pval_threshold = 5e-08)

#EXTRACT SNPS
#Correlating data structure
X1SNP<- get_effs(Corrdata$X1,Corrdata$Y,Corrdata[,c(1:250)],xname="X1",yname="Y")
X2SNP<- get_effs(Corrdata$X2,Corrdata$Y,Corrdata[,c(1:250)],xname="X2",yname="Y")
mr(subset(X1SNP, pval.exposure < 5e-8))
mr(subset(X2SNP, pval.exposure < 5e-8))


mvdat <- make_mvdat(list(Corrdata$X1, Corrdata$X2), Corrdata$Y,Corrdata[,c(2:251)] )
# Perform MV MR
mv_multiple(mvdat,pval_threshold = 5e-08)

#EXTRACT SNPS
#Mediating data structure
X1SNP<- get_effs(Meddata$X1,Meddata$Y,Meddata[,c(1:250)],xname="X1",yname="Y")
X2SNP<- get_effs(Meddata$X2,Meddata$Y,Meddata[,c(1:250)],xname="X2",yname="Y")
mr(subset(X1SNP, pval.exposure < 5e-8))
mr(subset(X2SNP, pval.exposure < 5e-8))


mvdat <- make_mvdat(list(Meddata$X1, Meddata$X2), Meddata$Y,Meddata[,c(1:250)] )
# Perform MV MR
mv_multiple(mvdat,pval_threshold = 5e-08)


#For loop

out_start=252
out_end= 252
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=2
exp_end=251
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

library(lme4)
  for (j in exp_start:exp_end){
    exposure = colnames(Condata)[j]
    model <- lm(Condata$X1 ~ get(exposure) ,
                  na.action = na.exclude,
                  data=Condata)
    
    Vcov <- vcov(model, useScale = FALSE)
    beta <- coef(model)
    se <- sqrt(diag(Vcov))
    zval <- beta / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)

    
    exp_beta[number] = as.numeric(beta[2])
    exp_se[number] = as.numeric(se[2])
    exp_pvalue[number] = as.numeric(pval[2])
    exp_variable[number] = exposure
    number = number + 1
  }


outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)


library(tidyverse)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
  )
allx2 = rbind(outcome, exposure)
allx2 = na.omit(allx2)
allx2<-tibble::rowid_to_column(allx2, "ID")

# X1

out_start=253
out_end= 253
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=2
exp_end=251
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

library(lme4)
for (j in exp_start:exp_end){
  exposure = colnames(Condata)[j]
  model <- lm(Condata$X2 ~ get(exposure) ,
              na.action = na.exclude,
              data=Condata)
  
  Vcov <- vcov(model, useScale = FALSE)
  beta <- coef(model)
  se <- sqrt(diag(Vcov))
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  
  exp_beta[number] = as.numeric(beta[2])
  exp_se[number] = as.numeric(se[2])
  exp_pvalue[number] = as.numeric(pval[2])
  exp_variable[number] = exposure
  number = number + 1
}


outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)


library(tidyverse)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
  )
allx1 = rbind(outcome, exposure)
allx1 = na.omit(allx1)
allx1<-tibble::rowid_to_column(allx1, "ID")
#Y


out_start=504
out_end= 504
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=254
exp_end=503
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

library(lme4)
for (j in exp_start:exp_end){
  exposure = colnames(Condata)[j]
  model <- lm(Condata$Y ~ get(exposure) ,
              na.action = na.exclude,
              data=Condata)
  
  Vcov <- vcov(model, useScale = FALSE)
  beta <- coef(model)
  se <- sqrt(diag(Vcov))
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  
  exp_beta[number] = as.numeric(beta[2])
  exp_se[number] = as.numeric(se[2])
  exp_pvalue[number] = as.numeric(pval[2])
  exp_variable[number] = exposure
  number = number + 1
}


outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)


library(tidyverse)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
  )
ally = rbind(outcome, exposure)
ally = na.omit(ally)
ally<-tibble::rowid_to_column(ally, "ID")

#x1+x2 

corr12<-merge(allx1,allx2, by="variable")
names(corr12)[3]<-"beta.x1"
names(corr12)[4]<-"se.x1"
names(corr12)[5]<-"pvalue.x1"
names(corr12)[7]<-"beta.x2"
names(corr12)[8]<-"se.x2"
names(corr12)[9]<-"pvalue.x2"
names(corr12)[9]<-"pvalue.x2"
names(corr12)[2]<-"ID"
corr12y<-merge(corr12,ally, by="ID")
corr12y<-subset(corr12y, pvalue.x1<5*10^-8  | pvalue.x2 <5*10^-8)
#format
F.data<-format_mvmr(BXGs=corr12y[,c(3,7)],
BYG=corr12y[,11],
seBXGs=corr12y[,c(4,8)],
seBYG=corr12y[,12],
RSID=corr12y[,1])
head(F.data)
sres<-strength_mvmr(r_input=F.data,gencov=0)
pres<-pleiotropy_mvmr(r_input=F.data,gencov=0)
res<-ivw_mvmr(r_input=F.data)

#Mediated


out_start=252
out_end= 252
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=2
exp_end=251
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

library(lme4)
for (j in exp_start:exp_end){
  exposure = colnames(Meddata)[j]
  model <- lm(Meddata$X1 ~ get(exposure) ,
              na.action = na.exclude,
              data=Meddata)
  
  Vcov <- vcov(model, useScale = FALSE)
  beta <- coef(model)
  se <- sqrt(diag(Vcov))
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  
  exp_beta[number] = as.numeric(beta[2])
  exp_se[number] = as.numeric(se[2])
  exp_pvalue[number] = as.numeric(pval[2])
  exp_variable[number] = exposure
  number = number + 1
}


outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)


outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
  )
mx2 = rbind(outcome, exposure)
mx2 = na.omit(mx2)
mx2<-tibble::rowid_to_column(mx2, "ID")

# X1

out_start=253
out_end= 253
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=2
exp_end=251
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

library(lme4)
for (j in exp_start:exp_end){
  exposure = colnames(Meddata)[j]
  model <- lm(Meddata$X2 ~ get(exposure) ,
              na.action = na.exclude,
              data=Meddata)
  
  Vcov <- vcov(model, useScale = FALSE)
  beta <- coef(model)
  se <- sqrt(diag(Vcov))
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  
  exp_beta[number] = as.numeric(beta[2])
  exp_se[number] = as.numeric(se[2])
  exp_pvalue[number] = as.numeric(pval[2])
  exp_variable[number] = exposure
  number = number + 1
}


outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)


library(tidyverse)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
  )
mx1 = rbind(outcome, exposure)
mx1 = na.omit(mx1)
mx1<-tibble::rowid_to_column(mx1, "ID")
#Y


out_start=504
out_end= 504
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=254
exp_end=503
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

library(lme4)
for (j in exp_start:exp_end){
  exposure = colnames(Meddata)[j]
  model <- lm(Meddata$Y ~ get(exposure) ,
              na.action = na.exclude,
              data=Meddata)
  
  Vcov <- vcov(model, useScale = FALSE)
  beta <- coef(model)
  se <- sqrt(diag(Vcov))
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  
  exp_beta[number] = as.numeric(beta[2])
  exp_se[number] = as.numeric(se[2])
  exp_pvalue[number] = as.numeric(pval[2])
  exp_variable[number] = exposure
  number = number + 1
}


outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)


library(tidyverse)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
  )
my = rbind(outcome, exposure)
my = na.omit(my)
my<-tibble::rowid_to_column(my, "ID")

#x1+x2 

mx12<-merge(mx1,mx2, by="variable")
names(mx12)[3]<-"beta.x1"
names(mx12)[4]<-"se.x1"
names(mx12)[5]<-"pvalue.x1"
names(mx12)[7]<-"beta.x2"
names(mx12)[8]<-"se.x2"
names(mx12)[9]<-"pvalue.x2"
names(mx12)[9]<-"pvalue.x2"
names(mx12)[2]<-"ID"
mx12y<-merge(mx12,my, by="ID")
mx12y<-subset(mx12y, pvalue.x1<5*10^-8  | pvalue.x2 <5*10^-8)
#format
F.data<-format_mvmr(BXGs=mx12y[,c(3,7)],
                    BYG=mx12y[,11],
                    seBXGs=mx12y[,c(4,8)],
                    seBYG=mx12y[,12],
                    RSID=mx12y[,1])
head(F.data)
sres<-strength_mvmr(r_input=F.data,gencov=0)
pres<-pleiotropy_mvmr(r_input=F.data,gencov=0)
res<-ivw_mvmr(r_input=F.data)


# Correlation

out_start=252
out_end= 252
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=2
exp_end=251
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

library(lme4)
for (j in exp_start:exp_end){
  exposure = colnames(Corrdata)[j]
  model <- lm(Corrdata$X1 ~ get(exposure) ,
              na.action = na.exclude,
              data=Corrdata)
  
  Vcov <- vcov(model, useScale = FALSE)
  beta <- coef(model)
  se <- sqrt(diag(Vcov))
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  
  exp_beta[number] = as.numeric(beta[2])
  exp_se[number] = as.numeric(se[2])
  exp_pvalue[number] = as.numeric(pval[2])
  exp_variable[number] = exposure
  number = number + 1
}


outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)


library(tidyverse)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
  )
corrx2 = rbind(outcome, exposure)
corrx2 = na.omit(corrx2)
corrx2<-tibble::rowid_to_column(corrx2, "ID")

# X1

out_start=253
out_end= 253
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=2
exp_end=251
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

library(lme4)
for (j in exp_start:exp_end){
  exposure = colnames(Condata)[j]
  model <- lm(Condata$X2 ~ get(exposure) ,
              na.action = na.exclude,
              data=Condata)
  
  Vcov <- vcov(model, useScale = FALSE)
  beta <- coef(model)
  se <- sqrt(diag(Vcov))
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  
  exp_beta[number] = as.numeric(beta[2])
  exp_se[number] = as.numeric(se[2])
  exp_pvalue[number] = as.numeric(pval[2])
  exp_variable[number] = exposure
  number = number + 1
}


outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)


library(tidyverse)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
  )
corrx1 = rbind(outcome, exposure)
corrx1 = na.omit(corrx1)
corrx1<-tibble::rowid_to_column(corrx1, "ID")
#Y


out_start=504
out_end= 504
out_nvar=out_end-out_start+1

out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)

# exposure
exp_start=254
exp_end=503
exp_nvar=exp_end-exp_start+1

exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)

number=1

library(lme4)
for (j in exp_start:exp_end){
  exposure = colnames(Condata)[j]
  model <- lm(Condata$Y ~ get(exposure) ,
              na.action = na.exclude,
              data=Condata)
  
  Vcov <- vcov(model, useScale = FALSE)
  beta <- coef(model)
  se <- sqrt(diag(Vcov))
  zval <- beta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  
  
  exp_beta[number] = as.numeric(beta[2])
  exp_se[number] = as.numeric(se[2])
  exp_pvalue[number] = as.numeric(pval[2])
  exp_variable[number] = exposure
  number = number + 1
}


outcome = data.frame(out_variable, out_beta, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta, exp_se, exp_pvalue)


library(tidyverse)
outcome = outcome %>% 
  rename(
    variable = out_variable,
    beta = out_beta,
    se = out_se,
    pvalue = out_pvalue
  )
exposure = exposure %>% 
  rename(
    variable = exp_variable,
    beta = exp_beta,
    se = exp_se,
    pvalue = exp_pvalue
  )
corry = rbind(outcome, exposure)
corry = na.omit(corry)
corry<-tibble::rowid_to_column(corry, "ID")

#x1+x2 

x12<-merge(corrx1,corrx2, by="variable")
names(x12)[3]<-"beta.x1"
names(x12)[4]<-"se.x1"
names(x12)[5]<-"pvalue.x1"
names(x12)[7]<-"beta.x2"
names(x12)[8]<-"se.x2"
names(x12)[9]<-"pvalue.x2"
names(x12)[9]<-"pvalue.x2"
names(x12)[2]<-"ID"
x12y<-merge(x12,corry, by="ID")
x12y<-subset(x12y, pvalue.x1<5*10^-8  | pvalue.x2 <5*10^-8)
#format
F.data<-format_mvmr(BXGs=x12y[,c(3,7)],
                    BYG=x12y[,11],
                    seBXGs=x12y[,c(4,8)],
                    seBYG=x12y[,12],
                    RSID=x12y[,1])
head(F.data)
sres<-strength_mvmr(r_input=F.data,gencov=0)
pres<-pleiotropy_mvmr(r_input=F.data,gencov=0)
res<-ivw_mvmr(r_input=F.data)



