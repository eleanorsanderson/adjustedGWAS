library(stringr)
library(simulateGP)
library(MVMR)
library(dplyr)
set.seed(5)
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


reps = 1000

results = data.frame()
snps = 250
for(j in 1:reps){

dat <- datagen("correlated",250,0,0.5) 

MR_dat = data.frame()

G <- dat[,1:snps]
G2 <- dat[,(snps+3):((2*snps)+2)]

for(i in 1:snps){
  a <- summary(lm(dat$X2~G[,i]))
  MR_dat[i,"bmi_b"] <- a$coefficient[2,1]
  MR_dat[i,"bmi_se"] <- a$coefficient[2,2]
  MR_dat[i,"bmi_p"] <- a$coefficient[2,4]
  b<-summary(lm(dat$X1~G[,i]+dat$X2))
  MR_dat[i,"sbpaj_b"] <- b$coefficient[2,1]
  MR_dat[i,"sbpaj_se"] <- b$coefficient[2,2]
  MR_dat[i,"sbpaj_p"] <- b$coefficient[2,4]
  c<-summary(lm(dat$Y~G2[,i]))
  MR_dat[i,"out_b"] <- c$coefficient[2,1]
  MR_dat[i,"out_se"] <-c$coefficient[2,2]
  MR_dat[i,"out_p"] <- c$coefficient[2,4]
  d<-summary(lm(dat$X1~G[,i]))
  MR_dat[i,"sbp_b"] <- d$coefficient[2,1]
  MR_dat[i,"sbp_se"] <-d$coefficient[2,2]
  MR_dat[i,"sbp_p"] <- d$coefficient[2,4]
  
}

#select sps associated with each exposure

MR_dat.bmi <- subset(MR_dat, MR_dat$bmi_p < 5e-08)
MR_dat.sbp <- subset(MR_dat, MR_dat$sbp_p < 5e-08)
MR_dat.sbpaj <- subset(MR_dat, MR_dat$sbpaj_p < 5e-08)

MR_dat[,"min_p"] <- apply(MR_dat[,c("bmi_p", "sbp_p")],1,min)
MR_dat.mvmr <- subset(MR_dat, MR_dat$min_p < 5e-08)

MR_dat[,"min_paj"] <- apply(MR_dat[,c("bmi_p", "sbpaj_p")],1,min)
MR_dat.mvmraj <- subset(MR_dat, MR_dat$min_paj < 5e-08)

bmi_out <- summary(lm(out_b~-1+bmi_b,weights=1/(out_se^2), data = MR_dat.bmi))
sbp_out <- summary(lm(out_b~-1+sbp_b,weights=1/(out_se^2), data = MR_dat.sbp))
sbpaj_out <- summary(lm(out_b~-1+sbpaj_b,weights=1/(out_se^2), data = MR_dat.sbpaj))
mvmr_out <- summary(lm(out_b~-1+bmi_b+sbp_b,weights=1/(out_se^2), data = MR_dat.mvmr))
mvmraj_out <- summary(lm(out_b~-1+bmi_b+sbpaj_b,weights=1/(out_se^2), data = MR_dat.mvmraj))

results[j,"bmi_effect"] <- bmi_out$coefficients["bmi_b","Estimate"]
results[j,"bmi_se"] <- bmi_out$coefficients["bmi_b","Std. Error"]

results[j,"sbp_effect"] <- sbp_out$coefficients["sbp_b","Estimate"]
results[j,"sbp_se"] <- sbp_out$coefficients["sbp_b","Std. Error"]

results[j,"sbpaj_effect"] <- sbpaj_out$coefficients["sbpaj_b","Estimate"]
results[j,"sbpaj_se"] <- sbpaj_out$coefficients["sbpaj_b","Std. Error"]

results[j,"mvmr.bmi_effect"] <- mvmr_out$coefficients["bmi_b","Estimate"]
results[j,"mvmr.bmi_se"] <- mvmr_out$coefficients["bmi_b","Std. Error"]
results[j,"mvmr.sbp_effect"] <-mvmr_out$coefficients["sbp_b","Estimate"]
results[j,"mvmr.sbp_se"] <- mvmr_out$coefficients["sbp_b","Std. Error"]
results[j,"mvmr.nsnps"] <- length(MR_dat.mvmr$bmi_b)

results[j,"mvmraj.bmi_effect"] <- mvmraj_out$coefficients["bmi_b","Estimate"]
results[j,"mvmraj.bmi_se"] <- mvmraj_out$coefficients["bmi_b","Std. Error"]
results[j,"mvmraj.sbpaj_effect"] <-mvmraj_out$coefficients["sbpaj_b","Estimate"]
results[j,"mvmraj.sbpaj_se"] <- mvmraj_out$coefficients["sbpaj_b","Std. Error"]
results[j,"mvmraj.nsnps"] <- length(MR_dat.mvmraj$bmi_b)

}

results
write.csv(results, file = "results_med_nullexp_5.csv")
