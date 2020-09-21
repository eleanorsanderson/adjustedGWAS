#function to generate data under different scenarios
#input (confounder, mediator, correlated) - relates to the position of the adjustment variable X2 relative to X1 and the outcome. 
#requires functions simulateGP and dplyr

datagen <- function(model,beta1,beta2){

n=100000
nsnps = 500
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
