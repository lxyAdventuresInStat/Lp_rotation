rm(list=ls());
library(doParallel)
library(lavaan)
library(gtools)
library(zoo)
library(GPArotation)

#setwd('C:/Users/LIUX114/OneDrive - London School of Economics/0501')
source('cl510_oblique_rotation_function.R')
sim = 500
c_list<-c(seq(0,0.19,by=0.01),seq(0.2,2,by=0.1),3)
c_bic_list<-c(seq(0.02,0.4,by=0.02))
lambda_list<-c(seq(0.01,0.37,by=0.01),0.4,0.5,2)
t1=proc.time()

####30*5 SMALL--------------------------------------------------------------------------

  L<-matrix(c(0.71,0.00,0.00,0.00,0.00,
0.00,0.75,0.00,0.00,0.00,
0.00,0.00,0.83,0.00,0.00,
0.00,0.00,0.00,0.96,0.00,
0.00,0.00,0.00,0.00,0.68,
0.96,0.00,0.00,0.00,0.00,
0.00,0.98,0.00,0.00,0.00,
0.00,0.00,0.86,0.00,0.00,
0.00,0.00,0.00,0.85,0.00,
0.00,0.00,0.00,0.00,0.62,
0.68,0.00,0.00,0.00,0.00,
0.00,0.67,0.00,0.00,0.00,
0.00,0.00,0.87,0.00,0.00,
0.00,0.00,0.00,0.75,0.00,
0.00,0.00,0.00,0.00,0.91,
0.80,0.34,0.00,0.00,0.00,
0.00,0.89,0.38,0.00,0.00,
0.00,0.00,1.00,0.35,0.00,
0.00,0.00,0.00,0.75,0.26,
0.45,0.00,0.00,0.00,0.91,
0.97,0.00,0.40,0.00,0.00,
0.00,0.68,0.00,0.44,0.00,
0.00,0.00,0.86,0.00,0.23,
0.42,0.00,0.00,0.65,0.00,
0.00,0.32,0.00,0.00,0.71,
0.75,0.45,0.39,0.00,0.00,
0.00,0.61,0.43,0.37,0.00,
0.00,0.00,0.75,0.36,0.44,
0.34,0.00,0.00,0.95,0.21,
0.42,0.41,0.00,0.00,0.74),nrow<-30,ncol=5,byrow<-TRUE);

B<-matrix(c(1,0.08507226,0.429328611,0.1480261,0.24930445,
0,0.99637478,0.005627324,0.1368557,0.28840407,
0,0.00000000,0.903130819,0.4776413,0.01373041,
0,0.00000000,0.000000000,0.8551126,0.12494345,
0,0.00000000,0.000000000,0.0000000,0.91589901),5,5,byrow=TRUE)
#P<-c(0.23547415, 0.31635983, 0.45289140, 0.64616447, 0.18372218, 0.64100599, 0.66509501, 0.50729808, 0.48803634, 0.05995265, 0.18728802, 0.16259217,0.52296535, 0.32505280, 0.57088995, 0.40393009, 0.54093874, 0.68909202, 0.32210899, 0.57517706, 0.65995498, 0.19238947, 0.50178918, 0.11827633,0.23682605, 0.32650422, 0.01330147, 0.32381241, 0.62577309, 0.29293003)
P<-c(0.24,0.32,0.45,0.65,0.18,0.64,0.67,0.51,0.49,0.06,0.19,0.16,0.52,0.33,0.57,0.40,0.54,0.69,0.32,0.58,0.66,0.19,0.50,0.12,0.24,0.33,0.01,0.32,0.63,0.29)
#########################################################

n_row<- nrow(L)
n_col<- ncol(L)
Sig <- L%*%t(B)%*%B%*%t(L)+diag(exp(P))
(Sig)
(exp(P))
(t(B)%*%B)
L_old=L!=0
seed=0

############ Replicates
ncores=detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)

############ Replicates
ncores=detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
out = foreach(i = 1:sim, .packages = c("lavaan","gtools","zoo","GPArotation"), .errorhandling = "pass") %dopar%{
  set.seed(i)
  N = 400
  ## generate data
  SLP0=Gen.SLP(Sig,N,n_col)
  S=SLP0$S
  L_vari=SLP0$L_vari
  Psi_cfa=SLP0$Psi_cfa
  
  
  res=list()
  res$oblimin$mse= mean((oblimin(L_vari)[[1]]-L)^2)
  res$quartimin$mse= mean((quartimin(L_vari)[[1]]-L)^2)
  res$simplimax$mse= mean((simplimax(L_vari)[[1]]-L)^2)
  res$geominQ$mse= mean((geominQ(L_vari)[[1]]-L)^2)
  res$promax$mse=mean((promax(L_vari)[[1]]-L)^2)
  
  set.seed(i)
  N = 800
  ## generate data
  SLP0=Gen.SLP(Sig,N,n_col)
  S=SLP0$S
  L_vari=SLP0$L_vari
  Psi_cfa=SLP0$Psi_cfa
  
  
  
  res$oblimin$mse2= mean((oblimin(L_vari)[[1]]-L)^2)
  res$quartimin$mse2= mean((quartimin(L_vari)[[1]]-L)^2)
  res$simplimax$mse2= mean((simplimax(L_vari)[[1]]-L)^2)
  res$geominQ$mse2= mean((geominQ(L_vari)[[1]]-L)^2)
  res$promax$mse2=mean((promax(L_vari)[[1]]-L)^2)
  
  set.seed(i)
  N = 1600
  ## generate data
  SLP0=Gen.SLP(Sig,N,n_col)
  S=SLP0$S
  L_vari=SLP0$L_vari
  Psi_cfa=SLP0$Psi_cfa
  res$oblimin$mse3= mean((oblimin(L_vari)[[1]]-L)^2)
  res$quartimin$mse3= mean((quartimin(L_vari)[[1]]-L)^2)
  res$simplimax$mse3= mean((simplimax(L_vari)[[1]]-L)^2)
  res$geominQ$mse3= mean((geominQ(L_vari)[[1]]-L)^2)
  res$promax$mse3=mean((promax(L_vari)[[1]]-L)^2)
  return(res)
}
t=proc.time()[3]-t1[3]
#save(out,L,B,P,sim,N,seed, t,file = paste0("0510A.15L_Fa_CI_p",nrow(L),"k",n_col,"nsim", sim,"nsample", N, ".RData"))
stopCluster(cl)
line<-function(out,j,k){
  mse=mean(sapply(1:sim,function(i,out){out[[i]][[j]][[k]]},out=out))
}
res.obl=matrix(0,5,3)
rownames(res.obl) = c("oblimin","quartimin","simplimax", "geominQ", "promax")
colnames(res.obl) = c("N=400","N=800","N=1600") 
for ( j in 1:5){
  for (k in 1:3){
    res.obl[j,k]=line(out,j,k)
  }
}

write.csv(res.obl,file = paste0("res.obl30",".csv"))