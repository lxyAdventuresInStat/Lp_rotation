rm(list=ls());
library(doParallel)
library(lavaan)
library(gtools)
library(zoo)

#setwd('C:/Users/LIUX114/OneDrive - London School of Economics/0501')
source('cl510_oblique_rotation_function.R')
sim = 500
N = 3000
c_list<-c(seq(0,0.19,by=0.01),seq(0.2,2,by=0.1),3)
c_bic_list<-c(seq(0.02,0.4,by=0.02))
lambda_list<-c(seq(0.01,0.37,by=0.01),0.4,0.5,2)
t1=proc.time()




L<-matrix(c(0.5310173, 0.76007036, 0.0000000,
0.7442478,0.00000000,0.2158873
,0.0000000,1.86941046,1.4474219
,1.8164156,0.42428504,0.0000000
,0.4033639,0.00000000,1.6418926
,0.0000000,0.25111019,1.2941204
,1.8893505,0.53444134,0.0000000
,1.3215956,0.00000000,1.1060726
,0.0000000,0.02678067,1.0594392
,0.1235725,0.76477591,0.0000000
,0.4119491,0.00000000,0.0466624
,0.0000000,0.68069799,0.9544601
,1.3740457,0.96416023,0.0000000
,0.7682074,0.00000000,1.3854631
,0.0000000,0.98708261,0.9552392
,0.9953985,0.37243520,0.0000000
,1.4352370,0.00000000,0.8761942
,0.0000000,1.33693348,0.4895946),nrow=18,ncol=3,byrow<-TRUE)

B<-matrix(c(1,0.02079577,0.5024378,
            0,0.99978374,0.2635086,
            0,0,0.8234801),3,3,byrow=TRUE) #



P<-c(0.24, 0.32, 0.45, 0.65, 0.18, 0.64, 0.67, 0.51, 0.49, 0.06, 0.19, 0.16, 0.52, 0.33, 0.57,0.40,0.54,0.69)


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

out = foreach(i = 1:sim, .packages = c("lavaan","gtools","zoo"), .errorhandling = "pass") %dopar%{
  set.seed(i)
  
  ## generate data
  SLP0=Gen.SLP(Sig,N,n_col)
  S=SLP0$S
  L_vari=SLP0$L_vari
  Psi_cfa=SLP0$Psi_cfa
  
  
  res=list()
  L_vari=permfilp(L_vari,L)
  #res$vari$mse= mean((L_vari-L)^2)
  
  ## irls p=1
  ans_est1<-irls(1,L_vari)
  res$irls1$L=L_est=permfilp(ans_est1$L,L)
  res$irls1$t=ans_est1$t
  res$irls1$it=ans_est1$it
  res$irls1$mse= mean((L_est-L)^2)
  #hard-thresholding
  ans_hard=hard(c_list,L_est,L_old)
  res$irls1$TPR =ans_hard$TPR
  res$irls1$TNR =ans_hard$TNR
  #res$irls1$AUC =ans_hard$AUC
  #bic-accuracy
  bic=sapply(c_bic_list,refit.bic,L_irls=L_est,S=S,N=N)
  res$irls1$c= c = c_bic_list[which.min(bic)[1]]
  res$irls1$L_bic=L_bic=refit.L(c,L_est,S,N)
  res$irls1$L_bic.res=hard.each(c,L_est,L_old)

  
  #irls p=0.5
  ans_est<-irls(0.5,L_vari,T=ans_est1$T)
  res$irls0.5$L=L_est=permfilp(ans_est$L,L)
  res$irls0.5$t=ans_est$t
  res$irls0.5$it=ans_est$it
  res$irls0.5$mse= mean((L_est-L)^2)
  #hard-thresholding
  ans_hard=hard(c_list,L_est,L_old)
  res$irls0.5$TPR =ans_hard$TPR
  res$irls0.5$TNR =ans_hard$TNR
  #res$irls0.5$AUC =ans_hard$AUC
  #bic-accuracy
  bic=sapply(c_bic_list,refit.bic,L_irls=L_est,S=S,N=N)
  res$irls0.5$c= c = c_bic_list[which.min(bic)[1]]
  res$irls0.5$L_bic=L_bic=refit.L(c,L_est,S,N)
  res$irls0.5$L_bic.res=hard.each(c,L_est,L_old)
  

  
  return(res)
}
t=proc.time()[3]-t1[3]
save(out,L,B,P,sim,N,seed, t,file = paste0("0510Counter.15S_Fa_CI_p",nrow(L),"k",n_col,"nsim", sim,"nsample", N, ".RData"))

stopCluster(cl)