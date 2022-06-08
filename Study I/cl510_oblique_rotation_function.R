library(lavaan)
#------------------------------------- IRGP Algorithm -------------------------------------- 
Q<-function(w,L){return(sum(w*L*L))};
dQ<-function(w,L){return(2*w*L)}
rho<-function(x){return(x%*%diag(diag(t(x)%*%x)^(-1/2)))};
Q1<- function(L,p){return(sum(abs(L)^p))}

# Gradient Projection Algorithm
gpa <-function(A,TR,W,iter,stop_ep,alpha0=1){
  # Q is the weighted least square function
  alpha<-alpha0;
  Ti<-solve(TR);
  L<-A%*%t(Ti);
  ft<-Q(W,L);
  G<- -t(t(L)%*%dQ(W,L)%*%Ti);
  
  for (j in 1:iter){   
    alpha<-2*alpha;
    # find a proper step size
    for( i in 1:20 ){
      X<-TR-alpha*G;
      T_new<-rho(X);
      L<-A%*%solve(t(T_new));
      if (Q(W,L)<ft){
        ft<-Q(W,L);
        break;
      }else{
        alpha<-alpha*0.5;
      }
    }
    
    # update gradient and T
    TR<-   T_new;
    G<-  -t(t(L)%*%dQ(W,L)%*%solve(TR));
  }
  T_final<-TR;
  L_final<-A%*%solve(t(T_final));
  return(list(T_final,L_final,ft,j))
}



# iterative re-weighted least square algorithm
irls<-function(p,A,irls_ep0=1e-3,
               stop_ep=1e-5,
               irls_iter=10000,gpa_iter=5,
               T=diag(rep(1,ncol(A))),W=(abs(A)+mean(abs(A))*irls_ep0)^(p-2), 
               alpha0=1,power=1/2){
 
  j_sum<-0
  start_time<-proc.time()
  for (it in 1:irls_iter){

    # solve by gpa
    gpa_result <- gpa(A,T,W,gpa_iter,stop_ep,alpha0);
    T_new <-gpa_result[[1]]
    L <-gpa_result[[2]]
    ft <-gpa_result[[3]]
    j <-gpa_result[[4]] 
    j_sum<-j_sum+j
    
    #repeatedly update weights in the objective function
    W <- (abs(L)+mean(abs(L))*irls_ep0)^(p-2)
    if (max(abs(L-A%*%solve(t(T))))< stop_ep){
      break
    }
    T<-T_new;
    t<-proc.time()[3]-start_time[3]
}
  return(list(T=T,L=L,
              ft=ft,
              it=it,gpa_it=j_sum/it,
              t=t,obj=Q1(L,p)))
}

#------------------------------------- other functions -------------------------------------- 
rho <- function(x){return(x%*%diag(diag(t(x)%*%x)^(-1/2)))};# normalization function

hard.each<-function(c,L_irls,L_old){
  L_new=abs(L_irls)>c
  TPR= sum(L_new[L_old==1])/(sum(L_old))
  TNR= 1-sum(L_new[L_old==0])/(length(L_old)-sum(L_old))
  TR=sum(L_new==L_old)/(nrow(L_old)*ncol(L_old))
  TRall=all(L_new==L_old)
  return(list(TPR=TPR,TNR=TNR,TR=TR,TRall=TRall))
}

hard<-function(c_list,L_irls,L_old){

  ans_hard.each=sapply(c_list,hard.each,L_irls=L_irls,L_old=L_old)
  TPR=unlist(ans_hard.each[1,])
  TNR=unlist(ans_hard.each[2,])
  id <- 1:(length(TNR)+1)
  AUC <- sum(diff((1-c(1,TPR))[id])*rollmean(c(0,TNR)[id],2))
  return(list(TPR=TPR,TNR=TNR,AUC=AUC))
}

# Permutation and Sign-flip of the solution
permfilp <- function(L_irls,L_list){
  # Sign-filp
  n_col=ncol(L_irls)
  for( i in 1:n_col){
    if(sum(L_irls[,i])<0){
      L_irls[,i]<- -L_irls[,i]}}
  #Column-swapping
  #library(gtools)
  P<-permutations(n_col, n_col)
  minf<-norm(L_list-L_irls,'f')
  L_new0<-L_irls
  for (i in 1:nrow(P)){
    if (norm(L_list-L_irls[,P[i,]],'f')< minf){
      minf<-norm(L_list-L_irls[,P[i,]],'f')
      L_new0<-L_irls[,P[i,]]
    }
  }
  return(L_new0)
}

# Convert matrix to syntax in CFA model
modeltext<-function(L_tmp){
  n_row=nrow(L_tmp)
  n_col=ncol(L_tmp)
  main=''
  for (i in 1:n_col){
    t=0
    for (j in 1:n_row)
    {if(L_tmp[j,i]!=0){
      if(t==0){
        tmp=paste0('y',i,'=~','NA*','x',j)
      }else{
        tmp=paste0('+','x',j)
      }
      t=1 
      main=paste0(main, tmp)
      
    } 
    } 
    main=paste0(main,'\n')
  }
  main=paste0(main, '# factors with unit variance\n')
  for ( i in 1:n_col){
    if (!all(!L_tmp[,i])){
      main = paste0(main,'y',i,'~~','1*y',i,'\n')}
  }
  main=paste0(main, '# constrains\n')
  for ( j in 1:n_row){
    main = paste0(main,'x',j,'~~','c',j,'*','x',j,'\n')
    main = paste0(main,'c',j,'>0','\n')}
  #cat(main)
  return(main)
}
#Generate data for simulation
Gen.SLP<-function(Sig,N,n_col){
  n_row<-nrow(Sig)
  X <-t(chol(Sig)) %*% matrix(rnorm(N*n_row),n_row,N)
  S <- cov(t(X))
  var_nam=c()
  for( j in 1:n_row){
    var_nam<-c(var_nam,paste0('x',j))}
  colnames(S)[1:n_row] <- var_nam
  rownames(S)[1:n_row] <- var_nam
  
  L_tmp<-matrix(1,n_row,n_col)
  L_tmp[upper.tri(L_tmp)]<-0
  main<-modeltext(L_tmp)
  fit <-cfa(main,
            sample.cov=S,
            sample.nobs=N,
            orthogonal=T,
            se="none")
  ncoef.L<-(n_row*n_col-(n_col-1)*n_col/2)
  L_cfa<-matrix(0,n_row,n_col)
  L_cfa[lower.tri(L_cfa,diag=TRUE)]<-coef(fit)[1:ncoef.L]
  Psi_cfa<-coef(fit)[ncoef.L+1:n_row]
  L_vari<-L_cfa %*% varimax(L_cfa)$rotmat
  return(list(S=S,
              L_vari=L_vari,
              Psi_cfa=Psi_cfa))
}
#refit cfa model by q-matrix and calculate extended bic
refit.bic<-function(c,L_irls,S,N){
  L_zeros<-abs(L_irls)>c
  main<-modeltext(L_zeros)
  fit <-cfa(main,
            sample.cov=S,
            sample.nobs=N,
            se="none")
  n_row=nrow(L_zeros)
  n_col=ncol(L_zeros)
  bic_ex= log(choose(n_row*n_col,sum(L_zeros)))
  return(bic=BIC(fit)+bic_ex)
}
#refit cfa model by q-matrix and calculate loading matrix
refit.L<-function(c,L_irls,S,N){
  L_zeros<-abs(L_irls)>c
  main<-modeltext(L_zeros)
  fit <-cfa(main,
            sample.cov=S,
            sample.nobs=N,
            se="none")
  L_cfa<-matrix(0,n_row,n_col)
  L_cfa[which(L_zeros!=0)]<-coef(fit)[1:sum(L_zeros)]
  return(L_cfa)
}
#refit cfa model by q-matrix and calculate loading matrix and extended bic
refit.Lbic<-function(c,L_irls,S,N){
  L_zeros<-abs(L_irls)>c
  main<-modeltext(L_zeros)
  fit <-cfa(main,
            sample.cov=S,
            sample.nobs=N,
            se="none")
  L_cfa<-matrix(0,n_row,n_col)
  L_cfa[which(L_zeros!=0)]<-coef(fit)[1:sum(L_zeros)]
  n_row=nrow(L_zeros)
  n_col=ncol(L_zeros)
  bic_ex= log(choose(n_row*n_col,sum(L_zeros)))
  return(list(bic=BIC(fit)+bic_ex,L_bic=L_cfa))
}
#calculate each line of CI
CI_line<-function(k,L_zeros,S,N){
  L_tmp<-L_zeros
  L_tmp[k,]<-1
  main<-modeltext(L_tmp)
  fit <-cfa(main,sample.cov=S,sample.nobs=N) 
  L_uppertmp<-matrix(0,n_row,n_col)
  L_uppertmp[which(L_tmp!=0)]<-parameterEstimates(fit)$ci.upper[1:sum(L_tmp)]
  L_lowertmp<-matrix(0,n_row,n_col)
  L_lowertmp[which(L_tmp!=0)]<-parameterEstimates(fit)$ci.lower[1:sum(L_tmp)]
  return( c(L_lowertmp[k,],L_uppertmp[k,]))
}
#calculate CI and set na to inf
CI<-function(L_bic,S,N,L){
  n_row=nrow(L)
  n_col=ncol(L)
  bounds=sapply(1:n_row,CI_line,L_zeros=(L_bic!=0),S=S,N=N)
  L_lower=t(bounds[1:n_col,])
  L_upper=t(bounds[n_col+1:n_col,])
  L_lower[which(is.na(L_lower))]<--Inf
  L_upper[which(is.na(L_upper))]<-Inf
  L_class<-(L_lower<=L)&(L<=L_upper)
  accuracy<-sum(L_class)/length(L_lower)
  return(list(L_class=L_class,
              L_upper=L_upper, L_lower= L_lower,
              ci.accuracy=accuracy,na.flag=!all(!is.infinite(L_lower)) ))
}
#calculate lasso path solution-not used
lasso.path<-function(lambda_list,L_vari,B,Psi_cfa,S,N,L){
  nlambda<-length(lambda_list)
  TPR=TNR=rep(0,nlambda)
  L0_list=list()
  t=it=0
  L0=L_vari
  B0=B
  Psi0=Psi_cfa
  bic.min=Inf
  for (j in 1:nlambda){
    ans_est<-prox_grad(L0,B0,Psi0,S,lambda_list[j],1)
    L0_list[[j]]=L0=permfilp(ans_est$L,L)
    Psi0=ans_est$Psi
    t=ans_est$t+t
    it=ans_est$it+it
    
    #soft-thresholding
    ans_hard=hard.each(0,L0,L_old)
    TPR[j] =ans_hard$TPR
    TNR[j] =ans_hard$TNR
    
    if(lambda_list[j]){
      ans_bic=refit.Lbic(0,L0,S,N)
      bic=ans_bic$bic
      if(bic<bic.min){
        L_bic=ans_bic$L_bic
        #print(L_bic)
        l=lambda_list[j]
        bic.min=bic
      }
    }
  }
  L_bic.res=hard.each(0,L_bic,L_old)
  id <- 1:(length(TNR)+1)
  AUC <- sum(diff((1-c(1,TPR))[id])*rollmean(c(0,TNR)[id],2))
  return(list(L=L0_list,t=t,it=it, 
              TPR=TPR,TNR=TNR,AUC=AUC,
              c=l,L_bic=L_bic,L_bic.res=L_bic.res))
}
#------------------------------------- lasso functions -------------------------------------- 
## the main function f = g + h
f <- function(Sigma, S, Lambda,lambda){
  log(det(Sigma))+sum(diag(S %*% solve(Sigma))) + lambda*sum(abs(Lambda))
}

## smooth function g
g <- function(Sigma,S){
  log(det(2*pi*Sigma))+sum(diag(S %*% solve(Sigma)))
}

lkhd1 <- function(Lambda,B,Psi,S){
  Sigma <-Lambda%*%t(B)%*%B%*%t(Lambda)+diag(exp(Psi))
  log(det(2*pi*Sigma))+sum(diag(S %*% solve(Sigma)))
}

lkhd2 <- function(Lambda,T_i,Psi_cfa,S){
  Sigma <-Lambda%*%t(B)%*%B%*%t(Lambda)+diag(Psi_cfa)
  log(det(2*pi*Sigma))+sum(diag(S %*% solve(Sigma)))
}

## subgradient of smooth function
Subg_Lambda <- function(Q,Lambda,Phi){
  2*Q %*% Lambda %*% Phi
}
Subg_B <- function(Q,Lambda,B){
  grad_B <-2*B %*% t(Lambda) %*% Q %*% Lambda
  grad_B[lower.tri(grad_B)] <- 0
  return(grad_B)
}

Subg_Psi <- function(Q,Psi){
  diag(Q)*exp(Psi)
}

## proximal function of none smooth function h
prox_L1 <- function(Lambda, lambda){
  sign(Lambda) * pmax(abs(Lambda) - lambda, 0)
}



# proximal gradient descend 
prox_grad<-function(Lambda0,B0,Psi0,S,lambda,if_fixB=0,
                    maxiter=10000,stop_ep=10^(-6),t0=1){
  
  Lambda<-Lambda0
  B<-B0
  Psi<-Psi0
  Sigma <-Lambda%*%t(B)%*%B%*%t(Lambda)+diag(exp(Psi))
  
  # initilization
  t <- t0                  # step size
  beta0 <- 0.5             # line search factor for t paramter

  start_time<-proc.time()
  for(i in 1:maxiter){
    
    t<-2*t
    Sigma_inv <-solve(Sigma)
    Q=Sigma_inv-Sigma_inv%*%S%*%Sigma_inv
    Phi=t(B)%*%B
    
    # gradient
    grad_Lambda <- Subg_Lambda(Q,Lambda,Phi)
    grad_B <- Subg_B(Q,Lambda,B)
    grad_Psi <- Subg_Psi(Q,Psi)
    
    # proximal gradient descend stage and line search
    for( i_iter in 1:20 ){
      
      Lambda_new <- prox_L1(Lambda - t*grad_Lambda, t*lambda)
      if (if_fixB==1){
        B_new<-B0  
      }else{
        B_new <-B- t*grad_B
        B_new <-rho(B_new)
      }                       
      # here we do not want Psi goes to negative infinity and Sigma degenerated
      Psi_new<- Psi- t*grad_Psi
      Psi_new[Psi_new > 3] <- 3
      Psi_new[Psi_new < -3]<- -3
      Sigma_new <-Lambda_new%*%t(B_new)%*%B_new%*%t(Lambda_new)+diag(exp(Psi_new)) 
      # line search the step size
      if(f(Sigma_new, S, Lambda_new,lambda) < f(Sigma, S, Lambda,lambda)) {
        break}
      t <- beta0*t
    }
  
    # all model parameters needs to converge
     if(i > 1 && (max(abs(Lambda- Lambda_new)) < stop_ep)
        && (max(abs(B-B_new)) < stop_ep)&& (max(abs(Psi- Psi_new)) < stop_ep)) {
       break }

    
    # update
    Lambda <- Lambda_new
    B <- B_new
    Psi <- Psi_new
    Sigma <-Sigma_new
  }
  time=proc.time()[3]-start_time[3]

  return(list(L=Lambda,B=B,Psi=Psi,it=i,t=time))
}





