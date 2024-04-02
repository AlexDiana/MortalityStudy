
#set working directory
#setwd("G:/IDEA-LR/")

##load packages
library(rjags)
#install.packages("insight")
library(insight)
library(demography)

###############
##Some functions
##############

#percentile

percentile.fn<-function(x){
  quantile(x,c(0.05,0.5,0.95))
}

#DIC

DIC_jags2_fn<-function(mcmc_object,deaths_combined_vectorised,expo_combined_vectorised,A,T,p){
  n<-dim(mcmc_object)[1] #number of posterior samples
  rates_mat<-matrix(0,nrow=n,ncol=A*T)
  
  deviance_vec<-vector(length=n)
  deviance_mean<-NULL
  
  rates_mat<-mcmc_object[,startsWith(colnames(mcmc_object),"q")]
  rates_mean<-apply(rates_mat,2,mean)
  
  for (s in 1:n){
    deviance_vec[s]<--2*sum(dbinom(x=deaths_combined_vectorised,size=expo_combined_vectorised,prob=rates_mat[s,],log=TRUE))
  }
  
  deviance_mean<--2*sum(dbinom(x=deaths_combined_vectorised,size=expo_combined_vectorised,prob=rates_mean,log=TRUE))
  DIC<-2*mean(deviance_vec)-deviance_mean
  DIC
}

#auto model-selection JAGS model fitting

auto_fit_DIC_fn<-function(deaths_combined_array,expo_combined_array,n_iter=10000){
  
  p<-dim(deaths_combined_array)[1]
  A<-dim(deaths_combined_array)[2]
  T<-dim(deaths_combined_array)[3]
  
  deaths_combined_vectorised=as.vector(deaths_combined_array)
  expo_combined_vectorised=as.vector(expo_combined_array)
  
  ages<-as.numeric(dimnames(deaths_combined_array)[[2]])
  years<-as.numeric(dimnames(deaths_combined_array)[[3]])
  
  M=23
  DIC_table<-matrix(0,nrow=1,ncol=M)
  colnames(DIC_table)<-c("LC_M1A","LC_M1U","LC_M1M","LC_M2A1","LC_M2A2","LC_M2Y1","LC_M2Y2","LC_MLiLee","LC_MLiLee_modified","LC_all_separately","CBD_M3","CBD_M3_sep","CBD_M5","CBD_M6","CBD_M6_sep","CBD_M7","CBD_M7_sep","CBD_M8","CBD_M8_sep","APCI","APCI_sep","RH","RH_sep")
  iter=1
  
  #Model 1A
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_kappa<-rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_kappa<-solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1),c_p_rest=rep(0,p-1)))
  vars<-c("q","alpha","beta","kappa","c_p")
  logit_LC_M1A_jags<-jags.model("logit_LC_M1A.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M1A_jags<-coda.samples(logit_LC_M1A_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_M1A"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M1A_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: LC_M1A"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  #Model 1U
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_kappa<-rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_kappa<-solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1)))
  vars<-c("q","alpha","beta","kappa")
  logit_LC_M1U_jags<-jags.model("logit_LC_M1U.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M1U_jags<-coda.samples(logit_LC_M1U_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_M1U"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M1U_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: LC_M1U"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  #Model 1M
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_kappa<-rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_kappa<-solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1),c_p_rest=rep(0,p-1)))
  vars<-c("q","alpha","beta","kappa","c_p")
  logit_LC_M1M_jags<-jags.model("logit_LC_M1M.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M1M_jags<-coda.samples(logit_LC_M1M_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_M1M"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M1M_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: LC_M1M"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  #Model 2A1
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_kappa<-rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_kappa<-solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1),c_p_rest=rep(0,p-1)))
  vars<-c("q","alpha","beta","kappa","c_p")
  logit_LC_M2A1_jags<-jags.model("logit_LC_M2A1.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M2A1_jags<-coda.samples(logit_LC_M2A1_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_M2A1"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M2A1_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: LC_M2A1"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  #Model 2A2
  
  prior_mean_beta<-rep(1/(p*A),p*A-1)
  sigma2_beta<-0.001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,p*A-1))-1/(p*A)*(matrix(1,nrow=(p*A-1),ncol=(p*A-1)))))
  prior_mean_kappa<-rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_kappa<-solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/(p*A),p*A-1),kappa_rest=rep(0,T-1)))
  vars<-c("q","alpha","beta","kappa")
  logit_LC_M2A2_jags<-jags.model("logit_LC_M2A2.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M2A2_jags<-coda.samples(logit_LC_M2A2_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_M2A2"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M2A2_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: LC_M2A2"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  #Model 2Y1
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_kappa<-rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_kappa<-solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1),c_p_rest=rep(0,p-1)))
  vars<-c("q","alpha","beta","kappa","c_p")
  logit_LC_M2Y1_jags<-jags.model("logit_LC_M2Y1.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M2Y1_jags<-coda.samples(logit_LC_M2Y1_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_M2Y1"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M2Y1_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: LC_M2Y1"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  #Model 2Y2
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_kappa<-rep(0,p*T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,p*T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_kappa<-solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,p*T-1)))
  vars<-c("q","alpha","beta","kappa")
  logit_LC_M2Y2_jags<-jags.model("logit_LC_M2Y2.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M2Y2_jags<-coda.samples(logit_LC_M2Y2_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_M2Y2"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M2Y2_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: LC_M2Y2"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  #Model LiLee
  
  prior_mean_Beta=prior_mean_beta=rep(1/A,A-1)
  sigma2_beta<-0.001
  prior_prec_Beta=prior_prec_beta=solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_Kappa=prior_mean_kappa=rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_Kappa=prior_prec_kappa=solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_Beta=prior_mean_Beta,prior_prec_Beta=prior_prec_Beta,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_Kappa=prior_mean_Kappa,prior_prec_Kappa=prior_prec_Kappa,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),Beta_rest=rep(1/A,A-1),beta_rest_mat=matrix(1/A,nrow=p-1,ncol=A-1),Kappa_rest=rep(0,T-1),kappa_rest_mat=matrix(0,nrow=p-1,ncol=T-1)))
  vars<-c("q","alpha","beta","kappa","Beta","Kappa")
  logit_LC_MLiLee_jags<-jags.model("logit_LC_MLiLee.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_MLiLee_jags<-coda.samples(logit_LC_MLiLee_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_MLiLee"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_MLiLee_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: MLiLee"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##Model LiLee modified 
  
  prior_mean_Beta=prior_mean_beta=rep(1/A,A-1)
  sigma2_beta<-0.001
  prior_prec_Beta=prior_prec_beta=solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_Kappa=prior_mean_kappa=rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_Kappa=prior_prec_kappa=solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_Beta=prior_mean_Beta,prior_prec_Beta=prior_prec_Beta,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_Kappa=prior_mean_Kappa,prior_prec_Kappa=prior_prec_Kappa,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=rep(0,A),Beta_rest=rep(1/A,A-1),beta_rest_mat=matrix(0,nrow=p-1,ncol=A-1),Kappa_rest=rep(0,T-1),kappa_rest_mat=matrix(0,nrow=p-1,ncol=T-1)))
  vars<-c("q","alpha","beta","kappa","Beta","Kappa")
  logit_LC_MLiLee_modified_jags<-jags.model("logit_LC_MLiLee_modified.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_MLiLee_modified_jags<-coda.samples(logit_LC_MLiLee_modified_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_MLiLee_modified"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_MLiLee_modified_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: MLiLee modified"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##all LC separately
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  mean_kappa=rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_kappa=solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),beta_rest_mat=matrix(1/A,nrow=p,ncol=(A-1)),kappa_rest_mat=matrix(0,nrow=p,ncol=(T-1))))
  vars<-c("q","alpha","beta","kappa")
  logit_LC_all_separately_jags<-jags.model("logit_LC_all_separately.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_all_separately_jags<-coda.samples(logit_LC_all_separately_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"LC_all_separately"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_all_separately_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: all LC_sep"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##CBD M3
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  cohorts_count_mat<-matrix(NA,nrow=A,ncol=T)
  for (i in 1:A){
    for(j in 1:T){
      cohorts_count_mat[i,j]<-j-i+A
    }
  }
  cohorts_count<-as.vector(table(cohorts_count_mat))
  prior_mean_kappa=rep(0,T-1)
  sigma2_kappa<-1000
  matrix_kappa_A<-diag(rep(1,T));matrix_kappa_A[1,]<-1
  matrix_kappa_B<-matrix_kappa_A%*%t(matrix_kappa_A)
  prior_prec_kappa=solve((matrix_kappa_B[-1,-1]-1/matrix_kappa_B[1,1]*matrix_kappa_B[-1,1]%*%t(matrix_kappa_B[1,-1]))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-2)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C-1))
  matrix_cohort_A<-rbind(cohorts_count[-1],matrix_cohort_A[(1:(C-2)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[-1,-1]-1/(matrix_cohort_B[1,1])*matrix_cohort_B[-1,1]%*%t(matrix_cohort_B[1,-1]))*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,cohorts_count=cohorts_count)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),kappa_rest_mat=matrix(0,nrow=p,ncol=T-1),gamma_rest=rep(0,C-2)))
  vars<-c("q","alpha","kappa","gamma")
  logit_CBD_M3_jags<-jags.model("logit_CBD_M3.jags",data=data,inits=inits,n.chain=1)
  fit_logit_CBD_M3_jags<-coda.samples(logit_CBD_M3_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"CBD_M3"]<-DIC_jags2_fn(mcmc_object=fit_logit_CBD_M3_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: CBD_M3"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##CBD M3 (sep)
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  cohorts_count_mat<-matrix(NA,nrow=A,ncol=T)
  for (i in 1:A){
    for(j in 1:T){
      cohorts_count_mat[i,j]<-j-i+A
    }
  }
  cohorts_count<-as.vector(table(cohorts_count_mat))
  prior_mean_kappa=rep(0,T-1)
  sigma2_kappa<-1000
  matrix_kappa_A<-diag(rep(1,T));matrix_kappa_A[1,]<-1
  matrix_kappa_B<-matrix_kappa_A%*%t(matrix_kappa_A)
  prior_prec_kappa=solve((matrix_kappa_B[-1,-1]-1/matrix_kappa_B[1,1]*matrix_kappa_B[-1,1]%*%t(matrix_kappa_B[1,-1]))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-2)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C-1))
  matrix_cohort_A<-rbind(cohorts_count[-1],matrix_cohort_A[(1:(C-2)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[-1,-1]-1/(matrix_cohort_B[1,1])*matrix_cohort_B[-1,1]%*%t(matrix_cohort_B[1,-1]))*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,cohorts_count=cohorts_count)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),kappa_rest_mat=matrix(0,nrow=p,ncol=T-1),gamma_rest_mat=matrix(0,nrow=p,ncol=C-2)))
  vars<-c("q","alpha","kappa","gamma")
  logit_CBD_M3_sep_jags<-jags.model("logit_CBD_M3_sep.jags",data=data,inits=inits,n.chain=1)
  fit_logit_CBD_M3_sep_jags<-coda.samples(logit_CBD_M3_sep_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"CBD_M3_sep"]<-DIC_jags2_fn(mcmc_object=fit_logit_CBD_M3_sep_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: CBD_M3_sep"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##CBD M5
  
  prior_mean_kappa1=prior_mean_kappa2=rep(0,T)
  sigma2_kappa<-1000
  prior_prec_kappa1=prior_prec_kappa2=solve(diag(rep(1,T))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_kappa1=prior_mean_kappa1,prior_prec_kappa1=prior_prec_kappa1,prior_mean_kappa2=prior_mean_kappa2,prior_prec_kappa2=prior_prec_kappa2,x=ages,xbar=mean(ages))
  inits<-function() (list(kappa1=matrix(0,nrow=p,ncol=T),kappa2=matrix(0,nrow=p,ncol=T)))
  vars<-c("q","kappa1","kappa2")
  logit_CBD_M5_jags<-jags.model("logit_CBD_M5.jags",data=data,inits=inits,n.chain=1)
  fit_logit_CBD_M5_jags<-coda.samples(logit_CBD_M5_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"CBD_M5"]<-DIC_jags2_fn(mcmc_object=fit_logit_CBD_M5_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: CBD_M5"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##CBD M6
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  prior_mean_kappa1=prior_mean_kappa2=rep(0,T)
  sigma2_kappa<-1000
  prior_prec_kappa1=prior_prec_kappa2=solve(diag(rep(1,T))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-2)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A<-rbind(rep(1/C,C),(cohorts-mean(cohorts))/C,matrix_cohort_A[(2:(C-1)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[3:C,3:C]-matrix_cohort_B[3:C,1:2]%*%solve(matrix_cohort_B[1:2,1:2])%*%matrix_cohort_B[1:2,3:C])*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_kappa1=prior_mean_kappa1,prior_prec_kappa1=prior_prec_kappa1,prior_mean_kappa2=prior_mean_kappa2,prior_prec_kappa2=prior_prec_kappa2,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,x=ages,xbar=mean(ages),cohorts=cohorts,cohorts_rev=cohorts_rev)
  inits<-function() (list(kappa1=matrix(0,nrow=p,ncol=T),kappa2=matrix(0,nrow=p,ncol=T),gamma_rest=rep(0,C-2)))
  vars<-c("q","kappa1","kappa2","gamma")
  logit_CBD_M6_jags<-jags.model("logit_CBD_M6.jags",data=data,inits=inits,n.chain=1)
  fit_logit_CBD_M6_jags<-coda.samples(logit_CBD_M6_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"CBD_M6"]<-DIC_jags2_fn(mcmc_object=fit_logit_CBD_M6_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: CBD_M6"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##CBD M6 (sep)
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  prior_mean_kappa1=prior_mean_kappa2=rep(0,T)
  sigma2_kappa<-1000
  prior_prec_kappa1=prior_prec_kappa2=solve(diag(rep(1,T))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-2)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A<-rbind(rep(1/C,C),(cohorts-mean(cohorts))/C,matrix_cohort_A[(2:(C-1)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[3:C,3:C]-matrix_cohort_B[3:C,1:2]%*%solve(matrix_cohort_B[1:2,1:2])%*%matrix_cohort_B[1:2,3:C])*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_kappa1=prior_mean_kappa1,prior_prec_kappa1=prior_prec_kappa1,prior_mean_kappa2=prior_mean_kappa2,prior_prec_kappa2=prior_prec_kappa2,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,x=ages,xbar=mean(ages),cohorts=cohorts,cohorts_rev=cohorts_rev)
  inits<-function() (list(kappa1=matrix(0,nrow=p,ncol=T),kappa2=matrix(0,nrow=p,ncol=T),gamma_rest_mat=matrix(0,nrow=p,ncol=C-2)))
  vars<-c("q","kappa1","kappa2","gamma")
  logit_CBD_M6_sep_jags<-jags.model("logit_CBD_M6_sep.jags",data=data,inits=inits,n.chain=1)
  fit_logit_CBD_M6_sep_jags<-coda.samples(logit_CBD_M6_sep_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"CBD_M6_sep"]<-DIC_jags2_fn(mcmc_object=fit_logit_CBD_M6_sep_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: CBD_M6_sep"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##CBD M7
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  prior_mean_kappa1=prior_mean_kappa2=prior_mean_kappa3=rep(0,T)
  sigma2_kappa<-1000
  prior_prec_kappa1=prior_prec_kappa2=prior_prec_kappa3=solve(diag(rep(1,T))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-3)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A<-rbind(rep(1/C,C),(cohorts-mean(cohorts))/sd(cohorts),(cohorts^2-mean(cohorts^2))/sd(cohorts^2),matrix_cohort_A[(2:(C-2)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[4:C,4:C]-matrix_cohort_B[4:C,1:3]%*%solve(matrix_cohort_B[1:3,1:3])%*%matrix_cohort_B[1:3,4:C])*sigma2_cohort)
  var_x=mean((ages-mean(ages))^2)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_kappa1=prior_mean_kappa1,prior_prec_kappa1=prior_prec_kappa1,prior_mean_kappa2=prior_mean_kappa2,prior_prec_kappa2=prior_prec_kappa2,prior_mean_kappa3=prior_mean_kappa3,prior_prec_kappa3=prior_prec_kappa3,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,x=ages,xbar=mean(ages),cohorts=cohorts,cohorts_rev=cohorts_rev,var_x=var_x)
  inits<-function() (list(kappa1=matrix(0,nrow=p,ncol=T),kappa2=matrix(0,nrow=p,ncol=T),kappa3=matrix(0,nrow=p,ncol=T),gamma_rest=rep(0,C-3)))
  vars<-c("q","kappa1","kappa2","kappa3","gamma")
  logit_CBD_M7_jags<-jags.model("logit_CBD_M7.jags",data=data,inits=inits,n.chain=1)
  fit_logit_CBD_M7_jags<-coda.samples(logit_CBD_M7_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"CBD_M7"]<-DIC_jags2_fn(mcmc_object=fit_logit_CBD_M7_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: CBD_M7"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##CBD M7 (sep)
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  prior_mean_kappa1=prior_mean_kappa2=prior_mean_kappa3=rep(0,T)
  sigma2_kappa<-1000
  prior_prec_kappa1=prior_prec_kappa2=prior_prec_kappa3=solve(diag(rep(1,T))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-3)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A<-rbind(rep(1/C,C),(cohorts-mean(cohorts))/sd(cohorts),(cohorts^2-mean(cohorts^2))/sd(cohorts^2),matrix_cohort_A[(2:(C-2)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[4:C,4:C]-matrix_cohort_B[4:C,1:3]%*%solve(matrix_cohort_B[1:3,1:3])%*%matrix_cohort_B[1:3,4:C])*sigma2_cohort)
  var_x=mean((ages-mean(ages))^2)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_kappa1=prior_mean_kappa1,prior_prec_kappa1=prior_prec_kappa1,prior_mean_kappa2=prior_mean_kappa2,prior_prec_kappa2=prior_prec_kappa2,prior_mean_kappa3=prior_mean_kappa3,prior_prec_kappa3=prior_prec_kappa3,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,x=ages,xbar=mean(ages),cohorts=cohorts,cohorts_rev=cohorts_rev,var_x=var_x)
  inits<-function() (list(kappa1=matrix(0,nrow=p,ncol=T),kappa2=matrix(0,nrow=p,ncol=T),kappa3=matrix(0,nrow=p,ncol=T),gamma_rest_mat=matrix(0,nrow=p,ncol=(C-3))))
  vars<-c("q","kappa1","kappa2","kappa3","gamma")
  logit_CBD_M7_sep_jags<-jags.model("logit_CBD_M7_sep.jags",data=data,inits=inits,n.chain=1)
  fit_logit_CBD_M7_sep_jags<-coda.samples(logit_CBD_M7_sep_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"CBD_M7_sep"]<-DIC_jags2_fn(mcmc_object=fit_logit_CBD_M7_sep_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: CBD_M7_sep"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##CBD M8
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  cohorts_count_mat<-matrix(NA,nrow=A,ncol=T)
  for (i in 1:A){
    for(j in 1:T){
      cohorts_count_mat[i,j]<-j-i+A
    }
  }
  cohorts_count<-as.vector(table(cohorts_count_mat))
  prior_mean_kappa1=prior_mean_kappa2=rep(0,T)
  sigma2_kappa<-1000
  prior_prec_kappa1=prior_prec_kappa2=solve(diag(rep(1,T))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-1)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A[1,]<-cohorts_count
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[-1,-1]-1/(matrix_cohort_B[1,1])*matrix_cohort_B[-1,1]%*%t(matrix_cohort_B[1,-1]))*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_kappa1=prior_mean_kappa1,prior_prec_kappa1=prior_prec_kappa1,prior_mean_kappa2=prior_mean_kappa2,prior_prec_kappa2=prior_prec_kappa2,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,x=ages,xbar=mean(ages),cohorts_count=cohorts_count)
  inits<-function() (list(xc=rep(0,p),kappa1=matrix(0,nrow=p,ncol=T),kappa2=matrix(0,nrow=p,ncol=T),gamma_rest=rep(0,C-1)))
  vars<-c("q","kappa1","kappa2","gamma","xc")
  logit_CBD_M8_jags<-jags.model("logit_CBD_M8.jags",data=data,inits=inits,n.chain=1)
  fit_logit_CBD_M8_jags<-coda.samples(logit_CBD_M8_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"CBD_M8"]<-DIC_jags2_fn(mcmc_object=fit_logit_CBD_M8_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: CBD_M8"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##CBD M8 (sep)
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  cohorts_count_mat<-matrix(NA,nrow=A,ncol=T)
  for (i in 1:A){
    for(j in 1:T){
      cohorts_count_mat[i,j]<-j-i+A
    }
  }
  cohorts_count<-as.vector(table(cohorts_count_mat))
  prior_mean_kappa1=prior_mean_kappa2=rep(0,T)
  sigma2_kappa<-1000
  prior_prec_kappa1=prior_prec_kappa2=solve(diag(rep(1,T))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-1)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A[1,]<-cohorts_count
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[-1,-1]-1/(matrix_cohort_B[1,1])*matrix_cohort_B[-1,1]%*%t(matrix_cohort_B[1,-1]))*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_kappa1=prior_mean_kappa1,prior_prec_kappa1=prior_prec_kappa1,prior_mean_kappa2=prior_mean_kappa2,prior_prec_kappa2=prior_prec_kappa2,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,x=ages,xbar=mean(ages),cohorts_count=cohorts_count)
  inits<-function() (list(xc=rep(0,p),kappa1=matrix(0,nrow=p,ncol=T),kappa2=matrix(0,nrow=p,ncol=T),gamma_rest_mat=matrix(0,nrow=p,ncol=C-1)))
  vars<-c("q","kappa1","kappa2","gamma","xc")
  logit_CBD_M8_sep_jags<-jags.model("logit_CBD_M8_sep.jags",data=data,inits=inits,n.chain=1)
  fit_logit_CBD_M8_sep_jags<-coda.samples(logit_CBD_M8_sep_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"CBD_M8_sep"]<-DIC_jags2_fn(mcmc_object=fit_logit_CBD_M8_sep_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: CBD_M8_sep"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##APCI
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  t<-1:T;t_rev<-rev(t)
  prior_mean_beta=rep(0,A)
  sigma2_beta<-1
  prior_prec_beta<-solve(diag(rep(1,A))*sigma2_beta)
  prior_mean_kappa=rep(0,T-2)
  sigma2_kappa<-1000
  matrix_kappa_A<-diag(rep(1,T))
  matrix_kappa_A<-rbind(rep(1/T,T),((1:T)-mean((1:T)))/sd((1:T)),matrix_kappa_A[(2:(T-1)),])
  matrix_kappa_B<-matrix_kappa_A%*%t(matrix_kappa_A)
  prior_prec_kappa=solve((matrix_kappa_B[3:T,3:T]-matrix_kappa_B[3:T,1:2]%*%solve(matrix_kappa_B[1:2,1:2])%*%matrix_kappa_B[1:2,3:T])*sigma2_kappa)
  prior_mean_cohort=rep(0,C-3)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A<-rbind(rep(1/C,C),(cohorts-mean(cohorts))/sd(cohorts),(cohorts^2-mean(cohorts^2))/sd(cohorts^2),matrix_cohort_A[(2:(C-2)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[4:C,4:C]-matrix_cohort_B[4:C,1:3]%*%solve(matrix_cohort_B[1:3,1:3])%*%matrix_cohort_B[1:3,4:C])*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,t=t,tbar=mean(t),t_rev=t_rev,cohorts=cohorts,cohorts_rev=cohorts_rev)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),beta=matrix(0,nrow=p,ncol=A),kappa_rest_mat=matrix(0,nrow=p,ncol=T-2),gamma_rest=rep(0,C-3)))
  vars<-c("q","alpha","beta","kappa","gamma")
  logit_APCI_jags<-jags.model("logit_APCI.jags",data=data,inits=inits,n.chain=1)
  fit_logit_APCI_jags<-coda.samples(logit_APCI_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"APCI"]<-DIC_jags2_fn(mcmc_object=fit_logit_APCI_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: APCI"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##APCI (sep)
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  t<-1:T;t_rev<-rev(t)
  prior_mean_beta=rep(0,A)
  sigma2_beta<-1
  prior_prec_beta<-solve(diag(rep(1,A))*sigma2_beta)
  prior_mean_kappa=rep(0,T-2)
  sigma2_kappa<-1000
  matrix_kappa_A<-diag(rep(1,T))
  matrix_kappa_A<-rbind(rep(1/T,T),((1:T)-mean((1:T)))/sd((1:T)),matrix_kappa_A[(2:(T-1)),])
  matrix_kappa_B<-matrix_kappa_A%*%t(matrix_kappa_A)
  prior_prec_kappa=solve((matrix_kappa_B[3:T,3:T]-matrix_kappa_B[3:T,1:2]%*%solve(matrix_kappa_B[1:2,1:2])%*%matrix_kappa_B[1:2,3:T])*sigma2_kappa)
  prior_mean_cohort=rep(0,C-3)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A<-rbind(rep(1/C,C),(cohorts-mean(cohorts))/sd(cohorts),(cohorts^2-mean(cohorts^2))/sd(cohorts^2),matrix_cohort_A[(2:(C-2)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[4:C,4:C]-matrix_cohort_B[4:C,1:3]%*%solve(matrix_cohort_B[1:3,1:3])%*%matrix_cohort_B[1:3,4:C])*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,t=t,tbar=mean(t),t_rev=t_rev,cohorts=cohorts,cohorts_rev=cohorts_rev)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),beta=matrix(0,nrow=p,ncol=A),kappa_rest_mat=matrix(0,nrow=p,ncol=T-2),gamma_rest_mat=matrix(0,nrow=p,ncol=C-3)))
  vars<-c("q","alpha","beta","kappa","gamma")
  logit_APCI_sep_jags<-jags.model("logit_APCI_sep.jags",data=data,inits=inits,n.chain=1)
  fit_logit_APCI_sep_jags<-coda.samples(logit_APCI_sep_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"APCI_sep"]<-DIC_jags2_fn(mcmc_object=fit_logit_APCI_sep_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: APCI_sep"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##RH
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  prior_mean_beta=rep(1/A,A-1)
  sigma2_beta<-0.1
  prior_prec_beta=solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_kappa=rep(0,T-1)
  sigma2_kappa<-1000
  matrix_kappa_A<-diag(rep(1,T));matrix_kappa_A[1,]<-1
  matrix_kappa_B<-matrix_kappa_A%*%t(matrix_kappa_A)
  prior_prec_kappa=solve((matrix_kappa_B[-1,-1]-1/matrix_kappa_B[1,1]*matrix_kappa_B[-1,1]%*%t(matrix_kappa_B[1,-1]))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-2)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A<-rbind(rep(1/C,C),(cohorts-mean(cohorts))/C,matrix_cohort_A[(2:(C-1)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[3:C,3:C]-matrix_cohort_B[3:C,1:2]%*%solve(matrix_cohort_B[1:2,1:2])%*%matrix_cohort_B[1:2,3:C])*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,cohorts=cohorts,cohorts_rev=cohorts_rev)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),beta_rest_mat=matrix(0,nrow=p,ncol=A-1),kappa_rest_mat=matrix(0,nrow=p,ncol=T-1),gamma_rest=rep(0,C-2)))
  vars<-c("q","alpha","beta","kappa","gamma")
  logit_RH_jags<-jags.model("logit_RH.jags",data=data,inits=inits,n.chain=1)
  fit_logit_RH_jags<-coda.samples(logit_RH_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"RH"]<-DIC_jags2_fn(mcmc_object=fit_logit_RH_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: RH"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  ##RH (sep)
  
  C<-A+T-1
  cohorts<-1:C;cohorts_rev<-rev(cohorts)
  prior_mean_beta=rep(1/A,A-1)
  sigma2_beta<-0.1
  prior_prec_beta=solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_kappa=rep(0,T-1)
  sigma2_kappa<-1000
  matrix_kappa_A<-diag(rep(1,T));matrix_kappa_A[1,]<-1
  matrix_kappa_B<-matrix_kappa_A%*%t(matrix_kappa_A)
  prior_prec_kappa=solve((matrix_kappa_B[-1,-1]-1/matrix_kappa_B[1,1]*matrix_kappa_B[-1,1]%*%t(matrix_kappa_B[1,-1]))*sigma2_kappa)
  prior_mean_cohort=rep(0,C-2)
  sigma2_cohort<-1
  matrix_cohort_A<-diag(rep(1,C))
  matrix_cohort_A<-rbind(rep(1/C,C),(cohorts-mean(cohorts))/C,matrix_cohort_A[(2:(C-1)),])
  matrix_cohort_B<-matrix_cohort_A%*%t(matrix_cohort_A)
  prior_prec_cohort=solve((matrix_cohort_B[3:C,3:C]-matrix_cohort_B[3:C,1:2]%*%solve(matrix_cohort_B[1:2,1:2])%*%matrix_cohort_B[1:2,3:C])*sigma2_cohort)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,C=C,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa,prior_mean_cohort=prior_mean_cohort,prior_prec_cohort=prior_prec_cohort,cohorts=cohorts,cohorts_rev=cohorts_rev)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),beta_rest_mat=matrix(0,nrow=p,ncol=A-1),kappa_rest_mat=matrix(0,nrow=p,ncol=T-1),gamma_rest_mat=matrix(0,nrow=p,ncol=C-2)))
  vars<-c("q","alpha","beta","kappa","gamma")
  logit_RH_sep_jags<-jags.model("logit_RH_sep.jags",data=data,inits=inits,n.chain=1)
  fit_logit_RH_sep_jags<-coda.samples(logit_RH_sep_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"RH_sep"]<-DIC_jags2_fn(mcmc_object=fit_logit_RH_sep_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  print_colour(paste0("Completed: RH_sep"," (",iter,"/",M,")","\n"),"red")
  iter<-iter+1
  
  best_model<-colnames(DIC_table)[which.min(DIC_table)]
  worst_model<-colnames(DIC_table)[which.max(DIC_table)]
  
  result<-list(best=get(paste0("fit_logit_",best_model,"_jags"))[[1]],
               worst=get(paste0("fit_logit_",worst_model,"_jags"))[[1]])
  
  list(result=result,DIC=DIC_table,best_model=best_model,worst_model=worst_model)
  
}


##function to analyse(plot) output

analyse_output_fn<-function(mcmc_object,deaths_combined_array,expo_combined_array,ages,years){
  
  p<-dim(deaths_combined_array)[1]
  A<-dim(deaths_combined_array)[2]
  T<-dim(deaths_combined_array)[3]
  
  n<-dim(mcmc_object)[1]
  
  p_names<-dimnames(deaths_combined_array)[[1]]
  ages_names<-dimnames(deaths_combined_array)[[2]]
  years_names<-dimnames(deaths_combined_array)[[3]]
  
  crude_rates<-deaths_combined_array/expo_combined_array
  crude_rates<-provideDimnames(crude_rates,base=list(p_names,ages_names,years_names))
  
  rates_lower<-array(dim=c(p,A,T))
  rates_median<-array(dim=c(p,A,T))
  rates_upper<-array(dim=c(p,A,T))
  rates_lower<-provideDimnames(rates_lower,base=list(p_names,ages_names,years_names))
  rates_median<-provideDimnames(rates_median,base=list(p_names,ages_names,years_names))
  rates_upper<-provideDimnames(rates_upper,base=list(p_names,ages_names,years_names))
  interval_name<-c("lower","median","upper")
  
  rates_mat<-matrix(0,nrow=n,ncol=A*T)
  
  rates_mat<-mcmc_object[,startsWith(colnames(mcmc_object),"q")]
  
  rates_pn<-apply(rates_mat,2,percentile.fn)
  
  param_mat<-mcmc_object[,!(startsWith(colnames(mcmc_object),"q"))]
  #ToDo: plot param
  
  param_pn<-apply(param_mat,2,percentile.fn)
  
  #for (s in 1:3){
  #  for (i in 1:p){
  #    for(j in 1:A){
  #      for(k in 1:T){
  #        index<-paste0("q[",i,",",j,",",k,"]")
  #        assign(paste0("rates_",interval_name[s],"[",i,",",j,",",k,"]"),rates_pn[s,index])
  #        #get(paste0("rates_",interval_name[s]))[i,j,k]<-1
  #      }
  #    }
  #  }
  #}
  
  for (i in 1:p){
    for(j in 1:A){
      for(k in 1:T){
        index<-paste0("q[",i,",",j,",",k,"]")
        rates_lower[i,j,k]<-rates_pn[1,index]
        rates_median[i,j,k]<-rates_pn[2,index]
        rates_upper[i,j,k]<-rates_pn[3,index]
      }
    }
  }
  
  length_years<-length(years)
  if (length_years<=3){
  par(mfrow=c(1,length_years))}else if(length_years>3 & length_years<=6){
    par(mfrow=c(2,ceiling(length_years/2)))
  }else if(length_years>6 & length_years<=9){
    par(mfrow=c(3,3))
  }else{
    par(mfrow=c(3,3))
    print("WARNING: Too many years selected, only printing the first 9 years.")
    years<-years[1:9]
  }
  
  
  yrange_plot<-range(log(rates_pn))
  
  for (i in 1:length(years)){
    plot(NULL,xlim=range(ages),ylim=yrange_plot,main=years[i],xlab="age",ylab="log rates")
    for (j in 1:p){
    lines(ages,log(rates_median[j,,as.character(years[i])]),type="l",col=(j+1));lines(ages,log(rates_lower[j,,as.character(years[i])]),lty=2,col=(j+1));lines(ages,log(rates_upper[j,,as.character(years[i])]),lty=2,col=(j+1))
    points(ages,log(crude_rates[j,,as.character(years[i])]),col=(j+1),pch=19,cex=0.5)
    legend("bottomright",p_names,lty=1,col=((1:p)+1))
    
  }}
  
}

###############
##Data (products)
##############

load("data_summarised.rda")

head(data_summarised)
dim(data_summarised)

names(data_summarised)

table(data_summarised$Claim/data_summarised$Exposure==data_summarised$Qx)
table(data_summarised$Product)

data_summarised_ACI<-data_summarised[data_summarised$Product=="ACI",]
data_summarised_Annuities<-data_summarised[data_summarised$Product=="Annuities",]
data_summarised_DB<-data_summarised[data_summarised$Product=="DB",]
data_summarised_SCI<-data_summarised[data_summarised$Product=="SCI",]

table(data_summarised_ACI$Age)
table(data_summarised_Annuities$Age)
table(data_summarised_DB$Age)
table(data_summarised_SCI$Age)

#subset ages 35-65 (to avoid 0 number of claims)

data_subset_ACI<-data_summarised_ACI[(34<data_summarised_ACI$Age & data_summarised_ACI$Age<=65),c("Age","Year","Claim","Exposure")]
data_subset_DB<-data_summarised_DB[(34<data_summarised_DB$Age & data_summarised_DB$Age<=65),c("Age","Year","Claim","Exposure")]
data_subset_SCI<-data_summarised_SCI[(34<data_summarised_SCI$Age & data_summarised_SCI$Age<=65),c("Age","Year","Claim","Exposure")]
data_subset_Annuities<-data_summarised_Annuities[(34<data_summarised_Annuities$Age & data_summarised_Annuities$Age<=65),c("Age","Year","Claim","Exposure")]

#we round off the number of claims/deaths to be consistent with the model specification due to Poisson/Binomial later

data_subset_ACI$Qx<-round(data_subset_ACI$Claim)/data_subset_ACI$Exposure
data_subset_DB$Qx<-round(data_subset_DB$Claim)/data_subset_DB$Exposure
data_subset_SCI$Qx<-round(data_subset_SCI$Claim)/data_subset_SCI$Exposure
data_subset_Annuities$Qx<-round(data_subset_Annuities$Claim)/data_subset_Annuities$Exposure

deaths_subset_ACI_mat<-matrix(round(data_subset_ACI$Claim),nrow=length(unique(data_subset_ACI$Age)),ncol=length(unique(data_subset_ACI$Year)),byrow = TRUE)
expo_subset_ACI_mat<-matrix(round(data_subset_ACI$Exposure),nrow=length(unique(data_subset_ACI$Age)),ncol=length(unique(data_subset_ACI$Year)),byrow = TRUE)

deaths_subset_DB_mat<-matrix(round(data_subset_DB$Claim),nrow=length(unique(data_subset_DB$Age)),ncol=length(unique(data_subset_DB$Year)),byrow = TRUE)
expo_subset_DB_mat<-matrix(round(data_subset_DB$Exposure),nrow=length(unique(data_subset_DB$Age)),ncol=length(unique(data_subset_DB$Year)),byrow = TRUE)

deaths_subset_SCI_mat<-matrix(round(data_subset_SCI$Claim),nrow=length(unique(data_subset_SCI$Age)),ncol=length(unique(data_subset_SCI$Year)),byrow = TRUE)
expo_subset_SCI_mat<-matrix(round(data_subset_SCI$Exposure),nrow=length(unique(data_subset_SCI$Age)),ncol=length(unique(data_subset_SCI$Year)),byrow = TRUE)

p<-3 #3 products
A<-dim(deaths_subset_ACI_mat)[1];T<-dim(deaths_subset_ACI_mat)[2]
deaths_combined_array<-array(dim=c(p,A,T))
deaths_combined_array[1,,]<-deaths_subset_ACI_mat;deaths_combined_array[2,,]<-deaths_subset_DB_mat;deaths_combined_array[3,,]<-deaths_subset_SCI_mat
expo_combined_array<-array(dim=c(p,A,T))
expo_combined_array[1,,]<-expo_subset_ACI_mat;expo_combined_array[2,,]<-expo_subset_DB_mat;expo_combined_array[3,,]<-expo_subset_SCI_mat

ages=unique(data_subset_ACI$Age);years=unique(data_subset_ACI$Year)
deaths_combined_array<-provideDimnames(deaths_combined_array,base=list(c("ACI","DB","SCI"),as.character(ages),as.character(years)))
expo_combined_array<-provideDimnames(expo_combined_array,base=list(c("ACI","DB","SCI"),as.character(ages),as.character(years)))

##########
##Analyse the data (products)
##########

#run the auto model-selection function to find the best Bayesian model
result_auto_fit<-auto_fit_DIC_fn(deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,n_iter=10000)
result_auto_fit$DIC
par(mfrow=c(1,1),mar=c(10,4,4,2))
col_min<-rep(1,length(result_auto_fit$DIC));col_min[which.min(result_auto_fit$DIC[1,])]<-2
barplot(result_auto_fit$DIC[1,],ylim=c(0,4000),las=3,col=col_min)
result_auto_fit$best_model
result_auto_fit$worst_model

#plot the output(fitted rates)
#best model results
analyse_output_fn(mcmc_object=result_auto_fit$result$best,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=unique(data_subset_ACI$Year))
#worst model results
analyse_output_fn(mcmc_object=result_auto_fit$result$worst,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=unique(data_subset_ACI$Year))

#focusing on selected years
#best model results
analyse_output_fn(mcmc_object=result_auto_fit$result$best,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=c(2016,2019,2020))
#worst model results
analyse_output_fn(mcmc_object=result_auto_fit$result$worst,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=c(2016,2019,2020))
gc()

###############
##Data (gender/sex)
##############

load("data_sex_UK.rda")

dim(d) #1=male,2=female
dim(E)
head(d[,,1]);head(E[,,2])
tail(d[,,1]);tail(E[,,2]) #quite a few zero deaths for male
rownames(d[,,1]);rownames(E[,,1])
colnames(d[,,1]);colnames(E[,,1])

dxt<-array(dim=c(dim(d)[3],dim(d)[1],dim(d)[2]));Ext<-array(dim=c(dim(E)[3],dim(E)[1],dim(E)[2]))
dxt[1,,]<-d[,,1];Ext[1,,]<-E[,,1]
dxt[2,,]<-d[,,2];Ext[2,,]<-E[,,2]

#subset ages 0-100 (to avoid 0 number of claims), years 1961-2021
dxt_array_sex<-round(dxt[,1:101,121:181])  #round off
Ext_array_sex<-round(Ext[,1:101,121:181])

p<-dim(dxt_array_sex)[1];A<-dim(dxt_array_sex)[2];T<-dim(dxt_array_sex)[3]
ages<-0:100;years<-as.numeric(colnames(d[,,1]))[121:181]

dxt_array_sex<-provideDimnames(dxt_array_sex,base=list(c("male","female"),as.character(ages),as.character(years)))
Ext_array_sex<-provideDimnames(Ext_array_sex,base=list(c("male","female"),as.character(ages),as.character(years)))

##########
##Analyse the data (gender/sex)
##########

#for some reasons,the following cell has larger number of death than exposure?(Check!) so i just set a random larger number
#Ext_array_sex[1,101,113]<-30

#run the auto model-selection function to find the best Bayesian model
result_auto_fit_sex<-auto_fit_DIC_fn(deaths_combined_array=dxt_array_sex,expo_combined_array=Ext_array_sex,n_iter=5000)
result_auto_fit_sex$DIC
result_auto_fit_sex$best_model
result_auto_fit_sex$worst_model
par(mfrow=c(1,1),mar=c(10,4,4,2))
col_min<-rep(1,length(result_auto_fit_sex$DIC));col_min[which.min(result_auto_fit_sex$DIC[1,])]<-2
barplot(result_auto_fit_sex$DIC[1,],las=3,col=col_min)

#plot the output(fitted rates)
#best model results
#analyse_output_fn(mcmc_object=result_auto_fit_sex$result$best,deaths_combined_array=dxt_array_sex,expo_combined_array=Ext_array_sex,ages=as.numeric(dimnames(dxt_array_sex)[[2]]),years=years)
#worst model results
#analyse_output_fn(mcmc_object=result_auto_fit_sex$result$worst,deaths_combined_array=dxt_array_sex,expo_combined_array=Ext_array_sex,ages=as.numeric(dimnames(dxt_array_sex)[[2]]),years=years)

#focusing on selected years
#best model results
#analyse_output_fn(mcmc_object=result_auto_fit_sex$result$best,deaths_combined_array=dxt_array_sex,expo_combined_array=Ext_array_sex,ages=as.numeric(dimnames(dxt_array_sex)[[2]]),years=c(1961,1981,2001))
#worst model results
#analyse_output_fn(mcmc_object=result_auto_fit_sex$result$worst,deaths_combined_array=dxt_array_sex,expo_combined_array=Ext_array_sex,ages=as.numeric(dimnames(dxt_array_sex)[[2]]),years=c(1961,1981,2001))


###############
##Data (countries)
##############

load("data_country.rda")

p<-dim(dxt_array_country)[1];A<-dim(dxt_array_country)[2];T<-dim(dxt_array_country)[3]
ages<-as.numeric(dimnames(dxt_array_country)[[2]]);years<-as.numeric(dimnames(dxt_array_country)[[3]])

#library(demography)
#rates_AUS_male_demog<-demogdata((dxt_array_country/Ext_array_country)[1,,],pop=Ext_array_country[1,,],ages=ages,years=years,label="AUS_male",type="mortality",name="AUS_male")
#plot(rates_AUS_male_demog)
#plot(lca(rates_AUS_male_demog))
#rates_US_male_demog<-demogdata((dxt_array_country/Ext_array_country)[5,,],pop=Ext_array_country[5,,],ages=ages,years=years,label="US_male",type="mortality",name="US_male")
#plot(rates_US_male_demog)
#plot(lca(rates_US_male_demog))

##########
##Analyse the data (countries)
##########

#for some reasons,the following cell has larger number of death than exposure?(Check!) so i just set a random larger number
#Ext_array_sex[1,101,113]<-30

#run the auto model-selection function to find the best Bayesian model
result_auto_fit_country<-auto_fit_DIC_fn(deaths_combined_array=dxt_array_country,expo_combined_array=Ext_array_country,n_iter=5000)
result_auto_fit_country$DIC
result_auto_fit_country$best_model
result_auto_fit_country$worst_model
par(mfrow=c(1,1),mar=c(10,4,4,2))
col_min<-rep(1,length(result_auto_fit_country$DIC));col_min[which.min(result_auto_fit_country$DIC[1,])]<-2
barplot(result_auto_fit_country$DIC[1,],las=3,col=col_min)

#plot the output(fitted rates)
#best model results
analyse_output_fn(mcmc_object=result_auto_fit_country$result$best,deaths_combined_array=dxt_array_country,expo_combined_array=Ext_array_country,ages=as.numeric(dimnames(dxt_array_country)[[2]]),years=years)
#worst model results
analyse_output_fn(mcmc_object=result_auto_fit_country$result$worst,deaths_combined_array=dxt_array_country,expo_combined_array=Ext_array_country,ages=as.numeric(dimnames(dxt_array_country)[[2]]),years=years)

#focusing on selected years
#best model results
analyse_output_fn(mcmc_object=result_auto_fit_country$result$best,deaths_combined_array=dxt_array_country,expo_combined_array=Ext_array_country,ages=as.numeric(dimnames(dxt_array_country)[[2]]),years=c(2000))
#worst model results
analyse_output_fn(mcmc_object=result_auto_fit_country$result$worst,deaths_combined_array=dxt_array_country,expo_combined_array=Ext_array_country,ages=as.numeric(dimnames(dxt_array_country)[[2]]),years=c(1960,1981,2000))

#############
#saving plots
#############

#uncomment and run if needed

#product data
#pdf(file="auto_DIC_product_barplot.pdf",width=8,height=5)
#par(mfrow=c(1,1),mar=c(10,4,4,2))
#col_min<-rep(1,length(result_auto_fit$DIC));col_min[which.min(result_auto_fit$DIC[1,])]<-2
#barplot(result_auto_fit$DIC[1,],ylim=c(0,4000),las=3,col=col_min)
#dev.off()

#pdf(file="auto_DIC_product.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit$result$best,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=unique(data_subset_ACI$Year))
#dev.off()

#pdf(file="auto_DIC_product_selectedyears.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit$result$best,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=c(2016,2019,2020))
#dev.off()

#sex data
#pdf(file="auto_DIC_sex_barplot.pdf",width=8,height=5)
#par(mfrow=c(1,1),mar=c(10,4,4,2))
#col_min<-rep(1,length(result_auto_fit_sex$DIC));col_min[which.min(result_auto_fit_sex$DIC[1,])]<-2
#barplot(result_auto_fit_sex$DIC[1,],las=3,col=col_min,ylim=c(0,600000))
#dev.off()

#pdf(file="auto_DIC_sex.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit_sex$result$best,deaths_combined_array=dxt_array_sex,expo_combined_array=Ext_array_sex,ages=as.numeric(dimnames(dxt_array_sex)[[2]]),years=c(1961,1981,2001))
#dev.off()

#pdf(file="auto_DIC_sex_y1961.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit_sex$result$best,deaths_combined_array=dxt_array_sex,expo_combined_array=Ext_array_sex,ages=as.numeric(dimnames(dxt_array_sex)[[2]]),years=c(1961))
#dev.off()

#pdf(file="auto_DIC_sex_y1981.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit_sex$result$best,deaths_combined_array=dxt_array_sex,expo_combined_array=Ext_array_sex,ages=as.numeric(dimnames(dxt_array_sex)[[2]]),years=c(1981))
#dev.off()

#pdf(file="auto_DIC_sex_y2001.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit_sex$result$best,deaths_combined_array=dxt_array_sex,expo_combined_array=Ext_array_sex,ages=as.numeric(dimnames(dxt_array_sex)[[2]]),years=c(2001))
#dev.off()

#country data
#pdf(file="auto_DIC_country_barplot.pdf",width=8,height=5)
#par(mfrow=c(1,1),mar=c(10,4,4,2))
#col_min<-rep(1,length(result_auto_fit_country$DIC));col_min[which.min(result_auto_fit_country$DIC[1,])]<-2
#barplot(result_auto_fit_country$DIC[1,],las=3,col=col_min,ylim=c(0,1000000))
#dev.off()

#pdf(file="auto_DIC_country.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit_country$result$best,deaths_combined_array=dxt_array_country,expo_combined_array=Ext_array_country,ages=as.numeric(dimnames(dxt_array_country)[[2]]),years=c(1951,1971,1991))
#dev.off()

#pdf(file="auto_DIC_country_y1951.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit_country$result$best,deaths_combined_array=dxt_array_country,expo_combined_array=Ext_array_country,ages=as.numeric(dimnames(dxt_array_country)[[2]]),years=c(1951))
#dev.off()

#pdf(file="auto_DIC_country_y1971.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit_country$result$best,deaths_combined_array=dxt_array_country,expo_combined_array=Ext_array_country,ages=as.numeric(dimnames(dxt_array_country)[[2]]),years=c(1971))
#dev.off()

#pdf(file="auto_DIC_country_y1991.pdf",width=8,height=5)
#analyse_output_fn(mcmc_object=result_auto_fit_country$result$best,deaths_combined_array=dxt_array_country,expo_combined_array=Ext_array_country,ages=as.numeric(dimnames(dxt_array_country)[[2]]),years=c(1991))
#dev.off()



