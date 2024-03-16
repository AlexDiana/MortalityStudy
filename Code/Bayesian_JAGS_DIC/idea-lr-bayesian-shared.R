

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
  
  DIC_table<-matrix(0,nrow=1,ncol=10)
  colnames(DIC_table)<-c("M1A","M1U","M1M","M2A1","M2A2","M2Y1","M2Y2","MLiLee","MLiLee_modified","all_separately")
  
  #Model 1A
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.0001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1),c_p_rest=rep(0,p-1)))
  vars<-c("q","alpha","beta","kappa","c_p")
  logit_LC_M1A_jags<-jags.model("G:/IDEA-LR/logit_LC_M1A.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M1A_jags<-coda.samples(logit_LC_M1A_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"M1A"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M1A_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  #Model 1U
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.0001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1)))
  vars<-c("q","alpha","beta","kappa")
  logit_LC_M1U_jags<-jags.model("G:/IDEA-LR/logit_LC_M1U.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M1U_jags<-coda.samples(logit_LC_M1U_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"M1U"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M1U_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  #Model 1M
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.0001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1),c_p_rest=rep(0,p-1)))
  vars<-c("q","alpha","beta","kappa","c_p")
  logit_LC_M1M_jags<-jags.model("G:/IDEA-LR/logit_LC_M1M.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M1M_jags<-coda.samples(logit_LC_M1M_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"M1M"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M1M_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  #Model 2A1
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.0001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1),c_p_rest=rep(0,p-1)))
  vars<-c("q","alpha","beta","kappa","c_p")
  logit_LC_M2A1_jags<-jags.model("G:/IDEA-LR/logit_LC_M2A1.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M2A1_jags<-coda.samples(logit_LC_M2A1_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"M2A1"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M2A1_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  #Model 2A2
  
  prior_mean_beta<-rep(1/(p*A),p*A-1)
  sigma2_beta<-0.0001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,p*A-1))-1/(p*A)*(matrix(1,nrow=(p*A-1),ncol=(p*A-1)))))
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/(p*A),p*A-1),kappa_rest=rep(0,T-1)))
  vars<-c("q","alpha","beta","kappa")
  logit_LC_M2A2_jags<-jags.model("G:/IDEA-LR/logit_LC_M2A2.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M2A2_jags<-coda.samples(logit_LC_M2A2_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"M2A2"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M2A2_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  #Model 2Y1
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.0001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,T-1),c_p_rest=rep(0,p-1)))
  vars<-c("q","alpha","beta","kappa","c_p")
  logit_LC_M2Y1_jags<-jags.model("G:/IDEA-LR/logit_LC_M2Y1.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M2Y1_jags<-coda.samples(logit_LC_M2Y1_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"M2Y1"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M2Y1_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  #Model 2Y2
  
  prior_mean_beta<-rep(1/A,A-1)
  sigma2_beta<-0.0001
  prior_prec_beta<-solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_kappa<-rep(0,p*T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,p*T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_kappa<-solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_beta=prior_mean_beta,prior_prec_beta=prior_prec_beta,prior_mean_kappa=prior_mean_kappa,prior_prec_kappa=prior_prec_kappa)
  inits<-function() (list(alpha=rep(0,A),beta_rest=rep(1/A,A-1),kappa_rest=rep(0,p*T-1)))
  vars<-c("q","alpha","beta","kappa")
  logit_LC_M2Y2_jags<-jags.model("G:/IDEA-LR/logit_LC_M2Y2.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_M2Y2_jags<-coda.samples(logit_LC_M2Y2_jags,vars,n.iter=n_iter,thin=1)
  DIC_table[1,"M2Y2"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_M2Y2_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  #Model LiLee
  
  prior_mean_Beta=rep(1/A,A-1)
  sigma2_beta<-0.0001
  prior_prec_Beta=solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_Kappa=rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_Kappa=solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_Beta=prior_mean_Beta,prior_prec_Beta=prior_prec_Beta,prior_mean_Kappa=prior_mean_Kappa,prior_prec_Kappa=prior_prec_Kappa)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),Beta_rest=rep(1/A,A-1),beta_rest_mat=matrix(1/A,nrow=p-1,ncol=A-1),Kappa_rest=rep(0,T-1),kappa_rest_mat=matrix(0,nrow=p-1,ncol=T-1)))
  vars<-c("q","alpha","beta","kappa","Beta","Kappa")
  logit_LC_MLiLee_jags<-jags.model("G:/IDEA-LR/logit_LC_MLiLee.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_MLiLee_jags<-coda.samples(logit_LC_MLiLee_jags,vars,n.iter=10000,thin=1)
  DIC_table[1,"MLiLee"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_MLiLee_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  ##Model LiLee modified 
  
  prior_mean_Beta=rep(1/A,A-1)
  sigma2_beta<-0.0001
  prior_prec_Beta=solve(sigma2_beta*(diag(rep(1,A-1))-1/A*(matrix(1,nrow=A-1,ncol=A-1))))
  prior_mean_Kappa=rep(0,T-1)
  sigma2_kappa<-1000
  C<-diag(rep(1,T));C[1,]<-1
  B<-C%*%t(C)
  prior_prec_Kappa=solve((B[-1,-1]-1/B[1,1]*B[-1,1]%*%t(B[1,-1]))*sigma2_kappa)
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p,prior_mean_Beta=prior_mean_Beta,prior_prec_Beta=prior_prec_Beta,prior_mean_Kappa=prior_mean_Kappa,prior_prec_Kappa=prior_prec_Kappa)
  inits<-function() (list(alpha=rep(0,A),Beta_rest=rep(1/A,A-1),beta_rest_mat=matrix(0,nrow=p-1,ncol=A-1),Kappa_rest=rep(0,T-1),kappa_rest_mat=matrix(0,nrow=p-1,ncol=T-1)))
  vars<-c("q","alpha","beta","kappa","Beta","Kappa")
  logit_LC_MLiLee_modified_jags<-jags.model("G:/IDEA-LR/logit_LC_MLiLee_modified.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_MLiLee_modified_jags<-coda.samples(logit_LC_MLiLee_modified_jags,vars,n.iter=10000,thin=1)
  DIC_table[1,"MLiLee_modified"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_MLiLee_modified_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  ##all LC separately
  
  data<-list(dxt=deaths_combined_array,ext=expo_combined_array,A=A,T=T,p=p)
  inits<-function() (list(alpha=matrix(0,nrow=p,ncol=A),beta_rest_mat=matrix(1/A,nrow=p,ncol=(A-1)),kappa_rest_mat=matrix(0,nrow=p,ncol=(T-1))))
  vars<-c("q","alpha","beta","kappa")
  logit_LC_all_separately_jags<-jags.model("G:/IDEA-LR/logit_LC_all_separately.jags",data=data,inits=inits,n.chain=1)
  fit_logit_LC_all_separately_jags<-coda.samples(logit_LC_all_separately_jags,vars,n.iter=10000,thin=1)
  DIC_table[1,"all_separately"]<-DIC_jags2_fn(mcmc_object=fit_logit_LC_all_separately_jags[[1]],deaths_combined_vectorised=deaths_combined_vectorised,expo_combined_vectorised =expo_combined_vectorised,A=A,T=T,p=p)
  
  best_model<-colnames(DIC_table)[which.min(DIC_table)]
  worst_model<-colnames(DIC_table)[which.max(DIC_table)]
  
  result<-list(best=get(paste0("fit_logit_LC_",best_model,"_jags"))[[1]],
               worst=get(paste0("fit_logit_LC_",worst_model,"_jags"))[[1]])
  
  list(result=result,DIC=DIC_table,best_model=best_model,worst_model=worst_model)
  
}


##function to analyse(plot) output

analyse_output_fn<-function(mcmc_object,deaths_combined_array,expo_combined_array,ages,years){
  
  p<-dim(deaths_combined_array)[1]
  A<-dim(deaths_combined_array)[2]
  T<-dim(deaths_combined_array)[3]
  
  n<-dim(mcmc_object)[1]
  
  crude_rates<-deaths_combined_array/expo_combined_array
  
  rates_lower<-array(dim=c(p,A,T))
  rates_median<-array(dim=c(p,A,T))
  rates_upper<-array(dim=c(p,A,T))
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
    plot(ages,log(rates_median[1,,i]),type="l",col=2,ylim=yrange_plot,main=years[i],ylab="log rates");lines(ages,log(rates_lower[1,,i]),lty=2,col=2);lines(ages,log(rates_upper[1,,i]),lty=2,col=2)
    points(ages,log(crude_rates[1,,i]),col=2,pch=19)
    
    lines(ages,log(rates_median[2,,i]),type="l",col=3,ylim=yrange_plot,main=years[i],ylab="log rates");lines(ages,log(rates_lower[2,,i]),lty=2,col=3);lines(ages,log(rates_upper[2,,i]),lty=2,col=3)
    points(ages,log(crude_rates[2,,i]),col=3,pch=19)
    
    lines(ages,log(rates_median[3,,i]),type="l",col=4,ylim=yrange_plot,main=years[i],ylab="log rates");lines(ages,log(rates_lower[3,,i]),lty=2,col=4);lines(ages,log(rates_upper[3,,i]),lty=2,col=4)
    points(ages,log(crude_rates[3,,i]),col=4,pch=19)
    legend("bottomright",c("ACI","DB","SCI"),lty=1,col=2:4)
    
  }
  
}

###############
##Data
##############

load("G:/IDEA-LR/data_summarised.rda")

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

##load packages
library(rjags)

##########
##Analyse the data
##########

#run the auto model-selection function to find the best Bayesian model
result_auto_fit<-auto_fit_DIC_fn(deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,n_iter=10000)
result_auto_fit$DIC
result_auto_fit$best_model
result_auto_fit$worst_model

#plot the output(fitted rates)
analyse_output_fn(mcmc_object=result_auto_fit$result$best,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=unique(data_subset_ACI$Year))
analyse_output_fn(mcmc_object=result_auto_fit$result$worst,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=unique(data_subset_ACI$Year))

#focusing on selected years
analyse_output_fn(mcmc_object=result_auto_fit$result$best,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=c(2016,2019,2020))
analyse_output_fn(mcmc_object=result_auto_fit$result$worst,deaths_combined_array=deaths_combined_array,expo_combined_array=expo_combined_array,ages=unique(data_subset_ACI$Age),years=c(2016,2019,2020))
