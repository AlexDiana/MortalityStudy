


########################
##Model 1A
########################

iteration.methodB.LC.poisson_M1A<-function(deaths_combine,exposure_combine,m,beta.initial,A,T){
  beta.new<-beta.initial
  alpha.new<-NULL
  k.new<-NULL
  Deviance<-vector(length=m)
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
  model_bic<-NULL
  
  design.matrix.intercept<-kronecker(X=matrix(c(1,0,-1,0,1,-1),ncol=2),Y=matrix(rep(1,A*T),ncol=1))
  
  for (i in 1:m){
    design.matrix.kappa<-kronecker(X=diag(rep(1,T)),Y=beta.new)
    design.matrix.kappa<-design.matrix.kappa-design.matrix.kappa[,1]
    design.matrix.kappa<-design.matrix.kappa[,-1]
    design.matrix.1<-cbind(design.matrix.alpha,design.matrix.kappa)
    design.matrix.A<-cbind(rbind(design.matrix.1,design.matrix.1,design.matrix.1),design.matrix.intercept)
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    alpha.new<-fit$coef[1:A]
    k.new<-fit$coef[(A+1):(A+T-1)]
    k.new<-c(-sum(k.new),k.new)
    
    design.matrix.beta<-kronecker(X=k.new,Y=diag(rep(1,A)))
    design.matrix.beta<-design.matrix.beta-design.matrix.beta[,1]
    design.matrix.beta<-design.matrix.beta[,-1]
    design.matrix.2<-cbind(design.matrix.alpha,design.matrix.beta)
    design.matrix.B<-cbind(rbind(design.matrix.2,design.matrix.2,design.matrix.2),design.matrix.intercept)
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)),family=poisson)
    alpha.new<-fit2$coef[1:A]
    beta.new<-fit2$coef[(A+1):(2*A-1)]
    beta.new<-c(1-sum(beta.new),beta.new)
    intercept<-fit2$coef[-(1:(2*A-1))]
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha=alpha.new,beta=beta.new,kappa=k.new,intercept=intercept,deviance=Deviance,model_bic=model_bic)
}



########################
##Model 1U
########################

iteration.methodB.LC.poisson_M1U<-function(deaths_combine,exposure_combine,m,beta.initial,A,T){
  beta.new<-beta.initial
  alpha.new<-NULL
  k.new<-NULL
  Deviance<-vector(length=m)
  model_bic<-NULL
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
  design.matrix.alpha<-kronecker(X=diag(rep(1,3)),Y=design.matrix.alpha)
  
  for (i in 1:m){
    design.matrix.kappa<-kronecker(X=diag(rep(1,T)),Y=beta.new)
    design.matrix.kappa<-design.matrix.kappa-design.matrix.kappa[,1]
    design.matrix.kappa<-design.matrix.kappa[,-1]
    design.matrix.A<-cbind(design.matrix.alpha,rbind(design.matrix.kappa,design.matrix.kappa,design.matrix.kappa))
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    k.new<-fit$coef[(3*A+1):(3*A+T-1)]
    k.new<-c(-sum(k.new),k.new)
    
    design.matrix.beta<-kronecker(X=k.new,Y=diag(rep(1,A)))
    design.matrix.beta<-design.matrix.beta-design.matrix.beta[,1]
    design.matrix.beta<-design.matrix.beta[,-1]
    design.matrix.B<-cbind(design.matrix.alpha,rbind(design.matrix.beta,design.matrix.beta,design.matrix.beta))
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)),family=poisson)
    alpha.prod1.new<-fit2$coef[1:A]
    alpha.prod2.new<-fit2$coef[(A+1):(2*A)]
    alpha.prod3.new<-fit2$coef[(2*A+1):(3*A)]
    beta.new<-fit2$coef[(3*A+1):(4*A-1)]
    beta.new<-c(1-sum(beta.new),beta.new)
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha.prod1=alpha.prod1.new,alpha.prod2=alpha.prod2.new,alpha.prod3=alpha.prod3.new,beta=beta.new,kappa=k.new,deviance=Deviance,model_bic=model_bic)
}

########################
##Model 1M
########################


iteration.methodB.LC.poisson_M1M<-function(deaths_combine,exposure_combine,m,beta.initial,c.initial,A,T){
  beta.new<-beta.initial
  alpha.new<-NULL
  k.new<-NULL
  c.new<-c.initial
  Deviance<-vector(length=m)
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
  model_bic<-NULL
  
  for (i in 1:m){
    design.matrix.kappa<-kronecker(X=diag(rep(1,T)),Y=beta.new)
    design.matrix.kappa<-design.matrix.kappa-design.matrix.kappa[,1]
    design.matrix.kappa<-design.matrix.kappa[,-1]
    design.matrix.1<-cbind(design.matrix.alpha*(1+c.new[1]),design.matrix.kappa)
    design.matrix.2<-cbind(design.matrix.alpha*(1+c.new[2]),design.matrix.kappa)
    #design.matrix.3<-cbind(design.matrix.alpha*(1-c.new[1]-c.new[2]),design.matrix.kappa)
    design.matrix.3<-cbind(design.matrix.alpha*(1+c.new[3]),design.matrix.kappa)
    design.matrix.A<-rbind(design.matrix.1,design.matrix.2,design.matrix.3)
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    alpha.new<-fit$coef[1:A]
    k.new<-fit$coef[(A+1):(A+T-1)]
    k.new<-c(-sum(k.new),k.new)
    
    design.matrix.beta<-kronecker(X=k.new,Y=diag(rep(1,A)))
    design.matrix.beta<-design.matrix.beta-design.matrix.beta[,1]
    design.matrix.beta<-design.matrix.beta[,-1]
    
    design.vector.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=alpha.new)
    #design.matrix.c<-kronecker(X=matrix(c(1,0,-1,0,1,-1),ncol=2),Y=design.vector.alpha)
    design.matrix.c<-kronecker(X=diag(rep(1,3)),Y=design.vector.alpha)
    design.matrix.B<-cbind(design.matrix.c,rbind(design.matrix.2,design.matrix.2,design.matrix.2))
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)+rep(design.vector.alpha,3)),family=poisson)
    c.new<-fit2$coef[1:3]
    beta.new<-fit2$coef[4:(A+2)]
    beta.new<-c(1-sum(beta.new),beta.new)
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha=alpha.new,beta=beta.new,kappa=k.new,c.new=c.new,deviance=Deviance,model_bic=model_bic)
}

########################
##Model 1M2
########################

iteration.methodB.LC.poisson_M1M2<-function(deaths_combine,exposure_combine,m,beta.initial,c.initial,A,T){
  beta.new<-beta.initial
  alpha.new<-NULL
  k.new<-NULL
  c.new<-c.initial
  Deviance<-vector(length=m)
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  model_bic<-NULL
  
  for (i in 1:m){
    design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
    design.matrix.alpha<-kronecker(X=matrix(c(c.new[1],c.new[2],(1-c.new[1]-c.new[2])),ncol=1),Y=design.matrix.alpha)
    
    design.matrix.kappa<-kronecker(X=diag(rep(1,T)),Y=beta.new)
    design.matrix.kappa<-design.matrix.kappa-design.matrix.kappa[,1]
    design.matrix.kappa<-design.matrix.kappa[,-1]
    design.matrix.kappa<-kronecker(X=matrix(rep(1,3),ncol=1),Y=design.matrix.kappa)
    design.matrix.A<-cbind(design.matrix.alpha,design.matrix.kappa)
    
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    alpha.new<-fit$coef[1:A]
    k.new<-fit$coef[(A+1):(A+T-1)]
    k.new<-c(-sum(k.new),k.new)
    
    design.matrix.c<-kronecker(matrix(c(1,0,-1,0,1,-1),ncol=2),Y=kronecker(X=matrix(rep(1,T),ncol=1),Y=alpha.new))
    
    design.matrix.beta<-kronecker(X=k.new,Y=diag(rep(1,A)))
    design.matrix.beta<-design.matrix.beta-design.matrix.beta[,1]
    design.matrix.beta<-design.matrix.beta[,-1]
    design.matrix.beta<-kronecker(X=matrix(rep(1,3),ncol=1),Y=design.matrix.beta)
    
    design.matrix.B<-cbind(design.matrix.c,design.matrix.beta)
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)+kronecker(X=matrix(c(0,0,1),ncol=1),Y=kronecker(X=matrix(rep(1,T),ncol=1),Y=alpha.new))),family=poisson)
    c.new<-fit2$coef[1:2]
    c.new<-c(c.new,1-sum(c.new))
    beta.new<-fit2$coef[(3):(A+1)]
    beta.new<-c(1-sum(beta.new),beta.new)
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha=alpha.new,beta=beta.new,kappa=k.new,c.new=c.new,deviance=Deviance,model_bic=model_bic)
}

########################
##Model 2A1
########################

iteration.methodB.LC.poisson_M2A1<-function(deaths_combine,exposure_combine,m,beta.initial,c.initial,A,T){
  beta.new<-beta.initial
  alpha.new<-NULL
  k.new<-NULL
  c.new<-c.initial
  Deviance<-vector(length=m)
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
  model_bic<-NULL
  
  for (i in 1:m){
    design.matrix.kappa.1<-kronecker(X=diag(rep(1,T)),Y=(beta.new+c.new[1]))
    design.matrix.kappa.1<-design.matrix.kappa.1-design.matrix.kappa.1[,1]
    design.matrix.kappa.1<-design.matrix.kappa.1[,-1]
    design.matrix.kappa.2<-kronecker(X=diag(rep(1,T)),Y=(beta.new+c.new[2]))
    design.matrix.kappa.2<-design.matrix.kappa.2-design.matrix.kappa.2[,1]
    design.matrix.kappa.2<-design.matrix.kappa.2[,-1]
    design.matrix.kappa.3<-kronecker(X=diag(rep(1,T)),Y=(beta.new-c.new[1]-c.new[2]))
    design.matrix.kappa.3<-design.matrix.kappa.3-design.matrix.kappa.3[,1]
    design.matrix.kappa.3<-design.matrix.kappa.3[,-1]
    design.matrix.1<-cbind(design.matrix.alpha,design.matrix.kappa.1)
    design.matrix.2<-cbind(design.matrix.alpha,design.matrix.kappa.2)
    design.matrix.3<-cbind(design.matrix.alpha,design.matrix.kappa.3)
    design.matrix.A<-rbind(design.matrix.1,design.matrix.2,design.matrix.3)
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    alpha.new<-fit$coef[1:A]
    k.new<-fit$coef[(A+1):(A+T-1)]
    k.new<-c(-sum(k.new),k.new)
    
    design.vector.kappa<-kronecker(X=k.new,Y=matrix(rep(1,A),ncol=1))
    design.matrix.c<-kronecker(X=matrix(c(1,0,-1,0,1,-1),ncol=2),Y=design.vector.kappa)
    
    design.matrix.beta<-kronecker(X=k.new,Y=diag(rep(1,A)))
    design.matrix.beta<-design.matrix.beta-design.matrix.beta[,1]
    design.matrix.beta<-design.matrix.beta[,-1]
    design.matrix.2<-cbind(design.matrix.alpha,design.matrix.beta)
    design.matrix.B<-cbind(rbind(design.matrix.2,design.matrix.2,design.matrix.2),design.matrix.c)
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)),family=poisson)
    alpha.new<-fit2$coef[1:A]
    beta.new<-fit2$coef[(A+1):(2*A-1)]
    beta.new<-c(1-sum(beta.new),beta.new)
    c.new<-fit2$coef[-(1:(2*A-1))]
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha=alpha.new,beta=beta.new,kappa=k.new,c.new=c.new,deviance=Deviance,model_bic=model_bic)
}


########################
##Model 2Y1
########################

iteration.methodB.LC.poisson_M2Y1<-function(deaths_combine,exposure_combine,m,beta.initial,c.initial,A,T){
  beta.new<-beta.initial
  alpha.new<-NULL
  k.new<-NULL
  c.new<-c.initial
  Deviance<-vector(length=m)
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
  design.matrix.alpha<-kronecker(X=matrix(rep(1,3),ncol=1),Y=design.matrix.alpha)
  model_bic<-NULL
  
  for (i in 1:m){
    design.matrix.kappa.1<-kronecker(X=diag(rep(1,T)),Y=(beta.new))
    design.matrix.kappa.1<-design.matrix.kappa.1-design.matrix.kappa.1[,1]
    design.matrix.kappa.1<-design.matrix.kappa.1[,-1]
    design.matrix.kappa.2<-kronecker(X=diag(rep(1,T)),Y=(beta.new))
    design.matrix.kappa.2<-design.matrix.kappa.2-design.matrix.kappa.2[,1]
    design.matrix.kappa.2<-design.matrix.kappa.2[,-1]
    design.matrix.kappa.3<-kronecker(X=diag(rep(1,T)),Y=(beta.new))
    design.matrix.kappa.3<-design.matrix.kappa.3-design.matrix.kappa.3[,1]
    design.matrix.kappa.3<-design.matrix.kappa.3[,-1]
    design.matrix.kappa<-rbind(design.matrix.kappa.1,design.matrix.kappa.2,design.matrix.kappa.3)
    
    design.vector.beta<-kronecker(X=matrix(rep(1,T),ncol=1),Y=beta.new)
    design.matrix.c<-kronecker(X=matrix(c(1,0,-1,0,1,-1),ncol=2),Y=design.vector.beta)
    
    design.matrix.A<-cbind(design.matrix.alpha,design.matrix.kappa,design.matrix.c)
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    alpha.new<-fit$coef[1:A]
    k.new<-fit$coef[(A+1):(A+T-1)]
    k.new<-c(-sum(k.new),k.new)
    c.new<-fit$coef[-(1:(A+T-1))]
    
    design.matrix.beta1<-kronecker(X=(k.new+c.new[1]),Y=diag(rep(1,A)))
    design.matrix.beta1<-design.matrix.beta1-design.matrix.beta1[,1]
    design.matrix.beta1<-design.matrix.beta1[,-1]
    design.matrix.beta2<-kronecker(X=(k.new+c.new[2]),Y=diag(rep(1,A)))
    design.matrix.beta2<-design.matrix.beta2-design.matrix.beta2[,1]
    design.matrix.beta2<-design.matrix.beta2[,-1]
    design.matrix.beta3<-kronecker(X=(k.new-c.new[1]-c.new[2]),Y=diag(rep(1,A)))
    design.matrix.beta3<-design.matrix.beta3-design.matrix.beta3[,1]
    design.matrix.beta3<-design.matrix.beta3[,-1]
    design.matrix.beta<-rbind(design.matrix.beta1,design.matrix.beta2,design.matrix.beta3)
    design.matrix.B<-cbind(design.matrix.alpha,design.matrix.beta)
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+c(kronecker(X=(k.new+c.new[1]),Y=c(1,rep(0,A-1))),kronecker(X=(k.new+c.new[2]),Y=c(1,rep(0,A-1))),kronecker(X=(k.new-c.new[1]-c.new[2]),Y=c(1,rep(0,A-1))))),family=poisson)
    alpha.new<-fit2$coef[1:A]
    beta.new<-fit2$coef[(A+1):(2*A-1)]
    beta.new<-c(1-sum(beta.new),beta.new)
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha=alpha.new,beta=beta.new,kappa=k.new,c.new=c.new,deviance=Deviance,model_bic=model_bic)
}


########################
##Model 2A2
########################

iteration.methodB.LC.poisson_M2A2<-function(deaths_combine,exposure_combine,m,beta1.initial,beta2.initial,beta3.initial,A,T){
  beta1.new<-beta1.initial
  beta2.new<-beta2.initial
  beta3.new<-beta3.initial
  alpha.new<-NULL
  k.new<-NULL
  Deviance<-vector(length=m)
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
  design.matrix.alpha<-kronecker(X=matrix(rep(1,3),ncol=1),Y=design.matrix.alpha)
  model_bic<-NULL
  
  for (i in 1:m){
    design.matrix.kappa1<-kronecker(X=diag(rep(1,T)),Y=beta1.new)
    design.matrix.kappa1<-design.matrix.kappa1-design.matrix.kappa1[,1]
    design.matrix.kappa1<-design.matrix.kappa1[,-1]
    design.matrix.kappa2<-kronecker(X=diag(rep(1,T)),Y=beta2.new)
    design.matrix.kappa2<-design.matrix.kappa2-design.matrix.kappa2[,1]
    design.matrix.kappa2<-design.matrix.kappa2[,-1]
    design.matrix.kappa3<-kronecker(X=diag(rep(1,T)),Y=beta3.new)
    design.matrix.kappa3<-design.matrix.kappa3-design.matrix.kappa3[,1]
    design.matrix.kappa3<-design.matrix.kappa3[,-1]
    design.matrix.kappa<-rbind(design.matrix.kappa1,design.matrix.kappa2,design.matrix.kappa3)
    design.matrix.A<-cbind(design.matrix.alpha,design.matrix.kappa)
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    alpha.new<-fit$coef[1:A]
    k.new<-fit$coef[(A+1):(A+T-1)]
    k.new<-c(-sum(k.new),k.new)
    
    #constraint:sum(beta1)=sum(beta2)=sum(beta3)=1 [does not work!]
    #design.matrix.beta<-kronecker(X=k.new,Y=diag(rep(1,A)))
    #design.matrix.beta<-design.matrix.beta-design.matrix.beta[,1]
    #design.matrix.beta<-design.matrix.beta[,-1]
    #design.matrix.beta<-kronecker(X=diag(rep(1,3)),Y=design.matrix.beta)
    #design.matrix.B<-cbind(design.matrix.alpha,design.matrix.beta)
    
    #constraint:sum(beta1)=1 [does not work!]
    #design.matrix.beta1<-kronecker(X=k.new,Y=diag(rep(1,A)))
    #design.matrix.beta1<-design.matrix.beta1-design.matrix.beta1[,1]
    #design.matrix.beta1<-design.matrix.beta1[,-1]
    #design.matrix.beta1<-kronecker(X=matrix(c(1,0,0),ncol=1),Y=design.matrix.beta1)
    #design.matrix.beta2<-kronecker(X=matrix(c(0,1,0),ncol=1),Y=kronecker(X=k.new,Y=diag(rep(1,A))))
    #design.matrix.beta3<-kronecker(X=matrix(c(0,0,1),ncol=1),Y=kronecker(X=k.new,Y=diag(rep(1,A))))
    #design.matrix.B<-cbind(design.matrix.alpha,design.matrix.beta1,design.matrix.beta2,design.matrix.beta3)
    
    #constraint:sum(beta1,beta2,beta3)=1 [works!]
    design.matrix.beta<-kronecker(X=k.new,Y=diag(rep(1,A)))
    design.matrix.beta<-kronecker(X=diag(rep(1,3)),Y=design.matrix.beta)
    design.matrix.beta<-design.matrix.beta-design.matrix.beta[,1]
    design.matrix.beta<-design.matrix.beta[,-1]
    design.matrix.B<-cbind(design.matrix.alpha,design.matrix.beta)
    
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+kronecker(X=matrix(c(1,0,0),ncol=1),Y=c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))))),family=poisson)
    alpha.new<-fit2$coef[1:A]
    beta1.new<-fit2$coef[(A+1):(2*A-1)]
    beta2.new<-fit2$coef[(2*A):(3*A-1)]
    beta3.new<-fit2$coef[(3*A):(4*A-1)]
    beta1.new<-c(1-sum(beta1.new,beta2.new,beta3.new),beta1.new)
    
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha=alpha.new,beta1=beta1.new,beta2=beta2.new,beta3=beta3.new,kappa=k.new,deviance=Deviance,model_bic=model_bic)
}


########################
##Model 2Y2
########################

iteration.methodB.LC.poisson_M2Y2<-function(deaths_combine,exposure_combine,m,beta.initial,A,T){
  beta.new<-beta.initial
  alpha.new<-NULL
  k.new<-NULL
  Deviance<-vector(length=m)
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
  design.matrix.alpha<-kronecker(X=matrix(rep(1,3),ncol=1),Y=design.matrix.alpha)
  model_bic<-NULL
  
  for (i in 1:m){
    design.matrix.kappa<-kronecker(X=diag(rep(1,T)),Y=beta.new)
    design.matrix.kappa<-kronecker(X=diag(rep(1,3)),Y=design.matrix.kappa)
    design.matrix.kappa<-design.matrix.kappa-design.matrix.kappa[,1]
    design.matrix.kappa<-design.matrix.kappa[,-1]
    design.matrix.A<-cbind(design.matrix.alpha,design.matrix.kappa)
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    alpha.new<-fit$coef[1:A]
    k1.new<-fit$coef[(A+1):(A+T-1)]
    k2.new<-fit$coef[(A+T):(A+2*T-1)]
    k3.new<-fit$coef[(A+2*T):(A+3*T-1)]
    k1.new<-c(-sum(k1.new,k2.new,k3.new),k1.new)
    
    design.matrix.beta1<-kronecker(X=k1.new,Y=diag(rep(1,A)))
    design.matrix.beta1<-design.matrix.beta1-design.matrix.beta1[,1]
    design.matrix.beta1<-design.matrix.beta1[,-1]
    design.matrix.beta2<-kronecker(X=k2.new,Y=diag(rep(1,A)))
    design.matrix.beta2<-design.matrix.beta2-design.matrix.beta2[,1]
    design.matrix.beta2<-design.matrix.beta2[,-1]
    design.matrix.beta3<-kronecker(X=k3.new,Y=diag(rep(1,A)))
    design.matrix.beta3<-design.matrix.beta3-design.matrix.beta3[,1]
    design.matrix.beta3<-design.matrix.beta3[,-1]
    design.matrix.beta<-rbind(design.matrix.beta1,design.matrix.beta2,design.matrix.beta3)
    design.matrix.B<-cbind(design.matrix.alpha,design.matrix.beta)
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+c(kronecker(X=k1.new,Y=c(1,rep(0,A-1))),kronecker(X=k2.new,Y=c(1,rep(0,A-1))),kronecker(X=k3.new,Y=c(1,rep(0,A-1))))),family=poisson)
    alpha.new<-fit2$coef[1:A]
    beta.new<-fit2$coef[(A+1):(2*A-1)]
    beta.new<-c(1-sum(beta.new),beta.new)
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha=alpha.new,beta=beta.new,kappa1=k1.new,kappa2=k2.new,kappa3=k3.new,deviance=Deviance,model_bic=model_bic)
}


########################
##Model Li and Lee 2005
########################

iteration.methodB.LC.poisson_lilee<-function(deaths_combine,exposure_combine,m,beta1.initial,beta2.initial,beta3.initial,Beta.initial,A,T){
  beta1.new<-beta1.initial
  beta2.new<-beta2.initial
  beta3.new<-beta3.initial
  Beta.new<-Beta.initial
  alpha.new<-NULL
  k.new<-NULL
  Deviance<-vector(length=m)
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
  design.matrix.alpha<-kronecker(X=diag(rep(1,3)),Y=design.matrix.alpha)
  model_bic<-NULL
  
  
  for (i in 1:m){
    design.matrix.kt1<-kronecker(X=diag(rep(1,T)),Y=beta1.new)
    design.matrix.kt1<-design.matrix.kt1-design.matrix.kt1[,1]
    design.matrix.kt1<-design.matrix.kt1[,-1]
    design.matrix.kt1<-kronecker(matrix(c(1,0,0),ncol=1),design.matrix.kt1)
    
    design.matrix.kt2<-kronecker(X=diag(rep(1,T)),Y=beta2.new)
    design.matrix.kt2<-design.matrix.kt2-design.matrix.kt2[,1]
    design.matrix.kt2<-design.matrix.kt2[,-1]
    design.matrix.kt2<-kronecker(matrix(c(0,1,0),ncol=1),design.matrix.kt2)
    
    #design.matrix.kt3<-kronecker(X=diag(rep(1,T)),Y=beta3.new)
    #design.matrix.kt3<-design.matrix.kt3-design.matrix.kt3[,1]
    #design.matrix.kt3<-design.matrix.kt3[,-1]
    #design.matrix.kt3<-kronecker(matrix(c(0,0,1),ncol=1),design.matrix.kt3)
    #we need to make kt3 zero for estimability 
    
    design.matrix.Kt<-kronecker(X=diag(rep(1,T)),Y=Beta.new)
    design.matrix.Kt<-design.matrix.Kt-design.matrix.Kt[,1]
    design.matrix.Kt<-design.matrix.Kt[,-1]
    design.matrix.Kt<-rbind(design.matrix.Kt,design.matrix.Kt,design.matrix.Kt)
    
    #design.matrix.A<-cbind(design.matrix.alpha,design.matrix.kt1,design.matrix.kt2,design.matrix.kt3,design.matrix.Kt)
    design.matrix.A<-cbind(design.matrix.alpha,design.matrix.kt1,design.matrix.kt2,design.matrix.Kt)
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    alpha1.new<-fit$coef[1:(A)]
    alpha2.new<-fit$coef[(A+1):(2*A)]
    alpha3.new<-fit$coef[(2*A+1):(3*A)]
    k1.new<-fit$coef[(3*A+1):(3*A+T-1)]
    k1.new<-c(-sum(k1.new),k1.new)
    k2.new<-fit$coef[(3*A+T):(3*A+2*T-2)]
    k2.new<-c(-sum(k2.new),k2.new)
    #k3.new<-fit$coef[(3*A+2*T-1):(3*A+3*T-3)]
    #k3.new<-c(-sum(k3.new),k3.new)
    k3.new<-rep(0,T)
    K.new<-fit$coef[-(1:(3*A+2*T-2))]
    K.new<-c(-sum(K.new),K.new)
    
    design.matrix.beta1<-kronecker(X=k1.new,Y=diag(rep(1,A)))
    design.matrix.beta1<-design.matrix.beta1-design.matrix.beta1[,1]
    design.matrix.beta1<-design.matrix.beta1[,-1]
    design.matrix.beta1<-kronecker(X=matrix(c(1,0,0),ncol=1),Y=design.matrix.beta1)
    
    design.matrix.beta2<-kronecker(X=k2.new,Y=diag(rep(1,A)))
    design.matrix.beta2<-design.matrix.beta2-design.matrix.beta2[,1]
    design.matrix.beta2<-design.matrix.beta2[,-1]
    design.matrix.beta2<-kronecker(X=matrix(c(0,1,0),ncol=1),Y=design.matrix.beta2)
    
    #design.matrix.beta3<-kronecker(X=k3.new,Y=diag(rep(1,A)))
    #design.matrix.beta3<-design.matrix.beta3-design.matrix.beta3[,1]
    #design.matrix.beta3<-design.matrix.beta3[,-1]
    #design.matrix.beta3<-kronecker(X=matrix(c(0,0,1),ncol=1),Y=design.matrix.beta3)
    
    design.matrix.Beta<-kronecker(X=K.new,Y=diag(rep(1,A)))
    design.matrix.Beta<-design.matrix.Beta-design.matrix.Beta[,1]
    design.matrix.Beta<-design.matrix.Beta[,-1]
    design.matrix.Beta<-rbind(design.matrix.Beta,design.matrix.Beta,design.matrix.Beta)
    
    #design.matrix.B<-cbind(design.matrix.alpha,design.matrix.beta1,design.matrix.beta2,design.matrix.beta3,design.matrix.Beta)
    design.matrix.B<-cbind(design.matrix.alpha,design.matrix.beta1,design.matrix.beta2,design.matrix.Beta)
    
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    offset.vector.k1<-kronecker(X=matrix(c(1,0,0),ncol=1),Y=kronecker(X=k1.new,Y=c(1,rep(0,A-1))))
    offset.vector.k2<-kronecker(X=matrix(c(0,1,0),ncol=1),Y=kronecker(X=k2.new,Y=c(1,rep(0,A-1))))
    offset.vector.k3<-kronecker(X=matrix(c(0,0,1),ncol=1),Y=kronecker(X=k3.new,Y=c(1,rep(0,A-1))))
    offset.vector.K<-1*rep(c(kronecker(X=K.new,Y=c(1,rep(0,A-1)))),3)
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+offset.vector.k1+offset.vector.k2+offset.vector.k3+offset.vector.K),family=poisson)
    alpha1.new<-fit2$coef[1:(A)]
    alpha2.new<-fit2$coef[(A+1):(2*A)]
    alpha3.new<-fit2$coef[(2*A+1):(3*A)]
    beta1.new<-fit2$coef[(3*A+1):(4*A-1)]
    beta1.new<-c(1-sum(beta1.new),beta1.new)
    beta2.new<-fit2$coef[(4*A):(5*A-2)]
    beta2.new<-c(1-sum(beta2.new),beta2.new)
    #beta3.new<-fit2$coef[(5*A-1):(6*A-3)]
    #beta3.new<-c(1-sum(beta3.new),beta3.new)
    Beta.new<-fit2$coef[-(1:(5*A-2))]
    Beta.new<-c(1-sum(Beta.new),Beta.new)
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha1=alpha1.new,alpha2=alpha2.new,alpha3=alpha3.new,beta1=beta1.new,beta2=beta2.new,beta3=beta3.new,Beta=Beta.new,kappa1=k1.new,kappa2=k2.new,kappa3=k3.new,Kappa=K.new,deviance=Deviance,model_bic=model_bic)
}


########################
##Model Li and Lee 2005 (modified,with shared alpha across products)
########################

iteration.methodB.LC.poisson_lilee_modified<-function(deaths_combine,exposure_combine,m,beta1.initial,beta2.initial,beta3.initial,Beta.initial,A,T){
  beta1.new<-beta1.initial
  beta2.new<-beta2.initial
  beta3.new<-beta3.initial
  Beta.new<-Beta.initial
  alpha.new<-NULL
  k.new<-NULL
  Deviance<-vector(length=m)
  deaths.vec<-as.vector(deaths_combine)
  exposure.vec<-as.vector(exposure_combine)
  design.matrix.alpha<-kronecker(X=matrix(rep(1,T),ncol=1),Y=diag(rep(1,A)))
  design.matrix.alpha<-kronecker(X=matrix(rep(1,3),ncol=1),Y=design.matrix.alpha)
  model_bic<-NULL
  
  
  for (i in 1:m){
    design.matrix.kt1<-kronecker(X=diag(rep(1,T)),Y=beta1.new)
    design.matrix.kt1<-design.matrix.kt1-design.matrix.kt1[,1]
    design.matrix.kt1<-design.matrix.kt1[,-1]
    design.matrix.kt1<-kronecker(matrix(c(1,0,0),ncol=1),design.matrix.kt1)
    
    design.matrix.kt2<-kronecker(X=diag(rep(1,T)),Y=beta2.new)
    design.matrix.kt2<-design.matrix.kt2-design.matrix.kt2[,1]
    design.matrix.kt2<-design.matrix.kt2[,-1]
    design.matrix.kt2<-kronecker(matrix(c(0,1,0),ncol=1),design.matrix.kt2)
    
    #design.matrix.kt3<-kronecker(X=diag(rep(1,T)),Y=beta3.new)
    #design.matrix.kt3<-design.matrix.kt3-design.matrix.kt3[,1]
    #design.matrix.kt3<-design.matrix.kt3[,-1]
    #design.matrix.kt3<-kronecker(matrix(c(0,0,1),ncol=1),design.matrix.kt3)
    
    design.matrix.Kt<-kronecker(X=diag(rep(1,T)),Y=Beta.new)
    design.matrix.Kt<-design.matrix.Kt-design.matrix.Kt[,1]
    design.matrix.Kt<-design.matrix.Kt[,-1]
    design.matrix.Kt<-rbind(design.matrix.Kt,design.matrix.Kt,design.matrix.Kt)
    
    #design.matrix.A<-cbind(design.matrix.alpha,design.matrix.kt1,design.matrix.kt2,design.matrix.kt3,design.matrix.Kt)
    design.matrix.A<-cbind(design.matrix.alpha,design.matrix.kt1,design.matrix.kt2,design.matrix.Kt)
    fit<-glm.fit(x=design.matrix.A,y=deaths.vec,family=poisson(),offset=log(exposure.vec))
    alpha.new<-fit$coef[1:(A)]
    k1.new<-fit$coef[(A+1):(A+T-1)]
    k1.new<-c(-sum(k1.new),k1.new)
    k2.new<-fit$coef[(A+T):(A+2*T-2)]
    k2.new<-c(-sum(k2.new),k2.new)
    #k3.new<-fit$coef[(A+2*T-1):(A+3*T-3)]
    #k3.new<-c(-sum(k3.new),k3.new)
    k3.new<-rep(0,T)
    K.new<-fit$coef[-(1:(A+2*T-2))]
    K.new<-c(-sum(K.new),K.new)
    
    design.matrix.beta1<-kronecker(X=k1.new,Y=diag(rep(1,A)))
    design.matrix.beta1<-design.matrix.beta1-design.matrix.beta1[,1]
    design.matrix.beta1<-design.matrix.beta1[,-1]
    design.matrix.beta1<-kronecker(X=matrix(c(1,0,0),ncol=1),Y=design.matrix.beta1)
    
    design.matrix.beta2<-kronecker(X=k2.new,Y=diag(rep(1,A)))
    design.matrix.beta2<-design.matrix.beta2-design.matrix.beta2[,1]
    design.matrix.beta2<-design.matrix.beta2[,-1]
    design.matrix.beta2<-kronecker(X=matrix(c(0,1,0),ncol=1),Y=design.matrix.beta2)
    
    #design.matrix.beta3<-kronecker(X=k3.new,Y=diag(rep(1,A)))
    #design.matrix.beta3<-design.matrix.beta3-design.matrix.beta3[,1]
    #design.matrix.beta3<-design.matrix.beta3[,-1]
    #design.matrix.beta3<-kronecker(X=matrix(c(0,0,1),ncol=1),Y=design.matrix.beta3)
    
    design.matrix.Beta<-kronecker(X=K.new,Y=diag(rep(1,A)))
    design.matrix.Beta<-design.matrix.Beta-design.matrix.Beta[,1]
    design.matrix.Beta<-design.matrix.Beta[,-1]
    design.matrix.Beta<-rbind(design.matrix.Beta,design.matrix.Beta,design.matrix.Beta)
    
    #design.matrix.B<-cbind(design.matrix.alpha,design.matrix.beta1,design.matrix.beta2,design.matrix.beta3,design.matrix.Beta)
    design.matrix.B<-cbind(design.matrix.alpha,design.matrix.beta1,design.matrix.beta2,design.matrix.Beta)
    #fit2<-glm.fit(x=design.matrix.B,y=deaths.vec,family=poisson(),offset=(log(exposure.vec)+1*rep(c(kronecker(X=k.new,Y=c(1,rep(0,A-1)))),3)))
    #replaced with below so we can find the BIC
    offset.vector.k1<-kronecker(X=matrix(c(1,0,0),ncol=1),Y=kronecker(X=k1.new,Y=c(1,rep(0,A-1))))
    offset.vector.k2<-kronecker(X=matrix(c(0,1,0),ncol=1),Y=kronecker(X=k2.new,Y=c(1,rep(0,A-1))))
    offset.vector.k3<-kronecker(X=matrix(c(0,0,1),ncol=1),Y=kronecker(X=k3.new,Y=c(1,rep(0,A-1))))
    offset.vector.K<-1*rep(c(kronecker(X=K.new,Y=c(1,rep(0,A-1)))),3)
    fit2<-glm(deaths.vec~design.matrix.B-1+offset(log(exposure.vec)+offset.vector.k1+offset.vector.k2+offset.vector.k3+offset.vector.K),family=poisson)
    alpha.new<-fit2$coef[1:(A)]
    beta1.new<-fit2$coef[(A+1):(2*A-1)]
    beta1.new<-c(1-sum(beta1.new),beta1.new)
    beta2.new<-fit2$coef[(2*A):(3*A-2)]
    beta2.new<-c(1-sum(beta2.new),beta2.new)
    #beta3.new<-fit2$coef[(3*A-1):(4*A-3)]
    #beta3.new<-c(1-sum(beta3.new),beta3.new)
    #Beta.new<-fit2$coef[-(1:(4*A-3))]
    Beta.new<-fit2$coef[-(1:(3*A-2))]
    Beta.new<-c(1-sum(Beta.new),Beta.new)
    Deviance[i]<-fit2$deviance
    model_bic[i]<-BIC(fit2)
  }
  list(alpha=alpha.new,beta1=beta1.new,beta2=beta2.new,beta3=beta3.new,Beta=Beta.new,kappa1=k1.new,kappa2=k2.new,kappa3=k3.new,Kappa=K.new,deviance=Deviance,model_bic=model_bic)
}




