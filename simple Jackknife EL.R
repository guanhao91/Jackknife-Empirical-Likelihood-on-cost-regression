library(splines)
library(survival)
L=10 #time period frame,limit research time to 10
a0=410 #regression intercept
b0=570 #regression coefficient
q=qchisq(0.95, 2, lower.tail = T, log.p = F)#95% chisq quantile value with df=2
#generate function
I<-function(x,y){ifelse(x<=y,1,0)}
#m is the median parameter in months
###sim1(median survival, sample size, censoring parameter)####
sim2<-function (m,n,cp){
  
  r=log(2)/m #rate to compute survival time given median value
  Count=0 #count of less than chisq threshold
  St<-numeric(length=n)  #survival time
  M0=numeric(n)  #simulated initial cost
  M=numeric(n) #simulated total cost
  d=numeric(n) #terminal cost
  MR=numeric(n) #expected median regression cost given covariates
  B=numeric(n)
  X=numeric(n)
  delta=numeric(n)
  Zeta=numeric(n)
  
  #j is the number of count
  for (j in 1:1000){
    
    St<-rexp(n,rate=r) 
    Z=runif(n,0,10) #set covariates
    e=rlnorm(n, meanlog = log(50), sdlog = 0.245) #errors for initial cost
    ep=matrix(rlnorm(n*10,log(10),0.245),nrow=n,ncol=10,byrow=T) #errors for monthly cost
    epp=rlnorm(n,log(10),0.632) ##errors for terminal cost
    
    b=matrix(nrow=n,ncol=10) ##monthly cost for different people in different period
    
    for(i in 1:n){
      d[i]=20*Z[i]+epp[i] ##terminal cost
      M0[i]=50+50*Z[i]+e[i]  ##initial cost
      Sj=0 ##total monthly cost
      for (j in 1:10){
        b[i,j]=50+100*Z[i]+ep[i,j]
        Sj=Sj+b[i,j]*(max(min(St[i],j)-(j-1),0))
      }
      M[i]=M0[i]+Sj+d[i]*I(St[i],L) ##simulated total cost for each patient
    }
    
    
    for(i in 1:n){
      MR[i]=a0 +b0*Z[i] ##median regression cost
    }
    
    C=runif(n,0,cp) ##simulated censoring time with given CP(censoring parameter)
    
    
    
    #B function values
    for (i in 1:n){
      B[i]<-(I(M[i],a0+b0*Z[i])-1/2)*Z[i]
    }
    
    
    
    for (i in 1:n){
      X[i]=min(St[i],C[i])
      delta[i]=I(St[i],C[i])##death indicator(1 death, 0 censored)
    }
    
    
    S.surv=survfit(Surv(X,delta)~1)
    ##KM estimator for censored time K.surv
    K.surv=survfit(Surv(X,1-delta)~1)
    
    ##generate estimated survival
    S.est<-stepfun(S.surv$time, c(1, S.surv$surv)) ##survival of any given time
    K.est<-stepfun(K.surv$time, c(1, K.surv$surv)) 
    
    #define G function
    Gfun<-function(u){
      Expt=0
      for(i in 1:n) {
        Expt=Expt+B[i]*I(u,X[i])
      }
      return(Expt/n/S.est(u))
    }
    
    #hazard function of the censoring distribution
    K.haz<-stepfun(summary(K.surv)$time, c(summary(K.surv)$n.event/summary(K.surv)$n.risk,1))
    
    #in case there is 0 survival probability. gfun is B[i]-Gfun(u))/K.est(u)
    gfun<-function(u){
      ifelse(S.est(u)==0 | K.est(u)==0, return(0), return((B[i]-Gfun(u))/K.est(u)))
    }
    
    #derivative of matingale process function 
    Mic2<-function(u) {
      return(gfun(u)*K.haz(u)*I(u,X[i]))
    }
    
    #influnence function
    for(i in 1:n){
      Zeta[i]=B[i]-gfun(X[i])*I(delta[i],0)+integrate(Mic2, lower = 0, upper = L,subdivisions=2000)$value
    }
    
    # we have Zeta[i] from simple method, then we can work on Jackknife part.
    Tn=sum(Zeta)/n
    Tnj=numeric(n)
    for(i in 1:n){
      Tnj[i]=(sum(Zeta)-Zeta[i])/(n-1)
      J[i]=n*Tn-(n-1)*Tnj[i]
    }
    #total/n
    Lambda=(sum(J^2))^(-1)*sum(J)
    #solution 
    IFEL=2*sum(log(1+Lambda*J))
    Count=Count+I(IFEL,q)
  }
  return (Count/1000)
}
