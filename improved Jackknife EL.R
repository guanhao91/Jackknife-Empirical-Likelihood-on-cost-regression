#Improved method with coverage probability using J.EL method
library(splines)
library(survival)
BN=numeric(1)
r=log(2)/5##median survival time is 5 months
L=10
#if we set median time is 10 months, then L=20
#r=log(2)/10##median survival time is 5 months
#L=20
a0=410
b0=570
q=qchisq(0.95, 2, lower.tail = T, log.p = F)
#generate function
I<-function(x,y){ifelse(x<=y,1,0)}

sim2<-function (m,BN,n,cp){
  r=log(2)/m
  Count=0
  T<-numeric(length=n)
  M0=numeric(n)
  M=numeric(n)
  d=numeric(n)
  MR=numeric(n)
  B=numeric(n)
  X=numeric(n)
  delta=numeric(n)
  Zeta=numeric(n)
  D=numeric(n)
  
  
  
  for (j in 1:BN){
    
    T<-rexp(n,rate=r)#survival time
    Z=runif(n,0,10)
    e=rlnorm(n, meanlog = log(50), sdlog = 0.245)
    ep=matrix(rlnorm(n*10,log(10),0.245),nrow=n,ncol=10,byrow=T)
    epp=rlnorm(n,log(10),0.632)
    
    b=matrix(nrow=n,ncol=10)
    
    for(i in 1:n){
      d[i]=20*Z[i]+epp[i]
      M0[i]=50+50*Z[i]+e[i]  
      Sj=0
      for (j in 1:10){
        b[i,j]=50+100*Z[i]+ep[i,j]
        Sj=Sj+b[i,j]*(max(min(St[i],j)-(j-1),0))
      }
      M[i]=M0[i]+Sj+d[i]*I(St[i],L)
    }
    
    
    
    #real median cost values with parameter a0 and b0
    for(i in 1:n){
      MR[i]=a0 +b0*Z[i]
    }
    
    #simulate censored time
    C=runif(n,0,cp)
    
    for (i in 1:n){
      B[i]<-(I(M[i],a0+b0*Z[i])-1/2)*Z[i]
    }
    
    
    
    for (i in 1:n){
      X[i]=min(T[i],C[i])
      delta[i]=I(T[i],C[i])##death is 1
    }
    S.surv=survfit(Surv(X,delta)~1)
    ##KM estimator for censored time K.surv
    K.surv=survfit(Surv(X,1-delta)~1)
    
    ##generate estimated survival
    S.est<-stepfun(S.surv$time, c(1, S.surv$surv)) ##survival of any given time
    K.est<-stepfun(K.surv$time, c(1, K.surv$surv))
    
    
    #1.need to estimate the cost in subinterval (tj-1,min(tj-1,u))
    #define a fixed number of sub functions as 5
    J=5
    Mij<-function(i,j,u){
      if(u>(j-1)){
        return(b[i,j]*(min(u,j)-(j-1)))
      }
      else{
        return( 0 )
      }
    }
    
    #Gfun compared with B fun
    Gfun<-function(u){
      Expt=0
      for(i in 1:n) {
        Expt=Expt+B[i]*I(u,T[i])
      }
      return(Expt/n/S.est(u))
    }
    
    #Gfun2 compared with M fun
    Gfun2<-function(j,u){
      Expt=0
      for(i in 1:n) {
        Expt=Expt+delta[i]*Mij(i,j,u)*I(u,T[i])/K.est(T[i])
      }
      return(Expt/n/S.est(u))
    }
    
    K.haz<-stepfun(summary(K.surv)$time,
                   c(summary(K.surv)$n.event/summary(K.surv)$n.risk,1))
    
    gfun<-function(u){
      ifelse(S.est(u)==0 | K.est(u)==0, return(0), 
             return((B[i]-Gfun(u))/K.est(u)))
    }
    gfun2<-function(u){
      ifelse(S.est(u)==0 | K.est(u)==0, return(0), 
             return((Mij(i,j,u)-Gfun2(j,u))/K.est(u)))
    }
    
    
    #2. matrix cov(Yi, Wi) with column j as following;
    Covij<-function(i,j){
      sumpart<-0
      for(ic in 1:n){
        sumpart<-sumpart+delta[ic]/K.est(T[ic])*(B[ic]-Gfun(X[i]))*
          (Mij(i,j,X[i])-Gfun2(j,X[i]))*I(X[i],T[ic])
      }
      return(1/n*1/(K.est(X[i])^2)*(1/n*S.est(X[i]))*sumpart*I(delta[i],0))
      
    }
    #defind Cov matrix which has dimension 2*J
    #Cov for observation i will be as following
    Cov<-matrix(0,ncol=J,nrow=2)
    for(ic in 1:2){
      for (j in 1:J){
        Cov[ic,j]=Covij(i,j)[ic,1]
      }
    }
    #3. var(Wi) has (j,l) element as;
    Varijl<-function(i,j,l){
      sumpart<-0
      for(ic in 1:n){
        sumpart<-sumpart+delta[ic]/K.est(T[ic])*(Mij(X[i])-Gfun2(j,X[i]))*
          (Mij(ic,l,X[i])-Gfun2(l,X[i]))*I(X[i],T[ic])
      }
      return(1/n*1/(K.est(X[i])^2)*(1/n*S.est(X[i]))*sumpart*I(delta[i],0))
    }
    #defind Var matrix which has dimension J*J
    Var<-matrix(0,ncol=J,nrow=J)
    for(ic in 1:J){
      for (j in 1:J){
        for (l in 1:j){
          Var[ic,j]=Varijl(i,j,l)[ic,1]
        } 
      }
    }
    #4. small gamma(sg=cov/var)
    
    gammaopt<-Cov%*%solve(Var)
    
    Mic2<-function(u) {
      return(gfun(u)*K.haz(u)*I(u,X[i]))
    }
    
    Mic22<-function(u){
      return(gfun2(u)*K.haz(u)*I(u,X[i]))
    }
    
    #calculate influence function Zeta
    for(i in 1:n){
      Zeta[i]=B[i]-gfun(X[i])*I(delta[i],0)+
        integrate(Mic2, lower = 0, upper = L,subdivisions=2000)$value
    }
    
    
    #5. D[i]=zeta[i]+sum
    #calculate improved influence function D
    for(i in 1:n){
      sum=0
      for(j in 1:J){
        sum=sum+gfun2(X[i])*I(delta[i],0)+
          integrate(Mic22, lower = 0, upper = L,subdivisions=2000)$value
      }
      D[i]=zeta[i]+gammaopt*sum
    }
    
    # we have D[i] from improved method, then we can work on Jackknife part.
    Tn=sum(D)/n
    Tnj=numeric(n)
    for(i in 1:n){
      Tnj[i]=(sum(D)-D[i])/(n-1)
      J[i]=n*Tn-(n-1)*Tnj[i]
    }
    
    #total/n
    Lambda=(sum(J^2))^(-1)*sum(J)
    #solution 
    IFEL=2*sum(log(1+Lambda*J))
    Count=Count+I(IFEL,q)
  }
  return (Count/BN)
}
