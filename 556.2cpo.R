library(coda)
library(HDInterval)
library(mcmcse)

##
q=c(32+1,32+4,32+7,32+10,32+13,32+16,32+19,32+22)
p=c(33+1,33+4,33+7,33+10,33+13,33+16,33+19,33+22)
o=c(34+1,34+4,34+7,34+10,34+13,34+16,34+19,34+22)

### traceplot 
par(mfrow=c(3,2))
traceplot(nS[,33:60],main="Traceplot")
traceplot(nS[,33],main="Traceplot of pi1")

## acf auto corelation 
par(mfrow=c(2,3))
for(j in 33:ncol(S)){
  acf(nS[,j],lag.max = 50)
}

par(mfrow=c(1,1))
acf(nS[,33],lag.max = 50,main="ACF of pi1")

## NSE
est=matrix(0,nrow=3,ncol=l)
nse=matrix(0,nrow=3,ncol=l)
Spi1=nS[,q]
Spi2=nS[,p]
Sgma=nS[,o]
for(i in 1:l){
  est[,i]=c(mcse(Spi1[,i])$est,mcse(Spi2[,i])$est,mcse(Sgma[,i])$est)
  nse[,i]=c(mcse(Spi1[,i])$se,mcse(Spi2[,i])$se,mcse(Sgma[,i])$se)
}

## Gweke test
pvalue=2*(1-pnorm(abs(geweke.diag(nS,frac1=0.5,frac2=0.5)$z)))
Pval=rbind(pvalue[q],pvalue[p],pvalue[o])

## Effective Sample Size
EFF=rbind(effectiveSize(nS)[q],effectiveSize(nS)[p],effectiveSize(nS)[o])
nrow(nS)

round(nse,4)
round(Pval,4)
round(EFF,0)

##Diagnostic
Dig=matrix(0,nrow=3*3,ncol=l)
for(i in 1:3){
  Dig[((i-1)*3+1):((i-1)*3+3),]=rbind(nse[i,],Pval[i,],EFF[i,])
}

## posterior density & posterior mean & sd & HPD
for(i in 1:l){
  plot(density(nS[,q[i]]),xlim=c(0.01,1),main="Density Estimator of pi1")
}
for(i in 1:l){
  plot(density(nS[,p[i]]),xlim=c(0.01,1),main="Density Estimator of pi2")
}
for(i in 1:l){
  plot(density(nS[,o[i]]),xlim=c(0.01,1),main="Density Estimator of gma")
}
plot(density(nS[,q[2]]),xlim=c(0.01,1),main="Density Estimator of pi1")
plot(nS[,33])

## Mean Sd HPD
A=cbind(apply(nS,2,mean)[q],apply(nS,2,sd)[q],t(hdi(nS)[,q]))
B=cbind(apply(nS,2,mean)[p],apply(nS,2,sd)[p],t(hdi(nS)[,q]))
C=cbind(apply(nS,2,mean)[o],apply(nS,2,sd)[o],t(hdi(nS)[,o]))
round(A,4)
round(B,4)
round(C,4)














##define function
P_1=function(pi,gma){
  (gma*p1*pi[1] + gma*(1-p1)*pi[2]) / ( gma*p1*pi[1] + gma*(1-p1)*pi[2]  +  gma*p1*(1-pi[1]) + gma*(1-p1)*(1-pi[2]) )
}
P_2=function(pi,gma){
  ((1-gma)*pi[1]) / ((1-gma)*pi[1] + (1-gma)*(1-pi[1]))
}
P_3=function(pi,gma){
  (gma*p2*pi[1] + gma*(1-p2)*pi[2]) / ( gma*p2*pi[1] + gma*(1-p2)*pi[2]  +  gma*p2*(1-pi[1]) + gma*(1-p2)*(1-pi[2]) )
}
P_4=function(pi,gma){
  ((1-gma)*pi[1]) / ((1-gma)*pi[1] + (1-gma)*(1-pi[1]))
}
P_5=function(pi,gma){
  pi[2]
}



## CPO
for(i in 1:l){
  if(i==1){
Spi=nS[,c(32+(i-1)*3+1,32+(i-1)*3+2)]
Sgma=nS[,32+(i-1)*3+3]
  }else{
    Spi=cbind(Spi,nS[,c(32+(i-1)*3+1,32+(i-1)*3+2)])
    Sgma=cbind(Sgma,nS[,32+(i-1)*3+3])
  }
}

CPO=matrix(0,nrow=nrow(y),ncol=ncol(y))
inv=matrix(0,nrow=nrow(nS),ncol=l*ncol(y))

for(k in 1:ncol(y)){
for(r in 1:l){
for(j in 1:nrow(nS)){
  if(k==1){inv[j,8*(k-1)+r]=1/dbinom(y[r,k],n[r,k],P_1(Spi[j,((r-1)*2+1):((r-1)*2+2)],Sgma[j,r]))}
  if(k==2){inv[j,8*(k-1)+r]=1/dbinom(y[r,k],n[r,k],P_2(Spi[j,((r-1)*2+1):((r-1)*2+2)],Sgma[j,r]))}
  if(k==3){inv[j,8*(k-1)+r]=1/dbinom(y[r,k],n[r,k],P_3(Spi[j,((r-1)*2+1):((r-1)*2+2)],Sgma[j,r]))}
  if(k==4){inv[j,8*(k-1)+r]=1/dbinom(y[r,k],n[r,k],P_4(Spi[j,((r-1)*2+1):((r-1)*2+2)],Sgma[j,r]))}
  if(k==5){inv[j,8*(k-1)+r]=1/dbinom(y[r,k],n[r,k],Spi[j,(r-1)*2+2])}
  }
CPO[r,k]=(mean(inv[,8*(k-1)+r]))^(-1)
}
}
CPO
(1/CPO)>40
round(1/CPO,4)


## Posterior predictive p-value
H=ncol(y)
K=10

Y=matrix(0,nrow=l*nrow(nS),ncol=H*K)
E=matrix(0,nrow=l*nrow(nS),ncol=H*K)
V=matrix(0,nrow=l*nrow(nS),ncol=H*K)
Yk=matrix(0,nrow=l,ncol=ncol(y))
Ek=matrix(0,nrow=l,ncol=ncol(y))
Vk=matrix(0,nrow=l,ncol=ncol(y))
Trep=matrix(0,nrow=nrow(nS),ncol=K)
Tobs=matrix(0,nrow=nrow(nS),ncol=K)

for(j in 1:nrow(nS)){
  for(k in 1:K){
    for(i in 1:l){
      yi1=rbinom(1,n[i,1],P_1(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      yi2=rbinom(1,n[i,2],P_2(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      yi3=rbinom(1,n[i,3],P_3(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      yi4=rbinom(1,n[i,4],P_4(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      yi5=rbinom(1,n[i,5],P_5(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      Yk[i,]=c(yi1,yi2,yi3,yi4,yi5)
      
      Ei1=n[i,1]*P_1(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ei2=n[i,2]*P_2(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ei3=n[i,3]*P_3(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ei4=n[i,4]*P_4(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ei5=n[i,5]*P_5(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ek[i,]=c(Ei1,Ei2,Ei3,Ei4,Ei5)
      
      Vi1=n[i,1] * P_1(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_1(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      Vi2=n[i,2] * P_2(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_2(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      Vi3=n[i,3] * P_3(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_3(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      Vi4=n[i,4] * P_4(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_4(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      Vi5=n[i,5] * P_5(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_5(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      if(Vi1==0){Vi1=1}
      if(Vi2==0){Vi2=1}
      if(Vi3==0){Vi3=1}
      if(Vi4==0){Vi4=1}
      if(Vi5==0){Vi5=1}
      Vk[i,]=c(Vi1,Vi2,Vi3,Vi4,Vi5)
    }
    Trep[j,k]=sum((Yk-Ek)^2/Vk)
  }
}

for(j in 1:nrow(nS)){
  for(k in 1:K){
    for(i in 1:l){
      yi1=y[i,1]
      yi2=y[i,2]
      yi3=y[i,3]
      yi4=y[i,4]
      yi5=y[i,5]
      Yk[i,]=c(yi1,yi2,yi3,yi4,yi5)
      
      Ei1=n[i,1]*P_1(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ei2=n[i,2]*P_2(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ei3=n[i,3]*P_3(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ei4=n[i,4]*P_4(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ei5=n[i,5]*P_5(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i])
      Ek[i,]=c(Ei1,Ei2,Ei3,Ei4,Ei5)
      
      Vi1=n[i,1] * P_1(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_1(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      Vi2=n[i,2] * P_2(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_2(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      Vi3=n[i,3] * P_3(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_3(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      Vi4=n[i,4] * P_4(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_4(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      Vi5=n[i,5] * P_5(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]) * (1-P_5(Spi[j,((i-1)*2+1):((i-1)*2+2)],Sgma[j,i]))
      if(Vi1==0){Vi1=1}
      if(Vi2==0){Vi2=1}
      if(Vi3==0){Vi3=1}
      if(Vi4==0){Vi4=1}
      if(Vi5==0){Vi5=1}
      Vk[i,]=c(Vi1,Vi2,Vi3,Vi4,Vi5)
    }
    Tobs[j,k]=sum((Yk-Ek)^2/Vk)
  }
}


BBP=sum((Trep-Tobs)>=0)/(nrow(nS)*K)
BBP

## BPP in (0.1,0.9)
