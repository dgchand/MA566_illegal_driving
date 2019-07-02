rm(list=ls())

## X1=c(1,0)=c(M,F)
## X2=c(1,0)=c(>24,<=24)
## X3=c(1,0)=c(>3,<=3)

data=read.csv("data.csv")
data2=read.csv("data2.csv")
Group1=data[1:32,]
Group2=data[33:64,]
Group3=data2[1:16,]

## Toss 
p1=1/3
p2=2/3

## i: domians, i=1,2,3,4,5,6,7,8
## s: groups, s=1,2,3
l=8
y=matrix(0,nrow=l,ncol=3)
colnames(y)=c("y1","y2","y3")
n=matrix(0,nrow=l,ncol=3)
for(i in 1:l){
  y[i,1]=Group1[(i-1)*4+1,6]+Group1[(i-1)*4+3,6]
  y[i,2]=Group2[(i-1)*4+1,6]+Group2[(i-1)*4+3,6]
  y[i,3]=Group3[(i-1)*2+1,5]
  n[i,1]=sum(Group1[((i-1)*4+1):((i-1)*4+4),6])
  n[i,2]=sum(Group2[((i-1)*4+1):((i-1)*4+4),6])
  n[i,3]=sum(Group3[((i-1)*2+1):((i-1)*2+2),5])
}

## introduce latent varibale zi11,zi12,zi13,zi14,zi21,zi22,zi23,zi24
z=matrix(1,nrow=l,ncol=12)
colnames(z)=c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12")

gma=rep(0.5,l)
pi=matrix(0.5,nrow=l,ncol=2)
mu=rep(0.5,3)
tau=1

##Sample
M=10000
##define function
G=100
grid=seq(0,1,length=G+1)
mgrid=c()
for(i in 1:G){
  mgrid[i]=(grid[i]+grid[i+1])/2
}

logh_phi=function(phi,pi,mu,gma){
  tau=phi/(1-phi)
  sum=(0)
  for(i in 1:l){
    sum=sum+tau * (mu[1]*log(pi[i,1]) + (1-mu[1])*log(1-pi[i,1]) + mu[2]*log(pi[i,2]) + (1-mu[2])*log(1-pi[i,2]) +
                     mu[3]*log(gma[i]) + (1-mu[3])*log(1-gma[i])) 
  }
  sum=sum-l* log( beta(mu[1]*tau,(1-mu[1])*tau) *beta(mu[2]*tau,(1-mu[2])*tau) * beta(mu[3]*tau,(1-mu[3])*tau))
  sum
}

f_phi=function(phi,pi,mu,gma){
  exp(logh_phi(phi,pi,mu,gma))/sum(exp(logh_phi(mgrid,pi,mu,gma)))
}

logh_mu=function(mu,tau,pi){
  mu*tau*log(sum(pi)) + (1-mu)*tau*log(sum(1-pi)) - l * log(beta(mu*tau,(1-mu)*tau))
}

f_mu=function(mu,tau,pi){
  exp(logh_mu(mu,tau,pi))/sum(exp(logh_mu(mgrid,tau,pi)))
}

###Gibbs sampler
S=matrix(0,nrow=M,ncol=124)

for(j in 1:M){
for(i in 1:l){
s1=gma[i]*p1*pi[i,1] + gma[i]*(1-p1)*pi[i,2] + (1-gma[i])*pi[i,1]
z[i,1:3]=rmultinom(1,y[i,1],c( gma[i]*p1*pi[i,1]/s1 , gma[i]*(1-p1)*pi[i,2]/s1 , (1-gma[i])*pi[i,1]/s1 ))
s2=gma[i]*p1*(1-pi[i,1]) + gma[i]*(1-p1)*(1-pi[i,2]) + (1-gma[i])*(1-pi[i,1])
z[i,4:6]=rmultinom(1,n[i,1]-y[i,1],c( gma[i]*p1*(1-pi[i,1])/s2 , gma[i]*(1-p1)*(1-pi[i,2])/s2 , (1-gma[i])*(1-pi[i,1])/s2 ))
s3=gma[i]*p2*pi[i,1] + gma[i]*(1-p2)*pi[i,2] + (1-gma[i])*pi[i,1]
z[i,7:9]=rmultinom(1,y[i,2],c( gma[i]*p2*pi[i,1]/s3 , gma[i]*(1-p2)*pi[i,2]/s3 , (1-gma[i])*pi[i,1]/s3 ))
s4=gma[i]*p2*(1-pi[i,1]) + gma[i]*(1-p2)*(1-pi[i,2]) + (1-gma[i])*(1-pi[i,1])
z[i,10:12]=rmultinom(1,n[i,2]-y[i,2],c( gma[i]*p2*(1-pi[i,1])/s4 , gma[i]*(1-p2)*(1-pi[i,2])/s4 , (1-gma[i])*(1-pi[i,1])/s4 ))

S[j,((i-1)*12+1):((i-1)*12+12)]=z[i,]

pi[i,1]=rbeta( 1, sum(z[i,c(1,3,7,9)])+mu[1]*tau , sum(z[i,c(4,6,10,12)])+(1-mu[1])*tau)
if(pi[i,1]==1){pi[i,1]=0.99999}
S[j,96+(i-1)*3+1]=pi[i,1]
pi[i,2]=rbeta( 1, sum(z[i,c(2,8)])+y[i,3]+mu[2]*tau , sum(z[i,c(5,11)])+n[i,3]-y[i,3]+(1-mu[2])*tau)
if(pi[i,2]==1){pi[i,1]=0.99999}
S[j,96+(i-1)*3+2]=pi[i,2]
gma[i]=rbeta( 1, sum(z[i,c(1,2,4,5,7,8,10,11)])+mu[3]*tau , sum(z[i,c(3,6,9,12)])+(1-mu[3])*tau )
if(gma[i]==1){gma[i]=0.99999}
S[j,96+(i-1)*3+3]=gma[i]
}

d_mu1=f_mu(mgrid,tau,pi[,1])
mu[1]=mean(sample(mgrid,size=M/10,replace=T,prob=d_mu1))
S[j,121]=mu[1]
d_mu2=f_mu(mgrid,tau,pi[i,2])
mu[2]=mean(sample(mgrid,size=M/10,replace=T,prob=d_mu2))
S[j,122]=mu[2]
d_mu3=f_mu(mgrid,tau,gma)
mu[3]=mean(sample(mgrid,size=M/10,replace=T,prob=d_mu3))
S[j,123]=mu[3]
d_phi=f_phi(mgrid,pi,mu,gma)
phi=mean(sample(mgrid,size=M/10,replace=T,prob=d_phi))
tau=phi/(1-phi)
S[j,124]=tau
}

nS=S[-(1:M/2),][seq(1,M-M/2,5),]
nS=as.mcmc(nS)

