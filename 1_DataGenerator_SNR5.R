#Simulated trajectories for the case with SNR=5
#
load("Bases.RData")
#
library(mnormt)
#
# number of subjects
U=40
# number of trajectories
N_curv=50
# number of points per trajectory
Timep=150
# 
# noisy trajectories
Y_it<-as.list(replicate(U,matrix(0,N_curv,Timep),simplify=F))
# smooth trajectories
X_it<- as.list(replicate(U,matrix(0,N_curv,Timep),simplify=F))
#
### GROUP A ###
# Subjects 1&2------------
#
fscoresA1=fscoresA2=matrix(NA,nrow=2,ncol=N_curv)
#
for(s in 1:2){
# generate the scores for the first dimension
set.seed(2+s) 
g11<-as.vector(rnorm(20,7,2.5))
g12<-as.vector(rnorm(30,15,1))
g1<-as.vector(c(g11,g12))
#
# generate the scores for the second dimension
set.seed(3+s)
g2<-as.vector(rnorm(50,-6,1.5))
#
fscoresA=cbind(g1,g2)
fscoresA1[s,]<-g1
fscoresA2[s,]<-g2
rm(g1,g2,g11,g12)
#
# Obtain smooth trajectories
for(l in 1:N_curv){
  for(t in 1:Timep){ 
    X_it[[s]][l,t]<-(t(fscoresA[l,])%*%Bases[[1]][t,])
  }
}
#
# Add noise
#
means=rep(0,N_curv)
var= 0.15
Sigma_E=diag(x=var,N_curv,N_curv)
#
set.seed(2+s)
E_it<- t(rmnorm(n = Timep, mean = means, Sigma_E, sqrt=NULL))
#
# build the datasets
#
for(i in 1:N_curv){
  for(t in 1:Timep){
    Y_it[[s]][i,t]<-X_it[[s]][i,t]+E_it[i,t]
  }
} 
}
#
# Subjects 3 to 20----------
#
fscoresB1=fscoresB2=matrix(NA,nrow=18,ncol=N_curv)
#
for(s in 3:20){
#
set.seed(2+s)
g1<-rnorm(50,7,2.5)
# 
# generate the scores for the second dimension
set.seed(3+s)
g2<-rnorm(50,-6,1.5)
#
fscoresB=cbind(g1,g2)
fscoresB1[(s-2),]<-g1
fscoresB2[(s-2),]<-g2
rm(g1,g2)
#
# Obtain smooth trajectories
for(l in 1:N_curv){
  for(t in 1:Timep){ 
    X_it[[s]][l,t]<-(t(fscoresB[l,])%*%Bases[[1]][t,])
  }
}
#
# add noise
#
set.seed(3+s)
E_it<- t(rmnorm(n = Timep, mean = means, Sigma_E, sqrt=NULL))
#
# build the datasets
#
for(i in 1:N_curv){
  for(t in 1:Timep){
     Y_it[[s]][i,t]<-X_it[[s]][i,t]+E_it[i,t]
  }
} 
}
### GROUP B ###
# Subjects 21 to 38--------
#
fscoresC1=fscoresC2=matrix(NA,nrow=18,ncol=N_curv)
#
for(s in 21:38){
#  
# generate the scores for the first dimension
set.seed(2+s)
g1<-rnorm(50,-7,2.5)
# 
# generate the scores for the second dimension
set.seed(3+s)
g2<-rnorm(50,6,1.5)
#
fscoresC=cbind(g1,g2)
fscoresC1[(s-20),]<-g1
fscoresC2[(s-20),]<-g2
rm(g1,g2)
#
# Obtain smooth trajectories
for(l in 1:N_curv){
  for(t in 1:Timep){ 
    X_it[[s]][l,t]<-(t(fscoresC[l,])%*%Bases[[1]][t,])
  }
}
#
# add noise
#
set.seed(3+s)
E_it<- t(rmnorm(n = Timep, mean = means, Sigma_E, sqrt=NULL))
#
# build the datasets
#
for(i in 1:N_curv){
  for(t in 1:Timep){
    Y_it[[s]][i,t]<-X_it[[s]][i,t]+E_it[i,t]
  }
} 
}
#
# Subjects 39 & 40---------------
#
fscoresD1=fscoresD2=matrix(NA,nrow=2,ncol=N_curv)
for(s in 39:40){
#  
# generate the scores for the first dimension
set.seed(2+s) 
g11<-as.vector(rnorm(20,-7,2.5))
g12<-as.vector(rnorm(30,-15,1))
g1<-as.vector(c(g11,g12))
#    
# generate the scores for the second dimension
set.seed(3+s)
g2<-as.vector(rnorm(50,6,1.5))
# 
fscoresD=cbind(g1,g2)
fscoresD1[(s-38),]<-g1
fscoresD2[(s-38),]<-g2
rm(g1,g2,g11,g12)
#
# Obtain smooth trajectories
for(l in 1:N_curv){
  for(t in 1:Timep){ 
    X_it[[s]][l,t]<-(t(fscoresD[l,])%*%Bases[[1]][t,])
  }
}
#
# add noise
#
set.seed(2+s)
E_it<- t(rmnorm(n = Timep, mean = means, Sigma_E, sqrt=NULL))
#
# build the datasets
#
for(i in 1:N_curv){
  for(t in 1:Timep){
    Y_it[[s]][i,t]<-X_it[[s]][i,t]+E_it[i,t]
  }
} 
}
#
ScoresD2=rbind(fscoresA1,fscoresB1,fscoresC1,fscoresD1)
ScoresD1=rbind(fscoresA2,fscoresB2,fscoresC2,fscoresD2)
rm(l,s,t,fscoresA1,fscoresB1,fscoresC1,fscoresD1,
   fscoresA2,fscoresB2,fscoresC2,fscoresD2,
   fscoresA,fscoresB,fscoresC,fscoresD)
#
X_all=Reduce(rbind.data.frame, X_it)
#
#
save(Y_it,var,X_all,X_it,
     ScoresD1,ScoresD2,file="simu_data")
#