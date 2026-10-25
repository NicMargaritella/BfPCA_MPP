# Simulation study - running model with NIMBLE
#
#NOTE: with the current MCMC settings it will run
# for a while (see paper for comp time).
#
library(nimble)
library(igraph)
library(coda)
library(R6)
#
#
load("simu_data_fpca")
s=1 # set a seed for replicability  e.g. 1
#
### prepare data list 
# subjects
U=40
# trajectories
N_curv=50
# data points
Timep<-150
#
Wc<-fPCA_info$Wc_it
# transform data matrix into array
Wc1<-array(0,dim=c(U,N_curv,Timep))
for(i in 0:(U-1)){
Wc1[(i+1),,]<-Wc[((i*N_curv)+1):((i+1)*N_curv),]
}
Wc<-Wc1
rm(Wc1)
#
E<-fPCA_info$E #extract eigenfunctions
dim.space<-dim(E)[2]
#
# define the fixed groups
Glob=c(rep(2,15),rep(3,25))
Glob=cbind(Glob,Glob)
#
### hyperparameters
source("HyperVal.R")
#
clmax=20
C=diag(3)
#
dataList = list(Wc=Wc,E=E,C=C,delta=c(1,1,1))
#
consts <- list(N_curv = N_curv, U=U, Timep = Timep,
                 dim.space=dim.space,clmax=clmax,Glob=Glob,
                 glob_empir_Msd=glob_empir_Msd,
                 glob_empir_Mpre=glob_empir_Mpre,
                 glob_empir_sd=glob_empir_sd,
                 empir_mean=empir_mean,
                 empir_Mpre=empir_Mpre,
                 empir_sd=empir_sd,
                 empirSingle_Mpre=empirSingle_Mpre,
                 empirSingle_sd=empirSingle_sd)
#
### initial values
source("InitialVal.R")
myinits <- simu.inits()
multi_SimuInits <- list(simu.inits(),simu.inits())# 
seed=s# for replicability  
#
# Import the "sort" function in the NIMBLE environment
Rsort<-nimbleRcall(
  function(x=double(1)){}, Rfun = "sort",
  returnType = double(1), where = .GlobalEnv)
#
#
### Model statement----------------------------------- 
Multimodel <- nimbleCode({
#
for (u in 1:U) { 
  for (i in 1:N_curv) {
    for (t in 1:Timep) {
    Wc[u,i,t] ~ dnorm(m[u,i,t],tau_eps) 
    m[u,i,t]<- sum(xi[u,i,1:2]* E[t,1:2])
    }# close t
  }# close i
#
  for(k in 1:dim.space) {
    for(i in 1:N_curv) {
    xi[u,i,k]~ dnorm(muofclust[u,k,z[u,i,k]],tauofclust[u,k,z[u,i,k]])
    z[u,i,k]<- C[g[u,k],1]+ C[g[u,k],2]*Glob[u,k]+C[g[u,k],3]*(3+nu[u,i,k])
    nu[u,i,k]~dcat(pii[1:20,(Glob[u,k]-1),k]) 
    }# close i
#
  g[u,k]~ dcat(omega[1:3,k])
#
### prior clusters mean and var### 
# global:
    for(f in 1:3){ 
    muofclust[u,k,f]<-Glob_mu[f,k]     
    tauofclust[u,k,f]<-Glob_tau[f,k] 
    }#close f
#
# local:
    for(j in 4:(clmax+3)){
    muofclust[u,k,(20+3+4-j)]~dnorm(0, empirSingle_Mpre[(Glob[u,k]-1),k])
    tauofc[u,k,(20+3+4-j)]~dunif(0.01,2*empirSingle_sd[(Glob[u,k]-1),k])
    tauofclust[u,k,(20+3+4-j)]<-1/(tauofc[u,k,(20+3+4-j)]^2)
    }#close j
  }#close k
}#close u
#
for(k in 1:dim.space){
omega[1:3,k]~ddirch(delta[1:3])
#
Glob_mu[1,k]~T(dnorm(0,glob_empir_Mpre[1,k]),-3.5*glob_empir_Msd[1,k], 3.5*glob_empir_Msd[1,k])
Glob_t[1,k]~dunif(0.01,2*glob_empir_sd[1,k])
Glob_tau[1,k]<-1/(Glob_t[1,k]^2)
#
  for(f in 2:3){
  Glob_mu[f,k]~dnorm(empir_mean[(f-1),k],empir_Mpre[(f-1),k])
  Glob_t[f,k]~ dunif(0.01,2*empir_sd[(f-1),k])
  Glob_tau[f,k]<-1/(Glob_t[f,k]^2)
  }# close f
#
# stick-breaking prior for the local clusters
  for(l in 1:2){
    for(j in 1:(clmax-1)){
    v[j,l,k] ~ T(dbeta(1, 1),0.001,0.999)
    }# close j
  pi[1:clmax,l,k]<- stick_breaking(v[1:(clmax-1),l,k])
  pii[1:clmax,l,k]<- Rsort(pi[1:clmax,l,k])
  }# close l
}#close k
#
tau_eps~dgamma(0.001,0.001)# 
})
#---------------------------------------------------
#
## to build the model
SIM_1 <- nimbleModel(code = Multimodel, name = "multi_mod",
                    constants = consts,
                    data = dataList, inits<-myinits)
# To compile the model
SIM_2 <- compileNimble(SIM_1)
#
# set up the monitored random quantities.
SIM_3 <- configureMCMC(SIM_1, monitors = c("tau_eps","omega","g",
                                          "muofclust","tauofclust",
                                          "Glob_mu","Glob_tau","pii",
                                          "z"), enableWAIC = FALSE,
                                          print = FALSE) 
# build the MCMC algorithm
SIM_4 <- buildMCMC(SIM_3)
# compile the MCMC chain 
SIM_5 <- compileNimble(SIM_4, project = SIM_1)
#
set.seed(seed)
posterior <- runMCMC(SIM_5, niter = 19000, thin=3, nburnin=10000, 
                     summary = TRUE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = multi_SimuInits) 
#
save(posterior, myinits,multi_SimuInits,
    file=sprintf("SIMU_%s",seed))
#
#