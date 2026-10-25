#initial values
simu.inits <- function(){
  omega=Glob_mu=Glob_t=matrix(NA,3,2)
  g=matrix(1,nrow=40,ncol=2)
  muofclust=tauofc=array(NA,c(40,consts$dim.space,(consts$clmax+3)))
  #
  omega[1,1]=omega[2,1]=1/3
  omega[3,1]=1/3
  omega[1,2]=omega[2,2]=1/3
  omega[3,2]=1/3
  #
  #
  g[1:20,1]=sample(c(rep(1,7),rep(2,7),rep(3,6)),replace = FALSE)
  g[1:20,2]=sample(c(rep(1,7),rep(2,7),rep(3,6)),replace = FALSE)
  #
  g[21:40,1]=sample(c(rep(1,7),rep(2,7),rep(3,6)),replace = FALSE)
  g[21:40,2]=sample(c(rep(1,7),rep(2,7),rep(3,6)),replace = FALSE)
  #
  #
  Glob_mu[1,1]=rnorm(1,0,sqrt(1/glob_empir_Mpre[1,1])) 
  while(
    Glob_mu[1,1] < -3.5*sqrt(1/glob_empir_Mpre[1,1]) | Glob_mu[1,1] > 3.5*sqrt(1/glob_empir_Mpre[1,1])) 
  {
    Glob_mu[1,1]=rnorm(1,0,sqrt(1/glob_empir_Mpre[1,1])) 
  }
  Glob_t[1,1]=runif(1,0.01,2*glob_empir_sd[1,1])
  #
  Glob_mu[1,2]=rnorm(1,0,sqrt(1/glob_empir_Mpre[1,2])) 
  while(
    Glob_mu[1,2] < -3.5*sqrt(1/glob_empir_Mpre[1,2]) | Glob_mu[1,2] > 3.5*sqrt(1/glob_empir_Mpre[1,2]))
  {
    Glob_mu[1,2]=rnorm(1,0,sqrt(1/glob_empir_Mpre[1,2])) 
  }    
  Glob_t[1,2]=runif(1,0.01,2*glob_empir_sd[1,2])
  #
  Glob_mu[2,1]=rnorm(1,empir_mean[1,1],sqrt(1/empir_Mpre[1,1]))
  Glob_t[2,1]=runif(1,0.01,2*empir_sd[1,1])
  Glob_mu[2,2]=rnorm(1,empir_mean[1,2],sqrt(1/empir_Mpre[1,2]))
  Glob_t[2,2]=runif(1,0.01,2*empir_sd[1,2])
  Glob_mu[3,1]=rnorm(1,empir_mean[2,1],sqrt(1/empir_Mpre[2,1]))
  Glob_t[3,1]=runif(1,0.01,2*empir_sd[2,1])
  Glob_mu[3,2]=rnorm(1,empir_mean[2,2],sqrt(1/empir_Mpre[2,2]))
  Glob_t[3,2]=runif(1,0.01,2*empir_sd[2,2])
  #   
  v=array(NA,c((consts$clmax-1),2,consts$dim.space))
  for(i in 1:(consts$clmax-1)){
    for(l in 1:2){
      for(k in 1:consts$dim.space){
        v[i,l,k]= rbeta(1, 1, 1)
      }
    }
  }
  #    
  nu=array(NA,c(consts$U,consts$N_curv,consts$dim.space))
  for(u in 1:(consts$U)){
    for(i in 1:consts$N_curv){
      for(k in 1:consts$dim.space){
        nu[u,i,k]=sample(c(20,19,18,17,16), size = 1, replace = TRUE)
      }
    }
  }
  #
  for(u in 1:(consts$U)){
    for(k in 1:(consts$dim.space)){  
      for(j in 4:(consts$clmax+3)){
        muofclust[u,k,(clmax+3+4-j)]=rnorm(1,0,sqrt(1/consts$empirSingle_Mpre[(consts$Glob[u,k]-1),k]))
      } 
    }
  } 
  #
  for(u in 1:(consts$U)){
    for(k in 1:(consts$dim.space)){  
      for(j in 1:3){
        muofclust[u,k,j]=Glob_mu[j,k]
      } 
    }
  }
  #
  for(u in 1:(consts$U)){
    for(k in 1:(consts$dim.space)){  
      for(j in 4:(consts$clmax+3)){
        tauofc[u,k,(clmax+3+4-j)]=runif(1,0.01,2*empirSingle_sd[(Glob[u,k]-1),k])
      } 
    }
  } 
  #
  for(u in 1:(consts$U)){
    for(k in 1:(consts$dim.space)){  
      for(j in 1:3){
        tauofc[u,k,j]=Glob_t[j,k]
      } 
    }
  }    
  #
  tau_eps=NULL
  tau_eps=rgamma(1,1,6)
  #
  #
  list(omega=omega,Glob_mu=Glob_mu,
       Glob_t=Glob_t,v=v,nu=nu,
       g=g,muofclust=muofclust,tauofc=tauofc,
       tau_eps=tau_eps)
}