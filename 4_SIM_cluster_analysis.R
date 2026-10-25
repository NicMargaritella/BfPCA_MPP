#
# Cluster analysis using ARI, VI and Cluster balls
#
#
library(pdfCluster)
library(coda)
library(mcclust)
library(mcclust.ext)
#
s=1 # set same value as in (3) 
#
#
# set constants and 
# R objects to store results
#
N_subj=40
N_curv=50
dim.space=2
N_subj_spc=4
dim.space1=1
ARI1<-NULL
ARI2<-NULL
clustSubSpec_label<-array(0,c(1,N_subj_spc,N_curv))
clustOv_label<-array(0,c(1,N_subj,dim.space))
#
ARI2_SubSpec1=ARI2_SubSpec2=matrix(NA,1,6000)
#6000=size of all retained iterations from all MCMC chains used
ARI2_SubSpec39=ARI2_SubSpec40=matrix(NA,1,6000)
ARI2_SubSpec1_mode=ARI2_SubSpec2_mode=NULL
ARI2_SubSpec39_mode=ARI2_SubSpec40_mode=NULL
#
VI_results1=VI_results2=vector(mode = "list", length =1)
CB_results1=CB_results2=vector(mode = "list", length =1)
#
VI_s1=VI_s2=VI_s39=VI_s40=vector(mode = "list", length =1)
CB_s1=CB_s2=CB_s39=CB_s40=vector(mode = "list", length =1)
#
#
load(sprintf("SIMU_%s",s))
# 
combinedchains <- mcmc.list(posterior$samples$chain1,
                              posterior$samples$chain2)
#
combinedchains=as.matrix(combinedchains)
#
# overall clustering for the 2 dimensions
ClustOv_true1<-c(rep(2,20),rep(4,20))
ClustOv_true2<-c(rep(3,2),rep(2,18),rep(4,18),rep(3,2))
#
### Clustering with Adjusted Rand Index (ARI) 
#
ARR1_Ov<-array(0,c(1,N_subj,dim.space,(dim(combinedchains)[1])))
a<-NULL
for(l in 1:N_subj){
  for(k in 1:dim.space){
  ARR1_Ov[1,l,k,]<-combinedchains[,sprintf("g[%s, %s]",l,k)]
  a<-which.max(prop.table(table(ARR1_Ov[1,l,k,])))
  clustOv_label[s,l,k]<-max(as.integer(names(a)),0)
    if(l>20 & clustOv_label[s,l,k]==2){
    clustOv_label[s,l,k]=4
    }
  }
}
#
ARI1<-adj.rand.index(ClustOv_true1,clustOv_label[1,,1])
ARI2<-adj.rand.index(ClustOv_true2,clustOv_label[1,,2])
#
#
## Clustering of subject specific trajectories with ARI
#
#true clusters  
Clust_SubSpec_true1<-c(rep(1,20),rep(2,30))
Clust_SubSpec_true2<-c(rep(1,20),rep(2,30))
#
ARR2_SubSpec<-array(0,c(1,N_subj_spc,N_curv,(dim(combinedchains)[1])))
a<-NULL
id_label=c(1,2,39,40)
for(l in 1:N_subj_spc){
  for(i in 1:N_curv){
  ARR2_SubSpec[1,l,i,]<-combinedchains[,sprintf("z[%s, %s, 2]",id_label[l], i)]
  a<-which.max(prop.table(table(ARR2_SubSpec[1,l,i,])))
  clustSubSpec_label[1,l,i]<-max(as.integer(names(a)),0)
  }
}
#
for(t in 1:dim(combinedchains)[1]){
ARI2_SubSpec1[1,t]= adj.rand.index(Clust_SubSpec_true1,ARR2_SubSpec[1,1,,t])  
ARI2_SubSpec2[1,t]= adj.rand.index(Clust_SubSpec_true1,ARR2_SubSpec[1,2,,t]) 
ARI2_SubSpec39[1,t]= adj.rand.index(Clust_SubSpec_true2,ARR2_SubSpec[1,3,,t])     
ARI2_SubSpec40[1,t]= adj.rand.index(Clust_SubSpec_true2,ARR2_SubSpec[1,4,,t])
}
#
ARI2_SubSpec1_mode=max(as.numeric(names
                      (which.max(table(round(ARI2_SubSpec1[1,],3))))),0)
ARI2_SubSpec2_mode=max(as.numeric(names(
                      which.max(table(round(ARI2_SubSpec2[1,],3))))),0)
ARI2_SubSpec39_mode=max(as.numeric(names(
                       which.max(table(round(ARI2_SubSpec39[1,],3))))),0)
ARI2_SubSpec40_mode=max(as.numeric(names(
                       which.max(table(round(ARI2_SubSpec40[1,],3))))),0)
#
#
### Clustering with Variation of Information
### and Credible Balls ---
#
# dim 1--- relabelling
ARR1_cu<-matrix(0,(dim(combinedchains)[1]),N_subj)
for(l in 1:N_subj){
ARR1_cu[,l]<-combinedchains[,sprintf("g[%s, 1]",l)]
  for(j in 1:dim(combinedchains)[1]){ 
    if(l>20 & ARR1_cu[j,l]==2){
    ARR1_cu[j,l]=4
    }
    if(ARR1_cu[j,l]==3){
      ARR1_cu[j,l]=3*l
    }
  }
}
#
for(g in 1:dim(combinedchains)[1]){
  config_g=ARR1_cu[g,]
  label_g=unique(config_g)
  for(k in 1:length(label_g)){
    ARR1_cu[g,config_g==label_g[k]]=k
  }
}
#
ppmat1=comp.psm(ARR1_cu)# post prob matrix
#
VI_results1=minVI(ppmat1,cls.draw=ARR1_cu,method="all",
      include.greedy = T,suppress.comment = F)
#
CB_results1=credibleball(VI_results1$cl[1,],
                              ARR1_cu, alpha = 0.05)
#
# dim 2--- relabelling
ARR2_cu<-matrix(0,(dim(combinedchains)[1]),N_subj)
for(l in 1:N_subj){
ARR2_cu[,l]<-combinedchains[,sprintf("g[%s, 2]",l)]
  for(j in 1:dim(combinedchains)[1]){ 
    if(l>20 & ARR2_cu[j,l]==2){
    ARR2_cu[j,l]=4
    }
    if(ARR2_cu[j,l]==3){
    ARR2_cu[j,l]=3*l
    }  
  }
}
#
for(g in 1:dim(combinedchains)[1]){
  config_g=ARR2_cu[g,]
  label_g=unique(config_g)
  for(k in 1:length(label_g)){
    ARR2_cu[g,config_g==label_g[k]]=k
  }
}
#
ppmat2=comp.psm(ARR2_cu)# post prob matrix
#
VI_results2=minVI(ppmat2,cls.draw=ARR2_cu,method="all",
                  include.greedy = T,suppress.comment = F)
#
CB_results2=credibleball(VI_results2$cl[1,],
                              ARR2_cu,alpha = 0.05)
#
#
### Clustering within subject-specific deviations using VI
#
# s1----------
s1ARR2=t(as.matrix(ARR2_SubSpec[1,1,,]))#
#
for(g in 1:dim(combinedchains)[1]){
config_g=s1ARR2[g,]
label_g=unique(config_g)
  for(k in 1:length(label_g)){
  s1ARR2[g,config_g==label_g[k]]=k
  }
}
ppmat_s1=comp.psm(s1ARR2)# post prob matrix
#
VI_s1=minVI(ppmat_s1,cls.draw=s1ARR2,method="all",
                 include.greedy = T,suppress.comment = F)
CB_s1=credibleball(VI_s1$cl[1,],
                        s1ARR2, alpha = 0.05)
#
rm(ppmat_s1,s1ARR2, config_g,label_g)
# s2----------
s2ARR2=t(as.matrix(ARR2_SubSpec[1,2,,]))# N=50 curves, M=6000 mcmc iterations
#
for(g in 1:dim(combinedchains)[1]){
config_g=s2ARR2[g,]
label_g=unique(config_g)
 for(k in 1:length(label_g)){
 s2ARR2[g,config_g==label_g[k]]=k
 }
}
ppmat_s2=comp.psm(s2ARR2)# post prob matrix
#
VI_s2=minVI(ppmat_s2,cls.draw=s2ARR2,method="all",
                 include.greedy = T,suppress.comment = F)
CB_s2=credibleball(VI_s2$cl[1,],
                        s2ARR2, alpha = 0.05)
#
rm(ppmat_s2,s2ARR2, config_g,label_g)
# s39----------
s39ARR2=t(as.matrix(ARR2_SubSpec[1,3,,]))# N=50 curves, M=6000 mcmc iterations
#
for(g in 1:dim(combinedchains)[1]){
  config_g=s39ARR2[g,]
  label_g=unique(config_g)
  for(k in 1:length(label_g)){
    s39ARR2[g,config_g==label_g[k]]=k
  }
}
ppmat_s39=comp.psm(s39ARR2)# post prob matrix
#
VI_s39=minVI(ppmat_s39,cls.draw=s39ARR2,method="all",
                  include.greedy = T,suppress.comment = F)
CB_s39=credibleball(VI_s39$cl[1,],
                         s39ARR2, alpha = 0.05)
#
rm(ppmat_s39,s39ARR2, config_g,label_g)
# s40----------
s40ARR2=t(as.matrix(ARR2_SubSpec[1,4,,]))# N=50 curves, M=6000 mcmc iterations
#
for(g in 1:dim(combinedchains)[1]){
  config_g=s40ARR2[g,]
  label_g=unique(config_g)
  for(k in 1:length(label_g)){
    s40ARR2[g,config_g==label_g[k]]=k
  }
}
ppmat_s40=comp.psm(s40ARR2)# post prob matrix
#
VI_s40=minVI(ppmat_s40,cls.draw=s40ARR2,method="all",
                  include.greedy = T,suppress.comment = F)
CB_s40=credibleball(VI_s40$cl[1,],
                         s40ARR2, alpha = 0.05)
#
rm(ppmat_s40,s40ARR2, config_g,label_g)

#
# VI classification performance - group clusters
VI_classif=VI2_classif=NULL

VI_classif=adj.rand.index(VI_results1$cl[1,],
                                             c(rep(1,20),rep(2,20)))
VI2_classif=adj.rand.index(VI_results2$cl[1,],
                                c(3,5,rep(2,18),rep(4,18),7,9))
#
# VI classification performance - subject-specific clusters
VI_s1_res=VI_s2_res=VI_s39_res=VI_s40_res=NULL

VI_s1_res=adj.rand.index(VI_s1$cl[1,],
                              Clust_SubSpec_true1)
VI_s2_res=adj.rand.index(VI_s2$cl[1,],
                              Clust_SubSpec_true1)
VI_s39_res=adj.rand.index(VI_s39$cl[1,],
                               Clust_SubSpec_true2)
VI_s40_res=adj.rand.index(VI_s40$cl[1,],
                               Clust_SubSpec_true2)
#
rm(ARR2_SubSpec,ARR1_Ov,ARR2_cu,
   posterior,ppmat1,ppmat2)
#
#
### clusters uncertainty with cluster balls
#
alpha=0.05
### D1 horizontal bound 
uncert_VI_h=NULL
  if(isTRUE(dim(CB_results1$c.horiz)[1]==(1-alpha)*dim(combinedchains)[1])==FALSE){
    uncert_VI_h= rbind(uncert_VI_h, aggregate(numdup ~.,data=transform(CB_results1$c.horiz,
                                             numdup=1,ID=s),length))
  }
uncert_VI_h[,43]=round(uncert_VI_h[,42]/((1-alpha)*dim(combinedchains)[1]),2)
uncert_VI_h=uncert_VI_h[uncert_VI_h[,42]>=1000,]
#
### D1 upper vertical bound
uncert_VI_u=NULL
  if(isTRUE(dim(CB_results1$c.uppervert)[1]==(1-alpha)*dim(combinedchains)[1])==FALSE){
    uncert_VI_u= rbind(uncert_VI_u, aggregate(numdup ~.,data=transform(CB_results1$c.uppervert,
                                              numdup=1,ID=s),length))
}

uncert_VI_u[,43]=round(uncert_VI_u[,42]/((1-alpha)*dim(combinedchains)[1]),2)
uncert_VI_u=uncert_VI_u[uncert_VI_u[,42]>=1000,]#
#
### D1 lower vertical bound
#
uncert_VI_l=NULL
  if(isTRUE(dim(CB_results1$c.lowervert)[1]==(1-alpha)*dim(combinedchains)[1])==FALSE){
    uncert_VI_l= rbind(uncert_VI_l, aggregate(numdup ~.,data=transform(CB_results1$c.lowervert,
                                             numdup=1,ID=s),length))
  }

uncert_VI_l[,43]=round(uncert_VI_l[,42]/((1-alpha)*dim(combinedchains)[1]),2)
uncert_VI_l=uncert_VI_l[uncert_VI_l[,42]>=1000,]
#
#
### D2 horizontal bound
uncert_VI2_h=NULL
  if(isTRUE(dim(CB_results2$c.horiz)[1]==(1-alpha)*dim(combinedchains)[1])==FALSE){
    uncert_VI2_h= rbind(uncert_VI2_h, aggregate(numdup ~.,data=transform(CB_results2$c.horiz,                                                numdup=1,ID=s),length))
  }

uncert_VI2_h[,43]=round(uncert_VI2_h[,42]/((1-alpha)*dim(combinedchains)[1]),2)
uncert_VI2_h=uncert_VI2_h[uncert_VI2_h[,42]>=1000,]
#
### D2 upper vertical bound
#
uncert_VI2_u=NULL
  if(isTRUE(dim(CB_results2$c.uppervert)[1]==(1-alpha)*dim(combinedchains)[1])==FALSE){
    uncert_VI2_u= rbind(uncert_VI2_u, aggregate(numdup ~.,data=transform(CB_results2$c.uppervert,                                               numdup=1,ID=s),length))
  }

uncert_VI2_u[,43]=round(uncert_VI2_u[,42]/((1-alpha)*dim(combinedchains)[1]),2)
uncert_VI2_u=uncert_VI2_u[uncert_VI2_u[,42]>=1000,]
#
### D2 lower vertical bound
#
uncert_VI2_l=NULL
  if(isTRUE(dim(CB_results2$c.lowervert)[1]==(1-alpha)*dim(combinedchains)[1])==FALSE){
    uncert_VI2_l= rbind(uncert_VI2_l, aggregate(numdup ~.,data=transform(CB_results2$c.lowervert,                                               numdup=1,ID=s),length))
  }

uncert_VI2_l[,43]=round(uncert_VI2_l[,42]/((1-alpha)*dim(combinedchains)[1]),2)
uncert_VI2_l=uncert_VI2_l[uncert_VI2_l[,42]>=1000,]
#
rm(list=setdiff(ls(),
  c("ARI1","ARI2","ARI2_SubSpec1_mode","ARI2_SubSpec2_mode","
    ARI2_SubSpec39_mode","ARI2_SubSpec40_mode","CB_results1","
    CB_results2","VI_results1","VI_results2","VI_classif","VI2_classif","
    uncert_VI_h","uncert_VI2_h","uncert_VI_u","uncert_VI2_u","
    uncert_VI_l","uncert_VI2_l","VI_s1_res","VI_s2_res","
    VI_s39_res","VI_s40_res")))
#
save(ARI1,ARI2,ARI2_SubSpec1_mode,ARI2_SubSpec2_mode,
     ARI2_SubSpec39_mode,ARI2_SubSpec40_mode,CB_results1,
     CB_results2,VI_results1,VI_results2,VI_classif,VI2_classif,
     uncert_VI_h,uncert_VI2_h,uncert_VI_u,uncert_VI2_u,
     uncert_VI_l,uncert_VI2_l,VI_s1_res,VI_s2_res,
     VI_s39_res,VI_s40_res,
     file="SIM_clusteringResults")
#