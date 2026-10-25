### Data smoothing and fPCA 
#
load("simu_data")
#
library(fda)
Y_all<-Reduce( rbind.data.frame, Y_it[1:40] )
#
N_curv=dim(Y_all)[1]
Timep=dim(Y_all)[2]
W_it=t(as.matrix(Y_all))# matrix N_obs*N_subj 
#
### --Data Smoothing--
#
fdarange<-c(1,Timep)
fdatime<-seq(1,Timep,by=1)# observation points
fdabasis<-create.bspline.basis(fdarange,nbasis=Timep+2,norder=4)
#
#  Lambda-tuning with GCV criterion
#
loglamb=seq(0,4,by=0.01)
nlamb=length(loglamb)
dfsave=rep(NA,nlamb)
gcvsave=rep(NA,nlamb)
  for(i in 1:nlamb){
      lambda=10^loglamb[i]
      fdParobj= fdPar(fdabasis,2,lambda)
      smoothlist= smooth.basis(fdatime,W_it,fdParobj)
      dfsave[i]=smoothlist$df
      gcvsave[i]=sum(smoothlist$gcv)
  }
#
log_gcv<-loglamb[which.min(gcvsave)]
lambda=10^log_gcv
#
fdParobj= fdPar(fdabasis,2,lambda)
curve.fit=smooth.basis(fdatime,W_it,fdParobj)
curve.fd=curve.fit$fd# this contains the smoothed curves
#
#### ---FPCA----
#
# this center the curves and applies fPCA
cWsmooth_pca<-pca.fd(curve.fd,nharm=min(N_curv-1,Timep-1),centerfns = FALSE)
dim.space<-min(N_curv-1,Timep-1)
Wc_it<-t(W_it)
#    
max_E<-min(c(min(which(cumsum(cWsmooth_pca$varprop)>0.95)),length(which
                 (cWsmooth_pca$varprop>0.025))))
rotated<-varmx.pca.fd(cWsmooth_pca,nharm=max_E)
max_E<-min(c(min(which(cumsum(rotated$varprop)>0.95)),
                   length(which(rotated$varprop>0.05))))
Ei<-eval.fd(fdatime,rotated$harmonics,returnMatrix = FALSE)
#
E<-as.matrix(Ei[,1:(max_E)])# eigenfunctions matrix  
eigval<-rotated$values[1:(max_E)]# relative eigenvalues
varprop<-rotated$varprop[1:(max_E)]
empir_score<-rotated$scores[,1:max_E]#
#
fPCA_info<-list(max_E=max_E, log_gcv=log_gcv, eigval=eigval, E=E,varprop=varprop,
                       empir_score=empir_score,Wc_it=Wc_it, Timep = Timep, N_curv = N_curv,
                       dim.space=dim.space)
#
#
save(fPCA_info,var,X_all,X_it,
       file=("simu_data_fpca"))