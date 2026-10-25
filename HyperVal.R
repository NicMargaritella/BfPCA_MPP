### hyperprior values-------------------------------------
#                      
glob_empir_Mvar<-matrix(NA,1,2)
gl1=gl2=Mgl1=Mgl2=NULL
for(i in 1:1000){
  string=sample(1:40,40,replace = TRUE) 
  for(j in 1:40){
    gl1[((j*50)-49):(j*50)]=
      fPCA_info$empir_score[((string[j]*50)-49):(string[j]*50),1]
  }
  Mgl1[i]=mean(gl1)
  gl1=NULL
}
glob_empir_Mvar[1,1]<-var(Mgl1)
#
for(i in 1:1000){
  string=sample(1:40,40,replace = TRUE) 
  for(j in 1:40){ 
    gl2[((j*50)-49):(j*50)]=
      fPCA_info$empir_score[((string[j]*50)-49):(string[j]*50),2]
  }
  Mgl2[i]=mean(gl2)
  gl2=NULL
}
glob_empir_Mvar[1,2]<-var(Mgl2)
glob_empir_Msd=sqrt(glob_empir_Mvar)
glob_empir_Mpre=1/glob_empir_Mvar
#
glob_empir_sd<-matrix(NA,1,2) 
glob_empir_sd[1,1]<-sd(fPCA_info$empir_score[1:2000,1])
glob_empir_sd[1,2]<-sd(fPCA_info$empir_score[1:2000,2])
#  
empir_mean<-matrix(NA,2,2)
empir_mean[1,1]<-mean(fPCA_info$empir_score[1:1000,1])
empir_mean[1,2]<-mean(fPCA_info$empir_score[1:1000,2])
empir_mean[2,1]<-mean(fPCA_info$empir_score[1001:2000,1])
empir_mean[2,2]<-mean(fPCA_info$empir_score[1001:2000,2])
empir_Mvar<-matrix(NA,2,2)
#
e1=e2=f1=f2=NULL
Me1=Me2=Mf1=Mf2=NULL
for(i in 1:1000){
  string=sample(1:20,8,replace = TRUE) 
  for(j in 1:8){
    e1[((j*50)-49):(j*50)]=
      fPCA_info$empir_score[((string[j]*50)-49):(string[j]*50),1]
  }
  Me1[i]=mean(e1)
  e1=NULL
}
empir_Mvar[1,1]<-var(Me1)
#
for(i in 1:1000){
  string=sample(1:20,8,replace = TRUE) 
  for(j in 1:8){
    e2[((j*50)-49):(j*50)]=
      fPCA_info$empir_score[((string[j]*50)-49):(string[j]*50),2]
  }
  Me2[i]=mean(e2)
  e2=NULL
}
empir_Mvar[1,2]<-var(Me2)  
#
for(i in 1:1000){
  string=sample(21:40,8,replace = TRUE) 
  for(j in 1:8){
    f1[((j*50)-49):(j*50)]=
      fPCA_info$empir_score[((string[j]*50)-49):(string[j]*50),1]
  }
  Mf1[i]=mean(f1)
  f1=NULL
}
empir_Mvar[2,1]<-var(Mf1)
#
for(i in 1:1000){
  string=sample(21:40,8,replace = TRUE) 
  for(j in 1:8){
    f2[((j*50)-49):(j*50)]=
      fPCA_info$empir_score[((string[j]*50)-49):(string[j]*50),2]
  }
  Mf2[i]=mean(f2)
  f2=NULL
}
empir_Mvar[2,2]<-var(Mf2)  
#
empir_Mpre=1/empir_Mvar
# 
#
empir_sd<-matrix(NA,2,2)
empir_sd[1,1]<-sd(fPCA_info$empir_score[1:1000,1])
empir_sd[1,2]<-sd(fPCA_info$empir_score[1:1000,2])
empir_sd[2,1]<-sd(fPCA_info$empir_score[1001:2000,1])
empir_sd[2,2]<-sd(fPCA_info$empir_score[1001:2000,2])
#
#
empirSingle_Msd=matrix(NA,2,2)
empirSingle_Msd[1,1]=(range(fPCA_info$empir_score[1:1000,1])[2]-range(fPCA_info$empir_score[1:1000,1])[1])/2.5
empirSingle_Msd[1,2]=(range(fPCA_info$empir_score[1:1000,2])[2]-range(fPCA_info$empir_score[1:1000,2])[1])/2.5
empirSingle_Msd[2,1]=(range(fPCA_info$empir_score[1001:2000,1])[2]-range(fPCA_info$empir_score[1001:2000,1])[1])/2.5
empirSingle_Msd[2,2]=(range(fPCA_info$empir_score[1001:1000,2])[2]-range(fPCA_info$empir_score[1001:2000,2])[1])/2.5
#
empirSingle_Mpre=1/(empirSingle_Msd^2)
#
sds1=sds2=sds3=sds4=NULL
empirSingle_sd<- matrix(NA,2,2)
for(i in 1:20){
  sds1[i]=sd(fPCA_info$empir_score[((i*50)-49):(i*50),1])
}
empirSingle_sd[1,1]=mean(sds1)
#
for(i in 1:20){
  sds2[i]=sd(fPCA_info$empir_score[((i*50)-49):(i*50),2])
}
empirSingle_sd[1,2]=mean(sds2)
#
for(i in 1:20){
  lis=c(21:40)
  sds3[i]=sd(fPCA_info$empir_score[((lis[i]*50)-49):(lis[i]*50),1])
}
empirSingle_sd[2,1]=mean(sds3)
#
for(i in 1:20){
  lis=c(21:40)
  sds4[i]=sd(fPCA_info$empir_score[((lis[i]*50)-49):(lis[i]*50),2])
}
empirSingle_sd[2,2]=mean(sds4)
#