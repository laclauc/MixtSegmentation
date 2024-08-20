#########################
# analyse
# MAJ 17/02/2024 vincent.brault@univ-grenoble-alpes.fr
#########################

rm(list=ls())
Niter<-100
K=3
p=4
source("R/Class_MixSegTime.R")
library(wavelets)
library(reshape2)

mu_res<-matrix(0,nrow=0,ncol=5)
T_max = 32 ## A day with 32 hours
time = (1:T_max)/T_max

for (n in c(100,1000)){
  for (d in c(50,100)){
    t = 1:d
    LK=matrix(NA,K,d)
    for (k in 1:K){
      LK[k,]=sort(rep(1:(k+1),length=d)) ### Segmentation of days
    }
    TL_ref<-sapply(1:K,function(k){which(!duplicated(LK[k,]))-1})
    Cluster<-unlist(sapply(1:K,function(k){rep(k,p*length(TL_ref[[k]]))}))
    for (alpha in c(0.1,0.2,1)){
      ##### Ref
      X = array(NA, dim = c(K,d,T_max))
      for (k in 1:K){
        for (j in 1:d){
          X[k,j,] = (-1)^(k) * alpha *cos(2*pi*time*LK[k,j])
        }
      }
      Xproj = array(NA, dim = c(K,d,p))
      for (k in 1:K){
        for (j in 1:d){
          Xproj[k,j,] = dwt(X[k,j,], filter="haar", n.levels = 3)@V[[3]]
        }
      }
      mu_star<-lapply(1:K,function(k){
        t(sapply(1:length(TL_ref[[k]]),function(l){
          Xproj[k,TL_ref[[k]][l]+1,1:p]
        }))
      })

      #####
      # Path
      Nom_com<-paste0("Plan_t_fixe_dgp/Plan_n_",n,"_d_",d,"_alpha_",10*alpha,"/")
      if (!dir.exists(Nom_com)){
        dir.create(Nom_com)
      }
      #####
      # Mu_est
      dist_mu=matrix(NA,nrow=Niter,ncol=length(Cluster))
      for (it in 1:Niter){
        Nom_fichier<-paste0(Nom_com,"Essai_",it,".Rdata")
        if (file.exists(Nom_fichier)){
          load(Nom_fichier)
          if (length(res_Dyn@TL)==K){
            dist_mu[it,]<-unlist(sapply(1:K,function(k){
              c(res_Dyn@mu[[k]]-mu_star[[k]])
            }))
          }
        }
      }
      N_part<-sum(apply(!is.na(dist_mu),MARGIN = 1,all))
      mu_res<-rbind(mu_res,cbind(as.matrix(melt(t(dist_mu[1:N_part,])))[,-(1:2)],
                                 rep(n,N_part*length(Cluster)),
                                 rep(d,N_part*length(Cluster)),
                                 rep(alpha,N_part*length(Cluster)),
                                 rep(Cluster,N_part)))
    }
  }
}

#######
# Display estimation of TL
#######
library(ggplot2)
mu_res<-data.frame(mu_res)
names(mu_res)=c("dist","n","d","alpha","Cluster")
M<-max(abs(mu_res$dist))
mu_res$alpha<-as.factor(mu_res$alpha)
mu_res$Cluster<-as.factor(mu_res$Cluster)
for (n in c(100,1000)){
  for (d in c(50,100)){
    dat<-mu_res[(mu_res$n==n)&(mu_res$d==d),]
    p<-ggplot(dat,aes(x=alpha,y=dist,fill=alpha,col=Cluster))
    pdf(paste0("Plan_t_fixe_dgp/Im/Mu_n_",n,"_d_",d,".pdf"),width=7,height=6)
    print(p+geom_boxplot()+
            labs(title=paste0("n=",n," and d=",d))+scale_y_continuous(limits = c(-M,M)))
    dev.off()
  }
}
