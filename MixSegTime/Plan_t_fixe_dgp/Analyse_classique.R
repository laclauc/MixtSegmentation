#########################
# analyse
#
#########################

rm(list=ls())
Niter<-100
K=3
source("R/Class_MixSegTime.R")

library(bikm1)

TL_res<-matrix(0,nrow=0,ncol=4)
ARI_res<-matrix(0,nrow=0,ncol=4)
CE_res<-matrix(0,nrow=0,ncol=4)

for (n in c(100,1000)){
  for (d in c(50,100)){
    LK=matrix(NA,K,d)
    for (k in 1:K){
      LK[k,]=sort(rep(1:(k+1),length=d)) ### Segmentation of days
    }
    TL_ref<-sapply(1:K,function(k){which(!duplicated(LK[k,]))[-1]-1})
    for (alpha in c(0.1,0.2,1)){
      Nom_com<-paste0("Plan_t_fixe_dgp/Plan_n_",n,"_d_",d,"_alpha_",10*alpha,"/")
      if (!dir.exists(Nom_com)){
        dir.create(Nom_com)
      }
      dist_PG=rep(NA,Niter)
      ARI=rep(NA,Niter)
      CE=rep(NA,Niter)
      for (it in 1:Niter){
        Nom_fichier<-paste0(Nom_com,"Essai_",it,".Rdata")
        if (file.exists(Nom_fichier)){
          load(Nom_fichier)
          if (length(res_Dyn@TL)==K){
            dist_PG[it]<-max(sapply(1:K,function(k){
              max(abs(TL_ref[[k]]-res_Dyn@TL[[k]][2:(length(TL_ref[[k]])+1)]))}))
            ARI[it]<-ARI(v = Z,res_Dyn@z)$ari
            CE[it]<-NCE_simple(Z,res_Dyn@z)
          }
        }
      }
      N_part<-sum(!is.na(dist_PG))
      TL_res<-rbind(TL_res,
                    matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),dist_PG[!is.na(dist_PG)]/d),
                           nrow=N_part))
      ARI_res<-rbind(ARI_res,
                     matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),ARI[!is.na(dist_PG)]),
                            nrow=N_part))
      CE_res<-rbind(CE_res,
                    matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),CE[!is.na(dist_PG)]),
                           nrow=N_part))
    }
  }
}

#######
# Display estimation of TL
#######
library(ggplot2)
TL_res<-data.frame(TL_res)
names(TL_res)=c("n","d","alpha","dist")
M<-max(TL_res$dist)
for (n in c(100,1000)){
  for (d in c(50,100)){
    dat<-TL_res[(TL_res$n==n)&(TL_res$d==d),]
    dat$alpha<-as.factor(dat$alpha)
    p<-ggplot(dat,aes(x=alpha,y=dist,fill=alpha))
    pdf(paste0("Plan_t_fixe_dgp/Im/TL_n_",n,"_d_",d,".pdf"),width=8,height=7)
    print(p+geom_boxplot()+
            theme(legend.position = "none")+scale_y_continuous(limits = c(0,M)))
    dev.off()
  }
}

#########
# Table ARI
#########
ARI_tab<-matrix("",4,3)
it<-1
ARI_res<-data.frame(ARI_res)
names(ARI_res)=c("n","d","alpha","value")
for (n in 10^(2:3)){
  for (d in c(50,100)){
    for (j in 1:3){
      alpha<-c(0.1,0.2,1)[j]
      dat<-ARI_res[(ARI_res$n==n)&(ARI_res$d==d)&(ARI_res$alpha==alpha),]
      ARI_tab[it,j]<-paste0(signif(mean(dat$value),2)," (",signif(sd(dat$value),2),")")
    }
    it<-it+1
  }
}

row.names(ARI_tab)<-c(sapply(c(100,1000),function(n){
  sapply(c(50,100),function(d){
    paste0("\n (",n,",",d,")")
  })
}))
row.names(ARI_tab)[1]<-paste0("\\hline ",row.names(ARI_tab)[1])
colnames(ARI_tab)<-as.character(c(0.1,0.2,1))
write.table(ARI_tab,file = "Plan_t_fixe_dgp/Tab/ARI.tex",quote = FALSE,sep = "&",eol = "\\\\")

#########
# Table ARI
#########
CE_tab<-matrix("",4,3)
it<-1
CE_res<-data.frame(CE_res)
names(CE_res)=c("n","d","alpha","value")
for (n in 10^(2:3)){
  for (d in c(50,100)){
    for (j in 1:3){
      alpha<-c(0.1,0.2,1)[j]
      dat<-CE_res[(CE_res$n==n)&(CE_res$d==d)&(CE_res$alpha==alpha),]
      if (all(dat$value<10^(-15))){
        CE_tab[it,j]<-"0 (0)"
      }else{
        CE_tab[it,j]<-paste0(signif(mean(dat$value),2)," (",signif(sd(dat$value),2),")")
      }
    }
    it<-it+1
  }
}

row.names(CE_tab)<-c(sapply(c(100,1000),function(n){
  sapply(c(50,100),function(d){
    paste0("\n (",n,",",d,")")
  })
}))
row.names(CE_tab)[1]<-paste0("\\hline ",row.names(CE_tab)[1])
colnames(CE_tab)<-as.character(c(0.1,0.2,1))
write.table(CE_tab,file = "Plan_t_fixe_dgp/Tab/CE.tex",quote = FALSE,sep = "&",eol = "\\\\")

##############
# TCL Mu
##############

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
M<-max(abs(mu_res$dist),na.rm = TRUE)
mu_res$alpha<-as.factor(mu_res$alpha)
mu_res$Cluster<-as.factor(mu_res$Cluster)
for (n in c(100,1000)){
  for (d in c(50,100)){
    dat<-mu_res[(mu_res$n==n)&(mu_res$d==d),]
    p<-ggplot(dat,aes(x=alpha,y=dist,fill=Cluster))
    pdf(paste0("Plan_t_fixe_dgp/Im/Mu_n_",n,"_d_",d,".pdf"),width=8,height=7)
    print(p+geom_boxplot()+scale_y_continuous(limits = c(-M,M)))
    dev.off()
  }
}
