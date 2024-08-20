#########################
# analyse
# 17/02//2023 vincent.brault@univ-grenoble-alpes.fr
#########################

rm(list=ls())
Niter<-100
K=3
source("R/Class_MixSegTime.R")

######################
# Comparison with mixture model gaussian
######################

library(bikm1)

TL_res<-matrix(0,nrow=0,ncol=4)
TL_res_2<-matrix(0,nrow=0,ncol=4)
TL_res_3<-matrix(0,nrow=0,ncol=4)
TL_res_5<-matrix(0,nrow=0,ncol=4)
ARI_res<-matrix(0,nrow=0,ncol=4)
ARI_res_simple<-matrix(0,nrow=0,ncol=4)
CE_res<-matrix(0,nrow=0,ncol=4)
CE_res_simple<-matrix(0,nrow=0,ncol=4)

for (n in c(100,1000)){
  for (d in c(50,100)){
    LK=matrix(NA,K,d)
    for (k in 1:K){
      LK[k,]=sort(rep(1:(k+1),length=d)) ### Segmentation of days
    }
    TL_ref<-sapply(1:K,function(k){which(!duplicated(LK[k,]))[-1]-1})
    for (alpha in c(0.1,0.2,1)){
      Nom_com<-paste0("Plan_Comp/Plan_n_",n,"_d_",d,"_alpha_",10*alpha,"/")
      if (!dir.exists(Nom_com)){
        dir.create(Nom_com)
      }
      dist_PG<-rep(NA,Niter)
      dist_PG_2<-rep(NA,Niter)
      dist_PG_3<-rep(NA,Niter)
      dist_PG_5<-rep(NA,Niter)
      ARI<-rep(NA,Niter)
      ARI_simple<-rep(NA,Niter)
      CE<-rep(NA,Niter)
      CE_simple<-rep(NA,Niter)
      for (it in 1:Niter){
        Nom_fichier<-paste0(Nom_com,"Essai_",it,".Rdata")
        if (file.exists(Nom_fichier)){
          load(Nom_fichier)
          if (length(res_Dyn@TL)==K){
            ## Classical
            dist_PG[it]<-max(sapply(1:K,function(k){
              max(abs(TL_ref[[k]]-res_Dyn@TL[[k]][2:(length(TL_ref[[k]])+1)]))}))
            ## Comp
            dist_PG_2[it]<-max(abs(TL_ref[[2]]-res_K_1_L_2@TL[[1]][2:3]))
            dist_PG_3[it]<-max(abs(TL_ref[[3]]-res_K_1_L_3@TL[[1]][2:4]))
            dist_PG_5[it]<-max(abs(sort(c(TL_ref[[2]],TL_ref[[3]]))-res_K_1_L_5@TL[[1]][2:6]))
            ## Mixture
            ARI[it]<-ARI(v = Z,res_Dyn@z)$ari
            ARI_simple[it]<-ARI(v = Z,res_K_3_L_0@z)$ari
            CE[it]<-NCE_simple(Z,res_Dyn@z)
            CE_simple[it]<-NCE_simple(Z,res_K_3_L_0@z)
          }
        }
      }
      N_part<-sum(!is.na(dist_PG))
      TL_res<-rbind(TL_res,
                    matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),dist_PG[!is.na(dist_PG)]/d),
                           nrow=N_part))
      TL_res_2<-rbind(TL_res_2,
                      matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),dist_PG_2[!is.na(dist_PG)]/d),
                             nrow=N_part))
      TL_res_3<-rbind(TL_res_3,
                      matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),dist_PG_3[!is.na(dist_PG)]/d),
                             nrow=N_part))
      TL_res_5<-rbind(TL_res_5,
                      matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),dist_PG_5[!is.na(dist_PG)]/d),
                             nrow=N_part))
      ARI_res<-rbind(ARI_res,
                     matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),ARI[!is.na(dist_PG)]),
                            nrow=N_part))
      ARI_res_simple<-rbind(ARI_res_simple,
                            matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),ARI_simple[!is.na(dist_PG)]),
                                   nrow=N_part))
      CE_res<-rbind(CE_res,
                    matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),CE[!is.na(dist_PG)]),
                           nrow=N_part))
      CE_res_simple<-rbind(CE_res_simple,
                           matrix(c(rep(n,N_part),rep(d,N_part),rep(alpha,N_part),CE_simple[!is.na(dist_PG)]),
                                  nrow=N_part))
    }
  }
}

TL_res_bis<-TL_res

#######
# Display estimation of TL
#######
library(ggplot2)
TL_res<-data.frame(TL_res)
names(TL_res)=c("n","d","alpha","dist")
TL_res$Model<-rep("MixSeg",nrow(TL_res))
for (l in c(2,3,5)){
  temp<-data.frame(get(paste0("TL_res_",l)))
  names(temp)=c("n","d","alpha","dist")
  temp$Model<-rep(paste0(l," breaks"),nrow(temp))
  TL_res<-rbind(TL_res,temp)
}

M<-max(TL_res$dist)
for (n in c(100,1000)){
  for (d in c(50,100)){
    dat<-TL_res[(TL_res$n==n)&(TL_res$d==d),]
    dat$alpha<-as.factor(dat$alpha)
    p<-ggplot(dat,aes(x=alpha,y=dist,fill=Model))
    pdf(paste0("Plan_Comp/Im/Comp_TL_n_",n,"_d_",d,".pdf"),width=9,height=7)
    print(p+geom_boxplot()+
            theme()+scale_y_continuous(limits = c(0,M)))
    dev.off()
  }
}

#######
# Display estimation of TL
#######
TL_res<-data.frame(TL_res_bis)
names(TL_res)=c("n","d","alpha","dist")
TL_res$Model<-rep("MixSeg",nrow(TL_res))
## 5 breaks
temp<-data.frame(TL_res_5)
names(temp)=c("n","d","alpha","dist")
temp$Model<-rep("SimpleSeg",nrow(temp))
TL_res<-rbind(TL_res,temp)

M<-max(TL_res$dist)
for (n in c(100,1000)){
  for (d in c(50,100)){
    dat<-TL_res[(TL_res$n==n)&(TL_res$d==d),]
    dat$alpha<-as.factor(dat$alpha)
    p<-ggplot(dat,aes(x=alpha,y=dist,fill=Model))+geom_boxplot()+
      scale_y_continuous(limits = c(0,M))
    if ((n==1000)&(d==100)){
      p<-p+theme()
      pdf(paste0("Plan_Comp/Im/Comp_5_TL_n_",n,"_d_",d,".pdf"),width=7,height=5)
      print(p)
      dev.off()
    }else{
      p<-p+theme(legend.position = "none")
      pdf(paste0("Plan_Comp/Im/Comp_5_TL_n_",n,"_d_",d,".pdf"),width=6,height=5)
      print(p)
      dev.off()
    }
  }
}

#########
# Table ARI
#########
ARI_tab<-matrix("",8,3)
it<-1
ARI_res<-data.frame(ARI_res)
names(ARI_res)=c("n","d","alpha","value")
ARI_res_simple<-data.frame(ARI_res_simple)
names(ARI_res_simple)=c("n","d","alpha","value")
for (n in 10^(2:3)){
  for (d in c(50,100)){
    for (j in 1:3){
      alpha<-c(0.1,0.2,1)[j]
      dat<-ARI_res[(ARI_res$n==n)&(ARI_res$d==d)&(ARI_res$alpha==alpha),]
      dat_simple<-ARI_res_simple[(ARI_res_simple$n==n)&(ARI_res_simple$d==d)&(ARI_res_simple$alpha==alpha),]
      ARI_tab[it,j]<-paste0(signif(mean(dat$value),2)," (",signif(sd(dat$value),2),")")
      ARI_tab[it+1,j]<-paste0(signif(mean(dat_simple$value),2)," (",signif(sd(dat_simple$value),2),")")
    }
    it<-it+2
  }
}

row.names(ARI_tab)<-c(sapply(c(100,1000),function(n){
  sapply(c(50,100),function(d){
    c(paste0("\\hline\n \\multirow{2}{*}{(",n,",",d,")}&MixSeg"),"\n &Simple")
  })
}))
colnames(ARI_tab)<-c("&0.1",as.character(c(0.2,1)))
write.table(ARI_tab,file = "Plan_Comp/Tab/Comp_ARI.tex",quote = FALSE,sep = "&",eol = "\\\\")

#########
# Table ARI
#########
CE_tab<-matrix("",8,3)
it<-1
CE_res<-data.frame(CE_res)
names(CE_res)=c("n","d","alpha","value")
CE_res_simple<-data.frame(CE_res_simple)
names(CE_res_simple)=c("n","d","alpha","value")
for (n in 10^(2:3)){
  for (d in c(50,100)){
    for (j in 1:3){
      alpha<-c(0.1,0.2,1)[j]
      dat<-CE_res[(CE_res$n==n)&(CE_res$d==d)&(CE_res$alpha==alpha),]
      dat_simple<-CE_res_simple[(CE_res_simple$n==n)&(CE_res_simple$d==d)&(CE_res_simple$alpha==alpha),]
      if (all(dat$value<10^(-15))){
        CE_tab[it,j]<-"0 (0)"
      }else{
        CE_tab[it,j]<-paste0(signif(mean(dat$value),2)," (",signif(sd(dat$value),2),")")
      }
      if (all(dat_simple$value<10^(-15))){
        CE_tab[it+1,j]<-"0 (0)"
      }else{
        CE_tab[it+1,j]<-paste0(signif(mean(dat_simple$value),2)," (",signif(sd(dat_simple$value),2),")")
      }
    }
    it<-it+2
  }
}

row.names(CE_tab)<-c(sapply(c(100,1000),function(n){
  sapply(c(50,100),function(d){
    c(paste0("\\hline\n \\multirow{2}{*}{(",n,",",d,")}&MixSeg"),"\n &Simple")
  })
}))
colnames(CE_tab)<-c("&0.1",as.character(c(0.2,1)))
write.table(CE_tab,file = "Plan_Comp/Tab/Comp_CE.tex",quote = FALSE,sep = "&",eol = "\\\\")

