#########################
# analyse
# 24/11/2023 vincent.brault@univ-grenoble-alpes.fr
#########################

rm(list=ls())
Niter<-20
N_time<-5
K=3
source("R/Class_MixSegTime.R")

nd_list <- 50*(1:10)
Nom_com<-paste0("Plan_t_fixe_time/Save/")
########
# Comp classification
########

library(reshape2)

res<-data.frame(n=NULL,d=NULL,alpha=NULL,Time=NULL)

for (n in nd_list){
  for (d in nd_list){
    if (all(sapply(c(0.1,0.2,1),function(alpha){
      file.exists(paste0(Nom_com,"Essai_alpha_",10*alpha,"_d_",d,"_n_",n,".Rdata"))
    }))){
      for (alpha in c(0.1,0.2,1)){
        Nom_fichier<-paste0(Nom_com,"Essai_alpha_",10*alpha,"_d_",d,"_n_",n,".Rdata")
        load(Nom_fichier)
        temp <- melt(time_estim)[,3]/(60*10^(3))
        res<-rbind(res,data.frame(n=rep(n,length(temp)),
                                  d=rep(d,length(temp)),
                                  alpha=rep(alpha,length(temp)),
                                  Time=temp))
      }
    }
  }
}

########
# Plot
########
library(ggplot2)
library(grid)

res$alpha<-as.factor(res$alpha)

########
# By d
########

pdf("Plan_t_fixe_time/Im/time_d.pdf",width=12,height=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 1)))
for (d_it in 1:5){
  d<-50*d_it
  dat<-res[res$d==d,]
  dat$n<-as.factor(dat$n)
  print(ggplot(dat,aes(x=n,y=Time,fill=alpha))+geom_boxplot()+
          ylab(paste0("Time (minutes) for d=",d)),
        vp=viewport(layout.pos.row = d_it, layout.pos.col = 1))
}
dev.off()

########
# By n
########

pdf("Plan_t_fixe_time/Im/time_n.pdf",width=12,height=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 1)))
for (n_it in 1:5){
  n<-50*n_it
  dat<-res[res$n==n,]
  dat$d<-as.factor(dat$d)
  print(ggplot(dat,aes(x=d,y=Time,fill=alpha))+geom_boxplot()+
          ylab(paste0("Time (minutes) for n=",n)),
        vp=viewport(layout.pos.row = n_it, layout.pos.col = 1))
}
dev.off()

#########
# Mixte
#########

pdf("Plan_t_fixe_time/Im/time_nd.pdf",width=12,height=8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
for (it_alpha in 1:2){
  alpha <- c(0.1,1)[it_alpha]
  ##### n = 100
  n<-100
  dat<-res[(res$n==n)&(res$alpha==alpha),]
  dat$d<-as.factor(dat$d)
  q<-ggplot(dat,aes(x=d,y=Time,fill=alpha))+geom_boxplot()+
    ylab(paste0("Time (minutes) for n=",n))+
    theme(legend.position = "none")
  if (it_alpha==1){
    q<-q+ggtitle(label = expression(paste(alpha, "=0.1",sep="")))
  }else{
    q<-q+ggtitle(label = expression(paste(alpha, "=1",sep="")))
  }
  print(q,
        vp=viewport(layout.pos.row = 1, layout.pos.col = it_alpha))
  ##### d = 100
  d<-100
  dat<-res[(res$d==d)&(res$alpha==alpha),]
  dat$n<-as.factor(dat$n)
  q <- ggplot(dat,aes(x=n,y=Time,fill=alpha))+geom_boxplot()+
    ylab(paste0("Time (minutes) for d=",d))+
    theme(legend.position = "none")
  if (it_alpha==1){
    q<-q+ggtitle(label = expression(paste(alpha, "=0.1",sep="")))
  }else{
    q<-q+ggtitle(label = expression(paste(alpha, "=1",sep="")))
  }
  print(q,
        vp=viewport(layout.pos.row = 2, layout.pos.col = it_alpha))
}
dev.off()


