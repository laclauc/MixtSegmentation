#############################
# Structure Plan
# 19/11/2023
#############################


source("R/Class_MixSegTime.R")
source("R/ProgDyn.R")
source("R/Init.R")
library(wavelets)
library(microbenchmark)

#### Parameter
T_max = 32 ## A day with 32 hours
sigma=1
K=3
time = (1:T_max)/T_max


### Init
Nom_com<-paste0("Plan_t_fixe_time/Save/")

if (!dir.exists(Nom_com)){
  dir.create(Nom_com)
}

Niter<-20
N_time<-5
pb<-txtProgressBar(0,Niter,style=3)
N_max <- 10

for (iter in 1:N_max){
  for (d in 50*(1:iter)){
    t = 1:d    ##

    LK=matrix(NA,K,d)
    for (k in 1:K){
      LK[k,]=sort(rep(1:(k+1),length=d)) ### Segmentation of days
    }
    for (n in 50*(1:iter)){
      cat("\n d=",d," et n=",n,"\n",sep="")
      Nom_fichier<-paste0(Nom_com,"Essai_alpha_",10*alpha,"_d_",d,"_n_",n,".Rdata")
      if ((!file.exists(Nom_fichier))&((n==100)|(d==100))){
        #### Time
        time_estim <- matrix(NA,Niter,N_time)
        #### Loop
        for (it in 1:Niter){
          set.seed(seed = 16*n+72*d+80*alpha+24749*it,kind = "Mersenne-Twister")
          Z = sample(1:K,n, replace=TRUE) ### appartenance aux clusters
          X = array(NA, dim = c(n,d,T_max))
          for (i in 1:n){
            for (j in 1:d){
              X[i,j,] = (-1)^(Z[i]) * alpha *cos(2*pi*time*LK[Z[i],j]) + rnorm(T_max,0,sigma)
            }
          }
          Xproj = array(NA, dim = c(n,d,4)) ### n: le nombre d'observation, d: le nombre de jours, 4: le nombre de coeff
          for (i in 1:n){
            for (j in 1:d){
              Xproj[i,j,] = dwt(X[i,j,], filter="haar", n.levels = 3)@V[[3]]
              #    Xrecon[i,j,] = rep(Xproj[i,j,], each = d/length(Xproj[i,j,]))
            }
          }
          #### Estimation
          temp <- microbenchmark(EM_Dyn(Xproj,K=K,Init = 1,Essai_max = 50,Essai_Min = 1,
                                        L=1:K,verbatim = FALSE,
                                        niter=100,epsilon=10^(-5)),
                                 times=N_time)

          time_estim[it,]<- temp$time/(10^6)
          #### Progression
          setTxtProgressBar(pb,it)
        }
        #### Save
        save("time_estim",file=Nom_fichier)
      }
    }
  }
}
######### Pour le reste des calculs
for (iter in 1:N_max){
  for (d in 50*(1:iter)){
    t = 1:d    ##

    LK=matrix(NA,K,d)
    for (k in 1:K){
      LK[k,]=sort(rep(1:(k+1),length=d)) ### Segmentation of days
    }
    for (n in 50*(1:iter)){
      cat("\n d=",d," et n=",n,"\n",sep="")
      Nom_fichier<-paste0(Nom_com,"Essai_alpha_",10*alpha,"_d_",d,"_n_",n,".Rdata")
      if ((!file.exists(Nom_fichier))){
        #### Time
        time_estim <- matrix(NA,Niter,N_time)
        #### Loop
        for (it in 1:Niter){
          set.seed(seed = 16*n+72*d+80*alpha+24749*it,kind = "Mersenne-Twister")
          Z = sample(1:K,n, replace=TRUE) ### appartenance aux clusters
          X = array(NA, dim = c(n,d,T_max))
          for (i in 1:n){
            for (j in 1:d){
              X[i,j,] = (-1)^(Z[i]) * alpha *cos(2*pi*time*LK[Z[i],j]) + rnorm(T_max,0,sigma)
            }
          }
          Xproj = array(NA, dim = c(n,d,4)) ### n: le nombre d'observation, d: le nombre de jours, 4: le nombre de coeff
          for (i in 1:n){
            for (j in 1:d){
              Xproj[i,j,] = dwt(X[i,j,], filter="haar", n.levels = 3)@V[[3]]
              #    Xrecon[i,j,] = rep(Xproj[i,j,], each = d/length(Xproj[i,j,]))
            }
          }
          #### Estimation
          temp <- microbenchmark(EM_Dyn(Xproj,K=K,Init = 1,Essai_max = 50,Essai_Min = 1,
                                        L=1:K,verbatim = FALSE,
                                        niter=100,epsilon=10^(-5)),
                                 times=N_time)

          time_estim[it,]<- temp$time/(10^6)
          #### Progression
          setTxtProgressBar(pb,it)
        }
        #### Save
        save("time_estim",file=Nom_fichier)
      }
    }
  }
}
