#############################
# Structure Plan
# 23/11/2023
#############################


source("R/Class_MixSegTime.R")
source("R/ProgDyn.R")
source("R/Init.R")
library(wavelets)

#### Parameter
T_max = 32 ## A day with 32 hours
t = 1:d    ##
sigma=1
K=3

LK=matrix(NA,K,d)
for (k in 1:K){
  LK[k,]=sort(rep(1:(k+1),length=d)) ### Segmentation of days
}
time = (1:T_max)/T_max



### Init
Nom_com<-paste0("Plan_Comp/Plan_n_",n,"_d_",d,"_alpha_",10*alpha,"/")

if (!dir.exists(Nom_com)){
  dir.create(Nom_com)
}

Niter<-100
pb<-txtProgressBar(0,Niter,style=3)

for (it in 1:Niter){
  Nom_fichier<-paste0(Nom_com,"Essai_",it,".Rdata")
  if (!file.exists(Nom_fichier)){
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
    ############
    # Means
    ############
    Moy<-Calcul_Moy(x=Xproj,verbatim=FALSE)
    ############
    # True model
    ############
    res_Dyn <- EM_Dyn(Xproj,K=K,Moy = Moy,
              L=1:K,verbatim=FALSE,
              niter=100,epsilon=10^(-5))
    ############
    # K=1 (determinist)
    ############
    res_K_1_L_2 <- EM_Dyn(Xproj,K=1,Moy = Moy,Essai_Min = 1,
                        L=2,verbatim=FALSE,
                        niter=100,epsilon=10^(-5))
    res_K_1_L_3 <- EM_Dyn(Xproj,K=1,Moy = Moy,Essai_Min = 1,
                      L=3,verbatim=FALSE,
                      niter=100,epsilon=10^(-5))
    res_K_1_L_5 <- EM_Dyn(Xproj,K=1,Moy = Moy,Essai_Min = 1,
                        L=5,verbatim=FALSE,
                        niter=100,epsilon=10^(-5))
    ############
    # K=3 and L=0
    ############
    res_K_3_L_0 <- EM_Dyn(Xproj,K=K,Moy = Moy,
                        L=0,verbatim=FALSE,
                        niter=100,epsilon=10^(-5))
    save("Xproj","Z","res_Dyn",
         "res_K_1_L_2","res_K_1_L_3",
         "res_K_1_L_5","res_K_3_L_0",
         file=Nom_fichier)
  }
  setTxtProgressBar(pb,it)
}

