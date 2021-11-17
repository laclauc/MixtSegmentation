##################
# Pour les essais
##################

source("R/Class_MixSegTime.R")
source("R/ProgDyn.R")
source("R/LASSO.R")
source("R/Init.R")

curve<-simu_courbe(1000,50,p = 2,K = 3,L = 1:3)

##### LASSO
Y<-curve$x
K=3
niter=100
ite=1
lambda=0.1
L=1:3
epsilon=10^(-5)

#################
res.init = Init_Em(Y,K) ### TO BE IMPROVED
par(mfrow=c(1,1))
matplot(t(res.init[,,2]),type="l")








plot(new(Class="MixSegTime",mu=curve$mu,
         sigma=curve$sigma,TL=curve$TL,
         pik=curve$pik,method="EM",z=curve$z),curve$x)
res_EM=EM_Dyn(curve$x,K=curve$K,
       L=sapply(1:curve$K,function(k){length(curve$TL[[k]])-2}),
       niter=100,epsilon=10^(-5))
plot(res_EM,curve$x)

res_LASSO = EMalgorithmLasso_div(curve$x, K = 3)
plot(res_LASSO,curve$x)

