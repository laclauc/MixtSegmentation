rm(list=ls())
library(wavelets)

n = 100
d = 32
T = 48
t = 1:T
sigma = 0.5
alpha = 3

LK = list(rep(1:2,each = T/2),rep(1:3, each = T/3), rep(1:4, each=T/4) ) ### appartenance aux segments
Z = sample(1:3,n, replace=TRUE) ### appartenance aux clusters
time = (1:d)/d

X = array(NA, dim = c(n,T,d))
for (i in 1:n){
  for (j in 1:T){
    X[i,j,] = (-1)^(Z[i]) * alpha *cos(2*pi*time*LK[[Z[i]]][j]) + rnorm(d,0,sigma)
  }
}
### concaténation des dimensions pour voir la rupture temporelle
#Xsemaine = matrix(NA, nrow = n, ncol = T*d)
#for (i in 1:n){
#  for (ind in 1:T){
#    Xsemaine[i,((ind-1)*d+1):(ind*d)] = X[i,ind,]
#  }
#}

#Xrecon = X
Xproj = array(NA, dim = c(n,T,4)) ### n: le nombre d'observation, T: le nombre de jours, 4: le nombre de coeff
for (i in 1:n){
  for (j in 1:T){
    Xproj[i,j,] = dwt(X[i,j,], filter="haar", n.levels = 3)@V[[3]]
#    Xrecon[i,j,] = rep(Xproj[i,j,], each = d/length(Xproj[i,j,]))
  }
}
par(mfrow=c(1,3))
matplot(Xproj[3,,], pch=1, lty=1)
matplot(Xproj[1,,], pch=1, lty=1)
matplot(Xproj[4,,], pch=1, lty=1)
