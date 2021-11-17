################
# Pour le beamer
################

rm(list=ls())

n=30
if (!file.exists("ref.Rdata")){
  ref<-1.4
  mu<-c(rep(-ref,n/3),rep(ref,n/3),rep(-ref,n/3))
  x<-mu+rnorm(n)
  plot(x,type="o",pch=19,col=c(rep(1,n/3),rep(2,n/3),rep(3,n/3)))
  min(x[1:10+n/3])-max(x[c(1:(n/3),(2*n/3+1):n)])

  alea<-sample(1:n,n)
  plot(x[alea],type="o",pch=19,col=c(rep(1,n/3),rep(2,n/3),rep(3,n/3))[alea])

  save("x","alea",file="ref.Rdata")
}else{
  load("ref.Rdata")
}


plot(x,type="o",pch=19,col=c(rep(1,n/3),rep(2,n/3),rep(3,n/3)))
min(x[1:10+n/3])-max(x[c(1:(n/3),(2*n/3+1):n)])
summary(x)

dat<-matrix(c(x,rep(1,n/3),rep(2,n/3),rep(1,n/3)),ncol=2)
dat[alea,]


##### Matrix
library(Matrix)
if (!file.exists("ref_matrix.Rdata")){
  z<-c(rep(1,4),rep(2,4),rep(1,4),rep(3,4))
  w<-c(rep(1,4),rep(2,3),rep(1,4),rep(3,3))
  Mat<-matrix(0,16,14)
  ### Carre haut gauche
  Mat[order(z)[1:8],order(w)[1:8]]<-sample(0:1,8^2,replace = TRUE,prob = c(0.1,0.9))
  ### complement carre
  Mat[order(z)[9:12],order(w)[1:8]]<-sample(0:1,8*4,replace = TRUE,prob = c(0.9,0.1))
  Mat[order(z)[1:8],order(w)[9:11]]<-sample(0:1,8*3,replace = TRUE,prob = c(0.9,0.2))
  Mat[order(z)[9:12],order(w)[9:11]]<-sample(0:1,4*3,replace = TRUE,prob = c(0.5,0.5))
  ### La fin
  Mat[order(z)[1:8],order(w)[12:14]]<-sample(0:1,8*3,replace = TRUE,prob = c(0.6,0.4))
  Mat[order(z)[9:12],order(w)[12:14]]<-sample(0:1,4*3,replace = TRUE,prob = c(0.4,0.6))
  Mat[order(z)[13:16],order(w)[12:14]]<-sample(0:1,4*3,replace = TRUE,prob = c(0.1,0.9))
  Mat[order(z)[13:16],order(w)[9:11]]<-sample(0:1,4*3,replace = TRUE,prob = c(0.5,0.5))
  Mat[order(z)[13:16],order(w)[1:8]]<-sample(0:1,4*8,replace = TRUE,prob = c(0.7,0.3))
  save("Mat","z","w",file="ref_matrix.Rdata")
}else{
  load("ref_matrix.Rdata")
}



