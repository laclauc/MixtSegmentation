##########################
# Algorithme EM
# MAJ 18/02/2021
# vincent.brault@univ-grenoble-alpes.fr
###########################

Init_Em<-function(bY,K){
  if (class(bY)=="array"){
    Y<-sapply(1:dim(bY)[1],function(i){c(bY[i,,])})
  }else{
    Y<-bY
  }
  n=dim(Y)[1]
  d=dim(Y)[2]
  ind<-sample(x = 1:n,size = 1)
  if (K>1){
    D<-rep(0,n)
    for (k in 2:K){
      D<-D+apply((Y-matrix(1,n,1)%*%Y[ind[k-1],])^2,MARGIN = 1,sum)
      ind<-c(ind,
             sample(x = (1:n)[-ind],size = 1,prob = D[-ind]/sum(D[-ind])))
    }
  }
  if (class(bY)=="array"){
    bY[ind,,]
  }else{
    bY[ind,]
  }
}


