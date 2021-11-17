##########################
# Package
# MAJ 12/02/2020
# vincent.brault@univ-grenoble-alpes.fr
###########################


setClass(
  Class="MixSegTime",
  representation=representation(
    mu="list",  ### The values of each means
    sigma="matrix", ### The variance
    TL="list",  ### Each group of Lk break points
    pik="numeric",
    method="character",
    z="numeric",
    Q="numeric")
)


setMethod(
  f="plot",
  signature="MixSegTime",
  definition=function(x,y,png=FALSE){

    if (!is.array(y)){
      stop("y must be the original data")
    }

    n=dim(y)[1]
    d=dim(y)[2]
    p=dim(y)[3]

    K=length(x@pik)

    #x11(width=60, height=30)
    if (png){
      png("Im/Affichage.png",width=200*(K+1),height=200*p)
    }
    par_temp=graphics::par()
    graphics::par(mfrow=c(p,K+1),oma=c(0,0,3,0))
    for (r in 1:p){
      M=1.1*max(y[,,r])
      m=0.9*min(y[,,r])
      for (k in 0:K){
        if (k==0){
          if (r==1){
            matplot(t(y[,,r]),col=x@z,main="Toutes les courbes",ylab=paste0("Dimension ",r),type="l",ylim=c(m,M))
          }else{
            matplot(t(y[,,r]),col=x@z,ylab=paste0("Dimension ",r),type="l",ylim=c(m,M))
          }
        }else{
          nk=sum(x@z==k)
          if (nk==0){
            if (r==1){
              plot(0,0,xaxt="n",pch="",main=paste0("Cluster ",k))
            }else{
              plot(0,0,xaxt="n",pch="")
            }
            text("Cluster vide",x = 0,y=0)
          }else{
            if (r==1){
              if (nk==1){
                plot(y[x@z==k,,r],main=paste0("Cluster ",k," avec ",nk," courbes"),type="l",col=k,ylab="",ylim=c(m,M))
              }else{
                matplot(t(y[x@z==k,,r]),main=paste0("Cluster ",k," avec ",nk," courbes"),type="l",col=k,ylab="",ylim=c(m,M))
              }
            }else{
              if (nk==1){
                plot(y[x@z==k,,r],type="l",col=k,ylab="",ylim=c(m,M))
              }else{
                matplot(t(y[x@z==k,,r]),type="l",col=k,ylab="",ylim=c(m,M))
              }
            }
            # Affichage moyennes
            if (k==1){
              col="red"
            }else{
              col="black"
            }
            lines(unlist(lapply(1:(length(x@TL[[k]])-1), function(l){
              rep(x@mu[[k]][l,r],x@TL[[k]][l+1]-x@TL[[k]][l])
            })),col=col,lty=2,lwd=2)
            # Affichage ruptures
            if (k==4){
              col="black"
            }else{
              col="blue"
            }
            abline(v=x@TL[[k]][-c(1,length(x@TL[[k]]))],lty=3,lwd=2,col=col)
          }
        }
      }
    }
    if (png){
      dev.off()
    }
  })

simu_courbe=function(n,d,p,K=NULL,L=NULL,sigma=NULL,mu_max=5){
  if (is.null(K)){
    K=rpois(1,2)+2
  }
  if (is.null(L)){
    TL=lapply(1:K,function(k){
      c(0,sort(sample(1:(d-1),size = rpois(1,3)+1,replace = FALSE)),d)
    })
  }else{
    TL=lapply(1:K,function(k){
      Temp<-c(0,sort(sample(1:(d-1),size = L[k],replace = FALSE)),d)
      while (any((Temp[2:length(Temp)]-Temp[1:(length(Temp)-1)])<5)){
        Temp<-c(0,sort(sample(1:(d-1),size = L[k],replace = FALSE)),d)
      }
      Temp
    })
  }
  if (is.null(sigma)){
    sigma=runif(n = K,min = 0.1,max = 2)
  }
  mu<-lapply(1:K,function(k){
    temp<-matrix(runif(p*(length(TL[[k]])-1)),ncol=p)
    temp<-temp-matrix(1,length(TL[[k]])-1,1)%*%apply(temp,2,min)
    temp<-temp/matrix(1,length(TL[[k]])-1,1)%*%apply(temp,2,max)*mu_max
  })
  x=array(0,c(n,d,p))
  pik=runif(K)
  pik=pik/sum(pik)
  z=sample(1:K,size = n,prob = pik,replace=TRUE)
  for (k in 1:K){
    nk=sum(z==k)
    for (r in 1:p){
      x[z==k,,r]=matrix(rep(unlist(lapply(1:(length(TL[[k]])-1), function(l){
        rep(mu[[k]][l,r],(TL[[k]][l+1]-TL[[k]][l]))
      })),nk),nrow=nk,byrow=TRUE)+rnorm(d*nk)*sigma[k]
    }
  }
  list(x=x,z=z,pik=pik,TL=TL,mu=mu,K=K,sigma=as.matrix(sigma))
}
simu_courbe_equ=function(n,d,p,K=NULL,L=NULL,sigma=NULL,mu_max=5){
  if (is.null(K)){
    K=rpois(1,2)+2
  }
  if (is.null(L)){
    TL=lapply(1:K,function(k){
      c(0,sort(sample(1:(d-1),size = rpois(1,3)+1,replace = FALSE)),d)
    })
  }else{
    TL=lapply(1:K,function(k){
      Temp<-c(0,sort(sample(1:(d-1),size = L[k],replace = FALSE)),d)
      while (any((Temp[2:length(Temp)]-Temp[1:(length(Temp)-1)])<5)){
        Temp<-c(0,sort(sample(1:(d-1),size = L[k],replace = FALSE)),d)
      }
      Temp
    })
  }
  if (is.null(sigma)){
    sigma=runif(n = K,min = 0.1,max = 2)
  }
  mu<-lapply(1:K,function(k){
    matrix(sapply(1:p,function(r){
      temp<-rep(c(0,mu_max)[sample(1:2,2)],(length(TL[[k]])+1)/2)[1:(length(TL[[k]]))]
    }),ncol=p)
  })
  x=array(0,c(n,d,p))
  pik=rbeta(K,4,4)
  pik=pik/sum(pik)
  z=sample(1:K,size = n,prob = pik,replace=TRUE)
  for (k in 1:K){
    nk=sum(z==k)
    for (r in 1:p){
      x[z==k,,r]=matrix(rep(unlist(lapply(1:(length(TL[[k]])-1), function(l){
        rep(mu[[k]][l,r],(TL[[k]][l+1]-TL[[k]][l]))
      })),nk),nrow=nk,byrow=TRUE)+rnorm(d*nk)*sigma[k]
    }
  }
  list(x=x,z=z,pik=pik,TL=TL,mu=mu,K=K,sigma=as.matrix(sigma))
}
