#if (!require(grplasso)) install.packages('grplasso')
library(gglasso)
#if (!require(glmnet)) install.packages('glmnet')
library(glmnet)
#EMalgorithmLasso(x[,,1], 2, 10, rep(0.0001,2))

###### p=1 #####

EMalgorithmLasso_1 = function(Y, K, lambda=NULL, niter=100,epsilon=10^(-5)){
  ### lambda= vector of length K
  n = dim(Y)[1]
  d = dim(Y)[2]

  ### Initialization
  res.init = Init_Em(Y,K) ### TO BE IMPROVED
  mu = res.init
  z = apply(sapply(1:K,function(k){
    sapply(1:n,function(i){
      sum((Y[i,]-mu[k,])^2)
    })
  }),1,which.min)
  s<-matrix(sapply(1:n,function(i){
    temp<-rep(0,K)
    temp[z[i]]<-1
    temp
  }),ncol=K,byrow=TRUE)
  Sigma<-sapply(1:K,function(k){
    S_temp<-sum(s[,k])
    unlist(sapply(1:d,function(l){
      sum(s[,k] %*% ((Y[,l] - mu[k,l])^2))
    }))/S_temp
  })
  ## Affichage
  # matplot(t(Y),type="l",col=apply(s,1,which.max))
  # for (k in 1:K){
  #   lines(mu[k,],col=3+k,lwd=2,lty=2)
  # }
  loc = rbeta(K,4,4)
  pi_k =  loc / sum(loc)
  ## Necessary
  Ld = matrix(0,nrow = d, ncol = d)
  Ld[lower.tri(Ld, diag=TRUE)] = 1
  num = rep(0,K)
  W_new<-Inf
  pb<-txtProgressBar(0,niter,style=3)
  Ycal=c(Y)
  lambda<-numeric(niter)
  for (ite in 1:niter){
    setTxtProgressBar(pb = pb,ite)
    #### old
    mu_old<-mu
    pi_k_old<-pi_k
    Sigma_old<-Sigma
    W_old<-W_new
    #### Estep ####
    s<-sapply(1:K,function(k){
      sapply(1:n,function(i){
        -1/2 * t(Y[i,] - mu[k,]) %*% solve(diag(Sigma[,k])) %*% (Y[i,] - mu[k,])
      })+log(pi_k[k]) -1/2 * sum(log(Sigma[,k]))
    })
    s<-exp(s-apply(s,1,max)%*%matrix(1,1,K))
    s<-s/(apply(s,1,sum)%*%matrix(1,1,K))

    #Z = apply(s,1,which.max)

    #### Mstep ####
    ## pik
    pi_k = apply(s,2,mean)
    ## moyenne
    Xcal = kronecker(Ld,s)
    #res<-glmnet(Xcal, Ycal, intercept = FALSE)
    lambda[ite]<-cv.glmnet(Xcal, Ycal, intercept = FALSE)$lambda.1se
    res<-glmnet(Xcal, Ycal, intercept = FALSE,
                lambda = mean(lambda[1:ite]))
    B<-matrix(res$beta,nrow=K)
    mu<-B%*%t(Ld)
    TL<-sapply(1:K,function(k){unique(c(0,which(B[k,]!=0)-1,d))})
    ## Affichage
    matplot(t(Y),type="l",col=apply(s,1,which.max),main=paste0("lambda=",mean(lambda[1:ite])))
    for (k in 1:K){
      lines(mu[k,],col=3+k,lwd=2,lty=2)
      abline(v=TL[[k]],col=3+k,lty=2,lwd=2)
    }
    ## variance
    Sigma<-sqrt(sapply(1:K,function(k){
      L<-length(TL[[k]])-1
      S_temp<-sum(s[,k])
      if (S_temp!=0){
        unlist(sapply(1:L,function(l){
          taille<-(TL[[k]][l]+1):TL[[k]][l+1]
          rep(sum(s[,k] %*% ((Y[,taille]
                              - matrix(rep(mu[k,taille],n),n,byrow=TRUE))^2))
              ,length(taille))
        }))/S_temp
      }else{
        rchisq(d, df = max(Yvec)-min(Yvec))
      }
    }))

    ## Break ?
    W_new<-mean(abs(mu_old-mu))+mean(abs(pi_k_old-pi_k))+
      mean(abs(Sigma_old-Sigma))
    if ((abs(W_new-W_old)/W_new)<epsilon){
      break
    }
  }
  res = list(mu = mu, Sigma = Sigma, pi = pi_k)
  return(res)
}

####### p > 1 ######
EMalgorithmLasso = function(Y, K, niter=100, lambda=1, package =  'gglasso',epsilon=10^(-5)){
  ### lambda= vector of length K
  n = dim(Y)[1]
  dold = dim(Y)[2]
  pold = dim(Y)[3]
  Yvec = matrix(c(Y), nrow = n) #### Mets les lignes bouts a bout qu'importe les dimensions
  d = dim(Yvec)[2]
  ### Initialization
  res.init = Init_Em(Yvec,K) ### TO BE IMPROVED
  mu = res.init
  z = apply(sapply(1:K,function(k){
    sapply(1:n,function(i){
      sum((Yvec[i,]-mu[k,])^2)
    })
  }),1,which.min)
  s<-matrix(sapply(1:n,function(i){
    temp<-rep(0,K)
    temp[z[i]]<-1
    temp
  }),ncol=K,byrow=TRUE)
  Sigma<-sapply(1:K,function(k){
    S_temp<-sum(s[,k])
    unlist(sapply(1:d,function(l){
      sum(s[,k] %*% ((Yvec[,l] - mu[k,l])^2))
    }))/S_temp
  })
  s_col<-colSums(s)
  if (any(s_col<=1)){
    for (s_it in 1:length(s_col)){
      if (s_col[s_it]<=1){
        Sigma[,s_it]<-rchisq(d, df = max(Yvec)-min(Yvec))
      }
    }
  }
  #res.kmeans = kmeans(Xvec, K)
  #mu = t(res.kmeans$centers)
  #Sigma = matrix(rchisq(d*K, df = K*d-1), ncol = K) ### Pourquoi simuler si on peut faire le calcul exact ?
  #loc = runif(K)
  #pi_k =  loc / sum(loc)
  #pi_k = c(0.5,0.5)
  pi_k = apply(s,2,mean)
  #s = matrix(nrow = n, ncol = K)
  ## Necessary
  Ld = matrix(0,nrow = dold, ncol = dold)
  Ld[lower.tri(Ld, diag=TRUE)] = 1
  Ip<-matrix(0,pold,pold)
  diag(Ip)<-1
  Ld<-kronecker(Ip,Ld)
  W_new<-0
  pb<-txtProgressBar(0,niter,style=3)
  Ycal=c(Yvec)
  lambda<-numeric(niter)
  for (ite in 1:niter){
    setTxtProgressBar(pb = pb,ite)
    #### old
    mu_old<-mu
    pi_k_old<-pi_k
    Sigma_old<-Sigma
    W_old<-W_new
    #### Estep ####
    s<-sapply(1:K,function(k){
      sapply(1:n,function(i){
        -1/2 * t(Yvec[i,] - mu[k,]) %*% solve(diag(Sigma[,k])) %*% (Yvec[i,] - mu[k,])
      })+log(pi_k[k]) -1/2 * sum(log(Sigma[,k]))
    })
    s<-exp(s-apply(s,1,max)%*%matrix(1,1,K))
    s<-s/(apply(s,1,sum)%*%matrix(1,1,K))

    #Z = apply(s,1,which.max)

    #### Mstep ####
    ## pik
    pi_k = apply(s,2,mean)
    ## moyenne
    Xcal = kronecker(Ld,s)
    #res<-glmnet(Xcal, Ycal, intercept = FALSE)
    # temp<-cv.gglasso(Xcal, Ycal, intercept = FALSE,group = rep(1:(K*dold),pold)) ### Il va falloir trouver mieux...
    lambda<-cv.gglasso(Xcal, Ycal, intercept = FALSE,group = rep(1:(K*dold),pold))$lambda.1se
    reg.gplasso = gglasso(Xcal, Ycal, lambda = lambda, intercept = FALSE,
                          group = rep(1:(K*dold),pold))


    #res<-glmnet(Xcal, Ycal, intercept = FALSE,lambda = mean(lambda[1:ite]))
    B<-matrix(reg.gplasso$beta,nrow=K)
    mu<-B%*%t(Ld)
    TL<-lapply(1:K,function(k){unique(c(0,which(B[k,]!=0)-1,d))})
    ### Sigma
    ## variance
    Sigma<-sqrt(sapply(1:K,function(k){
      L<-length(TL[[k]])-1
      S_temp<-sum(s[,k])
      if (S_temp!=0){
        unlist(sapply(1:L,function(l){
          taille<-(TL[[k]][l]+1):TL[[k]][l+1]
          rep(sum(s[,k] %*% ((Yvec[,taille]
                              - matrix(rep(mu[k,taille],n),n,byrow=TRUE))^2))
              ,length(taille))
        }))/S_temp
      }else{
        rchisq(d, df = max(Yvec)-min(Yvec))
      }
    }))

    W_new<-mean(abs(mu_old-mu))+mean(abs(pi_k_old-pi_k))+
      mean(abs(Sigma_old-Sigma))
    if (W_new<epsilon){
      break
    }
  }

  TL<-lapply(1:K,function(k){unique(sort(c(0,(TL[[k]]-1)%%dold+1,dold)))})
  ord<-order(sapply(1:K,function(k){length(TL[[k]])}))
  ### On ordonne par ordre croissant
  TL<-lapply(1:K,function(k){TL[[ord[k]]]})
  mu<-mu[ord,]
  Sigma<-Sigma[,ord]
  pi_k<-pi_k[ord]
  ### Vrais
  Q<-sum(sapply(1:K,function(k){
    sum(sapply(1:n,function(i){
      -1/2 * t(Yvec[i,] - mu[k,]) %*% solve(diag(Sigma[,k])) %*% (Yvec[i,] - mu[k,])
    })+log(pi_k[k]) -1/2 * sum(log(Sigma[,k])))
  }))
  ### On met en forme
  mu<-lapply(1:K,function(k){
    sapply(1:pold,function(r){mu[k,TL[[k]][-length(TL[[k]])]+1+(r-1)*dold]})
  })

  ### Renvoi
  new(Class="MixSegTime",mu=mu,
      sigma=Sigma,TL=TL,
      pik=pi_k,method="LASSO",z=apply(s[,ord],1,which.max),Q=Q)
}

####### L defini ######
EMalgorithmLasso_L = function(Y, K, L,lambda=1, niter=100, package =  'gglasso'){
  ### Verif
  if (is.numeric(L)){
    if (length(L)==1){
      L<-rep(L,K)
    }else if((any(L<=0))|(any(floor(L)!=L))|(length(L)!=K)){
      stop("L must be an positive integer or a positive vector integer of size K")
    }
  }else{
    stop("L must be an positive integer or a positive vector integer of size K")
  }
  ### lambda= vector of length K
  n = dim(Y)[1]
  dold = dim(Y)[2]
  pold = dim(Y)[3]
  Yvec <-matrix(sapply(1:n,function(i){c(Y[i,,])}),nrow=n,byrow=TRUE)

  matplot(t(Yvec),type="l",col=curve$z)
  d = dim(Yvec)[2]
  ### Initialization
  res.init = Init_Em(Yvec,K) ### TO BE IMPROVED
  mu = res.init
  z = apply(sapply(1:K,function(k){
    sapply(1:n,function(i){
      sum((Yvec[i,]-mu[k,])^2)
    })
  }),1,which.min)
  s<-matrix(sapply(1:n,function(i){
    temp<-rep(0,K)
    temp[z[i]]<-1
    temp
  }),ncol=K,byrow=TRUE)
  Sigma<-sapply(1:K,function(k){
    S_temp<-sum(s[,k])
    unlist(sapply(1:d,function(l){
      sum(s[,k] %*% ((Yvec[,l] - mu[k,l])^2))
    }))/S_temp
  })
  #res.kmeans = kmeans(Xvec, K)
  #mu = t(res.kmeans$centers)
  #Sigma = matrix(rchisq(d*K, df = K*d-1), ncol = K) ### Pourquoi simuler si on peut faire le calcul exact ?
  #loc = runif(K)
  #pi_k =  loc / sum(loc)
  #pi_k = c(0.5,0.5)
  pi_k = apply(s,2,mean)
  #s = matrix(nrow = n, ncol = K)
  ## Necessary
  Ld = matrix(0,nrow = dold, ncol = dold)
  Ld[lower.tri(Ld, diag=TRUE)] = 1
  Ip<-matrix(0,pold,pold)
  diag(Ip)<-1
  Ld<-kronecker(Ip,Ld)
  W_new<-Inf
  pb<-txtProgressBar(0,niter,style=3)
  Ycal=c(Yvec)
  lambda_temp<-lambda
  for (ite in 1:niter){
    setTxtProgressBar(pb = pb,ite)
    #### old
    mu_old<-mu
    pi_k_old<-pi_k
    Sigma_old<-Sigma
    W_old<-W_new
    #### Estep ####
    s<-sapply(1:K,function(k){
      sapply(1:n,function(i){
        -1/2 * t(Yvec[i,] - mu[k,]) %*% solve(diag(Sigma[,k])) %*% (Yvec[i,] - mu[k,])
      })+log(pi_k[k]) -1/2 * sum(log(Sigma[,k]))
    })
    s<-exp(s-apply(s,1,max)%*%matrix(1,1,K))
    s<-s/(apply(s,1,sum)%*%matrix(1,1,K))

    #Z = apply(s,1,which.max)

    #### Mstep ####
    ## pik
    pi_k = apply(s,2,mean)
    ## moyenne
    Xcal = kronecker(Ld,s)
    ## Search lamnda opt
    Born_int<-c(0,+Inf)
    mult<-1
    Bad_L<-TRUE
    while(Bad_L){
      reg.gplasso = gglasso(Xcal, Ycal, lambda = lambda_temp, intercept = FALSE,
                            group = rep(1:(K*dold),pold))

      B<-matrix(reg.gplasso$beta,nrow=K)
      mu<-B%*%Ld
      TL_double<-lapply(1:K,function(k){unique(c(0,which(B[k,]!=0)-1,d))})
      ## Affichage
      matplot(t(Yvec),type="l",col=apply(s,1,which.max),main=paste0("lambda=",lambda_temp))
      for (k in 1:K){
        lines(mu[k,],col=3+k,lwd=2,lty=2)
        abline(v=TL_double[[k]],col=3+k,lty=2,lwd=2)
      }

      TL<-lapply(1:K,function(k){unique(c(0,sort(unique(TL_double[[k]]%%dold)),dold))})
      L_temp<-sort(sapply(1:K,function(k){length(TL[[k]][-c(1,length(TL[[k]]))])}))
      if (any(sort(L_temp)<sort(L))){
        mult<-mult+1
        Born_int[2]<-lambda_temp
        lambda_temp<-mean(Born_int)
      }else if (any(((mult*sort(L))<sort(L_temp)))){
        Born_int[1]<-lambda_temp
        if (Born_int[2]==+Inf){
          lambda_temp<-2*lambda_temp
        }else{
          lambda_temp<-mean(Born_int)
        }
      }else{
        Bad_L<-FALSE
      }
    }
    W_new<-mean(abs(mu_old-mu))+mean(abs(pi_k_old-pi_k))+
      mean(abs(Sigma_old-Sigma))
    if ((abs(W_new-W_old)/W_new)<epsilon){
      break
    }

  }

  Sigma.tab = array(c(Sigma), dim = c(dold, pold, K))
  res = list(mu = mu.tab, Sigma = Sigma.tab, pi = pi_k, Z = Z)
  return(res)
}

####### L defini ######
EMalgorithmLasso_L_div = function(Y, K, L,lambda=1, niter=100, package =  'gglasso'){
  ### Verif
  if (is.numeric(L)){
    if (length(L)==1){
      L<-rep(L,K)
    }else if((any(L<=0))|(any(floor(L)!=L))|(length(L)!=K)){
      stop("L must be an positive integer or a positive vector integer of size K")
    }
  }else{
    stop("L must be an positive integer or a positive vector integer of size K")
  }
  ### lambda= vector of length K
  n = dim(Y)[1]
  dold = dim(Y)[2]
  pold = dim(Y)[3]
  Yvec <-matrix(sapply(1:n,function(i){c(Y[i,,])}),nrow=n,byrow=TRUE)

  matplot(t(Yvec),type="l",col=curve$z)
  d = dim(Yvec)[2]
  ### Initialization
  res.init = Init_Em(Yvec,K) ### TO BE IMPROVED
  mu = res.init
  z = apply(sapply(1:K,function(k){
    sapply(1:n,function(i){
      sum((Yvec[i,]-mu[k,])^2)
    })
  }),1,which.min)
  s<-matrix(sapply(1:n,function(i){
    temp<-rep(0,K)
    temp[z[i]]<-1
    temp
  }),ncol=K,byrow=TRUE)
  Sigma<-sapply(1:K,function(k){
    S_temp<-sum(s[,k])
    unlist(sapply(1:d,function(l){
      sum(s[,k] %*% ((Yvec[,l] - mu[k,l])^2))
    }))/S_temp
  })
  #res.kmeans = kmeans(Xvec, K)
  #mu = t(res.kmeans$centers)
  #Sigma = matrix(rchisq(d*K, df = K*d-1), ncol = K) ### Pourquoi simuler si on peut faire le calcul exact ?
  #loc = runif(K)
  #pi_k =  loc / sum(loc)
  #pi_k = c(0.5,0.5)
  pi_k = apply(s,2,mean)
  #s = matrix(nrow = n, ncol = K)
  ## Necessary
  Ld = matrix(0,nrow = dold, ncol = dold)
  Ld[lower.tri(Ld, diag=TRUE)] = 1
  Ip<-matrix(0,pold,pold)
  diag(Ip)<-1
  Ld<-kronecker(Ip,Ld)
  W_new<-Inf
  pb<-txtProgressBar(0,niter,style=3)
  Ycal=c(Yvec)
  lambda_temp<-lambda
  for (ite in 1:niter){
    setTxtProgressBar(pb = pb,ite)
    #### old
    mu_old<-mu
    pi_k_old<-pi_k
    Sigma_old<-Sigma
    W_old<-W_new
    #### Estep ####
    s<-sapply(1:K,function(k){
      sapply(1:n,function(i){
        -1/2 * t(Yvec[i,] - mu[k,]) %*% solve(diag(Sigma[,k])) %*% (Yvec[i,] - mu[k,])
      })+log(pi_k[k]) -1/2 * sum(log(Sigma[,k]))
    })
    s<-exp(s-apply(s,1,max)%*%matrix(1,1,K))
    s<-s/(apply(s,1,sum)%*%matrix(1,1,K))

    #Z = apply(s,1,which.max)

    #### Mstep ####
    ## pik
    pi_k = apply(s,2,mean)
    ## moyenne
    Xcal = kronecker(Ld,s)
    ## Search lamnda opt
    Born_int<-c(0,+Inf)
    mult<-1
    Bad_L<-TRUE
    while(Bad_L){
      reg.gplasso = gglasso(Xcal, Ycal, lambda = lambda_temp, intercept = FALSE,
                            group = rep(1:(K*dold),pold))

      B<-matrix(reg.gplasso$beta,nrow=K)
      mu<-B%*%Ld
      TL_double<-lapply(1:K,function(k){unique(c(0,which(B[k,]!=0)-1,d))})
      ## Affichage
      matplot(t(Yvec),type="l",col=apply(s,1,which.max),main=paste0("lambda=",lambda_temp))
      for (k in 1:K){
        lines(mu[k,],col=3+k,lwd=2,lty=2)
        #abline(v=TL_double[[k]],col=3+k,lty=2,lwd=2)
      }

      TL<-lapply(1:K,function(k){unique(c(0,sort(unique(TL_double[[k]]%%dold)),dold))})
      L_temp<-sort(sapply(1:K,function(k){length(TL[[k]][-c(1,length(TL[[k]]))])}))
      if (any(sort(L_temp)<sort(L))){
        mult<-mult+1
        Born_int[2]<-lambda_temp
        lambda_temp<-mean(Born_int)
      }else if (any(((mult*sort(L))<sort(L_temp)))){
        Born_int[1]<-lambda_temp
        if (Born_int[2]==+Inf){
          lambda_temp<-2*lambda_temp
        }else{
          lambda_temp<-mean(Born_int)
        }
      }else{
        Bad_L<-FALSE
      }
    }
    W_new<-mean(abs(mu_old-mu))+mean(abs(pi_k_old-pi_k))+
      mean(abs(Sigma_old-Sigma))
    if ((abs(W_new-W_old)/W_new)<epsilon){
      break
    }

  }

  Sigma.tab = array(c(Sigma), dim = c(dold, pold, K))
  res = list(mu = mu.tab, Sigma = Sigma.tab, pi = pi_k, Z = Z)
  return(res)
}

####### p > 1 ######
EMalgorithmLasso_div = function(Y, K, niter=100, lambda=1, package =  'gglasso',epsilon=10^(-5)){
  ### lambda= vector of length K
  n = dim(Y)[1]
  dold = dim(Y)[2]
  pold = dim(Y)[3]
  Yvec = matrix(c(Y), nrow = n) #### Mets les lignes bouts a bout qu'importe les dimensions
  d = dim(Yvec)[2]
  ### Initialization
  res.init = Init_Em(Yvec,K) ### TO BE IMPROVED
  mu = res.init
  z = apply(sapply(1:K,function(k){
    sapply(1:n,function(i){
      sum((Yvec[i,]-mu[k,])^2)
    })
  }),1,which.min)
  s<-matrix(sapply(1:n,function(i){
    temp<-rep(0,K)
    temp[z[i]]<-1
    temp
  }),ncol=K,byrow=TRUE)
  Sigma<-sapply(1:K,function(k){
    S_temp<-sum(s[,k])
    unlist(sapply(1:d,function(l){
      sum(s[,k] %*% ((Yvec[,l] - mu[k,l])^2))
    }))/S_temp
  })
  s_col<-colSums(s)
  if (any(s_col<=1)){
    for (s_it in 1:length(s_col)){
      if (s_col[s_it]<=1){
        Sigma[,s_it]<-rchisq(d, df = max(Yvec)-min(Yvec))
      }
    }
  }


  #res.kmeans = kmeans(Xvec, K)
  #mu = t(res.kmeans$centers)
  #Sigma = matrix(rchisq(d*K, df = K*d-1), ncol = K) ### Pourquoi simuler si on peut faire le calcul exact ?
  #loc = runif(K)
  #pi_k =  loc / sum(loc)
  #pi_k = c(0.5,0.5)
  pi_k = apply(s,2,mean)
  #s = matrix(nrow = n, ncol = K)
  ## Necessary
  Ip<-matrix(0,pold,pold)
  diag(Ip)<-1
  Ld = matrix(0,nrow = dold, ncol = dold)
  Ld[lower.tri(Ld, diag=TRUE)] = 1
  Ld<-kronecker(Ip,Ld)
  W_new<-0
  pb<-txtProgressBar(0,niter,style=3)
  Ycal=c(Yvec)
  lambda<-numeric(niter)
  B<-matrix(0,nrow=K,ncol=d)
  lambda_temp<-rep(10^(-255),K)
  for (ite in 1:niter){
    setTxtProgressBar(pb = pb,ite)
    #### old
    mu_old<-mu
    pi_k_old<-pi_k
    Sigma_old<-Sigma
    W_old<-W_new
    #### Estep ####
    s<-sapply(1:K,function(k){
      sapply(1:n,function(i){
        -1/2 * t(Yvec[i,] - mu[k,]) %*% solve(diag(Sigma[,k])) %*% (Yvec[i,] - mu[k,])
      })+log(pi_k[k]) -1/2 * sum(log(Sigma[,k]))
    })
    s<-exp(s-apply(s,1,max)%*%matrix(1,1,K))
    s<-s/(apply(s,1,sum)%*%matrix(1,1,K))

    #Z = apply(s,1,which.max)

    #### Mstep ####
    ## pik
    pi_k = apply(s,2,mean)
    ## moyenne
    for (k in 1:K){
      Ytemp<-(s[,k]%*%matrix(1,nrow=1,ncol=d))*Yvec
      Ytemp<-Ytemp[s[,k]>0.01,]
      Ycal<-c(Ytemp)
      Xcal = kronecker(Ld,s[s[,k]>0.01,k])
      #### lamda
      if (length(Ycal)>0){
        try_lambda<-try(cv.gglasso(Xcal, Ycal, intercept = FALSE,group = rep(1:(dold),pold))$lambda.1se,TRUE)
        if (class(try_lambda)!="try-error"){
          lambda_temp[k]<-try_lambda
        }
        reg.gplasso = gglasso(Xcal, Ycal, lambda = lambda_temp[k], intercept = FALSE,
                              group = rep(1:(dold),pold))
        B[k,]<-matrix(reg.gplasso$beta,nrow=1)
      }else{
        B[k,]=0
      }
    }

    mu<-B%*%t(Ld)
    ## Affichage
    # matplot(t(Yvec),type="l",col=apply(s,1,which.max),main=paste0("lambda=",lambda_temp))
    # for (k in 1:K){
    #   lines(mu[k,],col=3+k,lwd=2,lty=2)
    #   #abline(v=TL_double[[k]],col=3+k,lty=2,lwd=2)
    # }
    ### TL
    TL<-lapply(1:K,function(k){unique(c(0,which(B[k,]!=0)-1,d))})
    ### Sigma
    ## variance
    Sigma<-sqrt(sapply(1:K,function(k){
      L<-length(TL[[k]])-1
      S_temp<-sum(s[,k])
      if (S_temp>0){
        unlist(sapply(1:L,function(l){
          taille<-(TL[[k]][l]+1):TL[[k]][l+1]
          rep(sum(s[,k] %*% ((Yvec[,taille]
                              - matrix(rep(mu[k,taille],n),n,byrow=TRUE))^2))
              ,length(taille))
        }))/S_temp
      }else{
        rchisq(d, df = max(Yvec)-min(Yvec))
      }
    }))

    W_new<-mean(abs(mu_old-mu))+mean(abs(pi_k_old-pi_k))+
      mean(abs(Sigma_old-Sigma))
    if (W_new<epsilon){
      break
    }
  }

  TL<-lapply(1:K,function(k){unique(sort(c(0,(TL[[k]]-1)%%dold+1,dold)))})
  ord<-order(sapply(1:K,function(k){length(TL[[k]])}))
  ### On ordonne par ordre croissant
  TL<-lapply(1:K,function(k){TL[[ord[k]]]})
  mu<-mu[ord,]
  Sigma<-Sigma[,ord]
  pi_k<-pi_k[ord]
  ### Vrais
  Q<-sum(sapply(1:K,function(k){
    sum(sapply(1:n,function(i){
      -1/2 * t(Yvec[i,] - mu[k,]) %*% solve(diag(Sigma[,k])) %*% (Yvec[i,] - mu[k,])
    })+log(pi_k[k]) -1/2 * sum(log(Sigma[,k])))
  }))
  ### On met en forme
  mu<-lapply(1:K,function(k){
    sapply(1:pold,function(r){mu[k,TL[[k]][-length(TL[[k]])]+1+(r-1)*dold]})
  })

  ### Renvoi
  new(Class="MixSegTime",mu=mu,
      sigma=Sigma,TL=TL,
      pik=pi_k,method="LASSO",z=apply(s[,ord],1,which.max),Q=Q)
}



