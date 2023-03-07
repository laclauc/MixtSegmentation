##########################
# Algorithme EM
# MAJ 05/01/2023
# vincent.brault@univ-grenoble-alpes.fr
###########################

Deltak<-function(x,sk,Moy){
  # Parametres
  n=dim(x)[1]
  d=dim(x)[2]
  p=dim(x)[3]
  temp<-matrix(0,dim(x)[2],dim(x)[2])
  splus=sum(sk)
  splus=splus*(splus>0)+10^(-250)*(splus<=0)
  ## Compute mukstbar for t<s Attention au decalage 1->0 pour t
  Carre<-array(0,c(p,d,d))
  for (r in 1:p){
    Carre[r,,]<-sapply(1:d,function(s){
      sapply(1:d,function(t){
        if (s<t){
          return(0)
        }else{
          return(Moy[r,,t,s]%*%sk/splus)
        }})
    })
  }
  ### Apres affichage carre semble coherent
  sapply(1:d,function(s){
    sapply(1:d,function(t){
      if (s<t){
        return(0)
      }else if (t==s){
        return(sum(sapply(1:p,function(r){
          sk%*%((x[,t,r]-Carre[r,t,t])^2)
        })))
      }else{
        return(sum(sapply(1:p,function(r){
          sk%*%(rowSums((x[,t:s,r]-Carre[r,t,s])^2))
        })))
      }
    })
  })

}

ProgDyn<-function(Delta,L){
  d<-nrow(Delta)
  Ip<-matrix(0,L+1,d)
  ind<-matrix(0,L+1,d)
  Ip[1,]<-Delta[1,]
  for (l in 2:(L+1)){
    ind[l,l:d]<-sapply(l:d,function(q){which.min(Ip[l-1,l:(q-1)]+Delta[(l+1):q,q])})+l-1
    Ip[l,l:d]<-sapply(l:d,function(q){(Ip[l-1,ind[l,q]]+Delta[ind[l,q]+1,q])})
  }
  #### Recuperation des resultats
  # Ruptures
  TL=numeric(L+1)
  TL[L+1]=ind[L+1,d]
  for (l in L:1){
    TL[l]=ind[l,TL[l+1]]
  }
  list(TL=c(TL,d),S=Ip[L+1,d])
}


estim_mu<-function(Moy,TL,sk){
  p<-dim(Moy)[1]
  splus=sum(sk)
  sapply(1:p,function(r){sapply(1:(length(TL)-1),function(l){
    Moy[r,,TL[l]+1,TL[l+1]]%*%sk/splus
  })})
}

Calcul_Moy<-function(x){
  # Verification et init
  if (!is.array(x)){
    stop("x must be an array")
  }
  n=dim(x)[1]
  d=dim(x)[2]
  p=dim(x)[3]
  # Calculs
  cat("Calcul des valeurs moyennes (Chronophage)\n")
  Moy<-array(0,c(p,n,d,d))
  pb<-txtProgressBar(0,n*p,style = 3)
  for (r in 1:p){
    for (i in 1:n){
      setTxtProgressBar(pb,i+(r-1)*n)
      for (t in 1:d){
        for (s in t:d){
          Moy[r,i,t,s]<-mean(x[i,t:s,r])
        }
      }
    }
  }
  Moy
  #### Complexity : pnd^2
}

EM_Dyn<-function(x,K,L,niter=100,epsilon=10^(-5),Moy=NULL,mc.cores=1,Init=NULL,Essai_max=100,
             Essai_Min=10){
  #############
  # Initialisation
  #############
  # Verification
  if (!is.array(x)){
    stop("x must be an array")
  }
  if ((!is.numeric(K))|(K<=0)|(floor(K)!=K)){
    stop("K must be an positive integer")
  }
  if (is.numeric(L)){
    if (length(L)==1){
      L<-rep(L,K)
    }else if((any(L<=0))|(any(floor(L)!=L))|(length(L)!=K)){
      stop("L must be an positive integer or a positive vector integer of size K")
    }
  }else{
    stop("L must be an positive integer or a positive vector integer of size K")
  }
  if (!is.numeric(mc.cores)){
    stop("max.break must be an integer between 1 and n")
  } else if ((mc.cores<=0)||(length(mc.cores)!=1)||(floor(mc.cores)!=mc.cores)){
    stop("max.break must be a positive integer")
  }

  ## =============================================================
  ## INITIALIZATION & PARAMETERS RECOVERY
  if (mc.cores>1){
    if (Sys.info()[['sysname']] == "Windows") {
      warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
      mc.cores <- 1
    }
  }
  # Parametres
  n=dim(x)[1]
  d=dim(x)[2]
  p=dim(x)[3]
  # Les valeurs moyennes
  if (is.null(Moy)){
    Moy<-Calcul_Moy(x)
  }
  ### Au cas ou classes vides
  ref=Inf
  Vide=TRUE
  Essai=1
  while(((Vide)&(Essai<=Essai_max))|(Essai<Essai_Min)){
    # Init des parametres par etape M
    pik<-rbeta(K,4,4)
    pik<-pik/sum(pik)
    if (is.null(Init)||(Essai>1)){
      z<-sample(1:K,size = n,replace = TRUE,prob = pik)
    }else{
      z<-Init
    }
    sik=matrix(0,n,K)
    for (k in unique(z)){
      sik[z==k,k]=1
    }
    if (mc.cores>1){
      temp<- parallel::mclapply(1:K,function(k){ProgDyn(Deltak(x,sik[,k],Moy),L[k])}, mc.cores=mc.cores)
      mu<- parallel::mclapply(1:K,function(k){estim_mu(Moy,temp[[k]]$TL,sik[,k])}, mc.cores=mc.cores)
    } else {
      temp<-lapply(1:K,function(k){ProgDyn(Deltak(x,sik[,k],Moy),L[k])})
      mu<-lapply(1:K,function(k){estim_mu(Moy,temp[[k]]$TL,sik[,k])})
    }

    colS<-colSums(sik)
    sigma<-sapply(1:K,function(k){temp[[k]]$S/(p*d)})/(colS+10^(-250)*(colS<=0))
    pik<-colS/n
    ###########################
    # Debut iterations
    ###########################
    cat("\n Lancement numero ",Essai," de la boucle : si 100% -> toutes les iterations faites sans sortir\n")
    pb<-txtProgressBar(0,niter,style = 3)
    old=sum(sapply(1:K,function(k){temp[[k]]$S}))+sum(colSums(sik)*log(pik+10^(-250)*(pik<=0)))
    for (iter in 1:niter){
      setTxtProgressBar(pb,iter)
      #########################
      # Step E
      #########################
      lsik<-sapply(1:K,function(k){
        Nb<-diff(temp[[k]]$TL)
        if (sigma[k]>0){
        -apply((x-array(unlist(sapply(1:p,function(r){
          sapply(1:length(Nb),function(l){
            rep(mu[[k]][l,r],Nb[l]*n)})})),
          c(n,d,p)))^2,1,sum)/(sigma[k])
        }else{
         rep( -10^(250),n)
        }
      })-matrix(1,nrow=n,1)%*%(log(sigma*(sigma>0)+10^(-250)*(sigma<=0)))### 1/2 et log 2pi constants donc enleves
      # substraction of the max to have at least one positive value after exp
      lsik<-lsik-apply(lsik,MARGIN = 1,max)%*%matrix(1,nrow=1,K)
      # exponentielle
      sik<-exp(lsik)
      # normalisation
      sik<-sik/(rowSums(sik)%*%matrix(1,nrow=1,K))
      ##########################
      # Step M
      #########################
      if (mc.cores>1){
        temp<- parallel::mclapply(1:K,function(k){ProgDyn(Deltak(x,sik[,k],Moy),L[k])}, mc.cores=mc.cores)
        mu<- parallel::mclapply(1:K,function(k){estim_mu(Moy,temp[[k]]$TL,sik[,k])}, mc.cores=mc.cores)
      } else {
        temp<-lapply(1:K,function(k){ProgDyn(Deltak(x,sik[,k],Moy),L[k])})
        mu<-lapply(1:K,function(k){estim_mu(Moy,temp[[k]]$TL,sik[,k])})
      }

      colS<-colSums(sik)
      sigma<-sapply(1:K,function(k){temp[[k]]$S/(p*d)})/(colS+10^(-250)*(colS<=0))
      pik<-colS/n
      ##########################
      # Break ?
      ##########################
      new<-sum(sapply(1:K,function(k){temp[[k]]$S}))+sum(colSums(sik)*log(pik+10^(-250)*(pik<=0)))
      if (abs(new-old)/old<epsilon){
        break
      }else{
        old<-new
      }
    }
    z<-apply(sik,1,which.max)
    cat("\n Q: ",new," with (",paste0(names(table(z)),collapse =","),")=(",paste0(table(z),collapse =","),")\n")
    if (length(unique(z))<K){ ### Empty cluster
      if ((ref>new)&(Vide)){ ### If there exist a partition without empty cluster (vide==FALSE), we forgot thise proposition
        K_part<-length(unique(z))
        cluster_part<-sort(unique(z))
        for (kbis in 1:K_part){
          mu[[kbis]]<-mu[[cluster_part[kbis]]]
          sigma[kbis]<-sigma[cluster_part[kbis]]
          temp[[kbis]]<-temp[[cluster_part[kbis]]]
          pik[kbis]<-pik[cluster_part[kbis]]
          z[z==cluster_part[kbis]]<-kbis
        }
        mu<-lapply(1:K_part,function(it){mu[[it]]})
        sigma<-sigma[1:K_part]
        temp<-lapply(1:K_part,function(it){temp[[it]]})
        pik<-pik[1:K_part]
        pik<-pik/sum(pik)
        res_part=methods::new(Class="MixSegTime",mu=mu,
                         sigma=as.matrix(sigma),TL=lapply(1:length(temp),function(k){temp[[k]]$TL}),
                         pik=pik,method="EM-Dyn",z=z,Q=new)
        ref=new
      }
    }else{ ### Any empty cluster
      if (Vide){ ### It is the first case without empty cluster
        res_part<-methods::new(Class="MixSegTime",mu=mu,
                               sigma=as.matrix(sigma),TL=lapply(1:length(temp),function(k){temp[[k]]$TL}),
                               pik=pik,method="EM-Dyn",z=z,Q=new)
        ref=new
      }
      Vide=FALSE ### We have a solution without empty cluster
      if (ref>new){ ### If it is not the first proposition
        res_part<-methods::new(Class="MixSegTime",mu=mu,
                               sigma=as.matrix(sigma),TL=lapply(1:length(temp),function(k){temp[[k]]$TL}),
                               pik=pik,method="EM-Dyn",z=z,Q=new)
        ref=new
      }
    }
    Essai=Essai+1

  }


  if (Vide){
    warning("Empty cluster")
  }
  return(res_part)
}







