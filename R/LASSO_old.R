library(fields)
library(grplasso)
library(tictoc)
library(gglasso)
library(glmnet)
#EMalgorithmLasso(x[,,1], 2, 10, rep(0.0001,2))

###### p=1 #####

EMalgorithmLasso = function(X, K, lambda=NULL, niter=10){
  ### lambda= vector of length K
  n = dim(X)[1]
  d = dim(X)[2]
  #p = dim(X)[3]
  p=1

  ### Initialization
  res.kmeans = kmeans(X, K) ### TO BE IMPROVED
  mu = t(res.kmeans$centers)
  Sigma = matrix(rchisq(d*K, df = K*d-1), ncol = K)
  loc = rep(1,K)
  pi_k =  loc / sum(loc)
  s = matrix(nrow = n, ncol = K)
  num = rep(0,K)
  ite = 0
  while ((ite < niter)){
    print(ite)
    ite = ite +1
    #### Estep ####
    for (i in 1:n){
      for (k in 1:K){
        num[k] = log(pi_k[k]) -1/2 * sum(log(Sigma[,k])) -1/2 * t(X[i,] - mu[,k]) %*% solve(diag(Sigma[,k])) %*% (X[i,] - mu[,k])
      }
      num = num - max(num)
      num = exp(num)
      s[i,] = num / sum(num)
    }
    Z = apply(s,1,which.max)

    #### Mstep ####
    ## pik
    pi_k = apply(s,2,mean)
    ## variance
    for (k in 1:K){
      Sigma[,k] = s[,k] %*% ((X[,] - mu[,k])^2) / sum(s[,k])
    }
    ## moyenne
    covariates = matrix(0,nrow = d, ncol = d)
    covariates[upper.tri(covariates, diag=TRUE)] = 1
    for (k in 1:K){
      nk = sum(Z==k)
      kro.covariates = kronecker(covariates,matrix(1,nrow=  1, ncol = nk))
      Xk = X[Z==k,]
      reg.lasso = glmnet(t(kro.covariates), c(Xk), lambda = lambda[k], intercept = FALSE)
      mu[,k] = cumsum(reg.lasso$beta)
    }
    matplot(mu, type='l')
  }
  res = list(mu = mu, Sigma = Sigma, pi = pi_k)
  return(res)
}

####### p > 1 ######
EMalgorithmLasso = function(X, K, niter, lambda, package =  'gglasso'){
  ### lambda= vector of length K
  n = dim(X)[1]
  dold = dim(X)[2]
  pold = dim(X)[3]
  Xvec = matrix(c(X), nrow = n) #### Mets les lignes bouts a bout qu'importe les dimensions
  d = dim(Xvec)[2]
  ### Initialization
  res.kmeans = kmeans(Xvec, K)
  mu = t(res.kmeans$centers)
  Sigma = matrix(rchisq(d*K, df = K*d-1), ncol = K) ### Pourquoi simuler si on peut faire le calcul exact ?
  #loc = runif(K)
  #pi_k =  loc / sum(loc)
  #pi_k = c(0.5,0.5)
  pi_k = sapply(1:K, function(k) {
    mean(res.kmeans$cluster == k)
  })
  s = matrix(nrow = n, ncol = K)
  num = rep(0,K)
  ite = 0
  mu.tab = array(0, dim = c(dold, pold, K))
  while ((ite < niter)){
    print(ite)
    ite = ite +1
    #### Estep ####
    for (i in 1:n){
      num = sapply(1:K, function(k) {
        log(pi_k[k]) -1/2 * sum(log(Sigma[,k])) -1/2 * t(Xvec[i,] - mu[,k]) %*% solve(diag(Sigma[,k])) %*% (Xvec[i,] - mu[,k]) })
      num = num - max(num)
      num = exp(num)
      s[i,] = num / sum(num)
    }
    Z = apply(s,1,which.max)
    print(Z)
    print('E step')
    #### Mstep ####
    ## pik
    pi_k = apply(s,2,mean)
    ## variance ### TO BE MODIFIED si on veut que ce soit constant par segment !!
    Sigma = sapply(1:K, function(k){ s[,k] %*% ((Xvec - mu[,k])^2) / sum(s[,k])})
    print('pik et variance done')
    ## moyenne
    covariates = matrix(0,nrow = d, ncol = d)
    covariates[upper.tri(covariates, diag=TRUE)] = 1
    for (k in 1:K){
      nk = sum(Z==k)
      kro.covariates = kronecker(covariates,matrix(1,nrow=  1, ncol = nk))
      Xveck = Xvec[Z==k,]
      if ( package ==  'grpLasso'){
        reg.gplasso = grplasso(x = t(kro.covariates), y = c(Xveck), lambda = lambda[k], intercept = FALSE, index = rep(1:dold, pold), model = LinReg())
        param = reg.gplasso$coefficients
      } else {
        reg.gplasso = gglasso(t(kro.covariates), c(Xveck), lambda = lambda[k], intercept = FALSE,
                              group = rep(1:pold, dold))
        param = reg.gplasso$beta
      }
      mu.tab[,,k] = matrix(param, ncol = pold, byrow = TRUE)
      var.sec =  which(mu.tab[,1,k] != 0)
      mu.tab[,,k] = apply(mu.tab[,,k], 2, cumsum)
      mu[,k] = c(mu.tab[,,k])
    }
    print('moyenne done')
  }

  Sigma.tab = array(c(Sigma), dim = c(dold, pold, K))
  res = list(mu = mu.tab, Sigma = Sigma.tab, pi = pi_k, Z = Z)
  return(res)
}
test = FALSE
if (test == TRUE){

  source('~/ownCloud/Recherche/_ Projets en cours/Mixture of Segmentations/Simulation.R')
  dim(x)
  x = x[1:40,1:50,]
  res = EMalgorithmLasso(x, 2, 3, rep(0.1,2), package =  'grpLasso')
  res = EMalgorithmLasso(x, 2, 3, rep(0.01,2), package =  'gglasso')
  matplot((res$mu[,,1]), type='l')
  matplot(res$mu[,,2], type='l')

  ### test gglasso
  p = 20000
  n = 500
  library(mvtnorm)
  X = matrix(rnorm(n*p, 0,1), nrow = n)
  beta = c(rep(1,5), rep(0, p-5))
  Y =  X %*% beta + rnorm(n)
  reg.gplasso = gglasso(X, Y, lambda = lambda[k], intercept = FALSE, group = rep(1:5, p/5))
}
#### il manque: le choix de lambda, le choix de K,
#### comparaison des différents packages de group lasso pour trouver le plus rapide,
####
# n = 40
# d = 50
# dvec = d * 2
# covariates = matrix(0,nrow = dvec, ncol = dvec)
# covariates[upper.tri(covariates, diag=TRUE)] = 1
# data = simul_data(n, d, 2, 2, rep(2,2), 'facile')
# X = data$x
# kro.covariates = kronecker(covariates,matrix(1,nrow=  1, ncol = n))
# Xvec = matrix(c(X),nrow = n)
# reg.gplasso = gglasso(t(kro.covariates), c(Xvec), lambda = 0.01, intercept = FALSE,
#                       group = rep(1:d, 2))
# param = reg.gplasso$beta
# tic()
# par(mfrow=c(1,2))
# res = EMalgorithmLasso(X, 2, 3, rep(3,2), package =  'grpLasso')
# toc()
#
