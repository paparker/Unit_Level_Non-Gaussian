library(BayesLogit)
library(Matrix)
library(mvtnorm)
library(LaplacesDemon)

lse <- function(x){
  m <- max(x)
  m + log(sum(exp(x-m)))
}

softmax <- function(x){
  exp(x - lse(x))
}

PLMBmcmc <- function(X, Psi, Y, sig2b=1000, iter=1000, burn=500, wgt=NULL){
  ## This function fits the PL-MB model via Gibbs sampling
    # X is the input covariate matrix, one row per sample unit
    # Psi is the matrix of spatial basis functions, one row per sample unit
    # Y is a vector of binary survey responses
    # sig2B is the prior variance for Beta
    # iter is the number of iterations
    # burn is the length of burn in
    # wgt is the vector of survey weights
  p <- ncol(X)
  r <- ncol(Psi)
  n <- length(Y)
  if(is.null(wgt)) wgt <- rep(1,n)
  Binv <- (1/sig2b)*Diagonal(p)
  kappa <- wgt*(Y - 0.5)
  w <- rep(1, n)
  beta <- rep(0, p)
  eta <- rep(0, r)
  sig2E <- 1
  betaOut <- matrix(NA, nrow=iter, ncol=p)
  etaOut <- matrix(NA, nrow=iter, ncol=r)
  sig2Eout <- rep(NA, iter)
  print(paste0("Starting ",iter, " iterations."))
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## Sample fixed effects
    precBeta <- solve(t(X)%*%Diagonal(length(w),w)%*%X + Binv)
    meanBeta <- t(X)%*%Diagonal(length(w),w)%*%(kappa/w - Psi%*%eta )
    beta <- betaOut[i,]  <- as.numeric(rmvnorm(1, mean=precBeta%*%meanBeta, sigma=as.matrix(precBeta)))
    
    ## Sample random effects
    Einv <- (1/sig2E)*Diagonal(r)
    precEta <- solve(t(Psi)%*%Diagonal(length(w),w)%*%Psi + Einv)
    meanEta <- t(Psi)%*%Diagonal(length(w),w)%*%(kappa/w - X%*%beta )
    eT <- as.numeric(rmvnorm(1, mean=precEta%*%meanEta, sigma=as.matrix(precEta)))
    eta <- etaOut[i,]  <- eT - mean(eT)
    
    
    ## Sample RE variance
    sig2E <- sig2Eout[i] <- 1/rgamma(1,
                                     shape=0.5 + r/2,
                                     (0.5 + 0.5*t(eta)%*%(eta)))
    
    ## Sample latent PG variables
    w <- rpg(n, wgt, as.numeric(X%*%beta + Psi%*%eta))
    setTxtProgressBar(pb, i)
  }
  return(list(Beta=betaOut[-c(1:burn),], Eta=etaOut[-c(1:burn),], sig2E=sig2Eout[-c(1:burn)]))
}

PLMMmcmc <- function(X, Psi, Y, N=NULL, sig2b=1000,  iter=1000, burn=500, wgt=NULL, predX=NULL, predPsi=NULL){
  ## This function fits the PL-MM model via variational Bayes and multinomial stick breaking
  # X is the input covariate matrix, one row per sample unit
  # Psi is the matrix of spatial basis functions, one row per sample unit
  # Y is a vector of binary survey responses
  # N is a vector of multinomial sizes, defaults to 1 when left as NULL
  # sig2B is the prior variance for Beta
  # iter is the number of iterations
  # burn is the length of burn in
  # wgt is the vector of survey weights
  # predX is the covariate matrix for predictions (poststrat cells in our case)
  # predPsi is the matrix of spatial basis functions for predictions
  # nsamp is the number of samples from variational posterior for inference
  
  n <- nrow(Y)
  if(is.null(N)) N <- rep(1, n)
  if(is.null(wgt)) wgt <- rep(1, n)
  
  
  K <- ncol(Y)
  Nstar <- cbind(1, N - t(apply(Y, 1, cumsum)))[,1:(K-1)]
  wgt <- rep(wgt,K-1)
  preds <- array(NA, dim=c(nrow(predX), K-1, iter-burn))
  for(kk in 1:(K-1)){
    nn <- Nstar[,kk]
    keep <- which(nn>0)
    yy <- Y[keep, kk]
    xx <- X[keep, ]
    zz <- Psi[keep, ]
    ww <- wgt[keep]
    n2 <- length(yy)
    
    tt <- PLMBmcmc(X=xx, Psi=zz, Y=yy, sig2b=sig2b, iter=iter, burn=burn, wgt=ww)
    preds[,kk,] <- plogis(as.matrix(cbind(predX, predPsi)%*%t(cbind(tt$Beta, tt$Eta))))
  }
  preds <- aperm(apply(preds,c(1,3), Stick), c(2,1,3))
  return(list(Preds=preds))
}


PLMBvb <- function(X, Psi, Y, sig2b=1000,  eps=0.001, wgt=NULL){
  ## This function fits the PL-MB model via variational Bayes
    # X is the input covariate matrix, one row per sample unit
    # Psi is the matrix of spatial basis functions, one row per sample unit
    # Y is a vector of binary survey responses
    # sig2B is the prior variance for Beta
    # eps is convergence tolerance
    # wgt is the vector of survey weights
  Au <- Bu <- 0.5
  n <- length(Y)
  if(is.null(wgt)) wgt <- rep(1, n)
  p <- ncol(X)
  r <- ncol(Psi)
  C <- cbind(X,Psi)
  Bsig2u <- 1
  MUbu <- rep(1,p+r)
  SIGbu <- Diagonal(p+r)
  checkOld <- Inf
  checkVec <- c()
  iter <- 1
  repeat{
    ## Latent Variables
    xi <- as.numeric(sqrt(colSums(t(C)*(SIGbu%*%t(C))) + (C%*%MUbu)^2))
    
    ## Regression coefficients
    Zbar <- Diagonal(n,(wgt*0.5/xi)*tanh(0.5*xi))
    SIGbu <- solve(bdiag((1/sig2b)*Diagonal(p), 
                         as.numeric((Au+r/2)/(Bsig2u))*Diagonal(r)) +
                     t(C)%*%Zbar%*%C)
    MUbu <-  SIGbu %*% t(C) %*% (wgt*(Y - 0.5))
    
    
    
    ## RE Variance
    Bsig2u <- Bu + 0.5*(t(MUbu[-c(1:p)])%*%(MUbu[-c(1:p)]) + sum(diag(SIGbu[-c(1:p),-c(1:p)])))
    
    PPrec <- bdiag((1/sig2b)*Diagonal(p), 
                   as.numeric((Au+r/2)/(Bsig2u))*Diagonal(r))
    
    ## Check for convergence
    checkNew  <- 0.5*(p+r) + 0.5*determinant(SIGbu, logarithm = T)$modulus  - 
      0.5*t(MUbu)%*%PPrec%*%(MUbu) + sum(wgt*(Y-0.5)*as.numeric(C%*%MUbu) + wgt*log(plogis(xi)) - 0.5*wgt*xi) - 
      0.5*sum(diag(PPrec %*% SIGbu)) - log(Bsig2u)
    checkVec <- c(checkVec, as.numeric(checkNew))
    
    
    if (as.logical(abs(checkOld - checkNew) < eps)) break
    iter <- iter + 1
    checkOld <- checkNew
    print(checkNew)
  }
  return(list(iter=iter, BU=list(type="Normal", mean=MUbu, Cov=SIGbu, p=p, r=r),
              Bsig2u=list(type="Inv. Gamma", A=Au + r/2, B=Bsig2u), Elbo=checkVec))
}

PLMMvb <- function(X, Psi, Y, N=NULL, sig2b=1000,  eps=0.001, wgt=NULL, predX=NULL, predPsi=NULL, nsamp=1000){
  ## This function fits the PL-MM model via variational Bayes and multinomial stick breaking
    # X is the input covariate matrix, one row per sample unit
    # Psi is the matrix of spatial basis functions, one row per sample unit
    # Y is a vector of binary survey responses
    # N is a vector of multinomial sizes, defaults to 1 when left as NULL
    # sig2B is the prior variance for Beta
    # eps is convergence tolerance
    # wgt is the vector of survey weights
    # predX is the covariate matrix for predictions (poststrat cells in our case)
    # predPsi is the matrix of spatial basis functions for predictions
    # nsamp is the number of samples from variational posterior for inference
  
  n <- nrow(Y)
  if(is.null(N)) N <- rep(1, n)
  if(is.null(wgt)) wgt <- rep(1, n)
  
  
  K <- ncol(Y)
  Nstar <- cbind(1, N - t(apply(Y, 1, cumsum)))[,1:(K-1)]
  wgt <- rep(wgt,K-1)
  preds <- array(NA, dim=c(nrow(predX), K-1, nsamp))
  for(kk in 1:(K-1)){
    nn <- Nstar[,kk]
    keep <- which(nn>0)
    yy <- Y[keep, kk]
    xx <- X[keep, ]
    zz <- Psi[keep, ]
    ww <- wgt[keep]
    n2 <- length(yy)
    
    tt <- PLMBvb(X=xx, Psi=zz, Y=yy, sig2b=sig2b, eps=eps/K, wgt=ww)
    predTheta <- t(rmvnorm(nsamp, mean=tt$BU$mean, sigma=as.matrix(tt$BU$Cov)))
    preds[,kk,] <- plogis(as.matrix(cbind(predX, predPsi)%*%predTheta))
  }
  preds <- aperm(apply(preds,c(1,3), Stick), c(2,1,3))
  return(list(Preds=preds))
}