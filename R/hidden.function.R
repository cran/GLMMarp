########################################################################
# File: hidden.function.R                                                            
# Author: Xun Pang                                                             
# Creation Date: 17/09/2009                                                            
# Description: 	functions to support the MCMC simulation file GLMMARp.R
########################################################################

 # load(library)
 library(bayesSurv); library(panel);
 library(MCMCpack);  library(mvtnorm);  library(lattice)
 library(msm); library(MASS);  library(Matrix)
 library(accuracy); library(kinship)

 #########  Functions #########
 AcceptRate0 <- function(rho, orho,precision=1, unit, AY,
                         AX, AW, AS, Abeta, Ab, AT, Ac, Entry){
   fphi <- 0
   ofphi <- 0
   
   for (n in 1: max(unit)){
     roww <- Entry[which(unit==n)]
     T<-AT[roww]
     t<-min(T)
     P<-roww[1]
     ystari<-AY[P]
     Xi<- t(as.matrix(AX[P,]))
     Wi<- t(as.matrix(AW[P,]))
     Si <- t(as.matrix(AS[P,]))
     bbbn<- as.matrix(Ab[n,])
     ccn <- as.matrix(Ac[t,])
     rant <- Si%*%cbind(ccn)

     fphi <- fphi+1/2*log(precision*(1-rho^2))-precision*(1-rho^2)/2*(ystari-Xi%*%cbind(Abeta)-Wi%*%cbind(bbbn)-rant)^2
     
     ofphi <- ofphi+1/2*log(precision*(1-orho^2))-precision*(1-orho^2)/2*(ystari-Xi%*%cbind(Abeta)-Wi%*%cbind(bbbn)-rant)^2
   }
     ppalpha<-fphi-ofphi
     palpha<- exp(ppalpha)
     alpha <- min(palpha, 1)
     return(alpha)
 }

  ###################################
  AcceptRateARpT <- function(nrho, orho, Osig, unit, AY, AX, AW, AS, Abeta, 
                             precision=1, Ab, AT, Ac, Entry, f, nlag){
    II <- diag(nlag^2)
    hh <- rep(0, (nlag-1))
    e <- c(1, hh)
    newphi<-nrho
    oldphi<-orho
    pG<-matrix(NA, ncol=nlag, nrow=nlag)
    pG[2:nlag, 1:(nlag-1)] <- diag(nlag-1)
    pG[1,] <- newphi
    pG[-1, nlag] <- 0
     
    ll<-II-kronecker(pG,pG)
    pd<-as.vector(solve(ll)%*%as.vector(e%*%t(e)))
   
    newsig<-1/precision* matrix(pd, nrow=nlag)
    nSig <- det(newsig)
    oSig <- det(Osig)
    fphi1 <- 0
    ofphi1 <-0
    fphi <- 0
    ofphi <-0
    for (n in 1: max(unit)){
    roww <- Entry[which(unit==n)]
    T<-AT[roww]
    t<-min(T)
    Ti <- length(roww)
    if (Ti==1){
      fphi <- fphi+1/2*log(newsig[1,1])-(AY[roww]-AX[roww,]%*%Abeta-AW[roww,]%*%Ab[n,]-AS[roww,]%*%Ac[t,])^2/(2*newsig[1,1])
      ofphi <- ofphi + 1/2*log(Osig[1,1])-(AY[roww]-AX[roww,]%*% Abeta-AW[roww,]%*%Ab[n,]-AS[roww,]%*%Ac[t,])^2/(2*Osig[1,1])
    }else{
      Ji<-min(Ti, nlag)
      P<-roww[c(1:Ji)]
      ystari<-AY[P]
      Xi<- as.matrix(AX[P,])
      Wi<- as.matrix(AW[P,])
      Si <- as.matrix(AS[P,])
      bbbn<-as.matrix(Ab[n,])
      nnSig <- newsig[1:Ji, 1:Ji]
      ooSig <- Osig[1:Ji, 1:Ji]
      dnSig <- det(nnSig)
      doSig <- det(ooSig)
      bbbn<-as.matrix(Ab[n,])
      ccn <- as.matrix(cbind(Ac[ c(t:(t+Ji-1)),]))
      rant <- c()
      for (k in 1:Ji){
        rant[k] <- Si[k,]%*%cbind(ccn[k,])
      }
      fphi <- fphi-1/2*log(dnSig)+ss(ff=nnSig, Yb=ystari, Xb=cbind(Xi, Wi, rant), bb=c(Abeta, bbbn, 1))
      ofphi <- ofphi -1/2*log(doSig)+ss(ff=ooSig, Yb=ystari, Xb=cbind(Xi, Wi, rant), bb=c(Abeta, bbbn, 1))
    }
  }
  ppalpha<-fphi-ofphi
  palpha<- exp(ppalpha)
  alpha <- min(palpha, 1)
  return(alpha)
}

 #######################################
 BTnormlog <- function(x, mu, sigma=1){
   if(x==0) {
      u <- runif(1)
      xstar <- qnorm(log(u)+pnorm(0, mu, sigma, log=TRUE), mu, sigma, log=TRUE)
   } else{
      u <- runif(1)
      xstar <- -qnorm(log(u)+pnorm(0, -mu, sigma, log=TRUE), -mu, sigma, log=TRUE)
   }
   xstar
 }

 ######################################
  check.mcmc.parameters <- function(burnin, mcmc, thin) {
  
    if(mcmc %% thin != 0) {
      cat("Error: MCMC iterations not evenly divisible by thinning interval.\n")
      stop("Please respecify and call GLMMARp.Binary() again.",
         call.=FALSE)
    }
    if(mcmc < 0) {
      cat("Error: MCMC iterations negative.\n")
      stop("Please respecify and call GLMMARp.Binary() again.",
         call.=FALSE) 
    }
    if(burnin < 0) {
      cat("Error: Burnin iterations negative.\n")
      stop("Please respecify and call GLMMARp.Binary() again.",
         call.=FALSE)
    }
    if(thin < 0) {
      cat("Error: Thinning interval negative.\n")
      stop("Please respecify and call GLMMARp.Binary() again.",
         call.=FALSE)
    }    
    return(0)
  }

 #####################################
 Cov.AR1 <- function(rho, col, sigma){
   cov.vector <- sigma
   for (k in 2: col){
    cov.vector[k] <- rho*cov.vector[k-1]
   }
   covariance <- matrix(NA, nrow=col, ncol=col)
   covariance[1,1] <- sigma
   for (i in 1:col){
     for (j in 1:col){
      covariance[i, j] <- cov.vector[abs(i-j)+1]
     }
   }
   covariance
 }

 ####################################
 Cov.AR2 <- function(phig, col){
   II<-diag(4)
   e<-matrix(c(1,0), nrow=2)
   G<-matrix(c(phig[1],1, phig[2], 0), ncol=2)
   d<-as.vector(pseudoinverse(II-kronecker(G,G))%*%as.vector(e%*%t(e)))
   Sig<-matrix(d, nrow=2)
   var <- Sig[1,1]
   cov.vector <- c(Sig[1,1], Sig[1,2])
   for (k in 3: col){
     cov.vector[k] <- phig[1]*cov.vector[k-1]+phig[2]*cov.vector[k-2]
   }
   covariance <- matrix(NA, nrow=col, ncol=col)
   covariance[1:2, 1:2] <- Sig
   for (i in 1:col){
     for (j in 1:col){
      covariance[i, j] <- cov.vector[abs(i-j)+1]
     }
   }
   return(list(Sig, covariance))
 }

 ######################################
 Cov.ARp <- function(phig, col){
   ord <- length(phig)
   II<-diag(ord^2)
   om <- rep(0, (ord-1))
   e<-matrix(c(1,om), nrow=ord)
   G<-matrix(NA, ncol=ord, nrow=ord)
   G[2:ord, 1:(ord-1)] <- diag(ord-1)
   G[1,] <- phig
   G[-1, ord] <- 0

   d<-as.vector(pseudoinverse(II-kronecker(G,G))%*%as.vector(e%*%t(e)))
   Sig<-matrix(d, nrow=ord)
   var <- Sig[1,1]
   cov.vector <- as.vector(Sig[,1])
   for (k in (ord+1): col){
     cov.vector[k] <- rev(phig) %*% cov.vector[(k-ord): (k-1)]
   }
   covariance <- matrix(NA, nrow=col, ncol=col)
   covariance[1:ord, 1:ord] <- Sig
   for (i in 1:col){
     for (j in 1:col){
       covariance[i, j] <- cov.vector[abs(i-j)+1]
     }
   }
   return(list(Sig, covariance))
 }
  
 ######################################
 Covariance.update <- function(Covariate, prior.freedom, prior.scale,
                               dimension, obs.length, current.value){
   if (ncol(Covariate)==1){                                    
     Ei<-rgamma(1,(prior.freedom+obs.length)/2,
                rate=(prior.scale+sum(current.value^2))/2)
     Ei <- as.matrix(Ei)
     return(list(Ei=Ei))
   }else{
     OEin <- solve(prior.scale)+sumtrans(Mat=current.value)
     sE<-rWishart (1, prior.freedom+obs.length, solve(OEin))
     Ei <- Wishart(element=sE, col=dimension)
  }
     return(list(Ei=Ei, sE=sE))
 }

 #####################################
 CutPointUpdate <- function(obs, latent, oldcut, cut0){
   ntau<-c()
   catnum <- length(oldcut)+2       
   if (catnum==3){
     al<-max(max(latent[which(obs==2)]), cut0)
     bl<-min(latent[which(obs==3)])
     ntau<-runif(1, al, bl)
   }
   if (catnum==4){
     a1<-max(max(latent[which(obs==2)]), cut0) 
     b1<-min(min(latent[which(obs==3)]), oldcut[2])
     ntau[1]<-runif(1, a1, b1)
     al<-max(max(latent[which(obs==3)]), ntau[1])
     bl<-min(latent[which(obs==4)])
     ntau[2]<-runif(1, al, bl)
   }
   if (catnum>4){
      a1<-max(max(latent[which(obs==2)]), cut0) 
      b1<-min(min(latent[which(obs==3)]), oldcut[2])
      ntau[1]<-runif(1, a1, b1)
      for (i in 3:(catnum-2) ){
        aaa<-max(max(latent[which(obs==(i))]), ntau[i-2])
        bbb<-min(min(latent[which(obs==(i+1))]), oldcut[i])
        ntau[i-1]<-runif(1, aaa, bbb)
      }
      al<-max(max(latent[which(obs==(catnum-1))]), ntau[catnum-3])
      bl<-min(latent[which(obs==catnum)])
      ntau[catnum-2]<-runif(1, al, bl)
    }
    ntau
 }

  ###############################
  Index<-function(Names){
   uniq <- unique(na.omit(Names))                # get units
   J <- length(uniq)                          # how many units 
   index <- rep (NA, J)                      # how many observations for each unit
   for (i in 1:J) index[Names==uniq[i]] <- i
   index
 }

 ##############################
 OrdinalAugment <- function(obs, maxcat, mincat, cutpoint, mu, sigma){
   if(length(cutpoint)!=(maxcat-2)){
    warning("numbers of categories and cutpoints do not match")
        return(NULL)
    }
   ystar <- c()
   if (obs==mincat){
     ystar<-rtnorm(1, mu, sigma,  upper=0)
   }else{
     if  (obs==(mincat+1)){
       ystar<-rtnorm(1, mu, sigma, lower=0,  upper=cutpoint[(obs-1)])
     }else{
       if (obs==maxcat){
         ystar<-rtnorm(1, mu, sigma, lower=cutpoint[length(cutpoint)])
       }else{
         ystar <- rtnorm(1, mu, sigma, lower=cutpoint[(obs-2)], upper=cutpoint[(obs-1)])
       }
    }
   }
   ystar
  }

 ################################

 outerself<-function(x) outer(x,x)

 ########################################
 PhiTerm0 <- function(Entry,  unit, NN, yystar,newb,
                      PX, PW, PS, Ti, CC, newbeta){
   Ymtotal<-0
   Ymtotal2<-0
   for (n in 1: NN){
     rows<-Entry[which(unit==n)]
     upZ <- yystar[rows]
     upX <- as.matrix(PX[rows,])
     upW <- as.matrix(PW[rows,])
     upS <- as.matrix(PS[rows,])
     TT<-Ti[rows]
     upccin <- as.matrix(CC[TT,])
     JP <- length(rows)
     rant <- c()
     for (d in 1: JP){
       rant[d] <- upS[d,]%*%cbind(upccin[d,])
     }
     marki <- upZ-upX%*%cbind(newbeta)-upW%*%cbind(newb[n,])-rant
         
     if (length(rows)>1){
       sumt2<-crossprod(marki)-marki[length(marki)]^2
       that<-0
       for (i in 2: length(marki)){
         that<-that+marki[i-1]*marki[i]
       }
       Ymtotal <- Ymtotal + sumt2
       Ymtotal2 <-Ymtotal2+ that
      }
    }
    return (c(Ymtotal, Ymtotal2))
  }

  ###################################
  PhiPropARp <- function(Ymt, Ymt2, scale=1){
    HATPhi <- pseudoinverse(Ymt)*scale
    HATphi <- HATPhi %*%Ymt2
    unitcir=2
    while(any(unitcir>1)){
      phi.draw <- mvrnorm(1, mu=HATphi, Sigma=HATPhi)
      ord=length(phi.draw)
      GG<-matrix(NA, ncol=ord, nrow=ord)
      GG[2:ord, 1:(ord-1)] <- diag(ord-1)
      GG[1,] <- phi.draw
      GG[-1, ord] <- 0
      eigGG <- eigen(GG)$values
      unitcir <- abs(eigGG)
    }
    phi.draw
  }

 ##################################

   PhiTermARpT <- function(Entry, unit, NN, yystar,newb, PX,
                           PW, PS, Ti, CC, newbeta, nlag){
    Ymtotal<-matrix(0, nrow=nlag, ncol=nlag)
    Ymtotal2<-0
    for (n in 1: NN){
      rows<-Entry[which(unit==n)]
      upZ <- yystar[rows]
      upX <- as.matrix(PX[rows,])
      upW <- as.matrix(PW[rows,])
      upS <- as.matrix(PS[rows,])
      TT<-Ti[rows]
      upccin <- as.matrix(CC[TT, ])
      JP <- length(rows)
      rant <- c()
      for (d in 1: JP){
       rant[d] <- upS[d,]%*%cbind(upccin[d,])
      }
      marki <- upZ-upX%*%cbind(newbeta)-upW%*%cbind(newb[n,])-rant
      if (length(rows)>nlag){
        for(q in (nlag+1):length(marki)){
          trans <- c()
          for (w in 1:nlag){
            trans[w] <-marki[q-w]
          }
          Ymtotal<-Ymtotal+cbind(trans)%*%rbind(trans)
          Ymtotal2<-Ymtotal2+cbind(trans)*marki[q]
        }
      }
    }
    return (cbind(Ymtotal, Ymtotal2))
  }


  ###################

 rWishart <- function (n, df, S) 
{
    thispackage <- "bayesSurv"
    if (is.null(dim(S))) 
        wdim <- 1
    else {
        wdim <- nrow(S)
        if (ncol(S) != wdim) 
            stop("S must be a square matrix")
    }
    lW <- (wdim * (wdim + 1))/2
    if (df <= wdim - 1) 
        stop(paste("df must be > ", wdim - 1, sep = ""))
    Si <- chol(S)
    Si <- chol2inv(Si)
    Sitri <- Si[lower.tri(Si, diag = TRUE)]
    SAMPLE <- .C("rwishartR3", W = double(n * lW), work = double(2 * 
        wdim * wdim), df = as.double(df), invS = as.double(Sitri), 
        dim = as.integer(wdim), nrandom = as.integer(n), PACKAGE = thispackage)
    Imat <- diag(wdim)
    rowsI <- row(Imat)[lower.tri(row(Imat), diag = TRUE)]
    colsI <- col(Imat)[lower.tri(col(Imat), diag = TRUE)]
    naam <- paste("(", rowsI, ".", colsI, ")", sep = "")
    if (n == 1) {
        if (wdim > 1) {
            names(SAMPLE$W) <- naam
        }
    }
    else {
        if (wdim == 1) {
            names(SAMPLE$W) <- 1:n
        }
        else {
            SAMPLE$W <- matrix(SAMPLE$W, byrow = TRUE, ncol = lW, 
                nrow = n)
            rownames(SAMPLE$W) <- 1:n
            colnames(SAMPLE$W) <- naam
        }
    }
    return(SAMPLE$W)
}
  ################
  ss<-function(ff, Yb, Xb, bb){
   dds<--1/2*(t(Yb-Xb%*%cbind(bb))%*%pseudoinverse(ff)%*%(Yb-Xb%*%cbind(bb))) 
   dds<-as.vector(dds)
  }

  ################
  sumtrans<-function(Mat){
    total<-0
    for (r in 1: nrow(Mat)){
      total<-total+outerself(Mat[r,])
    }
    total
  }

  ################
  Wishart <- function (element, col){
    mm <- array(0, c(col,col))
    mm[lower.tri(mm, diag = TRUE)] <-element
    res <- mm + t(mm)
    diag(res) <- diag(res)/2
    return (res)
  }

  ################
  V.matrix <- function(Omega, col){
    kap <- min(eddcmp(Omega)$evalues)/2
    J <- Omega-diag(kap, col)
    V <- t(chol(J))
    return (list(kap, V))
  }

  ###################################
  ## Function for Data Arrangement ##
  ###################################
  RegressionData <- function(yob, Unit, Time, Fixed, UnitRandom, TimeRandom,
                             UnitPred, TimePred){
    if(!is.numeric(UnitPred)& !is.numeric(TimePred)){
      NewData <- DataArrange.NP(id1=Unit, id2=Time, y=yob, Fixed=Fixed,
                                UnitRandom=UnitRandom, TimeRandom=TimeRandom)
    } else{
        if (!is.numeric(UnitPred) &is.numeric(TimePred)){
           NewData <- DataArrange.TP(id1=Unit, id2=Time, y=yob, Fixed=Fixed,
                                     UnitRandom=UnitRandom, TimeRandom=TimeRandom,
                                     TimePred=TimePred)
        } else{
           if (is.numeric(UnitPred) &!is.numeric(TimePred)){
              NewData <- DataArrange.UP(id1=Unit, id2=Time, y=yob, Fixed=Fixed,
                                        UnitRandom=UnitRandom, TimeRandom=TimeRandom,
                                        UnitPred=UnitPred)
             }else{
               NewData <- DataArrange(id1=Unit, id2=Time, y=yob, Fixed=Fixed,
                                      UnitRandom=UnitRandom, TimeRandom=TimeRandom,
                                      UnitPred=UnitPred, TimePred=TimePred)
             }
        }
   }
    return(NewData)
  }

 ###################################
 DataArrange <- function(id1, id2, y, Fixed, UnitRandom,
                         TimeRandom, UnitPred, TimePred){
  # Make sure the design matrices are matrix formats
   if(!is.matrix(Fixed)){Fixed <- as.matrix(Fixed)}
   if(!is.matrix(UnitRandom)){UnitRandom <- as.matrix(UnitRandom)}
   if(!is.matrix(TimeRandom)){TimeRandom <- as.matrix(TimeRandom)}
   if(!is.matrix(UnitPred)){UnitPred <- as.matrix(UnitPred)}
   if(!is.matrix(TimePred)){TimePred <- as.matrix(TimePred)}

   J<-length(unique(na.omit(id1)))  # how many units, in case NAs
   len <- 1:length(na.omit(id1))    # how many observations in the dataset

   #_____  Create group-level predictors____________#

   # create a design matrix  at the group level of units
   UnitP <-GroupPred(X=as.matrix(UnitPred), index=id1)
   A <- as.matrix(UnitP)      
   # create a design matrix  at the group level of units
   UnitQ <-GroupPred(X=as.matrix(TimePred), index=id2) 
   B <- as.matrix(UnitQ)   

   #_______ Interacting terms in the individual level ____#

   # interacting terms of unit-level variables and individual
   # variables--- they are with fixed effects in the reduced model
   RanFix <- random(Var=UnitRandom, Pred=A, index=id1)

   # interacting terms of time-level variables and individual
   #  variables--- they are with fixed effects in the reduced model
   TimFix <- randomt2(Var=TimeRandom, Pred=B,
                      index1=id1, index2=id2, entry=len)

   #________ Get the dimensions ___________________________#

   RD <- ncol(RanFix)
   TD <- ncol(TimFix)
   FD <- ncol(Fixed)
   UD <- ncol(UnitRandom)
   UT <- ncol(TimeRandom)

   #________ Creat the dataset ______________________________#

   # conbined the data with two indices
   RawData<-cbind(id1, id2, y, Fixed, RanFix, TimFix, UnitRandom, TimeRandom) 
   Data2 <- remove.NA(RawData)

   # get rid of entries with missing data---no drop-out:
   # THIS NEEDS FURTHER THINKING---HOW TO DEAL WITH MISSING DATA

   DD <- ncol(Data2)
   #__________ Recale the data and get right covariates ___________#
   newid1 <- Data2[,1]
   newid2 <- Data2[,2]
   Y<-Data2[,3]
   X <- as.matrix(cbind(Data2[, c(4: (3+FD+TD+RD))]))
   W <- as.matrix(Data2[,  c((4 +FD+TD+RD):(3+FD+TD+RD+UD))])
   S <- as.matrix(Data2[, (4+FD+TD+RD+UD):(3+FD+TD+RD+UD+UT)])

   return(list(newid1, newid2,Y, X, W, S))
  }

 ###################################
 DataArrange.UP <- function(id1, id2, y, Fixed,
                             UnitRandom, TimeRandom, UnitPred){
   # Make sure the design matrices are matrix formats
   if(!is.matrix(Fixed)){Fixed <- as.matrix(Fixed)}
   if(!is.matrix(UnitRandom)){UnitRandom <- as.matrix(UnitRandom)}
   if(!is.matrix(TimeRandom)){TimeRandom <- as.matrix(TimeRandom)}
   if(!is.matrix(UnitPred)){UnitPred <- as.matrix(UnitPred)}

   J<-length(unique(na.omit(id1)))  # how many units, in case NAs
   len <- 1:length(na.omit(id1))    # how many observations in the dataset

   #______ Create group-level predictors______________#

   # create a design matrix  at the group level of units
   UnitP <-GroupPred(X=as.matrix(UnitPred), index=id1)
   A <- as.matrix(UnitP)      
   #_____ Interacting terms in the individual level _____#
   # interacting terms of unit-level variables and individual
   # variables--they are with fixed effects in the reduced model
   RanFix <- random(Var=UnitRandom, Pred=A, index=id1)
   #_______ Get the dimensions _________________________#

   RD <- ncol(RanFix)
   FD <- ncol(Fixed)
   UD <- ncol(UnitRandom)
   UT <- ncol(TimeRandom)

   #__________ Creat the dataset _______________________#
   # conbined the data with two indices
   if (ncol(TimeRandom)==1&length(unique(TimeRandom))==1){
     RawData<-cbind(id1, id2, y, Fixed,  RanFix, UnitRandom)
   }else{
     RawData<-cbind(id1, id2, y, Fixed, RanFix, UnitRandom, TimeRandom)
   }
   Data2 <- remove.NA(RawData)

   # get rid of entries with missing data---no drop-out:
   # THIS NEEDS FURTHER THINKING---HOW TO DEAL WITH MISSING DATA 
   DD <- ncol(Data2)
   #______ Recale the data and get right covariates ___#

   newid1 <- Data2[,1]
   newid2 <- Data2[,2]
   Y<-Data2[,3]
   X <- as.matrix(cbind(Data2[, c(4: (3+FD+RD))]))
   W <- as.matrix(Data2[, c((4+FD+RD):(3+FD+RD+UD))])
   if (ncol(TimeRandom)==1&length(unique(TimeRandom))==1){
     S <- as.matrix(rep(1, nrow(X)))
   }else{
     S <- as.matrix(Data2[, c((4+FD+RD+UD):(3+FD+RD+UD+UT))])
   }
   return(list(newid1, newid2,Y, X, W, S))
 }

 ###################################
 DataArrange.TP <- function(id1, id2, y, Fixed, UnitRandom,
                            TimeRandom,  TimePred){
   # Make sure the design matrices are matrix formats
   if(!is.matrix(Fixed)){Fixed <- as.matrix(Fixed)}
   if(!is.matrix(UnitRandom)){UnitRandom <- as.matrix(UnitRandom)}
   if(!is.matrix(TimeRandom)){TimeRandom <- as.matrix(TimeRandom)}
   if(!is.matrix(TimePred)){TimePred <- as.matrix(TimePred)}

   J<-length(unique(na.omit(id1)))  # how many units, in case NAs
   len <- 1:length(na.omit(id1))    # how many observations in the dataset

   #__________ Create group-level predictors____________#
   # create a design matrix  at the group level of units
   UnitQ <-GroupPred(X=as.matrix(TimePred), index=id2) 
   B <- as.matrix(UnitQ)   

   #______ Interacting terms in the individual level _____#
   # interacting terms of time-level variables and individual
   # variables--- they are with fixed effects in the reduced model
   TimFix <- randomt2(Var=TimeRandom, Pred=B,
                      index1=id1, index2=id2, entry=len)

   #__________ Get the dimensions _________________#

   TD <- ncol(TimFix)
   FD <- ncol(Fixed)
   UD <- ncol(UnitRandom)
   UT <- ncol(TimeRandom)

   #________ Creat the dataset ____________________#

   # conbined the data with two indices
     RawData<-cbind(id1, id2, y, Fixed, UnitRandom, TimeRandom, TimFix)
   Data2 <- remove.NA(RawData)

   # get rid of entries with missing data---no drop-out:
   # THIS NEEDS FURTHER THINKING---HOW TO DEAL WITH MISSING DATA
   DDT <- ncol(Data2)
   #__________ Recale the data and get right covariates _____#
   newid1 <- Data2[,1]
   newid2 <- Data2[,2]
   Y<-Data2[,3]
   X <- as.matrix(cbind(Data2[, -c(c(1:3), c((DDT-UT-UD-TD+1):(DDT-TD)))]))
   if (ncol(UnitRandom)==1&length(unique(UnitRandom))==1){
     W <- as.matrix(rep(1, nrow(X)))
   }else{
     W <-  as.matrix(Data2[, (4+FD):(3+FD+UD)])
   }
   S <- as.matrix(Data2[,  c((DDT-UT+1):DDT)])
   return(list(newid1, newid2,Y, X, W, S))
 }

 # Function: there is no group level predictors, if only varing-intercepts
 # are included, the universal intercept should not be included.
 DataArrange.NP <- function(id1, id2, y, Fixed, UnitRandom, TimeRandom){
   if(!is.matrix(Fixed)){Fixed <- as.matrix(Fixed)}
   if(!is.matrix(UnitRandom)){UnitRandom <- as.matrix(UnitRandom)}
   if(!is.matrix(TimeRandom)){TimeRandom <- as.matrix(TimeRandom)}

   # Make sure the design matrices are matrix formats
   J<-length(unique(na.omit(id1)))  # how many units, in case NAs
   len <- 1:length(na.omit(id1))    # how many obsedrvations in the dataset

   FD <- ncol(Fixed)
   UD <- ncol(UnitRandom)
   UT <- ncol(TimeRandom)
   
   #_________ Creat the dataset ____________#
   # conbined the data with two indices
   if ((UD==1&length(unique(UnitRandom))==1)&(UT==1&length(unique(TimeRandom))==1)){
     RawData<-cbind(id1, id2, y, Fixed)
     Data <- remove.NA(RawData)
     newid1 <- Data[,1]
     newid2 <- Data[,2]
     Y<-Data[,3]
     X <- as.matrix(Data[, -c(1:3)])
     W <- as.matrix(rep(1,nrow(X)))
     S <- as.matrix(rep(1,nrow(X)))
   }else{
     if ((UD==1&length(unique(UnitRandom))==1)&(UT!=1|length(unique(TimeRandom))!=1)){
       RawData<-cbind(id1, id2, y, Fixed,  TimeRandom)
       Data <- remove.NA(RawData)
       RR <- ncol(Data)
       newid1 <- Data[,1]
       newid2 <- Data[,2]
       Y<-Data[,3]
       X <- as.matrix(Data[, -c(1:3)])
       W <- as.matrix(rep(1,nrow(X)))
       S <- as.matrix(Data[, (RR-UT+1):RR])
      }else{
        if ((UD!=1|length(unique(UnitRandom))!=1)&(UT==1&length(unique(TimeRandom))==1)){
          RawData<-cbind(id1, id2, y, Fixed, UnitRandom)
          Data <- remove.NA(RawData)
          RR <- ncol(Data)
          newid1 <- Data[,1]
          newid2 <- Data[,2]
          Y <- Data[,3]
          X <- as.matrix(Data[, -c(1:3)])
          W <- as.matrix(Data[, (RR-UD+1):RR])
          S <- as.matrix(rep(1,nrow(X)))
         }else{
          RawData<-cbind(id1, id2, y, Fixed, UnitRandom, TimeRandom)
          Data <- remove.NA(RawData)
          RR <- ncol(Data)
          newid1 <- Data[,1]
          newid2 <- Data[,2]
          Y<-Data[,3]
          X <- as.matrix(Data[, -c(1:3)])
          W <- as.matrix(Data[, (RR-UD -UT+1):(RR-UT)])
          S <- as.matrix(Data[, (RR-UT+1):RR])
        }
      }
   }
   return(list(newid1, newid2,Y, X, W, S))
 }

 ###################################
fast.svd <- function (m, tol) {
    n = dim(m)[1]
    p = dim(m)[2]
    EDGE.RATIO = 2
    if (n > EDGE.RATIO * p) {
        return(psmall.svd(m, tol))
    }
    else if (EDGE.RATIO * n < p) {
        return(nsmall.svd(m, tol))
    }
    else {
        return(positive.svd(m, tol))
    }
}

 ###################################
 GroupPred <- function(X, index){
   b<-ncol(X)
   J <- length(unique(na.omit(index)))
   XPred<-matrix(0, nrow=J, ncol=b)
   for (i in 1:b){
     Y<-X[,i]
     a <- length(Y)
     Z <- is.na(X[,i])==TRUE
     m <- which(Z==TRUE)
     s <- length(Z[which(Z==FALSE)])
     if (s==a){
       XP <- tapply(Y, index, mean)
       XPred[,i]<-XP
     }else{
       XP <- tapply(X[-m, i], index[-m], mean)
       XPred[,i]<-XP
    }
    XGroup<-XPred
   }
 }

  ###################################
  GroupPred2 <- function(X, index){
    b<-ncol(X)
    J <- length(unique(na.omit(index)))
    XPred<-matrix(0, nrow=J, ncol=b)
    for (i in 1:b){
      if(length(which(is.na(X[,i])))>0){
        average <- tapply(na.omit(X[,i]), index[-which(is.na(X[,i]))], mean)
        XPred[, i] <- average
      }else{
        average <- tapply(X[,i], index, mean)
        XPred[, i] <- average
      }
    }
  }

  ###################################
make.positive.definite <- function (m, tol) {
    if (!is.matrix(m)) 
        m = as.matrix(m)
    d = dim(m)[1]
    if (dim(m)[2] != d) 
        stop("Input matrix is not square!")
    es = eigen(m)
    esv = es$values
    if (missing(tol)) 
        tol = d * max(abs(esv)) * .Machine$double.eps
    delta = 2 * tol
    tau = pmax(0, delta - esv)
    dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
    return(m + dm)
}

  ###################################
nsmall.svd <- function (m, tol) 
{
    B = m %*% t(m)
    s = svd(B, nv = 0)
    if (missing(tol)) 
        tol = dim(B)[1] * max(s$d) * .Machine$double.eps
    Positive = s$d > tol
    d = sqrt(s$d[Positive])
    u = s$u[, Positive, drop = FALSE]
    v = crossprod(m, u) %*% diag(1/d, nrow = length(d))
    return(list(d = d, u = u, v = v))
}

  ###################################
  random <- function(Var, Pred, index){
    J<-length(unique(na.omit(index)))
    RanFix <- matrix(0, nrow=1, ncol=ncol(Var)*ncol(Pred))
    for (j in 1:J){
      W <- as.matrix(Var[which(index==j),])
      a <- kronecker(Diagonal(ncol(Var)), t(as.matrix(Pred[j,])))
      RanFix<-rbind(RanFix, as.matrix(W%*%a))
    }
    TTT<-RanFix[-1,]
  }

  ###################################
  Random.Coef<- function(Cova, random, coef.pred, RBB){
    random <- as.matrix(random)
    coef.pred <- as.matrix(coef.pred)
    N.unit <- nrow(Cova)
    N.coef <- RBB
    beta.random <-matrix(NA, ncol=N.coef, nrow=N.unit)
    for (n in 1: N.unit){
      beta1 <- c()
      for (j in 1:N.coef){
        beta1[j] <- Cova[n,]%*%cbind(coef.pred[j,])+random[n,j]
      }
      beta.random[n,] <- beta1
    }
    return(as.vector(t(beta.random)))
 }

##################################
psmall.svd <- function (m, tol) 
{
    B = crossprod(m)
    s = svd(B, nu = 0)
    if (missing(tol)) 
        tol = dim(B)[1] * max(s$d) * .Machine$double.eps
    Positive = s$d > tol
    d = sqrt(s$d[Positive])
    v = s$v[, Positive, drop = FALSE]
    u = m %*% v %*% diag(1/d, nrow = length(d))
    return(list(d = d, u = u, v = v))
}
 ################################
positive.svd <- function (m, tol) 
{
    s = svd(m)
    if (missing(tol)) 
        tol = max(dim(m)) * max(s$d) * .Machine$double.eps
    Positive = s$d > tol
    return(list(d = s$d[Positive], u = s$u[, Positive, drop = FALSE], 
        v = s$v[, Positive, drop = FALSE]))
}
################################
 randomt2 <- function(Var, Pred, index1, index2, entry){
   index <- na.omit(index1)
   index2 <- na.omit(index2)
   index2 <- index2-min(index2)+1
   Var <- Var
   Pred <- na.omit(Pred)
   J1<-length(unique(na.omit(index)))
   WWW <- matrix(NA, nrow=length(entry), ncol=ncol(Var)*ncol(Pred))
   for (j in index){
     rows <- entry[which(index==j)]
     TT <- index2[rows]
     W <- as.matrix(Var[rows,])
     PW <- as.matrix(Pred[TT,])
     Ee <- length(rows)
     Tt <- matrix(NA, ncol=ncol(Var)*ncol(Pred), nrow=Ee)
     for (s in 1:Ee){
       Tt[s,] <- as.vector(t(cbind(W[s,])%*%rbind(PW[s,])))
     }
     WWW[rows,] <- Tt    
   }
   WWW
 }

############################
remove.NA <- function (data) 
{
    n1 <- dim(data)[1]
    n2 <- dim(data)[2]
    newdata <- as.matrix(data[, 3:n2])
    n2 <- n2 - 2
    index <- rep(NA, n1)
    for (i in 1:n1) {
        if (sum(is.finite(newdata[i, ])) == n2) 
            index[i] <- i
    }
    index <- index[is.finite(index) == TRUE]
    data[index, ]
}

############################
dmnorm <- function (x, mean = rep(0, d), varcov, log = FALSE) 
{
    d <- if (is.matrix(varcov)) 
        ncol(varcov)
    else 1
    if (d > 1 & is.vector(x)) 
        x <- matrix(x, 1, d)
    n <- if (d == 1) 
        length(x)
    else nrow(x)
    X <- t(matrix(x, nrow = n, ncol = d)) - mean
    Q <- apply((solve(varcov) %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
    logPDF <- as.vector(Q + d * logb(2 * pi) + logDet)/(-2)
    if (log) 
        logPDF
    else exp(logPDF)
}

################
pseudoinverse <- function (m, tol) 
{
    msvd = fast.svd(m, tol)
    if (length(msvd$d) == 0) {
        return(array(0, dim(m)[2:1]))
    }
    else {
        return(msvd$v %*% (1/msvd$d * t(msvd$u)))
    }
}

################################
 VaryingCoef <- function(KK, Betai, Residual, m, RB,
                         GU, unitint.add){
   betaii <- matrix(NA, nrow=1, ncol=RB*GU)
   for (i in 1: m){
     # if there is a random intercept
     if (unitint.add==1){
        Res1 <- matrix(Residual[i,], ncol=(RB+1))
        Res <- as.matrix(Res1[, -1])
     }else{
        Res <- matrix(Residual[i,], ncol=(RB))
     }
     Betaii <- matrix(Betai[i,], nrow=RB)
     beta2 <- Random.Coef(Cova=KK, random=Res, coef.pred=Betaii, RBB=RB)
     betaii <- rbind(betaii, as.vector(beta2))
   }
   betai <- betaii[-1,]
   return(betai)
 }
