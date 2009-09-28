#******************************************************************************
# Marginal Likelihood  ---  Program to compute Marginal Likelihood by using the
#                            methods of Chib (1995) and Chib (2001)
#                           This program is only for models with two clusterings;
#                           models with only one clustering should be navegated
#                           to another separate file
#                                                                             
# Warning: in this file, variance (covariance) for {b_i} and {c_t} are matrices
#                                                                              
# Author: Xun Pang                                                             
#                                                                              
# Usage: This program contains 5 reduced runs for binary model and more for    
#         ordinal model                                                        
#         Thin the MCMC output, and use 1,000 posteriors draws to estimate     
#         since for storage reason, the data augmentation and auxiliary para-  
#         meter are stored every 100 iteration, the output has to been handled 
#         carefully                                                            
#                                                                              
# Note: Now this program is only for binary outcomes, but for ordinal cases    
#        the cutpoint part can be added as the last ordinate and multiple      
#        reduced runs are needed, depending on how many cutpoints we have      
#
# Date of Created: 12/31/2008
# Date of Modified: 1/3/2009
#******************************************************************************
# call functions in other files

 MargLike.Binary.Builtin <- function(Y,X, W, S, id1, id2, mcmc.output, reduced.mcmc, reduced.burnin, 
                              nlag, beta0,B0, D0, d0,E0, e0, tracking){

   m <- reduced.mcmc
   N <- length(Y)
   # Dimension of all fixed parameter vector beta 
   FB <- ncol(X)
   # Dimension of random fixed parameter vector bi
   RB <- ncol(W)
   # Dimension of random fixed parameter vector ci
   TB <- ncol(S)
 
  # create indices for two non-nested groupings ##########
   entry <- seq(1: N)                   # the global length of observations
   GU <- max(id1)                       # how many units
   GT <- max(id2)                       # how many time periods in total
   rho.pos <- mcmc.output[[1]]
   bi.pos <- mcmc.output[[3]]
   ct.pos <- mcmc.output[[4]]
   ystar.pos <- mcmc.output[[8]]
   beta.pos <- mcmc.output[[2]]     
   if (nlag==1){
     rho.pos.mean <- mean(mcmc.output[[1]])
   }
   if (nlag>1){
    rho.pos.mean <- apply(mcmc.output[[1]], 2, mean)       # mean of posterior distribution of rho
   }
   beta.pos.mean <- apply(mcmc.output[[2]], 2, mean)
   bi.pos.mean <- apply(mcmc.output[[3]], 2, mean)
   ct.pos.mean <- apply(mcmc.output[[4]], 2, mean)
   u.pos.mean <- apply(mcmc.output[[7]], 2, mean)
   ystar.pos.mean <- apply(mcmc.output[[8]], 2, mean)
   if (ncol(W)>1){
    D.pos.mean<- apply(mcmc.output[[6]], 2, mean)
    D.matrix <- Wishart(D.pos.mean, RB)
    }else{
    D.pos.mean <- mean(mcmc.output[[6]])
    D.matrix <- as.matrix(D.pos.mean)
  }
  if (ncol(S)>1){
    E.pos.mean <- apply(mcmc.output[[5]], 2, mean)
    E.matrix <- Wishart(E.pos.mean, TB)
  }else{
    E.pos.mean <- mean(mcmc.output[[5]])
    E.matrix <- as.matrix(E.pos.mean)
  }

   bi.matrix <- matrix(bi.pos.mean, ncol=RB)
   ct.matrix <- matrix(ct.pos.mean, ncol=TB)
   mcmc.pos <- list(rho=rho.pos.mean, beta=beta.pos.mean, bi= bi.pos.mean, ct= ct.pos.mean,
                    u=u.pos.mean, ystar=ystar.pos.mean, D=D.pos.mean, E=E.pos.mean)

   ######################### Stage Two: Reduced Runs ######################

  # RR1(reduced run I): fixed rho. The output is named RR1.OUT
  monitored.parameter <- c( "beta", "bi", "ct", "D", "E", "ystar", "u", "q.rho")
  fixed.parameter <- "rho"
  RR1.OUT <- Reduced.Binary(Y=Y,X=X, W=W, S=S,unit.id=id1, time.id=id2, m=reduced.mcmc,
                            burnin=reduced.burnin, mcmc.pos=mcmc.pos, 
                            beta0=beta0,B0=B0, D0=D0 , d0=d0, E0=E0, e0=e0, nlag=nlag,
                            tracking=tracking,  monitor=monitored.parameter, fixed=fixed.parameter)                            
  # RR2: fixed rho, E, and {c_t}. The output is named RR2.OUT                              
  monitored.parameter <- c( "beta",  "ystar", "D")
  fixed.parameter <- c( "E", "ct", "rho")
  RR2.OUT <- Reduced.Binary(Y=Y,X=X, W=W, S=S,unit.id=id1, time.id=id2, m=reduced.mcmc,
                            burnin=reduced.burnin, mcmc.pos=mcmc.pos, 
                            beta0=beta0,B0=B0, D0=D0 , d0=d0, E0=E0, e0=e0, nlag=nlag,
                            tracking=tracking,  monitor=monitored.parameter, fixed=fixed.parameter)                           
  # RR3: fixed rho, E, {c_t}, {b_i}, and D. The output is named RR3.OUT                             
  monitored.parameter <- c( "ystar")
  fixed.parameter <- c( "E","bi","D", "ct", "rho")
  RR3.OUT <- Reduced.Binary(Y=Y,X=X, W=W, S=S,unit.id=id1, time.id=id2, m=reduced.mcmc,
                            burnin=reduced.burnin, mcmc.pos=mcmc.pos, 
                            beta0=beta0,B0=B0, D0=D0 , d0=d0, E0=E0, e0=e0, nlag=nlag,
                            tracking=tracking,  monitor=monitored.parameter, fixed=fixed.parameter)
     
  # RR4: fixed rho, E, {c_t}, {b_i}, D, and beta. The output is named RR4.OUT
  monitored.parameter <- "u"
  fixed.parameter <- c( "E","bi","D", "ct", "rho")
  RR4.OUT <-Reduced.Binary(Y=Y,X=X, W=W, S=S,unit.id=id1, time.id=id2, m=reduced.mcmc,
                            burnin=reduced.burnin, mcmc.pos=mcmc.pos, 
                            beta0=beta0,B0=B0, D0=D0 , d0=d0, E0=E0, e0=e0, nlag=nlag,
                            tracking=tracking,  monitor=monitored.parameter, fixed=fixed.parameter)  
  #_____________________________ Stage Three ___________________________________#
  # Compute each ordinate in the logrithm scale
  #@@@@@@@@@@@@@@@@@@@
  # The prior ordinate
  # beta
  beta.prior <- dmvnorm(beta.pos.mean, mean=beta0, sigma=B0, log=TRUE)
                              
  # bi and D
  bi.prior <- 0
  for (p in 1: max(id1)){
    bi.prior <- bi.prior+ dmvnorm(bi.matrix[p,], mean=rep(0, RB),
                                  sigma=solve(D.matrix), log=TRUE)
  }
  D.prior <- log(dwish(W=D.matrix, v=d0, S=solve(D0))) # check here: inverse or not

  # ct and E
  ct.prior <- 0
  for (q in 1:max(id2)){
   ct.prior <-  ct.prior+dmvnorm(ct.matrix[q,], mean=rep(0, TB),
                                   sigma=solve(E.matrix), log=TRUE)
 }
 if (ncol(S)==1){                     # if S only contains 1 element (most likely a constant)
   E.prior<-log(dgamma(E.matrix, e0/2, rate=E0/2))
 }else{
   E.prior <- log(dwish(W=E.matrix, v=e0, S=solve(E0))) # check here: inverse or not
 }                     
 
  #@@@@@@@@@@@@@@@@@@@
  # The likelihood ordinate (integrate u out)
     u.likelihood <- RR4.OUT[[8]]
     likelihood.pre <- rep(NA, reduced.mcmc)
     for (i in 1:reduced.mcmc){
       u.poss <- u.likelihood[i,]                           
       if (nlag==1){
         var.phig <- 1/(1- rho.pos.mean^2)                          
         TCov <- Cov.AR1(col=GT, rho= rho.pos.mean,sigma=var.phig)
       }else{
         CCov <- Cov.ARp(phig=rho.pos.mean, col=GT)        
         TCov <- CCov[[2]]
         Sig <- TCov[1:nlag, 1:nlag]
         SSig<-solve(Sig) 
         var <- Sig[1,1]
      }
      likelihood <- c()
      for (n in 1: max(id1)){                              # loop for each unit
        rows <- entry[which(id1==n)]                # the obserbations of unit n
        nT <- length(rows)
        yi <- Y[rows]
        Xi <- as.matrix(X[rows,])
        Wi <- as.matrix(W[rows,])
        Si <- as.matrix(S[rows,])
        bbbn <- bi.matrix[n,]
        TT <- id2[rows]
        Tccin <- as.matrix(ct.matrix[TT,])
        # u.term <- u.pos.mean[rows]
        u.term <- u.poss[rows]
      
        rant <- c()
        for (l in 1:nT){
         rant[l] <- Si[l,]%*%cbind(Tccin[l,])
        }
        nCov <- TCov[1:nT, 1:nT]
        nVK <- V.matrix(Omega=nCov, col=nT)
        nV <- nVK[[2]]                          
        nK <- nVK[[1]]
        u.decom <- nV%*%cbind(u.term)
        like <- rep(NA, nT)
        mean.like <-  (Xi%*% cbind(beta.pos.mean)+Wi%*%cbind(bbbn)+rant+u.decom)/sqrt(nK)
        for (j in 1:nT){
          individual <- pnorm(mean.like[j])^(yi[j])*(1-pnorm(mean.like[j]))^(1-yi[j])
          if (individual==0){
            individual <- 0.00001
          }
        like[j] <- individual
       }
        likelihood[rows] <- like
      }
      likelihood.pre[i] <- prod(likelihood)
     }
       likelihood.term <- log(mean(likelihood.pre))
      cat("likelihood", likelihood.term, "\n")
  #@@@@@@@@@@@@@@@@@@@
  # The posterior ordinate
  #++++++++++++++++++++++
  # rho|y
  # Get the results from the reduced run 1
   beta.rho <- RR1.OUT[[3]]
   b.rho <- RR1.OUT[[4]]
   c.rho <- RR1.OUT[[5]]
   ystar.rho <- RR1.OUT[[9]]
   q.rho <- RR1.OUT[[2]]
   # Compute the numoninator
   upstairs <- rep(NA, reduced.mcmc)
   downstairs <- rep(NA, reduced.mcmc)
   for (i in 1: reduced.mcmc){
     # by using the mcmc output
     bi.porho <- matrix(bi.pos[i,], ncol=RB)
     ct.porho <- matrix(ct.pos[i,], ncol=TB)
     # if nlag=1, call function "AcceptRate0" and "PhiTerm0"
     if (nlag==1){
       alpha <- AcceptRate0(rho=rho.pos.mean, orho=rho.pos[i], unit=id1, precision=1,
                          AY=ystar.pos[i,], AX=X, AS=S, AW=W, Abeta=beta.pos[i,],
                          Ab=bi.porho,AT=id2,  Ac=ct.porho, Entry=entry)
     PhiPr <- PhiTerm0(Entry=entry, unit=id1, NN=max(id1), yystar=ystar.pos[i,],newb=bi.porho,
                       PX=X, PW=W, PS=S, Ti=id2, CC=ct.porho, newbeta=beta.pos[i,])
     PhiSig <- PhiPr[1]
     PhiMean <- PhiPr[2] 
     hatPhi<-1/PhiSig
     hatphi<-hatPhi*PhiMean
     q.term <- dnorm(rho.pos.mean, hatphi, sqrt(hatPhi))
     upstairs[i] <- alpha*q.term
    }else{
     CCov <- Cov.ARp(phig=rho.pos.mean, col=(nlag+1))        
     TCov <- CCov[[2]]
     Sig <- TCov[1:nlag, 1:nlag]
     alpha<-AcceptRateARpT(nrho=rho.pos.mean, orho=rho.pos[i,], Osig=Sig , unit=id1,  AY=ystar.pos[i,],
                           AX=X, AS=S, AW=W, Abeta=beta.pos[i,], Ab=bi.porho, Ac=ct.porho,
                           Entry=entry,  AT=id2,  f=i, nlag=nlag)
     
     PhiPr <- PhiTermARpT(Entry=entry, unit=id1, NN=max(id1), yystar=ystar.pos[i,], newb=bi.porho,
                          PX=X, PW=W, PS=S, Ti=id2, CC=ct.porho, newbeta=beta.pos[i,], nlag=nlag)
     
     PhiSig <- PhiPr[,-(nlag+1)]
     PhiMean <- PhiPr[,(nlag+1)]
     HATPhi <- solve(PhiSig)
     HATphi <- HATPhi %*%PhiMean
     HATPhi <- make.positive.definite(HATPhi)
     q.term <- dmnorm(rho.pos.mean, mean=HATphi, HATPhi)
     upstairs[i]<- alpha*q.term
    }

    # compute the denominator
    bi.rhoi <- matrix(b.rho[i,], ncol=RB)
    ct.rhoi <- matrix(c.rho[i,], ncol=TB)
    if (nlag==1){
       alpha.q <- AcceptRate0(orho=rho.pos.mean, rho=q.rho[i], unit=id1, precision=1, AY=ystar.rho[i,],
                              AX=X, AS=S, AW=W, Abeta=beta.rho[i,], Ab=bi.rhoi,AT=id2, Ac=ct.rhoi, Entry=entry)
       downstairs[i] <- alpha.q
     } else {
       alpha.q <- AcceptRateARpT(nrho=q.rho[i,], orho=rho.pos.mean, Osig=Sig , unit=id1,
                                 AY=ystar.rho[i,], AX=X, AS=S, AW=W, Abeta=beta.rho[i,], Ab=bi.rhoi,
                                 Ac=ct.rhoi, Entry=entry,  AT=id2,  f=i, nlag=nlag)
       
       downstairs[i] <- alpha.q
     }
    } 
    rho.term.pos <- log(mean(upstairs)) - log(mean(downstairs))
   #++++++++++++++++++++++
   # {c_t}|rho^*, y: directly use reduced No.1
   u.rho <- RR1.OUT[[8]]
   E.rho <- as.matrix(RR1.OUT[[6]])
   # first get the term of C'u
   if (nlag==1){
     var.phig <- 1/(1- rho.pos.mean^2)                       
     TCov <- Cov.AR1(col=GT, rho= rho.pos.mean,sigma=var.phig)
   }else{
     CCov <- Cov.ARp(phig=rho.pos.mean, col=GT)        
     TCov <- CCov[[2]]
   }
       
    ct.term.pos.sum <- 0
    cttt <- c()
    
    for (i in 1: m){
     # get its convariance matrix
      Ei <- Wishart(E.rho[i,], TB)
      u.rhoi <- u.rho[i,]
      kap <- rep(NA, GU)
      qu <- rep(NA,N)
      
      for (n in 1: max(id1)){                              # loop for each unit
        rows <- entry[which(id1==n)]                # the obserbations of unit n
        nT <- length(rows)
        u.draw <- u.rhoi[rows]
        nCov <- TCov[1:nT, 1:nT]
        nVK <- V.matrix(Omega=nCov, col=nT)
        nV <- nVK[[2]]                          
        nK <- nVK[[1]]
        u.decom <- nV%*%cbind(u.draw)
        qu[rows] <- u.decom
        kap[n] <- nK
      }
      bi.ct <- matrix(b.rho[i,], ncol=RB)
      ystart <- ystar.rho[i,]
      ct.term.posp <- 1
      for (j in 1: max(id2)){
         row1<-entry[which(id2==j)]  # Those entries with t=j
         H <- length(row1)           #  #observations at t period
         iden1<-id1[row1]            # The units at t period
         tttn <- as.matrix(bi.ct[iden1,])         # The random unit-specific effects of
         yt1<- ystart[row1]
         Xt1<-X[row1,]
         Wt1<-as.matrix(W[row1,])
         St1 <- as.matrix(S[row1,])
         Ut1 <- qu[row1]
         kapt1 <- kap[iden1]
            
         term7 <- c()
         for (q in 1:H){
           term7[q] <- Wt1[q,]%*% cbind(tttn[q,])
         }
            
         PM <- matrix(0, nrow=H, ncol=H)
         for (c in 1:H){
           PM[c,c] <- kapt1[c]
         }
            
         cov.S <- solve(Ei+t(St1)%*%solve(PM)%*%St1)
         mean.S <- cov.S%*% t(St1) %*%solve(PM)%*%(yt1-term7-Xt1%*%cbind(beta.rho[i,])-Ut1)
         ct.term.posp <- ct.term.posp*dmvnorm(ct.matrix[j,], mean=mean.S, sigma=cov.S)
       }
       cttt[i] <- ct.term.posp
       ct.term.pos.sum <- ct.term.pos.sum+ct.term.posp
      }
      ct.term.pos <- log(ct.term.pos.sum/m)

  #++++++++++++++++++++++
  # E|{c_t^*}, rho^*, y: no reduced run needed
    
   if (ncol(S)==1){                     # if S only contains 1 element (most likely a constant)
     E.term.pos <- log(dgamma(1,(e0+GT)/2, rate=(E0+sum(ct.matrix^2))/2))
   }else{
     OEin <- solve(E0)+sumtrans(Mat=ct.matrix)
     E.Mean.Matrix <- Wishart(E.pos.mean, TB)
     E.term.pos <- log(dwish (W=E.Mean.Matrix, v=e0+GT, S=solve(OEin)))  # check whether there should be "solve" or not solve
     }
     E.term.pos

  #++++++++++++++++++++++
  # {b_i}|E^*, {c_t^*}, rho^*, y: by using reduced Run No. 2
  beta.bi <- RR2.OUT[[3]]
  D.bi <- as.matrix(RR2.OUT[[7]])
  ystar.bi <- RR2.OUT[[9]]
                            
  if (nlag==1){
    var.phig <- 1/(1- rho.pos.mean^2)                       
    TCov <- Cov.AR1(col=GT, rho= rho.pos.mean,sigma=var.phig)
  }else{
    CCov <- Cov.ARp(phig=rho.pos.mean, col=GT)        
    TCov <- CCov[[2]]
  }

  bi.term.post <- 0
  for (i in 1:m){
     Di <- Wishart(D.bi[i,], RB)
     ystarit <- ystar.bi[i,]
     betait <- beta.bi[i,]
     bi.term <- 1
     for (n in 1: max(id1)){                                # loop for each unit
       rows <- entry[which(id1==n)]                    # the obserbations of unit n
       nT <- length(rows)
       nCov <- TCov[1:nT, 1:nT]
       ystari <- ystarit[rows]
       Xi <- X[rows,]
       Wi <- W[rows,]
       Si <- as.matrix(S[rows,])
       b.this <- bi.matrix[n,]
       TT <- id2[rows]
       Tccin <- as.matrix(ct.matrix[TT,])
    
       rant <- c()
       for (l in 1:nT){
         rant[l] <- Si[l,]%*%cbind(Tccin[l,])
       }
     
       Cov.b <- pseudoinverse(Di + t(Wi)%*%solve(nCov)%*%Wi)
       mean.b <- Cov.b%*% t(Wi)%*%solve(nCov)%*%(ystari-Xi%*% cbind(betait)-rant)
       Cov.b <- make.positive.definite(Cov.b)
       bi.term <-bi.term* dmvnorm(b.this, mean=t(mean.b), sigma=Cov.b)
     }
     bi.term.post <- bi.term.post+bi.term
   }
   bi.term.pos <- log(bi.term.post/m)

  #++++++++++++++++++++++
  # D|{b_i}^*, E^*, {c_t^*}, rho^*, y
                            
  ODin <- solve(D0)+sumtrans(Mat=bi.matrix)
  D.Mean.Matrix <- Wishart(D.pos.mean, RB)
  D.term.pos <- log(dwish (W=D.Mean.Matrix, v=e0+GU, S=solve(ODin)))

  #++++++++++++++++++++++
  # beta|D^*, {b_i}^*, E^*, {c_t^*}, rho^*, y

  ystar.beta <- RR3.OUT[[9]]
   if (nlag==1){
     var.phig <- 1/(1- rho.pos.mean^2)                       
     TCov <- Cov.AR1(col=GT, rho= rho.pos.mean,sigma=var.phig)
   }else{
     CCov <- Cov.ARp(phig=rho.pos.mean, col=GT)        
     TCov <- CCov[[2]]
  }
  beta.term.post <- 0
  for (i in 1:m){                        
    betaterm1 <- matrix(0, nrow=FB, ncol=FB)  
    betaterm2 <- rep(0, FB)                   
    ystarbeta <- ystar.beta[i,]
    for (n in 1: GU){                              # loop for each unit
      rows <- entry[which(id1==n)]                # the obserbations of unit n
      nT <- length(rows)
      yi <- Y[rows]
      ystari <- ystarbeta[rows]
      Xi <- X[rows,]
      Wi <- W[rows,]
      Si <- as.matrix(S[rows,])
      TT <- id2[rows]
      Tccin <- as.matrix(ct.matrix[TT,])
     
      rant <- c()
      for (l in 1:nT){
        rant[l] <- Si[l,]%*%cbind(Tccin[l,])
      }
      
      nCov <- TCov[1:nT, 1:nT]
      var.beta <- pseudoinverse(nCov+Wi%*%solve(D.matrix)%*%t(Wi))
      betaterm1 <- betaterm1+t(Xi) %*% var.beta %*%Xi
      betaterm2 <- betaterm2+t(Xi)%*% var.beta %*% (ystari-rant)
    }
    B1<-pseudoinverse(betaterm1+solve(B0))
    beta1<-t(B1%*%(betaterm2+solve(B0)%*%beta0))
    B1 <- make.positive.definite(B1)
    beta.term <-dmvnorm(beta.pos.mean, beta1, B1)
    beta.term.post <- beta.term.post+beta.term
  }                     
  beta.term.pos <- log(beta.term.post/i)
     
   # Sum up the terms in logrithm scale: output Pos.Ord
   logMargLik <- likelihood.term+beta.prior+bi.prior+D.prior+ct.prior+E.prior-rho.term.pos-
                  beta.term.pos-bi.term.pos-D.term.pos-ct.term.pos-E.term.pos
   return(logMargLik)
 }
