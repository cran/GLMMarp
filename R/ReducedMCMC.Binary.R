#******************************************************************************
# Marginal Likelihood  ---  Program to compute Marginal Likelihood by using the
#                            methods of Chib (1995) and Chib (2001)
#                           This program is only for models with two clusterings;
#                           models with only one clustering should be navegated
#                           to another separate file
#                                                                             
# Note: 
#                                                                              
# Author: Xun Pang                                                             
#                                                                                                                 #                                                                              
# Note:  In this file, variance (covariance) for {b_i} and {c_t} are matrices
#        Now this program is only for binary outcomes, but for ordinal cases    
#        the cutpoint part can be added as the last ordinate and multiple      
#        reduced runs are needed, depending on how many cutpoints we have      
#
# Date of Created: 12/31/2008
# Date of Modified: 24/09/2009
#******************************************************************************

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Note: Fit in the function according to the reduced form
#
# Usage:        Reduced.Binary(Y,X, W, S,  ...)
# Arguments:
#         Y     response variable, dichotomous (0/1)
#         X     matrix of covariates with fixed effects
#         W     matrix of covariates with subject-varying effects
#         S     matrix of covariates with time-varying effects
#   unit.id     index of subjects
#   time.id     index of time periods
#         m     number of iterations to be returned
#    burnin     number of initial iterations which are discarded
#  mcmc.pos     mcmc output from GLMMARp.Binary()
#  monitor      a string contains the names of parameters whose MCMC output are
#               will be returned. The string has to be a subset of
#               ("rho", "beta", "bi", "ct", "D", "E", "ystar", "u") which is the default
#     fixed     a string specifying which parameters are fixed in the reduced run
#    priors     Priors have to be those which are used in the MCMC simulation
#               beta0--- mean vector of the multivariate normal distribution of beta,
#                        all the fixed effects (beta in the reduced form of the model)
#               B0--- variance matrix of the multivariate normal distribution of beta,
#                        all the fixed effects (beta in the reduced form of the model)
#               d0, D0 --- degree of freedom and scale matrix of the Wishart prior
#                          on b_i, the subject-level residual
#               e0, E0 --- degree of freedom and scale matrix of the Wishart prior
#                          on c_t, the time-level residual
#     nlag      the lag order p, integer >= 1
# tracking      every "tacking" iterations, return the information about how many
#               iterations in total have been done   
######################################################################################
Reduced.Binary <- function(Y,X, W, S, unit.id, time.id, m=1000, burnin=500, mcmc.pos,
                            beta0,B0, D0, d0, E0, e0, nlag,tracking,
                            monitor=c("rho","q.rho", "beta", "bi", "ct", "D", "E", "ystar", "u"), fixed="NULL"){

#****************** check whetehr illegal values in monitor and fixed ***************#
   if (any(monitor%in%fixed)){
       warning("parameters in monitor should not be in fixed at the same time")
        return(NULL)
   }

   N <- length(Y)
   # Dimension of all fixed parameter vector beta 
   FB <- ncol(X)
   PFB <- length(beta0)
    if (FB!=PFB){
       cat("Length of the prior on fixed effects is incorrect.\n")
       cat("Length of fixed effect should be ", FB, "\n", sep=" ")
       stop("Respecify beta0 and B0, call GLMMARp.Binary() again. \n") 
    }
   # Dimension of random fixed parameter vector bi
   RB <- ncol(W)
   if(is.matrix(D0)){ 
     PRB <- ncol(D0)
   }else{
     PRB <- length(D0)
   }
   if (RB!=PRB){
       cat("Dimension of the prior on unit random effects is incorrect.\n")
       cat("Length of D0 should be ", RB, "*", RB,  "\n", sep=" ")
       stop("Respecify D0, call GLMMARp.Binary() again. \n") 
    }
   # Dimension of random fixed parameter vector ci
   TB <- ncol(S)
   if(is.matrix(E0)){
      PTB <- ncol(E0)
    }else{
      PTB <- length(E0)
    }
   if (TB!=PTB){
       cat("Dimension of the prior on time random effects is incorrect.\n")
       cat("Length of E0 should be ", TB, "*", TB,  "\n", sep=" ")
       stop("Respecify E0, call GLMMARp.Binary() again. \n") 
    }
 
  # create indices for two non-nested groupings ##########
   entry <- seq(1: N)                   # the global length of observations
   GU <- max(unit.id)                       # how many units
   GT <- max(time.id)                       # how many time periods in total

   # initial values
   init.phi <- mcmc.pos$rho      # mean of posterior distribution of rho
   init.beta <- mcmc.pos$beta
   init.bi <- mcmc.pos$bi
   init.ct <- mcmc.pos$ct
   init.u <- mcmc.pos$u
   init.ystar <- mcmc.pos$ystar
   init.D <- mcmc.pos$D
   init.E <- mcmc.pos$E

#****************************** Dimensions again ************************************#
   
   # number of all observations
   N <- length(Y)
   # Dimension of all fixed parameter vector beta 
   FB <- ncol(X)
   # Dimension of random fixed parameter vector bi
   RB <- ncol(W)
   # Dimension of random fixed parameter vector ci
   TB <- ncol(S)
   
#**************************** initial values **************************************#
# in the reduced runs, all the initial values are the means or medians, or
# other specific values of the posterior distributions
   total.return <- m+burnin
   if ("rho"%in% monitor){ 
     if (nlag==1){
       phi <- rep(NA, total.return)
       phi[1] <- init.phi
     }else{
       phi <- matrix(NA, nrow=total.return, ncol=nlag)
       phi[1,] <- matrix(init.phi, nrow=1, ncol=nlag)      
     }
   }else{
     if (nlag==1){
       phi <- init.phi
     }else{
       phi <- matrix(init.phi, nrow=1, ncol=nlag)
     }
   }

  if ("q.rho"%in% monitor){
     if (nlag==1){
       q.rho <- rep(NA, total.return)
       q.rho[1] <- init.phi
     }else{
       q.rho <- matrix(NA, nrow=total.return, ncol=nlag)
       q.rho[1,] <- matrix(init.phi, nrow=1, ncol=nlag)      
     }     
   }else{
     if (nlag==1){
       q.rho <- init.phi
     }else{
       q.rho <- matrix(init.phi, nrow=1, ncol=nlag)
     }
   } 
   
   if ("beta"%in% monitor){ 
     beta <- matrix(NA, nrow=total.return, ncol=FB)
     beta[1,] <- init.beta
   }else{
     beta <- matrix(init.beta, nrow=1, ncol=FB)
   }
   
   if ("bi"%in% monitor){ 
     bbb <- matrix(NA, nrow=total.return, ncol=RB*GU)
     bbb[1,] <- matrix(init.bi, nrow=1, ncol=RB*GU)
   }else{
     bbb <- matrix(init.bi, nrow=1, ncol=RB*GU)
   }
   
   if ("ct"%in% monitor){ 
    ccc <- matrix(NA, nrow=total.return, ncol=GT*TB)
    ccc[1,] <- matrix(init.ct, nrow=1, ncol=GT*TB)
   }else{
    ccc <- matrix(init.ct, nrow=1, ncol=GT*TB)
   }

   EDim <- (TB^2-TB)/2+TB   
   if ("E"%in% monitor){
     EEE <- matrix(NA, nrow=total.return, ncol=EDim)
     EEE[1,] <- matrix(0, nrow=1, ncol=EDim)
   }else{
     EEE <- matrix(0, nrow=1, ncol=EDim)
   }
   
   DDim <- (RB^2-RB)/2+RB
   if ("D"%in% monitor){
     DDD <- matrix(NA, nrow=total.return, ncol=DDim)
     DDD[1,] <- matrix(0, nrow=1, ncol=DDim)
   }else{
     DDD <- matrix(0, nrow=1, ncol=DDim)
   }
   
   if ("u"%in% monitor){
     u.aux <- matrix(NA, nrow=total.return, ncol=N)
     u.aux[1,] <- matrix(rnorm(N), nrow=1, ncol=N)
   }else{
     u.aux <- matrix(rnorm(N), nrow=1, ncol=N)
   }

   if ("ystar"%in% monitor){ 
     ystarm <- matrix(NA, nrow=total.return, ncol=N)
     ystarm[1,] <- matrix(init.ystar, nrow=1, ncol=N)
   }else{
     ystarm <- matrix(init.ystar, nrow=1, ncol=N)
   }
   
   
############### loop begins ################################################
 #Tracking the process of MCMC simulation
    mcmc <- m+burnin

    for(g in 2:mcmc) {
     if (tracking>0){
       if ((g%% tracking ==0)| g==m){
          cat("g=",g,"\n")
        }
     }
       
       
# Get the values of the pervious iterations, used in the transition

  if (g==2){
     current.bi <- matrix(bbb[1,], ncol=RB)           # current values of bi
     current.ct <- matrix(ccc[1,], ncol=TB)           # current values of ct
     current.beta <- beta[1,]                         # current values of beta
     current.ystar <- ystarm[1,]                      # the current value of ystar
     current.phi <- init.phi
   }else{
     current.bi <- new.bi                            # current values of bi
     current.ct <- new.ct                            # current values of ct
     current.beta <- new.beta                        # current values of beta
     current.ystar <- new.ystar                      # the current value of ystar
     current.phi <- new.phi
   }
  
#_______________Block 1: variance Matrix of E _____________________________
       
     if ("E"%in% fixed){
       if (ncol(S)==1){
         Ei <- init.E
       }else{
         Ei <- Wishart(element=init.E, col=TB)
       }
      }else{
       Eupdate <- Covariance.update(Covariate=S, prior.freedom=e0, prior.scale=E0,
                              dimension=TB, obs.length=GT, current.value=current.ct)
       Ei <- Eupdate[[1]]
      }
       
#_______________Block 2: variance Matrix of D _____________________________
       
   if ("D"%in% fixed){
        if (ncol(W)==1){ 
          Di <- init.D
        }else{
          Di <- Wishart(element=init.D, col=RB)
        }
    }else{
       Dupdate <- Covariance.update(Covariate=W, prior.freedom=d0, prior.scale=D0,
                              dimension=RB, obs.length=GU, current.value=current.bi)
       Di <- Dupdate[[1]]
    }
       
 #____ Compute the covariance matrix for the maximum periods of time__________
  
    if (nlag==1){
      var.phig <- 1/(1- current.phi^2)                          
      TCov <- Cov.AR1(col=GT, rho= current.phi,sigma=var.phig)
    }else{
      CCov <- Cov.ARp(phig=current.phi, col=GT)        
      TCov <- CCov[[2]]
      Sig <- TCov[1:nlag, 1:nlag]
      SSig<-solve(Sig) 
      var <- Sig[1,1]
    }

##### Variables for outputs from the loop ################################

    betaterm1 <- matrix(0, nrow=FB, ncol=FB)     # for posterior covariance of beta
    betaterm2 <- rep(0, FB)                      # for posterior mean of beta
    new.u<- rep(NA, N)
    kap <- rep(NA, GU)                           # Store the kappa for variance
    qu <- rep(NA,N)                              # Store the output of the auxiliary term
    new.ystar <- rep(NA, N)
       
# Loop in the dimension of unit
   for (n in 1: GU){                               # loop for each unit
      rows <- entry[which(unit.id==n)]                 # the obserbations of unit n
      nT <- length(rows)
      yi <- Y[rows]
      ystari <- current.ystar[rows]
      Xi <- as.matrix(X[rows,])
      Wi <- as.matrix(W[rows,])
      Si <- as.matrix(S[rows,])
      bbbn <- current.bi[n,]
      TT <- time.id[rows]
      Tccin <- as.matrix(current.ct[TT,])
      rant <- c()
      for (l in 1:nT){
        rant[l] <- Si[l,]%*%cbind(Tccin[l,])
      }

   # Updata U 
      nCov <- as.matrix(TCov[1:nT, 1:nT])
      nVK <- V.matrix(Omega=nCov, col=nT)
      nV <- nVK[[2]]                          
      nK <- nVK[[1]]
      if ("u"%in% fixed){
        u.fixed <- init.u[rows]
        u.decom <- nV%*%cbind(u.fixed)
       }else{
        Cov.U <- solve(diag(nT) + 1/nK*t(nV)%*%nV)   
        mean.U <- Cov.U %*% t(nV) %*% (ystari-Xi%*%cbind(current.beta)-Wi%*%cbind(bbbn)-rant)/nK
        u.draw <- mvrnorm(1, mu=mean.U, Sigma=Cov.U)
        new.u[rows] <- u.draw
        u.decom <- nV%*%cbind(u.draw)
      }
      qu[rows] <- u.decom
      kap[n] <- nK
  #__________________________Block 4: Updata ystar ________________________

   if ("ystar"%in% fixed){
     new.ystar <- init.ystar
    }else{
       mean.ystar <- Xi%*% cbind(current.beta)+Wi%*%cbind(bbbn)+rant+u.decom
       ystar <- sapply(c(1:nT),
                      function(r){BTnormlog(x=yi[r], mu=mean.ystar[r], sigma=sqrt(nK))})
       new.ystar[rows] <- ystar
    }
        
   #_________________ calculate BB1 and BB2 for beta __________________________
      
    if (!"beta"%in% fixed){
      var.beta <- pseudoinverse(nCov+Wi%*%solve(Di)%*%t(Wi))
      betaterm1 <- betaterm1+t(Xi) %*% var.beta %*%Xi
      betaterm2 <- betaterm2+t(Xi)%*% var.beta %*% (ystar-rant)
     }
   }
     
 #_______________Block 5 Beta and {b_i}: ______________________________________
       # 1) \beta
       
   if ("beta"%in% fixed){
     new.beta <- init.beta
   }else{
    B1<-pseudoinverse(betaterm1+solve(B0))
    beta1<-t(B1%*%(betaterm2+solve(B0)%*%beta0))
    new.beta<-mvrnorm(1, beta1, B1)
   }

   # 2) {b_i}
   if ("bi"%in% fixed){
     new.bi <- current.bi
   }else{
       new.bi <- matrix(NA, nrow=GU, ncol=RB)
   for (n in 1: GU){                                # loop for each unit
     rows <- entry[which(unit.id==n)]                    # the obserbations of unit n
     nT <- length(rows)
     nCov <- as.matrix(TCov[1:nT, 1:nT])
     ystari <- new.ystar[rows]
     Xi <- as.matrix(X[rows,])
     Wi <- as.matrix(W[rows,])
     Si <- as.matrix(S[rows,])
     TT <- time.id[rows]
     Tccin <- as.matrix(current.ct[TT,])
     rant <- c()
     for (l in 1:nT){
       rant[l] <- Si[l,]%*%cbind(Tccin[l,])
     }
     
     Cov.b <- solve(Di + t(Wi)%*%solve(nCov)%*%Wi)
     mean.b <- Cov.b%*% t(Wi)%*%solve(nCov)%*%(ystari-Xi%*% cbind(new.beta)-rant)
     if (length(mean.b)==1){
       b.drawn <- as.matrix(rnorm(1, mean.b, sqrt(Cov.b)))
     }else{
       b.drawn <- mvrnorm(1, mu=mean.b, Sigma=Cov.b)
     }
     new.bi[n,] <-  t(b.drawn)
    }
   }

  #___________________ Block 6: Time-Specific Random Effect ________________
  if ("ct"%in% fixed){
    new.ct <- current.ct
  }else{
    new.ct <-matrix(NA, nrow=GT, ncol=TB)
    for (j in 1: GT){
      row1 <- entry[which(time.id==j)]  # Those entries with t=j
      H <- length(row1)           #  #observations at t period
      iden1<-unit.id[row1]            # The units at t period
      tttn <- as.matrix(new.bi[iden1,])         # The random unit-specific effects of
      yt1 <- new.ystar[row1]
      Xt1 <- as.matrix(X[row1,])
      Wt1<- as.matrix(W[row1,])
      St1 <- as.matrix(S[row1,])
      Ut1 <- qu[row1]
      kapt1 <- kap[iden1]
      mean1 <- function(WF, TT){
        TERM<- WF%*% cbind(TT)
      }
      term7 <-  sapply(c(1:H),
             function(q){mean1(WF=Wt1[q,], TT=tttn[q,])})
      PM <- matrix(0, nrow=H, ncol=H)
      for (c in 1:H){
        PM[c,c] <- kapt1[c]
      }
      cov.S <- solve(Ei+t(St1)%*%solve(PM)%*%St1)
      mean.S <- cov.S%*% t(St1) %*%solve(PM)%*%(yt1-term7-Xt1%*%cbind(new.beta)-Ut1)
      if (length(mean.S)==1){
         supd <- as.matrix(rnorm(1, mean.S, sqrt(cov.S)))
      } else {
         supd <- rmvnorm(1, mean.S, cov.S)
      }
         new.ct[j,] <- supd
    }
  }
          
#__________________Block 7: Autoregressive Coef.s_________________________
  if (nlag==1){
     PhiPr <- PhiTerm0(Entry=entry, unit=unit.id, NN=GU, yystar=new.ystar,newb=new.bi,
                       PX=X, PW=W, PS=S, Ti=time.id, CC=new.ct, newbeta=new.beta)
     PhiSig <- PhiPr[1]
     PhiMean <- PhiPr[2]
     hatPhi<-1/PhiSig
     hatphi<-hatPhi*PhiMean
     phi.draw <- rtnorm(1, hatphi, sqrt(hatPhi), lower=-1, upper=1)
     # cat("phi.draw ", phi.draw, "\n", sep=" ")
     ####### Accept and Reject ########
     alpha<-AcceptRate0(rho=phi.draw, orho=current.phi, unit=unit.id, precision=1,
                        AY=new.ystar, AX=X, AW=W, AS=S, Abeta=new.beta, Ab=new.bi,
                        AT=time.id, Ac=new.ct, Entry=entry)

     if( runif(1) <= alpha ) {
        new.phi<-phi.draw
     }else{
        new.phi<-current.phi
     }
   }
   if (nlag>1){
      PhiPr <- PhiTermARpT(Entry=entry, unit=unit.id, NN=GU, yystar=new.ystar,
                           newb=new.bi, PX=X, PW=W, PS=S, Ti=time.id, CC=new.ct,
                           newbeta=new.beta, nlag=nlag)
      PhiSig <- PhiPr[,-(nlag+1)]
      PhiMean <- PhiPr[,(nlag+1)]
      phigg<-PhiPropARp(Ymt=PhiSig, Ymt2=PhiMean,scale=1)
         
      alpha<-AcceptRateARpT(nrho=phigg, orho=current.phi, Osig=Sig , unit=unit.id,
                            AY=new.ystar, AX=X, AW=W, AS=S, Abeta=new.beta, Ab=new.bi,
                            AT=time.id, Ac=new.ct, Entry=entry, f=g, nlag=nlag)
      #   cat("alpha ", alpha, "\n", sep=" ")
      u <- runif(1)
      if( u <= alpha ) {
        new.phi <- phigg
      } else {
        new.phi<- current.phi
      }
    }
   if (tracking>0){
    if ((g%% tracking ==0)| g==m){
      cat("phi = ", new.phi, "\n", sep=" ")
    }
   }   
    if ("rho"%in% fixed){
      new.phi <- init.phi
    } 
   
 #___________________ Input the new iteration ______________________#
         if ("ystar"%in% monitor){
           ystarm[g,] <- new.ystar
         }
         
         if ("u"%in% monitor){
           u.aux[g,] <- new.u
         }
       
         if ("bi" %in% monitor){
           bbb[g,] <- as.vector(new.bi)
         }
        
         if ("ct"%in% monitor){
           timetime <- as.vector(new.ct)
           ccc[g,] <- timetime
         }
         
         if ("beta"%in% monitor){
           beta[g,] <- new.beta
         }
       
         if ("rho"%in% monitor){
           if (nlag==1){
             phi[g] <- new.phi
           }else{
             phi[g,] <- new.phi
           }
         }
         
         if ("q.rho"%in% monitor){
           if (nlag==1){
             q.rho[g] <- phi.draw
           }else{
             q.rho[g,] <- phigg
           }
         }
       
         if ("D"%in% monitor){
           if (ncol(W)==1){
             DDD[g,] <-  as.vector(Di)
           }else{
             DDD[g,] <-  Dupdate[[2]]
           }
          }
       
         if ("E"%in% monitor){
           if (ncol(S)==1){
           EEE[g,] <- as.vector(Ei)
         }else{
           EEE[g,] <- Eupdate[[2]]
         }
         }
     }
#__________________________Returning MCMC output_______________________________  
  if ("rho"%in% monitor){ 
    if (nlag==1){
      posterior.phi <- phi[-c(1:burnin)]
     }else{
      posterior.phi <- phi[-c(1:burnin),]
     }
  }else{
        posterior.phi <- phi
  }
  if ("q.rho"%in% monitor){ 
    if (nlag==1){
      proposal.phi <- q.rho[-c(1:burnin)]
    }else{
      proposal.phi <- q.rho[-c(1:burnin),]
     }
  }else{
     proposal.phi<- q.rho
  }
   
  if ("beta"%in% monitor){ 
     posterior.beta <- beta[-c(1:burnin),]
  }else{
     posterior.beta <- beta
   }
  if ("bi"%in% monitor){ 
     posterior.unit.res <- bbb[-c(1:burnin),]
  }else{
     posterior.unit.res <- bbb
  }
 if ("ct"%in% monitor){ 
     posterior.time.res <- ccc[-c(1:burnin),]
  }else{
     posterior.time.res <- ccc
  }   
  if ("E"%in% monitor){ 
     posterior.time.cov <- EEE[-c(1:burnin),]
  }else{
     posterior.time.cov <- EEE
  }    
  if ("D"%in% monitor){ 
     posterior.unit.cov <- DDD[-c(1:burnin),]
  }else{
     posterior.unit.cov <- DDD
  }
   if ("u"%in% monitor){
      posterior.u.aux <- u.aux[-c(1:burnin),]
  }else{
      posterior.u.aux <- u.aux
  }    
   if ("ystar"%in% monitor){ 
     posterior.data <- ystarm[-c(1:burnin),]
  }else{
     posterior.data <- ystarm
  }

   return(list(Phi=posterior.phi, proposal=proposal.phi, Beta=posterior.beta,  bi=posterior.unit.res,
               ct=posterior.time.res, var.c=posterior.time.cov, Cov.betai=posterior.unit.cov,
               auxiliary.u=posterior.u.aux, data.augmented=posterior.data))
 }

