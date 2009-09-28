#######################################################################################
# File: GLMMARp.R                                                             
# Author: Xun Pang                                                             
# Creation Date: 1/1/2009                                                            
# Description: 	MCMC simulation for estimating binary GLMM-AR(p) model
#               based on Xun Pang, 2009, Modeling Intertemporal and Contemporal 
#               Dependence in Binary TSCS Data: A Bayesian Model with AR(p) Errors and 
#               Non-nested Clustering1, http://xpang.wustl.edu/Dissertation.html
#
# Usage:        GLMMARp.Binary(y,X1, Wi, St, Ai="NULL", Ft="NULL", ...)
# Arguments:
#         y     response variable, dichotomous (0/1)
#        X1     matrix of covariates with fixed effects
#        Wi     matrix of covariates with subject-varying effects
#        St     matrix of covariates with time-varying effects
#        Ai     matrix of covariates explaining the subject-varying effects
#               with the same length of y (the values of time-invariant variables
#               have to be repeated over the same times of each subject,
#               for time-varying covariate in Ai, the function will automatically
#               use the within-subject mean). The default is "NULL", no group-level
#               covariates
#        Ft     matrix of covariates explaining the time-varying effects
#               with the same length of y (the function will automatically
#               use the within-time-period mean). The default is "NULL", no group-level
#               covariates
# Unit.Index    index of subjects
# Time.Index    index of time periods
# timeint.add   add an time-specific intercept? "TRUE" or "FALSE". Default is "FALSE"
# unitint.add   add an unit-specific intercept? "TRUE" or "FALSE". Default is "FALSE"
#         m     number of iterations to be returned
#    burnin     number of initial iterations which are discarded
#    priors     beta0--- mean vector of the multivariate normal distribution of beta,
#                        all the fixed effects (beta in the reduced form of the model)
#               B0--- variance matrix of the multivariate normal distribution of beta,
#                        all the fixed effects (beta in the reduced form of the model)
#               d0, D0 --- degree of freedom and scale matrix of the Wishart prior
#                          on b_i, the subject-level residual
#               e0, E0 --- degree of freedom and scale matrix of the Wishart prior
#                          on c_t, the time-level residual
#     nlag      the lag order p, integer >= 1
# initial values
#               init.phi --- vector of initial values of the autoregressive coef.
#               init.beta, init.b, and init.c --- initial values of the fixed and
#                   random coefficients
#  
#      pgm      apply parameter expantion method? Default is TRUE
# tracking      every "tacking" iterations, return the information about how many
#               iterations in total have been done
# marginal.likelihood
#               return the marginal likelihood, default is FALSE
#  monitor      a string contains the names of parameters whose MCMC output are
#               will be returned. The string has to be a subset of
#               ("rho", "beta", "bi", "ct", "D", "E", "ystar", "u") which is the default
#     thin      thin the chain by recording each "thin" iterations
#                                       
######################################################################################
  
 GLMMARp.Binary <- function(y,X1, Wi, St, Ai="NULL", Ft="NULL", Unit.Index, Time.Index,
                            timeint.add=FALSE, unitint.add=FALSE, m=10000, burnin=5000,
                            beta0,B0, D0, d0, E0, e0, nlag, init.phi,  init.beta=0,
                            init.b=0, init.c=0, pgm=TRUE, tracking=0,
                            marginal.likelihood=FALSE, reduced.mcmc=1000, reduced.burnin=100,
                            monitor=c("rho", "beta", "bi", "ct", "D", "E", "ystar", "u"),
                            thin=1){
    
    ## y \in {0, 1} error checking
     yy <- y[!is.na(y)]
     if (sum(yy!=0 & yy!=1) > 0) {
       cat("Elements of Y equal to something other than 0 or 1.\n")
       stop("Check data and call GLMMARp.Binary() again. \n") 
     }
     #check.mcmc.parameters (burnin=burnin, mcmc=m, thin=thin)
     
    #****************** Arrange Data for MCMC updating  **********************#
     NewData <- RegressionData(yob=y, Unit=Unit.Index, Time=Time.Index, Fixed=X1,
                               UnitRandom=Wi, TimeRandom=St,
                               UnitPred=Ai, TimePred=Ft)

    #Get results from the approprate function
    newid1 <- NewData[[1]]
    newid2 <- NewData[[2]]
    Y <- NewData[[3]]
    X <- as.matrix(NewData[[4]])
    W <- as.matrix(NewData[[5]])
    if (unitint.add==TRUE){
     W <- cbind(1, W)
    }
    S <- as.matrix(NewData[[6]])
     if (timeint.add==TRUE){
      S <- cbind(1, S)
    } 
   # number of all observations
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
   id1 <- Index(newid1)                 # unit id repeating t times t=# observations
   id2 <- (newid2)-min(newid2)+1        # year id repeating n times n=# unit in each time period
   GU <- max(id1)                       # how many units
   GT <- max(id2)                       # how many time periods in total
   #****************** Storage and Initial Values ********************
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
   
   if ("beta"%in% monitor){ 
     beta <- matrix(NA, nrow=total.return, ncol=FB)
     beta[1,] <- init.beta
   }else{
     beta <- matrix(init.beta, nrow=1, ncol=FB)
   }
   
   if ("bi"%in% monitor){ 
     bbb <- matrix(NA, nrow=total.return, ncol=RB*GU)
     bbb[1,] <- matrix(init.b, nrow=1, ncol=RB*GU)
   }else{
     bbb <- matrix(init.b, nrow=1, ncol=RB*GU)
   }
   
   if ("ct"%in% monitor){ 
    ccc <- matrix(NA, nrow=total.return, ncol=GT*TB)
    ccc[1,] <- matrix(init.c, nrow=1, ncol=GT*TB)
   }else{
    ccc <- matrix(init.c, nrow=1, ncol=GT*TB)
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
     ystarm[1,] <- matrix(0, nrow=1, ncol=N)
   }else{
     ystarm <- matrix(0, nrow=1, ncol=N)
   }
   
############### MCMC Updating ######################################

 #Tracking the process of MCMC simulation
    mcmc <- (m+burnin)*thin
    for(g in 2:(mcmc)) {
      if (tracking>0){
       if ((g%% tracking ==0)| g==(m*thin)){
          cat("g=",g,"\n")
        }
     }
       
# The current values of those parameters

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

  # MCMC Updating 1: Covariance Matrices (random errors)
     Eupdate <- Covariance.update(Covariate=S, prior.freedom=e0, prior.scale=E0,
                              dimension=TB, obs.length=GT, current.value=current.ct)
     Ei <- Eupdate[[1]]
     Dupdate <- Covariance.update(Covariate=W, prior.freedom=d0, prior.scale=D0,
                              dimension=RB, obs.length=GU, current.value=current.bi)
     Di <- Dupdate[[1]]
  # Compute the covariance-matrix of xi
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

  # Variables for outputs from the loop      
    betaterm1 <- matrix(0, nrow=FB, ncol=FB)     # for posterior covariance of beta
    betaterm2 <- rep(0, FB)                      # for posterior mean of beta
    new.u<- rep(NA, N)
    kap <- rep(NA, GU)                           # Store the kappa for variance
    qu <- rep(NA,N)                              # Store the output of the auxiliary term
    new.ystar <- rep(NA, N)

    for (n in 1: GU){                               # loop for each unit
      rows <- entry[which(id1==n)]                 # the obserbations of unit n
      nT <- length(rows)
      yi <- Y[rows]
      ystari <- current.ystar[rows]
      Xi <- as.matrix(X[rows,])
      Wi <- as.matrix(W[rows,])
      Si <- as.matrix(S[rows,])
      bbbn <- current.bi[n,]
      TT <- id2[rows]
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
      Cov.U <- solve(diag(nT) + 1/nK*t(nV)%*%nV)   
      mean.U <- Cov.U %*% t(nV) %*% (ystari-Xi%*%cbind(current.beta)-Wi%*%cbind(bbbn)-rant)/nK
      u.draw <- mvrnorm(1, mu=mean.U, Sigma=Cov.U)
      new.u[rows] <- u.draw
      u.decom <- nV%*%cbind(u.draw)
      qu[rows] <- u.decom
      kap[n] <- nK
      
  # Updata ystar
    mean.ystar <- Xi%*% cbind(current.beta)+Wi%*%cbind(bbbn)+rant+u.decom
    ystar <- sapply(c(1:nT),
             function(r){BTnormlog(x=yi[r], mu=mean.ystar[r], sigma=sqrt(nK))})
 
    new.ystar[rows] <- ystar

   # calculate BB1 and BB2 for beta
    var.beta <- pseudoinverse(nCov+Wi%*%solve(Di)%*%t(Wi))
    betaterm1 <- betaterm1+t(Xi) %*% var.beta %*%Xi
    betaterm2 <- betaterm2+t(Xi)%*% var.beta %*% (ystar-rant)
  }
       
 # Update Beta and {b_i} in one block
       
   # 1) \beta
   B1<-pseudoinverse(betaterm1+solve(B0))
   beta1<-t(B1%*%(betaterm2+solve(B0)%*%beta0))
   new.beta<-mvrnorm(1, beta1, B1)
   # 2) {b_i}
   new.bi <- matrix(NA, nrow=GU, ncol=RB)
   for (n in 1: GU){                                # loop for each unit
     rows <- entry[which(id1==n)]                    # the obserbations of unit n
     nT <- length(rows)
     nCov <- as.matrix(TCov[1:nT, 1:nT])
     ystari <- new.ystar[rows]
     Xi <- as.matrix(X[rows,])
     Wi <- as.matrix(W[rows,])
     Si <- as.matrix(S[rows,])
     TT <- id2[rows]
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
       
  # Time-Specific Random Effect 
    new.ct <-matrix(NA, nrow=GT, ncol=TB)
    for (j in 1: GT){
      row1 <- entry[which(id2==j)]  # Those entries with t=j
      H <- length(row1)           #  #observations at t period
      iden1<-id1[row1]            # The units at t period
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

       # End of Block 6
   # Update: Autoregressive Coef.s

   if (nlag==1){
     PhiPr <- PhiTerm0(Entry=entry, unit=id1, NN=GU, yystar=new.ystar,newb=new.bi,
                       PX=X, PW=W, PS=S, Ti=id2, CC=new.ct, newbeta=new.beta)
     PhiSig <- PhiPr[1]
     PhiMean <- PhiPr[2]
     hatPhi<-1/PhiSig
     hatphi<-hatPhi*PhiMean
     phi.draw <- rtnorm(1, hatphi, sqrt(hatPhi), lower=-1, upper=1)
     # cat("phi.draw ", phi.draw, "\n", sep=" ")
     ####### Accept and Reject ########
     alpha<-AcceptRate0(rho=phi.draw, orho=current.phi, unit=id1, precision=1,
                        AY=new.ystar, AX=X, AW=W, AS=S, Abeta=new.beta, Ab=new.bi,
                        AT=id2, Ac=new.ct, Entry=entry)

     if( runif(1) <= alpha ) {
        new.phi<-phi.draw
     }else{
        new.phi<-current.phi
     }
   }
   if (nlag>1){
      PhiPr <- PhiTermARpT(Entry=entry, unit=id1, NN=GU, yystar=new.ystar,
                           newb=new.bi, PX=X, PW=W, PS=S, Ti=id2, CC=new.ct,
                           newbeta=new.beta, nlag=nlag)
      PhiSig <- PhiPr[,-(nlag+1)]
      PhiMean <- PhiPr[,(nlag+1)]
      phigg<-PhiPropARp(Ymt=PhiSig, Ymt2=PhiMean,scale=1)

         
      alpha<-AcceptRateARpT(nrho=phigg, orho=current.phi, Osig=Sig , unit=id1,
                            AY=new.ystar, AX=X, AW=W, AS=S, Abeta=new.beta, Ab=new.bi,
                            AT=id2, Ac=new.ct, Entry=entry, f=g, nlag=nlag)
      #   cat("alpha ", alpha, "\n", sep=" ")
      u <- runif(1)
      if( u <= alpha ) {
        new.phi <- phigg
      } else {
        new.phi<- current.phi
      }
    }
    if (tracking>0){
     if ((g%% tracking ==0)| g==(m*thin)){
       cat("phi = ", new.phi, "\n", sep=" ")
     }
    }
       
 #  PGM Updating
       
   if (pgm==TRUE){
      if (any(beta0!=0)){
         warning("priors have to be centered at 0 for using PGM-MGMC updating")
         return(NULL)
      }
      phig <-new.phi   
    if (nlag==1){
      Phig <- 1/(1-phig^2)
      TCov <- Cov.AR1(col=GT, rho=phig,sigma=Phig)
    }else{
      CCov <- Cov.ARp(phig=phig, col=GT)        
      TCov <- CCov[[2]]
    }
           
    qq1 <- 0
    qq2 <- 0
    for (f in 1: GU){
      row2 <- entry[which(id1==f)]       # Those entries with n=i
      H <- length(row2)               
      iden3 <-id2[row2]                  # The units at t period
      tttn <- new.bi[f,]                 # The random unit-specific effects
      yt1 <- new.ystar[row2]
      Xt1 <- as.matrix(X[row2,])
      Wt1 <- as.matrix(W[row2,])
      St1 <- as.matrix(S[row2,])
      nCov <- as.matrix(TCov[1:H, 1:H])
      timec <- as.matrix(new.ct[iden3,])
      rrant <- c()
      for (w in 1:H){
        rrant[w] <- St1[w,]%*% cbind(timec[w,])
      }
      termm <- yt1-Xt1%*%as.matrix(new.beta)-Wt1%*%as.matrix(tttn)-rrant
      qq1 <- qq1+t(termm)%*%solve(nCov)%*%termm
      qq2 <- qq2+rbind(tttn)%*%Di%*%cbind(tttn)
     }
           
     qq3 <- rbind(new.beta)%*% solve(B0)%*% cbind(new.beta)
     qq4 <- 0
     for (t in 1: GT){
       qq4 <- qq4+rbind(new.ct[t,])%*%Ei%*%cbind(new.ct[t,])
     }
     bb <- (qq1+qq2+qq3+qq4)/2
     aa <- (N+ ncol(beta)+ncol(bbb)+ncol(ccc)+1)/2
     gamm <- sqrt(rgamma(1, shape=aa, scale=1/bb))
     new.ystar <- new.ystar*gamm
     new.beta <- new.beta*gamm
     new.bi <- new.bi*gamm
     new.ct <- new.ct*gamm

     if (tracking>0){
      if ((g%% tracking ==0)| g==(m*thin)){ 
        cat("gamm ",gamm, "\n", sep=" ")
      }
    }
   }
   #*********************** Update other parameters ************************************

      if ((g%%thin ==0)| g==mcmc){
         gg <- g%/%thin
         if ("ystar"%in% monitor){
           ystarm[gg,] <- new.ystar
         }
         
         if ("u"%in% monitor){
           u.aux[gg,] <- new.u
         }
       
         if ("bi" %in% monitor){
           bbb[gg,] <- as.vector(new.bi)
         }
        
         if ("ct"%in% monitor){
           timetime <- as.vector(new.ct)
           ccc[gg,] <- timetime
         }
         
         if ("beta"%in% monitor){
           beta[gg,] <- new.beta
         }
       
         if ("rho"%in% monitor){
           if (nlag==1){
             phi[gg] <- new.phi
           }else{
             phi[gg,] <- new.phi
           }
         }
       
         if ("D"%in% monitor){
           if (ncol(W)==1){
             DDD[gg,] <-  as.vector(Di)
           }else{
             DDD[gg,] <-  Dupdate[[2]]
           }
          }
       
         if ("E"%in% monitor){
           if (ncol(S)==1){
           EEE[gg,] <- as.vector(Ei)
         }else{
           EEE[gg,] <- Eupdate[[2]]
           }
         }
       }
     }
#__________________________Returning MCMC output_______________________________
 cat("MCMC simulation completed", "\n", sep=" ")
 cat("Returning MCMC output...", "\n", sep=" ")  
  if ("rho"%in% monitor){ 
    if (nlag==1){posterior.phi <- phi[-c(1:burnin)]
             }else{
        posterior.phi <- phi[-c(1:burnin),]
     }
  }else{
        posterior.phi <- phi
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
      posterior.u.aux <- u.aux[-c(1:(burnin)),]
  }else{
      posterior.u.aux <- u.aux
  }    
   if ("ystar"%in% monitor){ 
     posterior.data <- ystarm[-c(1:(burnin)),]
  }else{
     posterior.data <- ystarm
  }
  mcmc.output <- list(Phi=posterior.phi, Beta=posterior.beta,  bi=posterior.unit.res,
                      ct=posterior.time.res, var.c=posterior.time.cov, Cov.betai=posterior.unit.cov,
                      auxiliary.u=posterior.u.aux,data.augmented=posterior.data)
     
  if  (marginal.likelihood==TRUE){
    cat("Now computing marginal likelihood...", "\n", sep=" ")
    Log.Marg.Lik <- MargLike.Binary.Builtin(Y=Y,X=X, W=W, S=S,id1=id1,id2=id2, reduced.mcmc=reduced.mcmc, 
                                     reduced.burnin=reduced.burnin, mcmc.output=mcmc.output, nlag=nlag, beta0=beta0,
                                     B0=B0, D0=D0, d0=d0, E0=E0, e0=e0, tracking=tracking)
    return(list(mcmc.output=mcmc.output, Log.Marginal.Likelihood=Log.Marg.Lik))
   }else{
    return(mcmc.output)
   }
 }
    

########################### END ###################################