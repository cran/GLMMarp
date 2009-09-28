#******************************************************************************
# File: predtion.R                                                           
# Author: Xun Pang                                                             
# Creation Date: 1/1/2009                                                            
# Description: 	R function for compute predictive probability
#
# Usage:        Prediction.Binary(y,X1, Wi, St, Ai="NULL", Ft="NULL", ...)
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
# mcmc.output   mcmc output from GLMMARp.Binary()
#     nlag      the lag order p, integer >= 1
#                                       
######################################################################################

Prediction.Binary <- function(y,X1, Wi, St, Ai="NULL", Ft="NULL", Unit.Index, Time.Index,
                              timeint.add=FALSE,unitint.add=FALSE, mcmc.output,
                              nlag){

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
  # Dimension of random fixed parameter vector bi
   RB <- ncol(W)
   # Dimension of random fixed parameter vector ci
   TB <- ncol(S)
  # create indices for two non-nested groupings ##########
   entry <- seq(1: N)                   # the global length of observations
   id1 <- Index(newid1)                 # unit id repeating t times t=# observations
   id2 <- (newid2)-min(newid2)+1        # year id repeating n times n=# unit in each time period
   GU <- max(id1)                       # how many units
   GT <- max(id2)                       # how many time periods in total
   phi <- as.matrix(mcmc.output[[1]])
   bi <- as.matrix(mcmc.output[[3]])
   ct <- as.matrix(mcmc.output[[4]])
   ystar <- as.matrix(mcmc.output[[8]])
   u <- as.matrix(mcmc.output[[7]])
   beta <- as.matrix(mcmc.output[[2]])
 
   m <- nrow(beta)

   prob.ystarm <- matrix(NA, nrow=1, ncol=N)
   for (j in 1:m){
    current.phi <- as.matrix(phi[j,])
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

     current.beta <- beta[j,]
     current.bi <- matrix(bi[j,], ncol=RB)           # current values of bi
     current.ct <- matrix(ct[j,], ncol=TB)  
     current.u <- u[j,]
     prob.ystar <- rep(NA, N)
     
     for (n in 1: GU){                               # loop for each unit
       rows <- entry[which(id1==n)]                 # the obserbations of unit n
       nT <- length(rows)
       yi <- Y[rows]
       Xi <- as.matrix(X[rows,])
       Wi <- as.matrix(W[rows,])
       Si <- as.matrix(S[rows,])
       bbbn <- current.bi[n,]
       u.draw <- current.u[rows]
       TT <- id2[rows]
       Tccin <- as.matrix(current.ct[TT,])
     
       rant <- c()
       for (l in 1:nT){
         rant[l] <- Si[l,]%*%cbind(Tccin[l,])
       }
     
   #__________________________Block 3:Updata U _____________________________
    
       colcov <- rows-min(rows)+1 
       nCov <- as.matrix(TCov[colcov, colcov])
       nVK <- V.matrix(Omega=nCov, col=nT)
       nV <- nVK[[2]]                          
       nK <- nVK[[1]]                  
       u.decom <- nV%*%cbind(u.draw)
      
  #__________________________Block 4: Updata ystar ________________________
       mean.ystar <- Xi%*% cbind(current.beta)+Wi%*%cbind(bbbn)+rant+u.decom
       proystar <- c()
       sigma <- sqrt(nK)
       prob.yy <- function(sigma, mu){
         ystar <- rnorm(1, mu, sigma)
         pro.ystar <- pnorm(ystar, 0, sigma)
       }
       proystar <- sapply(c(1:nT),
                          function(r){prob.yy(sigma=sigma, mu=mean.ystar[r])})       
       prob.ystar[rows] <- proystar
      }
       prob.ystarm <- rbind( prob.ystarm, prob.ystar)
     }
     return(prob.ystarm[-1,])
   }

 
