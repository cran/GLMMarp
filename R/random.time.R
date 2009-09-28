#######################################################################################
# File: RandomCoef2                                                            
# Author: Xun Pang                                                             
# Creation Date: 1/1/2009                                                            
# Description: 	Function to compute the random coefficients in the structural form
#               i.e., beta_{2i} and beta_{3t_i}
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
#    timeint    add an time-specific intercept? "TRUE" or "FALSE". Default is "FALSE"
#      model     mcmc output returned by GLMMARp.Binary() function
######################################################################################
  random.time <- function(y, X1, Wi, St, Ai, Ft, Unit.Index, Time.Index, timeint, model)
 {
  NewData <- RegressionData(yob=y, Unit=Unit.Index, Time=Time.Index, Fixed=X1,
                               UnitRandom=Wi, TimeRandom=St,
                               UnitPred=Ai, TimePred=Ft)

  #Get results from the approprate function
  newid1 <- NewData[[1]]
  newid2 <- NewData[[2]]
  id1=Unit.Index; id2=Time.Index
  # to compute the unit-specific effects
   UnitP <-GroupPred(X=as.matrix(Ai), index=id1)
   A <- as.matrix(UnitP)
   AA <- A[sort(unique(newid1)),]
  # to compute the time-specific effects
   UnitQ <-GroupPred(X=as.matrix(Ft), index=id2) 
   B <- as.matrix(UnitQ)
   B <- as.matrix(B[unique(id2)%in%unique(newid2),])
   BB <- as.matrix(B[sort(unique(Index(newid2))),])
   betat <- as.matrix(model[[2]][, (ncol(X1)+ncol(Wi)*ncol(AA)+1):ncol(model[[2]])])
   ct <- model[[4]]
   beta.time <- VaryingCoef(KK=BB, Betai=betat, Residual=ct, m=nrow(betat),
                            RB=ncol(St), GU=length(unique(newid2)),
                            unitint.add=timeint)
  return(beta.time)
 }
