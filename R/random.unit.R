#######################################################################################
# File: RandomCoef                                                            
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
#    unitint     add an unit-specific intercept? "TRUE" or "FALSE". Default is "FALSE"
#      model     mcmc output returned by GLMMARp.Binary() function
######################################################################################
 random.unit <- function(y, X1, Wi, St, Ai, Ft, Unit.Index, Time.Index, unitint,  model)
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
   W <- as.matrix(NewData[[5]])
   RB <- ncol(W)
   betai <- as.matrix(model[[2]][, (ncol(X1)+1):(ncol(X1)+ncol(Wi)*ncol(AA))])
   bi <- model[[3]]
   beta.unit <- VaryingCoef(KK=AA, Betai=betai, Residual=bi, m=nrow(betai),
                            RB=RB, GU=length(unique(newid1)),
                            unitint.add=unitint)
  return(beta.unit)
 }


 


