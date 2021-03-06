\name{Prediction.Binary}
\alias{Prediction.Binary}

\title{Within-Sample and Out-of-Sample Prediction Based on the Posteriors Returned from 
       GLMMARp.Binary}
\description{
This functions conducts Bayes prediction by using the posterior distributions of 
the parameters simulated with the Markov Chain Monte Carlo method.
}
\usage{
Prediction.Binary(y,X1, Wi, St, Ai="NULL", Ft="NULL", Unit.Index, Time.Index,
                 timeint.add=FALSE,unitint.add=FALSE, mcmc.output,nlag)
}

\arguments{
    \item{y}{A vector of the dichotomous responses, taking two unique values of 0 and 1.}

    \item{X1}{A matrix of covariates with fixed effects.}

    \item{Wi}{A matrix of covariates with subject-varying effects.}

    \item{St}{A matrix of covariates with time-varying effects.}

    \item{Ai}{A matrix of group-level predictors explaining the subject-varying effects.
               The number of rows of this matrix is equal to the length of y. If some or 
               all predictors are time-varying, the function will automatically use the 
               within-subject mean; for time-invariant predictors, they should be repeated 
               over the same times of each subject. The default is "NULL", and there are 
               no group-level covariates.}
   
    \item{Ft}{A matrix of group-level predictors explaining the time-varying effects.
	       The number of rows of this matrix is equal to the length of y. If some or
               all predictors are subject-varying, the function will automatically use the
               within-time-period mean; for subject-invariant predictors, they should be repeated
               over the same number of subjects in the specific time period.
               The default is "NULL", and there are no group-level covariates.}

    \item{Unit.Index}{A vector of the subject index,i.e., the i's. Note: the number of
                observations of each unit should be larger than the lag order, nlag. Those
                units which have fewer than or equal to nlag observations should be taken
                out of the sample in order to use the function.}

    \item{Time.Index}{A vector of the time index, i.e., the t's.  Note: no missing observations
               in the middle of the sample time periods of a unit are allowed. In other words,
               unbalanced data structures are allowed, but no broken data structure.}

    \item{timeint.add}{Should a time-specific intercept be added into the model? 
              It takes two values: TRUE or FALSE with default as FALSE. Be sure that 
              there is no time-specific interceot already in the reduced form of the 
              model before adding an intercept here. }

    \item{unitint.add}{Should a subject-specific intercept be added into the model? 
              It takes two values: TRUE or FALSE with default as FALSE. Be sure that
              there is no subject-specific interceot already in the reduced form of the
              model before adding an intercept here.}

    \item{mcmc.output}{The MCMC output generated by the full MCMC by using the function 
                of \code{GLMMARp.Binary}. takes two values: TRUE or FALSE with default as FALSE.}

    \item{nlag}{A scalar of the lag order p, which should be an integer equal to 
                or greater than 1.}

}

\value{
   A matrix with the row dimension equal to the number of iterations returned by 
   \code{GLMMARp.Binary}; with the column dimension equal to the length of the response
   variable. Therefore, each column forms a predictive distribution of one observation.   
   The user can convert the output into an mcmc output by using the \code{coda} package, and 
   summarize the predicted distribution with various \code{coda} functions.
}

\examples{
\dontrun{
## Example
 require(panel)
 require(bayesSurv)
 data(StateFailure)
 y <- StateFailure$failure
 Unit <- StateFailure$country
 Time <- StateFailure$year
 Fixed <- cbind(StateFailure$poldemoc,  StateFailure$bnkv123, StateFailure$bnkv117,
              StateFailure$poldurab, StateFailure$faocalry, StateFailure$pwtcgdp,
              StateFailure$macnac,StateFailure$macnciv,  StateFailure$wdiinfmt,
              StateFailure$ptsamnes, StateFailure$dis1, StateFailure$bnkv81,
              StateFailure$change.demo)

 UnitRandom <- cbind(log(StateFailure$pwtopen))
 TimeRandom <- as.matrix(rep( 1, length(y)))
 UnitPred <- cbind(StateFailure$macnac, StateFailure$poldemoc)
 TimePred <- "NULL"
 
 data(StateAR2)
 prediction <- Prediction.Binary(y=y, X1=Fixed, Wi=UnitRandom, St=TimeRandom, Ai=UnitPred,
                           Ft=TimePred, Unit.Index=Unit, Time.Index=Time, timeint.add=0, 
                           unitint.add=1, mcmc.output=StateAR2, nlag=2)
 }
}
\keyword{manip}
