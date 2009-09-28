\name{Marg.Like.Binary}
\alias{Marg.Like.Binary}

\title{Compute the Marginal Likelihood of the GLMM-AR(p) Model by Using the MCMC Output and Reduced Samples}
\description{
This function estimates the marginal likelihood for model comparison by using the Bayes factor. It is built in the \code{GLMMARp.Binary} function, but can also be used independently, and the way to fit in the function is similar to how to fit in the \code{GLMMARp.Binary} function.
}

\usage{
 Marg.Like.Binary(y,X1, Wi, St, Ai="NULL", Ft="NULL", Unit.Index, 
                 Time.Index, timeint.add=FALSE,unitint.add=FALSE, 
                 mcmc.output, reduced.mcmc, reduced.burnin, nlag,
                 beta0,B0, D0 , d0, E0, e0, tracking)
}

\arguments{
    \item{y}{A vector of response variable, dichotomous (0/1).}

    \item{X1}{A matrix of covariates with fixed effects.}

    \item{Wi}{A matrix of covariates with subject-varying effects.}

    \item{St}{A matrix of covariates with time-varying effects.}

    \item{Ai}{A matrix of covariates explaining the subject-varying effects
               with the same length of y (the values of time-invariant variables
               have to be repeated over the same times of each subject,
               for time-varying covariate in Ai, the function will automatically
               use the within-subject mean). The default is "NULL", no group-level
               covariates.}
   
    \item{Ft}{matrix of covariates explaining the time-varying effects
               with the same length of y (the function will automatically
               use the within-time-period mean). The default is "NULL", 
               no group-level covariates.}

    \item{Unit.Index}{A vector of subject index, i's. Note: the number of
                observations of each unit should be larger than the lag order, nlag. Those
                units which have fewer than or equal to nlag observations should be taken
                out of the sample in order to use the function.}

    \item{Time.Index}{A vector of time index, t's. Note: no missing observations
               in the middle of the sample time periods of a unit are allowed. In other words,
               unbalanced data structures are allowed, but no broken data structure.}

    \item{timeint.add}{Should a time-specific intercept be added into the model? 
              It takes two values: TRUE or FALSE with default as FALSE.}

    \item{unitint.add}{Should a subject-specific intercept be added into the model? 
              It takes two values: TRUE or FALSE with default as FALSE.}

    \item{reduced.mcmc}{The number of iterations to return in the reduced mcmc simulations.}

    \item{reduced.burnin}{The number of burn-in iterations for the sampler in the reduced
                     mcmc simulations.}
 
    \item{mcmc.output}{The mcmc output from the full simulation of GLMM-AR(p) model.
                     The format has to be the same as in the GLMMARp.Binary() function.}

    \item{beta0}{The prior mean of \eqn{\beta}{beta}.  This should be a 
          vector with dimension equal to the number of fixed effect parameters 
          in the reduced form. Since the dimension is difficult for the user to
          compute when the model contains multiple random coefficients and multiple 
          group-level predictors, the function will provide the correct dimension 
          in the error message if the dimension of beta0 is specified incorrectly, 
          and the user can respecify beta0 with this information and recall the 
          function. No default is provided.}
  
    \item{B0}{The prior covariance matrix of \eqn{\beta}{beta}.  This should be a 
          positive definite matrix with dimensions equal to the number of betas in 
          the reduced form of GLMM-AR(p). No default is provided.}
  
    \item{d0, D0}{The degree of freedom and scale matrix of the Wishart prior
                  on \eqn{b_i} which is the subject-level residual. D0 should 
                  not be too defuse, otherwise it may takes a long time for 
                  the chain to converge. Recommended values is a*I, where a is 
                  between 1.67 and 10. No default is provided.}

    \item{e0, E0}{The degree of freedom and scale matrix of the Wishart prior
                  on \eqn{c_t} which is the time-level residual. E0 should 
                  not be too defuse, otherwise it may takes a long time for 
                  the chain to converge. Recommended values is a*I, where a is 
                  between 1.67 and 10. No default is provided (priors have to 
                  be the same as are used in the full mcmc. }

    \item{nlag}{A scalar of the lag order p, which should be an integer equal to 
                or greater than 1. In this version, the function does not support 
                the model with p=0, which can be estimated by using BUGs or JAGs.}

    \item{tracking}{The tracking interval used in the simulation. Every "tacking" 
              iterations, the function will return the information about how many 
              iterations in total have been done.}
  
}

\value{
  A scalar which is the marginal likelihood (on a natural logarithm scale).   
}

\examples{
\dontrun{
require(panel)
require(bayesSurv)
data(StateAR3)
data(StateAR2)
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
 TimePred <- cbind(StateFailure$dis1)

TimePred <- "NULL"
marginalLikely <- Marg.Like.Binary(y=y,X1=Fixed, Wi=UnitRandom, St=TimeRandom, 
                      Ai=UnitPred, Ft=TimePred, Unit.Index=Unit, Time.Index=Time, 
                      timeint.add=0, unitint.add=1, mcmc.output=StateAR2, reduced.mcmc=100,                     
                      reduced.burnin=50, nlag=2, beta0=rep(0,15),
                      B0=diag(10, 15), D0=diag(1.67, 2) , d0=9,  E0=6, e0=6, 
                      tracking=100)
}
}

\concept{Bayes factor}
\concept{model comparison}

\seealso{\code{\link{GLMMARp.Binary}}}
\keyword{models}
