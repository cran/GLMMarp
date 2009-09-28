\name{GLMMARp.Binary}
\alias{GLMMARp.Binary}
\title{Markov Chain Monte Carlo for a Probit Multilevel Model with a pth-order Autoregressive Error Process}

\description{
  This function generates a sample from the posterior distribution of
  a probit multilevel model with a pth-order autoregressive error process for 
  modeling binary Time-Series Cross-Sectional responses. The model allows observed 
  and unobserved heterogeneities in both the time and spatial dimensions as well as 
  group-level predictors to explain the group-level variations. Serial correlation 
  is corrected with a pth-order autoregressive error process, and the lag order p 
  should be chosen through model comparison with the Bayes factor. 
  The lag order p must be equal or greater than 1, and the algorithm is based on Pang (2009). 
  The user supplies data and priors, and takes the responsibility of specifying
  different groups of regressors in the structural form of the GLMM-AR(p) model (not the reduced 
  form). The user can choose a subset of parameters to be monitored in the 
  MCMC process. But for the Bayes factor be calculated, all parameters (including 
  the augmented data and the posteriors of the auxiliary variable \eqn{u}{u}) should 
  be monitored. A sample of the posterior distribution is returned as a list 
  containing the individual parameter chain. The output can be easily converted 
  into mcmc output and then analyzed with functions provided in the 
  \code{coda} package. The function does not provide convergence diagnostics,
  and it is the user's responsibility to check convergence 
  by using the tools provided in the \code{coda} package.
  }
  
\usage{
 GLMMARp.Binary(y,X1, Wi, St, Ai="NULL", Ft="NULL", Unit.Index, Time.Index, 
     timeint.add=FALSE, unitint.add=FALSE, m=10000, burnin=5000,
     beta0,B0, D0, d0, E0, e0, nlag, init.phi,  init.beta=0,
     init.b=0, init.c=0,  pgm=TRUE, tracking=0,
     marginal.likelihood=FALSE, reduced.mcmc=1000, reduced.burnin=100,
     monitor=c("rho", "beta", "bi", "ct", "D", "E", "ystar", "u"),thin=1)
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

    \item{Time.Index}{A vector of the time index, i.e., the t's. Note: no missing observations 
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

    \item{m}{The number of iterations after burn-in and to be returned. If the chain is 
             thinned, thinning does not affect the number of iterations to be returned.}

    \item{burnin}{The number of burn-in iterations for the sampler. If the 
                chain is thinned, the number of burn-in interations will 
                be burnin*thin.}

    \item{thin}{The thinning interval used in the simulation. It does not affect the 
                total number of iterations set by the argument "m", but will affect 
                the total number of iterations in the simulation.}

    \item{beta0}{The prior mean of \eqn{\beta}{beta}.  This should be a 
          vector with dimension equal to the number of fixed-effect parameters 
          in the reduced form. Since the dimension is difficult for the user to
          compute when the model contains multiple random coefficients and multiple 
          group-level predictors, the function will provide the correct dimension 
          in the error message if the dimension of beta0 is incorrectly specified, 
          and the user can respecify beta0 with this information and recall the 
          function.}
  
    \item{B0}{The prior covariance matrix of \eqn{\beta}{beta}.  This should be a 
          positive definite matrix with dimension equal to the number of betas in 
          the reduced form of GLMM-AR(p).}
  
    \item{d0, D0}{The degree of freedom and scale matrix of the Wishart prior
                  on \eqn{b_i} which are the subject-level residuals. D0 should 
                  not be too defuse, otherwise it may takes a long time for 
                  the chain to converge. Recommended values is a*I, where a is 
                  between 1.67 and 10. No default is provided.}

    \item{e0, E0}{The degree of freedom and scale matrix of the Wishart prior
                  on \eqn{c_t} which are the time-level residuals. E0 should 
                  not be too defuse, otherwise it may takes a long time for 
                  the chain to converge. Recommended values is a*I, where a is 
                  between 1.67 and 10. No default is provided.}

    \item{nlag}{A scalar of the lag order p, which should be an integer equal to 
                or greater than 1. At this stage, the function does not support the 
                GLMM-AR(0) moel, which can be estimated by using \code{BUGS} or \code{JAGS}.}

    \item{pgm}{Should the parameter expantion method be applied? It takes two values: 
               TRUE or FALSE, and the default is TRUE.}

    \item{tracking}{The tracking interval in the simulation. Every "tracking" 
              iterations, the function will return the information about how many 
              iterations in total have been done. The default is 0 and no tracking
              information will be given during the simulation.}

    \item{init.phi}{A vector of starting values of the autoregressive coefficients.}
 
    \item{init.beta}{A vector of starting values of \eqn{\beta}{beta}. The default is 0.}

    \item{init.b}{A vector of starting values of \eqn{b_i}. The default is 0.}

    \item{init.c}{A vector of starting values of \eqn{c_t}. The default is 0.}

    \item{monitor}{A string contains the names of parameters whose MCMC output are
              to be returned. The string should be a subset of
               ("rho", "beta", "bi", "ct", "D", "E", "ystar", "u") which is also the default.
               the meanings of the parameter names are: "rho" -- \eqn{\rho}{rho}, the 
               autoregressive coefficients; "beta" -- \eqn{\beta}{\beta}, the fixed-effect 
               coefficients; "bi" -- \eqn{b_i}{b_i}, the subject-level residuals; "ct" --
               \eqn{c_t}{c_t}, the time-level errors residuals; "D" -- the covariance matrix 
               of \eqn{b_i}{b_i}; "E" -- the covariance matrix of \eqn{c_t}{c_t}; "ystar" -- 
               the augmented data; and "u" -- the auxiliary variable. }
   
    \item{marginal.likelihood}{Should the marginal likelihood be calculated. The default is FALSE.}
    
    \item{reduced.mcmc}{The number of iterations to return in the reduced mcmc simulations.}

    \item{reduced.burnin}{The number of burn-in iterations for the sampler in the reduced
                     mcmc simulations.}

}

\value{
   A list contains the posterior sample of all monitored parameters. This object can 
   be converted into mcmc object by using the \code{coda} function \code{as.mcmc}, then be 
   summarized by functions provided by the \code{coda} package.
}

\details{
\code{GLMMARp.Binary} simulates from the posterior distribution of a probit
  multilevel model with a pth-order autoregressive error process by using 
  the Cholesky-decomposition-and-auxiliary-variable approach. Please consult
  the coda documentation for a comprehensive list of functions that can be
  used to analyze the posterior sample.

  The model has a structural form and a reduced form. The user fits the data 
  in the function according to the structural form, and the function will estimate 
  the reduced form. Therefore, the posteriors of the fixed-effect coefficients 
  include the fixed-effect coefficients at both the individual and group levels. 
  The function only returns the group-level residuals, and the random-effect 
  coefficients at the individual level have to be recovered by using the function 
  \code{random.time} and \code{random.unit}. 

  The structural form of the GLMM-AR(p) model is as following: 
  \deqn{y_{i,t_i}=\mathbf{I}(z_{i,t_i}>0)}{y_{i,t_i}=I(z_{i,t_i}>0)}

  \deqn{z_{i,t_i}=\mathbf{x}_{1i,t_i}'\mathbf{\beta}_1+\mathbf{w}_{i,t_i}'\mathbf{\beta}_{2i}+
   \mathbf{s}_{i,t_i}'\mathbf{\beta}_{3t_i}+\xi_{i,t_i}}{z_{i,t_i}=x_{1i,t_i}'
   \beta_1+w_{i,t_i}'\beta_{2i}+
   s_{i,t_i}'\beta_{3t_i}+\xi_{i,t_i}}

   \deqn{\mathbf{\beta}_{2i}=\mathbf{A}_{i}\mathbf{\beta}_2+\mathbf{b}_{i}}{\beta_{2i}
    =A_{i}\beta_2+b_{i}}
   
  \deqn{\mathbf{\beta}_{3t_i}=\mathbf{F}_{t_i} \mathbf{\beta}_3+ \mathbf{c}_{t_i}}{\beta_{3t_i}
    =F_{t_i}\beta_3+c_{t_i}}


  \deqn{\xi_{i,t_i}=\rho_1 \xi_{i, t_i-1}+...+\rho_p \xi_{i, t_i-p}+\epsilon_{i,t_i}}{\xi_{i,t_i}
    =\rho_1 \xi_{i, t_i-1}+...+\rho_p \xi_{i, t_i-p}+\epsilon_{i,t_i}}


   and the reduced form of the model at the latent levels is

  \deqn{z_{i,t_i}=\mathbf{x}_{i,t_i}'\mathbf{\beta}+\mathbf{w}'_{i,t_i}\mathbf{b}_{i}
   +\mathbf{s}'_{i,t_i}\mathbf{c}_{t_i}+\xi_{i,t_i}}{z_{i,t_i}=x_{i,t_i}'\beta+w'_{i,t_i}b_{i}
   +s'_{i,t_i}c_{t_i}+\xi_{i,t_i}}

  \deqn{\xi_{i,t_i}=\rho_1 \xi_{i, t_i-1}+...+\rho_p \xi_{i, t_i-p}+\epsilon_{i,t_i}}{\xi_{i,t_i}=
        \rho_1 \xi_{i, t_i-1}+...+\rho_p \xi_{i, t_i-p}+\epsilon_{i,t_i}}
 }
  
\references{
   Pang, Xun. 2009. "Intertemporal and Contemporal Dependence in Binary Time-Series Cross-Sectional Data: Bayesian Hierarchical Model with AR(p) Errors and Non-nested Clustering."\url{http://xpang.wustl.edu/Dissertation.html}.

   Albert, James A. and Siddhartha Chib. 1993. "Bayesian Analysis of Binary and Polychotomous Response Data." \emph{Journal of the American Statistical Association.} 88: 669-679.

   Chib, Siddhartha. 1995. "Marginal Likelihood from the Gibbs Output." \emph{Journal of the American Statistical Association.} 90: 1313-1321. 

  Chib, Siddhartha, and Ivan Jeliazkov. 2001. "Marginal Likelihood from the Metropolis-Hastings Output." \emph{Journal of the American Statistical Association. } 96: 270-281.
}

\examples{
 \dontrun{

## Example 1: GLMM-AR(p) model only with subject-specific and time-specific random intercepts
##            and serial correlation; no group-level predictors
 
 data(StateFailure)
 require(bayesSurv)
 require(panel)
 y <- StateFailure$failure
 Unit <- StateFailure$country
 Time <- StateFailure$year
 Fixed <- cbind(StateFailure$poldemoc,  StateFailure$bnkv123, StateFailure$bnkv117, 
              StateFailure$poldurab, StateFailure$faocalry, StateFailure$pwtcgdp, 
              StateFailure$macnac,StateFailure$macnciv,  StateFailure$wdiinfmt, 
              StateFailure$ptsamnes, StateFailure$dis1, StateFailure$bnkv81, 
              StateFailure$change.demo)
  nphi <- 0.5
  nbeta <- 0
  nb <- 0
  nc <- 0

 UnitRandom <- as.matrix(rep( 1, length(y)))
 TimeRandom <- as.matrix(rep( 1, length(y)))
 UnitPred <- "NULL"
 TimePred <- "NULL"

 State0AR1 <- GLMMARp.Binary(y=y,X1=Fixed, Wi=UnitRandom, St=TimeRandom, Ai=UnitPred, Ft=TimePred,
                        Unit.Index=Unit, Time.Index=Time, timeint.add=0, unitint.add=0, m=1000,
                        burnin=500, init.phi=nphi, init.b=nb, init.c=nc, init.beta=nbeta, 
                        beta0=rep(0,13),B0=diag(10, 13),
                        D0=1.67 ,E0=6, e0=5, d0=5,  pgm=TRUE, nlag=1, tracking=5,
                        monitor=c("rho", "beta", "bi", "ct", "D", "E", "ystar", "u"), 
                        marginal.likelihood=TRUE, reduced.mcmc=500, reduced.burnin=100, thin=1)


## Example 2: GLMM-AR(p) model with a subject-specific coefficient, random intercepts
##            and serial correlation; no group-level predictors

 UnitRandom <- cbind(log(StateFailure$pwtopen))
 TimeRandom <- as.matrix(rep( 1, length(y)))
 UnitPred <- "NULL"
 TimePred <- "NULL"


StateAR1 <- GLMMARp.Binary(y=y,X1=Fixed, Wi=UnitRandom, St=TimeRandom, Ai=UnitPred, Ft=TimePred,
                        Unit.Index=Unit, Time.Index=Time, timeint.add=0, unitint.add=1, m=1000,
                        burnin=500, init.phi=nphi, init.b=nb,
                        init.c=nc, init.beta=nbeta,beta0=rep(0,14),B0=diag(10, 14),
                        D0=diag(1.67,2) ,E0=6,e0=3, d0=4,  pgm=TRUE, nlag=1, tracking=5,
                        monitor=c("rho", "beta", "bi", "ct", "D", "E", "ystar", "u"), 
                        marginal.likelihood=TRUE, reduced.mcmc=500, reduced.burnin=100,thin=1)

## Example 3: GLMM-AR(p) model with a subject-specific coefficient, random intercepts
##            and serial correlation; only subject-level predictors.

 TimeRandom <- as.matrix(rep( 1, length(y)))
 UnitPred <- cbind(StateFailure$macnac, StateFailure$poldemoc)
 TimePred <- "NULL"
 nphi <- c(0.2, 0.1)

 StateAR2 <- GLMMARp.Binary(y=y,X1=Fixed, Wi=UnitRandom, St=TimeRandom, Ai=UnitPred, Ft=TimePred,
                        Unit.Index=Unit, Time.Index=Time, timeint.add=0, unitint.add=1, m=1000,
                        burnin=500, init.phi=nphi, init.b=nb,
                        init.c=nc, init.beta=nbeta,beta0=rep(0,15),B0=diag(10, 15),
                        D0=diag(1.67,2) ,E0=6,e0=5, d0=6, pgm=TRUE, nlag=2, tracking=5,
                        monitor=c("rho", "beta", "bi", "ct", "D", "E", "ystar", "u"), 
                        marginal.likelihood=TRUE, reduced.mcmc=500, reduced.burnin=100, thin=1)


## Example 3: GLMM-AR(p) model with a subject-specific coefficient, random intercepts
##            and serial correlation; with both subject- and time-level predictors.

 TimePred <- cbind(StateFailure$dis1)
 nphi <- c(0.3, 0.1, 0.1)

StateAR3 <- GLMMARp.Binary(y=y,X1=Fixed, Wi=UnitRandom, St=TimeRandom, Ai=UnitPred, Ft=TimePred,
                        Unit.Index=Unit, Time.Index=Time, timeint.add=1, unitint.add=1, m=1000,
                        burnin=500, init.phi=nphi, init.b=nb,
                        init.c=nc, init.beta=nbeta,beta0=rep(0,16),B0=diag(10, 16),
                        D0=diag(1.67,2) ,E0=diag(1.67,2),e0=5, d0=6, pgm=TRUE, nlag=3, tracking=5,
                        monitor=c("rho", "beta", "bi", "ct", "D", "E", "ystar", "u"), 
                        marginal.likelihood=TRUE, reduced.mcmc=500, reduced.burnin=100, thin=1)

}
}

\keyword{models}

