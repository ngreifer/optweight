\name{optweight.svy.fit}
\alias{optweight.svy.fit}
\title{Fitting Function for Optweight for Survey Weights}
\usage{
optweight.svy.fit(covs,
                  tols = 0,
                  targets,
                  s.weights = NULL,
                  norm = "l2",
                  std.binary = FALSE,
                  std.cont = TRUE,
                  min.w = 1E-8,
                  verbose = FALSE,
                  ...)
}
\description{
\code{optweight.svy.fit} performs the optimization (via \pkg{rosqp}; Anderson, 2018) for \code{optweight.svy} and should, in most coses, not be used directly. No processing of inputs is performed, so they must be given exactly as described below.
}
\arguments{
  \item{covs}{
A matrix of covariates to be targeted. Should must be numeric but does not have to be full rank.
}
  \item{tols}{
A vector of target balance tolerance values.
}
  \item{targets}{
A vector of target populaton mean values for each covariate. The resulting weights will yield sample means within \code{tols} units of the target values for each covariate. If any target values are \code{NA}, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance. To ensure the weighted mean for a covairate is equal to its unweighted mean (i.e., so that its orginal mean is its target mean), its original mean must be supplied as a target.
}
  \item{s.weights}{
A vector of sampling weights. Optimization occurs on the product of the sampling weights and the estimated weights.
}
  \item{norm}{
A string containing the name of the norm corresponding to the objective function to minimize. The options are \code{"l1"} for the L1 norm, \code{"l2"} for the L2 norm (the default), and \code{"linf"} for the L\eqn{\infty} norm. The L1 norm minimizes the average absolute distance between each weight and the mean of the weights; the L2 norm minimizes the variance of the weights; the L\eqn{\infty} norm minimizes the largest weight. The L2 norm has a direct correspondence with the effective sample size, making it ideal if this is your criterion of interest.
}
  \item{std.binary, std.cont}{
\code{logical}; whether the tolerances are in standardized mean units (\code{TRUE}) or raw units (\code{FALSE}) for binary variables and continuous variables, respectively. The default is \code{FALSE} for \code{std.binary} because raw proportion differences make more sense than standardized mean difference for binary variables.
}
  \item{min.w}{
A single \code{numeric} value between 0 and 1 for the smallest allowable weight. Some analyses require nonzero weights for all units, so a small, nonzero minimum may be desirable. Doing so will likely (slightly) increase the variance of the resulting weights depending on the magntiude of the minimum. The default is 1e-8, which does not materially change the properties of the weights from a minimum of 0 but prevents warnings in some packages that use weights to estimate treatment effects.
}
  \item{verbose}{
Whether information on the optimization problem solution should be printed. This information contains how many iterations it took to estimate the weights and whether the solution is optimal.
}
  \item{\dots}{
Options that are passed to \code{\link[rosqp]{osqpSettings}} for use in the \code{par} arguments of \code{\link[rosqp]{solve_osqp}}.
}

}
\value{
  An \code{optweight.svy.fit} object with the following elements:
  \item{w}{The estimated weights, one for each unit.}
  \item{duals}{A data.frame containing the dual variables for each covariate. See Zubizarreta (2015) for interpretation of these values.}
  \item{info}{The \code{info} component of the output of \code{\link[rosqp]{solve_osqp}}, which contains information on the performance of the optimization at termination.}

}
\details{
\code{optweight.svy.fit} transforms the inputs into the required inputs for \code{\link[rosqp]{solve_osqp}}, which are (sparse) matrices and vectors, and then supplies the outputs (the weights, duals variables, and convergence information) back to \code{optweight.svy}. No processing of inputs is performed, as this is normally handled by \code{optweight.svy}.
}

\author{
Noah Greifer
}

\references{
Anderson, E. (2018). rosqp: Quadratic Programming Solver using the 'OSQP' Library. R package version 0.1.0. \url{https://CRAN.R-project.org/package=rosqp}

Wang, Y., & Zubizarreta, J. R. (2017). Approximate Balancing Weights: Characterizations from a Shrinkage Estimation Perspective. ArXiv:1705.00998 [Math, Stat]. Retrieved from \url{http://arxiv.org/abs/1705.00998}

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
}

\seealso{
\code{\link{optweight.svy}} which you should use for estimating the balancing weights, unless you know better.

\url{https://osqp.org/docs/index.html} for more information on \pkg{rosqp}, the underlying solver, and the options for \code{\link[rosqp]{solve_osqp}}.

\code{\link[rosqp]{osqpSettings}} for details on options for \code{solve_osqp}.
}

\examples{
library("cobalt")
data("lalonde", package = "cobalt")

covs <- splitfactor(lalonde[c("age", "educ", "race",
                  "married", "nodegree")],
                  drop.first = FALSE)

targets <- c(23, 9, .3, .3, .4, .2, .5)

tols <- rep(0, 7)

ows.fit <- optweight.svy.fit(covs,
                             tols = tols,
                             targets = targets,
                             norm = "l2")

#Unweighted means
apply(covs, 2, mean)

#Weighted means; same as targets
apply(covs, 2, weighted.mean, w = ows.fit$w)
}
