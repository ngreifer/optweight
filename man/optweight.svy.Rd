\name{optweight.svy}
\alias{optweight.svy}
\alias{print.optweight.svy}
\title{Estimate Targeting Weights Using Optimization}
\usage{
optweight.svy(formula,
              data = NULL,
              tols = 0,
              targets = NULL,
              s.weights = NULL,
              verbose = FALSE,
              ...)

\method{print}{optweight.svy}(x, ...)

}
\description{
Estimate targeting weights for covariates specified in \code{formula}. The target means are specified with \code{targets} and the maximum distance between each weighted covariate mean and the corresponding target mean is specified by \code{tols}. See Zubizarreta (2015) for details of the properties of the weights and the methods used to fit them.
}
\arguments{
  \item{formula}{
A formula with nothing on the left hand side and the covariates to be targeted on the right hand side. See \code{\link{glm}} for more details. Interactions and functions of covariates are allowed.
}
  \item{data}{
An optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{tols}{
A vector of target balance tolerance values for each covariate. The resulting weighted covariate means will be no further away from the targets than the specified values. If only one value is supplied, it will be applied to all covariates. Can also be the output of a call to \code{\link{check.tols}}. See Details.
}
  \item{targets}{
A vector of target populaton mean values for each covariate. The resulting weights will yield sample means within \code{tols} units of the target values for each covariate. If any target values are \code{NA}, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance. To ensure the weighted mean for a covairate is equal to its unweighted mean (i.e., so that its orginal mean is its target mean), its original mean must be supplied as a target.
}
  \item{s.weights}{
A vector of sampling weights or the name of a variable in \code{data} that contains sampling weights. Optimization occurs on the product of the sampling weights and the estimated weights.
}
  \item{verbose}{
Whether information on the optimization problem solution should be printed. This information contains how many iterations it took to estimate the weights and whether the solution is optimal.
}
  \item{\dots}{
For \code{optweight.svy}, arguments passed to \code{\link{optweight.svy.fit}}. Ignored otherwise.
}

  \item{x}{
An \code{optweight.svy} object; the output of a call to \code{optweight.svy()}.
}

}
\value{
  An \code{optweight.svy} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{s.weights}{The provided sampling weights.}
  \item{call}{The function call.}
  \item{tols}{The tolerance values for each covariate.}
  \item{duals}{A data.frame containing the dual variables for each covariate. See Details for interpretation of these values.}
  \item{info}{The \code{info} component of the output of \code{\link[osqp]{solve_osqp}}, which contains information on the performance of the optimization at termination.}

}
\details{
The optimization is performed by the lower-level function \code{\link{optweight.svy.fit}} using \code{\link[osqp]{solve_osqp}} in the \pkg{osqp} package, which provides a straightforward interface to specifying the constraints and objective function for quadratic optimization problems and uses a fast and flexible solving algorithm.

Weights are estimated so that the standardized differences between the weighted covariate means and the corresponding targets are within the given tolerance thresholds (unless \code{std.binary} or \code{std.cont} are \code{FALSE}, in which case unstandardized mean differences are considered for binary and continuous variables, respectively). For a covariate \eqn{x} with specified tolerance \eqn{\delta}, the weighted mean will be within \eqn{\delta} of the target. If standardized tolerance values are requested, the standardization factor is the standard deviation of the covariate in the whole sample. The standardization factor is always unweighted.

See the \code{\link{optweight}} help page for information on interpreting dual variables and solving convergence failure.
}
}

\author{
Noah Greifer
}

\references{
Anderson, E. (2018). osqp: Quadratic Programming Solver using the 'OSQP' Library. R package version 0.1.0. \url{https://CRAN.R-project.org/package=osqp}

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
}

\seealso{
\url{https://osqp.org/docs/index.html} for more information on \pkg{osqp}, the underlying solver, and the options for \code{\link[osqp]{solve_osqp}}.

\code{\link[osqp]{osqpSettings}} for details on options for \code{solve_osqp}.

\code{\link{optweight.svy.fit}}, the lower-level function that performs the fitting.

\code{\link{optweight}} for estimating weights that balance treatment groups.
}

\examples{
library("cobalt")
data("lalonde", package = "cobalt")

cov.formula <- ~ age + educ + race + married +
                      nodegree

targets <- check.targets(cov.formula, data = lalonde,
                        targets = c(23, 9, .3, .3, .4,
                                    .2, .5))

tols <- check.tols(cov.formula, data = lalonde,
                   tols = 0)

ows <- optweight.svy(cov.formula,
                     data = lalonde,
                     tols = tols,
                     targets = targets)
ows

covs <- splitfactor(lalonde[c("age", "educ", "race",
                  "married", "nodegree")],
                  drop.first = FALSE)
#Unweighted means
apply(covs, 2, mean)

#Weighted means; same as targets
apply(covs, 2, weighted.mean, w = ows$weights)
}
