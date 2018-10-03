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
A vector of target populaton mean values for each covariate. The resulting weights will yield sample means within \code{tols} units of the target values for each covariate. If any target values are \code{NA}, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance. To ensure the weighted mean for a covairate is equal to its unweighted mean (i.e., so that its orginal mean is its target mean), its original mean must be supplied as a target. See Details.
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
  \item{info}{The \code{info} component of the output of \code{\link[rosqp]{solve_osqp}}, which contains information on the performance of the optimization at termination.}

}
\details{
The optimization is performed by the lower-level function \code{\link{optweight.svy.fit}} using \code{\link[rosqp]{solve_osqp}} in the \pkg{rosqp} package, which provides a straightforward interface to specifying the constraints and objective function for quadratic optimization problems and uses a fast and flexible solving algorithm.

Weights are estimated so that the standardized differences between the weighted covariate means and the corresponding targets are within the given tolerance thresholds (unless \code{std.binary} or \code{std.cont} are \code{FALSE}, in which case unstandardized mean differences are considered for binary and continuous variables, respectively). For a covariate \eqn{x} with specified tolerance \eqn{\delta}, the weighted mean will be within \eqn{\delta} of teh target. If standardized tolerance values are requested, the standardization factor is the standard deviation of the covariate in the whole sample. The standardization factor is always unweighted.

\subsection{Dual Variables}{
Each covariate is associated with a target constraint. In the \code{duals} component of the output, each covariate has a dual variable for the constraint placed on it. The dual variable for each constraint is the instantaneous rate of change of the objective function at the optimum due to a change in the constraint. Because this relationship is not linear, large changes in the constraint will not exactly map onto corresponding changes in the objective function at the optimum, but will be close for small changes in the constraint. For example, for a covariate with a target balance constraint of .01 and a corresponding dual variable of .4, increasing (i.e., relaxing) the constraint to .025 will decrease the value of the objective function at the optimum by approximately (.025 - .01) * .4 = .006. When the L2 norm is used, this change corresponds to a change in the variance of the weights, which directly affects the effective sample size (though the magntidue of this effect depends on the original value of the effective sample size).

For factor variables, \code{optweight} takes the sum of the absolute dual variables for the constraints for all levels and reports it as the the single dual variable for the variable itself. This summed dual variable works the same way as dual variables for continuous variables do.
}

\subsection{Solving Convergence Failure}{

Sometimes the opimization will fail to converge at a solution. There are a variety of reasons why this might happen, which include that the constraints are nearly impossible to satisfy or that the optimization surface is relatively flat. It can be hard to know the exact cause or how to solve it, but this section offers some solutions one might try.

Rarely is the problem too few iterations, though this is possible. Most problems can be solved in the default 200,000 iterations, but sometimes it can help to increase this number with the \code{max_iter} argument. Usually, though, this just ends up taking more time without a solution found.

If the problem is that the constraints are too tight, it can be helpful to loosen the constraints. Sometimes examining the dual variables of a solution that has failed to converge can reveal which constraints are causing the problem.

Sometimes a suboptimal solution is possible; such a solution does not satisfy the constraints exactly but will come pretty close. To allow these solutions, the arguments \code{eps_abs} and \code{eps_rel} can be increased from 1E-9 to larger values. These should be adjusted together since they both must be satisfied for converge to occur.
}
}

\author{
Noah Greifer
}

\references{
Anderson, E. (2018). rosqp: Quadratic Programming Solver using the 'OSQP' Library. R package version 0.1.0. \url{https://CRAN.R-project.org/package=rosqp}

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
}

\seealso{
\url{https://osqp.org/docs/index.html} for more information on \pkg{rosqp}, the underlying solver, and the options for \code{\link[rosqp]{solve_osqp}}.

\code{\link[rosqp]{osqpSettings}} for details on options for \code{solve_osqp}.

\code{\link{optweight.svy.fit}}, the lower-level function that performs the fitting.
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
