\name{optweight.fit}
\alias{optweight.fit}
\title{Fitting Function for Optweight}
\usage{
optweight.fit(treat.list,
              covs.list,
              tols,
              estimand = "ATE",
              targets = NULL,
              s.weights = NULL,
              focal = NULL,
              norm = "l2",
              std.binary = FALSE,
              std.cont = TRUE,
              min.w = 1E-8,
              verbose = FALSE,
              force = FALSE,
              ...)
}
\description{
\code{optweight.fit} performs the optimization (via \pkg{rosqp}; Anderson, 2018) for \code{optweight} and should, in most coses, not be used directly. No processing of inputs is performed, so they must be given exactly as described below.
}
\arguments{
  \item{treat.list}{
A list containing one vector of treatment statuses for each time point. Non-numeric (i.e., factor or character) vectors are allowed.
}
  \item{covs.list}{
A list containing one matrix of covariates to be balanced for each time point. All matrices must be numeric but do not have to be full rank.
}
  \item{tols}{
A list containing one vector of balance tolerance values for each time point.
}
  \item{estimand}{
The desired estimand, which determines the target population. For binary treatments, can be "ATE", "ATT", "ATC", or \code{NULL}. For multinomial treatments, can be "ATE", "ATT", or \code{NULL}. For continuous treatments, can be "ATE" or \code{NULL}. The default for both is "ATE". For longitudinal treatments, only "ATE" is supported. \code{estimand} is ignored when \code{targets} is non-\code{NULL}. If both \code{estimand} and \code{targets} are \code{NULL}, no targeting will take place. See Details.
}
  \item{targets}{
A vector of target populaton mean values for each baseline covariate. The resulting weights will yield sample means within \code{tols}/2 units of the target values for each covariate. If \code{NULL} or all \code{NA}, \code{estimand} will be used to determine targets. Otherwise, \code{estimand} is ignored. If any target values are \code{NA}, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance.
}
  \item{s.weights}{
A vector of sampling weights. Optimization occurs on the product of the sampling weights and the estimated weights.
}

  \item{focal}{
When multinomial treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group.
}
  \item{norm}{
A string containing the name of the norm corresponding to the objective function to minimize. The options are \code{"l1"} for the L1 norm, \code{"l2"} for the L2 norm (the default), and \code{"linf"} for the L\eqn{\infty} norm. The L1 norm minimizes the average absolute distance between each weight and the mean of the weights; the L2 norm minimizes the variance of the weights; the L\eqn{\infty} norm minimizes the largest weight. The L2 norm has a direct correspondence with the effective sample size, making it ideal if this is your criterion of interest.
}
  \item{std.binary, std.cont}{
\code{logical}; whether the tolerances are in standardized mean units (\code{TRUE}) or raw units (\code{FALSE}) for binary variables and continuous variables, respectively. The default is \code{FALSE} for \code{std.binary} because raw proportion differences make more sense than standardized mean difference for binary variables. These arguments are analogous to the \code{binary} and \code{continuous} arguments in \code{\link[cobalt]{bal.tab}} in \pkg{cobalt}.
}
  \item{min.w}{
A single \code{numeric} value between 0 and 1 for the smallest allowable weight. Some analyses require nonzero weights for all units, so a small, nonzero minimum may be desirable. Doing so will likely (slightly) increase the variance of the resulting weights depending on the magntiude of the minimum. The default is 1e-8, which does not materially change the properties of the weights from a minimum of 0 but prevents warnings in some packages that use weights to estimate treatment effects.
}
  \item{verbose}{
Whether information on the optimization problem solution should be printed. This information contains how many iterations it took to estimate the weights and whether the solution is optimal.
}
  \item{force}{
optweights are currently not valid for use with longitudinal treatments, and will produce an error message if attempted. Set to \code{TRUE} to bypass this error message.
}
  \item{\dots}{
Options that are passed to \code{\link[rosqp]{osqpSettings}} for use in the \code{par} arguments of \code{\link[rosqp]{solve_osqp}}. See Details for defaults.
}

}
\value{
  An \code{optweight.fit} object with the following elements:
  \item{w}{The estimated weights, one for each unit.}
  \item{duals}{A data.frame containing the dual variables for each covariate, or a list thereof. See Zubizarreta (2015) for interpretation of these values.}
  \item{info}{The \code{info} component of the output of \code{\link[rosqp]{solve_osqp}}, which contains information on the performance of the optimization at termination.}

}
\details{
\code{optweight.fit} transforms the inputs into the required inputs for \code{\link[rosqp]{solve_osqp}}, which are (sparse) matrices and vectors, and then supplies the outputs (the weights, duals variables, and convergence information) back to \code{optweight}. No processing of inputs is performed, as this is normally handled by \code{optweight}.

The default values for some of the parameters sent to \code{solve_osqp} are not the same as those in \code{osqpSettings}. The following are the differences: \code{max_iter} is set to 20000 and \code{eps_abs} and \code{eps_rel} are set to 1E-8 (i.e., 10^-8). All other values are the same.

Note that optweights with longitudinal treatments are not valid and should not be used until further research is done. \code{optweight.fit} will not provide any error messages, but
}

\author{
Noah Greifer
}

\references{
Anderson, E. (2018). rosqp: Quadratic Programming Solver using the 'OSQP' Library. R package version 0.1.0. \url{https://CRAN.R-project.org/package=rosqp}

Wang, Y., & Zubizarreta, J. R. (2017). Approximate Balancing Weights: Characterizations from a Shrinkage Estimation Perspective. ArXiv:1705.00998 [Math, Stat]. Retrieved from \url{http://arxiv.org/abs/1705.00998}

Yiu, S., & Su, L. (2018). Covariate association eliminating weights: a unified weighting framework for causal effect estimation. Biometrika. \doi{10.1093/biomet/asy015}

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
}

\seealso{
\code{\link{optweight}} which you should use for estimating the balancing weights, unless you know better.

\url{https://osqp.org/docs/index.html} for more information on \pkg{rosqp}, the underlying solver, and the options for \code{\link[rosqp]{solve_osqp}}.

\code{\link[rosqp]{osqpSettings}} for details on options for \code{solve_osqp}.
}

\examples{
library("cobalt")
data("lalonde", package = "cobalt")

treat.list <- list(lalonde$treat)
covs.list <- list(splitfactor(lalonde[2:8], drop.first = "if2"))
tols.list <- list(rep(.01, ncol(covs.list[[1]])))

ow.fit <- optweight.fit(treat.list,
                        covs.list,
                        tols = tols.list,
                        estimand = "ATE",
                        norm = "l2")
}
