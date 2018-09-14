\name{optweight}
\alias{optweight}
\alias{optweight.fit}
\alias{print.optweight}
\alias{print.optweightMSM}
\title{Estimate Balancing Weights Using Optimization}
\usage{
optweight(formula,
          data = NULL,
          tols = 0,
          estimand = "ATE",
          targets = NULL,
          s.weights = NULL,
          focal = NULL,
          verbose = FALSE,
          ...)

\method{print}{optweight}(x, ...)

\method{print}{optweightMSM}(x, ...)

optweight.fit(treat,
              covs,
              tols = 0,
              estimand = "ATE",
              targets = NULL,
              s.weights = NULL,
              focal = NULL,
              norm = "l2",
              std.binary = FALSE,
              std.cont = TRUE,
              min.w = 0,
              verbose = FALSE,
              ...)
}
\description{
Estimate balancing weights for treatments and covariates specified in \code{formula}. The degree of balance for each covariate is specified by \code{tols} and the target population can be specified with \code{targets} or \code{estimand}. \code{optweight.fit} performs the optimization (via \pkg{rosqp}; Anderson, 2018) and outputs a vector of weights, and \code{optweight} prepares input for processing and outputs an object that includes the estimated weights and other values. Users should use \code{optweight} for most purposes. See Zubizarreta (2015), Wang & Zubizarreta (2017), and Yiu & Su (2018) for details in the properties of these weights and the methods used to fit them.
}
\arguments{
  \item{formula}{
A formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side, or a list thereof. See \code{\link{glm}} for more details. Interactions and functions of covariates are allowed.
}
  \item{data}{
An optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{tols}{
A vector of balance tolerance values for each covariate, or a list thereof. The resulting weighted standardized balance statistics will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. Can also be the output of a call to \code{\link{check.tols}} for point treatments. See Details.
}
  \item{estimand}{
The desired estimand, which determines the target population. For binary treatments, can be "ATE", "ATT", "ATC", or \code{NULL}. For multinomial treatments, can be "ATE", "ATT", or \code{NULL}. For continuous treatments, can be "ATE" or \code{NULL}. The default for both is "ATE". For longitudinal treatments, only "ATE" is supported. \code{estimand} is ignored when \code{targets} is non-\code{NULL}. If both \code{estimand} and \code{targets} are \code{NULL}, no targeting will take place. See Details.
}
  \item{targets}{
A vector of target populaton mean values for each baseline covariate. The resulting weights will yield sample means within \code{tols}/2 units of the target values for each covariate. If \code{NULL}, \code{estimand} will be used to determine targets. Otherwise, \code{estimand} is ignored. If any target values are \code{NA}, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance. Can also be the output of a call to \code{\link{check.targets}}. See Details.
}
  \item{s.weights}{
A vector of sampling weights or the name of a variable in \code{data} that contains sampling weights.
}
  \item{focal}{
When multinomial treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group.
}
  \item{verbose}{
Whether information on the optimization problem solution should be printed. This information contains how many iterations it took to estimate the weights and whether the solution is optimal.
}
  \item{\dots}{
For \code{optweight}, arguments passed to \code{optweight.fit}. For \code{optweight.fit}, options that are passed to \code{\link[rosqp]{osqpSettings}} for use in the \code{par} arguments of \code{\link[rosqp]{solve_osqp}}.
}
  \item{treat}{
A vector of treatment statuses for each unit or a list thereof. Non-numeric (i.e., factor or character) values are allowed.
}
  \item{covs}{
A data frame or matrix of covariates to be balanced. All variables must be numeric. Does not have to be full rank.
}

  \item{x}{
An \code{optweight} or \code{optweightMSM} object; the output of a call to \code{optweight()}.
}
  \item{norm}{
A string containing the name of the norm to minimize. The options are \code{"l1"} for the L1 norm, \code{"l2"} for the L2 norm (the default), and \code{"linf"} for the L\eqn{\infty} norm. The L1 norm minimizes the average absolute distance between each weight and the mean of the weights. The L2 norm minimizes the variance of the weights. The L\eqn{\infty} norm minimizes the largest weight.
}
  \item{std.binary, std.cont}{
\code{logical}; whether the tolerances are in standardized mean units (\code{TRUE}) or raw units (\code{FALSE}) for binary variables and contoinous variables, respectively. The default is \code{FALSE} for \code{std.binary} because raw proportion differences make more sense than standardized mean difference for binary variables. These arguments are analogous to the \code{binary} and \code{continuous} arguments in \code{\link[cobalt]{bal.tab}} in \pkg{cobalt}.
}
  \item{min.w}{
A single \code{numeric} value between 0 and 1 for the smallest allowable weight. Some analyses require nonzero weights for all units, so a nonzero minimum (e.g., \code{1E-8} can be specified. Doing so will likely (slightly) increase the variance of the resulting weights depending on the magntiude of the minimum.
}

}
\value{
  If only one time point is specified, an \code{optweight} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat}{The values of the treatment variable.}
  \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{s.weights}{The provided sampling weights.}
  \item{estimand}{The estimand requested.}
  \item{focal}{The focal variable if the ATT was requested with a multinomial treatment.}
  \item{call}{The function call.}
  \item{tols}{The tolerance values for each covariate.}
  \item{duals}{A data.frame containing the dual variables for each covariate. See Zubizarreta (2015) for interpretation of these values.}
  \item{info}{The \code{info} component of the output of \code{\link[rosqp]{solve_osqp}}, which contains information on the performance of the optimization at termination.}

  Otherwise, if multiple time points are specified, an \code{optmatchMSM} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat.list}{A list of the values of the treatment variables at each time point.}
  \item{covs.list}{A list of the covariates at each time point used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{s.weights}{The provided sampling weights.}
  \item{estimand}{The estimand requested; "ATE" for longtiduinal treatments.}
  \item{call}{The function call.}
  \item{tols}{A list of tolerance values for each covariate at each time point.}
  \item{duals}{A list of data.frames containing the dual variables for each covariate at each time point. See Zubizarreta (2015) for interpretation of these values.}
  \item{info}{The \code{info} component of the output of \code{\link[rosqp]{solve_osqp}}, which contains information on the performance of the optimization at termination.}

}
\details{
The optimization is fit using \code{\link[rosqp]{solve_osqp}} in the \pkg{rosqp} package, which provides a straightforward interface to specifying the constraints and objective function for quadratic optimization problems and uses a fast and flexible solving algorithm.

For binary and multinomial treatments, weights are estimated so that the standardized mean differences of the covariates are within the given tolerance thresholds (unless \code{std.binary} or \code{std.cont} are \code{FALSE}, in which case unstandardized mean differences are considered for binary and continuous variables, respectively). For a covariate \eqn{x} with specified tolerance \eqn{\delta}, the weighted means of each each group will be within \eqn{\delta} of each other. Additionally, when the ATE is specified as the estimand or a target population is specified, the weighted means of each group will each be within \eqn{\delta/2} of the target means; this ensures generalizability to the same population from which the original sample was drawn. If standardized tolerance values are requested, the standardization factor corresponds to the estimand requested: when the ATE is requested or a target population specified, the standardization factor is the square root of the average variance for that covariate across treatment groups, and when the ATT or ATC are requested, the standardization factor is the standard deviation of the covariate in the focal group. The standardization factor is always unweighted.

For continuous treatments, weights are estimated so that the weighted correlation between the treatment and each covariate is within the specified tolerance threshold. If the ATE is requested or a target population is specified, the means of the weighted covariates and treatment are restricted to be equal to those of the target population to ensure generalizability to the desired target population. The weighted correlation is computed as the weighted covariance divided by the product of the \emph{unweighted} standard deviations. The means used to center the variables in computing the covariance are those specified in the target population.

For longitudinal treatments, only "wide" data sets, where each row corresponds to a unit's entire variable history, are supported. You can use \code{\link{reshape}} or other functions to transform your data into this format; see example in the documentation for \code{weightitMSM} in the \pkg{WeightIt} package.

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

Wang, Y., & Zubizarreta, J. R. (2017). Approximate Balancing Weights: Characterizations from a Shrinkage Estimation Perspective. ArXiv:1705.00998 [Math, Stat]. Retrieved from \url{http://arxiv.org/abs/1705.00998}

Yiu, S., & Su, L. (2018). Covariate association eliminating weights: a unified weighting framework for causal effect estimation. Biometrika. \doi{10.1093/biomet/asy015}

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
}

\seealso{
\url{https://osqp.org/docs/index.html} for more information on \pkg{rosqp}, the underlying solver, and the options for \code{\link[rosqp]{solve_osqp}}.

\code{\link[rosqp]{osqpSettings}} for details on options for \code{solve_osqp}.
}

\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(ow1 <- optweight(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = c(.01, .02, .03, .04, .05),
                estimand = "ATE"))
bal.tab(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                weights = ow1$weights)

#Exactly alancing covariates with respect to race (multinomial)
(ow2 <- optweight(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                exact = TRUE, estimand = "ATT", focal = "black"))
bal.tab(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                weights = ow2$weights, focal = "black")

#Balancng covariates with longitudinal treatments
library("twang")
data("iptwExWide")

##Weighting more recent covariates more strictly
(ow3 <- optweight(list(tx1 ~ use0 + gender + age,
                      tx2 ~ tx1 + use1 + use0  + gender +
                            age,
                      tx3 ~ tx2 + use2 + tx1 + use1 +
                            use0 + gender + age),
                 data = iptwExWide,
                 tols = list(c(.001, .001, .001),
                             c(.001, .001, .01, .01, .01),
                             c(.001, .001, .01, .01,
                               .1, .1, .1))))
bal.tab(ow3)

#Balancing covariates between treatment groups (binary)
#and requesting a specified target population
(ow4a <- optweight(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = 0,
                targets = c(26, 12, .4, .5, 1000),
                estimand = NULL))
bal.tab(ow4a, disp.means = TRUE)

#Balancing covariates between treatment groups (binary)
#and not requesting a target population
(ow4b <- optweight(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = 0,
                targets = NULL,
                estimand = NULL))
bal.tab(ow4b, disp.means = TRUE)
}
