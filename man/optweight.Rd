\name{optweight}
\alias{optweight}
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
          b.weights = NULL,
          focal = NULL,
          verbose = FALSE,
          force = FALSE,
          ...)

\method{print}{optweight}(x, ...)

\method{print}{optweightMSM}(x, ...)
}
\description{
Estimate balancing weights for treatments and covariates specified in \code{formula}. The degree of balance for each covariate is specified by \code{tols} and the target population can be specified with \code{targets} or \code{estimand}. See Zubizarreta (2015), Wang & Zubizarreta (2017), and Yiu & Su (2018) for details of the properties of the weights and the methods used to fit them.
}
\arguments{
  \item{formula}{
A formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side, or a list thereof. See \code{\link{glm}} for more details. Interactions and functions of covariates are allowed.
}
  \item{data}{
An optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{tols}{
A vector of balance tolerance values for each covariate, or a list thereof. The resulting weighted balance statistics will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. Can also be the output of a call to \code{\link{check.tols}} for point treatments. See Details.
}
  \item{estimand}{
The desired estimand, which determines the target population. For binary treatments, can be "ATE", "ATT", "ATC", or \code{NULL}. For multinomial treatments, can be "ATE", "ATT", or \code{NULL}. For continuous treatments, can be "ATE" or \code{NULL}. The default for both is "ATE". For longitudinal treatments, only "ATE" is supported. \code{estimand} is ignored when \code{targets} is non-\code{NULL}. If both \code{estimand} and \code{targets} are \code{NULL}, no targeting will take place. See Details.
}
  \item{targets}{
A vector of target populaton mean values for each baseline covariate. The resulting weights will yield sample means within \code{tols}/2 units of the target values for each covariate. If \code{NULL} or all \code{NA}, \code{estimand} will be used to determine targets. Otherwise, \code{estimand} is ignored. If any target values are \code{NA}, the corresponding variable will not be targeted and its weighted mean will be wherever the weights yield the smallest variance. Can also be the output of a call to \code{\link{check.targets}}. See Details.
}
  \item{s.weights}{
A vector of sampling weights or the name of a variable in \code{data} that contains sampling weights. Optimization occurs on the product of the sampling weights and the estimated weights.
}
  \item{b.weights}{
A vector of base weights or the name of a variable in \code{data} that contains base weights. If supplied, the desired norm of the distance between the estimated weights and the base weights is minimized. Currently only supported with \code{norm = "l2"} (the default norm; see \code{\link{optweight.fit}}).
}
  \item{focal}{
When multinomial treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group. If specified, \code{estimand} will automatically be set to \code{"ATT"}.
}
  \item{verbose}{
Whether information on the optimization problem solution should be printed. This information contains how many iterations it took to estimate the weights and whether the solution is optimal.
}
  \item{force}{
optweights are currently not valid for use with longitudinal treatments, and will produce an error message if attempted. Set to \code{TRUE} to bypass this error message.
}
  \item{\dots}{
For \code{optweight}, arguments passed to \code{\link{optweight.fit}}. Ignored otherwise.
}

  \item{x}{
An \code{optweight} or \code{optweightMSM} object; the output of a call to \code{optweight()}.
}

}
\value{
  If only one time point is specified, an \code{optweight} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat}{The values of the treatment variable.}
  \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{s.weights}{The provided sampling weights.}
  \item{b.weights}{The provided base weights.}
  \item{estimand}{The estimand requested.}
  \item{focal}{The focal variable if the ATT was requested with a multinomial treatment.}
  \item{call}{The function call.}
  \item{tols}{The tolerance values for each covariate.}
  \item{duals}{A data.frame containing the dual variables for each covariate. See Details for interpretation of these values.}
  \item{info}{The \code{info} component of the output of \code{\link[osqp]{solve_osqp}}, which contains information on the performance of the optimization at termination.}

  Otherwise, if multiple time points are specified, an \code{optmatchMSM} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat.list}{A list of the values of the treatment variables at each time point.}
  \item{covs.list}{A list of the covariates at each time point used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{s.weights}{The provided sampling weights.}
  \item{b.weights}{The provided base weights.}
  \item{call}{The function call.}
  \item{tols}{A list of tolerance values for each covariate at each time point.}
  \item{duals}{A list of data.frames containing the dual variables for each covariate at each time point. See Details for interpretation of these values.}
  \item{info}{The \code{info} component of the output of \code{\link[osqp]{solve_osqp}}, which contains information on the performance of the optimization at termination.}

}
\details{
The optimization is performed by the lower-level function \code{\link{optweight.fit}} using \code{\link[osqp]{solve_osqp}} in the \pkg{osqp} package, which provides a straightforward interface to specifying the constraints and objective function for quadratic optimization problems and uses a fast and flexible solving algorithm.

For binary and multinomial treatments, weights are estimated so that the weighted mean differences of the covariates are within the given tolerance thresholds (unless \code{std.binary} or \code{std.cont} are \code{TRUE}, in which case standardized mean differences are considered for binary and continuous variables, respectively). For a covariate \eqn{x} with specified tolerance \eqn{\delta}, the weighted means of each each group will be within \eqn{\delta} of each other. Additionally, when the ATE is specified as the estimand or a target population is specified, the weighted means of each group will each be within \eqn{\delta/2} of the target means; this ensures generalizability to the same population from which the original sample was drawn.

If standardized tolerance values are requested, the standardization factor corresponds to the estimand requested: when the ATE is requested or a target population specified, the standardization factor is the square root of the average variance for that covariate across treatment groups, and when the ATT or ATC are requested, the standardization factor is the standard deviation of the covariate in the focal group. The standardization factor is always unweighted.

For continuous treatments, weights are estimated so that the weighted correlation between the treatment and each covariate is within the specified tolerance threshold. If the ATE is requested or a target population is specified, the means of the weighted covariates and treatment are restricted to be equal to those of the target population to ensure generalizability to the desired target population. The weighted correlation is computed as the weighted covariance divided by the product of the \emph{unweighted} standard deviations. The means used to center the variables in computing the covariance are those specified in the target population.

For longitudinal treatments, only "wide" data sets, where each row corresponds to a unit's entire variable history, are supported. You can use \code{\link{reshape}} or other functions to transform your data into this format; see example in the documentation for \code{weightitMSM} in the \pkg{WeightIt} package. Currently, longtiduinal treatments are not recommended as optweight's use with them has not been validated.

\subsection{Dual Variables}{
Two types of constriants may be associated with each covariate: target constraints and balance constraints. Target constraints require the mean of the covariate to be at (or near) a specific target value in each treatment group (or for the whole group when treatment is continuous). Balance constraints require the means of the covariate in pairs of treatments to be near each other. For binary and multinomial treatments, balance constraints are redundant if target constraints are provided for a variable. For continuous variables, balance constraints refer to the correlation between treatment and the covariate and are not redundant with target constraints. In the \code{duals} component of the output, each covariate has a dual variable for each nonredundant constraint placed on it.

The dual variable for each constraint is the instantaneous rate of change of the objective function at the optimum due to a change in the constraint. Because this relationship is not linear, large changes in the constraint will not exactly map onto corresponding changes in the objective function at the optimum, but will be close for small changes in the constraint. For example, for a covariate with a balance constraint of .01 and a corresponding dual variable of .4, increasing (i.e., relaxing) the constraint to .025 will decrease the value of the objective function at the optimum by approximately (.025 - .01) * .4 = .006. When the L2 norm is used, this change corresponds to a change in the variance of the weights, which directly affects the effective sample size (though the magnitude of this effect depends on the original value of the effective sample size).

For factor variables, \code{optweight} takes the sum of the absolute dual variables for the constraints for all levels and reports it as the the single dual variable for the variable itself. This summed dual variable works the same way as dual variables for continuous variables do.

}

\subsection{Solving Convergence Failure}{

Sometimes the optimization will fail to converge at a solution. There are a variety of reasons why this might happen, which include that the constraints are nearly impossible to satisfy or that the optimization surface is relatively flat. It can be hard to know the exact cause or how to solve it, but this section offers some solutions one might try.

Rarely is the problem too few iterations, though this is possible. Most problems can be solved in the default 200,000 iterations, but sometimes it can help to increase this number with the \code{max_iter} argument. Usually, though, this just ends up taking more time without a solution found.

If the problem is that the constraints are too tight, it can be helpful to loosen the constraints. Sometimes examining the dual variables of a solution that has failed to converge can reveal which constraints are causing the problem.

Sometimes a suboptimal solution is possible; such a solution does not satisfy the constraints exactly but will come pretty close. To allow these solutions, the arguments \code{eps_abs} and \code{eps_rel} can be increased from 1E-8 to larger values. These should be adjusted together since they both must be satisfied for convergence to occur; this can be done easily using the shortcut argument \code{eps}, which changes both \code{eps_abs} and \code{eps_rel} to the set value.

With continuous treatments, solutions that failed to converge may still be useable. Make sure to assess balance and examine the weights even after a optimal solution is not found, because the solution that is found may be good enough.
}
}

\author{
Noah Greifer
}

\references{
Anderson, E. (2018). osqp: Quadratic Programming Solver using the 'OSQP' Library. R package version 0.1.0. \url{https://CRAN.R-project.org/package=osqp}

Wang, Y., & Zubizarreta, J. R. (2017). Approximate Balancing Weights: Characterizations from a Shrinkage Estimation Perspective. ArXiv:1705.00998 [Math, Stat]. Retrieved from \url{http://arxiv.org/abs/1705.00998}

Yiu, S., & Su, L. (2018). Covariate association eliminating weights: a unified weighting framework for causal effect estimation. Biometrika. \doi{10.1093/biomet/asy015}

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
}

\seealso{
\url{https://osqp.org/docs/index.html} for more information on \pkg{osqp}, the underlying solver, and the options for \code{\link[osqp]{solve_osqp}}.

\code{\link[osqp]{osqpSettings}} for details on options for \code{solve_osqp}.

\code{\link{optweight.fit}}, the lower-level function that performs the fitting.

The package \pkg{sbw}, which was the inspiration for this package and provides additional functionality for binary treatments.
}

\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(ow1 <- optweight(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = c(.01, .02, .03, .04, .05),
                estimand = "ATE"))
bal.tab(ow1)

#Exactly alancing covariates with respect to race (multinomial)
(ow2 <- optweight(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = 0, estimand = "ATT", focal = "black"))
bal.tab(ow2)

# #Balancing covariates with longitudinal treatments
# #NOT VALID; DO NOT DO THIS.
# library("twang")
# data("iptwExWide")
#
# ##Weighting more recent covariates more strictly
# (ow3 <- optweight(list(tx1 ~ use0 + gender + age,
#                        tx2 ~ tx1 + use1 + use0  + gender +
#                          age,
#                        tx3 ~ tx2 + use2 + tx1 + use1 +
#                          use0 + gender + age),
#                   data = iptwExWide,
#                   tols = list(c(.001, .001, .001),
#                               c(.001, .001, .01, .01, .01),
#                               c(.001, .001, .01, .01,
#                                 .1, .1, .1))))
# bal.tab(ow3)

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
