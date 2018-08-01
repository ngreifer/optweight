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
          s.weights = NULL,
          focal = NULL,
          verbose = FALSE,
          ...)

\method{print}{optweight}(x, ...)

\method{print}{optweightMSM}(x, ...)

optweight.fit(treat,
              covs,
              tols = .001,
              estimand = "ATE",
              s.weights = NULL,
              focal = NULL,
              std.binary = FALSE,
              exact = FALSE,
              verbose = FALSE,
              ...)
}
\description{
Estimate balancing weights for treatments and covariates specified in \code{formula}. The degree of balance for each covariate is specified by \code{tols}. \code{optweight.fit} performs the optimization (via \pkg{rosqp}; Anderson, 2018) and outputs a vector of weights, and \code{optweight} prepares input for processing and outputs an object that includes the estimated weights and other values. Users should use \code{optweight} for most purposes. See Zubizarreta (2015), Wang & Zubizarreta (2017), and Yiu & Su (2018) for details in the properties of these weights and the methods used to fit them.

}
\arguments{
  \item{formula}{
a formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side, or a list thereof. See \code{\link{glm}} for more details. Interactions and functions of covariates are allowed.
}
  \item{data}{
an optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{tols}{
a vector of balance tolerance values in standardized mean difference units for each covariate, or a list thereof. The resulting weighted standardized mean differences will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. For binary variables, the values correspond to raw difference in proportion rather than standardized mean differences. Can also but the output of a call to \code{\link{check.tols}} for point treatments. For \code{optweight}, if \code{0}, exact balancing will be performed with a slight speedup relative to setting 0 for each covariate individually.
}
  \item{estimand}{
the desired estimand. For binary treatments, can be "ATE", "ATT", "ATC". For multinomial treatments, can be "ATE" or "ATT". The default for both is "ATE". For longitudinal treatments, only "ATE" is supported.
}
  \item{s.weights}{
a vector of sampling weights or the name of a variable in \code{data} that contains sampling weights.
}
  \item{focal}{
when multinomial treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group.
}
  \item{exact}{
\code{logical}; whether exact rather than approximate balancing should occur. If \code{TRUE}, \code{tols} is ignored. Essentially equivalent to setting \code{tols = 0}, but code might run a bit faster and less memory is allocated.
}
  \item{verbose}{
whether information on the optimization problem solution should be printed. This information contains how many iterations it took to estimate the weights and whether the solution is optimal.
}
  \item{\dots}{
for \code{optweight} and \code{optweight.fit}, options that are passed to code{\link[rosqp]{osqpSettings}} for use in the \code{par} arguments of code{\link[rosqp]{solve_osqp}}. Of particular note are \code{adaptive_rho}, which when set to \code{TRUE} dramatically speeds up the estimation of the weights at the cost of a slightly submoptimal solution (default \code{FALSE}), and \code{max_iter}, which controls the maximum number of iterations of the solver, and can be set higher if the problem is not solved in the default number of iterations (200000).
}
  \item{treat}{
a vector of treatment statuses for each unit or a list thereof. Non-numeric (i.e., factor or character) values are allowed.
}
  \item{covs}{
a data frame or matrix of covariates to be balanced. All variables must be numeric. Does not have to be full rank.
}
  \item{std.binary}{
whether the tolerances are in standardized mean units (\code{TRUE}) or raw units (\code{FALSE}) for binary variables. The default is \code{FALSE} because raw proportion differences make more sense than standardized mean difference for binary variables. This argument is analogous to the \code{binary} argument in \code{\link[cobalt]{bal.tab}} in \pkg{cobalt}.
}
  \item{x}{
an \code{optweight} or \code{optweightMSM} object; the output of a call to \code{optweight()}.
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

  Otherwise, if multiple time points are specified, an \code{optmatchMSM} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat.list}{A list of the values of the treatment variables at each time point.}
  \item{covs.list}{A list of the covariates at each time point used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{s.weights}{The provided sampling weights.}
  \item{estimand}{The estimand requested; "ATE" for longtiduinal treatments.}
  \item{call}{The function call.}

}
\details{
The optimization is fit using \code{\link[rosqp]{solve_osqp}} in the \pkg{rosqp} package, which provides a straghtforward interface to specifying the constraints and objective function for quadratic optimization problems and uses a fast and flexible solving algorithm.

For longitudinal treatments, only "wide" data sets, where each row corresponds to a unit's entire variable history, are supported. You can use \code{\link{reshape}} or other functions to transform your data into this format; see example in the documentation for \code{weightitMSM} in the \pkg{WeightIt} package.
}

\author{
Noah Greifer
}

\references{
Eric Anderson (2018). rosqp: Quadratic Programming Solver using the 'OSQP' Library. R package version 0.1.0. https://CRAN.R-project.org/package=rosqp

Wang, Y., & Zubizarreta, J. R. (2017). Approximate Balancing Weights: Characterizations from a Shrinkage Estimation Perspective. ArXiv:1705.00998 [Math, Stat]. Retrieved from \url{http://arxiv.org/abs/1705.00998}

Yiu, S., & Su, L. (2018). Covariate association eliminating weights: a unified weighting framework for causal effect estimation. Biometrika. \doi{10.1093/biomet/asy015}

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}

\seealso{
\url{https://osqp.org/docs/index.html} for more information on \pkg{rosqp}, the underlying solver, and the options for \code{\link[rosqp]{solve_osqp}}.

code{\link[rosqp]{osqpSettings}} for details on options for \code{solve_osqp}.
}

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
ow3 <- optweight(list(tx1 ~ use0 + gender + age,
                         tx2 ~ tx1 + use1 + use0  + gender + age,
                         tx3 ~ tx2 + use2 + tx1 + use1 + use0 + gender + age),
                    data = iptwExWide,
                    tols = list(c(.001, .001, .001),
                                c(.001, .001, .01, .01, .01),
                                c(.001, .001, .01, .01, .1, .1, .1)))

}
