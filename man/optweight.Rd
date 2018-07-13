\name{optweight}
\alias{optweight}
\alias{print.optweight}
\alias{print.optweightMSM}
\title{Estimate Balancing Weights Using Optimization}
\usage{
optweight(formula,
          data = NULL,
          tols = .0001,
          estimand = "ATE",
          s.weights = NULL,
          focal = NULL)

\method{print}{optweight}(x, ...)

\method{print}{optweightMSM}(x, ...)

optweight.fit(treat,
              covs,
              tols = .001,
              estimand = "ATE",
              s.weights = NULL,
              focal = NULL,
              std.binary = FALSE)
}
\description{
Estimate balancing weights for treatments and covariates specified in \code{formula}. The degree of balance for each covariate is specified by \code{tols}. \code{optweight.fit} performs the optimization (via \pkg{limSolve}) and outputs a vector of weights, and \code{optweight} prepares input for processing and outputs an object that includes the estimated weights and other values. Users should use \code{optweight} for most purposes.

}
\arguments{
  \item{formula}{
a formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side, or a list thereof. See \code{\link{glm}} for more details. Interactions and functions of covariates are allowed.
}
  \item{data}{
an optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{tols}{
a vector of balance tolerance values in standardized mean difference units for each covariate, or a list thereof. The resulting weighted standardized mean differences will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. For binary variables, the values correspond to raw difference in proportion rather than standardized mean differences.
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
  \item{treat}{
a vector of treatment statuses for each unit or a list thereof. Non-numeric (i.e., factor or character) values are allowed.
}
  \item{covs}{
a data frame or matrix of covariates to be balanced. All variables must be numeric. Does not have to be full rank.
}
  \item{std.binary}{
whether the tolerances are in standardized mean units (\code{TRUE}) or raw units (\code{FALSE}) for binary variables. The default is \code{FALSE} because raw proportion differences make more sense than standardized mean difference for binary variables. This argument is analogous to the \code{binary} argument in \code{bal.tab} in \pkg{cobalt}.
}
  \item{x}{
an \code{optweight} or \code{optweightMSM} object; the output of a call to \code{optweight()}.
}
  \item{\dots}{
ignored.
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
The optimization is fit using \code{\link{lsei}} in the \pkg{limSolve} package, which provides a simple interface to specifying the constraints and objective function for quadratic optimization problems.

For longitudinal treatments, only "wide" data sets, where each row corresponds to a unit's entire variable history, are supported. You can use \code{\link{reshape}} or other functions to transform your data into this format; see example in the documentation for \code{weightitMSM} in the \pkg{WeightIt} package.
}

\author{
Noah Greifer
}

\references{

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

#Balancing covariates with respect to race (multinomial)
(ow2 <- optweight(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = .01, estimand = "ATT", focal = "black"))
bal.tab(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                weights = ow2$weights, focal = "black")

#Balancng covariates with longtidinal treatments
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
