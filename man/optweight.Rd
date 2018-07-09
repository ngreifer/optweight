\name{optweight}
\alias{optweight}
\title{Estimate Balancing Weights Using Optimization}
\usage{
optweight(formula,
         data = NULL,
         tols = .0001,
         estimand = "ATE",
         s.weights = NULL,
         focal = NULL,
         full.output = TRUE)

}
\description{
Estimate balancing weights for treatments and covariates specified in \code{formula}. The degree of balance for each covariate is specified by \code{tols}.
}
\arguments{
  \item{formula}{
a formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side. See \code{\link{glm}} for more details. Interactions and functions of covariates are allowed.
}
  \item{data}{
an optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{tols}{
a vector of balance tolerance values in standardized mean difference units for each covariate. The resulting weighted standardized mean differences will be at least as small as these values. If only one value is supplied, it will be applied to all covariates. For binary variables, the values correspond to raw difference in proportion rather than standardized mean differences.
}
  \item{estimand}{
the desired estimand. For binary treatments, can be "ATE", "ATT", "ATC". For multinomial treatments, can be "ATE" or "ATT". The default for both is "ATE".
}
  \item{s.weights}{
a vector of sampling weights or the name of a variable in \code{data} that contains sampling weights.
}
  \item{focal}{
when multinomial treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group.
}
}
\value{
  An \code{optweight} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat}{The values of the treatment variable.}
  \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{s.weights}{The provided sampling weights.}
  \item{estimand}{The estimand requested.}
  \item{focal}{The focal variable if the ATT was requested with a multinomial treatment.}
  \item{call}{The function call.}
}
\details{
The optimization is fit using \code{\link{lsei}} in the \pkg{limSolve} package, which provides a simple interface to specifying the constraints and objective function for quadratic optimization problems.
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

}
