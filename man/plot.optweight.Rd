\name{plot.optweight}
\alias{plot.optweight}
\alias{plot.optweight.svy}
\title{
Plot Dual Variables for Assessing Balance Constraints
}
\description{
Plots the dual variables resulting from \code{optweight} in a way similar to figure 2 of Zubizarreta (2015), which explained how to interpret these values. These represent the cost of changing the constraint on the variance of the resulting weights. For covariates with large values of the dual variable, tightening the constraint will increase the variability of the weights, and loosening the constraint will decrease the variability of the weights, both to a greater extent than would doing the same for covariate with small values of the dual variable.
}
\usage{
\method{plot}{optweight}(x, which.time = 1, ...)

\method{plot}{optweight.svy}(x, ...)
}
\arguments{
  \item{x}{
An \code{optweight} or \code{optweight.svy} object; the output of a call to \code{\link{optweight}} or \code{\link{optweight.svy}}.
}
  \item{which.time}{
For longitudinal treatments, which time period to display. Only one may be displayed at a time.
}
  \item{\dots}{
Ignored.
}
}
\value{
A \code{ggplot} object that can be used with other \pkg{ggplot2} functions.
}
\references{
Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
}
\author{
Noah Greifer
}

\seealso{
\code{\link{optweight}} or \code{\link{optweight.svy}} to estimate the weights and the dual variables
\code{\link{plot.summary.optweight}} for plots of the distribution of weights
}
\examples{
# library("cobalt")
# data("lalonde", package = "cobalt")

data("lalonde", package = "MatchIt")
lalonde$race <- factor(ifelse(lalonde$black == 1, "black",
                              ifelse(lalonde$hispan == 1, "hispan",
                                     "white")))

#Balancing covariates between treatment groups (binary)
ow1 <- optweight(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = c(.1, .1, .1, .1, .1),
                estimand = "ATT")

summary(ow1) # Note the coefficient of variation
             # and effective sample size (ESS)

plot(ow1) # age has a low value, married is high

ow2 <- optweight(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = c(0, .1, .1, .1, .1),
                estimand = "ATT")

summary(ow2) # Notice that tightening the constraint
             # on age had a negligible effect on the
             # variability of the weights and ESS

ow3 <- optweight(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = c(.1, .1, 0, .1, .1),
                estimand = "ATT")

summary(ow3) # In contrast, tightening the constraint
             # on married had a large effect on the
             # variability of the weights, shrinking
             # the ESS
}
