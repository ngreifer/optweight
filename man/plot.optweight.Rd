\name{plot.optweight}
\alias{plot.optweight}
\title{
Plot Dual Variables for Assessing Balance Constraints
}
\description{
Plots the dual variables resulting from \code{optweight} in a way similar to figure 2 of Zubizarreta (2015), which explained how to interpret these values. These can be thought of as the cost of tightening the constraint on the variance of the resulting weights. For covariates with large values of the dual variable, tightening the constraint will increase the variability of the weights, and loosening the constraint will decrease the variability fo the weights, both to a greater extent than would doing the same for covariate with small values of the dual variable.
}
\usage{
\method{plot}{optweight}(x, which.time = 1, ...)
}
\arguments{
  \item{x}{
An \code{optweight} object; the output of a call to \code{\link{optweight}}.
}
  \item{which.time}{
For longitudinal treatments, which time period to display. Only one may be displayed at a time.
}
  \item{\dots}{
Ignored.
}
}
\details{
If binary or multinomial treatments are used, there will be a seperate plot for each treatment (other than the focal group if the ATT or ATC are requested). Because constraints are specified simutaneously for each treatment, the dual variables for all treatments should be considered when adjusting the tolerance for a covariate. For example, if the dual variable for a covariate is near 0 for the treated group but large for the control group, its tolerance should not be tightened because one of the constraints for the covariate is active at the optimum.

As recommended by Zubizarreta (2015), the dual variables are standardized to the value of the objective function at the optimum to enhance interpretability. A consequence of this is that the values of the standardized dual variables cannot be compared across different tolerance specifications because the value of the objective function changes.
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
\code{\link{optweight}} to estimate the weights and the dual variables
\code{\link{plot.summary.optweight}} for plots of the distribution of weights
}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

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
