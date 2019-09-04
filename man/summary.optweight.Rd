\name{summary.optweight}
\alias{summary.optweight}
\alias{summary.optweightMSM}
\alias{summary.optweight.svy}
\alias{print.summary.optweight}
\alias{print.summary.optweightMSM}
\alias{print.summary.optweight.svy}
\alias{plot.summary.optweight}

\title{
Summarize, print, and plot information about estimated weights
}
\description{
These functions summarize the weights resulting from a call to \code{\link{optweight}} or \code{\link{optweight.svy}}. \code{summary} produces summary statistics on the distribution of weights, including their range and variability, and the effective sample size of the weighted sample (computing using the formula in McCaffrey, Rudgeway, & Morral, 2004). \code{plot} creates a histogram of the weights.
}
\usage{
\method{summary}{optweight}(object, top = 5, ignore.s.weights = FALSE, ...)

\method{summary}{optweightMSM}(object, top = 5, ignore.s.weights = FALSE, ...)

\method{summary}{optweight.svy}(object, top = 5, ignore.s.weights = FALSE, ...)

\method{print}{summary.optweight}(x, ...)

\method{print}{summary.optweightMSM}(x, ...)

\method{print}{summary.optweight.svy}(x, ...)

\method{plot}{summary.optweight}(x, ...)
}

\arguments{
  \item{object}{
An \code{optweight}, \code{optweightMSM}, or \code{optweight.svy} object; the output of a call to \code{\link{optweight}} or \code{\link{optweight.svy}}.
}
  \item{top}{
How many of the largest and smallest weights to display. Default is 5.
}
  \item{ignore.s.weights}{
Whether or not to ignore sampling weights when computing the weight summary. If \code{FALSE}, the default, the estimated weights will be multiplied by the sampling weights (if any) before values are computed.
}
  \item{x}{
A \code{summary.optweight}, \code{summary.optweightMSM}, or \code{summary.optweight.svy} object; the output of a call to \code{summary.optweight}, \code{summary.optweightMSM}, or \code{summary.optweight.svy}.
}
  \item{...}{
Additional arguments. For \code{plot}, additional arguments passed to \code{\link[graphics]{hist}} to determine the number of bins, though \code{\link[ggplot2]{geom_histogram}} from \pkg{ggplot2} is actually used to create the plot.
  }
}

\value{
For point treatments (i.e., \code{optweight} objects), \code{summary} returns a \code{summary.optweight} object with the following elements:
\item{weight.range}{The range (minimum and maximum) weight for each treatment group.}
\item{weight.top}{The units with the greatest weights in each treatment group; how many are included is determined by \code{top}.}
\item{coef.of.var}{The coefficient of variation (standard deviation divided by mean) of the weights in each treatment group and overall. When no sampling weights are used, this is simply the standard deviation of the weights.}
\item{mean.abs.dev}{The mean absolute deviation of the weights in each treatment group and overall.}
\item{effective.sample.size}{The effective sample size for each treatment group before and after weighting.}

For longitudinal treatments (i.e., \code{optweightMSM} objects), a list of the above elements for each treatment period.

For \code{optweight.svy} objects, a list of the above elements but with no treatment group divisions.

\code{plot} returns a \code{ggplot} object with a histogram displaying the distribution of the estimated weights. If the estimand is the ATT or ATC, only the weights for the non-focal group(s) will be displayed (since the weights for the focal group are all 1). A dotted line is displayed at the mean of the weights (usually 1).
}
\references{
McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004). Propensity Score Estimation With Boosted Regression for Evaluating Causal Effects in Observational Studies. Psychological Methods, 9(4), 403–425. \doi{10.1037/1082-989X.9.4.403}
}
\author{
Noah Greifer
}

\seealso{
\code{\link{plot.optweight}} for plotting the values of the dual variables.
}
\examples{
# library("cobalt")
# data("lalonde", package = "cobalt")

data("lalonde", package = "MatchIt")
lalonde$race <- factor(ifelse(lalonde$black == 1, "black",
                              ifelse(lalonde$hispan == 1, "hispan",
                                     "white")))

#Balancing covariates between treatment groups (binary)
(ow1 <- optweight(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = .001,
                estimand = "ATT"))

(s <- summary(ow1))

plot(s, breaks = 12)
}
