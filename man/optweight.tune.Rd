\name{optweight.tune}
\alias{optweight.tune}
\title{Tuning the tolerance}
\usage{
optweight.tune(formula,
               data = NULL,
               candidate.tols = list(.1, .01, .001, .0001),
               nboot = 100,
               bal.tab.control = list(),
               ...)
}
\description{
Implements the tuning algorithm described in Wang & Zubizarreta (2017) for determining an optimal tolerance value to be used for all covariates. The algorithm works by estimating stable balancing weights at each candidate tolerance level, and finding the set of weights that best balances bootstrap samples drawn from the original sample. Wang & Zubizarreta (2017) claim that this tolerance value will yield effect estimates with the lowest RMSE.

\code{optweight.tune} works by running \code{optweight} on the provided data using each of the candidate tolerance values, drawing \code{nboot} bootstrap samples form the original data, computing an overall balance statistic for each bootstrapped sample, taking the average of those balance statistics across the bootstrapped samples, and choosing the weights with lowest value of these averages. The \code{optweight} object corresponding to weights estimated at the selected best tolerance value is returned.

}
\arguments{
  \item{formula}{
A formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side, or a list thereof. See \code{\link{glm}} for more details. Interactions and functions of covariates are allowed.
}
  \item{data}{
An optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{candidate.tols}{
A list of candidate tolerance values to compare. \code{optweight} will be run once for each entry provided. Each entry should be a suitable argument to \code{tols} in \code{optweight}, i.e., a single number or a numeric vector with length equal to the number of covariates.
}
  \item{nboot}{
The number of bootstrap samples to draw. Should be fairly large to get consistent results. Only balance computation (using \code{bal.tab} in \pkg{cobalt}) is performed, not optimization, so this number can be set somewhat high without too much of a computational burden.
}
  \item{bal.tab.control}{
A list of arguments to \code{bal.tab} in \pkg{cobalt} used for assessing balance on the bootstrapped samples. These might include arguments like \code{s.d.denom} or \code{binary} that control how balance statistics are computed.
}
  \item{\dots}{
Other arguments to \code{\link{optweight}}.
}

}
\value{
  The \code{optweight} object resulting from the candiate tolerance entry that yields the best balance averaged across the bootstrapped samples. The \code{tols} entry can be checked to see which tolerance entry was selected.

}
\details{
Within each bootstrapped sample, the balance statistic is computed as the square root of the mean of the squared primary balance measures. The primary balance measure is the standardized mean difference for binary treatments, the maximum standardized mean difference across treatment pairs for multinomial treatments, and the treatment-covairtae correlation for continuous treatments. How these values are ocmputed depends on the options in \code{bal.tab}, which are controlled by \code{bal.tab.control}. The relevant arguments are as follows:

\bold{Binary treatments}:
\describe{
\item{\code{s.d.denom}}{how the denominator of the standardized mean differences should be computed. The default is is to use the standard deviation of the target group. See \code{bal.tab} for more details.
}
\item{\code{binary}}{whether the balance measure for binary covariates should be the standardized mean difference or the raw difference in proportion. The default is to use the raw difference in proportion. If \code{std.binary} is set to \code{TRUE}, \code{binary} should be set to \code{"std"} in \code{bal.tab.control}. If it isn't, a warning will appear.
}
\item{\code{continuous}}{whether the balance measure for continuous covariates should be the standardized mean difference or the raw mean difference. The default is to use the standardized mean difference. If \code{std.cont} is set to \code{FALSE}, \code{continuous} should be set to \code{"raw"} in \code{bal.tab.control}. If it isn't, a warning will appear.
}
}

\bold{Multinomial treatments}:
\describe{
\item{The same as for binary treatments, plus:}{}
\item{\code{pairwise}}{whether the balance statistcis should be computed pairwise, i.e., between pairs of treatments, or not, i.e., between each treatment and the other groups combined, before the mximum is taken. The default is to use pairwise comparisons.
}
}

\bold{Continuous treatments}:
\describe{
\item{None.}{}
}

Other arguments will be ignored.

}

\author{
Noah Greifer
}

\references{
Wang, Y., & Zubizarreta, J. R. (2017). Minimal Approximately Balancing Weights: Asymptotic Properties and Practical Considerations. ArXiv:1705.00998v2 [Math, Stat]. Retrieved from \url{http://arxiv.org/abs/1705.00998v2}
}

\seealso{
\code{\link{optweight}}
}

\examples{
library("cobalt")
data("lalonde", package = "cobalt")
set.seed(999)

#Balancing covariates between treatment groups (binary)
(ow <- optweight.tune(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                candidate.tols = list(.01, .001, .0001),
                nboot = 100,
                estimand = "ATE"))
ow$tols
bal.tab(ow)
}
