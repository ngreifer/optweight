\name{check.tols}
\alias{check.tols}
\alias{print.optweight.tols}
\title{Construct and Check Tolerance Input}
\usage{
check.tols(formula,
           data = NULL,
           tols,
           stop = FALSE)

\method{print}{optweight.tols}(x, internal = FALSE, digits = 5, ...)
}
\description{
Checks whether proposed tolerance values for \code{tols} are suitable in number and order for submission to \code{optweight}. Users should include one value per item in \code{formula}. The output can also be used as an input to \code{tols} in \code{optweight}.

}
\arguments{
  \item{formula}{
a formula with the covariates to be balanced with \code{optweight} on the right hand side. See \code{\link{glm}} for more details. Interactions and functions of covariates are allowed. Lists of formulas are not allows; multiple formulas must be checked one at a time.
}
  \item{data}{
an optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{tols}{
a vector of balance tolerance values in standardized mean difference units for each covariate. These should be in the order corresponding to the order of the corresponding variable in \code{formula}, except for interactions, which will appear after all lower-order terms. If only one value is supplied, it will be applied to all covariates.
}
  \item{stop}{
\code{logical}; if \code{TRUE}, an error will be thrown if the number of values in \code{tols} is not equal to the correct number of covariates in \code{formula}, and no messages will be displayed if the \code{tols} input is satisfactory. If \code{FALSE}, a message will be displayed if the numbe of values in \code{tols} is not equal to the correct number of covariates in \code{formula}, and other messages will be displayed.
}
  \item{x}{
an \code{optweight.tols} object; the output of a call to \code{check.tols}.
}
  \item{internal}{
\code{logical}; whether to print the tolerance values that are to be used internally by \code{optweight}. See Value section.
}
  \item{digits}{
how many digits to print.
}
  \item{\dots}{
ignored.
}
}
\value{
An \code{optweight.tols} object, which contains the following items:
\item{tols}{A named vector of tolerance values, one for each variable specified in \code{formula}. The should be used as user inputs to \code{optweight}.}
  \item{internal.tols}{The tolerance values to be used internally by \code{optweight}. These will differ from \code{tols} when there are factor variables that are split up; the user only needs to submit one tolerance per factor variable, but seperate tolerance values are produced for each new dummy created.}

}
\details{
The purpose of \code{check.tols} is to allow users to ensure that their proposed input to \code{tols} in \code{optweight} is correct both in the number of entries and their order. This is especially important when factor variables and interactions are included in the formula because factor variables are split into several dummies and interactions are moved to the end of the variable list, both of which can cause some confusion and potential error when entering \code{tols} values.

Factor variables are internally split into a dummy variable for each level, but the user only needs to specify one tolerance value per original variable; \code{check.tols} automatically expands the \code{tols} input to match the newly created variables.

Interactions (e.g., \code{a:b} or \code{a*b} in the \code{formula} input) are always sent to the end of the variable list even if they are specified elsewhere in the \code{formula}. It is important to run \code{check.tols} to ensure the order of the proposed \code{tols} corresponds to the represented order of covariates used in \code{optweight}. You can run \code{check.tols} with no \code{tols} input to see the order of covariates that is required.

\code{check.tols} was designed to be used primarily for its message printing and \code{print} method, but you can also assign its output to an object for use as an input to \code{tols} in \code{optweight}.

Note that only one formula and vector of tolerance values can be assessed at a time; for multiple treatment periods, each formula and tolerance vector must be entered seperately.
}

\author{
Noah Greifer
}

\references{

}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Checking if the correct number of entries are include:
check.tols(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                tols = c(.01, .02, .03, .04))

#Checking the order of interactions; notice they go
#at the end even if specified at the beginning. The
#.09 values are where the interactions might be expected
#to be, but they are in fact not.
c <- check.tols(treat ~ age:educ + married*race +
                nodegree + re74, data = lalonde,
                tols = c(.09, .01, .01, .09, .01, .01))

print(c, internal = TRUE)

}
