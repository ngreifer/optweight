
<!-- README.md is generated from README.Rmd. Please edit that file -->

# optweight

[![CRAN\_Status\_Badge](http://r-pkg.org/badges/version-last-release/optweight?color=0047ab)](https://cran.r-project.org/package=optweight)
[![CRAN\_Downloads\_Badge](http://cranlogs.r-pkg.org/badges/optweight?color=0047ab)](https://cran.r-project.org/package=optweight)

I’ll provide a better description later but for now, `optweight`
contains the function `optweight()` which generates weights to balance
treatments to given balance thresholds. It solves a quadratic
programming problem to minimize the variance of the weights using
`lsei()` in the `limSolve` package. This is the method described in
Zubizarreta (2015).

In addition to binary point treatments, `optmatch` can handle
multinomial treatments; longitudinal treatments are coming soon.

If you are using a recent version of `WeightIt` that supports
user-defined methods, you can enter `method = optmatch.fit` to use the
method with binary, multinomial, and longitudinal treatments. Please let
me know if you find any bugs or would like to see additional features.
