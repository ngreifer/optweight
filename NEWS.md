optweight News and Updates
======

# optweight (development version)

* Added `b.weights` argument to supply base weights. When supplied, rather than minimizing the variance of the weights, the squared distance from each base weight is minimized, mirroring the functionality of the `base.weights` argument in `ebalance` for entropy balancing.

* Formula interfaces now accept `poly(x, .)` and other matrix-generating functions of variables, including the `rms`-class-generating functions from the `rms` package (e.g., `pol()`, `rcs()`, etc.) (the `rms` package must be loaded to use these latter ones) and the `basis`-class-generating functions from the `splines` package (i.e., `bs()` and `ns()`). A bug in an early version of this was found by @ahinton-mmc.

* Negative values are now allowed for the `min.w` argument of `optweight.fit()`.

# optweight 0.2.5

* Reverting back to using `osqp` instead of `rosqp` now that `osqp` works. `cobalt` is back.

# optweight 0.2.4

* Reverting back to using `rosqp` instead of `osqp` due to package failure. Also removed reliance on `cobalt` in favor of `MatchIt` for data. Both changes are temporary. 

# optweight 0.2.3

* The `rosqp` package is now `osqp`, and is faster with fewer bugs.

* If `focal` is set, the estimand is automatically changed to `"ATT"`. In the past, `focal` would be ignored unless `estimand = "ATT"`.

* Fixed some bugs with processing formula inputs. In particular, functions can be used inside `lapply` loops and nested functions more gracefully.

* Other bugs fixes and small changes.

# optweight 0.2.2

* Fixed bug with duals displaying improperly when factor variables are present.

# optweight 0.2.1

* Changed default `min.w` in `optweight.fit()` and `optweight.svy.fit()` to 1E-8 from 0. This ensures all weights are nonzero, which can reduce bugs in other functions that require nonzero weights (e.g, `summ()` in `jtools` and `svyglm()` in survey`).

* Fixed warning that would occur when interactions were present in the model formula in `optweight()`.

* optweights have been discovered to be invalid for longitudinal treatments, so attempting to use `optweight()`or `optweight.fit()` with longitudinal treatments will now produce an error. This can be overridden by setting `force = TRUE`, though this is not recommended until further research is done.

# optweight 0.2.0

* Added `optweight.svy` and associated methods and functions for estimating survey weights using optimization. These weights when applied to the sample will yield a sample whose covariate means are equal (within a specified tolerance) to given target values.

* Minor changes to `check.targets`. It will now produce the covariate means when the `targets` argument is empty and will produce the previous empty output, a named vector of `NA`s, when `targets = NULL`.

* Changes to how dual variables are processed and displayed. Now, each dual variable coming from `optweight` represents the change in the objective function corresponding to a 1-unit change in `tols`. The reported duals are the sum of all the duals affected by the constraint, so you can now reliably predict the change in the objective function from a change in `tols` (it was obscured and error-prone previously). The distinction between targeting duals and balance duals is maintained.

# optweight 0.1.0

* First version!
