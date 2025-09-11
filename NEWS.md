optweight News and Updates
======

# optweight (development version)

* Moved some documentation around; `optweight()` and `optweight.fit()` are now documented on the same page, `optweightMV()` and `optweightMV.fit()` are now documented on the same page, and `optweight.svy()` and `optweight.svy.fit()` are now documented on the same page.

* `min.w` is now a visible argument to `optweight()`, `optweightMV()`, and `optweight.svy()`.

* Messages are now a little prettier thanks to *crayon*, which is a new dependency.

# optweight 1.0.0

The version involved a full rewrite and may not be backward compatible with prior versions. Basic functionality has not changed, but some of the more advanced functionality has changed.

* Added a new function, `optweightMV()`, for performing weighting for multivariate treatments. This takes the place of the old functionality of `optweight()` to take in a list of formulas. Documentation and other supporting text has been updated to emphasize that the purpose of this function is for multivariate treatments, not longitudinal treatments as previously documented. Note that `optweight()` now throws an error if a list of formulas is supplied. A lower-level version, `optweightMV.fit()`, is also available.

* `check.tols()` and `check.targets()` have been renamed `process_tols()` and `process_targets()`, respectively, and have some new functionality. The `stop` argument has been removed. `process_targets()` now accepts `s.weights` to compute sampling-weighted target means.

* `norm` can be set to `"entropy"` in `optweight()`, etc., which minimizes the negative entropy of the weights. This is equivalent to entropy balancing, which is implemented more efficiently in *WeightIt*, but has support for inexact balance and multivariate treatments.

* `norm` can also be set to `"log"` in `optweight()`, etc., which minimizes the sum of the negative log of the weights. This is equivalent to nonparametric covariate balancing propensity score (npCBPS) weighting, a version of which is implemented in *CBPS*. The implementation here has support for inexact balance and multivariate treatments.

* Negative values are now allowed for the `min.w` argument.

* Different solvers can be used by supplying an argument to `solver`. See `?optweight.fit()` for defaults and allowable options.

* `optweight.svy()`, `process_tols()`, and `process_targets()` can all be supplied with without a formula, which uses all variables in the supplied data.

* Documentation is now in *roxygen2*.

* Added a new vignette (`vignette("optweight")`).

* Added `b.weights` argument to supply base weights. When supplied, rather than minimizing the variance of the weights, the squared distance from each base weight is minimized, mirroring the functionality of the `base.weights` argument in *ebalance* for entropy balancing. All norms now support base weights. Omitting base weights is equivalent to setting them equal to 1.

* Optimization proceeds differently with `s.weights` supplied. First, the L$\infinity$ cannot be used with sampling weights. Second, the norm minimized is the weighted norm of the difference between the estimated weights and base weights, with the estimated weights not incorporating the sampling weights. That is, the L2 norm minimizes $\sum_i s_i(w_i-b_i)^2$, and the L1 norm minimizes $\sum_i s_i|w_i-b_i|$, where $s_i$ is the sampling weight for unit $i$, $b_i$ is the base weight (1 by default), and $w_i$ are the weights to be estimated. The weights used in the balance constraints (and ultimately in effect estimation) are $w^*_i=s_i w_i$. An implication of this is that the ESS of the $w^*_i$ is not maximized with the L2 norm. This also ensures that the weighted bootstrap correctly accounts for estimation of the weights.

* `summary()` now displays the L2, L1, and L$\infinity$ norms and the relative entropy between the estimated weights and the base weights, and the number of weights estimated to be 0. The L2 and L1 norms and relative entropy are weighted by the `s.weights` if present.

* `polish` is now `TRUE` by default for norms other than `"entropy"`; this slightly improves estimation.

* Some default arguments to the solvers have changed.

* Formula interfaces now accept `poly(x, .)` and other matrix-generating functions of variables, including the `rms`-class-generating functions from the *rms* package (e.g., `pol()`, `rcs()`, etc.) (the *rms* package must be loaded to use these latter ones) and the `basis`-class-generating functions from the *splines* package (i.e., `bs()` and `ns()`). A bug in an early version of this was found by @ahinton-mmc.

* The returned covariates are now those without any transformations.

* Updated the README.

* Added a new logo.

# optweight 0.2.5

* Reverting back to using *osqp* instead of *rosqp* now that *osqp* works. *cobalt* is back.

# optweight 0.2.4

* Reverting back to using *rosqp* instead of *osqp* due to package failure. Also removed reliance on *cobalt* in favor of *MatchIt* for data. Both changes are temporary. 

# optweight 0.2.3

* The *rosqp* package is now *osqp*, and is faster with fewer bugs.

* If `focal` is set, the estimand is automatically changed to `"ATT"`. In the past, `focal` would be ignored unless `estimand = "ATT"`.

* Fixed some bugs with processing formula inputs. In particular, functions can be used inside `lapply()` loops and nested functions more gracefully.

* Other bugs fixes and small changes.

# optweight 0.2.2

* Fixed bug with duals displaying improperly when factor variables are present.

# optweight 0.2.1

* Changed default `min.w` in `optweight.fit()` and `optweight.svy.fit()` to 1E-8 from 0. This ensures all weights are nonzero, which can reduce bugs in other functions that require nonzero weights (e.g, `jtools::summ()` and `survey::svyglm()`).

* Fixed warning that would occur when interactions were present in the model formula in `optweight()`.

* Stable balancing weights have been discovered to be invalid for longitudinal treatments, so attempting to use `optweight()`or `optweight.fit()` with longitudinal treatments will now produce an error. This can be overridden by setting `force = TRUE`, though this is not recommended until further research is done.

# optweight 0.2.0

* Added `optweight.svy()` and associated methods and functions for estimating survey weights using optimization. These weights when applied to the sample will yield a sample whose covariate means are equal (within a specified tolerance) to given target values.

* Minor changes to `check.targets()`. It will now produce the covariate means when the `targets` argument is empty and will produce the previous empty output, a named vector of `NA`s, when `targets = NULL`.

* Changes to how dual variables are processed and displayed. Now, each dual variable coming from `optweight()` represents the change in the objective function corresponding to a 1-unit change in `tols`. The reported duals are the sum of all the duals affected by the constraint, so you can now reliably predict the change in the objective function from a change in `tols` (it was obscured and error-prone previously). The distinction between targeting duals and balance duals is maintained.

# optweight 0.1.0

* First version!
