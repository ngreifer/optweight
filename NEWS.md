optweight News and Updates
======

Version 0.2.1

* Changed default `min.w` in `optweight.fit()` and `optweight.svy.fit()` to 1E-8 from 0. This ensures all weights are nonzero, which can reduce bugs in other functions that require nonzero weights (e.g, `summ()` in `jtools` and `svyglm()` in survey`).

Version 0.2.0

* Added `optweight.svy` and associated methods and functions for estimating survey weights using optimization. These weights when applied to the sample will yield a sample whose covariate means are equal (within a specified tolerance) to given target values.

* Minor changes to `check.targets`. It will now produce the covariate means when the `targets` argument is empty and will produce the previous empty output, a named vector of `NA`s, when `targets = NULL`.

* Changes to how dual variables are processed and displayed. Now, each dual variable coming from `optweight` represents the change in the objective function corresponding to a 1-unit change in `tols`. The reported duals are the sum of all the duals affected by the constraint, so you can now reliably predict the change in the objective function from a change in `tols` (it was obscured and error-prone previously). The distinction between targeting duals and balance duals is maintained.

Version 0.1.0

* First version!
