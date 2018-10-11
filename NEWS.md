optweight News and Updates
======

Version 0.2.0

* Added `optweight.svy` and associated methods and functions for estimating survey weights using optimization. These weights when applied to the sample will yield a sample whose covariate means are equal (within a specified tolerance) to given target values.

* Minor changes to `check.targets`. It will now produce the covariate means when the `targets` argument is empty and will produce the previous empty output, a named vector of `NA`s, when `targets = NULL`.

Version 0.1.0

* First version!
