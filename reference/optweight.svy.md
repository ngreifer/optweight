# Stable Balancing Weights for Generalization

Estimates stable balancing weights to generalize a sample characterized
by supplied covariates to a given target population. The target means
are specified with `targets` and the maximum distance between each
weighted covariate mean. See Jackson et al. (2021) for details of the
properties of the weights and the methods used to fit them.

## Usage

``` r
optweight.svy(
  formula,
  data = NULL,
  tols = 0,
  targets = NULL,
  s.weights = NULL,
  b.weights = NULL,
  norm = "l2",
  min.w = 1e-08,
  verbose = FALSE,
  ...
)

optweight.svy.fit(
  covs,
  targets,
  tols = 0,
  s.weights = NULL,
  b.weights = NULL,
  norm = "l2",
  std.binary = FALSE,
  std.cont = TRUE,
  min.w = 1e-08,
  verbose = FALSE,
  solver = NULL,
  ...
)
```

## Arguments

- formula:

  a formula with nothing on the left hand side and the covariates to be
  targeted on the right hand side. Interactions and functions of
  covariates are allowed. Can be omitted, in which case all variables in
  `data` are assumed targeted. If `data` is `NULL` and `formula` is a
  data.frame, `data` will be replaced with `formula`.

- data:

  an optional data set in the form of a data frame that contains the
  variables in `formula`.

- tols:

  a vector of target balance tolerance values for each covariate. The
  resulting weighted covariate means will be no further away from the
  targets than the specified values. If only one value is supplied, it
  will be applied to all covariates. Can also be the output of a call to
  [`process_tols()`](https://ngreifer.github.io/optweight/reference/process_tols.md).
  Default is 0 for all covariates.

- targets:

  a vector of target population mean values for each covariate. The
  resulting weights will yield sample means within `tols` units of the
  target values for each covariate. If any target values are `NA`, the
  corresponding variable will not be targeted and its weighted mean will
  be wherever the weights yield the smallest variance. To ensure the
  weighted mean for a covariate is equal to its unweighted mean (i.e.,
  so that its original mean is its target mean), its original mean must
  be supplied as a target. For factor variables, a target value must be
  specified for each level of the factor, and these values must add up
  to 1. Can also be the output of a call to
  [`process_targets()`](https://ngreifer.github.io/optweight/reference/process_targets.md).

- s.weights:

  a vector of sampling weights. For
  [`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md),
  can also be the name of a variable in `data` that contains sampling
  weights.

- b.weights:

  a vector of base weights. If supplied, the desired norm of the
  distance between the estimated weights and the base weights is
  minimized. For
  [`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md),
  can also the name of a variable in `data` that contains base weights.

- norm:

  `character`; a string containing the name of the norm corresponding to
  the objective function to minimize. Allowable options include `"l1"`
  for the \\L_1\\ norm, `"l2"` for the \\L_2\\ norm (the default),
  `"linf"` for the \\L\_\infty\\ norm, `"entropy"` for the relative
  entropy, and `"log"` for the sum of the negative logs. See Details.

- min.w:

  `numeric`; a single value less than 1 for the smallest allowable
  weight. Some analyses require nonzero weights for all units, so a
  small, nonzero minimum may be desirable. The default is `1e-8`
  (\\10^{-8}\\), which does not materially change the properties of the
  weights from a minimum of 0 but prevents warnings in some packages
  that use weights in model fitting. When `norm` is `"entropy"` or
  `"log"` and `min.w <= 0`, `min.w` will be set to the smallest nonzero
  value.

- verbose:

  `logical`; whether information on the optimization problem solution
  should be printed. Default is `FALSE`.

- ...:

  for
  [`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md),
  additional arguments passed to
  [`optweight.fit()`](https://ngreifer.github.io/optweight/reference/optweight.md),
  including options that are passed to the settings function
  corresponding to `solver`.

- covs:

  a numeric matrix of covariates to be targeted.

- std.binary, std.cont:

  `logical`; whether the tolerances are in standardized mean units
  (`TRUE`) or raw units (`FALSE`) for binary variables and continuous
  variables, respectively. The default is `FALSE` for `std.binary`
  because raw proportion differences make more sense than standardized
  mean difference for binary variables. These arguments are analogous to
  the `binary` and `continuous` arguments in
  [`cobalt::bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html).

- solver:

  string; the name of the optimization solver to use. Allowable options
  depend on `norm`. Default is to use whichever eligible solver is
  installed, if any, or the default solver for the corresponding `norm`.
  See Details for information.

## Value

For `optweight.svy()`, an `optweight.svy` object with the following
elements:

- weights:

  The estimated weights, one for each unit.

- covs:

  The covariates used in the fitting. Only includes the raw covariates,
  which may have been altered in the fitting process.

- s.weights:

  The provided sampling weights.

- call:

  The function call.

- tols:

  The tolerance values for each covariate.

- duals:

  A data.frame containing the dual variables for each covariate. See
  [`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md)
  for interpretation of these values.

- info:

  A list containing information about the performance of the
  optimization at termination.

- norm:

  The `norm` used.

- solver:

  The `solver` used.

For `optweight.svy.fit()`, an `optweight.svy.fit` object with the
following elements:

- w:

  The estimated weights, one for each unit.

- duals:

  A data.frame containing the dual variables for each covariate.

- info:

  A list containing information about the performance of the
  optimization at termination.

- norm:

  The `norm` used.

- solver:

  The `solver` used.

## Details

`optweight.svy()` is the primary user-facing function for estimating
stable balancing weights for generalization to a target population. The
optimization is performed by the lower-level function
`optweight.svy.fit()`, which transforms the inputs into the required
inputs for the optimization functions and then supplies the outputs (the
weights, dual variables, and convergence information) back to
`optweight.svy()`. Little processing of inputs is performed by
`optweight.svy.fit()`, as this is normally handled by `optweight.svy()`.

Weights are estimated so that the standardized differences between the
weighted covariate means and the corresponding targets are within the
given tolerance thresholds (unless `std.binary` or `std.cont` are
`FALSE`, in which case unstandardized mean differences are considered
for binary and continuous variables, respectively). For a covariate
\\x\\ with specified tolerance \\\delta\\, the weighted mean will be
within \\\delta\\ of the target. If standardized tolerance values are
requested, the standardization factor is the standard deviation of the
covariate in the whole sample. The standardization factor is always
unweighted.

Target constraints are applied to the product of the estimated weights
and the sampling weights. In addition, sum of the product of the
estimated weights and the sampling weights is constrained to be equal to
the sum of the product of the base weights and sampling weights.

See
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md)
for information on `norm`, `solver`, and convergence failure.

## References

Jackson, D., Rhodes, K., & Ouwens, M. (2021). Alternative weighting
schemes when performing matching-adjusted indirect comparisons.
*Research Synthesis Methods*, 12(3), 333–346.
[doi:10.1002/jrsm.1466](https://doi.org/10.1002/jrsm.1466)

## See also

[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md)
for estimating weights that balance treatment groups.

[`process_targets()`](https://ngreifer.github.io/optweight/reference/process_targets.md)
for specifying the covariate target means supplied to `targets`.

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

cov.names <- c("age", "educ", "race",
               "married", "nodegree")

targets <- c(age = 23,
             educ = 9,
             race_black = .3,
             race_hispan = .3,
             race_white = .4,
             married = .2,
             nodegree = .5)

ows <- optweight.svy(lalonde[cov.names],
                     targets = targets)
ows
#> An optweight.svy object
#>  - number of obs.: 614
#>  - norm minimized: "l2"
#>  - sampling weights: present
#>  - base weights: present
#>  - covariates: age, educ, race, married, nodegree

# Unweighted means
col_w_mean(lalonde[cov.names])
#>         age        educ  race_black race_hispan  race_white     married 
#>  27.3631922  10.2687296   0.3957655   0.1172638   0.4869707   0.4153094 
#>    nodegree 
#>   0.6302932 

# Weighted means; same as targets
col_w_mean(lalonde[cov.names],
           w = ows$weights)
#>         age        educ  race_black race_hispan  race_white     married 
#>        23.0         9.0         0.3         0.3         0.4         0.2 
#>    nodegree 
#>         0.5 
```
