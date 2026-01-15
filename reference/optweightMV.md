# Stable Balancing Weights for Multivariate Treatments

Estimates stable balancing weights for the supplied multivariate (i.e.,
multiple) treatments and covariates. The degree of balance for each
covariate is specified by `tols.list`. See Zubizarreta (2015) and Wang &
Zubizarreta (2020) for details of the properties of the weights and the
methods used to fit them.

## Usage

``` r
optweightMV(
  formula.list,
  data = NULL,
  tols.list = list(0),
  estimand = "ATE",
  targets = NULL,
  target.tols.list = list(0),
  s.weights = NULL,
  b.weights = NULL,
  norm = "l2",
  min.w = 1e-08,
  verbose = FALSE,
  ...
)

optweightMV.fit(
  covs.list,
  treat.list,
  tols.list = list(0),
  estimand = "ATE",
  targets = NULL,
  target.tols.list = list(0),
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

- formula.list:

  a list of formulas, each with a treatment variable on the left hand
  side and the covariates to be balanced on the right hand side.

- data:

  an optional data set in the form of a data frame that contains the
  variables in `formula.list`.

- tols.list:

  a list of vectors of balance tolerance values for each covariate for
  each treatment. The resulting weighted balance statistics will be at
  least as small as these values. If only one value is supplied, it will
  be applied to all covariates. See Details. Default is 0 for all
  covariates.

- estimand:

  the desired estimand, which determines the target population. Only
  "ATE" or `NULL` are supported. `estimand` is ignored when `targets` is
  non-`NULL`. If both `estimand` and `targets` are `NULL`, no targeting
  will take place.

- targets:

  an optional vector of target population mean values for each
  covariate. The resulting weights ensure the midpoint between group
  means are within `target.tols` units of the target values for each
  covariate. If `NULL` or all `NA`, `estimand` will be used to determine
  targets. Otherwise, `estimand` is ignored. If any target values are
  `NA`, the corresponding variable will not be targeted and its weighted
  mean will be wherever the weights yield the smallest value of the
  objective function; this is only allowed if all treatments are binary
  or multi-category. Can also be the output of a call to
  [`process_targets()`](https://ngreifer.github.io/optweight/reference/process_targets.md).
  See Details.

- target.tols.list:

  a list of vectors of target balance tolerance values for each
  covariate for each treatment. For binary and multi-category
  treatments, the average of each pair of means will be at most as far
  from the target means as these values. Can also be the output of a
  call to
  [`process_tols()`](https://ngreifer.github.io/optweight/reference/process_tols.md).
  See Details. Default is 0 for all covariates. Ignored with continuous
  treatments.

- s.weights:

  a vector of sampling weights. For `optweightMV()`, can also be the
  name of a variable in `data` that contains sampling weights.

- b.weights:

  a vector of base weights. If supplied, the desired norm of the
  distance between the estimated weights and the base weights is
  minimized. For `optweightMV()`, can also the name of a variable in
  `data` that contains base weights.

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

  for `optweightMV()`, additional arguments passed to
  `optweightMV.fit()`, including options that are passed to the settings
  function corresponding to `solver`.

- covs.list:

  a list containing one numeric matrix of covariates to be balanced for
  each treatment.

- treat.list:

  a list containing one vector of treatment statuses for each treatment.

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

For `optweightMV()`, an `optweightMV` object with the following
elements:

- weights:

  The estimated weights, one for each unit.

- treat.list:

  A list of the values of the treatment variables.

- covs.list:

  A list of the covariates for each treatment used in the fitting. Only
  includes the raw covariates, which may have been altered in the
  fitting process.

- s.weights:

  The provided sampling weights.

- b.weights:

  The provided base weights.

- call:

  The function call.

- tols:

  A list of tolerance values for each covariate for each treatment.

- duals:

  A list of data.frames containing the dual variables for each covariate
  for each treatment. See
  [`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md)
  for interpretation of these values.

- info:

  A list containing information about the performance of the
  optimization at termination.

- norm:

  The `norm` used.

- solver:

  The `solver` used.

For `optweightMV.fit()`, an `optweightMV.fit` object with the following
elements:

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

`optweightMV()` is the primary user-facing function for estimating
stable balancing weights for multivariate treatments. The optimization
is performed by the lower-level function `optweightMV.fit()`, which
transforms the inputs into the required inputs for the optimization
functions and then supplies the outputs (the weights, dual variables,
and convergence information) back to `optweightMV()`. Little processing
of inputs is performed by `optweightMV.fit()`, as this is normally
handled by `optweightMV()`.

See
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md)
for more information about balance tolerances (i.e., those specified in
`tols.list`), `targets`, `norm`, `solver`, and convergence failure.

## References

Chattopadhyay, A., Cohn, E. R., & Zubizarreta, J. R. (2024). One-Step
Weighting to Generalize and Transport Treatment Effect Estimates to a
Target Population. *The American Statistician*, 78(3), 280–289.
[doi:10.1080/00031305.2023.2267598](https://doi.org/10.1080/00031305.2023.2267598)

Källberg, D., & Waernbaum, I. (2023). Large Sample Properties of Entropy
Balancing Estimators of Average Causal Effects. *Econometrics and
Statistics*.
[doi:10.1016/j.ecosta.2023.11.004](https://doi.org/10.1016/j.ecosta.2023.11.004)

Wang, Y., & Zubizarreta, J. R. (2020). Minimal dispersion approximately
balancing weights: Asymptotic properties and practical considerations.
*Biometrika*, 107(1), 93–105.
[doi:10.1093/biomet/asz050](https://doi.org/10.1093/biomet/asz050)

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for
Estimation With Incomplete Outcome Data. *Journal of the American
Statistical Association*, 110(511), 910–922.
[doi:10.1080/01621459.2015.1023805](https://doi.org/10.1080/01621459.2015.1023805)

## See also

[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md)
for more information on the optimization, specifications, and options.

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

# Balancing two treatments
(ow1 <- optweightMV(list(treat ~ age + educ + race + re74,
                         re75 ~ age + educ + race + re74),
                    data = lalonde))
#> An optweightMV object
#>  - number of obs.: 614
#>  - norm minimized: "l2"
#>  - sampling weights: present
#>  - base weights: present
#>  - number of treatments: 2
#>     treat: 2-category
#>     re75: continuous
#>  - covariates: 
#>     + for treat: age, educ, race, re74
#>     + for re75: age, educ, race, re74

summary(ow1)
#>                   Summary of weights
#> ─────────────────── Treatment 1 ───────────────────
#> - Weight ranges:
#> 
#>         Min                                 Max
#> treated   0 |---------------------------| 8.786
#> control   0 |--------------------|        6.704
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>            179   166   162   124    23
#>  treated 5.682 5.913   6.4 6.819 8.786
#>            300    48    26    19    15
#>  control 3.854 3.955 3.991 5.418 6.704
#> 
#> 
#> - Weight statistics:
#> 
#>            L2    L1    L∞ Rel Ent # Zeros
#> treated 1.783 1.35  7.786   1.276       0
#> control 0.955 0.743 5.704   0.477       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.    185.  
#> Weighted    224.31   44.26
#> ─────────────────── Treatment 2 ───────────────────
#> - Weight ranges:
#> 
#>     Min                                 Max
#> all   0 |---------------------------| 8.786
#> 
#> - Units with the 5 most extreme weights:
#>                                 
#>        200 179   166   124    23
#>  all 5.913 6.4 6.704 6.819 8.786
#> 
#> 
#> - Weight statistics:
#> 
#>        L2    L1    L∞ Rel Ent # Zeros
#> all 1.263 0.926 7.786   0.718       0
#> 
#> - Effective Sample Sizes:
#> 
#>             Total
#> Unweighted 614.  
#> Weighted   236.54

bal.tab(ow1)
#> Balance by Time Point
#> 
#>  - - - Time: 1 - - - 
#> Balance Measures
#>                Type Diff.Adj
#> age         Contin.       -0
#> educ        Contin.       -0
#> race_black   Binary       -0
#> race_hispan  Binary        0
#> race_white   Binary        0
#> re74        Contin.       -0
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    224.31   44.26
#> 
#>  - - - Time: 2 - - - 
#> Balance Measures
#>                Type Corr.Adj Diff.Target.Adj
#> age         Contin.       -0               0
#> educ        Contin.       -0               0
#> race_black   Binary       -0               0
#> race_hispan  Binary       -0              -0
#> race_white   Binary        0              -0
#> re74        Contin.       -0               0
#> 
#> Effective sample sizes
#>             Total
#> Unadjusted 614.  
#> Adjusted   236.54
#>  - - - - - - - - - - - 
#> 
```
