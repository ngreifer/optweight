# Stable Balancing Weights

Estimates stable balancing weights for the supplied treatments and
covariates. The degree of balance for each covariate is specified by
`tols` and the target population can be specified with `targets` or
`estimand`. See Zubizarreta (2015) and Wang & Zubizarreta (2020) for
details of the properties of the weights and the methods used to fit
them.

## Usage

``` r
optweight(
  formula,
  data = NULL,
  tols = 0,
  estimand = "ATE",
  targets = NULL,
  target.tols = 0,
  s.weights = NULL,
  b.weights = NULL,
  focal = NULL,
  norm = "l2",
  min.w = 1e-08,
  verbose = FALSE,
  ...
)

optweight.fit(
  covs,
  treat,
  tols = 0,
  estimand = "ATE",
  targets = NULL,
  target.tols = 0,
  s.weights = NULL,
  b.weights = NULL,
  focal = NULL,
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

  a formula with a treatment variable on the left hand side and the
  covariates to be balanced on the right hand side, or a list thereof.
  Interactions and functions of covariates are allowed.

- data:

  an optional data set in the form of a data frame that contains the
  variables in `formula`.

- tols:

  a vector of balance tolerance values for each covariate. The resulting
  weighted balance statistics will be at least as small as these values.
  If only one value is supplied, it will be applied to all covariates.
  Can also be the output of a call to
  [`process_tols()`](https://ngreifer.github.io/optweight/reference/process_tols.md).
  See Details. Default is 0 for all covariates.

- estimand:

  a string containing the desired estimand, which determines the target
  population. For binary treatments, can be "ATE", "ATT", "ATC", or
  `NULL`. For multi-category treatments, can be "ATE", "ATT", or `NULL`.
  For continuous treatments, can be "ATE" or `NULL`. The default for
  both is "ATE". `estimand` is ignored when `targets` is non-`NULL`. If
  both `estimand` and `targets` are `NULL`, no targeting will take
  place. See Details.

- targets:

  an optional vector of target population mean values for each
  covariate. The resulting weights ensure the midpoint between group
  means are within `target.tols` units of the target values for each
  covariate. If `NULL` or all `NA`, `estimand` will be used to determine
  targets. Otherwise, `estimand` is ignored. If any target values are
  `NA`, the corresponding variable will not be targeted and its weighted
  mean will be wherever the weights yield the smallest value of the
  objective function; this is only allowed for binary and multi-category
  treatments. Can also be the output of a call to
  [`process_targets()`](https://ngreifer.github.io/optweight/reference/process_targets.md).
  See Details.

- target.tols:

  a vector of target balance tolerance values for each covariate. For
  binary and multi-category treatments, the average of each pair of
  means will be at most as far from the target means as these values.
  Can also be the output of a call to
  [`process_tols()`](https://ngreifer.github.io/optweight/reference/process_tols.md).
  See Details. Default is 0 for all covariates. Ignored with continuous
  treatments and when `estimand` is `"ATT"` or `"ATC"`.

- s.weights:

  a vector of sampling weights. For `optweight()`, can also be the name
  of a variable in `data` that contains sampling weights.

- b.weights:

  a vector of base weights. If supplied, the desired norm of the
  distance between the estimated weights and the base weights is
  minimized. For `optweight()`, can also the name of a variable in
  `data` that contains base weights.

- focal:

  when multi-category treatments are used and `estimand = "ATT"`, which
  group to consider the "treated" or focal group. This group will not be
  weighted, and the other groups will be weighted to be more like the
  focal group. If specified, `estimand` will automatically be set to
  `"ATT"`.

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

  for `optweight()`, additional arguments passed to `optweight.fit()`,
  including options that are passed to the settings function
  corresponding to `solver`.

- covs:

  a numeric matrix of covariates to be balanced.

- treat:

  a vector of treatment statuses. Non-numeric (i.e., factor or
  character) vectors are allowed.

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

For `optweight()`, an `optweight` object with the following elements:

- weights:

  The estimated weights, one for each unit.

- treat:

  The values of the treatment variable.

- covs:

  The covariates used in the fitting. Only includes the raw covariates,
  which may have been altered in the fitting process.

- s.weights:

  The provided sampling weights.

- b.weights:

  The provided base weights.

- estimand:

  The estimand requested.

- focal:

  The focal variable if the ATT was requested with a multi-category
  treatment.

- call:

  The function call.

- tols:

  The balance tolerance values for each covariate.

- target.tols:

  The target balance tolerance values for each covariate.

- duals:

  A data.frame containing the dual variables for each covariate. See
  Details for interpretation of these values.

- info:

  A list containing information about the performance of the
  optimization at termination.

- norm:

  The `norm` used.

- solver:

  The `solver` used.

For `optweight.fit()`, an `optweight.fit` object with the following
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

`optweight()` is the primary user-facing function for estimating stable
balancing weights. The optimization is performed by the lower-level
function `optweight.fit()`, which transforms the inputs into the
required inputs for the optimization functions and then supplies the
outputs (the weights, dual variables, and convergence information) back
to `optweight()`. Little processing of inputs is performed by
`optweight.fit()`, as this is normally handled by `optweight()`.

For binary and multi-category treatments, weights are estimated so that
the weighted mean differences of the covariates are within the given
tolerance thresholds controlled by `tols` and `target.tols` (unless
`std.binary` or `std.cont` are `TRUE`, in which case standardized mean
differences are considered for binary and continuous variables,
respectively). For a covariate \\x\\ with specified balance tolerance
\\\delta\\ and target tolerance \\\varepsilon\\, the weighted means of
each each group will be within \\\delta\\ of each other, and the
midpoint between the weighted group means will be with \\\varepsilon\\
of the target means. More specifically, the constraints are specified as
follows: \$\$ \left\| \bar{x}^w_1 - \bar{x}^w_0 \right\| \le \delta \\
\left\| \frac{\bar{x}^w_1 + \bar{x}^w_0}{2} - \bar{x}^\* \right\| \le
\varepsilon \$\$ where \\\bar{x}^w_1\\ and \\\bar{x}^w_0\\ are the
weighted means of covariate \\x\\ for treatment groups 1 and 0,
respectively, and \\\bar{x}^\*\\ is the target mean for that covariate.
\\\delta\\ corresponds to `tols`, and \\\varepsilon\\ corresponds to
`target.tols`. Setting a covariate's value of `target.tols` to `Inf` or
its `target` to `NA` both serve to remove the second constraint, as is
done in Barnard et al. (2025).

If standardized tolerance values are requested, the standardization
factor corresponds to the estimand requested: when the ATE is requested
or a target population specified, the standardization factor is the
square root of the average variance for that covariate across treatment
groups, and when the ATT or ATC are requested, the standardization
factor is the standard deviation of the covariate in the focal group.
The standardization factor is computed accounting for `s.weights`.

Target and balance constraints are applied to the product of the
estimated weights and the sampling weights. In addition, the sum of the
product of the estimated weights and the sampling weights is constrained
to be equal to the sum of the product of the base weights and sampling
weights. For binary and multi-category treatments, these constraints
apply within each treatment group.

### Continuous treatments

For continuous treatments, weights are estimated so that the weighted
correlation between the treatment and each covariate is within the
specified tolerance threshold. The means of the weighted covariates and
treatment are restricted to be exactly equal to those of the target
population to ensure generalizability to the desired target population,
regardless of `tols` or `target.tols`. The weighted correlation is
computed as the weighted covariance divided by the product of the
*unweighted* standard deviations. The means used to center the variables
in computing the covariance are those specified in the target
population.

### `norm`

The objective function for the optimization problem is
\\f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right)\\, where
\\\mathbf{w}=\\w_1, \dots, w_n\\\\ are the estimated weights,
\\\mathbf{s}=\\s_1, \dots, s_n\\\\ are sampling weights (supplied by
`s.weights`), and \\\mathbf{b}=\\b_1, \dots, b_n\\\\ are base weights
(supplied by `b.weights`). The `norm` argument determines \\f(.,.,.)\\,
as detailed below:

- when `norm = "l2"`, \\f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right)
  = \frac{1}{n} \sum_i {s_i(w_i - b_i)^2}\\

- when `norm = "l1"`, \\f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right)
  = \frac{1}{n} \sum_i {s_i \vert w_i - b_i \vert}\\

- when `norm = "linf"`, \\f\left(\mathbf{w},
  \mathbf{b},\mathbf{s}\right) = \max_i {\vert w_i - b_i \vert}\\

- when `norm = "entropy"`, \\f\left(\mathbf{w},
  \mathbf{b},\mathbf{s}\right) = \frac{1}{n} \sum_i {s_i w_i \log
  \frac{w_i}{b_i}}\\

- when `norm = "log"`, \\f\left(\mathbf{w}, \mathbf{b},\mathbf{s}\right)
  = \frac{1}{n} \sum_i {-s_i \log \frac{w_i}{b_i}}\\

By default, `s.weights` and `b.weights` are set to 1 for all units
unless supplied. `b.weights` must be positive when `norm` is `"entropy"`
or `"log"`, and `norm = "linf"` cannot be used when `s.weights` are
supplied.

When `norm = "l2"` and both `s.weights` and `b.weights` are `NULL`,
weights are estimated to maximize the effective sample size. When
`norm = "entropy"`, the estimated weights are equivalent to entropy
balancing weights (Källberg & Waernbaum, 2023). When `norm = "log"`,
`b.weights` are ignored in the optimization, as they do not affect the
estimated weights.

### Dual Variables

Two types of constraints may be associated with each covariate: target
constraints and balance constraints, controlled by `target.tols` and
`tols`, respectively. In the `duals` component of the output, each
covariate has a dual variable for each constraint placed on it. The dual
variable for each constraint is the instantaneous rate of change of the
objective function at the optimum corresponding to a change in the
constraint. Because this relationship is not linear, large changes in
the constraint will not exactly map onto corresponding changes in the
objective function at the optimum, but will be close for small changes
in the constraint. For example, for a covariate with a balance
constraint of .01 and a corresponding dual variable of 40, increasing
(i.e., relaxing) the constraint to .025 will decrease the value of the
objective function at the optimum by approximately \\(.025 - .01) \* 40
= .6\\.

For factor variables, `optweight()` takes the sum of the absolute dual
variables for the constraints for all levels and reports it as the the
single dual variable for the variable itself. This summed dual variable
works the same way as dual variables for continuous variables do.

An addition dual variable is computed for the constraint on the range of
the weights, controlled by `min.w`. A high dual variable for this
constraint implies that decreasing `min.w` will decrease the value of
the objective function at the optimum.

### `solver`

The `solver` argument controls which optimization solver is used.
Different solvers are compatible with each `norm`. See the table below
for allowable options, which package they require, which function does
the solving, and which function controls the settings.

|              |                          |                                                         |                                                                                              |                                                                                                                                                                                               |
|--------------|--------------------------|---------------------------------------------------------|----------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `solver`     | `norm`                   | Package                                                 | Solver function                                                                              | Settings function                                                                                                                                                                             |
| `"osqp"`     | `"l2"`, `"l1"`, `"linf"` | [osqp](https://CRAN.R-project.org/package=osqp)         | [`osqp::solve_osqp()`](https://rdrr.io/pkg/osqp/man/solve_osqp.html)                         | [`osqp::osqpSettings()`](https://rdrr.io/pkg/osqp/man/osqpSettings.html)                                                                                                                      |
| `"highs"`    | `"l2"`, `"l1"`, `"linf"` | [highs](https://CRAN.R-project.org/package=highs)       | [`highs::highs_solve()`](https://rdrr.io/pkg/highs/man/highs_solve.html)                     | [`highs::highs_control()`](https://rdrr.io/pkg/highs/man/highs_control.html) / [`highs::highs_available_solver_options()`](https://rdrr.io/pkg/highs/man/highs_available_solver_options.html) |
| `"lpsolve"`  | `"l1"`, `"linf"`         | [lpSolve](https://CRAN.R-project.org/package=lpSolve)   | [`lpSolve::lp()`](https://rdrr.io/pkg/lpSolve/man/lp.html)                                   | .                                                                                                                                                                                             |
| `"scs"`      | `"entropy"`, `"log"`     | [scs](https://CRAN.R-project.org/package=scs)           | [`scs::scs()`](https://rdrr.io/pkg/scs/man/scs.html)                                         | [`scs::scs_control()`](https://rdrr.io/pkg/scs/man/scs_control.html)                                                                                                                          |
| `"clarabel"` | `"entropy"`, `"log"`     | [clarabel](https://CRAN.R-project.org/package=clarabel) | [`clarabel::clarabel()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.html) | [`clarabel::clarabel_control()`](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_control.html)                                                                                  |

Note that `"lpsolve"` can only be used when `min.w` is nonnegative.

The default `solver` for each `norm` is as follows:

|             |                  |
|-------------|------------------|
| `norm`      | Default `solver` |
| `"l2"`      | `"osqp"`         |
| `"l1"`      | `"highs"`        |
| `"linf"`    | `"highs"`        |
| `"entropy"` | `"scs"`          |
| `"log"`     | `"scs"`          |

If the package corresponding to a default `solver` is not installed but
the package for a different eligible solver is, that will be used.
Otherwise, you will be asked to install the required package. osqp is
required for optweight, and so will be the default for the `"l1"` and
`"linf"` norms if highs is not installed. The default package is the one
has shown good performance for the given norm in informal testing;
generally, all eligible solvers perform about equally well in terms of
accuracy but differ in time taken.

### Solving Convergence Failure

Sometimes the optimization will fail to converge at a solution. There
are a variety of reasons why this might happen, which include that the
constraints are nearly impossible to satisfy or that the optimization
surface is relatively flat. It can be hard to know the exact cause or
how to solve it, but this section offers some solutions one might try.
Typically, solutions can be found most easily when using the `"l2"`
norm; other norms, especially `"linf"` and `"l1"`, are more likely to
see problems.

Rarely is the problem too few iterations, though this is possible. Most
problems can be solved in the default 200,000 iterations, but sometimes
it can help to increase this number with the `max_iter` argument.
Usually, though, this just ends up taking more time without a solution
found.

If the problem is that the constraints are too tight, it can be helpful
to loosen the constraints. Sometimes examining the dual variables of a
solution that has failed to converge can reveal which constraints are
causing the problem. An extreme value of a dual variable typically
suggests that its corresponding constraint is one cause of the failure
to converge.

Sometimes a suboptimal solution is possible; such a solution does not
satisfy the constraints exactly but will come pretty close. To allow
these solutions, the argument `eps` can be increased to larger values.
This is more likely to occur when `s.weights` are supplied.

Sometimes using a different solver can improve performance. Using the
default `solver` for each `norm`, as described above, can reduce the
probability of convergence failures.

## References

Barnard, M., Huling, J. D., & Wolfson, J. (2025). Partially Retargeted
Balancing Weights for Causal Effect Estimation Under Positivity
Violations (No. arXiv:2510.22072). arXiv.
[doi:10.48550/arXiv.2510.22072](https://doi.org/10.48550/arXiv.2510.22072)

Chattopadhyay, A., Cohn, E. R., & Zubizarreta, J. R. (2024). One-Step
Weighting to Generalize and Transport Treatment Effect Estimates to a
Target Population. *The American Statistician*, 78(3), 280–289.
[doi:10.1080/00031305.2023.2267598](https://doi.org/10.1080/00031305.2023.2267598)

de los Angeles Resa, M., & Zubizarreta, J. R. (2020). Direct and Stable
Weight Adjustment in Non-Experimental Studies With Multivalued
Treatments: Analysis of the Effect of an Earthquake on Post-Traumatic
Stress. *Journal of the Royal Statistical Society Series A: Statistics
in Society*, 183(4), 1387–1410.
[doi:10.1111/rssa.12561](https://doi.org/10.1111/rssa.12561)

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

[`optweightMV()`](https://ngreifer.github.io/optweight/reference/optweightMV.md)
for estimating stable balancing weights for multivariate (i.e.,
multiple) treatments simultaneously.

[sbw](https://CRAN.R-project.org/package=sbw), which was the inspiration
for this package and provides some additional functionality for binary
treatments.

[WeightIt](https://CRAN.R-project.org/package=WeightIt), which provides
a simplified interface to `optweight()` and a more efficient
implementation of entropy balancing.

## Examples

``` r
library("cobalt")
#>  cobalt (Version 4.6.1, Build Date: 2025-08-20)
data("lalonde", package = "cobalt")

# Balancing covariates between treatment groups (binary)
(ow1 <- optweight(treat ~ age + educ + married +
                    nodegree + re74, data = lalonde,
                  tols = c(.01, .02, .03, .04, .05),
                  estimand = "ATE"))
#> An optweight object
#>  - number of obs.: 614
#>  - norm minimized: "l2"
#>  - sampling weights: present
#>  - base weights: present
#>  - treatment: 2-category
#>  - estimand: ATE
#>  - covariates: age, educ, married, nodegree, re74
bal.tab(ow1)
#> Balance Measures
#>             Type Diff.Adj
#> age      Contin.    -0.01
#> educ     Contin.     0.02
#> married   Binary    -0.03
#> nodegree  Binary     0.04
#> re74     Contin.    -0.05
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.     185. 
#> Adjusted    415.09   125.3

# Exactly alancing covariates with respect to race (multi-category)
(ow2 <- optweight(race ~ age + educ + married +
                    nodegree + re74, data = lalonde,
                  tols = 0, estimand = "ATT",
                  focal = "black"))
#> An optweight object
#>  - number of obs.: 614
#>  - norm minimized: "l2"
#>  - sampling weights: present
#>  - base weights: present
#>  - treatment: 3-category (black, hispan, white)
#>  - estimand: ATT (focal: black)
#>  - covariates: age, educ, married, nodegree, re74
bal.tab(ow2)
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.            0
#> educ     Contin.            0
#> married   Binary            0
#> nodegree  Binary            0
#> re74     Contin.            0
#> 
#> Effective sample sizes
#>            hispan  white black
#> Unadjusted  72.   299.     243
#> Adjusted    45.96 181.39   243

# Balancing covariates between treatment groups (binary)
# and requesting a specified target population
targets <- process_targets(~ age + educ + married +
                             nodegree + re74,
                           data = lalonde,
                           targets = c(26, 12, .4, .5,
                                       1000))

(ow3a <- optweight(treat ~ age + educ + married +
                     nodegree + re74, data = lalonde,
                   targets = targets,
                   estimand = NULL))
#> An optweight object
#>  - number of obs.: 614
#>  - norm minimized: "l2"
#>  - sampling weights: present
#>  - base weights: present
#>  - treatment: 2-category
#>  - estimand: targets
#>  - covariates: age, educ, married, nodegree, re74

bal.tab(ow3a, disp.means = TRUE)
#> Note: `s.d.denom` not specified; assuming "pooled".
#> Balance Measures
#>             Type M.0.Adj M.1.Adj Diff.Adj
#> age      Contin.    26.0    26.0       -0
#> educ     Contin.    12.0    12.0       -0
#> married   Binary     0.4     0.4       -0
#> nodegree  Binary     0.5     0.5       -0
#> re74     Contin.  1000.0  1000.0       -0
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    158.04   64.09

# Balancing covariates between treatment groups (binary)
# and not requesting a target population
(ow3b <- optweight(treat ~ age + educ + married +
                     nodegree + re74, data = lalonde,
                   targets = NULL,
                   estimand = NULL))
#> An optweight object
#>  - number of obs.: 614
#>  - norm minimized: "l2"
#>  - sampling weights: present
#>  - base weights: present
#>  - treatment: 2-category
#>  - estimand: targets
#>  - covariates: age, educ, married, nodegree, re74

bal.tab(ow3b, disp.means = TRUE)
#> Note: `s.d.denom` not specified; assuming "pooled".
#> Balance Measures
#>             Type   M.0.Adj   M.1.Adj Diff.Adj
#> age      Contin.   26.4160   26.4160       -0
#> educ     Contin.   10.3547   10.3547       -0
#> married   Binary    0.3615    0.3615       -0
#> nodegree  Binary    0.6305    0.6305       -0
#> re74     Contin. 3908.9059 3908.9059       -0
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    382.74  139.23

# Using a different norm
(ow1b <- optweight(treat ~ age + educ + married +
                    nodegree + re74, data = lalonde,
                  tols = c(.01, .02, .03, .04, .05),
                  estimand = "ATE",
                  norm = "l1"))
#> An optweight object
#>  - number of obs.: 614
#>  - norm minimized: "l1"
#>  - sampling weights: present
#>  - base weights: present
#>  - treatment: 2-category
#>  - estimand: ATE
#>  - covariates: age, educ, married, nodegree, re74

summary(ow1b, weight.range = FALSE)
#>                   Summary of weights
#> 
#> - Weight statistics:
#> 
#>            L2    L1     L∞ Rel Ent # Zeros
#> treated 1.637 0.422 16.956   0.567       0
#> control 1.03  0.165 16.058   0.222       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.    185.  
#> Weighted    208.09   50.29
summary(ow1, weight.range = FALSE)
#>                   Summary of weights
#> 
#> - Weight statistics:
#> 
#>            L2    L1    L∞ Rel Ent # Zeros
#> treated 0.69  0.536 3.419   0.198       0
#> control 0.183 0.166 0.424   0.017       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.     185. 
#> Weighted    415.09   125.3

# Allowing for negative weights
ow4 <- optweight(treat ~ age + educ + married + race +
                   nodegree + re74 + re75,
                 data = lalonde,
                 estimand = "ATE",
                 min.w = -Inf)

summary(ow4)
#>                   Summary of weights
#> - Weight ranges:
#> 
#>            Min                                 Max
#> treated -0.987 |---------------------------| 7.254
#> control  0.407     |-----|                   2.17 
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>            137   124    68    23    10
#>  treated 5.193 5.206 6.116 6.205 7.254
#>            388   375   226   196   118
#>  control 2.109  2.11 2.111 2.133  2.17
#> 
#> 
#> - Weight statistics:
#> 
#>            L2    L1    L∞ # Zeros
#> treated 1.608 1.216 6.254       0
#> control 0.499 0.39  1.17        0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.    185.  
#> Weighted    343.49   51.57

# Using `optweight.fit()`
treat <- lalonde$treat
covs <- splitfactor(lalonde[2:8], drop.first = "if2")

ow.fit <- optweight.fit(covs,
                        treat = treat,
                        tols = .02,
                        estimand = "ATE",
                        norm = "l2")
```
