# Construct and Check Targets Input

Checks whether proposed target population means values for `targets` are
suitable in number and order for submission to
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md),
[`optweightMV()`](https://ngreifer.github.io/optweight/reference/optweightMV.md),
and
[`optweight.svy()`](https://ngreifer.github.io/optweight/reference/optweight.svy.md),
and returns an object that can supplied to the `targets` argument of
these functions.

## Usage

``` r
process_targets(formula, data = NULL, targets = NULL, s.weights = NULL)

check.targets(...)

# S3 method for class 'optweight.targets'
print(x, digits = 5, ...)
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

- targets:

  a vector of target population mean values for each covariate. These
  should be in the order corresponding to the order of the corresponding
  variable in `formula`, except for interactions, which will appear
  after all lower-order terms. For factor variables, a target value must
  be specified for each level of the factor, and these values must add
  up to 1. If `NULL`, the current sample means will be produced
  (weighted by `s.weights`). If `NA`, an `NA` vector named with the
  covariate names will be produced.

- s.weights:

  a vector of sampling weights. For
  [`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md),
  can also be the name of a variable in `data` that contains sampling
  weights.

- ...:

  for
  [`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md),
  additional arguments passed to
  [`optweight.fit()`](https://ngreifer.github.io/optweight/reference/optweight.md),
  including options that are passed to the settings function
  corresponding to `solver`.

- x:

  an `optweight.targets` object; the output of a call to
  `process_targets()`.

- digits:

  how many digits to print.

## Value

An `optweight.targets` object, which is a named vector of target
population mean values, one for each (expanded) covariate specified in
`formula`. This should be used as an input to the `targets` argument of
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md),
[`optweightMV()`](https://ngreifer.github.io/optweight/reference/optweightMV.md),
and
[`optweight.svy()`](https://ngreifer.github.io/optweight/reference/optweight.svy.md).

## Details

The purpose of `process_targets()` is to allow users to ensure that
their proposed input to `targets` in
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md),
[`optweightMV()`](https://ngreifer.github.io/optweight/reference/optweightMV.md),
and
[`optweight.svy()`](https://ngreifer.github.io/optweight/reference/optweight.svy.md)
is correct both in the number of entries and their order. This is
especially important when factor variables and interactions are included
in the formula because factor variables are split into several dummies
and interactions are moved to the end of the variable list, both of
which can cause some confusion and potential error when entering
`targets` values.

Factor variables are internally split into a dummy variable for each
level, so the user must specify a target population mean value for each
level of the factor. These must add up to 1, and an error will be
displayed if they do not. These values represent the proportion of units
in the target population with each factor level.

Interactions (e.g., `a:b` or `a*b` in the `formula` input) are always
sent to the end of the variable list even if they are specified
elsewhere in the `formula`. It is important to run `process_targets()`
to ensure the order of the proposed `targets` corresponds to the
represented order of covariates used in the formula. You can run
`process_targets(., targets = NA)` to see the order of covariates that
is required without specifying any targets.

## See also

[`process_tols()`](https://ngreifer.github.io/optweight/reference/process_tols.md)

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

# Generating targets; means by default
targets <- process_targets(~ age + race + married +
                             nodegree + re74,
                           data = lalonde)

# Notice race is split into three values
targets
#>  - targets:
#>         age  race_black race_hispan  race_white     married    nodegree 
#>    27.36319     0.39577     0.11726     0.48697     0.41531     0.63029 
#>        re74 
#>  4557.54657 

# Generating targets; NA by default
targets <- process_targets(~ age + race + married +
                             nodegree + re74,
                           data = lalonde,
                           targets = NA)
targets
#>  - variables:
#>  age   race_black   race_hispan   race_white   married   nodegree   re74

# Can also supply just a dataset
covs <- lalonde |>
  subset(select = c(age, race, married,
                    nodegree, re74))

targets <- process_targets(covs)

targets
#>  - targets:
#>         age  race_black race_hispan  race_white     married    nodegree 
#>    27.36319     0.39577     0.11726     0.48697     0.41531     0.63029 
#>        re74 
#>  4557.54657 
```
