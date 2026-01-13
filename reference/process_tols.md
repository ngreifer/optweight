# Construct and Check Tolerance Input

Checks whether proposed tolerance values for `tols` are suitable in
number and order for submission to
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md)
and
[`optweight.svy()`](https://ngreifer.github.io/optweight/reference/optweight.svy.md),
and returns an object that can supplied to the `tols` argument of these
functions.

## Usage

``` r
process_tols(formula, data = NULL, tols = 0)

check.tols(...)

# S3 method for class 'optweight.tols'
print(x, internal = FALSE, digits = 5, ...)
```

## Arguments

- formula:

  a formula with the covariates to be balanced on the right-hand side.
  Interactions and functions of covariates are allowed. Lists of
  formulas are not allowed; multiple formulas must be checked one at a
  time.

- data:

  an optional data set in the form of a data frame that contains the
  variables in `formula`.

- tols:

  a vector of balance tolerance values in standardized mean difference
  units for each covariate. These should be in the order corresponding
  to the order of the corresponding variable in `formula`, except for
  interactions, which will appear after all lower-order terms. If only
  one value is supplied, it will be applied to all covariates.

- ...:

  ignored.

- x:

  an `optweight.tols` object; the output of a call to `process_tols()`.

- internal:

  `logical`; whether to print the tolerance values that are to be used
  internally by
  [`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md).
  See Value section.

- digits:

  how many digits to print.

## Value

An `optweight.tols` object, which is a named vector of tolerance values,
one for each variable specified in `formula`. This should be used as an
input to the `tols` argument of
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md).
The `"internal.tols"` attribute contains the tolerance values to be used
internally by
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md).
These will differ from the vector values when there are factor variables
that are split up; the user only needs to submit one tolerance per
factor variable, but separate tolerance values are produced for each new
dummy created.

## Details

The purpose of `process_tols()` is to allow users to ensure that their
proposed input to `tols` in
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md)
is correct both in the number of entries and their order. This is
especially important when factor variables and interactions are included
in the formula because factor variables are split into several dummies
and interactions are moved to the end of the variable list, both of
which can cause some confusion and potential error when entering `tols`
values.

Factor variables are internally split into a dummy variable for each
level, but the user only needs to specify one tolerance value per
original variable; `process_tols()` automatically expands the `tols`
input to match the newly created variables.

Interactions (e.g., `a:b` or `a*b` in the `formula` input) are always
sent to the end of the variable list even if they are specified
elsewhere in the `formula`. It is important to run `process_tols()` to
ensure the order of the proposed `tols` corresponds to the represented
order of covariates used in
[`optweight()`](https://ngreifer.github.io/optweight/reference/optweight.md).
You can run `process_tols()` with no `tols` input to see the order of
covariates that is required.

Note that only one formula and vector of tolerance values can be
assessed at a time; for multiple treatments, each formula and tolerance
vector must be entered separately.

## See also

[`process_targets()`](https://ngreifer.github.io/optweight/reference/process_targets.md)

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

# Generating tols; 0 by default
tols <- process_tols(treat ~ age + educ + married +
                       nodegree + re74,
                     data = lalonde)

tols
#>  - tols:
#>      age     educ  married nodegree     re74 
#>        0        0        0        0        0 

tols <- process_tols(treat ~ age + educ + married +
                       nodegree + re74,
                     data = lalonde,
                     tols = .05)

tols
#>  - tols:
#>      age     educ  married nodegree     re74 
#>     0.05     0.05     0.05     0.05     0.05 

# Checking the order of interactions; notice they go
# at the end even if specified at the beginning.
tols <- process_tols(treat ~ age:educ + married*race +
                       nodegree + re74,
                     data = lalonde,
                     tols = .05)

tols
#>  - tols:
#>      married         race     nodegree         re74     age:educ married:race 
#>         0.05         0.05         0.05         0.05         0.05         0.05 

# Internal tolerances for expanded covariates
print(tols, internal = TRUE)
#>  - tols:
#>      married         race     nodegree         re74     age:educ married:race 
#>         0.05         0.05         0.05         0.05         0.05         0.05 
#> 
#>  - tols used internally by optweight:
#>             married          race_black         race_hispan          race_white 
#>                0.05                0.05                0.05                0.05 
#>            nodegree                re74            age:educ  married:race_black 
#>                0.05                0.05                0.05                0.05 
#> married:race_hispan  married:race_white 
#>                0.05                0.05 
```
