process.focal.and.estimand <- function(focal, estimand, targets, treat, treat.type) {
  if ((is_null(targets) || all(is.na(targets))) && is_not_null(estimand)) {
    chk::chk_string(estimand)
    estimand_ <- toupper(estimand)

    #Allowable estimands
    AE <- list(binary =  c("ATT", "ATC", "ATE"),
               `multi-category` = c("ATT", "ATE"),
               continuous = "ATE")

    if (estimand_ %nin% AE[[treat.type]]) {
      .err(sprintf("%s is not an allowable estimand with %s treatments. Only %s allowed",
                   add_quotes(estimand),
                   treat.type,
                   word_list(AE[[treat.type]], quotes = TRUE, and.or = "and", is.are = TRUE)))
    }

    reported.estimand <- estimand_
  }
  else {
    if (is_not_null(estimand)) {
      .wrn("`targets` are not NULL; ignoring `estimand`")
    }

    estimand <- NULL
    reported.estimand <- "targets"
    estimand_ <- NULL
  }

  #Check focal
  if (treat.type %in% c("binary", "multi-category")) {
    if (is_null(estimand)) { #Targets were supplied
      if (is_not_null(focal)) {
        .wrn("only `estimand = \"ATT\"` is compatible with `focal`. Ignoring `focal`")
        focal <- NULL
      }
    }
    else if (estimand_ == "ATT") {
      if (is_null(focal)) {
        if (treat.type == "multi-category") {
          .err("when `estimand = \"ATT\"` for multi-category treatments, an argument must be supplied to `focal`")
        }
      }
      else if (length(focal) > 1L || !is.atomic(focal) || focal %nin% treat) {
        .err("the argument supplied to `focal` must be the name of a level of treatment")
      }
    }
    else if (is_not_null(focal)) {
      .wrn(sprintf("`estimand = %s` is not compatible with `focal`. Setting `estimand` to \"ATT\"",
                   add_quotes(estimand_)))
      estimand_ <- "ATT"
    }

    #Get focal, estimand, and reported estimand
    if (identical(treat.type, "binary")) {
      unique.treat <- unique(treat, nmax = 2L)
      unique.treat.bin <- unique(binarize(treat), nmax = 2L)

      if (is_not_null(estimand)) {
        if (estimand_ == "ATT") {
          if (is_null(focal)) {
            focal <- unique.treat[unique.treat.bin == 1]
          }
          else if (focal == unique.treat[unique.treat.bin == 0]){
            reported.estimand <- "ATC"
          }
        }
        else if (estimand_ == "ATC") {
          focal <- unique.treat[unique.treat.bin == 0]
          estimand_ <- "ATT"
        }
      }
    }
  }

  list(focal = focal,
       estimand = estimand_,
       reported.estimand = reported.estimand)
}

process.s.weights <- function(s.weights, data = NULL) {
  #Process s.weights
  if (is_null(s.weights)) {
    return(NULL)
  }

  if (chk::vld_numeric(s.weights)) {
    return(as.numeric(s.weights))
  }

  if (!chk::vld_string(s.weights)) {
    .err("the argument to `s.weights` must be a vector or data frame of sampling weights or the (quoted) name of the variable in `data` that contains sampling weights")
  }

  if (is_null(data)) {
    .err("`s.weights` was specified as a string but there was no argument to `data`")
  }

  if (s.weights %nin% names(data)) {
    .err("the name supplied to `s.weights` is not the name of a variable in `data`")
  }

  as.numeric(data[[s.weights]])
}
process.b.weights <- function(b.weights, data = NULL) {
  #Process b.weights
  if (is_null(b.weights)) {
    return(NULL)
  }

  if (chk::vld_numeric(b.weights)) {
    return(as.numeric(b.weights))
  }

  if (!chk::vld_string(b.weights)) {
    .err("the argument to `b.weights` must be a vector or data frame of base weights or the (quoted) name of the variable in `data` that contains base weights")
  }

  if (is_null(data)) {
    .err("`b.weights` was specified as a string but there was no argument to `data`")
  }

  if (b.weights %nin% names(data)) {
    .err("the name supplied to `b.weights` is not the name of a variable in `data`")
  }

  as.numeric(data[[b.weights]])
}

process.bin.vars <- function(bin.vars, mat) {
  if (missing(bin.vars)) {
    return(apply(mat, 2L, is_0_1))
  }

  if (is_null(bin.vars)) {
    return(rep.int(FALSE, ncol(mat)))
  }

  if (is.logical(bin.vars)) {

    if (length(bin.vars) != ncol(mat)) {
      .err("if `bin.vars` is logical, it must have length equal to the number of columns of `mat`")
    }

    bin.vars[is.na(bin.vars)] <- FALSE
  }
  else if (is.numeric(bin.vars)) {
    bin.vars <- bin.vars[!is.na(bin.vars) & bin.vars != 0]

    if (any(bin.vars < 0) && any(bin.vars > 0)) {
      .err("positive and negative indices cannot be mixed with `bin.vars`")
    }

    if (any(abs(bin.vars) > ncol(mat))) {
      .err("if `bin.vars` is numeric, none of its values can exceed the number of columns of `mat`")
    }

    logical.bin.vars <- rep.int(any(bin.vars < 0), ncol(mat))
    logical.bin.vars[abs(bin.vars)] <- !logical.bin.vars[abs(bin.vars)]
    bin.vars <- logical.bin.vars
  }
  else if (is.character(bin.vars)) {
    bin.vars <- bin.vars[!is.na(bin.vars) & nzchar(bin.vars)]

    if (is_null(colnames(mat))) {
      .err("if `bin.vars` is character, `mat` must have column names")
    }

    if (!all(bin.vars %in% colnames(mat))) {
      .err("if `bin.vars` is character, all its values must be column names of `mat`")
    }

    bin.vars <- colnames(mat) %in% bin.vars
  }
  else {
    .err("`bin.vars` must be a logical, numeric, or character vector")
  }

  bin.vars
}

check_missing_covs <- function(covs) {
  k <- ncol(covs)
  for (i in seq_len(k)) {
    if (anyNA(covs[[i]]) || (is.numeric(covs[[i]]) && !all(is.finite(covs[[i]])))) {
      covariates.with.missingness <- names(covs)[i:k][vapply(i:k, function(j) anyNA(covs[[j]]) ||
                                                               (is.numeric(covs[[j]]) && !all(is.finite(covs[[j]]))),
                                                             logical(1L))]
      .err(sprintf("Missing and non-finite values are not allowed in the covariates. Covariates with missingness or non-finite values:\n\t%s",
                  toString(covariates.with.missingness)), tidy = FALSE)
    }
  }
}


#To pass CRAN checks:
# utils::globalVariables(c("covs", "dual", "treat", "constraint"))
