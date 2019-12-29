process.focal.and.estimand <- function(focal, estimand, targets, treat, treat.type) {
  if ((is_null(targets) || all(is.na(targets))) && is_not_null(estimand)) {
    if (!(length(estimand) == 1 && is.character(estimand))) {
      stop("estimand must be a character vector of length 1.", call. = FALSE)
    }
    estimand_ <- toupper(estimand)[[1]]

    #Allowable estimands
    AE <- list(binary =  c("ATT", "ATC", "ATE"),
               multinomial = c("ATT", "ATE"),
               continuous = "ATE")

    if (estimand_ %nin% AE[[treat.type]]) {
      stop(paste0("\"", estimand, "\" is not an allowable estimand with ", treat.type, " treatments. Only ", word_list(AE[[treat.type]], quotes = TRUE, and.or = "and", is.are = TRUE),
                  " allowed."), call. = FALSE)
    }

    reported.estimand <- estimand_
  }
  else {
    if (is_not_null(estimand)) warning("targets are not NULL; ignoring estimand.", call. = FALSE, immediate. = TRUE)
    estimand <- NULL
    reported.estimand <- "targets"
    estimand_ <- NULL
  }

  #Check focal
  if (treat.type %in% c("binary", "multinomial")) {
    if (is_null(estimand)) { #Targets were supplied
      if (is_not_null(focal)) {
        warning(paste("Only estimand = \"ATT\" is compatible with focal. Ignoring focal."), call. = FALSE)
        focal <- NULL
      }
    }
    else if (estimand_ == "ATT") {
      if (is_null(focal)) {
        if (treat.type == "multinomial") {
          stop("When estimand = \"ATT\" for multinomial treatments, an argument must be supplied to focal.", call. = FALSE)
        }
      }
      else if (length(focal) > 1L || !is.atomic(focal) || !any(unique(treat) == focal)) {
        stop("The argument supplied to focal must be the name of a level of treat.", call. = FALSE)
      }
    }
    else {
      if (is_not_null(focal)) {
        warning(paste(estimand_, "is not compatible with focal. Setting estimand to \"ATT\"."), call. = FALSE)
        estimand_ <- "ATT"
      }
    }
  }

  #Get focal, estimand, and reported estimand
  if (isTRUE(treat.type == "binary")) {
    unique.treat <- unique(treat, nmax = 2)
    unique.treat.bin <- unique(binarize(treat), nmax = 2)
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
  return(list(focal = focal,
              estimand = estimand_,
              reported.estimand = reported.estimand))
}
process.b.weights <- function(b.weights, data = NULL) {
  #Process b.weights
  if (is_not_null(b.weights)) {
    if (!(is.character(b.weights) && length(b.weights) == 1) && !is.numeric(b.weights)) {
      stop("The argument to b.weights must be a vector or data frame of base weights or the (quoted) name of the variable in data that contains base weights.", call. = FALSE)
    }
    if (is.character(b.weights) && length(b.weights)==1) {
      if (is_null(data)) {
        stop("b.weights was specified as a string but there was no argument to data.", call. = FALSE)
      }
      else if (b.weights %in% names(data)) {
        b.weights <- data[[b.weights]]
      }
      else stop("The name supplied to b.weights is not the name of a variable in data.", call. = FALSE)
    }
  }
  else b.weights <- NULL
  return(b.weights)
}

#To pass CRAN checks:
utils::globalVariables(c("covs", "dual", "treat", "constraint"))
