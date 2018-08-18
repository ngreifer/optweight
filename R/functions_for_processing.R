get.covs.and.treat.from.formula <- function(f, data = NULL, env = .GlobalEnv, ...) {
  A <- list(...)

  tt <- terms(f, data = data)

  #Check if data exists
  if (is_not_null(data) && is.data.frame(data)) {
    data.specified <- TRUE
  }
  else data.specified <- FALSE

  #Check if response exists
  if (is.formula(tt, 2)) {
    resp.vars.mentioned <- as.character(tt)[2]
    resp.vars.failed <- vapply(resp.vars.mentioned, function(v) {
      is_null_or_error(try(eval(parse(text = v), c(data, env)), silent = TRUE))
    }, logical(1L))

    if (any(resp.vars.failed)) {
      if (is_null(A[["treat"]])) stop(paste0("The given response variable, \"", as.character(tt)[2], "\", is not a variable in ", word.list(c("data", "the global environment")[c(data.specified, TRUE)], "or"), "."), call. = FALSE)
      tt <- delete.response(tt)
    }
  }
  else resp.vars.failed <- TRUE

  if (any(!resp.vars.failed)) {
    treat.name <- resp.vars.mentioned[!resp.vars.failed][1]
    tt.treat <- terms(as.formula(paste0(treat.name, " ~ 1")))
    mf.treat <- quote(stats::model.frame(tt.treat, data,
                                         drop.unused.levels = TRUE,
                                         na.action = "na.pass"))

    tryCatch({mf.treat <- eval(mf.treat, c(data, env))},
             error = function(e) {stop(conditionMessage(e), call. = FALSE)})
    treat <- model.response(mf.treat)
  }
  else {
    treat <- A[["treat"]]
    treat.name <- NULL
  }

  #Check if RHS variables exist
  tt.covs <- delete.response(tt)
  rhs.vars.mentioned.lang <- attr(tt.covs, "variables")[-1]
  rhs.vars.mentioned <- vapply(rhs.vars.mentioned.lang, deparse, character(1L))
  rhs.vars.failed <- vapply(rhs.vars.mentioned.lang, function(v) {
    is_null_or_error(try(eval(v, c(data, env)), silent = TRUE))
  }, logical(1L))

  if (any(rhs.vars.failed)) {
    stop(paste0(c("All variables in formula must be variables in data or objects in the global environment.\nMissing variables: ",
                  paste(rhs.vars.mentioned[rhs.vars.failed], collapse=", "))), call. = FALSE)

  }

  rhs.term.labels <- attr(tt.covs, "term.labels")
  rhs.term.orders <- attr(tt.covs, "order")

  rhs.df <- vapply(rhs.vars.mentioned.lang, function(v) {
    is.data.frame(try(eval(v, c(data, env)), silent = TRUE))
  }, logical(1L))

  if (any(rhs.df)) {
    if (any(rhs.vars.mentioned[rhs.df] %in% unlist(sapply(rhs.term.labels[rhs.term.orders > 1], function(x) strsplit(x, ":", fixed = TRUE))))) {
      stop("Interactions with data.frames are not allowed in the input formula.", call. = FALSE)
    }
    addl.dfs <- setNames(lapply(rhs.vars.mentioned.lang[rhs.df], function(x) {eval(x, env)}),
                         rhs.vars.mentioned[rhs.df])

    for (i in rhs.term.labels[rhs.term.labels %in% rhs.vars.mentioned[rhs.df]]) {
      ind <- which(rhs.term.labels == i)
      rhs.term.labels <- append(rhs.term.labels[-ind],
                                values = names(addl.dfs[[i]]),
                                after = ind - 1)
    }
    new.form <- as.formula(paste("~", paste(rhs.term.labels, collapse = " + ")))

    tt.covs <- terms(new.form)
    if (is_not_null(data)) data <- do.call("cbind", unname(c(addl.dfs, list(data))))
    else data <- do.call("cbind", unname(addl.dfs))
  }

  #Get model.frame, report error
  mf.covs <- quote(stats::model.frame(tt.covs, data,
                                      drop.unused.levels = TRUE,
                                      na.action = "na.pass"))
  tryCatch({covs <- eval(mf.covs, c(data, env))},
           error = function(e) {stop(conditionMessage(e), call. = FALSE)})

  if (is_not_null(treat.name) && treat.name %in% names(covs)) stop("The variable on the left side of the formula appears on the right side too.", call. = FALSE)

  if (is_null(rhs.vars.mentioned)) {
    covs <- data.frame(Intercept = rep(1, if (is_null(treat)) 1 else length(treat)))
  }
  else attr(tt.covs, "intercept") <- 0

  #Get full model matrix with interactions too
  covs.matrix <- model.matrix(tt.covs, data = covs,
                              contrasts.arg = lapply(Filter(is.factor, covs),
                                                     contrasts, contrasts=FALSE))
  #attr(covs, "terms") <- NULL

  return(list(reported.covs = covs,
              model.covs = covs.matrix,
              treat = treat,
              treat.name = treat.name))
}
get.treat.type <- function(treat) {
  #Returns treat with treat.type attribute
  nunique.treat <- nunique(treat)
  if (nunique.treat == 2) {
    treat.type <- "binary"
  }
  else if (nunique.treat < 2) {
    stop("The treatment must have at least two unique values.", call. = FALSE)
  }
  else if (is.factor(treat) || is.character(treat)) {
    treat.type <- "multinomial"
    treat <- factor(treat)
  }
  else {
    treat.type <- "continuous"
  }
  attr(treat, "treat.type") <- treat.type
  return(treat)
}
process.focal.and.estimand <- function(focal, estimand, treat, treat.type) {
  if (length(estimand) != 1 || !is.character(estimand)) {
    stop("estimand must be a character vector of length 1.", call. = FALSE)
  }
  estimand <- toupper(estimand)[[1]]

  #Allowable estimands
  AE <- list(binary =  c("ATT", "ATC", "ATE"),
             multinomial = c("ATT", "ATE"))

  if (treat.type != "continuous" && estimand %nin% AE[[treat.type]]) {
    stop(paste0("\"", estimand, "\" is not an allowable estimand with ", treat.type, " treatments. Only ", word.list(AE[[treat.type]], quotes = TRUE, and.or = "and", is.are = TRUE),
                " allowed."), call. = FALSE)
  }

  reported.estimand <- estimand

  #Check focal
  if (treat.type %in% c("binary", "multinomial")) {
    if (estimand == "ATT") {
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
        warning(paste(estimand, "is not compatible with focal. Ignoring focal."), call. = FALSE)
        focal <- NULL
      }
    }
  }

  #Get focal, estimand, and reported estimand
  if (isTRUE(treat.type == "binary")) {
    unique.treat <- unique(treat, nmax = 2)
    unique.treat.bin <- unique(binarize(treat), nmax = 2)
    if (estimand == "ATT") {
      if (is_null(focal)) {
        focal <- unique.treat[unique.treat.bin == 1]
      }
      else if (focal == unique.treat[unique.treat.bin == 0]){
        reported.estimand <- "ATC"
      }
    }
    else if (estimand == "ATC") {
      focal <- unique.treat[unique.treat.bin == 0]
      estimand <- "ATT"
    }
  }
  return(list(focal = focal,
              estimand = estimand,
              reported.estimand = reported.estimand))
}
process.s.weights <- function(s.weights, data = NULL) {
  #Process s.weights
  if (is_not_null(s.weights)) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (is_null(data)) {
        stop("s.weights was specified as a string but there was no argument to data.", call. = FALSE)
      }
      else if (s.weights %in% names(data)) {
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
  }
  return(s.weights)
}
nunique <- function(x, nmax = NA, na.rm = TRUE) {
  if (is_null(x)) return(0)
  else {
    if (na.rm) x <- x[!is.na(x)]
    if (is.factor(x)) return(nlevels(x))
    else return(length(unique(x, nmax = nmax)))
  }

}
nunique.gt <- function(x, n, na.rm = TRUE) {
  if (missing(n)) stop("n must be supplied.")
  if (n < 0) stop("n must be non-negative.")
  if (is_null(x)) FALSE
  else {
    if (na.rm) x <- x[!is.na(x)]
    if (n == 1 && is.numeric(x)) !check_if_zero(max(x) - min(x))
    else if (length(x) < 2000) nunique(x) > n
    else tryCatch(nunique(x, nmax = n) > n, error = function(e) TRUE)
  }
}
is_binary <- function(x) !nunique.gt(x, 2)
all_the_same <- function(x) !nunique.gt(x, 1)
check_if_zero <- function(x) {
  # this is the default tolerance used in all.equal
  tolerance <- .Machine$double.eps^0.5
  # If the absolute deviation between the number and zero is less than
  # the tolerance of the floating point arithmetic, then return TRUE.
  # This means, to me, that I can treat the number as 0 rather than
  # -3.20469e-16 or some such.
  abs(x - 0) < tolerance
}
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
is_null_or_error <- function(x) {is_null(x) || class(x) == "try-error"}
binarize <- function(variable) {
  nas <- is.na(variable)
  if (!is_binary(variable[!nas])) stop(paste0("Cannot binarize ", deparse(substitute(variable)), ": more than two levels."))
  if (is.character(variable)) variable <- factor(variable)
  variable.numeric <- as.numeric(variable)
  if (!is.na(match(0, unique(variable.numeric)))) zero <- 0
  else zero <- min(unique(variable.numeric), na.rm = TRUE)
  newvar <- setNames(ifelse(!nas & variable.numeric==zero, 0, 1), names(variable))
  newvar[nas] <- NA
  return(newvar)
}
col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
  if (is_null(w)) {
    w <- 1
    w.sum <- apply(mat, 2, function(x) sum(!is.na(x)))
  }
  else {
    w.sum <- apply(mat, 2, function(x) sum(w, na.rm = na.rm))
  }
  return(colSums(mat*w, na.rm = na.rm)/w.sum)
}
w.cov.scale <- function(w) {
  (sum(w, na.rm = TRUE)^2 - sum(w^2, na.rm = TRUE)) / sum(w, na.rm = TRUE)
}
col.w.v <- function(mat, w = NULL, na.rm = TRUE) {
  if (is_null(w)) {
    w <- rep(1, nrow(mat))
  }
  return(colSums(t((t(mat) - col.w.m(mat, w, na.rm = na.rm))^2) * w, na.rm = na.rm) / w.cov.scale(w))
}
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
is.formula <- function(f, sides = NULL) {
  res <- is.name(f[[1]])  && deparse(f[[1]]) %in% c( '~', '!') &&
    length(f) >= 2
  if (is_not_null(sides) && is.numeric(sides) && sides %in% c(1,2)) {
    res <- res && length(f) == sides + 1
  }
  return(res)
}
word.list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  if (quotes) word.list <- vapply(word.list, function(x) paste0("\"", x, "\""), character(1L))
  if (L == 0) {
    out <- ""
    attr(out, "plural") = FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") = FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") = FALSE
    }
    else {
      and.or <- match.arg(and.or)
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or," "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "),
                     word.list[L], sep = paste0(", ", and.or," "))

      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") = TRUE
    }


  }
  return(out)
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  nas <- is.na(df)
  if (!is.data.frame(df)) df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  rn <- rownames(df)
  cn <- colnames(df)
  df <- as.data.frame(lapply(df, function(col) {
    if (suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
      as.numeric(as.character(col))
    } else {
      col
    }
  }), stringsAsFactors = FALSE)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  o.negs <- sapply(1:ncol(df), function(x) if (nums[x]) df[[x]] < 0 else rep(FALSE, length(df[[x]])))
  df[nums] <- round(df[nums], digits = digits)
  df[nas] <- ""

  df <- as.data.frame(lapply(df, format, scientific = FALSE, justify = "none"), stringsAsFactors = FALSE)

  for (i in which(nums)) {
    if (any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      lengths <- lengths(s)
      digits.r.of.. <- vapply(seq_along(s), function(x) {
        if (lengths[x] > 1) nchar(s[[x]][lengths[x]])
        else 0 }, numeric(1L))
      df[[i]] <- sapply(seq_along(df[[i]]), function(x) {
        if (df[[i]][x] == "") ""
        else if (lengths[x] <= 1) {
          paste0(c(df[[i]][x], rep(".", pad == 0), rep(pad, max(digits.r.of..) - digits.r.of..[x] + as.numeric(pad != 0))),
                 collapse = "")
        }
        else paste0(c(df[[i]][x], rep(pad, max(digits.r.of..) - digits.r.of..[x])),
                    collapse = "")
      })
    }
  }

  df[o.negs & df == 0] <- paste0("-", df[o.negs & df == 0])

  # Insert NA placeholders
  df[nas] <- na_vals

  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn

  return(df)
}
text.box.plot <- function(range.list, width = 12) {
  full.range <- range(unlist(range.list))
  ratio = diff(full.range)/(width+1)
  rescaled.range.list <- lapply(range.list, function(x) round(x/ratio))
  rescaled.full.range <- round(full.range/ratio)
  d <- as.data.frame(matrix(NA_character_, ncol = 3, nrow = length(range.list),
                            dimnames = list(names(range.list), c("Min", paste(rep(" ", width + 1), collapse = ""), "Max"))),
                     stringsAsFactors = FALSE)
  d[,"Min"] <- vapply(range.list, function(x) x[1], numeric(1L))
  d[,"Max"] <- vapply(range.list, function(x) x[2], numeric(1L))
  for (i in seq_len(nrow(d))) {
    spaces1 <- rescaled.range.list[[i]][1] - rescaled.full.range[1]
    #|
    dashes <- max(0, diff(rescaled.range.list[[i]]) - 2)
    #|
    spaces2 <- max(0, diff(rescaled.full.range) - (spaces1 + 1 + dashes + 1))

    d[i, 2] <- paste0(paste(rep(" ", spaces1), collapse = ""), "|", paste(rep("-", dashes), collapse = ""), "|", paste(rep(" ", spaces2), collapse = ""))
  }
  return(d)
}

#To pass CRAN checks:
utils::globalVariables(c("covs", "dual", "treat"))
