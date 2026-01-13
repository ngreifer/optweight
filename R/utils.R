#May use collapse

#Strings
word_list <- function(word.list = NULL, and.or = "and", is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately

  word.list <- setdiff(word.list, c(NA_character_, ""))

  if (is_null(word.list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word.list <- add_quotes(word.list, quotes)

  L <- length(word.list)

  if (L == 1L) {
    out <- word.list
    if (is.are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is_null(and.or) || isFALSE(and.or)) {
    out <- toString(word.list)
  }
  else {
    and.or <- match_arg(and.or, c("and", "or"))

    if (L == 2L) {
      out <- sprintf("%s %s %s",
                     word.list[1L],
                     and.or,
                     word.list[2L])
    }
    else {
      out <- sprintf("%s, %s %s",
                     toString(word.list[-L]),
                     and.or,
                     word.list[L])
    }
  }

  if (is.are) out <- sprintf("%s are", out)

  attr(out, "plural") <- TRUE

  out
}
add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }

  if (isTRUE(quotes)) {
    quotes <- '"'
  }

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, str_rev(quotes)))
  }

  if (!chk::vld_count(quotes) || quotes > 2L) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1L) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }

  x
}
ordinal <- function(x) {
  if (is_null(x) || !is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (length(x) > 1L) {
    out <- setNames(vapply(x, ordinal, character(1L)), names(x))
    return(out)
  }

  x0 <- abs(x)
  out <- paste0(x0, switch(substring(x0, nchar(x0), nchar(x0)),
                           "1" = "st",
                           "2" = "nd",
                           "3" = "rd",
                           "th"))
  if (x < 0) out <- sprintf("-%s", out)

  setNames(out, names(x))
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  if (NROW(df) == 0L || NCOL(df) == 0L) {
    return(df)
  }

  if (!is.data.frame(df)) {
    df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  }

  rn <- rownames(df)
  cn <- colnames(df)

  infs <- o.negs <- array(FALSE, dim = dim(df))
  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1L))
  for (i in which(nums)) {
    infs[, i] <- !nas[, i] & !is.finite(df[[i]])
  }

  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }

  o.negs[, nums] <- !nas[, nums] & df[nums] < 0 & round(df[nums], digits) == 0
  df[nums] <- round(df[nums], digits = digits)

  for (i in which(nums)) {
    df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                      drop0trailing = !identical(as.character(pad), "0"))

    if (!identical(as.character(pad), "0") && any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      s_lengths <- lengths(s)
      digits.r.of.. <- rep.int(0, NROW(df))
      digits.r.of..[s_lengths > 1L] <- nchar(vapply(s[s_lengths > 1L], `[[`, character(1L), 2L))
      max.dig <- max(digits.r.of..)

      dots <- ifelse(s_lengths > 1L, "", if (nzchar(as.character(pad))) "." else pad)
      pads <- vapply(max.dig - digits.r.of.., function(n) strrep(pad, n), character(1L))

      df[[i]] <- paste0(df[[i]], dots, pads)
    }
  }

  df[o.negs] <- paste0("-", df[o.negs])

  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"

  if (is_not_null(rn)) rownames(df) <- rn
  if (is_not_null(cn)) names(df) <- cn

  df
}
text_box_plot <- function(range.list, width = 12L) {
  full.range <- range(unlist(range.list))
  if (all_the_same(full.range)) {
    for (i in seq_along(range.list)) {
      range.list[[i]][1L] <- range.list[[i]][1L] - 1e-6
      range.list[[i]][2L] <- range.list[[i]][2L] + 1e-6
    }
    full.range <- range(unlist(range.list))
  }
  ratio <- diff(full.range) / (width + 1)
  rescaled.range.list <- lapply(range.list, function(x) round(x / ratio))
  rescaled.full.range <- round(full.range / ratio)
  d <- make_df(c("Min", space(width + 1L), "Max"),
               names(range.list),
               "character")
  d[["Min"]] <- vapply(range.list, function(x) x[1L], numeric(1L))
  d[["Max"]] <- vapply(range.list, function(x) x[2L], numeric(1L))
  for (i in seq_row(d)) {
    spaces1 <- rescaled.range.list[[i]][1L] - rescaled.full.range[1L]
    dashes <- max(c(0L, diff(rescaled.range.list[[i]]) - 2L))
    spaces2 <- max(c(0L, diff(rescaled.full.range) - (spaces1 + 1L + dashes + 1L)))

    d[i, 2L] <- sprintf("%s|%s|%s",
                        space(spaces1),
                        strrep("-", dashes),
                        space(spaces2))
  }

  d
}
can_str2num <- function(x) {
  if (is.numeric(x) || is.logical(x)) {
    return(TRUE)
  }

  nas <- is.na(x)
  x_num <- suppressWarnings(as.numeric(as.character(x[!nas])))

  !anyNA(x_num)
}
str2num <- function(x) {
  nas <- is.na(x)
  if (!is.numeric(x) && !is.logical(x)) x <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x))
  is.na(x_num)[nas] <- TRUE
  x_num
}
cat0 <- function(..., file = "", sep = "", fill = FALSE, labels = NULL,
                 append = FALSE) {
  cat(..., file = file, sep = sep, fill = fill, labels = labels,
      append = append)
}
space <- function(n) {
  strrep(" ", n)
}
txtbar <- function(n) {
  strrep("\u2500", n)
}
str_rev <- function(x) {
  strsplit(x, NULL) |>
    lapply(rev) |>
    vapply(paste, character(1L), collapse = "")
}
safe_str2expression <- function(text) {
  expr <- try(str2expression(text), silent = TRUE)

  if (null_or_error(expr)) {
    expr <- str2expression(add_quotes(text, "`"))
  }

  expr
}

#Numbers
check_if_zero <- function(x, tol = sqrt(.Machine$double.eps)) {
  # this is the default tolerance used in all.equal
  abs(x) < tol
}
squish <- function(p, lo = 1e-6, hi = 1 - lo) {
  if (lo > -Inf)
    p[p < lo] <- lo

  if (hi < Inf)
    p[p > hi] <- hi

  p
}

#Statistics
ESS <- function(w) {
  sum(w)^2 / sum(w^2)
}
col.w.v <- function(mat, g = NULL, w = NULL, bin.vars = NULL) {
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any_apply(mat, chk::vld_character_or_factor)) {
        stop("'mat' must be a numeric matrix.")
      }
    }
    else if (is.numeric(mat)) {
      mat <- matrix(mat, ncol = 1L)
    }
    else {
      stop("'mat' must be a numeric matrix.")
    }
  }

  if (is_null(bin.vars)) {
    bin.vars <- alloc(FALSE, ncol(mat))
  }
  else if (length(bin.vars) != ncol(mat) || anyNA(as.logical(bin.vars))) {
    stop("'bin.vars' must be a logical vector with length equal to the number of columns of 'mat'.")
  }

  if (!any(bin.vars)) {
    return(fvar(mat, g = g, w = w))
  }

  .f <- function(z) {z * (1 - z)}

  if (all(bin.vars)) {
    return(.f(fmean(mat, g = g, w = w)))
  }

  vars <- .f(fmean(mat, g = g, w = w))

  if (is.matrix(vars)) {
    vars[,!bin.vars] <- fvar(ss(mat, j = !bin.vars), g = g, w = w)
  }
  else {
    vars[!bin.vars] <- fvar(ss(mat, j = !bin.vars), g = g, w = w)
  }

  vars
}
col.w.sd <- function(mat, g = NULL, w = NULL, bin.vars = NULL) {
  sqrt(col.w.v(mat, g = g, w = w, bin.vars = bin.vars))
}
col.w.cov <- function(mat, y, g = NULL, w = NULL) {
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any_apply(mat, chk::vld_character_or_factor)) {
        stop("'mat' must be a numeric matrix.")
      }
    }
    else if (is.numeric(mat)) {
      mat <- matrix(mat, ncol = 1L)
    }
    else {
      stop("'mat' must be a numeric matrix.")
    }
  }

  if (is_not_null(g)) {
    g <- GRP(g)
  }

  if (is_null(w)) {
    mult <- 1 + 1 / (fnobs(y, g = g) - 1)
  }
  else {
    mult <- 1 / (1 - fsum(w^2, g = g) / fsum(w, g = g)^2)
  }

  fmean(fwithin(mat, g = g, w = w) * fwithin(y, g = g, w = w), g = g, w = w) * mult
}

bw.nrd <- function(x) {
  #R's bw.nrd doesn't always work, but bw.nrd0 does
  bw.nrd0(x) * 1.06 / .9
}

rms_dev <- function(x, bw = NULL, sw = NULL) {
  if (is_null(bw)) {
    bw <- fmean(x)
  }

  if (is_null(sw)) {
    sw <- 1
  }

  sqrt(fmean(sw * (x - bw)^2))
}
mean_abs_dev <- function(x, bw = NULL, sw = NULL) {
  if (is_null(bw)) {
    bw <- mean(x)
  }

  if (is_null(sw)) {
    sw <- 1
  }

  mean(sw * abs(x - bw))
}
max_abs_dev <- function(x, bw = NULL, sw = NULL) {
  if (is_null(bw)) {
    bw <- mean(x)
  }

  if (is_null(sw)) {
    sw <- 1
  }

  max(sw * abs(x - bw))
}
rel_ent <- function(x, bw = NULL, sw = NULL) {
  if (is_null(bw)) {
    bw <- mean(x)
  }

  if (is_null(sw)) {
    sw <- 1
  }

  mean(sw * x * (log(x) - log(bw)))
}

#Formulas
get_varnames <- function(expr) {
  recurse <- function(e) {
    if (is.symbol(e)) {
      # bare variable like age
      return(as.character(e))
    }

    if (!is.call(e)) {
      return(NULL)
    }

    # keep as-is for $, [[, and [
    fn <- e[[1L]]
    if (fn == as.name("$") || fn == as.name("[[") || fn == as.name("[")) {
      return(deparse1(e))
    }

    # strip outer function, recurse into arguments
    lapply(as.list(e)[-1L], recurse) |>
      unlist()
  }

  recurse(expr)
}

#treat/covs
get_covs_and_treat_from_formula2 <- function(f, data = NULL, sep = "", ...) {

  if (!rlang::is_formula(f)) {
    .err("{.arg formula} must be a formula")
  }

  chk::chk_string(sep)

  env <- environment(f)

  #Check if data exists
  if (is_null(data)) {
    data <- env
    data.specified <- FALSE
  }
  else if (is.data.frame(data)) {
    data.specified <- TRUE
  }
  else {
    .wrn("the argument supplied to {.arg data} is not a data frame object. This may causes errors or unexpected results")
    data <- env
    data.specified <- FALSE
  }

  tryCatch({
    tt <- terms(f, data = data)
  },
  error = function(e) {
    msg <- {
      if (identical(conditionMessage(e), "'.' in formula and no 'data' argument"))
        "`.` is not allowed in formulas"
      else
        conditionMessage(e)
    }
    .err(msg, tidy = FALSE, cli = FALSE)
  })

  treat <- ...get("treat")
  treat.name <- NULL
  m <- NULL

  #Check if response exists
  if (rlang::is_formula(tt, lhs = TRUE)) {
    resp.var.mentioned <- attr(tt, "variables")[[2L]]
    resp.var.mentioned.char <- deparse1(resp.var.mentioned)

    test <- tryCatch(eval(resp.var.mentioned, data, env),
                     error = identity)

    if (inherits(test, "simpleError")) {
      m <- conditionMessage(test)
      if (!startsWith(m, "object '") || !endsWith(m, "' not found")) {
        .err(m, tidy = FALSE)
      }

      resp.var.failed <- TRUE
    }
    else {
      resp.var.failed <- is_null(test)
    }

    if (!resp.var.failed) {
      treat.name <- resp.var.mentioned.char
      treat <- test
    }
    else if (is_not_null(treat)) {
      tt <- delete.response(tt)
    }
    else {
      var_name <- utils::strcapture("object '(.*)' not found", m, character(1L))[[1L]]
      env_name <- c("data the supplied dataset"[data.specified], "the global environment")
      .err("the given response variable, {.var {var_name}}, is not a variable in {.or {env_name}}")
    }
  }

  #Check if RHS variables exist
  tt.covs <- delete.response(tt)

  rhs.term.labels <- attr(tt.covs, "term.labels")

  if (is_null(rhs.term.labels)) {
    new.form <- as.formula("~ 0")
    tt.covs <- terms(new.form)

    if (is_null(treat)) {
      covs <- make_df(ncol = 0L, nrow = 1)
      covs.matrix <- model.matrix(tt.covs, data = covs)
    }
    else {
      covs <- make_df(ncol = 0L, nrow = length(treat))
      covs.matrix <- model.matrix(tt.covs, data = covs)

      class(treat) <- unique(c("treat", class(treat)))
      attr(treat, "treat.name") <- treat.name
    }

    return(list(reported.covs = covs,
                model.covs = covs.matrix,
                simple.covs = covs,
                treat = treat))
  }

  rhs.term.orders <- attr(tt.covs, "order")
  rhs.vars.mentioned <- attr(tt.covs, "variables")[-1L]
  rhs.vars.mentioned.char <- vapply(rhs.vars.mentioned, deparse1, character(1L))

  simple.covs <- rhs.vars.mentioned |>
    lapply(get_varnames) |>
    unlist() |>
    unique() |>
    lapply(function(i) {
      iexp <- safe_str2expression(i)

      test <- tryCatch(eval(iexp, data, env),
                       error = identity)

      if (inherits(test, "simpleError")) {
        m <- conditionMessage(test)

        if (!startsWith(m, "object '") || !endsWith(m, "' not found")) {
          .err(m, tidy = FALSE, cli = FALSE)
        }

        return(NULL)
      }

      if (length(dim(test)) == 2L) {
        out <- as.data.frame(test)

        if (is_null(colnames(test))) {
          names(out) <- paste(i, seq_col(out), sep = sep)
        }

        return(as.list(out))
      }

      list(test) |> setNames(i)
    }) |>
    clear_null() |>
    unlist(recursive = FALSE) |>
    list2DF()

  rhs.vars.failed <- rhs.df <- rep_with(FALSE, rhs.vars.mentioned.char)
  addl.dfs <- make_list(rhs.vars.mentioned.char)
  terms.with.interactions <- unlist(lapply(rhs.term.labels[rhs.term.orders > 1], strsplit, ":", fixed = TRUE))

  for (i in seq_along(rhs.vars.mentioned.char)) {
    test <- tryCatch(eval(rhs.vars.mentioned[[i]], data, env),
                     error = identity)

    if (inherits(test, "simpleError")) {
      m <- conditionMessage(test)
      if (!startsWith(m, "object '") || !endsWith(m, "' not found")) {
        .err(m, tidy = FALSE, cli = FALSE)
      }

      rhs.vars.failed[i] <- TRUE
      rhs.vars.mentioned.char[i] <- utils::strcapture("object '(.*)' not found", m, character(1L))[[1L]]
    }

    if (any(rhs.vars.failed)) {
      next
    }

    rhs.vars.failed[i] <- is_null(test)

    if (length(dim(test)) == 2L) {
      rhs.df[i] <- TRUE

      if (rhs.vars.mentioned.char[i] %in% terms.with.interactions) {
        .err("interactions with data frames are not allowed in the input formula")
      }

      if (inherits(test, "rms")) {
        class(test) <- "matrix"
        test <- setNames(as.data.frame(as.matrix(test)),
                         attr(test, "colnames"))
      }
      else if (is_not_null(colnames(test))) {
        colnames(test) <- paste(rhs.vars.mentioned.char[i], colnames(test), sep = sep)
      }
      else {
        colnames(test) <- paste(rhs.vars.mentioned.char[i], seq_col(test), sep = sep)
      }

      addl.dfs[[i]] <- as.data.frame(test)
    }
  }

  if (any(rhs.vars.failed)) {
    .err("all variables in {.arg formula} must be variables in {.arg data} or objects in the global environment.\n
          Missing variables: {.var rhs.vars.mentioned.char[rhs.vars.failed]}")

  }

  rhs.term.labels.list <- setNames(as.list(rhs.term.labels), rhs.term.labels)

  if (any(rhs.df)) {
    for (i in intersect(rhs.term.labels, rhs.vars.mentioned.char[rhs.df])) {
      ind <- match(i, rhs.term.labels)
      rhs.term.labels <- append(rhs.term.labels[-ind],
                                values = names(addl.dfs[[i]]),
                                after = ind - 1L)
      rhs.term.labels.list[[i]] <- names(addl.dfs[[i]])
    }

    data <- {
      if (data.specified) do.call("cbind", unname(c(clear_null(addl.dfs), list(data))))
      else do.call("cbind", unname(clear_null(addl.dfs)))
    }
  }

  new.form <- sprintf("~ %s", paste(vapply(names(rhs.term.labels.list), function(x) {
    if (x %in% rhs.vars.mentioned.char[rhs.df]) paste(add_quotes(rhs.term.labels.list[[x]], "`"), collapse = " + ")
    else rhs.term.labels.list[[x]]
  } , character(1L)), collapse = " + ")) |>
    as.formula()

  tt.covs <- terms(update(new.form,  ~ . - 1))

  #Get model.frame, report error
  mf.covs <- quote(stats::model.frame(tt.covs, data,
                                      drop.unused.levels = TRUE,
                                      na.action = "na.pass"))

  covs <- tryCatch(eval(mf.covs),
                   error = function(e) {
                     .err(conditionMessage(e), tidy = FALSE, cli = FALSE)
                   })

  if (is_not_null(treat.name) && utils::hasName(covs, treat.name)) {
    .err("the variable on the left side of the formula appears on the right side too")
  }

  s <- nzchar(sep)

  if (s) {
    original.covs.levels <- make_list(names(covs))
  }

  for (i in names(covs)) {
    if (is.character(covs[[i]])) {
      covs[[i]] <- factor(covs[[i]])
    }
    else if (!is.factor(covs[[i]])) {
      next
    }

    if (length(unique(covs[[i]])) == 1L) {
      covs[[i]] <- 1
    }
    else if (s) {
      original.covs.levels[[i]] <- levels(covs[[i]])
      levels(covs[[i]]) <- paste0(sep, original.covs.levels[[i]])
    }
  }

  #Get full model matrix with interactions too
  covs.matrix <- model.matrix(tt.covs, data = covs,
                              contrasts.arg = lapply(Filter(is.factor, covs),
                                                     contrasts, contrasts = FALSE))

  if (s) {
    for (i in names(covs)[vapply(covs, is.factor, logical(1L))]) {
      levels(covs[[i]]) <- original.covs.levels[[i]]
    }
  }

  if (is_not_null(treat)) {
    class(treat) <- unique(c("treat", class(treat)))
    attr(treat, "treat.name") <- treat.name
  }

  list(reported.covs = covs,
       model.covs = covs.matrix,
       simple.covs = simple.covs,
       treat = treat)
}
assign_treat_type <- function(treat, use.multi = FALSE) {
  #Returns treat with treat.type attribute
  nunique.treat <- fnunique(treat)

  if (nunique.treat < 2L) {
    .err("the treatment must have at least two unique values")
  }

  if (!use.multi && nunique.treat == 2L) {
    treat.type <- "binary"
  }
  else if (use.multi || chk::vld_character_or_factor(treat)) {
    treat.type <- "multi-category"
    if (!inherits(treat, "processed.treat")) treat <- factor(treat)
  }
  else {
    treat.type <- "continuous"
  }

  attr(treat, "treat.type") <- treat.type

  treat
}
get_treat_type <- function(treat) {
  attr(treat, "treat.type")
}
has_treat_type <- function(treat) {
  is_not_null(get_treat_type(treat))
}

#Uniqueness
all_the_same <- function(x, na.rm = TRUE) {
  if (anyNA(x)) {
    x <- na_rm(x)

    if (!na.rm) {
      return(is_null(x))
    }
  }

  if (length(x) == 1L) {
    return(TRUE)
  }

  if (is.numeric(x)) {
    return(check_if_zero(diff1(.range(x))))
  }

  allv(x, x[1L])
}
is_binary <- function(x, na.rm = TRUE) {
  if (na.rm) x <- na_rm(x)

  !all_the_same(x) && all_the_same(x[x != x[1L]])
}
is_binary_col <- function(dat, na.rm = TRUE) {
  if (length(dim(dat)) != 2L) {
    stop("is_binary_col() cannot be used with objects that don't have 2 dimensions.")
  }

  apply(dat, 2L, is_binary)
}

#R Processing
make_list <- function(n) {
  if (length(n) == 1L && is.numeric(n)) {
    vector("list", as.integer(n))
  }
  else if (is_not_null(n) && is.atomic(n)) {
    setNames(vector("list", length(n)), as.character(n))
  }
  else {
    stop("'n' must be an integer(ish) scalar or an atomic variable.")
  }
}
make_df <- function(ncol, nrow = 0L, types = "numeric") {
  if (is_null(ncol)) {
    ncol <- 0L
  }

  if (length(ncol) == 1L && is.numeric(ncol)) {
    col_names <- NULL
    ncol <- as.integer(ncol)
  }
  else if (is.atomic(ncol)) {
    col_names <- as.character(ncol)
    ncol <- length(ncol)
  }

  if (is_null(nrow)) nrow <- 0L

  if (length(nrow) == 1L && is.numeric(nrow)) {
    row_names <- NULL
    nrow <- as.integer(nrow)
  }
  else if (is.atomic(nrow)) {
    row_names <- as.character(nrow)
    nrow <- length(nrow)
  }

  df <- as.data.frame.matrix(matrix(NA_real_, nrow = nrow, ncol = ncol))

  names(df) <- col_names
  rownames(df) <- row_names

  if (is_null(types)) {
    return(df)
  }

  if (length(types) %nin% c(1L, ncol)) {
    stop("'types' must be equal to the number of columns.")
  }

  if (!is.character(types) ||
      !all(types %in% c("numeric", "integer", "logical", "character", NA))) {
    stop("'types' must be an acceptable type. For factors, use NA.")
  }

  if (length(types) == 1L) {
    types <- rep.int(types, ncol)
  }

  for (i in which(!is.na(types))) {
    if (types[i] != "numeric") {
      df[[i]] <- get(types[i])(nrow)
    }
  }

  df
}
rep_with <- function(x, y) {
  #Helper function to fill named vectors with x and given names of y
  rep.int(x, length(y)) |>
    setNames(names(y))
}
is_null <- function(x) {length(x) == 0L}
is_not_null <- function(x) {!is_null(x)}
clear_null <- function(x) {
  x[lengths(x) == 0L] <- NULL
  x
}
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
`%or%` <- function(x, y) {
  # like `%||%` but works for non-NULL length 0 objects
  if (is_null(x)) y else x
}
null_or_error <- function(x) {is_null(x) || inherits(x, "try-error")}
.attr <- function(x, which, exact = TRUE) {
  attr(x, which, exact = exact)
}

# Needed to be able to use collapse::ss() on Matrix objects
ss <- function(x, i, j, check = TRUE) {
  if (inherits(x, "Matrix")) {
    return(if (missing(j)) x[i, , drop = FALSE]
           else if (missing(i)) x[, j, drop = FALSE]
           else x[i, j, drop = FALSE])
  }

  collapse::ss(x, i, j, check)
}

#More informative and cleaner version of base::match.arg(). Uses chk and cli.
match_arg <- function(arg, choices, several.ok = FALSE, context = NULL) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg)) {
    .err("no argument was supplied to {.fn match_arg} (this is a bug)")
  }

  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)

  if (missing(choices)) {
    sysP <- sys.parent()
    formal.args <- formals(sys.function(sysP))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is_null(arg)) {
    return(choices[1L])
  }

  if (several.ok) {
    chk::chk_character(arg, x_name = add_quotes(arg.name, "`"))
  }
  else {
    chk::chk_string(arg, x_name = add_quotes(arg.name, "`"))

    if (identical(arg, choices)) {
      return(arg[1L])
    }
  }

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)

  if (allv(i, 0L)) {
    one_of <- {
      if (length(choices) <= 1L) NULL
      else if (several.ok) "at least one of"
      else "one of"
    }

    .err("{(context)} the argument to {.arg {arg.name}} should be {one_of} {.or {add_quotes(choices)}}")
  }

  i <- i[i > 0L]

  choices[i]
}

grab <- function(x, what) {
  lapply(x, function(z) z[[what]])
}
cbind_distinct <- function(x) {
  nm <- lapply(x, colnames)

  for (i in seq_along(x)[-1L]) {
    x[[i]] <- x[[i]][, setdiff(colnames(x[[i]]), unlist(nm[seq_len(i)])), drop = FALSE]
  }

  do.call("cbind", clear_null(x))
}
diff1 <- function(x) {
  x[-1L] - x[-length(x)]
}

#Extract variables from ..., similar to ...elt(), by name without evaluating list(...)
...get <- function(x, ifnotfound = NULL) {
  expr <- quote({
    .m1 <- match(.x, ...names())
    if (anyNA(.m1)) {
      .ifnotfound
    }
    else {
      ...elt(.m1[1L]) %or% .ifnotfound
    }
  })

  eval(expr,
       pairlist(.x = x[1L], .ifnotfound = ifnotfound),
       parent.frame(1L))
}
...mget <- function(x) {
  found <- match(x, eval(quote(...names()), parent.frame(1L)))

  not_found <- is.na(found)

  if (all(not_found)) {
    return(list())
  }

  lapply(found[!not_found], function(z) {
    eval(quote(...elt(.z)),
         pairlist(.z = z),
         parent.frame(3L))
  }) |>
    setNames(x[!not_found])
}

any_apply <- function(X, FUN, ...) {
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }

  for (x in X) {
    if (isTRUE(FUN(x, ...))) {
      return(TRUE)
    }
  }

  FALSE
}
all_apply <- function(X, FUN, ...) {
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }

  for (x in X) {
    if (isFALSE(FUN(x, ...))) {
      return(FALSE)
    }
  }

  TRUE
}

#-------#cli utilities-------
.it <- function(...) cli::style_italic(...)
.ul <- function(...) cli::style_underline(...)
.st <- function(...) cli::style_strikethrough(...)
