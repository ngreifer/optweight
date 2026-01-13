#chk utilities
pkg_caller_call <- function() {
  pn <- utils::packageName()
  package.funs <- c(getNamespaceExports(pn),
                    .getNamespaceInfo(asNamespace(pn), "S3methods")[, 3L])

  for (i in seq_len(sys.nframe())) {
    e <- sys.call(i)

    n <- rlang::call_name(e)

    if (is_not_null(n) && n %in% package.funs) {
      return(e)
    }
  }

  NULL
}

## Use cli
.err <- function(m, n = NULL, tidy = TRUE, cli = TRUE) {
  if (cli) {
    m <- eval.parent(substitute(cli::format_inline(.m), list(.m = m)))
  }

  chk::message_chk(m, n = n, tidy = tidy) |>
    cli::ansi_strwrap() |>
    paste(collapse = "\n") |>
    rlang::abort(call = pkg_caller_call())
}
.wrn <- function(m, n = NULL, tidy = TRUE, immediate = TRUE, cli = TRUE) {
  if (cli) {
    m <- eval.parent(substitute(cli::format_inline(.m), list(.m = m)))
  }

  m <- chk::message_chk(m, n = n, tidy = tidy)

  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
    rlang::with_options({
      m |>
        cli::ansi_strwrap() |>
        paste(collapse = "\n") |>
        rlang::warn()
    }, warn = 1)
  }
  else {
    m |>
      cli::ansi_strwrap() |>
      paste(collapse = "\n") |>
      rlang::warn()
  }
}
.msg <- function(m, n = NULL, tidy = TRUE, cli = TRUE) {
  if (cli) {
    m <- eval.parent(substitute(cli::format_inline(.m), list(.m = m)))
  }

  chk::message_chk(m, n = n, tidy = tidy) |>
    cli::ansi_strwrap() |>
    paste(collapse = "\n") |>
    rlang::inform(tidy = FALSE)
}

.chk_is <- function(x, class, x_name = NULL) {
  if (!inherits(x, class)) {
    if (is_null(x_name)) {
      x_name <- chk::deparse_backtick_chk(substitute(x))
    }

    .err("{.arg {chk::unbacktick_chk(x_name)}} must inherit from class {.or {add_quotes(class, 1)}}")
  }
}
