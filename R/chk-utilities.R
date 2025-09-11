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

.err <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::abort(paste(strwrap(m), collapse = "\n"),
               call = pkg_caller_call())
}
.wrn <- function(..., n = NULL, tidy = TRUE, immediate = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)

  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
    rlang::with_options({
      rlang::warn(paste(strwrap(m), collapse = "\n"))
    }, warn = 1)
  }
  else {
    rlang::warn(paste(strwrap(m), collapse = "\n"))
  }
}
.msg <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::inform(paste(strwrap(m), collapse = "\n"), tidy = FALSE)
}
