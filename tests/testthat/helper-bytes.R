# Compare two raw vectors (e.g. serialize() output) and, on mismatch, report
# only the lengths and the first differing offset instead of a full diff,
# which can run to hundreds of thousands of lines and hide other failures.
expect_same_bytes <- function(actual, expected) {
  same <- identical(actual, expected)
  if (!same) {
    n <- min(length(actual), length(expected))
    diff_at <- which(actual[seq_len(n)] != expected[seq_len(n)])
    first <- if (length(diff_at)) diff_at[[1L]] else n + 1L
    msg <- sprintf(
      "bytes differ: lengths %d vs %d, first difference at offset %d (%d differing bytes in common prefix)",
      length(actual), length(expected), first, length(diff_at)
    )
  } else {
    msg <- ""
  }
  expect(same, msg)
  invisible(actual)
}

# When a package is installed with source references kept
# (R_KEEP_PKG_SOURCE=yes, as CI does), every package closure's srcref points
# at a srcfile environment whose bindings are lazy-load promises. Reading
# source text anywhere in the session forces them, which changes serialize()
# output for any object holding such a closure without changing any value.
# Force them up front so byte-identity tests measure real mutations only.
materialize_srcfiles <- function(x) {
  seen <- new.env(parent = emptyenv())
  force_env <- function(e) {
    if (!is.environment(e) || isNamespace(e) ||
        identical(e, globalenv()) || identical(e, baseenv()) ||
        identical(e, emptyenv())) {
      return(invisible(NULL))
    }
    key <- format(e)
    if (exists(key, envir = seen, inherits = FALSE)) {
      return(invisible(NULL))
    }
    assign(key, TRUE, envir = seen)
    for (nm in ls(e, all.names = TRUE)) {
      value <- tryCatch(get(nm, envir = e), error = function(err) NULL)
      walk(value)
    }
    invisible(NULL)
  }
  walk <- function(v) {
    at <- attributes(v)
    for (a in at) {
      if (is.environment(a)) force_env(a) else if (length(attributes(a))) walk(a)
    }
    if (is.function(v)) {
      force_env(environment(v))
    } else if (is.environment(v)) {
      force_env(v)
    } else if (is.list(v)) {
      for (item in v) walk(item)
    }
    invisible(NULL)
  }
  walk(x)
  invisible(x)
}
