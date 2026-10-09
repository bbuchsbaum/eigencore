# Oracle-based differential sweep (helper-oracle.R, docs/test-assurance.md).
#
# Level from EIGENCORE_ORACLE_LEVEL (cran / ci / extended); default "cran"
# under R CMD check on CRAN (40 small cases, in-process, a few seconds) and
# "ci" when NOT_CRAN = "true" (400 cases in callr subprocesses). A failure
# prints, per violating case, a one-line reproducer:
#   source("tests/testthat/helper-oracle.R"); library(eigencore)
#   str(oracle_run_case(<id>))

test_that("oracle sweep: no crash, no certified-but-wrong result, API invariants hold", {
  level <- oracle_level()
  ids <- oracle_level_ids(level)
  subprocess <- level != "cran" && oracle_subprocess_available()
  workers <- switch(level, extended = 3L, ci = 2L, 1L)
  recs <- oracle_run_ids(ids, subprocess = subprocess, workers = workers)
  df <- oracle_record_frame(recs)
  expect_identical(nrow(df), length(ids))
  expect_setequal(df$id, ids)

  hard <- recs[nzchar(df$hard)]
  if (length(hard)) {
    msg <- paste(vapply(hard, oracle_reproducer, ""), collapse = "\n")
    fail(sprintf("%d of %d oracle cases violate a hard invariant (level %s):\n%s",
                 length(hard), length(ids), level, msg))
  } else {
    succeed()
  }

  # Report (never fail on) quality signals: uncertified rate per family.
  ok <- df[df$status == "ok", , drop = FALSE]
  rates <- tapply(!(ok$certified %in% TRUE), ok$family, mean)
  if (!identical(Sys.getenv("EIGENCORE_ORACLE_QUIET"), "true")) {
    message(sprintf(
      "oracle sweep (%s): %d cases, %d solved, %d certified, %d expected errors; uncertified %%: %s",
      level, nrow(df), nrow(ok), sum(ok$certified %in% TRUE),
      sum(df$status == "error_expected"),
      paste(sprintf("%s=%.0f", names(rates), 100 * rates), collapse = " ")))
  }
})

test_that("oracle cases are reproducible from their id alone", {
  a <- oracle_case(17L)
  b <- oracle_case(17L)
  expect_identical(a, b)
  pa <- oracle_build(a)
  pb <- oracle_build(b)
  expect_identical(pa$A, pb$A)
  # building a case never touches the caller's RNG stream
  set.seed(99)
  before <- .Random.seed
  invisible(oracle_build(oracle_case(5L)))
  expect_identical(.Random.seed, before)
})

test_that("the harness detects a certified-but-wrong result", {
  # Self-test: feed oracle_check() a doctored result and make sure every
  # class of hard violation is caught.
  case <- oracle_case(1L)
  case$api <- "eig"
  case$structure <- "herm"
  case$target <- "largest"
  case$k <- 2L
  case$tol <- 1e-8
  case$shim <- NULL
  A <- diag(c(5, 4, 3, 2, 1))
  prob <- list(A = A, V = NULL, jordan = 1L)
  truth <- oracle_truth(case, prob)
  good <- list(values = c(5, 4), vectors = diag(5)[, 1:2],
               certificate = list(passed = TRUE, backward_error = c(0, 0),
                                  max_backward_error = 0,
                                  target_completeness = "probed"))
  expect_length(oracle_check(case, prob, truth, 0, good)$hard, 0L)

  wrong_set <- good
  wrong_set$values <- c(5, 3)
  wrong_set$vectors <- diag(5)[, c(1, 3)]
  expect_match(oracle_check(case, prob, truth, 0, wrong_set)$hard,
               "differs from oracle target set", all = FALSE)
  # without a completeness claim it is a soft finding only
  wrong_set$certificate$target_completeness <- "not_checked"
  chk <- oracle_check(case, prob, truth, 0, wrong_set)
  expect_length(chk$hard, 0L)
  expect_match(chk$soft, "differs from oracle target set", all = FALSE)

  bad_pair <- good
  bad_pair$vectors[, 2] <- (diag(5)[, 2] + diag(5)[, 3]) / sqrt(2)
  bad_pair$certificate$target_completeness <- "not_checked"
  expect_match(oracle_check(case, prob, truth, 0, bad_pair)$hard,
               "true backward error", all = FALSE)

  bad_order <- good
  bad_order$values <- c(4, 5)
  bad_order$vectors <- diag(5)[, c(2, 1)]
  expect_match(oracle_check(case, prob, truth, 0, bad_order)$hard,
               "not ordered", all = FALSE)

  short <- good
  short$values <- 5
  short$vectors <- diag(5)[, 1, drop = FALSE]
  expect_match(oracle_check(case, prob, truth, 0, short)$hard,
               "certified with 1 of 2", all = FALSE)

  under <- good
  under$vectors[, 2] <- (diag(5)[, 2] + 1e-3 * diag(5)[, 3]) / sqrt(1 + 1e-6)
  under$certificate$passed <- FALSE
  under$certificate$backward_error <- c(0, 1e-12)
  expect_match(oracle_check(case, prob, truth, 0, under)$hard,
               "reported backward error", all = FALSE)
})
