# Tranche 3 workstream B: C10 (C++ cleanup before R errors), P16 (user
# interrupts in native loops) and the C28 remainder (CSC structure validation
# at the block Lanczos / LOBPCG / scalar Krylov / Arnoldi entry points).

unwind_selftest <- function(mode) {
  .Call("eigencore_unwind_selftest", mode, PACKAGE = "eigencore")
}

test_that("C10: every failure kind becomes an R condition after C++ cleanup", {
  expect_identical(unwind_selftest("live"), 0L)
  expect_identical(unwind_selftest("ok"), 1L)  # one tracked object while running
  expect_identical(unwind_selftest("live"), 0L)

  expect_error(unwind_selftest("error"),
               "eigencore selftest: error with 1 live C\\+\\+ object\\(s\\)")
  expect_identical(unwind_selftest("live"), 0L)

  expect_error(unwind_selftest("bad_alloc"),
               "memory allocation failed (std::bad_alloc)", fixed = TRUE)
  expect_identical(unwind_selftest("live"), 0L)

  expect_error(unwind_selftest("length_error"))
  expect_identical(unwind_selftest("live"), 0L)

  # An unsatisfiable std::vector size is an R error, not std::terminate.
  expect_error(unwind_selftest("huge_vector"))
  expect_identical(unwind_selftest("live"), 0L)

  # An R error raised inside unwind-protected R API code (a longjmp) becomes a
  # C++ exception, the C++ objects are destroyed, then R's unwind continues.
  expect_error(unwind_selftest("r_stop"), "eigencore selftest: R-level stop")
  expect_identical(unwind_selftest("live"), 0L)

  # The same through the remapped allocVector: R refuses an over-long vector.
  # The wording differs across R versions ("vector is too large" in older R,
  # "cannot allocate vector of length ..." in R 4.6).
  expect_error(unwind_selftest("r_alloc_failure"),
               "vector is too large|cannot allocate vector")
  expect_identical(unwind_selftest("live"), 0L)

  expect_error(unwind_selftest("no-such-mode"), "unknown selftest mode")
  expect_identical(unwind_selftest("live"), 0L)
})

test_that("P16: interrupt plumbing signals an R interrupt condition after cleanup", {
  # No interrupt is pending, so the R_ToplevelExec probe returns normally.
  expect_true(unwind_selftest("probe"))
  expect_identical(unwind_selftest("live"), 0L)

  caught <- tryCatch(unwind_selftest("interrupt"),
                     interrupt = function(cond) cond)
  expect_s3_class(caught, "interrupt")
  expect_s3_class(caught, "condition")
  expect_identical(unwind_selftest("live"), 0L)

  # Without an interrupt handler the condition falls through to an R error.
  expect_error(unwind_selftest("interrupt"), "interrupted by user")
  expect_identical(unwind_selftest("live"), 0L)
})

rss_mb <- function() {
  status <- tryCatch(readLines("/proc/self/status"), error = function(e) NULL,
                     warning = function(w) NULL)
  line <- grep("^VmRSS:", status, value = TRUE)
  if (length(line) != 1L) {
    return(NA_real_)
  }
  as.numeric(gsub("[^0-9]", "", line)) / 1024
}

test_that("C10: the dense shift-invert singular-sigma error path does not leak", {
  # AddressSanitizer keeps freed blocks in a quarantine (256 MB by default),
  # so RSS grows by the freed factors; LeakSanitizer checks leaks there.
  skip_if(running_under_asan(), "RSS growth is not meaningful under ASan")
  n <- 300L
  A <- diag(as.numeric(seq_len(n)))
  sigma <- 5  # exactly an eigenvalue: A - sigma I is singular
  call_si <- function() {
    .Call("eigencore_shift_invert_lanczos_dense", A, sigma, 20L, rep(1, n),
          2L, 0L, 1e-8, PACKAGE = "eigencore")
  }
  expect_error(call_si(), "perturb sigma")
  for (rep in 1:20) try(call_si(), silent = TRUE)  # warm up the allocator
  invisible(gc())
  before <- rss_mb()
  for (rep in 1:200) {
    expect_error(call_si(), "perturb sigma")
  }
  invisible(gc())
  after <- rss_mb()
  # Each failed call used to leak the n x n LDL^T factor (0.7 MB): 200 calls
  # leaked ~140 MB. Allow generous allocator noise.
  if (is.finite(before) && is.finite(after)) {
    expect_lt(after - before, 40)
  }
  expect_identical(unwind_selftest("live"), 0L)
})

test_that("C10: absurd native sizes give an R error instead of terminating R", {
  A <- crossprod(matrix(stats::rnorm(400), 20L))
  expect_error(
    .Call("eigencore_lanczos_dense", A, .Machine$integer.max, rep(1, 20L),
          2L, 0L, 1e-8, PACKAGE = "eigencore")
  )
  expect_error(unwind_selftest("huge_vector"))
  # R is still fully functional afterwards.
  expect_equal(sum(1:10), 55L)
})

test_that("C28: block Lanczos, LOBPCG, scalar Krylov and Arnoldi CSC entries reject malformed slots", {
  skip_if_not_installed("Matrix")
  set.seed(28)
  S <- Matrix::rsparsematrix(12L, 12L, density = 0.4)
  S <- as(Matrix::forceSymmetric(S + Matrix::t(S)), "generalMatrix")
  i <- S@i
  p <- S@p
  x <- S@x
  dim <- S@Dim
  n <- 12L
  start1 <- rep(1, n)
  start2 <- matrix(1, n, 2L)
  bdiag <- rep(1, n)

  # Each closure passes the (possibly malformed) slots as the CSC operand and
  # well-typed placeholders elsewhere; validation runs before other arguments
  # are interpreted.
  entries <- list(
    block_lanczos_csc = function(i, p, x) .Call(
      "eigencore_block_lanczos_csc", i, p, x, dim, 2L, 6L, 2L, 0L, 1e-8, start2,
      PACKAGE = "eigencore"),
    block_thick_restart_lanczos_csc = function(i, p, x) .Call(
      "eigencore_block_thick_restart_lanczos_csc", i, p, x, dim, 2L, 6L, 2L, 0L,
      1e-8, 10L, 1, start2, 0L, PACKAGE = "eigencore"),
    normal_thick_restart_lanczos_csc = function(i, p, x) .Call(
      "eigencore_normal_thick_restart_lanczos_csc", i, p, x, dim, 0L, 2L, 6L, 2L,
      0L, 1e-8, 10L, 1, start2, PACKAGE = "eigencore"),
    lobpcg_dense_csc_b = function(i, p, x) .Call(
      "eigencore_lobpcg_dense_csc_b", as.matrix(S), i, p, x, dim, 2L, 10L, 0L,
      1e-8, start2, numeric(), numeric(), numeric(), NULL, PACKAGE = "eigencore"),
    lobpcg_csc_diagonal_b = function(i, p, x) .Call(
      "eigencore_lobpcg_csc_diagonal_b", i, p, x, dim, bdiag, FALSE, 2L, 10L, 0L,
      1e-8, start2, numeric(), numeric(), numeric(), NULL, PACKAGE = "eigencore"),
    lobpcg_csc_csc_b_A = function(i, p, x) .Call(
      "eigencore_lobpcg_csc_csc_b", i, p, x, dim, S@i, S@p, S@x, dim, 2L, 10L, 0L,
      1e-8, start2, numeric(), numeric(), numeric(), NULL, PACKAGE = "eigencore"),
    lobpcg_csc_csc_b_B = function(i, p, x) .Call(
      "eigencore_lobpcg_csc_csc_b", S@i, S@p, S@x, dim, i, p, x, dim, 2L, 10L, 0L,
      1e-8, start2, numeric(), numeric(), numeric(), NULL, PACKAGE = "eigencore"),
    lobpcg_csc_operator_b = function(i, p, x) .Call(
      "eigencore_lobpcg_csc_operator_b", i, p, x, dim, function(X) X, 2L, 10L, 0L,
      1e-8, start2, numeric(), numeric(), numeric(), NULL, PACKAGE = "eigencore"),
    lobpcg_csc = function(i, p, x) .Call(
      "eigencore_lobpcg_csc", i, p, x, dim, 2L, 10L, 0L, 1e-8, start2,
      numeric(), numeric(), numeric(), NULL, PACKAGE = "eigencore"),
    lobpcg_csc_shifted_tridiagonal = function(i, p, x) .Call(
      "eigencore_lobpcg_csc_shifted_tridiagonal", i, p, x, dim, 2L, 10L, 0L, 1e-8,
      start2, 0, PACKAGE = "eigencore"),
    lanczos_csc = function(i, p, x) .Call(
      "eigencore_lanczos_csc", i, p, x, dim, 6L, start1, 2L, 0L, 1e-8,
      PACKAGE = "eigencore"),
    golub_kahan_csc = function(i, p, x) .Call(
      "eigencore_golub_kahan_csc", i, p, x, dim, 6L, start1, 2L, 0L, 1e-8, FALSE,
      PACKAGE = "eigencore"),
    golub_kahan_centered_scaled_csc = function(i, p, x) .Call(
      "eigencore_golub_kahan_centered_scaled_csc", i, p, x, dim, rep(0, n),
      rep(1, n), 6L, start1, 2L, 0L, 1e-8, FALSE, TRUE, TRUE,
      PACKAGE = "eigencore"),
    arnoldi_csc_cycle = function(i, p, x) .Call(
      "eigencore_arnoldi_csc_cycle", i, p, x, dim, start1, 6L,
      PACKAGE = "eigencore"),
    arnoldi_ks_csc = function(i, p, x) .Call(
      "eigencore_arnoldi_ks_csc", i, p, x, dim, FALSE, start1, 2L, 8L, "LM",
      numeric(), 1e-8, 10L, PACKAGE = "eigencore")
  )

  bad_p_len <- p[-1L]
  bad_p_mono <- p
  bad_p_mono[3L] <- bad_p_mono[5L] + 1L
  bad_p_end <- p
  bad_p_end[length(p)] <- length(x) + 5L
  bad_i_high <- i
  bad_i_high[1L] <- n
  bad_i_neg <- i
  bad_i_neg[1L] <- -1L
  bad_i_na <- i
  bad_i_na[2L] <- NA_integer_

  for (name in names(entries)) {
    fn <- entries[[name]]
    info <- paste("entry:", name)
    expect_error(fn(i, bad_p_len, x), "invalid CSC structure", info = info)
    expect_error(fn(i, bad_p_mono, x), "invalid CSC structure", info = info)
    expect_error(fn(i, bad_p_end, x), "invalid CSC structure", info = info)
    expect_error(fn(bad_i_high, p, x), "invalid CSC structure", info = info)
    expect_error(fn(bad_i_neg, p, x), "invalid CSC structure", info = info)
    expect_error(fn(bad_i_na, p, x), "invalid CSC structure", info = info)
    expect_error(fn(i, p, x[-1L]), "invalid CSC structure", info = info)
  }
})

test_that("C28: well-formed CSC input still solves through the validated entries", {
  skip_if_not_installed("Matrix")
  set.seed(29)
  S <- Matrix::rsparsematrix(60L, 60L, density = 0.1)
  S <- as(Matrix::forceSymmetric(S + Matrix::t(S)), "generalMatrix")
  ref <- eigen(as.matrix(S), symmetric = TRUE, only.values = TRUE)$values
  fit <- eigs(S, 3L)
  expect_equal(sort(abs(fit$values), decreasing = TRUE),
               sort(abs(ref), decreasing = TRUE)[1:3], tolerance = 1e-6)
  A <- Matrix::rsparsematrix(50L, 30L, density = 0.2)
  sv <- svds(A, 3L)
  expect_equal(sv$d, svd(as.matrix(A))$d[1:3], tolerance = 1e-6)
})
