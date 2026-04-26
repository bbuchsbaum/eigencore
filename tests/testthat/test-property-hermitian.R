test_that("Hermitian eig_partial certificate contract holds across spectra", {
  skip_on_cran()
  patterns <- c("uniform", "clustered", "exponential", "two_cluster")
  sizes    <- c(80L, 200L)
  seeds    <- 1:4
  k        <- 5L
  tol      <- 1e-8

  for (n in sizes) {
    for (pat in patterns) {
      for (seed in seeds) {
        A    <- random_symmetric_with_spectrum(n, pattern = pat, seed = seed)
        true <- sort(spectrum_pattern(pat, n), decreasing = TRUE)[seq_len(k)]
        info <- sprintf("n=%d, pattern=%s, seed=%d", n, pat, seed)

        fit  <- eig_partial(A, k = k, target = largest())
        cert <- fit$certificate
        expect_true(all(is.finite(cert$backward_error)), info = info)
        expect_lte(cert$max_backward_error, tol, label = paste0(info, "/backward"))
        expect_true(cert$passed,                       label = paste0(info, "/passed"))

        # eigenvalues must agree to a backward-error-scaled tolerance
        gap <- abs(fit$values - true) / pmax(abs(true), 1)
        expect_lt(max(gap), 1e-6, label = paste0(info, "/values"))
      }
    }
  }
})

test_that("smallest-target Hermitian eig_partial certificate contract holds", {
  skip_on_cran()
  patterns <- c("uniform", "exponential")
  sizes    <- c(80L, 200L)
  seeds    <- 1:3
  k        <- 4L
  tol      <- 1e-8

  for (n in sizes) {
    for (pat in patterns) {
      for (seed in seeds) {
        A    <- random_symmetric_with_spectrum(n, pattern = pat, seed = seed)
        true <- sort(spectrum_pattern(pat, n))[seq_len(k)]
        info <- sprintf("n=%d, pattern=%s, seed=%d", n, pat, seed)

        fit  <- eig_partial(A, k = k, target = smallest())
        cert <- fit$certificate
        expect_lte(cert$max_backward_error, tol, label = paste0(info, "/backward"))

        gap <- abs(fit$values - true) / pmax(abs(true), 1)
        expect_lt(max(gap), 1e-6, label = paste0(info, "/values"))
      }
    }
  }
})
