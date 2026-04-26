test_that("svd_partial certificate contract holds across SV patterns", {
  skip_on_cran()
  patterns <- c("uniform", "clustered", "exponential", "geometric")
  shapes   <- list(c(80L, 30L), c(120L, 60L))
  seeds    <- 1:3
  rank     <- 5L
  tol      <- 1e-8

  for (shape in shapes) {
    m <- shape[[1]]; n <- shape[[2]]
    for (pat in patterns) {
      for (seed in seeds) {
        full <- min(m, n)
        M    <- random_rectangular_with_pattern(m, n, rank = full,
                                                pattern = pat, seed = seed)
        true <- sort(spectrum_pattern(pat, full), decreasing = TRUE)[seq_len(rank)]
        info <- sprintf("m=%d, n=%d, pattern=%s, seed=%d", m, n, pat, seed)

        fit  <- svd_partial(M, rank = rank, target = largest())
        cert <- fit$certificate
        expect_true(all(is.finite(cert$backward_error)), info = info)
        expect_lte(cert$max_backward_error, tol,
                   label = paste0(info, "/backward"))

        gap <- abs(fit$d - true) / pmax(abs(true), 1)
        expect_lt(max(gap), 1e-6, label = paste0(info, "/values"))
      }
    }
  }
})

test_that("svd_partial handles a rank-deficient rectangular matrix", {
  skip_on_cran()
  m <- 100L; n <- 60L
  for (seed in 1:3) {
    M <- random_rectangular_with_pattern(m, n, rank = 30L,
                                         pattern = "exponential", seed = seed)
    info <- sprintf("seed=%d", seed)
    fit  <- svd_partial(M, rank = 5L, target = largest())
    expect_true(all(is.finite(fit$d)), info = info)
    expect_true(all(fit$d > 0), info = info)
    expect_true(!any(is.na(fit$certificate$backward_error)), info = info)
  }
})
