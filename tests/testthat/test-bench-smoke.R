test_that("benchmark harness produces certificate-inclusive rows", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- path_laplacian(30)
  eigen_rows <- benchmark_eigen_case(
    A,
    k = 2,
    target = smallest(),
    methods = "eigencore",
    iterations = 1,
    seed = 300
  )

  X <- tall_skinny_sparse(40, 12, density = 0.08, seed = 301)
  svd_rows <- benchmark_svd_case(
    X,
    rank = 2,
    methods = "eigencore",
    iterations = 1,
    seed = 301
  )

  required <- c(
    "method", "median", "mem_alloc", "max_residual", "max_backward_error",
    "orthogonality_loss", "certificate_passed", "certificate_type",
    "norm_bound_type", "scale_is_estimate", "nconv", "preconditioner_kind",
    "preconditioner_native", "preconditioner_calls", "seed", "pkg_version"
  )
  expect_true(all(required %in% names(eigen_rows)))
  expect_true(all(required %in% names(svd_rows)))
  expect_true(all(is.finite(eigen_rows$median)))
  expect_true(all(is.finite(svd_rows$median)))
  expect_true(all(is.finite(eigen_rows$max_backward_error)))
  expect_true(all(is.finite(svd_rows$max_backward_error)))
  expect_true(all(eigen_rows$nconv >= 0))
  expect_true(all(svd_rows$nconv >= 0))
})

test_that("benchmark rows expose typed preconditioner diagnostics", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- path_laplacian(30)
  rows <- benchmark_eigen_case(
    A,
    k = 2,
    target = smallest(),
    methods = "eigencore_lobpcg_tridiagonal",
    iterations = 1,
    seed = 303
  )

  expect_equal(rows$preconditioner_kind, "shifted_tridiagonal")
  expect_true(rows$preconditioner_native)
  expect_gte(rows$preconditioner_calls, 1L)
})

test_that("native Hermitian benchmark gate reports current pass/fail state", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")
  skip_if_not_installed("RSpectra")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- path_laplacian(30)
  out <- benchmark_native_hermitian_gate(
    A,
    k = 2,
    target = smallest(),
    methods = c("eigencore", "RSpectra"),
    iterations = 1,
    seed = 702
  )

  expect_true(all(c("rows", "gate") %in% names(out)))
  expect_true(all(c("eigencore_certified", "speed_ratio_vs_best_reference", "passed") %in% names(out$gate)))
  expect_equal(out$gate$requested, 2)
  expect_true(is.logical(out$gate$passed))
  expect_true(all(c("eigencore", "RSpectra") %in% out$rows$method))
})

test_that("performance baseline helpers cover release regimes", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)
  if (!exists("dense_low_rank_spd", mode = "function")) {
    source(test_path("../../inst/benchmarks/_helpers.R"))
  }

  dense <- dense_low_rank_spd(20, rank = 4, seed = 11)
  expect_equal(dim(dense), c(20L, 20L))
  expect_true(isSymmetric(dense))

  slow <- slow_decay_svd_matrix(18, 9, seed = 12)
  expect_equal(dim(slow), c(18L, 9L))

  pair <- generalized_spd_pair(12, sparse = FALSE, seed = 13)
  expect_equal(dim(pair$A), c(12L, 12L))
  expect_equal(dim(pair$B), c(12L, 12L))

  rows <- data.frame(
    method = c("eigencore", "reference"),
    median = c(2, 1),
    mem_alloc = c(4, 2),
    certificate_passed = c(TRUE, TRUE),
    nconv = c(2L, 2L)
  )
  gate <- evaluate_reference_gate(rows, requested = 2L)
  expect_false(gate$passed)
  expect_equal(gate$speed_ratio_vs_best_reference, 0.5)
  expect_equal(gate$memory_ratio_vs_best_reference, 0.5)
  expect_equal(release_speed_gate("hermitian"), 1.25)
  expect_equal(release_speed_gate("svd"), 1.5)
  expect_equal(release_speed_gate("randomized_svd"), 2.0)
})
