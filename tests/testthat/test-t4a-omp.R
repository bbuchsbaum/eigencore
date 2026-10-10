# Tranche 4a: P8 (OpenMP CSC kernels, eigencore.threads), P11 (block-apply
# output allocation) and C36 (OpenMP reorthogonalisation during sparse
# solves).

with_options <- function(opts, code) {
  old <- options(opts)
  on.exit(options(old), add = TRUE)
  force(code)
}

thread_info <- function() eigencore:::native_thread_info()
has_openmp <- function() isTRUE(thread_info()[["openmp"]] == 1L)

ref_apply <- function(A, X, transpose, alpha, beta, Y, col_means = NULL,
                      weights = NULL) {
  M <- as.matrix(A)
  if (!is.null(col_means)) {
    M <- sweep(sweep(M, 2L, col_means, "-"), 2L, weights, "*")
  }
  alpha * (if (transpose) crossprod(M, X) else M %*% X) + beta * Y
}

# Unsorted row indices inside each column (a valid but unusual dgCMatrix).
reverse_rows_within_columns <- function(A) {
  ord <- unlist(lapply(seq_len(ncol(A)), function(j) {
    if (A@p[j + 1L] > A@p[j]) rev((A@p[j] + 1L):A@p[j + 1L]) else integer()
  }))
  A@i <- A@i[ord]
  A@x <- A@x[ord]
  A
}

test_that("default thread count follows the CRAN / OMP_NUM_THREADS policy", {
  # cores / quota = NA: no physical-core or cgroup information.
  f <- function(env, processors, cores = NA, quota = NA) {
    eigencore:::eigencore_default_threads(env, processors = processors,
                                          cores = cores, quota = quota)
  }
  expect_identical(f(c(`_R_CHECK_LIMIT_CORES_` = "TRUE"), processors = 16L), 1L)
  expect_identical(f(c(`_R_CHECK_LIMIT_CORES_` = "warn"), processors = 16L), 1L)
  expect_identical(f(c(`_R_CHECK_LIMIT_CORES_` = "false"), processors = 16L), 8L)
  expect_identical(f(c(`_R_CHECK_PACKAGE_NAME_` = "eigencore"), processors = 16L), 1L)
  expect_identical(f(c(OMP_NUM_THREADS = "3,2"), processors = 16L), 3L)
  expect_identical(f(c(OMP_NUM_THREADS = "junk"), processors = 2L), 2L)
  expect_identical(f(character(), processors = 32L), 8L)
  expect_identical(f(character(), processors = 3L), 3L)
  expect_identical(f(character(), processors = NA), 1L)
  # P18: physical cores and a container CPU quota cap the default.
  expect_identical(f(character(), processors = 16L, cores = 6L), 6L)
  expect_identical(f(character(), processors = 16L, cores = 12L), 8L)
  expect_identical(f(character(), processors = 16L, quota = 2.5), 2L)
  expect_identical(f(character(), processors = 4L, quota = 0.5), 1L)
  expect_identical(f(character(), processors = 4L, cores = 8L, quota = 16), 4L)
  expect_identical(f(c(OMP_NUM_THREADS = "6"), processors = 4L, quota = 1), 6L)
})

test_that("P18: cgroup CPU quotas are read from cgroup v2 and v1 files", {
  quota <- eigencore:::eigencore_cgroup_cpu_quota
  root <- tempfile("cgroup")
  self <- tempfile("self")
  on.exit(unlink(c(root, self), recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "pod"), recursive = TRUE)
  writeLines("0::/pod", self)
  writeLines("150000 100000", file.path(root, "pod", "cpu.max"))
  expect_equal(quota(root, self), 1.5)
  writeLines("max 100000", file.path(root, "pod", "cpu.max"))
  expect_identical(quota(root, self), NA_real_)
  v1 <- file.path(root, "cpu,cpuacct", "docker")
  dir.create(v1, recursive = TRUE)
  unlink(file.path(root, "pod", "cpu.max"))
  writeLines(c("4:cpu,cpuacct:/docker", "0::/"), self)
  writeLines("200000", file.path(v1, "cpu.cfs_quota_us"))
  writeLines("100000", file.path(v1, "cpu.cfs_period_us"))
  expect_equal(quota(root, self), 2)
  writeLines("-1", file.path(v1, "cpu.cfs_quota_us"))
  expect_identical(quota(root, self), NA_real_)
  expect_identical(quota(file.path(root, "missing"), tempfile()), NA_real_)
})

test_that("eigencore.threads is read at every native call", {
  info <- thread_info()
  expect_named(info, c("openmp", "processors", "default", "current",
                       "blas_kind", "blas_threads", "effective", "reductions"))
  default <- info[["default"]]
  expect_gte(default, 1L)
  if (!has_openmp()) {
    expect_identical(default, 1L)
  }
  expect_identical(
    with_options(list(eigencore.threads = 3), thread_info()[["current"]]),
    if (has_openmp()) 3L else 1L
  )
  expect_identical(
    with_options(list(eigencore.threads = 2L), as.integer(eigencore:::eigencore_threads())),
    if (has_openmp()) 2L else 1L
  )
  # Invalid values fall back to the load-time default; huge ones are capped.
  for (bad in list(0, -1, NA_real_, NA_integer_, "4", Inf)) {
    expect_identical(
      with_options(list(eigencore.threads = bad), thread_info()[["current"]]),
      default
    )
  }
  expect_identical(
    with_options(list(eigencore.threads = 1e6), thread_info()[["current"]]),
    if (has_openmp()) 256L else 1L
  )
  expect_identical(with_options(list(eigencore.threads = NULL),
                                thread_info()[["current"]]), default)
  expect_error(eigencore:::eigencore_threads(0), "positive integer")
  expect_error(eigencore:::eigencore_threads(1.5), "positive integer")
  expect_error(eigencore:::eigencore_threads(c(1, 2)), "positive integer")
})

test_that("P8: cached parallel CSC kernels are bitwise thread-count invariant", {
  set.seed(401)
  shapes <- list(
    square = Matrix::rsparsematrix(3000, 3000, density = 0.006),  # CSR gather
    tall = Matrix::rsparsematrix(12000, 600, density = 0.008),    # row slabs
    wide = Matrix::rsparsematrix(600, 9000, density = 0.008)
  )
  for (name in names(shapes)) {
    A <- shapes[[name]]
    mu <- rnorm(ncol(A))
    w <- runif(ncol(A), 0.5, 2)
    for (b in c(1L, 2L, 3L, 10L, 13L)) {
      for (transpose in c(FALSE, TRUE)) {
        X <- matrix(rnorm((if (transpose) nrow(A) else ncol(A)) * b), ncol = b)
        X[sample(length(X), length(X) %/% 4L)] <- 0
        if (b >= 3L) X[, 2L] <- 0
        Y0 <- matrix(rnorm((if (transpose) ncol(A) else nrow(A)) * b), ncol = b)
        for (scaled in c(FALSE, TRUE)) {
          run <- function(threads, cache_mb = 4096) {
            with_options(
              list(eigencore.threads = threads, eigencore.csr_cache_mb = cache_mb),
              eigencore:::csc_apply_repeat(
                A, X, alpha = 0.7, beta = -0.3, Y = Y0, transpose = transpose,
                reps = 3L, col_means = if (scaled) mu, weights = if (scaled) w
              )$Y
            )
          }
          y1 <- run(1L)
          label <- sprintf("%s b=%d transpose=%s scaled=%s", name, b, transpose, scaled)
          expect_equal(
            y1,
            ref_apply(A, X, transpose, 0.7, -0.3, Y0,
                      if (scaled) mu, if (scaled) w),
            tolerance = 1e-12, ignore_attr = TRUE, label = label
          )
          expect_identical(run(2L), y1, label = label)
          expect_identical(run(4L), y1, label = label)
          # No room for the CSR copy / slabs: serial or column-chunk fallback.
          expect_identical(run(4L, cache_mb = 0), y1, label = label)
        }
      }
    }
  }
})

test_that("P8: unsorted row indices and per-call applies stay exact", {
  set.seed(402)
  A <- Matrix::rsparsematrix(8000, 1500, density = 0.01)
  U <- reverse_rows_within_columns(A)
  expect_false(identical(U@i, A@i))
  for (b in c(1L, 4L, 12L)) {
    X <- matrix(rnorm(ncol(A) * b), ncol = b)
    ref <- as.matrix(A %*% X)
    for (threads in c(1L, 4L)) {
      out <- with_options(list(eigencore.threads = threads),
                          eigencore:::csc_apply_repeat(U, X, reps = 3L)$Y)
      expect_equal(out, ref, tolerance = 1e-12, ignore_attr = TRUE)
      per_call <- with_options(list(eigencore.threads = threads),
                               eigencore:::csc_block_apply(A, X))
      expect_equal(per_call, ref, tolerance = 1e-12, ignore_attr = TRUE)
    }
    expect_identical(
      with_options(list(eigencore.threads = 4L), eigencore:::csc_block_apply(A, X)),
      with_options(list(eigencore.threads = 1L), eigencore:::csc_block_apply(A, X))
    )
    Xt <- matrix(rnorm(nrow(A) * b), ncol = b)
    expect_identical(
      with_options(list(eigencore.threads = 4L),
                   eigencore:::csc_block_apply(A, Xt, transpose = TRUE)),
      with_options(list(eigencore.threads = 1L),
                   eigencore:::csc_block_apply(A, Xt, transpose = TRUE))
    )
  }
})

test_that("P11: block applies allocate their output once and accept Y = NULL", {
  set.seed(403)
  A <- Matrix::rsparsematrix(40, 30, density = 0.2)
  X <- matrix(rnorm(30 * 3), 30)
  Y <- matrix(NaN, 40, 3, dimnames = list(NULL, c("a", "b", "c")))
  ref <- as.matrix(A %*% X)
  # beta = 0 never reads Y (NaN does not propagate) and keeps its dimnames.
  out <- eigencore:::csc_block_apply(A, X, alpha = 1, beta = 0, Y = Y)
  expect_equal(unname(out), unname(ref), tolerance = 1e-14)
  expect_identical(colnames(out), c("a", "b", "c"))
  expect_true(all(is.nan(Y)))
  Y2 <- matrix(1, 40, 3)
  out2 <- eigencore:::csc_block_apply(A, X, alpha = 2, beta = 0.5, Y = Y2)
  expect_equal(out2, 2 * ref + 0.5, tolerance = 1e-14, ignore_attr = TRUE)
  expect_identical(Y2, matrix(1, 40, 3))
  out_null <- .Call("eigencore_csc_block_apply", A@i, A@p, A@x, A@Dim, X, 1, 3,
                    NULL, FALSE, PACKAGE = "eigencore")
  expect_equal(out_null, ref, tolerance = 1e-14, ignore_attr = TRUE)
  D <- matrix(rnorm(40 * 30), 40)
  out_dense <- .Call("eigencore_dense_block_apply", D, X, 1, 5, NULL, FALSE,
                     PACKAGE = "eigencore")
  expect_equal(out_dense, D %*% X, tolerance = 1e-13)
  mu <- rnorm(30)
  w <- runif(30)
  out_cs <- .Call("eigencore_csc_centered_scaled_block_apply", A@i, A@p, A@x,
                  A@Dim, mu, w, X, 1, 7, NULL, FALSE, PACKAGE = "eigencore")
  expect_equal(out_cs, ref_apply(A, X, FALSE, 1, 0, 0, mu, w),
               tolerance = 1e-13, ignore_attr = TRUE)
})

test_that("C36: sparse solves agree across thread counts and restore BLAS threads", {
  skip_if_not(has_openmp(), "built without OpenMP")
  set.seed(404)
  B <- Matrix::rsparsematrix(4000, 4000, density = 1.5e-3)
  S <- methods::as(methods::as(Matrix::forceSymmetric(B + Matrix::t(B)),
                               "generalMatrix"), "CsparseMatrix")
  X <- Matrix::rsparsematrix(9000, 1200, density = 0.006)
  blas_before <- thread_info()[["blas_threads"]]
  solve_all <- function(threads) {
    with_options(list(eigencore.threads = threads), {
      set.seed(1)
      e <- eigs_sym(S, 6, which = "LA")
      set.seed(2)
      s <- svds(X, 6)
      list(e = e, s = s)
    })
  }
  r1 <- solve_all(1L)
  r2 <- solve_all(2L)
  r4 <- solve_all(4L)
  expect_identical(thread_info()[["blas_threads"]], blas_before)
  # OpenMP kernels are thread-count invariant: 2 and 4 threads agree exactly.
  expect_identical(r2$e$values, r4$e$values)
  expect_identical(r2$e$vectors, r4$e$vectors)
  expect_identical(r2$s$d, r4$s$d)
  expect_identical(r2$s$u, r4$s$u)
  # One thread uses BLAS for the reorthogonalisation: equal up to rounding.
  expect_equal(r1$e$values, r4$e$values, tolerance = 1e-12)
  expect_equal(r1$s$d, r4$s$d, tolerance = 1e-12)
  align <- function(U, V) sweep(U, 2L, sign(colSums(U * V)), "*")
  expect_equal(align(r1$e$vectors, r4$e$vectors), r4$e$vectors, tolerance = 1e-8)
  expect_equal(align(r1$s$v, r4$s$v), r4$s$v, tolerance = 1e-8)
})

test_that("P18: the parallel-efficiency governor changes team sizes, not results", {
  skip_if_not(has_openmp(), "built without OpenMP")
  governor <- function(cap) eigencore:::native_thread_info_governor(cap)
  on.exit(governor(0), add = TRUE)
  set.seed(405)
  A <- Matrix::rsparsematrix(3000, 3000, density = 0.006)
  tall <- Matrix::rsparsematrix(12000, 600, density = 0.008)
  B <- Matrix::rsparsematrix(4000, 4000, density = 1.5e-3)
  S <- methods::as(methods::as(Matrix::forceSymmetric(B + Matrix::t(B)),
                               "generalMatrix"), "CsparseMatrix")
  run <- function(cap) {
    with_options(list(eigencore.threads = 4L), {
      info <- governor(cap)
      expect_identical(info[["effective"]], if (cap == 0) 4L else as.integer(cap))
      out <- list()
      for (M in list(A, tall)) {
        for (transpose in c(FALSE, TRUE)) {
          X <- matrix(rnorm((if (transpose) nrow(M) else ncol(M)) * 12), ncol = 12)
          out[[length(out) + 1L]] <- eigencore:::csc_apply_repeat(
            M, X, transpose = transpose, reps = 3L)$Y
        }
      }
      set.seed(1)
      e <- eigs_sym(S, 6, which = "LA")
      set.seed(2)
      s <- svds(tall, 5)
      c(out, list(e$values, e$vectors, s$d, s$u, s$v))
    })
  }
  set.seed(7)
  full <- run(0)
  set.seed(7)
  one <- run(1)
  set.seed(7)
  two <- run(2)
  expect_identical(one, full)
  expect_identical(two, full)
  # Disabled governor: the configured count is used.
  with_options(list(eigencore.threads = 4L, eigencore.adaptive_threads = FALSE), {
    governor(1)
    expect_identical(thread_info()[["effective"]], 4L)
  })
})

test_that("C62: completeness gates use the smaller of elapsed and CPU solve time", {
  clock <- eigencore:::completeness_clock()
  expect_named(clock, c("elapsed", "cpu"))
  since <- eigencore:::completeness_seconds_since
  # Starved solve: 10 s elapsed, 1 s of CPU -> about 1 s.
  expect_equal(since(clock - c(elapsed = 10, cpu = 1)), 1, tolerance = 0.5)
  # Parallel solve: 1 s elapsed, 4 s of CPU -> about 1 s.
  expect_equal(since(clock - c(elapsed = 1, cpu = 4)), 1, tolerance = 0.5)
})
