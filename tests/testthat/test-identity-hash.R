# C45: native structural identity hash (src/identity_hash.cpp) and the
# identity-format guard on persisted plans and restart states.

h <- function(x) eigencore:::stable_raw_hash(x)

test_that("C45: identity hash is a fixed, platform-independent function of values", {
  digest <- h(list(a = 1:3, b = "x", c = c(1.5, NA, NaN)))
  expect_match(digest, "^[0-9a-f]{32}$")
  # Pinned value: a change here changes every persisted identity and token,
  # so it must come with a new identity_hash_format().
  expect_identical(digest, "cfda38219a438b067edd10ed2d61ea20")
  expect_identical(eigencore:::identity_hash_format(), "eigencore-identity-hash-v2")
})

test_that("C45: identity hash follows identical() semantics on values and structure", {
  expect_identical(h(c(0, 1)), h(c(-0, 1)))
  expect_identical(h(NaN), h(-NaN))
  expect_false(identical(h(NaN), h(NA_real_)))
  expect_identical(h(complex(real = -0, imaginary = 1)), h(complex(real = 0, imaginary = 1)))
  expect_false(identical(h(1L), h(1)))
  expect_false(identical(h(TRUE), h(1L)))
  expect_false(identical(h(list(1)), h(1)))
  expect_false(identical(h(matrix(1:4, 2)), h(1:4)))
  expect_false(identical(h(matrix(1:4, 2)), h(matrix(1:4, 1))))
  expect_false(identical(h(c("ab", "c")), h(c("a", "bc"))))
  expect_false(identical(h(NA_character_), h("NA")))
  expect_false(identical(h(1:3), h(1:4)))
  # Attribute order does not matter; attribute values and names do.
  expect_identical(h(structure(1, a = 1, b = 2)), h(structure(1, b = 2, a = 1)))
  expect_false(identical(h(structure(1, a = 1)), h(structure(1, a = 2))))
  expect_false(identical(h(structure(1, a = 1)), h(structure(1, b = 1))))
  # Dimnames are part of the identity (as under the serialised digest).
  m <- matrix(1, 2, 2)
  named <- m
  dimnames(named) <- list(c("a", "b"), NULL)
  expect_false(identical(h(m), h(named)))
  # UTF-8 and native encodings of the same string hash equal.
  s <- "café"
  latin <- iconv(s, "UTF-8", "latin1")
  expect_identical(Encoding(latin), "latin1")
  expect_identical(h(s), h(latin))
  # Odd-length integer and byte tails are padded unambiguously.
  expect_false(identical(h(c(1L, 0L)), h(1L)))
  expect_false(identical(h(as.raw(c(1, 0))), h(as.raw(1))))
  # Objects without a value fast path still hash deterministically.
  f <- function(x) x + 1
  expect_identical(h(list(f = f)), h(list(f = f)))
  expect_match(h(new.env()), "^[0-9a-f]{32}$")
})

test_that("C45: one changed entry changes the operator identity (dense and sparse)", {
  set.seed(7)
  D <- crossprod(matrix(rnorm(40 * 40), 40))
  id <- operator_identity(as_operator(D))$revision
  expect_identical(operator_identity(as_operator(D + 0))$revision, id)
  D2 <- D
  D2[3, 5] <- D2[3, 5] + 1e-13
  expect_false(identical(operator_identity(as_operator(D2))$revision, id))

  S <- Matrix::rsparsematrix(200, 200, density = 0.02)
  S <- S + Matrix::t(S)
  ids <- operator_identity(as_operator(S))$revision
  S2 <- S
  S2@x[1] <- S2@x[1] * 2
  expect_false(identical(operator_identity(as_operator(S2))$revision, ids))
  # Same values, moved to another position (only i changes).
  S3 <- S
  S3@i[1] <- if (S3@i[1] == 0L) 1L else 0L
  expect_false(identical(h(S3), h(S)))
})

test_that("C45: identity is identical across separate R sessions", {
  skip_on_cran()
  rscript <- file.path(R.home("bin"), "Rscript")
  skip_if_not(file.exists(rscript))
  code <- paste(
    "suppressMessages(library(eigencore));",
    "set.seed(3); D <- crossprod(matrix(rnorm(400), 20));",
    "S <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, -0, 2.5), dims = c(3, 3));",
    "cat(operator_identity(as_operator(D))$revision,",
    "operator_identity(as_operator(S))$revision,",
    "plan_solver(eigen_problem(D), k = 2)$serialization$operator_identity_token)"
  )
  run <- function() {
    out <- suppressWarnings(system2(
      rscript, c("-e", shQuote(code)), stdout = TRUE, stderr = FALSE,
      env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))
    ))
    if (!is.null(attr(out, "status"))) return(NULL)
    out
  }
  first <- run()
  second <- run()
  skip_if(is.null(first) || is.null(second), "child R session unavailable")
  expect_identical(first, second)
  set.seed(3)
  D <- crossprod(matrix(rnorm(400), 20))
  S <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, -0, 2.5), dims = c(3, 3))
  here <- paste(
    operator_identity(as_operator(D))$revision,
    operator_identity(as_operator(S))$revision,
    plan_solver(eigen_problem(D), k = 2)$serialization$operator_identity_token
  )
  expect_identical(first, here)
})

test_that("C45: plans and restart states from an older identity format ask for a re-plan", {
  plan <- plan_solver(eigen_problem(diag(c(7, 5, 3, 1))), k = 2L)
  expect_identical(plan$serialization$hash_format, "eigencore-identity-hash-v2")
  old <- plan
  old$serialization$hash_format <- NULL
  err <- tryCatch(solve(old), error = identity)
  expect_s3_class(err, "eigencore_plan_error")
  expect_identical(err$code, "identity_format_changed")
  expect_match(conditionMessage(err), "identity format changed")
  expect_match(conditionMessage(err), "Re-plan")

  fit <- solve(plan_solver(eigen_problem(diag(seq(30, 1))), k = 3L))
  state <- restart_state(fit)
  expect_identical(state$serialization$hash_format, "eigencore-identity-hash-v2")
  stale <- state
  stale$serialization$hash_format <- NULL
  err <- tryCatch(restart_state(stale), error = identity)
  expect_s3_class(err, "eigencore_restart_state_error")
  expect_identical(err$code, "identity_format_changed")
  expect_match(conditionMessage(err), "identity format changed")
})
