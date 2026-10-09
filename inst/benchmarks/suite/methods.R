# Benchmark suite: method adapters.
#
# suite_adapter(method, case, prob, tol) returns NULL when the method does not
# apply to the task (e.g. irlba for an eigenproblem), otherwise a list with
#   skip     NULL, or a reason string (e.g. base R on a 20000 x 20000 matrix)
#   call     human-readable call recorded in results.csv
#   run      zero-argument function: the timed region
#   extract  function(raw) -> standardised result (see suite_std())
#
# Matvec conventions (column "matvecs"): the number of operator applications
# to single vectors, counting products with A and with A' separately, as
# reported by each package itself:
#   eigencore  work(fit)$operator_columns + adjoint_columns (solve phase;
#              certification columns are reported separately). NA when the
#              work record is incomplete (e.g. explicit-Gram SVD paths).
#   RSpectra   nops (eigs/eigs_sym: applications of op, which for sigma
#              modes are linear solves; svds: A and A' products, verified with
#              counting callbacks).
#   irlba      mprod (A and A' products).
#   PRIMME     stats$numMatvecs.
#   base       NA (dense LAPACK, no operator applications).

SUITE_METHODS <- c("eigencore", "RSpectra", "irlba", "base", "PRIMME")

suite_method_available <- function(method) {
  switch(method,
         eigencore = requireNamespace("eigencore", quietly = TRUE),
         base = TRUE,
         requireNamespace(method, quietly = TRUE))
}

suite_std <- function(values, vectors = NULL, left = NULL, u = NULL, v = NULL,
                      label = NA_character_, certified = NA, own_backward_error = NA_real_,
                      matvecs = NA_real_, matvec_source = NA_character_,
                      ec_operator_columns = NA_real_, ec_adjoint_columns = NA_real_,
                      ec_cert_columns = NA_real_, iterations = NA_real_,
                      warnings = character()) {
  list(values = values, vectors = vectors, left = left, u = u, v = v,
       label = label, certified = certified, own_backward_error = own_backward_error,
       matvecs = matvecs, matvec_source = matvec_source,
       ec_operator_columns = ec_operator_columns, ec_adjoint_columns = ec_adjoint_columns,
       ec_cert_columns = ec_cert_columns, iterations = iterations,
       warnings = warnings)
}

suite_limits <- list(base_sym_n = 4000L, base_nonsym_n = 2000L,
                     base_svd_cells = 5e6, base_gen_n = 3000L)

suite_select <- function(values, target, k, sigma = NA_real_) {
  ord <- switch(target,
                LA = order(Re(values), decreasing = TRUE),
                SA = order(Re(values)),
                LM = order(-Mod(values), Im(values)),
                near = order(abs(values - sigma)),
                top = order(values, decreasing = TRUE),
                all = order(Re(values)))
  ord[seq_len(min(k, length(ord)))]
}

suite_ec_target <- function(case) {
  switch(case$target,
         LA = eigencore::largest(), SA = eigencore::smallest(),
         LM = eigencore::largest_magnitude(), near = eigencore::nearest(case$sigma),
         top = eigencore::largest())
}

# eigencore's own verdict; NA when no certificate was computed (values only).
suite_ec_certified <- function(cert) {
  if (is.null(cert) || identical(cert$certificate_type, "uncomputed")) NA else isTRUE(cert$passed)
}

suite_ec_std <- function(fit, task) {
  w <- tryCatch(eigencore::work(fit), error = function(e) NULL)
  cert <- fit$certificate
  op <- if (!is.null(w)) w$operator_columns else NA
  adj <- if (!is.null(w)) w$adjoint_columns else NA
  certc <- if (!is.null(w)) sum(c(w$certification_operator_columns,
                                  w$certification_adjoint_columns)) else NA
  complete <- !is.null(w) && isTRUE(w$complete)
  mv <- if (complete) op + adj else NA_real_
  src <- if (complete) "eigencore work(): operator_columns + adjoint_columns" else
    sprintf("eigencore work() incomplete (legacy_matvecs=%s)",
            if (is.null(w)) "NA" else format(w$legacy_matvecs))
  base <- list(label = fit$method %||% NA_character_,
               certified = suite_ec_certified(cert),
               own_backward_error = if (is.null(cert)) NA_real_ else
                 suppressWarnings(max(cert$backward_error, na.rm = TRUE)),
               matvecs = mv, matvec_source = src,
               ec_operator_columns = op %||% NA, ec_adjoint_columns = adj %||% NA,
               ec_cert_columns = certc %||% NA,
               iterations = if (!is.null(w)) w$iterations else NA)
  if (task == "svd") {
    out <- suite_std(values = eigencore::values(fit),
                     u = eigencore::left_vectors(fit), v = eigencore::right_vectors(fit))
  } else {
    left <- if (task == "nonsym") tryCatch(eigencore::left_vectors(fit), error = function(e) NULL) else NULL
    vec <- tryCatch(eigencore::vectors(fit), error = function(e) NULL)
    out <- suite_std(values = eigencore::values(fit), vectors = vec, left = left)
  }
  utils::modifyList(out, base)
}

suite_adapter <- function(method, case, prob, tol) {
  fn <- switch(method,
    eigencore = suite_adapter_eigencore,
    RSpectra = suite_adapter_rspectra,
    irlba = suite_adapter_irlba,
    PRIMME = suite_adapter_primme,
    base = suite_adapter_base,
    stop("unknown method ", method))
  fn(case, prob, tol)
}

suite_adapter_eigencore <- function(case, prob, tol) {
  A <- prob$A
  k <- case$k
  seed <- case$seed
  switch(case$task,
    sym = list(
      call = sprintf("eig_partial(A, k=%d, target=%s, tol=%g, seed=%d)", k,
                     if (case$target == "near") sprintf("nearest(%g)", case$sigma) else
                       c(LA = "largest()", SA = "smallest()")[[case$target]], tol, seed),
      run = function() eigencore::eig_partial(A, k = k, target = suite_ec_target(case), tol = tol, seed = seed),
      extract = function(fit) suite_ec_std(fit, "sym")),
    nonsym = list(
      call = sprintf("eig_partial(A, k=%d, target=largest_magnitude(), tol=%g, seed=%d, left_vectors=\"none\")", k, tol, seed),
      run = function() eigencore::eig_partial(A, k = k, target = eigencore::largest_magnitude(),
                                              tol = tol, seed = seed, left_vectors = "none"),
      extract = function(fit) suite_ec_std(fit, "nonsym")),
    gen = list(
      call = sprintf("eig_partial(A, k=%d, target=smallest(), B=B, tol=%g, seed=%d)", k, tol, seed),
      run = function() eigencore::eig_partial(A, k = k, target = suite_ec_target(case), B = prob$B,
                                              tol = tol, seed = seed),
      extract = function(fit) suite_ec_std(fit, "gen")),
    svd = if (isTRUE(case$center)) list(
      call = sprintf("svd_partial(center(A), rank=%d, tol=%g, seed=%d)", k, tol, seed),
      run = function() eigencore::svd_partial(eigencore::center(A), rank = k, tol = tol, seed = seed),
      extract = function(fit) suite_ec_std(fit, "svd")) else list(
      call = sprintf("svd_partial(A, rank=%d, tol=%g, seed=%d)", k, tol, seed),
      run = function() eigencore::svd_partial(A, rank = k, tol = tol, seed = seed),
      extract = function(fit) suite_ec_std(fit, "svd")),
    full = list(
      call = sprintf("eig_full(A, vectors=%s)", case$vectors),
      run = function() eigencore::eig_full(A, vectors = case$vectors),
      extract = function(fit) {
        vals <- eigencore::values(fit)
        vec <- if (case$vectors) eigencore::vectors(fit) else NULL
        utils::modifyList(suite_std(values = vals, vectors = vec),
                          list(label = fit$method %||% "eig_full",
                               certified = suite_ec_certified(fit$certificate)))
      }),
    NULL)
}

suite_rs_counts <- function(res, label) {
  list(matvecs = res$nops %||% NA_real_, matvec_source = "RSpectra nops",
       iterations = res$niter %||% NA_real_, label = label)
}

suite_adapter_rspectra <- function(case, prob, tol) {
  A <- prob$A
  k <- case$k
  seed <- case$seed
  opts <- list(tol = tol, maxitr = 1000L)
  switch(case$task,
    sym = {
      which <- if (case$target == "near") "LM" else case$target
      sigma <- if (case$target == "near") case$sigma else NULL
      list(
        call = sprintf("eigs_sym(A, k=%d, which=\"%s\"%s, opts=list(tol=%g, maxitr=1000))", k, which,
                       if (is.null(sigma)) "" else sprintf(", sigma=%g", sigma), tol),
        run = function() { suite_set_seed(seed); RSpectra::eigs_sym(A, k, which = which, sigma = sigma, opts = opts) },
        extract = function(res) utils::modifyList(
          suite_std(values = res$values, vectors = res$vectors),
          suite_rs_counts(res, if (is.null(sigma)) "eigs_sym" else "eigs_sym shift-invert (sparse LU)")))
    },
    nonsym = list(
      call = sprintf("eigs(A, k=%d, which=\"LM\", opts=list(tol=%g, maxitr=1000))", k, tol),
      run = function() { suite_set_seed(seed); RSpectra::eigs(A, k, which = "LM", opts = opts) },
      extract = function(res) utils::modifyList(
        suite_std(values = res$values, vectors = res$vectors), suite_rs_counts(res, "eigs"))),
    gen = list(
      # RSpectra has no generalized mode: the standard user route is a sparse
      # Cholesky B = R'R and the symmetric operator R^-T A R^-1 (factorisation
      # included in the timing), then x = R^-1 y.
      call = sprintf("R <- chol(B); eigs_sym(function(y) R^-T A R^-1 y, k=%d, which=\"SA\", opts=list(tol=%g, maxitr=1000))", k, tol),
      run = function() {
        suite_set_seed(seed)
        R <- Matrix::chol(Matrix::forceSymmetric(prob$B))
        Rt <- Matrix::t(R)
        op <- function(y, args) as.numeric(Matrix::solve(Rt, A %*% Matrix::solve(R, y)))
        res <- RSpectra::eigs_sym(op, k, which = "SA", n = nrow(A), opts = opts)
        res$vectors <- as.matrix(Matrix::solve(R, res$vectors))
        res
      },
      extract = function(res) utils::modifyList(
        suite_std(values = res$values, vectors = res$vectors),
        suite_rs_counts(res, "eigs_sym on chol-transformed operator"))),
    svd = if (isTRUE(case$center)) {
      mu <- prob$mu
      list(
        call = sprintf("svds(function(x) A x - 1 (mu'x), Atrans=function(y) A'y - mu sum(y), k=%d, opts=list(tol=%g))", k, tol),
        run = function() {
          suite_set_seed(seed)
          f <- function(x, args) as.numeric(A %*% x) - sum(mu * x)
          ft <- function(y, args) as.numeric(Matrix::crossprod(A, y)) - mu * sum(y)
          RSpectra::svds(f, k, nu = k, nv = k, Atrans = ft, dim = dim(A), opts = opts)
        },
        extract = function(res) utils::modifyList(
          suite_std(values = res$d, u = res$u, v = res$v), suite_rs_counts(res, "svds (function operator)")))
    } else list(
      call = sprintf("svds(A, k=%d, opts=list(tol=%g, maxitr=1000))", k, tol),
      run = function() { suite_set_seed(seed); RSpectra::svds(A, k, nu = k, nv = k, opts = opts) },
      extract = function(res) utils::modifyList(
        suite_std(values = res$d, u = res$u, v = res$v), suite_rs_counts(res, "svds"))),
    NULL)
}

suite_adapter_irlba <- function(case, prob, tol) {
  if (case$task != "svd") return(NULL)
  A <- prob$A
  k <- case$k
  seed <- case$seed
  ctr <- if (isTRUE(case$center)) prob$mu else NULL
  list(
    call = sprintf("irlba(A, nv=%d, nu=%d, tol=%g%s)", k, k, tol, if (is.null(ctr)) "" else ", center=colMeans(A)"),
    run = function() {
      suite_set_seed(seed)
      if (is.null(ctr)) irlba::irlba(A, nv = k, nu = k, tol = tol) else
        irlba::irlba(A, nv = k, nu = k, tol = tol, center = ctr)
    },
    extract = function(res) utils::modifyList(
      suite_std(values = res$d, u = res$u, v = res$v),
      list(matvecs = res$mprod %||% NA_real_, matvec_source = "irlba mprod",
           iterations = res$iter %||% NA_real_, label = "irlba")))
}

suite_adapter_primme <- function(case, prob, tol) {
  A <- prob$A
  k <- case$k
  seed <- case$seed
  pstats <- function(res, label) list(matvecs = res$stats$numMatvecs %||% NA_real_,
                                      matvec_source = "PRIMME stats$numMatvecs",
                                      iterations = res$stats$numOuterIterations %||% NA_real_,
                                      label = label)
  switch(case$task,
    sym = {
      which <- if (case$target == "near") case$sigma else case$target
      list(call = sprintf("PRIMME::eigs_sym(A, NEig=%d, which=%s, tol=%g)", k, format(which), tol),
           run = function() { suite_set_seed(seed); PRIMME::eigs_sym(A, NEig = k, which = which, tol = tol) },
           extract = function(res) utils::modifyList(suite_std(values = res$values, vectors = res$vectors),
                                                     pstats(res, "PRIMME eigs_sym")))
    },
    gen = list(call = sprintf("PRIMME::eigs_sym(A, NEig=%d, which=\"SA\", B=B, tol=%g)", k, tol),
               run = function() { suite_set_seed(seed); PRIMME::eigs_sym(A, NEig = k, which = "SA", B = prob$B, tol = tol) },
               extract = function(res) utils::modifyList(suite_std(values = res$values, vectors = res$vectors),
                                                         pstats(res, "PRIMME eigs_sym (generalized)"))),
    svd = if (isTRUE(case$center)) NULL else
      list(call = sprintf("PRIMME::svds(A, NSvals=%d, which=\"L\", tol=%g)", k, tol),
           run = function() { suite_set_seed(seed); PRIMME::svds(A, NSvals = k, which = "L", tol = tol) },
           extract = function(res) utils::modifyList(suite_std(values = res$d, u = res$u, v = res$v),
                                                     pstats(res, "PRIMME svds"))),
    NULL)
}

suite_adapter_base <- function(case, prob, tol) {
  A <- prob$A
  k <- case$k
  n <- ncol(A)
  lim <- suite_limits
  too_big <- function(reason) list(skip = reason, call = NA_character_, run = NULL, extract = NULL)
  dense_sel <- function(vals, vecs) {
    idx <- suite_select(vals, case$target, k, case$sigma)
    suite_std(values = vals[idx], vectors = vecs[, idx, drop = FALSE], label = "LAPACK")
  }
  switch(case$task,
    sym = if (n > lim$base_sym_n) too_big(sprintf("dense eigen() skipped: n=%d > %d", n, lim$base_sym_n)) else
      list(call = "eigen(as.matrix(A), symmetric=TRUE) then select",
           run = function() eigen(as.matrix(A), symmetric = TRUE),
           extract = function(res) dense_sel(res$values, res$vectors)),
    nonsym = if (n > lim$base_nonsym_n) too_big(sprintf("dense eigen() skipped: n=%d > %d", n, lim$base_nonsym_n)) else
      list(call = "eigen(as.matrix(A)) then select",
           run = function() eigen(as.matrix(A)),
           extract = function(res) dense_sel(res$values, res$vectors)),
    gen = if (n > lim$base_gen_n) too_big(sprintf("dense generalized eigen skipped: n=%d > %d", n, lim$base_gen_n)) else
      list(call = "R <- chol(B); eigen(R^-T A R^-1, symmetric=TRUE); x = R^-1 y",
           run = function() {
             R <- chol(as.matrix(prob$B))
             C <- backsolve(R, t(backsolve(R, as.matrix(A), transpose = TRUE)), transpose = TRUE)
             e <- eigen((C + t(C)) / 2, symmetric = TRUE)
             e$vectors <- backsolve(R, e$vectors)
             e
           },
           extract = function(res) dense_sel(res$values, res$vectors)),
    svd = {
      cells <- as.numeric(nrow(A)) * ncol(A)
      if (cells > lim$base_svd_cells) too_big(sprintf("dense svd() skipped: %d x %d", nrow(A), ncol(A))) else
        list(call = sprintf("svd(%s, nu=%d, nv=%d)", if (isTRUE(case$center)) "scale(as.matrix(A), scale=FALSE)" else "as.matrix(A)", k, k),
             run = function() {
               X <- as.matrix(A)
               if (isTRUE(case$center)) X <- sweep(X, 2L, prob$mu)
               svd(X, nu = k, nv = k)
             },
             extract = function(res) suite_std(values = res$d[seq_len(k)], u = res$u, v = res$v, label = "LAPACK"))
    },
    full = list(call = sprintf("eigen(A, symmetric=TRUE, only.values=%s)", !case$vectors),
                run = function() eigen(A, symmetric = TRUE, only.values = !case$vectors),
                extract = function(res) suite_std(values = res$values, vectors = res$vectors, label = "LAPACK")),
    NULL)
}
