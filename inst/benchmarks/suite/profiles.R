# Benchmark suite: case catalogue and profiles.
#
# A case is a list with
#   id        unique, stable identifier (part of the reference cache key)
#   family    problem family (see generators.R), used by --families
#   group     "core" for fixed cases, or the sweep name for scaling curves
#   task      "sym" | "nonsym" | "svd" | "gen" | "full"
#   target    "LA" | "SA" | "LM" | "near" | "top" | "all"
#   k         number of wanted pairs / triplets (n for task "full")
#   sigma     shift for target "near"
#   center    TRUE for column-centred SVD (PCA)
#   vectors   FALSE for values-only "full" cases
#   seed      generator seed (also the solver seed)
#   sweep, sweep_value   scaling variable (NA for core cases)
#   gen       zero-argument function returning the generator list
#
# Profiles:
#   quick     CI smoke (a few minutes): scaled-down versions of every family
#   standard  the problem sizes quoted in README / docs/review-2026-10.md
#   scaling   sweeps over n, nnz and k for a few families

suite_case <- function(id, family, task, target, k, gen, seed, sigma = NA_real_,
                       center = FALSE, vectors = TRUE, group = "core",
                       sweep = NA_character_, sweep_value = NA_real_) {
  list(id = id, family = family, group = group, task = task, target = target,
       k = as.integer(k), sigma = sigma, center = center, vectors = vectors,
       seed = as.integer(seed), sweep = sweep, sweep_value = sweep_value,
       gen = gen)
}

suite_core_cases <- function(scale = c("standard", "quick")) {
  scale <- match.arg(scale)
  q <- scale == "quick"
  n_sp <- if (q) 1500L else 20000L
  m_svd <- if (q) 5000L else 50000L
  n_svd1 <- if (q) 300L else 2000L
  n_svd2 <- if (q) 1000L else 20000L
  n_dense <- if (q) 300L else 1500L
  g_si <- if (q) 40L else 100L
  g_gen <- if (q) 20L else 50L
  g_2d <- if (q) 40L else 150L
  pca_m <- if (q) 5000L else 50000L
  pca_n <- if (q) 200L else 1000L
  n_band <- if (q) 3000L else 20000L
  n_pl <- if (q) 400L else 2000L
  lr <- if (q) c(1000L, 200L) else c(5000L, 500L)

  list(
    suite_case(sprintf("sym_sparse_n%d_LA_k10", n_sp), "sym_sparse", "sym", "LA", 10,
               function() gen_sparse_sym(n_sp, 5, seed = 101L), seed = 101L),
    suite_case(sprintf("sym_sparse_n%d_SA_k10", n_sp), "sym_sparse", "sym", "SA", 10,
               function() gen_sparse_sym(n_sp, 5, seed = 101L), seed = 101L),
    suite_case(sprintf("nonsym_sparse_n%d_LM_k6", n_sp), "nonsym_sparse", "nonsym", "LM", 6,
               function() gen_sparse_nonsym(n_sp, 5, seed = 102L), seed = 102L),
    suite_case(sprintf("svd_sparse_%dx%d_k20", m_svd, n_svd1), "svd_sparse", "svd", "top", 20,
               function() gen_sparse_rect(m_svd, n_svd1, density = 0.002 * (if (q) 5 else 1), seed = 103L),
               seed = 103L),
    suite_case(sprintf("svd_sparse_%dx%d_k20", m_svd, n_svd2), "svd_sparse", "svd", "top", 20,
               function() gen_sparse_rect(m_svd, n_svd2, density = 0.0005 * (if (q) 5 else 1), seed = 104L),
               seed = 104L),
    suite_case(sprintf("dense_sym_n%d_LA_k10", n_dense), "dense_sym", "sym", "LA", 10,
               function() gen_dense_wishart(n_dense, seed = 105L), seed = 105L),
    suite_case(sprintf("eig_full_n%d_vectors", n_dense), "eig_full", "full", "all", n_dense,
               function() gen_dense_wishart(n_dense, seed = 105L), seed = 105L),
    suite_case(sprintf("eig_full_n%d_values", n_dense), "eig_full", "full", "all", n_dense,
               function() gen_dense_wishart(n_dense, seed = 105L), seed = 105L, vectors = FALSE),
    suite_case(sprintf("shift_invert_lap2d_g%d_sigma4.01_k6", g_si), "shift_invert", "sym", "near", 6,
               function() gen_laplacian_2d(g_si), seed = 106L, sigma = 4.01),
    suite_case(sprintf("generalized_fem2d_g%d_SA_k6", g_gen), "generalized", "gen", "SA", 6,
               function() gen_fem_pencil_2d(g_gen), seed = 107L),
    suite_case(sprintf("pca_centered_%dx%d_k10", pca_m, pca_n), "pca", "svd", "top", 10,
               function() gen_sparse_rect(pca_m, pca_n, density = 0.01, seed = 108L, positive = TRUE),
               seed = 108L, center = TRUE),
    suite_case(sprintf("banded_lap1d_n%d_SA_k8", n_band), "banded", "sym", "SA", 8,
               function() gen_laplacian_1d(n_band), seed = 109L),
    suite_case(sprintf("laplacian2d_g%d_SA_k10", g_2d), "laplacian2d", "sym", "SA", 10,
               function() gen_laplacian_2d(g_2d), seed = 110L),
    suite_case(sprintf("powerlaw_dense_n%d_LA_k10", n_pl), "powerlaw", "sym", "LA", 10,
               function() gen_powerlaw_dense(n_pl, alpha = 1, seed = 111L), seed = 111L),
    suite_case(sprintf("clustered_dense_n%d_LA_k10", n_pl), "clustered", "sym", "LA", 10,
               function() gen_clustered_dense(n_pl, seed = 112L), seed = 112L),
    suite_case(sprintf("lowrank_noise_%dx%d_k10", lr[1], lr[2]), "lowrank", "svd", "top", 10,
               function() gen_lowrank_noise(lr[1], lr[2], rank = 12L, noise = 1e-2, seed = 113L),
               seed = 113L)
  )
}

suite_scaling_cases <- function() {
  out <- list()
  add <- function(x) out[[length(out) + 1L]] <<- x
  # n sweep, random sparse symmetric, 5 nnz/row, LA k=10
  # (Matrix::rsparsematrix(symmetric = TRUE) allocates O(n^2) workspace, so the
  # sweep stops at 40000.)
  for (n in c(2500L, 5000L, 10000L, 20000L, 40000L)) {
    add(suite_case(sprintf("scale_n_sym_sparse_n%d_LA_k10", n), "sym_sparse", "sym", "LA", 10,
                   local({ nn <- n; function() gen_sparse_sym(nn, 5, seed = 201L) }), seed = 201L,
                   group = "scaling_n", sweep = "n", sweep_value = n))
  }
  # nnz sweep at n = 20000
  for (r in c(2, 5, 10, 20, 50)) {
    add(suite_case(sprintf("scale_nnz_sym_sparse_n20000_r%g_LA_k10", r), "sym_sparse", "sym", "LA", 10,
                   local({ rr <- r; function() gen_sparse_sym(20000L, rr, seed = 202L) }), seed = 202L,
                   group = "scaling_nnz", sweep = "nnz_per_row", sweep_value = r))
  }
  # k sweep at n = 20000
  for (k in c(2L, 5L, 10L, 20L, 40L)) {
    add(suite_case(sprintf("scale_k_sym_sparse_n20000_LA_k%d", k), "sym_sparse", "sym", "LA", k,
                   function() gen_sparse_sym(20000L, 5, seed = 203L), seed = 203L,
                   group = "scaling_k", sweep = "k", sweep_value = k))
  }
  # tall sparse SVD: long side sweep, 2000 columns, density 0.002, k=20
  for (m in c(12500L, 25000L, 50000L, 100000L, 200000L)) {
    add(suite_case(sprintf("scale_m_svd_sparse_%dx2000_k20", m), "svd_sparse", "svd", "top", 20,
                   local({ mm <- m; function() gen_sparse_rect(mm, 2000L, density = 0.002, seed = 204L) }),
                   seed = 204L, group = "scaling_svd_m", sweep = "m", sweep_value = m))
  }
  # banded smallest: 1-D Laplacian n sweep (gaps shrink like 1/n^2)
  for (n in c(2500L, 5000L, 10000L, 20000L, 40000L)) {
    add(suite_case(sprintf("scale_n_banded_lap1d_n%d_SA_k8", n), "banded", "sym", "SA", 8,
                   local({ nn <- n; function() gen_laplacian_1d(nn) }), seed = 205L,
                   group = "scaling_banded_n", sweep = "n", sweep_value = n))
  }
  out
}

suite_profile_cases <- function(profile) {
  switch(profile,
         quick = suite_core_cases("quick"),
         standard = suite_core_cases("standard"),
         scaling = suite_scaling_cases(),
         stop("unknown profile: ", profile, call. = FALSE))
}

suite_profile_defaults <- function(profile) {
  switch(profile,
         quick = list(reps = 3L, threads = "1", budget = 30),
         standard = list(reps = 5L, threads = "1,4", budget = 120),
         scaling = list(reps = 3L, threads = "1,4", budget = 60),
         list(reps = 3L, threads = "1", budget = 60))
}

suite_filter_cases <- function(cases, families = NULL, ids = NULL) {
  if (length(families)) {
    keep <- vapply(cases, function(cs) cs$family %in% families || cs$group %in% families, logical(1))
    cases <- cases[keep]
  }
  if (length(ids)) {
    keep <- vapply(cases, function(cs) any(vapply(ids, grepl, logical(1), x = cs$id, fixed = TRUE)), logical(1))
    cases <- cases[keep]
  }
  cases
}
