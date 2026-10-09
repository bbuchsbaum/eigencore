# Benchmark suite: optional SuiteSparse Matrix Collection cases (--suitesparse).
#
# Matrices are downloaded once into a cache directory as Matrix Market
# tarballs from https://sparse.tamu.edu/MM/<group>/<name>.tar.gz. When the
# collection cannot be reached (offline CI, proxy) the cases are skipped with a
# message; the rest of the suite still runs.

SUITE_SUITESPARSE <- list(
  list(group = "HB", name = "bcsstk17", task = "sym", target = "LA", k = 10L,
       note = "structural stiffness, SPD, n=10974"),
  list(group = "HB", name = "1138_bus", task = "sym", target = "SA", k = 6L,
       note = "power network admittance, SPD, n=1138, ill-conditioned"),
  list(group = "SNAP", name = "ca-GrQc", task = "sym", target = "LA", k = 10L,
       note = "collaboration graph adjacency, n=5242"),
  list(group = "HB", name = "west2021", task = "nonsym", target = "LM", k = 6L,
       note = "chemical engineering, nonsymmetric, n=2021"),
  list(group = "HB", name = "well1850", task = "svd", target = "top", k = 10L,
       note = "least-squares, 1850 x 712")
)

suite_suitesparse_fetch <- function(group, name, cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  mtx <- file.path(cache_dir, name, paste0(name, ".mtx"))
  if (file.exists(mtx)) return(mtx)
  url <- sprintf("https://sparse.tamu.edu/MM/%s/%s.tar.gz", group, name)
  tgz <- file.path(cache_dir, paste0(name, ".tar.gz"))
  ok <- tryCatch({
    old <- options(timeout = max(120, getOption("timeout")))
    on.exit(options(old), add = TRUE)
    utils::download.file(url, tgz, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(tgz)) {
    unlink(tgz)
    return(NULL)
  }
  utils::untar(tgz, exdir = cache_dir)
  unlink(tgz)
  if (file.exists(mtx)) mtx else NULL
}

suite_suitesparse_cases <- function(cache_dir) {
  out <- list()
  for (s in SUITE_SUITESPARSE) {
    path <- suite_suitesparse_fetch(s$group, s$name, cache_dir)
    if (is.null(path)) {
      suite_log("SuiteSparse ", s$group, "/", s$name, ": unavailable (offline?), skipped")
      next
    }
    out[[length(out) + 1L]] <- local({
      p <- path
      ss <- s
      suite_case(sprintf("suitesparse_%s_%s_%s_k%d", ss$group, ss$name, ss$target, ss$k),
                 "suitesparse", ss$task, ss$target, ss$k,
                 function() {
                   A <- Matrix::readMM(p)
                   if (ss$task == "sym") A <- Matrix::forceSymmetric(A, uplo = "L")
                   list(A = as_dgc(A), description = sprintf("SuiteSparse %s/%s: %s", ss$group, ss$name, ss$note))
                 }, seed = 300L)
    })
  }
  out
}
