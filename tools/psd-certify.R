#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "manifest"

suppressPackageStartupMessages({
  library(eigencore)
  library(Matrix)
})

emit <- function(...) cat(paste0(...), "\n", sep = "")

profile_allocations <- function(expr) {
  path <- tempfile("eigencore-psd-rprofmem-", fileext = ".out")
  Rprofmem(path)
  on.exit(Rprofmem(NULL), add = TRUE)
  value <- force(expr)
  Rprofmem(NULL)
  lines <- readLines(path, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  list(
    value = value,
    allocation_bytes = sum(bytes[is.finite(bytes)]),
    allocation_events = sum(is.finite(bytes))
  )
}

emit_manifest <- function() {
  constructors <- list(
    identity = psd_identity(4),
    diagonal = psd_factor(c(4, 1, 0)),
    dense_spectral = psd_factor(diag(c(4, 1, 0))),
    gram_dense = psd_gram_factor(matrix(c(1, 0, 0, 1, 1, 1), 3, 2, byrow = TRUE)),
    gram_sparse = psd_gram_factor(as(
      Matrix(matrix(c(1, 0, 0, 1, 1, 1), 3, 2, byrow = TRUE), sparse = TRUE),
      "dgCMatrix"
    )),
    laplacian_sparse = psd_laplacian(forceSymmetric(sparseMatrix(
      i = c(1, 2, 2, 3, 3),
      j = c(1, 1, 2, 2, 3),
      x = c(1, -1, 2, -1, 1),
      dims = c(3, 3)
    ), uplo = "L"))
  )
  for (name in names(constructors)) {
    factor <- constructors[[name]]
    caps <- psd_capabilities(factor)
    available <- names(Filter(function(x) is.list(x) && isTRUE(x$available), caps))
    emit(
      "manifest representation=", name,
      " certificate=", certificate(factor)$passed,
      " coverage=", factor$evidence$spectrum_coverage,
      " fidelity=", factor$evidence$action_fidelity,
      " dense_n_by_n=", factor$materialization$dense_n_by_n,
      " retained_bytes=", as.numeric(retained_bytes(factor)),
      " capabilities=", paste(available, collapse = ",")
    )
  }
}

resource_sparse_gram <- function(n = 50000L, width = 32L, block_width = 4L) {
  factor_matrix <- sparseMatrix(
    i = seq_len(n),
    j = rep(seq_len(width), length.out = n),
    x = rep(c(1, -0.5, 0.25, -0.125), length.out = n),
    dims = c(n, width)
  )
  factor_matrix <- as(factor_matrix, "dgCMatrix")
  profiled <- profile_allocations(psd_gram_factor(factor_matrix))
  factor <- profiled$value
  block <- matrix(rep(seq_len(block_width), each = n), n, block_width)
  storage.mode(block) <- "double"
  elapsed <- system.time(result <- psd_apply(factor, block, "form"))[["elapsed"]]
  emit(
    "resource route=gram_sparse n=", n,
    " width=", width,
    " block_width=", block_width,
    " source_nnz=", length(factor_matrix@x),
    " source_bytes=", as.numeric(object.size(factor_matrix)),
    " retained_bytes=", as.numeric(retained_bytes(factor)),
    " allocation_bytes=", profiled$allocation_bytes,
    " allocation_events=", profiled$allocation_events,
    " dense_equivalent_bytes=", as.double(n) * as.double(n) * 8,
    " action_seconds=", format(elapsed, digits = 8),
    " checksum=", format(sum(result), digits = 17),
    " certificate=", certificate(factor)$passed,
    " dense_n_by_n=", factor$materialization$dense_n_by_n
  )
}

resource_laplacian <- function(n = 50000L, block_width = 4L) {
  laplacian <- sparseMatrix(
    i = c(seq_len(n), seq_len(n - 1L), seq_len(n - 1L) + 1L),
    j = c(seq_len(n), seq_len(n - 1L) + 1L, seq_len(n - 1L)),
    x = c(c(1, rep(2, n - 2L), 1), rep(-1, 2L * (n - 1L))),
    dims = c(n, n)
  )
  laplacian <- as(laplacian, "dgCMatrix")
  profiled <- profile_allocations(psd_laplacian(laplacian))
  factor <- profiled$value
  block <- matrix(rep(seq_len(block_width), each = n), n, block_width)
  storage.mode(block) <- "double"
  elapsed <- system.time(result <- psd_apply(factor, block, "form"))[["elapsed"]]
  emit(
    "resource route=laplacian_sparse n=", n,
    " block_width=", block_width,
    " source_nnz=", length(laplacian@x),
    " source_bytes=", as.numeric(object.size(laplacian)),
    " retained_bytes=", as.numeric(retained_bytes(factor)),
    " allocation_bytes=", profiled$allocation_bytes,
    " allocation_events=", profiled$allocation_events,
    " dense_equivalent_bytes=", as.double(n) * as.double(n) * 8,
    " action_seconds=", format(elapsed, digits = 8),
    " checksum=", format(sum(result), digits = 17),
    " algebraic_nullity=", psd_nullity(factor, "algebraic"),
    " certificate=", certificate(factor)$passed,
    " dense_n_by_n=", factor$materialization$dense_n_by_n
  )
}

benchmark_routes <- function() {
  sizes <- c(32L, 64L, 128L, 256L)
  repetitions <- 25L
  for (n in sizes) {
    values <- seq(2, 1, length.out = n)
    block <- matrix(sin(seq_len(n * 8L)), n, 8L)
    diagonal_setup <- system.time(diagonal <- psd_factor(values))[["elapsed"]]
    dense_source <- diag(values)
    dense_setup <- system.time(dense <- psd_factor(dense_source))[["elapsed"]]
    diagonal_action <- system.time(for (i in seq_len(repetitions)) {
      psd_apply(diagonal, block, "form")
    })[["elapsed"]]
    dense_action <- system.time(for (i in seq_len(repetitions)) {
      psd_apply(dense, block, "form")
    })[["elapsed"]]
    rebuild_action <- system.time(for (i in seq_len(repetitions)) {
      psd_apply(psd_factor(dense_source), block, "form")
    })[["elapsed"]]
    emit(
      "benchmark n=", n,
      " repetitions=", repetitions,
      " diagonal_setup_seconds=", format(diagonal_setup, digits = 8),
      " dense_setup_seconds=", format(dense_setup, digits = 8),
      " diagonal_action_seconds=", format(diagonal_action, digits = 8),
      " dense_action_seconds=", format(dense_action, digits = 8),
      " dense_rebuild_action_seconds=", format(rebuild_action, digits = 8),
      " diagonal_retained_bytes=", as.numeric(retained_bytes(diagonal)),
      " dense_retained_bytes=", as.numeric(retained_bytes(dense))
    )
  }
}

switch(
  mode,
  manifest = emit_manifest(),
  resource_sparse_gram = resource_sparse_gram(),
  resource_laplacian = resource_laplacian(),
  benchmark = benchmark_routes(),
  stop("unknown certification mode: ", mode, call. = FALSE)
)
