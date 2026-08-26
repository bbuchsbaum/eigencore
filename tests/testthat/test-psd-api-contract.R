psd_contract_text <- function() {
  path <- testthat::test_path("..", "..", "docs", "plans", "1.3-psd-api-snapshot.md")
  if (!file.exists(path)) {
    skip("repository PSD API contract is excluded from package tarballs")
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("PSD policy and tolerance schemas and defaults are frozen", {
  tolerance <- psd_tolerance()
  expect_s3_class(tolerance, "eigencore_psd_tolerance")
  expect_named(tolerance, c("schema_version", "abs", "rel"))
  expect_identical(tolerance$schema_version, 1L)
  expect_identical(tolerance$abs, 0)
  expect_identical(tolerance$rel, sqrt(.Machine$double.eps))

  policy <- psd_policy()
  expect_s3_class(policy, "eigencore_psd_policy")
  expect_named(
    policy,
    c(
      "schema_version", "scale", "symmetry", "positivity", "rank", "rhs",
      "symmetry_repair", "negative_repair", "structure_repair"
    )
  )
  expect_identical(policy$scale, "frobenius")
  expect_identical(policy$positivity$rel, 64 * .Machine$double.eps)
  expect_identical(policy$symmetry_repair, "average")
  expect_identical(policy$negative_repair, "clip")
  expect_identical(policy$structure_repair, "canonicalize")

  expect_psd_condition(
    psd_tolerance(abs = 0L),
    "eigencore_psd_invalid_input", "invalid_policy", "abs"
  )
  expect_psd_condition(
    psd_tolerance(rel = Inf),
    "eigencore_psd_invalid_input", "invalid_policy", "rel"
  )
  expect_psd_condition(
    psd_policy(scale = "spectral"),
    "eigencore_psd_invalid_input", "invalid_policy", "scale"
  )
  expect_psd_condition(
    psd_policy(negative_repair = "project"),
    "eigencore_psd_invalid_input", "invalid_policy", "negative_repair"
  )
})

test_that("the installed PSD export and S3 surface matches the frozen snapshot", {
  required_exports <- c(
    "psd_tolerance", "psd_policy", "psd_identity", "psd_factor",
    "psd_gram_factor", "psd_laplacian", "psd_capabilities", "psd_spectrum",
    "psd_rank", "psd_nullity", "psd_apply", "psd_operator", "psd_reduce",
    "psd_lift", "psd_solve", "psd_gram", "psd_orthonormalize",
    "psd_reduced_operator"
  )
  expect_setequal(
    intersect(getNamespaceExports("eigencore"), required_exports),
    required_exports
  )

  registrations <- getNamespaceInfo("eigencore", "S3methods")
  actual <- paste(registrations[, 1L], registrations[, 2L], sep = ".")
  expect_true(all(c(
    "as_operator.eigencore_psd_factor",
    "print.eigencore_psd_block_result",
    "print.eigencore_psd_capabilities",
    "print.eigencore_psd_certificate",
    "print.eigencore_psd_factor",
    "print.eigencore_psd_policy",
    "print.eigencore_psd_solve_result"
  ) %in% actual))

  contract <- psd_contract_text()
  for (name in required_exports) {
    expect_true(
      grepl(paste0(name, "("), contract, fixed = TRUE),
      info = paste("missing exact PSD signature for", name)
    )
  }
})

test_that("factor, spectrum, evidence, capability, and certificate fields are frozen", {
  factor <- psd_factor(c(9, 4, 0))
  expect_s3_class(factor, "eigencore_psd_factor")
  expect_named(
    factor,
    c(
      "schema_version", "dim", "dtype", "representation", "method", "policy",
      "scale", "thresholds", "spectrum", "classification", "rank", "nullity",
      "algebraic_nullity", "rank_bounds", "evidence", "capabilities",
      "operator_identity", "source_identity", "source_semantics",
      "materialization", "certificate", "work", "memory", "serialization",
      "warnings"
    )
  )
  expect_named(
    factor$spectrum,
    c(
      "coverage", "original", "repaired", "category", "source_index",
      "retained_indices", "exact_zero_indices", "accepted_negative_indices",
      "numerical_null_indices", "signed_zero_count", "lower_bound", "upper_bound"
    )
  )
  expect_named(
    factor$classification,
    c(
      "retained_positive", "exact_zero", "accepted_negative", "numerical_null",
      "materially_negative", "user_truncated"
    )
  )
  expect_named(
    factor$evidence,
    c(
      "schema_version", "spectrum_coverage", "validation", "action_fidelity",
      "source_semantics", "theorem", "bound_type", "details"
    )
  )
  capability_names <- c(
    "schema_version", "representation", "evidence", "form", "sqrt",
    "inverse_sqrt", "pseudoinverse", "image_projector", "null_projector",
    "reduction", "lift", "strict_solve", "gram", "orthonormalize",
    "reduced_operator", "numerical_rank", "numerical_nullity",
    "algebraic_nullity", "serialization", "cache_reuse"
  )
  expect_named(factor$capabilities, capability_names)
  for (name in capability_names[-seq_len(3L)]) {
    expect_named(
      factor$capabilities[[name]],
      c("available", "fidelity", "evidence_required", "materialization", "reason")
    )
  }

  cert <- certificate(factor)
  expect_s3_class(cert, "eigencore_psd_certificate")
  expect_s3_class(cert, "eigencore_certificate")
  expect_named(
    cert,
    c(
      "passed", "tolerance", "orthogonality_tolerance", "orthogonality_required",
      "certificate_type", "norm_bound_type", "scale_is_estimate",
      "max_backward_error", "max_residual", "max_orthogonality_loss",
      "orthogonality_passed", "failed_indices", "scale", "notes", "residuals",
      "backward_error", "orthogonality", "converged", "schema_version", "scope",
      "thresholds", "source_identity", "factor_identity", "representation",
      "evidence", "repair_applied", "symmetry_defect", "repair_defect",
      "source_action_defect", "original_spectrum", "repaired_spectrum",
      "classification", "rank", "nullity", "algebraic_nullity", "rank_bounds",
      "capabilities", "action_bounds"
    )
  )
  expect_false("exact" %in% names(factor))
  expect_false("source_exact" %in% names(cert))
})

test_that("strict-solve and block result fields are frozen", {
  factor <- psd_factor(c(9, 4, 0))
  solve_result <- psd_solve(factor, c(9, 8, 0))
  expect_named(
    solve_result,
    c(
      "schema_version", "solution", "compatibility_defect",
      "compatibility_threshold", "compatible", "factor_identity", "certificate",
      "work", "warnings"
    )
  )
  block_result <- psd_orthonormalize(factor, diag(3), required_rank = 2)
  expect_named(
    block_result,
    c(
      "schema_version", "basis", "rank", "required_rank", "condition", "dropped",
      "gram", "postcondition_error", "factor_identity", "certificate", "work",
      "warnings"
    )
  )
})
