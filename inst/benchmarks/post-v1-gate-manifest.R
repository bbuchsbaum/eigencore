post_v1_gate_manifest <- list(
  version = 1L,
  generated_on = "2026-06-06",
  purpose = paste(
    "Post-V1 benchmark truth surface for solver promotion decisions.",
    "Rows here are gates to run before broadening public planner labels beyond",
    "the scoped V1 surfaces."
  ),
  installed_package_prefix = paste(
    "R CMD INSTALL --library=/tmp/eigencore-bench-lib .",
    "R_LIBS=/tmp/eigencore-bench-lib Rscript"
  ),
  artifact_policy = list(
    directory = "inst/benchmarks/results",
    naming = "YYYYMMDD-<gate>-<artifact>.rds",
    required_artifacts = c("rows", "gates", "contracts", "memory"),
    release_record = "benchmarks/RELEASES.md"
  ),
  tier_profile = list(
    smoke = list(
      default_gate_ids = c("post_v1_operator_sidecars"),
      description = "Fast local/CI boundary-truth smoke. Validates the manifest and runs the matrix-free operator sidecar strict gate."
    ),
    strict = list(
      default_gate_ids = "all",
      description = "Installed-package strict release profile for promoted and candidate gates."
    ),
    long = list(
      default_gate_ids = "all",
      description = "Nightly/weekly long benchmark profile. Uses the strict command with expanded iteration counts unless a gate supplies a long_command."
    )
  ),
  required_metrics = c(
    "time_to_certified_answer",
    "memory",
    "max_residual",
    "max_backward_error",
    "orthogonality_loss",
    "planner_label",
    "certificate_type",
    "native_or_reference_boundary",
    "nconv"
  ),
  current_gate_owner_issue_ids = c(
    "bd-01KTE8G6RYE4RD5F6CN7SNKKC6",
    "bd-01KVWKVFC8QFZSTR6KFQ3MM3NH"
  ),
  gates = list(
    list(
      id = "post_v1_svd_hard_surface",
      owner_issue = "bd-01KTE8G6RYE4RD5F6CN7SNKKC6",
      surface = "sparse_and_hard_svd",
      script = "inst/benchmarks/bench-svd-surface.R",
      command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-svd-surface.R",
        "--iterations=3 --h-candidate",
        "--methods=eigencore,RSpectra,PRIMME,irlba",
        "--cases=tall_sparse,wide_sparse,rank_deficient_sparse,clustered_dense,slow_decay_dense,low_rank_sparse",
        "--subject=eigencore --strict --save"
      ),
      quick_smoke_command = paste(
        "Rscript inst/benchmarks/bench-svd-surface.R",
        "--quick --iterations=1 --h-candidate",
        "--methods=eigencore,RSpectra,irlba",
        "--cases=tall_sparse:600x90,wide_sparse:90x600,rank_deficient_sparse:160x50",
        "--subject=eigencore --strict"
      ),
      long_command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-svd-surface.R",
        "--iterations=10 --h-candidate",
        "--methods=eigencore,RSpectra,PRIMME,irlba",
        "--cases=tall_sparse,wide_sparse,rank_deficient_sparse,clustered_dense,slow_decay_dense,low_rank_sparse",
        "--subject=eigencore --strict --save"
      ),
      cases = c(
        "tall_sparse:100000x500",
        "wide_sparse:500x100000",
        "rank_deficient_sparse:5000x500",
        "clustered_dense:2000x500",
        "slow_decay_dense:2000x500",
        "low_rank_sparse:10000x500"
      ),
      baselines = c("RSpectra", "PRIMME", "irlba", "base_LAPACK_small"),
      artifacts = c(
        "post-v1-svd-hard-surface-rows.rds",
        "post-v1-svd-hard-surface-gates.rds",
        "post-v1-svd-hard-surface-memory.rds"
      ),
      thresholds = list(
        speed_ratio_min = 1.10,
        memory_ratio_min = 1.00,
        certificate_passed = TRUE,
        nconv_equals_requested = TRUE,
        max_backward_error_lte_tol = TRUE,
        materialized_normal_equation_allowed = FALSE,
        planner_label_must_not_match = c("reference", "prototype", "oracle")
      )
    ),
    list(
      id = "post_v1_operator_sidecars",
      owner_issue = "bd-01KTE8G6RYE4RD5F6CN7SNKKC6",
      surface = "matrix_free_boundary_truth",
      script = "inst/benchmarks/bench-post-v1-operator-sidecars.R",
      command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-post-v1-operator-sidecars.R",
        "--iterations=3 --strict --save"
      ),
      quick_smoke_command = paste(
        "Rscript inst/benchmarks/bench-nonsymmetric.R",
        "--quick --iterations=1 --strict"
      ),
      long_command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-post-v1-operator-sidecars.R",
        "--iterations=10 --strict --save"
      ),
      cases = c(
        "matrix_free_svd:80x5",
        "matrix_free_nonnormal:30",
        "matrix_free_b:diagonal6"
      ),
      baselines = c("current_planner_label", "exact_certificate"),
      artifacts = c("post-v1-operator-sidecars-rows.rds"),
      thresholds = list(
        certificate_passed = TRUE,
        max_backward_error_lte_tol = TRUE,
        planner_label_exact = TRUE,
        native_boundary_exact = TRUE
      )
    ),
    list(
      id = "post_v1_randomized_svd_hard_surface",
      owner_issue = "bd-01KTE8G6RYE4RD5F6CN7SNKKC6",
      surface = "randomized_svd_slow_decay_and_sparse",
      script = "inst/benchmarks/bench-randomized-rsvd.R",
      command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-randomized-rsvd.R",
        "--iterations=3 --methods=eigencore_randomized,rsvd,irlba",
        "--cases=exact_low_rank_dense,nearly_low_rank_dense,slow_decay_dense,low_rank_sparse",
        "--strict --save"
      ),
      quick_smoke_command = paste(
        "Rscript inst/benchmarks/bench-randomized-rsvd.R",
        "--quick --iterations=1 --methods=eigencore_randomized,rsvd",
        "--strict"
      ),
      long_command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-randomized-rsvd.R",
        "--iterations=10 --methods=eigencore_randomized,rsvd,irlba",
        "--cases=exact_low_rank_dense,nearly_low_rank_dense,slow_decay_dense,low_rank_sparse",
        "--strict --save"
      ),
      cases = c(
        "exact_low_rank_dense:2000x500",
        "nearly_low_rank_dense",
        "slow_decay_dense",
        "low_rank_sparse"
      ),
      baselines = c("rsvd", "irlba"),
      artifacts = c(
        "post-v1-randomized-rsvd-rows.rds",
        "post-v1-randomized-rsvd-gates.rds"
      ),
      thresholds = list(
        speed_ratio_min = 2.00,
        certificate_passed = TRUE,
        baseline_certified_required = TRUE,
        subspace_error_no_worse_than_baseline = TRUE,
        planner_label_must_not_match = c("reference-control")
      )
    ),
    list(
      id = "post_v1_generalized_preconditioned_surface",
      owner_issue = "bd-01KTE8G6RYE4RD5F6CN7SNKKC6",
      surface = "generalized_spd_preconditioned",
      script = "inst/benchmarks/bench-generalized-lobpcg.R",
      command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-generalized-lobpcg.R",
        "--iterations=3 --strict --save",
        "--cases=sparse_generalized_path_smallest:1000,sparse_generalized_path_largest:1000,adversarial_explicit_spd_matrix_free_b_smallest,adversarial_ill_conditioned_diagonal_b_smallest,adversarial_sparse_csc_b_largest",
        "--methods=eigencore_shifted_tridiagonal,eigencore,eigencore_lanczos_reference,base",
        "--subject=eigencore_shifted_tridiagonal"
      ),
      quick_smoke_command = paste(
        "Rscript inst/benchmarks/bench-generalized-lobpcg.R",
        "--quick --iterations=1 --strict"
      ),
      long_command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-generalized-lobpcg.R",
        "--iterations=10 --strict --save",
        "--cases=sparse_generalized_path_smallest:1000,sparse_generalized_path_largest:1000,adversarial_explicit_spd_matrix_free_b_smallest,adversarial_ill_conditioned_diagonal_b_smallest,adversarial_sparse_csc_b_largest",
        "--methods=eigencore_shifted_tridiagonal,eigencore,eigencore_lanczos_reference,base",
        "--subject=eigencore_shifted_tridiagonal"
      ),
      cases = c(
        "sparse_generalized_path_smallest:1000",
        "sparse_generalized_path_largest:1000",
        "adversarial_explicit_spd_matrix_free_b_smallest",
        "adversarial_ill_conditioned_diagonal_b_smallest",
        "adversarial_sparse_csc_b_largest"
      ),
      baselines = c("base", "eigencore_lanczos_reference"),
      artifacts = c(
        "post-v1-generalized-lobpcg-rows.rds",
        "post-v1-generalized-lobpcg-gates.rds",
        "post-v1-generalized-lobpcg-native-contracts.rds"
      ),
      thresholds = list(
        speed_ratio_min = 1.25,
        memory_ratio_min = 4.00,
        certificate_passed = TRUE,
        b_orthogonality_passed = TRUE,
        no_sparse_densification = TRUE,
        planner_label_must_not_match = c("reference", "prototype", "oracle")
      )
    ),
    list(
      id = "post_v1_generalized_eigen_surface",
      owner_issue = "bd-01KVWKVFC8QFZSTR6KFQ3MM3NH",
      surface = "generalized_eigen_replacement",
      script = "inst/benchmarks/bench-generalized-eigen.R",
      command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-generalized-eigen.R",
        "--iterations=3 --strict --save"
      ),
      quick_smoke_command = paste(
        "Rscript inst/benchmarks/bench-generalized-eigen.R",
        "--quick --iterations=1 --strict"
      ),
      long_command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-generalized-eigen.R",
        "--iterations=10 --strict --save"
      ),
      cases = c(
        "dense_full_generalized:spd_real",
        "dense_full_generalized:general_pencil_real",
        "sparse_general_pencil_partial:diagonal_B_transform",
        "sparse_general_pencil_partial:unsupported_boundary",
        "qz_dense:real_unsorted_and_sorted",
        "gsvd_dense:reconstruction"
      ),
      baselines = c("base::eigen", "dense LAPACK oracle"),
      artifacts = c(
        "post-v1-generalized-eigen-rows.rds",
        "post-v1-generalized-eigen-gates.rds"
      ),
      thresholds = list(
        certificate_passed = TRUE,
        oracle_match = TRUE,
        no_sparse_densification = TRUE,
        unsupported_boundary_fails_loud = TRUE,
        planner_label_must_not_match = c("reference-control")
      )
    ),
    list(
      id = "post_v1_shift_invert_boundaries",
      owner_issue = "bd-01KTE8G6RYE4RD5F6CN7SNKKC6",
      surface = "shift_invert_sparse_and_user_solve",
      script = "inst/benchmarks/bench-shift-invert.R",
      command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-shift-invert.R",
        "--iterations=3 --strict --save",
        "--cases=sparse_general_reference,sparse_general_diagonal_b_reference,matrix_free_user_solve_reference"
      ),
      quick_smoke_command = paste(
        "Rscript inst/benchmarks/bench-shift-invert.R",
        "--quick --iterations=1 --strict",
        "--cases=matrix_free_user_solve_reference"
      ),
      long_command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-shift-invert.R",
        "--iterations=10 --strict --save",
        "--cases=sparse_general_reference,sparse_general_diagonal_b_reference,matrix_free_user_solve_reference"
      ),
      cases = c(
        "sparse_general_reference",
        "sparse_general_diagonal_b_reference",
        "matrix_free_user_solve_reference"
      ),
      baselines = c("Matrix::lu", "user_solve", "dense_base_small"),
      artifacts = c(
        "post-v1-shift-invert-rows.rds",
        "post-v1-shift-invert-contracts.rds"
      ),
      thresholds = list(
        certificate_passed = TRUE,
        original_coordinate_certificate = TRUE,
        cache_provenance_present = TRUE,
        external_cache_contract_present = TRUE,
        native_label_required_for_promoted_rows = TRUE
      )
    ),
    list(
      id = "post_v1_nonsymmetric_matrix_free_surface",
      owner_issue = "bd-01KTE8G6RYE4RD5F6CN7SNKKC6",
      surface = "nonsymmetric_matrix_free",
      script = "inst/benchmarks/bench-nonsymmetric.R",
      command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-nonsymmetric.R",
        "--iterations=3 --strict --save"
      ),
      quick_smoke_command = paste(
        "Rscript inst/benchmarks/bench-post-v1-operator-sidecars.R",
        "--quick --iterations=1 --strict"
      ),
      long_command = paste(
        "R_LIBS=/tmp/eigencore-bench-lib Rscript",
        "inst/benchmarks/bench-nonsymmetric.R",
        "--iterations=10 --strict --save"
      ),
      cases = c(
        "matrix_free_nonnormal:30",
        "dense_native_arnoldi_lm",
        "sparse_native_arnoldi_lr"
      ),
      baselines = c("RSpectra", "base_LAPACK_small", "current_native_matrix_free_arnoldi"),
      artifacts = c(
        "post-v1-nonsymmetric-rows.rds",
        "post-v1-nonsymmetric-contracts.rds"
      ),
      thresholds = list(
        certificate_passed = TRUE,
        right_residual_certificate = TRUE,
        native_arnoldi_label_for_promoted_rows = TRUE,
        matrix_free_native_label_required = TRUE,
        restart_diagnostics_present = TRUE
      )
    )
  )
)
