#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" SEXP eigencore_dense_block_apply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_complex_block_apply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_randomized_apply(SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_randomized_sketch(SEXP, SEXP);
extern "C" SEXP eigencore_dense_randomized_project_transposed(SEXP, SEXP);
extern "C" SEXP eigencore_dense_randomized_svd_controller(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_block_apply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_randomized_apply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_randomized_sketch(SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_randomized_project_transposed(SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_randomized_svd_controller(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_centered_block_apply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_diagonal_block_apply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_native_apply_noalloc_check(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_apply_int_guard_check();
extern "C" SEXP eigencore_col_norms(SEXP);
extern "C" SEXP eigencore_mgs2(SEXP, SEXP);
extern "C" SEXP eigencore_cholqr2(SEXP);
extern "C" SEXP eigencore_b_cholqr2(SEXP, SEXP);
extern "C" SEXP eigencore_diagonal_b_cholqr2(SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_reorthogonalize_against(SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_basis_workspace_create(SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_basis_workspace_info(SEXP);
extern "C" SEXP eigencore_reorthogonalize_against_workspace(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lanczos_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_shift_invert_lanczos_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_shift_invert_lanczos_tridiagonal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_shift_invert_lanczos_tridiagonal_generalized(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_shift_invert_lanczos_dense_generalized(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lanczos_csc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_golub_kahan_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_golub_kahan_csc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_golub_kahan_r_operator(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_golub_kahan_dense_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_golub_kahan_csc_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_irlba_lbd_dense_retained(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_irlba_lbd_csc_retained(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_dense_basis(SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_dense_basis_cached(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_csc_basis(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_csc_basis_cached(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_dense_fit(SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_dense_fit_cached(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_csc_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_csc_fit_cached(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_dense_retained_cycle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_csc_retained_cycle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_lanczos_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_lanczos_csc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_thick_restart_lanczos_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_thick_restart_lanczos_csc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_block_thick_restart_lanczos_r_operator(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_normal_thick_restart_lanczos_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_normal_thick_restart_lanczos_csc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_dense_dense_b(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_dense_diagonal_b(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_dense_csc_b(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_csc_diagonal_b(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_csc_csc_b(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_diagonal_diagonal_b(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_dense_operator_b(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_csc_operator_b(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_diagonal_operator_b(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_csc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_lobpcg_csc_shifted_tridiagonal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_arnoldi_dense_cycle(SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_arnoldi_csc_cycle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_arnoldi_r_operator_cycle(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_arnoldi_ritz(SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_arnoldi_refined_ritz(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_orthogonality_loss(SEXP, SEXP);
extern "C" SEXP eigencore_dense_eigen_residuals(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_eigen_certificate(SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_svd_residuals(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_svd_certificate(SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_svd_certificate_cached_av(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_eigen_certificate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_diagonal_eigen_certificate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_tridiagonal_eigen_certificate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_svd_certificate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_svd_certificate_cached_av(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_diagonal_svd_certificate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_diagonal_svd_certificate_cached_av(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_rayleigh_ritz_symmetric(SEXP, SEXP);
extern "C" SEXP eigencore_tridiagonal_eigen(SEXP, SEXP);
extern "C" SEXP eigencore_tridiagonal_eigen_selected(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_bidiagonal_svd(SEXP, SEXP);
extern "C" SEXP eigencore_block_golub_kahan_ritz(SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_golub_kahan_ritz(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_is_symmetric(SEXP, SEXP);
extern "C" SEXP eigencore_dense_symmetric_eigen(SEXP);
extern "C" SEXP eigencore_dense_symmetric_eigen_dsyevd(SEXP);
extern "C" SEXP eigencore_dense_complex_hermitian_eigen(SEXP);
extern "C" SEXP eigencore_dense_complex_general_eigen(SEXP);
extern "C" SEXP eigencore_dense_generalized_pencil_eigen(SEXP, SEXP);
extern "C" SEXP eigencore_dense_complex_generalized_hpd_eigen(SEXP, SEXP);
extern "C" SEXP eigencore_dense_complex_generalized_pencil_eigen(SEXP, SEXP);
extern "C" SEXP eigencore_dense_generalized_schur(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_complex_generalized_schur(SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_symmetric_eigen_selected(SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_symmetric_eigen_dsyevx_selected(SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_left_gram_svd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_right_gram_svd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_left_gram_svd_fast_result(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_csc_right_gram_svd_fast_result(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP eigencore_dense_generalized_spd_eigen(SEXP, SEXP);
extern "C" SEXP eigencore_dense_svd(SEXP);
extern "C" SEXP eigencore_dense_complex_svd(SEXP);
extern "C" SEXP eigencore_dense_generalized_svd(SEXP, SEXP);
extern "C" SEXP eigencore_dense_complex_generalized_svd(SEXP, SEXP);
extern "C" SEXP eigencore_tridiagonal_solve(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"eigencore_dense_block_apply", (DL_FUNC) &eigencore_dense_block_apply, 6},
  {"eigencore_dense_complex_block_apply", (DL_FUNC) &eigencore_dense_complex_block_apply, 6},
  {"eigencore_dense_randomized_apply", (DL_FUNC) &eigencore_dense_randomized_apply, 3},
  {"eigencore_dense_randomized_sketch", (DL_FUNC) &eigencore_dense_randomized_sketch, 2},
  {"eigencore_dense_randomized_project_transposed", (DL_FUNC) &eigencore_dense_randomized_project_transposed, 2},
  {"eigencore_dense_randomized_svd_controller", (DL_FUNC) &eigencore_dense_randomized_svd_controller, 6},
  {"eigencore_csc_block_apply", (DL_FUNC) &eigencore_csc_block_apply, 9},
  {"eigencore_csc_randomized_apply", (DL_FUNC) &eigencore_csc_randomized_apply, 6},
  {"eigencore_csc_randomized_sketch", (DL_FUNC) &eigencore_csc_randomized_sketch, 5},
  {"eigencore_csc_randomized_project_transposed", (DL_FUNC) &eigencore_csc_randomized_project_transposed, 5},
  {"eigencore_csc_randomized_svd_controller", (DL_FUNC) &eigencore_csc_randomized_svd_controller, 9},
  {"eigencore_csc_centered_block_apply", (DL_FUNC) &eigencore_csc_centered_block_apply, 13},
  {"eigencore_diagonal_block_apply", (DL_FUNC) &eigencore_diagonal_block_apply, 7},
  {"eigencore_native_apply_noalloc_check", (DL_FUNC) &eigencore_native_apply_noalloc_check, 4},
  {"eigencore_dense_apply_int_guard_check", (DL_FUNC) &eigencore_dense_apply_int_guard_check, 0},
  {"eigencore_col_norms", (DL_FUNC) &eigencore_col_norms, 1},
  {"eigencore_mgs2", (DL_FUNC) &eigencore_mgs2, 2},
  {"eigencore_cholqr2", (DL_FUNC) &eigencore_cholqr2, 1},
  {"eigencore_b_cholqr2", (DL_FUNC) &eigencore_b_cholqr2, 2},
  {"eigencore_diagonal_b_cholqr2", (DL_FUNC) &eigencore_diagonal_b_cholqr2, 3},
  {"eigencore_reorthogonalize_against", (DL_FUNC) &eigencore_reorthogonalize_against, 3},
  {"eigencore_basis_workspace_create", (DL_FUNC) &eigencore_basis_workspace_create, 3},
  {"eigencore_basis_workspace_info", (DL_FUNC) &eigencore_basis_workspace_info, 1},
  {"eigencore_reorthogonalize_against_workspace", (DL_FUNC) &eigencore_reorthogonalize_against_workspace, 4},
  {"eigencore_lanczos_dense", (DL_FUNC) &eigencore_lanczos_dense, 6},
  {"eigencore_shift_invert_lanczos_dense", (DL_FUNC) &eigencore_shift_invert_lanczos_dense, 7},
  {"eigencore_shift_invert_lanczos_tridiagonal", (DL_FUNC) &eigencore_shift_invert_lanczos_tridiagonal, 8},
  {"eigencore_shift_invert_lanczos_tridiagonal_generalized", (DL_FUNC) &eigencore_shift_invert_lanczos_tridiagonal_generalized, 9},
  {"eigencore_shift_invert_lanczos_dense_generalized", (DL_FUNC) &eigencore_shift_invert_lanczos_dense_generalized, 8},
  {"eigencore_lanczos_csc", (DL_FUNC) &eigencore_lanczos_csc, 9},
  {"eigencore_golub_kahan_dense", (DL_FUNC) &eigencore_golub_kahan_dense, 7},
  {"eigencore_golub_kahan_csc", (DL_FUNC) &eigencore_golub_kahan_csc, 10},
  {"eigencore_golub_kahan_r_operator", (DL_FUNC) &eigencore_golub_kahan_r_operator, 11},
  {"eigencore_golub_kahan_dense_fit", (DL_FUNC) &eigencore_golub_kahan_dense_fit, 9},
  {"eigencore_golub_kahan_csc_fit", (DL_FUNC) &eigencore_golub_kahan_csc_fit, 12},
  {"eigencore_irlba_lbd_dense_retained", (DL_FUNC) &eigencore_irlba_lbd_dense_retained, 14},
  {"eigencore_irlba_lbd_csc_retained", (DL_FUNC) &eigencore_irlba_lbd_csc_retained, 17},
  {"eigencore_block_golub_kahan_dense_basis", (DL_FUNC) &eigencore_block_golub_kahan_dense_basis, 3},
  {"eigencore_block_golub_kahan_dense_basis_cached", (DL_FUNC) &eigencore_block_golub_kahan_dense_basis_cached, 4},
  {"eigencore_block_golub_kahan_csc_basis", (DL_FUNC) &eigencore_block_golub_kahan_csc_basis, 6},
  {"eigencore_block_golub_kahan_csc_basis_cached", (DL_FUNC) &eigencore_block_golub_kahan_csc_basis_cached, 7},
  {"eigencore_block_golub_kahan_dense_fit", (DL_FUNC) &eigencore_block_golub_kahan_dense_fit, 5},
  {"eigencore_block_golub_kahan_dense_fit_cached", (DL_FUNC) &eigencore_block_golub_kahan_dense_fit_cached, 6},
  {"eigencore_block_golub_kahan_csc_fit", (DL_FUNC) &eigencore_block_golub_kahan_csc_fit, 8},
  {"eigencore_block_golub_kahan_csc_fit_cached", (DL_FUNC) &eigencore_block_golub_kahan_csc_fit_cached, 9},
  {"eigencore_block_golub_kahan_dense_retained_cycle", (DL_FUNC) &eigencore_block_golub_kahan_dense_retained_cycle, 11},
  {"eigencore_block_golub_kahan_csc_retained_cycle", (DL_FUNC) &eigencore_block_golub_kahan_csc_retained_cycle, 14},
  {"eigencore_block_lanczos_dense", (DL_FUNC) &eigencore_block_lanczos_dense, 7},
  {"eigencore_block_lanczos_csc", (DL_FUNC) &eigencore_block_lanczos_csc, 10},
  {"eigencore_block_thick_restart_lanczos_dense", (DL_FUNC) &eigencore_block_thick_restart_lanczos_dense, 10},
  {"eigencore_block_thick_restart_lanczos_csc", (DL_FUNC) &eigencore_block_thick_restart_lanczos_csc, 13},
  {"eigencore_block_thick_restart_lanczos_r_operator", (DL_FUNC) &eigencore_block_thick_restart_lanczos_r_operator, 11},
  {"eigencore_normal_thick_restart_lanczos_dense", (DL_FUNC) &eigencore_normal_thick_restart_lanczos_dense, 10},
  {"eigencore_normal_thick_restart_lanczos_csc", (DL_FUNC) &eigencore_normal_thick_restart_lanczos_csc, 13},
  {"eigencore_lobpcg_dense", (DL_FUNC) &eigencore_lobpcg_dense, 10},
  {"eigencore_lobpcg_dense_dense_b", (DL_FUNC) &eigencore_lobpcg_dense_dense_b, 11},
  {"eigencore_lobpcg_dense_diagonal_b", (DL_FUNC) &eigencore_lobpcg_dense_diagonal_b, 12},
  {"eigencore_lobpcg_dense_csc_b", (DL_FUNC) &eigencore_lobpcg_dense_csc_b, 14},
  {"eigencore_lobpcg_csc_diagonal_b", (DL_FUNC) &eigencore_lobpcg_csc_diagonal_b, 15},
  {"eigencore_lobpcg_csc_csc_b", (DL_FUNC) &eigencore_lobpcg_csc_csc_b, 17},
  {"eigencore_lobpcg_diagonal_diagonal_b", (DL_FUNC) &eigencore_lobpcg_diagonal_diagonal_b, 14},
  {"eigencore_lobpcg_dense_operator_b", (DL_FUNC) &eigencore_lobpcg_dense_operator_b, 11},
  {"eigencore_lobpcg_csc_operator_b", (DL_FUNC) &eigencore_lobpcg_csc_operator_b, 14},
  {"eigencore_lobpcg_diagonal_operator_b", (DL_FUNC) &eigencore_lobpcg_diagonal_operator_b, 13},
  {"eigencore_lobpcg_csc", (DL_FUNC) &eigencore_lobpcg_csc, 13},
  {"eigencore_lobpcg_csc_shifted_tridiagonal", (DL_FUNC) &eigencore_lobpcg_csc_shifted_tridiagonal, 10},
  {"eigencore_arnoldi_dense_cycle", (DL_FUNC) &eigencore_arnoldi_dense_cycle, 3},
  {"eigencore_arnoldi_csc_cycle", (DL_FUNC) &eigencore_arnoldi_csc_cycle, 6},
  {"eigencore_arnoldi_r_operator_cycle", (DL_FUNC) &eigencore_arnoldi_r_operator_cycle, 4},
  {"eigencore_arnoldi_ritz", (DL_FUNC) &eigencore_arnoldi_ritz, 3},
  {"eigencore_arnoldi_refined_ritz", (DL_FUNC) &eigencore_arnoldi_refined_ritz, 4},
  {"eigencore_orthogonality_loss", (DL_FUNC) &eigencore_orthogonality_loss, 2},
  {"eigencore_dense_eigen_residuals", (DL_FUNC) &eigencore_dense_eigen_residuals, 4},
  {"eigencore_dense_eigen_certificate", (DL_FUNC) &eigencore_dense_eigen_certificate, 5},
  {"eigencore_dense_svd_residuals", (DL_FUNC) &eigencore_dense_svd_residuals, 4},
  {"eigencore_dense_svd_certificate", (DL_FUNC) &eigencore_dense_svd_certificate, 5},
  {"eigencore_dense_svd_certificate_cached_av", (DL_FUNC) &eigencore_dense_svd_certificate_cached_av, 6},
  {"eigencore_csc_eigen_certificate", (DL_FUNC) &eigencore_csc_eigen_certificate, 8},
  {"eigencore_diagonal_eigen_certificate", (DL_FUNC) &eigencore_diagonal_eigen_certificate, 7},
  {"eigencore_tridiagonal_eigen_certificate", (DL_FUNC) &eigencore_tridiagonal_eigen_certificate, 6},
  {"eigencore_csc_svd_certificate", (DL_FUNC) &eigencore_csc_svd_certificate, 9},
  {"eigencore_csc_svd_certificate_cached_av", (DL_FUNC) &eigencore_csc_svd_certificate_cached_av, 10},
  {"eigencore_diagonal_svd_certificate", (DL_FUNC) &eigencore_diagonal_svd_certificate, 8},
  {"eigencore_diagonal_svd_certificate_cached_av", (DL_FUNC) &eigencore_diagonal_svd_certificate_cached_av, 9},
  {"eigencore_rayleigh_ritz_symmetric", (DL_FUNC) &eigencore_rayleigh_ritz_symmetric, 2},
  {"eigencore_tridiagonal_eigen", (DL_FUNC) &eigencore_tridiagonal_eigen, 2},
  {"eigencore_tridiagonal_eigen_selected", (DL_FUNC) &eigencore_tridiagonal_eigen_selected, 4},
  {"eigencore_bidiagonal_svd", (DL_FUNC) &eigencore_bidiagonal_svd, 2},
  {"eigencore_block_golub_kahan_ritz", (DL_FUNC) &eigencore_block_golub_kahan_ritz, 5},
  {"eigencore_golub_kahan_ritz", (DL_FUNC) &eigencore_golub_kahan_ritz, 7},
  {"eigencore_dense_is_symmetric", (DL_FUNC) &eigencore_dense_is_symmetric, 2},
  {"eigencore_dense_symmetric_eigen", (DL_FUNC) &eigencore_dense_symmetric_eigen, 1},
  {"eigencore_dense_symmetric_eigen_dsyevd", (DL_FUNC) &eigencore_dense_symmetric_eigen_dsyevd, 1},
  {"eigencore_dense_complex_hermitian_eigen", (DL_FUNC) &eigencore_dense_complex_hermitian_eigen, 1},
  {"eigencore_dense_complex_general_eigen", (DL_FUNC) &eigencore_dense_complex_general_eigen, 1},
  {"eigencore_dense_generalized_pencil_eigen", (DL_FUNC) &eigencore_dense_generalized_pencil_eigen, 2},
  {"eigencore_dense_complex_generalized_hpd_eigen", (DL_FUNC) &eigencore_dense_complex_generalized_hpd_eigen, 2},
  {"eigencore_dense_complex_generalized_pencil_eigen", (DL_FUNC) &eigencore_dense_complex_generalized_pencil_eigen, 2},
  {"eigencore_dense_generalized_schur", (DL_FUNC) &eigencore_dense_generalized_schur, 4},
  {"eigencore_dense_complex_generalized_schur", (DL_FUNC) &eigencore_dense_complex_generalized_schur, 4},
  {"eigencore_dense_symmetric_eigen_selected", (DL_FUNC) &eigencore_dense_symmetric_eigen_selected, 3},
  {"eigencore_dense_symmetric_eigen_dsyevx_selected", (DL_FUNC) &eigencore_dense_symmetric_eigen_dsyevx_selected, 3},
  {"eigencore_csc_left_gram_svd", (DL_FUNC) &eigencore_csc_left_gram_svd, 6},
  {"eigencore_csc_right_gram_svd", (DL_FUNC) &eigencore_csc_right_gram_svd, 6},
  {"eigencore_csc_left_gram_svd_fast_result", (DL_FUNC) &eigencore_csc_left_gram_svd_fast_result, 6},
  {"eigencore_csc_right_gram_svd_fast_result", (DL_FUNC) &eigencore_csc_right_gram_svd_fast_result, 6},
  {"eigencore_dense_generalized_spd_eigen", (DL_FUNC) &eigencore_dense_generalized_spd_eigen, 2},
  {"eigencore_dense_svd", (DL_FUNC) &eigencore_dense_svd, 1},
  {"eigencore_dense_complex_svd", (DL_FUNC) &eigencore_dense_complex_svd, 1},
  {"eigencore_dense_generalized_svd", (DL_FUNC) &eigencore_dense_generalized_svd, 2},
  {"eigencore_dense_complex_generalized_svd", (DL_FUNC) &eigencore_dense_complex_generalized_svd, 2},
  {"eigencore_tridiagonal_solve", (DL_FUNC) &eigencore_tridiagonal_solve, 4},
  {NULL, NULL, 0}
};

extern "C" void R_init_eigencore(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
