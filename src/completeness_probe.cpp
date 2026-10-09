// Target-completeness probe (C50): a short block Lanczos process on a
// Hermitian operator compressed to the orthogonal complement of a returned
// eigenvector block V. See R/target_completeness.R for the contract.
//
// Each step applies the operator to an orthonormal block X (X orthogonal to V
// and to every earlier block), grows the projected matrix H = Q' A Q by one
// block column, and checks the Rayleigh-Ritz values of H against the returned
// edge. Every Ritz value of H is a Rayleigh quotient of a unit vector
// orthogonal to V, so one beyond the edge (by more than the caller's margin)
// proves that the operator has a more-preferred eigenvalue outside span(V).
// Start blocks come from a private splitmix64 stream, never from R's RNG.
#include <cmath>
#include <cfloat>
#include <cstdint>
#include <cstring>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "eigencore_lapack_compat.h"
#include "eigencore_common.h"
#include "native_operators.h"

namespace {

struct ProbeRng {
  uint64_t state;
  uint64_t next() {
    uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }
  double uniform() {  // in (0, 1)
    return (static_cast<double>(next() >> 11) + 0.5) * (1.0 / 9007199254740992.0);
  }
  double normal() {
    const double u1 = uniform();
    const double u2 = uniform();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
  }
};

// Serial skinny products for the probe's tiny blocks (b <= 8 columns against
// a basis of a few dozen columns). Plain loops avoid waking a threaded BLAS
// for each of the probe's many small products, which costs more than the
// arithmetic here.
// C (c x b) = B' X.
void skinny_tn(int n, int c, int b, const double* B, const double* X, double* C) {
  for (int col = 0; col < b; ++col) {
    const double* x = X + static_cast<int64_t>(col) * n;
    for (int j = 0; j < c; ++j) {
      const double* v = B + static_cast<int64_t>(j) * n;
      double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
      int r = 0;
      for (; r + 4 <= n; r += 4) {
        s0 += v[r] * x[r];
        s1 += v[r + 1] * x[r + 1];
        s2 += v[r + 2] * x[r + 2];
        s3 += v[r + 3] * x[r + 3];
      }
      for (; r < n; ++r) s0 += v[r] * x[r];
      C[j + static_cast<int64_t>(col) * c] = (s0 + s1) + (s2 + s3);
    }
  }
}

// X (n x b) -= B (n x c) C (c x b).
void skinny_nn_minus(int n, int c, int b, const double* B, const double* C,
                     double* X) {
  for (int col = 0; col < b; ++col) {
    double* x = X + static_cast<int64_t>(col) * n;
    for (int j = 0; j < c; ++j) {
      const double cj = C[j + static_cast<int64_t>(col) * c];
      if (cj == 0.0) continue;
      const double* v = B + static_cast<int64_t>(j) * n;
      for (int r = 0; r < n; ++r) x[r] -= cj * v[r];
    }
  }
}

// X (n x b) -= B (B' X), two passes. B is n x c column-major.
void project_out_block(double* X, int n, int b, const double* B, int c,
                       std::vector<double>* coef) {
  if (c <= 0 || b <= 0) return;
  coef->assign(static_cast<size_t>(c) * static_cast<size_t>(b), 0.0);
  for (int pass = 0; pass < 2; ++pass) {
    skinny_tn(n, c, b, B, X, coef->data());
    skinny_nn_minus(n, c, b, B, coef->data(), X);
  }
}

// Two-pass MGS within the block; columns whose norm falls below
// max(1e-10 * original norm, floor_abs) are dropped. Returns the kept count
// (kept columns are compacted to the front).
int orthonormalize_block(double* X, int n, int b, const double* original_norms,
                         double floor_abs) {
  int kept = 0;
  double coef = 0.0;
  for (int j = 0; j < b; ++j) {
    double* x = X + static_cast<int64_t>(j) * n;
    for (int pass = 0; pass < 2; ++pass) {
      for (int i = 0; i < kept; ++i) {
        const double* q = X + static_cast<int64_t>(i) * n;
        skinny_tn(n, 1, 1, q, x, &coef);
        skinny_nn_minus(n, 1, 1, q, &coef, x);
      }
    }
    double nx2 = 0.0;
    skinny_tn(n, 1, 1, x, x, &nx2);
    const double nx = std::sqrt(nx2);
    const double floor_rel = 1e-10 * original_norms[j];
    if (std::isfinite(nx) && nx > floor_rel && nx > floor_abs && nx > 0.0) {
      const double inv = 1.0 / nx;
      double* dst = X + static_cast<int64_t>(kept) * n;
      for (int r = 0; r < n; ++r) dst[r] = x[r] * inv;
      ++kept;
    }
  }
  return kept;
}

// Orthonormalize the block X (b columns) against V (k columns) and the first
// m columns of Q, then within itself. Returns the kept column count.
int prepare_block(double* X, int n, int b, const double* V, int k,
                  const double* Q, int m, double floor_abs,
                  std::vector<double>* coef) {
  std::vector<double> norms(static_cast<size_t>(b), 0.0);
  for (int j = 0; j < b; ++j) {
    const double* x = X + static_cast<int64_t>(j) * n;
    double s = 0.0;
    skinny_tn(n, 1, 1, x, x, &s);
    norms[static_cast<size_t>(j)] = std::sqrt(s);
  }
  project_out_block(X, n, b, V, k, coef);
  project_out_block(X, n, b, Q, m, coef);
  return orthonormalize_block(X, n, b, norms.data(), floor_abs);
}

// Next Lanczos block from AX, reusing hcol = Q[:, 0:m]' AX as the first
// projection pass against Q. The block is then projected against V (a
// second pass only when the first removed a significant part, DGKS-style:
// AX is nearly orthogonal to an almost-invariant V, so usually one pass
// suffices), given a second full pass against Q, and orthonormalized within
// itself. Returns the kept column count.
int prepare_next_block(double* X, const double* AX, int n, int b,
                       const double* V, int k, const double* Q, int m,
                       const double* hcol, double floor_abs,
                       std::vector<double>* coef) {
  std::memcpy(X, AX, sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(b));
  std::vector<double> norms(static_cast<size_t>(b), 0.0);
  for (int j = 0; j < b; ++j) {
    const double* x = X + static_cast<int64_t>(j) * n;
    double s = 0.0;
    skinny_tn(n, 1, 1, x, x, &s);
    norms[static_cast<size_t>(j)] = std::sqrt(s);
  }
  skinny_nn_minus(n, m, b, Q, hcol, X);
  if (k > 0) {
    coef->assign(static_cast<size_t>(k) * static_cast<size_t>(b), 0.0);
    std::vector<double> before(static_cast<size_t>(b), 0.0);
    for (int j = 0; j < b; ++j) {
      const double* x = X + static_cast<int64_t>(j) * n;
      double s = 0.0;
      skinny_tn(n, 1, 1, x, x, &s);
      before[static_cast<size_t>(j)] = s;
    }
    skinny_tn(n, k, b, V, X, coef->data());
    skinny_nn_minus(n, k, b, V, coef->data(), X);
    bool again = false;
    for (int j = 0; j < b && !again; ++j) {
      const double* x = X + static_cast<int64_t>(j) * n;
      double s = 0.0;
      skinny_tn(n, 1, 1, x, x, &s);
      if (s < 0.5 * before[static_cast<size_t>(j)]) again = true;
    }
    if (again) {
      skinny_tn(n, k, b, V, X, coef->data());
      skinny_nn_minus(n, k, b, V, coef->data(), X);
    }
  }
  if (m > 0) {
    coef->assign(static_cast<size_t>(m) * static_cast<size_t>(b), 0.0);
    skinny_tn(n, m, b, Q, X, coef->data());
    skinny_nn_minus(n, m, b, Q, coef->data(), X);
  }
  return orthonormalize_block(X, n, b, norms.data(), floor_abs);
}

bool beyond(int target_kind, double theta, double edge, double margin) {
  if (target_kind == 1) return theta > edge + margin;          // largest
  if (target_kind == 2) return theta < edge - margin;          // smallest
  return std::fabs(theta) > edge + margin;                     // magnitude
}

// Symmetric eigensolve of the leading m x m block of H (leading dimension
// ldh). Returns 0 on success; values ascending, vectors column-major m x m.
int small_symmetric_eigen(const double* H, int ldh, int m,
                          std::vector<double>* values,
                          std::vector<double>* vectors) {
  vectors->assign(static_cast<size_t>(m) * static_cast<size_t>(m), 0.0);
  values->assign(static_cast<size_t>(m), 0.0);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      (*vectors)[static_cast<size_t>(i) + static_cast<size_t>(j) * m] =
        0.5 * (H[i + static_cast<int64_t>(j) * ldh] +
               H[j + static_cast<int64_t>(i) * ldh]);
    }
  }
  const char jobz = 'V';
  const char uplo = 'U';
  int info = 0;
  int lwork = -1;
  double query = 0.0;
  F77_CALL(dsyev)(&jobz, &uplo, &m, vectors->data(), &m, values->data(),
                  &query, &lwork, &info FCONE FCONE);
  if (info != 0) return info;
  lwork = static_cast<int>(query);
  if (lwork < 3 * m) lwork = 3 * m;
  std::vector<double> work(static_cast<size_t>(lwork), 0.0);
  F77_CALL(dsyev)(&jobz, &uplo, &m, vectors->data(), &m, values->data(),
                  work.data(), &lwork, &info FCONE FCONE);
  return info;
}

SEXP completeness_probe_run(void* impl, EigencoreApplyFn apply, int n,
                            SEXP V_, SEXP params_) {
  if (!isReal(V_) || !isReal(params_) || LENGTH(params_) < 7) {
    error("invalid completeness probe inputs");
  }
  SEXP dimV = getAttrib(V_, R_DimSymbol);
  if (dimV == R_NilValue || INTEGER(dimV)[0] != n) {
    error("completeness probe: V must be an n x k matrix");
  }
  const int k = INTEGER(dimV)[1];
  const double* V = REAL(V_);
  const double* params = REAL(params_);
  const int block = static_cast<int>(params[0]);
  const int steps = static_cast<int>(params[1]);
  const int target_kind = static_cast<int>(params[2]);
  const double edge = params[3];
  const double margin = params[4];
  const uint64_t stream = static_cast<uint64_t>(params[5]);
  const double scale = (std::isfinite(params[6]) && params[6] > 0.0) ? params[6] : 1.0;
  if (block < 1 || steps < 1) error("completeness probe: block and steps must be >= 1");
  if (target_kind < 1 || target_kind > 3) error("completeness probe: unsupported target");

  const int nc = n - k;
  int cap = block * steps;
  if (cap > nc) cap = nc < 0 ? 0 : nc;

  std::vector<double> Q(static_cast<size_t>(n) * static_cast<size_t>(cap > 0 ? cap : 1), 0.0);
  std::vector<double> H(static_cast<size_t>(cap > 0 ? cap : 1) * static_cast<size_t>(cap > 0 ? cap : 1), 0.0);
  std::vector<double> X(static_cast<size_t>(n) * static_cast<size_t>(block), 0.0);
  std::vector<double> AX(static_cast<size_t>(n) * static_cast<size_t>(block), 0.0);
  std::vector<double> coef;
  std::vector<double> hcol;
  std::vector<double> theta;
  std::vector<double> U;

  ProbeRng rng{0x5EEDC0DEULL ^ (stream * 0x9E3779B97F4A7C15ULL)};
  auto seed_block = [&](double* dst, int cols) {
    for (int64_t i = 0; i < static_cast<int64_t>(n) * cols; ++i) dst[i] = rng.normal();
  };

  int m = 0;
  int used_steps = 0;
  int block_calls = 0;
  int columns = 0;
  int reseeds = 0;
  bool exhausted = (nc <= 0);
  bool intruder = false;
  int width = 0;
  const double floor_abs = 1e-12 * scale;

  if (!exhausted) {
    seed_block(X.data(), block);
    width = prepare_block(X.data(), n, block, V, k, Q.data(), 0, 0.0, &coef);
  }
  for (int step = 0; step < steps && !exhausted; ++step) {
    eigencore_check_interrupt();
    while (width == 0) {
      // Breakdown: the explored space is invariant; re-seed in the complement
      // of everything explored so far. A re-seed that finds nothing new (or a
      // full complement) means the complement is exhausted.
      if (m >= nc) {
        exhausted = true;
        break;
      }
      if (reseeds >= 3) break;
      ++reseeds;
      seed_block(X.data(), block);
      width = prepare_block(X.data(), n, block, V, k, Q.data(), m, 0.0, &coef);
      if (width == 0 && reseeds >= 3) {
        exhausted = true;
      }
    }
    if (exhausted || width == 0) break;
    if (width > cap - m) width = cap - m;
    if (width <= 0) break;

    std::fill(AX.begin(), AX.end(), 0.0);
    const int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, width, X.data(), n,
                         1.0, 0.0, AX.data(), n, nullptr);
    if (rc != 0) {
      error("completeness probe: operator apply failed with status=%d", rc);
    }
    ++block_calls;
    columns += width;
    used_steps = step + 1;
    std::memcpy(Q.data() + static_cast<int64_t>(m) * n, X.data(),
                sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(width));
    const int m_new = m + width;
    // New block column of H: Q[:, 0:m_new]' * AX.
    hcol.assign(static_cast<size_t>(m_new) * static_cast<size_t>(width), 0.0);
    skinny_tn(n, m_new, width, Q.data(), AX.data(), hcol.data());
    {
      for (int j = 0; j < width; ++j) {
        for (int i = 0; i < m_new; ++i) {
          const double h = hcol[static_cast<size_t>(i) + static_cast<size_t>(j) * m_new];
          H[static_cast<size_t>(i) + static_cast<size_t>(m + j) * cap] = h;
          H[static_cast<size_t>(m + j) + static_cast<size_t>(i) * cap] = h;
        }
      }
    }
    m = m_new;
    if (small_symmetric_eigen(H.data(), cap, m, &theta, &U) != 0) {
      error("completeness probe: projected eigensolve failed");
    }
    for (int i = 0; i < m; ++i) {
      if (beyond(target_kind, theta[static_cast<size_t>(i)], edge, margin)) {
        intruder = true;
        break;
      }
    }
    if (intruder) break;
    if (m >= nc) {
      exhausted = true;
      break;
    }
    if (m >= cap) break;
    width = prepare_next_block(X.data(), AX.data(), n, width, V, k, Q.data(),
                               m, hcol.data(), floor_abs, &coef);
  }

  // Ritz vectors Y = Q[:, 0:m] U.
  SEXP theta_ = PROTECT(allocVector(REALSXP, m));
  SEXP Y_ = PROTECT(allocMatrix(REALSXP, n, m));
  for (int i = 0; i < m; ++i) REAL(theta_)[i] = theta[static_cast<size_t>(i)];
  if (m > 0) {
    const char trans_N = 'N';
    const double one = 1.0;
    const double zero = 0.0;
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &m, &m, &one, Q.data(), &n,
                    U.data(), &m, &zero, REAL(Y_), &n FCONE FCONE);
  }
  const char* names[] = {"theta", "Y", "steps", "block_calls", "columns",
                         "exhausted", "intruder", "reseeds", ""};
  SEXP out = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, theta_);
  SET_VECTOR_ELT(out, 1, Y_);
  SET_VECTOR_ELT(out, 2, ScalarInteger(used_steps));
  SET_VECTOR_ELT(out, 3, ScalarInteger(block_calls));
  SET_VECTOR_ELT(out, 4, ScalarInteger(columns));
  SET_VECTOR_ELT(out, 5, ScalarLogical(exhausted ? 1 : 0));
  SET_VECTOR_ELT(out, 6, ScalarLogical(intruder ? 1 : 0));
  SET_VECTOR_ELT(out, 7, ScalarInteger(reseeds));
  UNPROTECT(3);
  return out;
}

}  // namespace

extern "C" SEXP eigencore_completeness_probe_dense(SEXP A_, SEXP V_, SEXP params_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) error("A must be a double matrix");
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue || INTEGER(dimA)[0] != INTEGER(dimA)[1]) {
    error("A must be a square matrix");
  }
  const int n = INTEGER(dimA)[0];
  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  return completeness_probe_run(&impl, eigencore_dense_apply, n, V_, params_);
  EIGENCORE_ENTRY_END
}

namespace {

// Serial CSC forward apply for the probe. The shared CSC kernel builds a
// per-operator CSR/slab cache on its second multithreaded apply, which a
// fresh 8-step probe operator would rebuild on every call (tens of ms at
// n = 20000); a serial column scatter costs far less for a handful of
// applies.
struct ProbeCsc {
  int n;
  const int* row_idx;
  const int* col_ptr;
  const double* values;
};

int probe_csc_apply(void* impl, EigencoreTranspose op, int64_t block_cols,
                    const double* X, int64_t ldx, double alpha, double beta,
                    double* Y, int64_t ldy, EigencoreWorkspace* workspace) {
  (void) workspace;
  if (op != EIGENCORE_TRANSPOSE_NONE) return -1;
  const ProbeCsc* A = static_cast<const ProbeCsc*>(impl);
  const int n = A->n;
  for (int64_t c = 0; c < block_cols; ++c) {
    const double* x = X + c * ldx;
    double* y = Y + c * ldy;
    if (beta == 0.0) {
      std::memset(y, 0, sizeof(double) * static_cast<size_t>(n));
    } else if (beta != 1.0) {
      for (int r = 0; r < n; ++r) y[r] *= beta;
    }
    for (int j = 0; j < n; ++j) {
      const double xj = alpha * x[j];
      if (xj == 0.0) continue;
      for (int p = A->col_ptr[j]; p < A->col_ptr[j + 1]; ++p) {
        y[A->row_idx[p]] += A->values[p] * xj;
      }
    }
  }
  return 0;
}

}  // namespace

extern "C" SEXP eigencore_completeness_probe_csc(SEXP i_, SEXP p_, SEXP x_,
                                                 SEXP dim_, SEXP V_,
                                                 SEXP params_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_)) {
    error("invalid CSC completeness probe inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "completeness probe");
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n) error("A must be square");
  ProbeCsc impl = {n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return completeness_probe_run(&impl, probe_csc_apply, n, V_, params_);
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_completeness_probe_r_operator(SEXP dim_, SEXP apply_,
                                                        SEXP V_, SEXP params_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(dim_) || LENGTH(dim_) != 2 || TYPEOF(apply_) != CLOSXP) {
    error("invalid matrix-free completeness probe inputs");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n) error("A must be square");
  RApplyOperator impl = {n, n, apply_, R_NilValue};
  return completeness_probe_run(&impl, eigencore_r_operator_apply, n, V_, params_);
  EIGENCORE_ENTRY_END
}
