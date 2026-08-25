#include <R.h>
#include <Rinternals.h>

#include <cinttypes>
#include <cstdint>
#include <cstdio>

extern "C" SEXP eigencore_stable_raw_hash(SEXP raw_) {
  if (TYPEOF(raw_) != RAWSXP) {
    Rf_error("eigencore_stable_raw_hash requires a raw vector");
  }

  const Rbyte* bytes = RAW(raw_);
  const R_xlen_t size = XLENGTH(raw_);
  std::uint64_t hash = UINT64_C(14695981039346656037);
  constexpr std::uint64_t prime = UINT64_C(1099511628211);
  for (R_xlen_t i = 0; i < size; ++i) {
    hash ^= static_cast<std::uint64_t>(bytes[i]);
    hash *= prime;
  }

  char output[17];
  std::snprintf(output, sizeof(output), "%016" PRIx64, hash);
  return Rf_mkString(output);
}
