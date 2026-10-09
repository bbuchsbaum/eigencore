// Native structural hash for operator identity and workflow tokens (C45).
//
// eigencore_identity_hash(x) walks an R object and hashes its type, length,
// attributes and values directly from the data buffers: no serialisation and
// no copy of the payload. It replaces stable_raw_hash(serialize(x)), whose
// serialisation of a dense source matrix dominated plan time.
//
// Design ("eigencore identity hash v2", 128-bit output, 32 hex characters):
//
// * The stream is a sequence of 64-bit words fed into four XXH64-style lanes
//   (32 bytes per round). Words are formed from *values*, never from raw
//   memory order: doubles via their IEEE-754 bit pattern as a uint64,
//   integers/logicals packed two per word as uint32 values, and bytes
//   (strings, raw vectors) packed in an explicit little-endian order with
//   shifts. The digest therefore does not depend on the host endianness.
// * Doubles are canonicalised so values R treats as identical() hash equal:
//   -0.0 hashes as 0.0, and every non-NA NaN payload hashes as one canonical
//   NaN. NA_real_ stays distinct from NaN (identical() distinguishes them).
// * Every node contributes a type tag and its length; attributes are hashed
//   sorted by name (identical() ignores attribute order), so dim, dimnames,
//   class and S4 slots (Matrix i, p, x, Dim, uplo, diag, ...) all count.
//   Dimnames are part of the identity, as they were under the serialised
//   digest.
// * Strings are hashed as UTF-8 bytes (Rf_translateCharUTF8) with their
//   length; NA_character_ has its own marker.
// * Nodes without a value-level fast path (closures, environments, external
//   pointers, bytecode, builtins, ...) are hashed through the old route:
//   base::serialize(node, NULL, version = 3) and the raw bytes.
// * Two different finalisations of the 256-bit lane state give the two
//   64-bit halves of the digest. This is a non-cryptographic hash: it guards
//   against accidental mismatch, not adversarial collision.

#include <algorithm>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <utility>
#include <vector>
#include "eigencore_common.h"


namespace {

constexpr std::uint64_t P1 = UINT64_C(0x9E3779B185EBCA87);
constexpr std::uint64_t P2 = UINT64_C(0xC2B2AE3D27D4EB4F);
constexpr std::uint64_t P3 = UINT64_C(0x165667B19E3779F9);
constexpr std::uint64_t P4 = UINT64_C(0x85EBCA77C2B2AE63);
constexpr std::uint64_t P5 = UINT64_C(0x27D4EB2F165667C5);

constexpr std::uint64_t kCanonicalNaN = UINT64_C(0x7FF8000000000000);
constexpr std::uint64_t kNAStringMarker = UINT64_C(0xFFFFFFFFFFFFFFFF);
constexpr std::uint64_t kNodeTag = UINT64_C(0xEC1D000000000000);
constexpr std::uint64_t kFormatVersion = 2;

inline std::uint64_t rotl64(std::uint64_t x, int r) {
  return (x << r) | (x >> (64 - r));
}

inline std::uint64_t lane_round(std::uint64_t acc, std::uint64_t input) {
  acc += input * P2;
  acc = rotl64(acc, 31);
  return acc * P1;
}

inline std::uint64_t merge_round(std::uint64_t acc, std::uint64_t val) {
  val = lane_round(0, val);
  acc ^= val;
  return acc * P1 + P4;
}

inline std::uint64_t avalanche(std::uint64_t h) {
  h ^= h >> 33;
  h *= P2;
  h ^= h >> 29;
  h *= P3;
  h ^= h >> 32;
  return h;
}

inline std::uint64_t canonical_double_bits(double value) {
  std::uint64_t bits;
  std::memcpy(&bits, &value, sizeof(bits));
  if ((bits << 1) == 0) {
    return 0;  // +0.0 and -0.0
  }
  if ((bits & UINT64_C(0x7FF0000000000000)) == UINT64_C(0x7FF0000000000000) &&
      (bits & UINT64_C(0x000FFFFFFFFFFFFF)) != 0) {
    // NaN: keep R's NA_real_ (low word 1954) distinct, canonicalise the rest.
    if (R_IsNA(value)) {
      return UINT64_C(0x7FF00000000007A2);
    }
    return kCanonicalNaN;
  }
  return bits;
}

class StreamHash {
 public:
  StreamHash() {
    v_[0] = P1 + P2;
    v_[1] = P2;
    v_[2] = 0;
    v_[3] = static_cast<std::uint64_t>(0) - P1;
  }

  inline void word(std::uint64_t w) {
    buf_[n_++] = w;
    if (n_ == 4) {
      v_[0] = lane_round(v_[0], buf_[0]);
      v_[1] = lane_round(v_[1], buf_[1]);
      v_[2] = lane_round(v_[2], buf_[2]);
      v_[3] = lane_round(v_[3], buf_[3]);
      n_ = 0;
    }
    ++total_;
  }

  void doubles(const double* x, R_xlen_t n) {
    for (R_xlen_t i = 0; i < n; ++i) {
      word(canonical_double_bits(x[i]));
    }
  }

  void int32s(const int* x, R_xlen_t n) {
    R_xlen_t i = 0;
    for (; i + 1 < n; i += 2) {
      word(static_cast<std::uint64_t>(static_cast<std::uint32_t>(x[i])) |
           (static_cast<std::uint64_t>(static_cast<std::uint32_t>(x[i + 1])) << 32));
    }
    if (i < n) {
      word(static_cast<std::uint64_t>(static_cast<std::uint32_t>(x[i])));
    }
  }

  void bytes(const unsigned char* x, std::size_t n) {
    word(static_cast<std::uint64_t>(n));
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
      std::uint64_t w = 0;
      for (int b = 0; b < 8; ++b) {
        w |= static_cast<std::uint64_t>(x[i + b]) << (8 * b);
      }
      word(w);
    }
    if (i < n) {
      std::uint64_t w = 0;
      for (int b = 0; i + b < n; ++b) {
        w |= static_cast<std::uint64_t>(x[i + b]) << (8 * b);
      }
      word(w);
    }
  }

  void finish(std::uint64_t out[2]) const {
    out[0] = finish_one(0);
    out[1] = finish_one(1);
  }

 private:
  std::uint64_t finish_one(int variant) const {
    std::uint64_t h;
    if (variant == 0) {
      h = rotl64(v_[0], 1) + rotl64(v_[1], 7) + rotl64(v_[2], 12) + rotl64(v_[3], 18);
      for (int i = 0; i < 4; ++i) h = merge_round(h, v_[i]);
    } else {
      h = rotl64(v_[3], 3) + rotl64(v_[2], 11) + rotl64(v_[1], 23) + rotl64(v_[0], 41) + P5;
      for (int i = 3; i >= 0; --i) h = merge_round(h, v_[i] ^ P3);
    }
    h += total_ * 8;
    for (int i = 0; i < n_; ++i) {
      std::uint64_t k = lane_round(variant == 0 ? 0 : P5, buf_[i]);
      h ^= k;
      h = rotl64(h, 27) * P1 + P4;
    }
    if (variant == 1) {
      h ^= rotl64(h, 17) * P3;
    }
    return avalanche(h);
  }

  std::uint64_t v_[4];
  std::uint64_t buf_[4] = {0, 0, 0, 0};
  int n_ = 0;
  std::uint64_t total_ = 0;
};

inline void hash_charsxp(StreamHash& h, SEXP s) {
  if (s == NA_STRING) {
    h.word(kNAStringMarker);
    return;
  }
  const char* str = Rf_translateCharUTF8(s);
  h.bytes(reinterpret_cast<const unsigned char*>(str), std::strlen(str));
}

void hash_node(StreamHash& h, SEXP x);

void hash_serialized(StreamHash& h, SEXP x) {
  SEXP version = PROTECT(Rf_ScalarInteger(3));
  SEXP call = PROTECT(Rf_lang4(Rf_install("serialize"), x, R_NilValue, version));
  SET_TAG(CDR(CDR(CDR(call))), Rf_install("version"));
  SEXP raw = PROTECT(eigencore_unwind_protect([&] { return Rf_eval(call, R_BaseEnv); }));
  h.bytes(RAW(raw), static_cast<std::size_t>(XLENGTH(raw)));
  UNPROTECT(3);
}

// Attributes are read through base::attributes() rather than ATTRIB(), which
// is not part of R's API (R >= 4.5 no longer exports it). The argument is
// quoted so language objects and symbols are not evaluated.
void hash_attributes(StreamHash& h, SEXP x) {
  SEXP quoted = PROTECT(Rf_lang2(Rf_install("quote"), x));
  SEXP call = PROTECT(Rf_lang2(Rf_install("attributes"), quoted));
  SEXP attrs = PROTECT(eigencore_unwind_protect([&] { return Rf_eval(call, R_BaseEnv); }));
  if (attrs == R_NilValue || XLENGTH(attrs) == 0) {
    UNPROTECT(3);
    h.word(0);
    return;
  }
  SEXP names = Rf_getAttrib(attrs, R_NamesSymbol);
  std::vector<std::pair<const char*, SEXP>> entries;
  for (R_xlen_t i = 0; i < XLENGTH(attrs); ++i) {
    entries.emplace_back(CHAR(STRING_ELT(names, i)), VECTOR_ELT(attrs, i));
  }
  std::sort(entries.begin(), entries.end(),
            [](const std::pair<const char*, SEXP>& a,
               const std::pair<const char*, SEXP>& b) {
              return std::strcmp(a.first, b.first) < 0;
            });
  h.word(static_cast<std::uint64_t>(entries.size()));
  for (const auto& entry : entries) {
    h.bytes(reinterpret_cast<const unsigned char*>(entry.first),
            std::strlen(entry.first));
    hash_node(h, entry.second);
  }
  UNPROTECT(3);
}

void hash_node(StreamHash& h, SEXP x) {
  const int type = TYPEOF(x);
  h.word(kNodeTag | (static_cast<std::uint64_t>(Rf_isS4(x) ? 1 : 0) << 8) |
         static_cast<std::uint64_t>(type));
  switch (type) {
    case NILSXP:
      return;
    case LGLSXP:
      h.word(static_cast<std::uint64_t>(XLENGTH(x)));
      h.int32s(LOGICAL_RO(x), XLENGTH(x));
      break;
    case INTSXP:
      h.word(static_cast<std::uint64_t>(XLENGTH(x)));
      h.int32s(INTEGER_RO(x), XLENGTH(x));
      break;
    case REALSXP:
      h.word(static_cast<std::uint64_t>(XLENGTH(x)));
      h.doubles(REAL_RO(x), XLENGTH(x));
      break;
    case CPLXSXP: {
      const R_xlen_t n = XLENGTH(x);
      h.word(static_cast<std::uint64_t>(n));
      const Rcomplex* z = COMPLEX_RO(x);
      for (R_xlen_t i = 0; i < n; ++i) {
        h.word(canonical_double_bits(z[i].r));
        h.word(canonical_double_bits(z[i].i));
      }
      break;
    }
    case RAWSXP:
      h.bytes(RAW_RO(x), static_cast<std::size_t>(XLENGTH(x)));
      break;
    case STRSXP: {
      const R_xlen_t n = XLENGTH(x);
      h.word(static_cast<std::uint64_t>(n));
      for (R_xlen_t i = 0; i < n; ++i) hash_charsxp(h, STRING_ELT(x, i));
      break;
    }
    case VECSXP:
    case EXPRSXP: {
      const R_xlen_t n = XLENGTH(x);
      h.word(static_cast<std::uint64_t>(n));
      for (R_xlen_t i = 0; i < n; ++i) hash_node(h, VECTOR_ELT(x, i));
      break;
    }
    case SYMSXP:
      hash_charsxp(h, PRINTNAME(x));
      return;  // symbols carry no attributes
    case LISTSXP:
    case LANGSXP: {
      std::uint64_t n = 0;
      for (SEXP node = x; node != R_NilValue; node = CDR(node)) ++n;
      h.word(n);
      for (SEXP node = x; node != R_NilValue; node = CDR(node)) {
        if (TAG(node) == R_NilValue) {
          h.word(0);
        } else {
          h.word(1);
          hash_charsxp(h, PRINTNAME(TAG(node)));
        }
        hash_node(h, CAR(node));
      }
      break;
    }
    case S4SXP:
      break;  // an S4 object's slots are its attributes
    default:
      // Closures, environments, external pointers, bytecode, promises,
      // builtins: no value-level fast path. Serialise this node only.
      hash_serialized(h, x);
      return;
  }
  hash_attributes(h, x);
}

}  // namespace

extern "C" SEXP eigencore_identity_hash(SEXP x) {
  EIGENCORE_ENTRY_BEGIN
  StreamHash h;
  h.word(kNodeTag | (kFormatVersion << 16));
  hash_node(h, x);
  std::uint64_t out[2];
  h.finish(out);
  char output[33];
  std::snprintf(output, sizeof(output), "%016" PRIx64 "%016" PRIx64, out[0], out[1]);
  return Rf_mkString(output);
  EIGENCORE_ENTRY_END
}
