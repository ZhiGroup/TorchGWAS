/* CPU BGEN Layout 2 decode: inflate, validate and unpack, in C.
 *
 * The Python path in bgen.py allocates several cohort-width temporaries for
 * every variant -- a padded copy of the probability block, a uint64 word per
 * sample, then shifts and masks over all of them -- and holds the GIL between
 * those numpy calls, so reader threads serialise. It is the slowest decode we
 * have by a wide margin.
 *
 * This does the same arithmetic in one pass with no allocation, writing
 * straight into the caller's chunk column. ctypes drops the GIL around the
 * call, so the reader threads that already exist actually overlap.
 *
 * Every validation below mirrors one in `_decode_record`, and the test
 * compares refusals as well as values: a record that is truncated,
 * non-diploid, phased, or carries probabilities summing past one must be
 * refused rather than decoded into plausible-looking numbers.
 */
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

#define TG_BGEN_CPU_OK 0
#define TG_BGEN_CPU_TRUNCATED (-1)
#define TG_BGEN_CPU_INFLATE (-2)
#define TG_BGEN_CPU_RAW_LENGTH (-3)
#define TG_BGEN_CPU_DIMENSIONS (-4)
#define TG_BGEN_CPU_PLOIDY (-5)
#define TG_BGEN_CPU_PHASED (-6)
#define TG_BGEN_CPU_BITS (-7)
#define TG_BGEN_CPU_PROBABILITY (-8)
#define TG_BGEN_CPU_COMPRESSION (-9)
#define TG_BGEN_CPU_SCRATCH (-10)
#define TG_BGEN_CPU_SAMPLE_INDEX (-11)

int tg_bgen_cpu_abi_version(void) { return 2; }

/* Transpose a (rows x columns) float matrix into (columns x rows).
 *
 * A chunk is decoded variant-major because writing each variant into its
 * column of a (samples, variants) buffer scatters one float every `variants`
 * floats and costs about four times the decode itself. That leaves one
 * transpose per chunk, and doing it in numpy put the cost back where it
 * started: `ascontiguousarray` holds the GIL, so every reader thread queued
 * behind it and concurrency stopped at about three cores of forty-eight.
 *
 * Here it is blocked so both sides stay in cache, and ctypes releases the GIL
 * for the call, so it overlaps across readers like the decode does.
 */
void tg_bgen_cpu_transpose(const float *source, int64_t rows, int64_t columns,
                           float *destination) {
  const int64_t block = 64;  /* 64 x 64 floats is 16 KB, comfortably in L1 */
  for (int64_t i0 = 0; i0 < rows; i0 += block) {
    const int64_t i1 = i0 + block < rows ? i0 + block : rows;
    for (int64_t j0 = 0; j0 < columns; j0 += block) {
      const int64_t j1 = j0 + block < columns ? j0 + block : columns;
      for (int64_t i = i0; i < i1; ++i) {
        const float *row = source + i * columns;
        for (int64_t j = j0; j < j1; ++j) destination[j * rows + i] = row[j];
      }
    }
  }
}

const char *tg_bgen_cpu_message(int code) {
  switch (code) {
    case TG_BGEN_CPU_OK: return "ok";
    case TG_BGEN_CPU_TRUNCATED: return "truncated BGEN probability block";
    case TG_BGEN_CPU_INFLATE: return "zlib inflate failed";
    case TG_BGEN_CPU_RAW_LENGTH: return "BGEN inflated length mismatch";
    case TG_BGEN_CPU_DIMENSIONS: return "unsupported BGEN dimensions";
    case TG_BGEN_CPU_PLOIDY: return "non-diploid selected sample";
    case TG_BGEN_CPU_PHASED: return "phased BGEN record";
    case TG_BGEN_CPU_BITS: return "invalid BGEN packed probability length";
    case TG_BGEN_CPU_PROBABILITY: return "invalid BGEN probabilities";
    case TG_BGEN_CPU_COMPRESSION: return "unsupported BGEN compression";
    case TG_BGEN_CPU_SCRATCH: return "inflate scratch buffer too small";
    case TG_BGEN_CPU_SAMPLE_INDEX: return "sample index outside cohort";
    default: return "unknown BGEN decode error";
  }
}

#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define TG_LITTLE_ENDIAN 1
#else
#define TG_LITTLE_ENDIAN 0
#endif

/* Assemble one value byte by byte. Correct anywhere and safe at the end of a
 * block, but about forty operations per sample -- only the tail uses it. */
static inline uint64_t read_bits_slow(const uint8_t *data, size_t length,
                                      uint64_t offset, uint32_t bits) {
  const uint64_t first = offset >> 3;
  const uint32_t shift = (uint32_t)(offset & 7u);
  const uint32_t needed = (bits + shift + 7u) >> 3;
  uint64_t word = 0;
  for (uint32_t i = 0; i < needed; ++i) {
    const uint64_t index = first + i;
    const uint64_t byte = index < length ? (uint64_t)data[index] : 0u;
    word |= byte << (8u * i);
  }
  const uint64_t mask = bits >= 64 ? ~(uint64_t)0 : (((uint64_t)1 << bits) - 1u);
  return (word >> shift) & mask;
}

/* Read `bits` bits starting at bit `offset`, little-endian within the stream.
 *
 * The byte-by-byte version above was the whole cost of this decoder: about
 * forty operations per sample, called twice per sample, which came to 1.5 ms
 * per variant at 35,365 samples -- around 43 ns to extract two integers. One
 * unaligned 64-bit load covers any depth up to 32 whatever the bit offset
 * (a 7-bit shift plus 32 bits needs 39), so the loop only exists for the last
 * seven bytes of a block, where a wide load would read past the end.
 *
 * The load is a memcpy so it stays defined for an unaligned address; every
 * compiler we build with turns it into one instruction. */
static inline uint64_t read_bits(const uint8_t *data, size_t length,
                                 uint64_t offset, uint32_t bits) {
  const uint64_t first = offset >> 3;
#if TG_LITTLE_ENDIAN
  if (bits <= 32 && first + 8 <= length) {
    uint64_t word;
    memcpy(&word, data + first, sizeof word);
    const uint64_t mask = ((uint64_t)1 << bits) - 1u;
    return (word >> (offset & 7u)) & mask;
  }
#endif
  return read_bits_slow(data, length, offset, bits);
}

/* `payload` is one record's probability block, compressed as the file header
 * declares; `raw_length` is the inflated size the record states. `scratch`
 * receives the inflated block. Results are written to
 * out[j * out_stride] for each selected sample j: the alternate-allele
 * dosage, or NaN where the record marks that sample missing. */
int tg_bgen_cpu_decode(const uint8_t *payload, size_t payload_length,
                       int compression, uint32_t raw_length,
                       uint32_t n_bgen_samples, const uint32_t *sample_indices,
                       uint32_t n_selected, uint8_t *scratch,
                       size_t scratch_length, float *out, int64_t out_stride) {
  const uint8_t *raw;
  size_t raw_size;

  if (compression == 0) {
    raw = payload;
    raw_size = payload_length;
  } else if (compression == 1) {
    if (scratch_length < (size_t)raw_length) return TG_BGEN_CPU_SCRATCH;
    uLongf produced = (uLongf)raw_length;
    if (uncompress((Bytef *)scratch, &produced, (const Bytef *)payload,
                   (uLong)payload_length) != Z_OK)
      return TG_BGEN_CPU_INFLATE;
    raw = scratch;
    raw_size = (size_t)produced;
  } else {
    /* zstd records stay on the Python path, as they do on the GPU path. */
    return TG_BGEN_CPU_COMPRESSION;
  }

  if (raw_size != (size_t)raw_length) return TG_BGEN_CPU_RAW_LENGTH;
  if (raw_size < (size_t)10 + n_bgen_samples) return TG_BGEN_CPU_RAW_LENGTH;

  uint32_t declared_samples;
  uint16_t alleles;
  memcpy(&declared_samples, raw, sizeof declared_samples);
  memcpy(&alleles, raw + 4, sizeof alleles);
  if (declared_samples != n_bgen_samples || alleles != 2 || raw[6] != 2 ||
      raw[7] != 2)
    return TG_BGEN_CPU_DIMENSIONS;
  if (raw[8 + n_bgen_samples] != 0) return TG_BGEN_CPU_PHASED;

  const uint32_t bits = raw[9 + n_bgen_samples];
  if (bits < 1 || bits > 32) return TG_BGEN_CPU_BITS;
  const uint64_t packed_bits = (uint64_t)n_bgen_samples * 2u * bits;
  const size_t expected =
      (size_t)10 + n_bgen_samples + (size_t)((packed_bits + 7u) >> 3);
  if (raw_size != expected) return TG_BGEN_CPU_BITS;

  const uint8_t *ploidy = raw + 8;
  const uint8_t *data = raw + 10 + n_bgen_samples;
  const size_t data_length = raw_size - (size_t)10 - n_bgen_samples;
  const uint64_t denominator = ((uint64_t)1 << bits) - 1u;
  const double inverse = 1.0 / (double)denominator;

  /* Ploidy and phase are properties of the record, so reject before writing
   * anything: a partially filled column would otherwise outlive the error. */
  for (uint32_t j = 0; j < n_selected; ++j) {
    const uint32_t sample = sample_indices[j];
    if (sample >= n_bgen_samples) return TG_BGEN_CPU_SAMPLE_INDEX;
    if ((ploidy[sample] & 63u) != 2u) return TG_BGEN_CPU_PLOIDY;
  }

  /* At a byte-aligned depth both probabilities sit in one 64-bit window, so
   * the pair costs a single load instead of two. Real files are overwhelmingly
   * 8- or 16-bit -- UK Biobank is 8 -- and `bits` is fixed for the record, so
   * this branch predicts perfectly.
   */
  const int paired = TG_LITTLE_ENDIAN && (bits % 8 == 0) && bits <= 32;
  const uint32_t pair_bytes = bits / 8 * 2;

  for (uint32_t j = 0; j < n_selected; ++j) {
    const uint32_t sample = sample_indices[j];
    if (ploidy[sample] & 128u) {
      out[(int64_t)j * out_stride] = NAN;
      continue;
    }
    uint64_t p0, p1;
    const uint64_t byte0 = (uint64_t)sample * pair_bytes;
    if (paired && byte0 + 8 <= data_length) {
      uint64_t word;
      memcpy(&word, data + byte0, sizeof word);
      p0 = word & denominator;
      p1 = (word >> bits) & denominator;
    } else {
      const uint64_t offset = (uint64_t)sample * 2u * bits;
      p0 = read_bits(data, data_length, offset, bits);
      p1 = read_bits(data, data_length, offset + bits, bits);
    }
    if (p0 + p1 > denominator) return TG_BGEN_CPU_PROBABILITY;
    /* Alternate dosage p1 + 2*p2 with p2 = 1 - p0 - p1, in the record's own
     * fixed-point units, divided once at the end as the Python path does. */
    out[(int64_t)j * out_stride] =
        (float)((double)(2u * denominator - 2u * p0 - p1) * inverse);
  }
  return TG_BGEN_CPU_OK;
}
