// Fused genotype-by-design product: decode, centre and multiply in one kernel.
//
// The unfused path writes the centred genotype block to global memory and then
// reads it back into a cuBLAS GEMM. At a 10,000-variant chunk and 22,250
// samples that is 890 MB written and 890 MB read per chunk -- sixteen times
// the packed genotypes themselves, because two bits per call become four
// bytes, and the largest single item in the scan's device traffic.
//
// This kernel never materialises it. A block owns a tile of variants and a
// tile of design columns, walks the samples, decodes two-bit calls into shared
// memory once per sample tile, centres them against a mean from an earlier
// cheap pass, and accumulates against a tile of the design staged beside them.
//
// Three parameters decide whether that is fast, and they pull in different
// directions, so the kernel is templated on all three and the caller sweeps
// them rather than the author guessing:
//
//   kSampleTile  how many samples are staged per pass. This is a *load*
//                parameter before it is a blocking parameter. Packed calls are
//                two bits, so one variant contributes kSampleTile/4 bytes per
//                pass: at 16 samples that is four bytes out of a 32-byte
//                sector. The first version of this kernel fixed it at 16 and
//                paid eight times the genotype read traffic for it.
//
//   kMicroVariants x kMicroColumns  the rectangle of the output each thread
//                keeps in registers. This sets the ratio that binds the inner
//                loop: a thread reads MV + MC values from shared and does
//                MV * MC multiply-adds, so 4x4 is two FMAs per read and 8x8 is
//                four. Bigger costs registers and so occupancy.
//
// Shared memory is laid out sample-major, because the inner loop walks
// variants and columns for one fixed sample and those must be adjacent. The
// genotype tile is padded by one float per staged sample: the staging loop
// writes with consecutive threads varying the sample, which unpadded is a
// stride of kVariantTile floats and lands every one of them in the same bank.
//
// The arithmetic matches the unfused path: same centring, same masked missing
// calls. It is not bit-identical, because the reduction is ordered by sample
// tile rather than however cuBLAS ordered it -- the same reason two GEMM
// tilings of one product differ in the last digits. Note that this file is
// built *with* FMA, unlike the statistics kernel: cuBLAS accumulates with
// FFMA, so forbidding it here would not make the comparison fairer, it would
// make this kernel issue two instructions where cuBLAS issues one.
#include <cuda_runtime.h>
#include <cuda_pipeline.h>
#include <cstdint>
#include <cmath>
#include <string>

namespace {
thread_local std::string last_error;

// A 16 x 16 thread grid over 256 threads. Every configuration below keeps this
// shape and varies only what each thread owns, so the staging loops are
// written once.
constexpr int kThreads = 256;
constexpr int kThreadRows = 16;
constexpr int kThreadColumns = 16;

enum PackedCodes { kPackedPgen = 1, kPackedPlink1 = 2 };

// THE TWO BRANCHES COUNT DIFFERENT ALLELES, AND EACH IS CORRECT FOR ITS FORMAT.
// PLINK1 decodes to A2 dosage (.bim column 6); PGEN decodes to ALT. Verified
// against the official PLINK example (cog-genomics .bed spec): bim "G A", PED
// GG/AA/missing/AA/AA/AA decodes to [0, 2, nan, 2, 2, 2] and the reader reports
// effect allele "A" -- self-consistent, counting and declaring the same allele.
//
// For a .bim written from a .pvar, A2 == REF and ALT == A1, so the SAME cohort
// read as .bed and as .pgen yields dosages summing to exactly 2.0 per variant
// and betas of opposite sign. That is a CONTRACT DIFFERENCE, not a decode bug:
// each reader names its own effect allele. Making them agree means flipping the
// dosage AND the declared effect allele together -- flipping only the dosage
// (attempted 2026-09-15) leaves the reader counting A1 while announcing A2,
// which is worse than either consistent choice.
__device__ __forceinline__ float decode_two_bit(unsigned code, int kind) {
  if (kind == kPackedPlink1) return code == 1u ? nanf("") : float((code + 1u) >> 1);
  return code == 3u ? nanf("") : float(code);
}

// products[variant][column] = sum_sample centred(variant, sample) * design[sample][column]
//
// A missing call is centred to exactly zero, so it drops out of every product
// without needing a mask tensor alongside.
// `kMinBlocks` is the second argument to __launch_bounds__: the occupancy the
// kernel promises to reach, which ptxas honours by capping registers. At 256
// threads the budget is 65536 / (256 * kMinBlocks) registers -- 128 at two
// blocks, 85 at three, 64 at four. The measured configurations sit at 65-118
// registers with *zero* bytes spilled, so there is headroom to spend: trading a
// few spills for another resident block is worth measuring on a loop that is
// latency-bound rather than issue-bound. kMinBlocks = 1 is effectively
// unconstrained, since its implied cap of 256 exceeds the hardware maximum of
// 255 registers per thread.
// `kAsyncDesign` stages the design tile with `cp.async` instead of routing it
// through registers. It only applies with double buffering, because the point
// is to have the copy in flight while the previous tile computes.
//
// The register-prefetch version of this lost, and the reason was occupancy: it
// doubled shared memory *and* held the values in registers. `cp.async` writes
// global-to-shared without the register round trip, and the occupancy
// arithmetic says the doubled shared memory is free here -- at 96 registers a
// block is limited to 2-3 per multiprocessor by registers, while two staging
// buffers come to about 51 KB against the 228 KB an H100 offers, so shared
// memory is not what binds.
//
// Copies are four bytes rather than sixteen: `columns` is a trait count plus
// covariates and is not generally a multiple of four, so a row of the design
// has no alignment guarantee and a wider `cp.async` would be unaligned.
// The body lives in a __device__ function so two entry points can wrap it with
// different launch bounds. `__launch_bounds__(kThreads)` and
// `__launch_bounds__(kThreads, 1)` are **not** the same thing: the second
// promises only one block per multiprocessor, which licences ptxas to spend up
// to 255 registers, and configurations that had peaked at 118 rebuilt at 168
// and ran slower. Without a second wrapper there is no way to express "no
// hint", and the occupancy sweep was scored against a regressed baseline.
template <int kCodes, int kSampleTile, int kMicroVariants, int kMicroColumns,
          bool kDoubleBuffer, bool kAsyncDesign>
__device__ void fused_products_body(
    const uint8_t* __restrict__ packed, int64_t variants, int64_t samples,
    int64_t row_stride, const float* __restrict__ design, int64_t columns,
    const float* __restrict__ means, float* __restrict__ products) {
  constexpr int kVariantTile = kThreadRows * kMicroVariants;
  constexpr int kColumnTile = kThreadColumns * kMicroColumns;
  // Padded by four rather than one. One float is enough to break the bank
  // conflict on the staging write, but the inner loop reads its operands as
  // float4, which needs each row to start 16-byte aligned; four floats gives
  // both. The micro-tile widths are multiples of four, so a thread's slice of
  // a row is aligned too.
  constexpr int kGenotypeStride = kVariantTile + 4;
  // Sixteen two-bit calls per 32-bit word. Every configuration's sample tile
  // is a multiple of sixteen, so a tile is a whole number of words.
  static_assert(kSampleTile % 16 == 0, "sample tile must be a multiple of 16");
  constexpr int kWordsPerVariant = kSampleTile / 16;
  constexpr int kStagingWords = kVariantTile * kWordsPerVariant;
  constexpr int kDesignElements = kSampleTile * kColumnTile;
  // Per-thread slots for the prefetch registers. Both staging loops are
  // strided by the block size, so a thread owns a fixed set of elements and
  // can hold the next tile's values between iterations.
  constexpr int kWordSlots = (kStagingWords + kThreads - 1) / kThreads;
  constexpr int kDesignSlots = (kDesignElements + kThreads - 1) / kThreads;
  constexpr int kBuffers = kDoubleBuffer ? 2 : 1;
  constexpr int kGenotypeSpan = kSampleTile * kGenotypeStride;

  extern __shared__ float shared[];
  float* genotypes = shared;
  float* design_tile = genotypes + kBuffers * kGenotypeSpan;
  float* mean_tile = design_tile + kBuffers * kDesignElements;

  const int64_t variant_base = int64_t(blockIdx.x) * kVariantTile;
  const int64_t column_base = int64_t(blockIdx.y) * kColumnTile;
  const int variant_span = int(min(int64_t(kVariantTile), variants - variant_base));
  const int column_span = int(min(int64_t(kColumnTile), columns - column_base));
  if (variant_span <= 0 || column_span <= 0) return;

  for (int v = threadIdx.x; v < kVariantTile; v += kThreads)
    mean_tile[v] = v < variant_span ? means[variant_base + v] : 0.0f;
  // Every thread reads means written by other threads during staging. The
  // single-buffered loop opens with a barrier and so was ordered by accident;
  // the pipelined one stages before its first barrier and was not, which
  // showed up as a product wrong by a fifth of the matrix scale, differently
  // on each run.
  __syncthreads();

  const int thread_variant = int(threadIdx.x) / kThreadColumns;
  const int thread_column = int(threadIdx.x) % kThreadColumns;

  float accumulator[kMicroVariants][kMicroColumns];
#pragma unroll
  for (int a = 0; a < kMicroVariants; ++a)
#pragma unroll
    for (int b = 0; b < kMicroColumns; ++b) accumulator[a][b] = 0.0f;

  // The next tile's raw bytes and design values, held in registers across the
  // current tile's arithmetic. Issuing these loads before the maths and
  // consuming them after is what lets a global fetch overlap with compute:
  // a register is not read until the store, so the warp does not stall at the
  // load.
  uint32_t words[kWordSlots];
  float designs[kDesignSlots];

  // One 32-bit load per slot, decoded into sixteen calls at store time.
  // Loading a byte per (sample, variant) instead means four threads contend
  // for the same byte and the kernel issues sixteen loads where one would do;
  // two bits per call is exactly what makes a wide load worth unpacking.
  //
  // Consecutive threads take consecutive words of one variant's row, so a warp
  // walks that row rather than touching 32 rows row_stride apart.
  auto fetch = [&](int64_t sample0) {
#pragma unroll
    for (int slot = 0; slot < kWordSlots; ++slot) {
      const int index = int(threadIdx.x) + slot * kThreads;
      uint32_t word = 0;
      if (index < kStagingWords) {
        const int v = index / kWordsPerVariant;
        if (v < variant_span) {
          const int w = index % kWordsPerVariant;
          const uint8_t* row = packed + (variant_base + v) * row_stride;
          const int64_t byte0 = (sample0 >> 2) + int64_t(w) * 4;
          if (byte0 + 4 <= row_stride) {
            word = __ldg(reinterpret_cast<const uint32_t*>(row + byte0));
          } else {
            // Only a cohort whose packed width is not a multiple of four can
            // land here; fall back rather than over-read the row.
            for (int b = 0; b < 4 && byte0 + b < row_stride; ++b)
              word |= uint32_t(__ldg(row + byte0 + b)) << (8 * b);
          }
        }
      }
      words[slot] = word;
    }
    if (kAsyncDesign) return;  // the design arrives via cp.async instead
#pragma unroll
    for (int slot = 0; slot < kDesignSlots; ++slot) {
      const int index = int(threadIdx.x) + slot * kThreads;
      float value = 0.0f;
      if (index < kDesignElements) {
        const int s = index / kColumnTile;
        const int c = index % kColumnTile;
        if (sample0 + s < samples && c < column_span)
          value = __ldg(design + (sample0 + s) * columns + column_base + c);
      }
      designs[slot] = value;
    }
  };

  // Issue the design tile's copies without waiting for them. Out-of-range
  // elements are zero-filled by the copy itself rather than branched around,
  // so every thread issues the same number of copies and the pipeline stays
  // uniform; the source pointer is clamped to a valid address in that case
  // because a zero-filled copy still names one.
  auto async_design = [&](int buffer, int64_t sample0) {
    float* design_buffer = design_tile + buffer * kDesignElements;
#pragma unroll
    for (int slot = 0; slot < kDesignSlots; ++slot) {
      const int index = int(threadIdx.x) + slot * kThreads;
      if (index >= kDesignElements) break;
      const int s = index / kColumnTile;
      const int c = index % kColumnTile;
      const bool present = (sample0 + s < samples) && (c < column_span);
      const int64_t offset =
          present ? (sample0 + s) * columns + column_base + c : 0;
      __pipeline_memcpy_async(design_buffer + index, design + offset,
                              sizeof(float), present ? 0 : sizeof(float));
    }
    __pipeline_commit();
  };

  auto store = [&](int buffer, int span) {
    float* genotype_buffer = genotypes + buffer * kGenotypeSpan;
    float* design_buffer = design_tile + buffer * kDesignElements;
#pragma unroll
    for (int slot = 0; slot < kWordSlots; ++slot) {
      const int index = int(threadIdx.x) + slot * kThreads;
      if (index >= kStagingWords) break;
      const int v = index / kWordsPerVariant;
      const int w = index % kWordsPerVariant;
      const uint32_t word = words[slot];
      const float mean = mean_tile[v];
#pragma unroll
      for (int k = 0; k < 16; ++k) {
        const int s = w * 16 + k;
        float centred = 0.0f;
        if (s < span && v < variant_span) {
          const float value = decode_two_bit((word >> (2 * k)) & 3u, kCodes);
          centred = isnan(value) ? 0.0f : value - mean;
        }
        genotype_buffer[s * kGenotypeStride + v] = centred;
      }
    }
    if (kAsyncDesign) return;  // cp.async has already written this buffer
#pragma unroll
    for (int slot = 0; slot < kDesignSlots; ++slot) {
      const int index = int(threadIdx.x) + slot * kThreads;
      if (index >= kDesignElements) break;
      design_buffer[index] = designs[slot];
    }
  };

  auto accumulate = [&](int buffer, int span) {
    const float* genotype_buffer = genotypes + buffer * kGenotypeSpan;
    const float* design_buffer = design_tile + buffer * kDesignElements;
    // Four floats per instruction. The inner loop was issuing one shared-memory
    // instruction per operand -- twelve of them against thirty-two multiply-
    // adds at an 8x4 tile -- which makes shared issue, not arithmetic, the
    // limit. Reading float4 cuts those twelve instructions to three and moves
    // the same bytes.
    static_assert(kMicroVariants % 4 == 0 && kMicroColumns % 4 == 0,
                  "micro-tile widths must be multiples of four for float4");
    for (int s = 0; s < span; ++s) {
      float local_genotypes[kMicroVariants];
      float local_design[kMicroColumns];
#pragma unroll
      for (int a = 0; a < kMicroVariants; a += 4) {
        const float4 quad = *reinterpret_cast<const float4*>(
            genotype_buffer + s * kGenotypeStride +
            thread_variant * kMicroVariants + a);
        local_genotypes[a + 0] = quad.x;
        local_genotypes[a + 1] = quad.y;
        local_genotypes[a + 2] = quad.z;
        local_genotypes[a + 3] = quad.w;
      }
#pragma unroll
      for (int b = 0; b < kMicroColumns; b += 4) {
        const float4 quad = *reinterpret_cast<const float4*>(
            design_buffer + s * kColumnTile +
            thread_column * kMicroColumns + b);
        local_design[b + 0] = quad.x;
        local_design[b + 1] = quad.y;
        local_design[b + 2] = quad.z;
        local_design[b + 3] = quad.w;
      }
#pragma unroll
      for (int a = 0; a < kMicroVariants; ++a)
#pragma unroll
        for (int b = 0; b < kMicroColumns; ++b)
          accumulator[a][b] += local_genotypes[a] * local_design[b];
    }
  };

  const int64_t tiles = (samples + kSampleTile - 1) / kSampleTile;
  auto tile_span = [&](int64_t tile) {
    return int(min(int64_t(kSampleTile), samples - tile * kSampleTile));
  };

  if (kDoubleBuffer) {
    fetch(0);
    if (kAsyncDesign) async_design(0, 0);
    store(0, tile_span(0));
    if (kAsyncDesign) __pipeline_wait_prior(0);
    __syncthreads();
    for (int64_t tile = 0; tile < tiles; ++tile) {
      // Issue the next tile's loads first: they are in flight throughout the
      // arithmetic below, which is the whole point of holding them in
      // registers rather than going straight to shared memory. With
      // `cp.async` the design copy is in flight there too, and without
      // occupying registers to do it.
      if (tile + 1 < tiles) {
        fetch((tile + 1) * kSampleTile);
        if (kAsyncDesign)
          async_design(int((tile + 1) & 1), (tile + 1) * kSampleTile);
      }
      accumulate(int(tile & 1), tile_span(tile));
      if (tile + 1 < tiles) {
        store(int((tile + 1) & 1), tile_span(tile + 1));
        // Drain before the barrier so the next iteration reads a complete tile.
        if (kAsyncDesign) __pipeline_wait_prior(0);
      }
      // Writing the other buffer, so this only has to order the store against
      // the next iteration's reads of it.
      __syncthreads();
    }
  } else {
    for (int64_t tile = 0; tile < tiles; ++tile) {
      const int span = tile_span(tile);
      __syncthreads();
      fetch(tile * kSampleTile);
      store(0, span);
      __syncthreads();
      accumulate(0, span);
    }
  }

#pragma unroll
  for (int a = 0; a < kMicroVariants; ++a) {
    const int v = thread_variant * kMicroVariants + a;
    if (v >= variant_span) break;
#pragma unroll
    for (int b = 0; b < kMicroColumns; ++b) {
      const int c = thread_column * kMicroColumns + b;
      if (c >= column_span) break;
      products[(variant_base + v) * columns + column_base + c] =
          accumulator[a][b];
    }
  }
}

// Entry point with no occupancy promise: ptxas picks its own register target.
template <int kCodes, int kSampleTile, int kMicroVariants, int kMicroColumns,
          bool kDoubleBuffer, bool kAsyncDesign>
__global__ __launch_bounds__(kThreads) void fused_products_unhinted(
    const uint8_t* __restrict__ packed, int64_t variants, int64_t samples,
    int64_t row_stride, const float* __restrict__ design, int64_t columns,
    const float* __restrict__ means, float* __restrict__ products) {
  fused_products_body<kCodes, kSampleTile, kMicroVariants, kMicroColumns,
                      kDoubleBuffer, kAsyncDesign>(
      packed, variants, samples, row_stride, design, columns, means, products);
}

// Entry point promising `kMinBlocks` blocks per multiprocessor, which caps
// registers at 65536 / (256 * kMinBlocks).
template <int kCodes, int kSampleTile, int kMicroVariants, int kMicroColumns,
          bool kDoubleBuffer, int kMinBlocks, bool kAsyncDesign>
__global__ __launch_bounds__(kThreads, kMinBlocks) void fused_products_kernel(
    const uint8_t* __restrict__ packed, int64_t variants, int64_t samples,
    int64_t row_stride, const float* __restrict__ design, int64_t columns,
    const float* __restrict__ means, float* __restrict__ products) {
  fused_products_body<kCodes, kSampleTile, kMicroVariants, kMicroColumns,
                      kDoubleBuffer, kAsyncDesign>(
      packed, variants, samples, row_stride, design, columns, means, products);
}

// Warp-specialised variant: dedicated staging warps, dedicated FMA warps.
//
// The block runs 384 threads. The first 256 keep the 16 x 16 consumer grid and
// do nothing but accumulate; the last 128 do nothing but stage the next tile.
// One barrier per tile separates them, so a consumer never executes a load
// instruction and a producer never executes an FMA.
//
// This is the third way of overlapping staging with compute that this file has
// tried, and the first two both lost: register double-buffering was slower than
// single-buffering in every shape, and cp.async was mixed. Both lost on
// occupancy, and this one spends more of it -- two staging buffers *and* half
// again as many threads per block. Built because the prediction is worth
// testing rather than asserting.
template <int kCodes, int kSampleTile, int kMicroVariants, int kMicroColumns,
          int kMinBlocks>
__global__ __launch_bounds__(384, kMinBlocks) void fused_products_warpspec(
    const uint8_t* __restrict__ packed, int64_t variants, int64_t samples,
    int64_t row_stride, const float* __restrict__ design, int64_t columns,
    const float* __restrict__ means, float* __restrict__ products) {
  constexpr int kConsumers = 256;
  constexpr int kProducers = 128;
  constexpr int kVariantTile = kThreadRows * kMicroVariants;
  constexpr int kColumnTile = kThreadColumns * kMicroColumns;
  constexpr int kGenotypeStride = kVariantTile + 4;
  static_assert(kSampleTile % 16 == 0, "sample tile must be a multiple of 16");
  constexpr int kWordsPerVariant = kSampleTile / 16;
  constexpr int kStagingWords = kVariantTile * kWordsPerVariant;
  constexpr int kDesignElements = kSampleTile * kColumnTile;
  constexpr int kGenotypeSpan = kSampleTile * kGenotypeStride;

  extern __shared__ float shared[];
  float* genotypes = shared;
  float* design_tile = genotypes + 2 * kGenotypeSpan;
  float* mean_tile = design_tile + 2 * kDesignElements;

  const int64_t variant_base = int64_t(blockIdx.x) * kVariantTile;
  const int64_t column_base = int64_t(blockIdx.y) * kColumnTile;
  const int variant_span = int(min(int64_t(kVariantTile), variants - variant_base));
  const int column_span = int(min(int64_t(kColumnTile), columns - column_base));
  if (variant_span <= 0 || column_span <= 0) return;

  for (int v = threadIdx.x; v < kVariantTile; v += kConsumers + kProducers)
    mean_tile[v] = v < variant_span ? means[variant_base + v] : 0.0f;
  __syncthreads();

  const bool producer = threadIdx.x >= kConsumers;
  const int producer_lane = int(threadIdx.x) - kConsumers;

  auto stage = [&](int buffer, int64_t sample0, int span) {
    float* genotype_buffer = genotypes + buffer * kGenotypeSpan;
    float* design_buffer = design_tile + buffer * kDesignElements;
    for (int index = producer_lane; index < kStagingWords; index += kProducers) {
      const int v = index / kWordsPerVariant;
      const int w = index % kWordsPerVariant;
      uint32_t word = 0;
      if (v < variant_span) {
        const uint8_t* row = packed + (variant_base + v) * row_stride;
        const int64_t byte0 = (sample0 >> 2) + int64_t(w) * 4;
        if (byte0 + 4 <= row_stride)
          word = __ldg(reinterpret_cast<const uint32_t*>(row + byte0));
        else
          for (int b = 0; b < 4 && byte0 + b < row_stride; ++b)
            word |= uint32_t(__ldg(row + byte0 + b)) << (8 * b);
      }
      const float mean = mean_tile[v];
#pragma unroll
      for (int k = 0; k < 16; ++k) {
        const int s = w * 16 + k;
        float centred = 0.0f;
        if (s < span && v < variant_span) {
          const float value = decode_two_bit((word >> (2 * k)) & 3u, kCodes);
          centred = isnan(value) ? 0.0f : value - mean;
        }
        genotype_buffer[s * kGenotypeStride + v] = centred;
      }
    }
    for (int index = producer_lane; index < kDesignElements; index += kProducers) {
      const int s = index / kColumnTile;
      const int c = index % kColumnTile;
      design_buffer[index] =
          (sample0 + s < samples && c < column_span)
              ? __ldg(design + (sample0 + s) * columns + column_base + c)
              : 0.0f;
    }
  };

  const int thread_variant = int(threadIdx.x) / kThreadColumns;
  const int thread_column = int(threadIdx.x) % kThreadColumns;
  float accumulator[kMicroVariants][kMicroColumns];
#pragma unroll
  for (int a = 0; a < kMicroVariants; ++a)
#pragma unroll
    for (int b = 0; b < kMicroColumns; ++b) accumulator[a][b] = 0.0f;

  const int64_t tiles = (samples + kSampleTile - 1) / kSampleTile;
  auto tile_span = [&](int64_t tile) {
    return int(min(int64_t(kSampleTile), samples - tile * kSampleTile));
  };

  if (producer) stage(0, 0, tile_span(0));
  __syncthreads();

  for (int64_t tile = 0; tile < tiles; ++tile) {
    if (producer) {
      if (tile + 1 < tiles)
        stage(int((tile + 1) & 1), (tile + 1) * kSampleTile, tile_span(tile + 1));
    } else {
      const float* genotype_buffer = genotypes + int(tile & 1) * kGenotypeSpan;
      const float* design_buffer = design_tile + int(tile & 1) * kDesignElements;
      const int span = tile_span(tile);
      for (int s = 0; s < span; ++s) {
        float local_genotypes[kMicroVariants];
        float local_design[kMicroColumns];
#pragma unroll
        for (int a = 0; a < kMicroVariants; a += 4) {
          const float4 quad = *reinterpret_cast<const float4*>(
              genotype_buffer + s * kGenotypeStride +
              thread_variant * kMicroVariants + a);
          local_genotypes[a + 0] = quad.x;
          local_genotypes[a + 1] = quad.y;
          local_genotypes[a + 2] = quad.z;
          local_genotypes[a + 3] = quad.w;
        }
#pragma unroll
        for (int b = 0; b < kMicroColumns; b += 4) {
          const float4 quad = *reinterpret_cast<const float4*>(
              design_buffer + s * kColumnTile + thread_column * kMicroColumns + b);
          local_design[b + 0] = quad.x;
          local_design[b + 1] = quad.y;
          local_design[b + 2] = quad.z;
          local_design[b + 3] = quad.w;
        }
#pragma unroll
        for (int a = 0; a < kMicroVariants; ++a)
#pragma unroll
          for (int b = 0; b < kMicroColumns; ++b)
            accumulator[a][b] += local_genotypes[a] * local_design[b];
      }
    }
    __syncthreads();
  }

  if (producer) return;
#pragma unroll
  for (int a = 0; a < kMicroVariants; ++a) {
    const int v = thread_variant * kMicroVariants + a;
    if (v >= variant_span) break;
#pragma unroll
    for (int b = 0; b < kMicroColumns; ++b) {
      const int c = thread_column * kMicroColumns + b;
      if (c >= column_span) break;
      products[(variant_base + v) * columns + column_base + c] =
          accumulator[a][b];
    }
  }
}

int checked_launch() {
  auto error = cudaGetLastError();
  if (error != cudaSuccess) { last_error = cudaGetErrorString(error); return -1; }
  return 0;
}

template <int kSampleTile, int kMicroVariants, int kMicroColumns, int kMinBlocks>
int launch_warpspec(const uint8_t* packed, int64_t variants, int64_t samples,
                    int64_t row_stride, int encoding, const float* design,
                    int64_t columns, const float* means, float* products,
                    cudaStream_t stream) {
  constexpr int kVariantTile = kThreadRows * kMicroVariants;
  constexpr int kColumnTile = kThreadColumns * kMicroColumns;
  constexpr size_t bytes =
      sizeof(float) * (2u * kSampleTile * (kVariantTile + 4) +
                       2u * size_t(kSampleTile) * kColumnTile + kVariantTile);
  dim3 blocks(unsigned((variants + kVariantTile - 1) / kVariantTile),
              unsigned((columns + kColumnTile - 1) / kColumnTile));
  auto prepare = [](const void* function) {
    if (bytes <= 48u * 1024u) return true;
    return cudaFuncSetAttribute(function,
                                cudaFuncAttributeMaxDynamicSharedMemorySize,
                                int(bytes)) == cudaSuccess;
  };
  if (encoding == kPackedPlink1) {
    auto kernel = fused_products_warpspec<kPackedPlink1, kSampleTile,
                                          kMicroVariants, kMicroColumns,
                                          kMinBlocks>;
    if (!prepare(reinterpret_cast<const void*>(kernel))) {
      last_error = "device refused the shared memory this configuration needs";
      return -1;
    }
    kernel<<<blocks, 384, bytes, stream>>>(packed, variants, samples, row_stride,
                                           design, columns, means, products);
  } else {
    auto kernel = fused_products_warpspec<kPackedPgen, kSampleTile,
                                          kMicroVariants, kMicroColumns,
                                          kMinBlocks>;
    if (!prepare(reinterpret_cast<const void*>(kernel))) {
      last_error = "device refused the shared memory this configuration needs";
      return -1;
    }
    kernel<<<blocks, 384, bytes, stream>>>(packed, variants, samples, row_stride,
                                           design, columns, means, products);
  }
  return checked_launch();
}

// kMinBlocks defaults to 0 = no occupancy hint, which is what the baseline
// configurations are supposed to measure. It was briefly 1, which is a real
// (and weaker) constraint rather than the absence of one, and that silently
// re-baselined every unhinted configuration in the sweep.
template <int kSampleTile, int kMicroVariants, int kMicroColumns,
          bool kDoubleBuffer, int kMinBlocks = 0, bool kAsyncDesign = false>
int launch(const uint8_t* packed, int64_t variants, int64_t samples,
           int64_t row_stride, int encoding, const float* design,
           int64_t columns, const float* means, float* products,
           cudaStream_t stream) {
  constexpr int kVariantTile = kThreadRows * kMicroVariants;
  constexpr int kColumnTile = kThreadColumns * kMicroColumns;
  constexpr int kBuffers = kDoubleBuffer ? 2 : 1;
  constexpr size_t bytes =
      sizeof(float) * (size_t(kBuffers) * kSampleTile * (kVariantTile + 4) +
                       size_t(kBuffers) * kSampleTile * kColumnTile +
                       kVariantTile);
  dim3 blocks(unsigned((variants + kVariantTile - 1) / kVariantTile),
              unsigned((columns + kColumnTile - 1) / kColumnTile));

  // Past 48 KB a kernel must ask for the larger dynamic allocation explicitly,
  // and the request can be refused on a device that does not have it.
  auto prepare = [](const void* function) {
    if (bytes <= 48u * 1024u) return true;
    return cudaFuncSetAttribute(function,
                                cudaFuncAttributeMaxDynamicSharedMemorySize,
                                int(bytes)) == cudaSuccess;
  };

  if (encoding == kPackedPlink1) {
    // kMinBlocks == 0 means "no occupancy promise", which needs the wrapper
    // declared without a second launch-bounds argument -- passing 1 is a
    // different and weaker constraint, not the absence of one.
    auto kernel = kMinBlocks == 0
        ? fused_products_unhinted<kPackedPlink1, kSampleTile, kMicroVariants,
                                  kMicroColumns, kDoubleBuffer, kAsyncDesign>
        : fused_products_kernel<kPackedPlink1, kSampleTile,
                                kMicroVariants, kMicroColumns,
                                kDoubleBuffer, kMinBlocks == 0 ? 1 : kMinBlocks,
                                kAsyncDesign>;
    if (!prepare(reinterpret_cast<const void*>(kernel))) {
      last_error = "device refused the shared memory this configuration needs";
      return -1;
    }
    kernel<<<blocks, kThreads, bytes, stream>>>(packed, variants, samples,
                                                row_stride, design, columns,
                                                means, products);
  } else {
    auto kernel = kMinBlocks == 0
        ? fused_products_unhinted<kPackedPgen, kSampleTile, kMicroVariants,
                                  kMicroColumns, kDoubleBuffer, kAsyncDesign>
        : fused_products_kernel<kPackedPgen, kSampleTile,
                                kMicroVariants, kMicroColumns,
                                kDoubleBuffer, kMinBlocks == 0 ? 1 : kMinBlocks,
                                kAsyncDesign>;
    if (!prepare(reinterpret_cast<const void*>(kernel))) {
      last_error = "device refused the shared memory this configuration needs";
      return -1;
    }
    kernel<<<blocks, kThreads, bytes, stream>>>(packed, variants, samples,
                                                row_stride, design, columns,
                                                means, products);
  }
  return checked_launch();
}

}  // namespace

extern "C" {

int tg_fused_abi_version() { return 2; }

const char* tg_fused_error() { return last_error.c_str(); }

// The configurations the sweep can reach, addressed by index rather than by a
// tuple so the harness and the library cannot disagree about what was run.
int tg_fused_config_count() { return 30; }

const char* tg_fused_config_name(int config) {
  switch (config) {
    case 0: return "sample16_micro4x4";
    case 1: return "sample32_micro4x4";
    case 2: return "sample64_micro4x4";
    case 3: return "sample128_micro4x4";
    case 4: return "sample32_micro8x4";
    case 5: return "sample64_micro8x4";
    case 6: return "sample32_micro8x8";
    case 7: return "sample64_micro8x8";
    // Double-buffered. These trade shared memory for overlap: two staging
    // buffers instead of one, so the sample tile that fits is smaller, and
    // whether the overlap is worth the occupancy is exactly what the sweep
    // is for.
    case 8: return "sample32_micro4x4_pipelined";
    case 9: return "sample64_micro4x4_pipelined";
    case 10: return "sample128_micro4x4_pipelined";
    case 11: return "sample32_micro8x4_pipelined";
    case 12: return "sample64_micro8x4_pipelined";
    case 13: return "sample32_micro8x8_pipelined";
    // Occupancy-forced variants of the configurations that won without a cap.
    // The suffix is the promised blocks per multiprocessor, which is what caps
    // the register budget.
    case 14: return "sample32_micro8x4_occ3";
    case 15: return "sample32_micro8x4_occ4";
    case 16: return "sample64_micro8x4_occ3";
    case 17: return "sample64_micro8x4_occ4";
    case 18: return "sample32_micro8x8_occ3";
    case 19: return "sample32_micro8x8_occ4";
    case 20: return "sample64_micro4x4_occ4";
    case 21: return "sample64_micro4x4_occ6";
    // cp.async staging of the design tile, double buffered.
    case 22: return "sample32_micro8x4_async";
    case 23: return "sample64_micro8x4_async";
    case 24: return "sample32_micro8x4_async_occ3";
    case 25: return "sample32_micro8x8_async";
    // Warp-specialised: 256 consumer threads, 128 dedicated staging threads.
    case 26: return "sample32_micro8x4_warpspec";
    case 27: return "sample64_micro8x4_warpspec";
    case 28: return "sample32_micro8x8_warpspec";
    case 29: return "sample32_micro8x4_warpspec_occ2";
    default: return "";
  }
}

int tg_fused_products_config(const uint8_t* packed, int64_t variants,
                             int64_t samples, int64_t row_stride, int encoding,
                             const float* design, int64_t columns,
                             const float* means, int config, float* products,
                             void* stream) {
  if (variants <= 0 || samples <= 0 || columns <= 0 ||
      row_stride < (samples + 3) / 4) {
    last_error = "invalid fused product dimensions";
    return -1;
  }
  cudaStream_t s = static_cast<cudaStream_t>(stream);
#define TG_CASE(index, st, mv, mc, db)                                      \
  case index:                                                               \
    return launch<st, mv, mc, db>(packed, variants, samples, row_stride,    \
                                  encoding, design, columns, means,         \
                                  products, s)
#define TG_CASE_OCC(index, st, mv, mc, db, occ)                             \
  case index:                                                               \
    return launch<st, mv, mc, db, occ>(packed, variants, samples,           \
                                       row_stride, encoding, design,        \
                                       columns, means, products, s)
  switch (config) {
    TG_CASE(0, 16, 4, 4, false);
    TG_CASE(1, 32, 4, 4, false);
    TG_CASE(2, 64, 4, 4, false);
    TG_CASE(3, 128, 4, 4, false);
    TG_CASE(4, 32, 8, 4, false);
    TG_CASE(5, 64, 8, 4, false);
    TG_CASE(6, 32, 8, 8, false);
    TG_CASE(7, 64, 8, 8, false);
    TG_CASE(8, 32, 4, 4, true);
    TG_CASE(9, 64, 4, 4, true);
    TG_CASE(10, 128, 4, 4, true);
    TG_CASE(11, 32, 8, 4, true);
    TG_CASE(12, 64, 8, 4, true);
    TG_CASE(13, 32, 8, 8, true);
    TG_CASE_OCC(14, 32, 8, 4, false, 3);
    TG_CASE_OCC(15, 32, 8, 4, false, 4);
    TG_CASE_OCC(16, 64, 8, 4, false, 3);
    TG_CASE_OCC(17, 64, 8, 4, false, 4);
    TG_CASE_OCC(18, 32, 8, 8, false, 3);
    TG_CASE_OCC(19, 32, 8, 8, false, 4);
    TG_CASE_OCC(20, 64, 4, 4, false, 4);
    TG_CASE_OCC(21, 64, 4, 4, false, 6);
#define TG_CASE_ASYNC(index, st, mv, mc, occ)                               \
  case index:                                                               \
    return launch<st, mv, mc, true, occ, true>(packed, variants, samples,   \
                                               row_stride, encoding, design,\
                                               columns, means, products, s)
    TG_CASE_ASYNC(22, 32, 8, 4, 0);
    TG_CASE_ASYNC(23, 64, 8, 4, 0);
    TG_CASE_ASYNC(24, 32, 8, 4, 3);
    TG_CASE_ASYNC(25, 32, 8, 8, 0);
#undef TG_CASE_ASYNC
#define TG_CASE_WARPSPEC(index, st, mv, mc, occ)                             \
  case index:                                                                \
    return launch_warpspec<st, mv, mc, occ>(packed, variants, samples,       \
                                            row_stride, encoding, design,    \
                                            columns, means, products, s)
    TG_CASE_WARPSPEC(26, 32, 8, 4, 1);
    TG_CASE_WARPSPEC(27, 64, 8, 4, 1);
    TG_CASE_WARPSPEC(28, 32, 8, 8, 1);
    TG_CASE_WARPSPEC(29, 32, 8, 4, 2);
#undef TG_CASE_WARPSPEC
    default:
      last_error = "unknown fused configuration";
      return -1;
  }
#undef TG_CASE
#undef TG_CASE_OCC
}

// The default the scan would use. A separate entry so callers that do not care
// about tuning need not name a configuration, and so the choice lives here
// rather than being repeated in each caller.
int tg_fused_products(const uint8_t* packed, int64_t variants, int64_t samples,
                      int64_t row_stride, int encoding, const float* design,
                      int64_t columns, const float* means, int sample_tile,
                      float* products, void* stream) {
  (void)sample_tile;  // superseded by the configuration index
  // Configuration 14, `sample32_micro8x4_occ3`: the fastest of the thirty on an
  // idle H100 at K=512 (26.83 TFLOP/s). The previous default was 5, which the
  // same sweep puts at 24.96. Chosen by measurement, not by shape.
  return tg_fused_products_config(packed, variants, samples, row_stride,
                                  encoding, design, columns, means, 14,
                                  products, stream);
}
}
