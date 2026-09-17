/* Direct PGEN hard-call record decoding.
 *
 * Decodes a contiguous run of PGEN variant records from a caller-supplied
 * byte buffer into packed two-bit genotype vectors, one padded row per
 * variant. This is the hot half of src/torchgwas/pgen_reader.py, which
 * remains the reference implementation and the thing this file is validated
 * against.
 *
 * File input stays on the Python side: the existing PGEN path already has
 * parallel positional reads and pinned buffers, and ctypes releases the GIL
 * for the duration of the call, so several threads can decode disjoint
 * variant runs concurrently.
 *
 * Record forms and the two undocumented details this relies on are described
 * in docs/direct_pgen_reader.md. In short: difflist per-group sizes exist for
 * skipping and must not be used to advance a sequential decode, and the
 * one-bit mode byte encodes an unordered category pair as 3*low + high.
 *
 * C ABI version 1. Returns 0 on success, or a negative code naming the
 * failure so the caller can distinguish "unsupported input" from "corrupt
 * input" rather than guessing.
 */

#include <stdint.h>
#include <string.h>

#define TORCHGWAS_PGEN_ABI 1

/* Error codes. Kept distinct so the Python wrapper can raise the right type. */
#define PGEN_OK 0
#define PGEN_ERR_TRUNCATED (-1)
#define PGEN_ERR_UNSUPPORTED_TRACK (-2)
#define PGEN_ERR_RESERVED_FORM (-3)
#define PGEN_ERR_SAMPLE_OUT_OF_RANGE (-4)
#define PGEN_ERR_NO_LD_BASE (-5)
#define PGEN_ERR_BAD_MODE_BYTE (-6)
#define PGEN_ERR_ULEB_OVERFLOW (-7)

#define TRACK_MULTIALLELIC 0x08u
#define TRACK_PHASE 0x30u
#define TRACK_DOSAGE 0xC0u
#define UNSUPPORTED_TRACKS (TRACK_MULTIALLELIC | TRACK_PHASE | TRACK_DOSAGE)

#define DIFFLIST_GROUP 64u

typedef struct {
    const uint8_t *data;
    uint64_t size;
    uint64_t pos;
} Cursor;

/* Expand the low eight bits of `bits` so that bit i occupies the whole of
 * two-bit slot i: 0b1011 becomes 0b00_11_00_11_11. The first three steps are
 * the standard Morton spread, and the multiply by three widens each isolated
 * bit to the pair above it, which cannot carry because the bits are isolated.
 *
 * This turns one byte of a one-bit record into two genovec bytes with no
 * branch and no per-sample read-modify-write. See the note in `decode_record`
 * for why that form was worth attacking. */
static inline uint16_t spread_bits8(uint8_t bits) {
    uint32_t x = bits;
    x = (x | (x << 4)) & 0x0F0Fu;
    x = (x | (x << 2)) & 0x3333u;
    x = (x | (x << 1)) & 0x5555u;
    return (uint16_t)(x * 3u);
}

static inline void set_category(uint8_t *row, uint64_t sample, uint8_t category) {
    const uint64_t byte = sample >> 2;
    const unsigned shift = (unsigned)((sample & 3u) * 2u);
    row[byte] = (uint8_t)((row[byte] & ~(0x03u << shift)) | ((category & 0x03u) << shift));
}

static inline uint8_t get_category(const uint8_t *row, uint64_t sample) {
    return (uint8_t)((row[sample >> 2] >> ((sample & 3u) * 2u)) & 0x03u);
}

/* Apply a sparse (sample, category) list onto an already-populated row.
 *
 * with_genotypes is always true for the hard-call forms; the flag exists so
 * the layout stays explicit rather than implied by the call site.
 *
 * The cursor is unpacked into locals for the duration of the call and written
 * back once on the way out, and the varint reader is spelled inline rather
 * than called. That is not a style choice. This library is built with
 * -fno-strict-aliasing, so the compiler must assume the two-bit store in
 * `set_category` can land on `cursor->pos`, `cursor->data` and `cursor->size`;
 * with the cursor reached through a pointer, every single difflist entry
 * reloaded all three from memory immediately after storing to the row. Keeping
 * the read position in a local whose address is never taken removes that, and
 * `restrict` on the row states what is already true -- the output block and
 * the record bytes are different allocations -- so the genotype and group-base
 * pointers stay in registers across the stores as well.
 */
static int apply_difflist(Cursor *cursor, uint8_t *restrict row, uint64_t sample_ct,
                          uint32_t id_bytes, int with_genotypes) {
    const uint8_t *const data = cursor->data;
    const uint8_t *const end = data + cursor->size;
    const uint8_t *p = data + cursor->pos;

/* Read one ULEB128 into `target`, or leave the function with the right error
 * after parking the read position. Deliberately a macro: passing `&p` to a
 * helper would let the address escape and reinstate exactly the reloads this
 * arrangement exists to remove.
 *
 * Peeling a single-byte fast path out of this loop was tried against the old
 * shape and reverted: 45.3 ms against 46.3 ms on a 3,000-variant block, inside
 * the noise. The deltas really are almost all one byte, but gcc already folds
 * the first iteration -- the shift is zero there and the overflow test is dead
 * -- so there was nothing left to peel. The cost was never the varint. */
#define DIFFLIST_ULEB128(target)                                              \
    do {                                                                      \
        unsigned shift_ = 0;                                                  \
        (target) = 0;                                                         \
        for (;;) {                                                            \
            if (p >= end) {                                                   \
                cursor->pos = (uint64_t)(p - data);                           \
                return PGEN_ERR_TRUNCATED;                                    \
            }                                                                 \
            const uint8_t byte_ = *p++;                                       \
            (target) |= (uint64_t)(byte_ & 0x7Fu) << shift_;                  \
            if ((byte_ & 0x80u) == 0u) {                                      \
                break;                                                        \
            }                                                                 \
            shift_ += 7u;                                                     \
            if (shift_ > 63u) {                                               \
                cursor->pos = (uint64_t)(p - data);                           \
                return PGEN_ERR_ULEB_OVERFLOW;                                \
            }                                                                 \
        }                                                                     \
    } while (0)

    uint64_t entry_ct;
    DIFFLIST_ULEB128(entry_ct);
    if (entry_ct == 0) {
        cursor->pos = (uint64_t)(p - data);
        return PGEN_OK;
    }
    if (entry_ct > sample_ct) {
        cursor->pos = (uint64_t)(p - data);
        return PGEN_ERR_SAMPLE_OUT_OF_RANGE;
    }
    const uint64_t group_ct = (entry_ct + DIFFLIST_GROUP - 1u) / DIFFLIST_GROUP;

    const uint8_t *const bases = p;
    if ((uint64_t)(end - p) < group_ct * (uint64_t)id_bytes) {
        cursor->pos = (uint64_t)(p - data);
        return PGEN_ERR_TRUNCATED;
    }
    p += group_ct * id_bytes;

    /* Per-group delta-run sizes: present for every group but the last, and
     * deliberately skipped. See the note at the top of this file. */
    if (group_ct > 1u) {
        if ((uint64_t)(end - p) < group_ct - 1u) {
            cursor->pos = (uint64_t)(p - data);
            return PGEN_ERR_TRUNCATED;
        }
        p += group_ct - 1u;
    }

    const uint8_t *genotypes = NULL;
    if (with_genotypes) {
        const uint64_t packed = (entry_ct + 3u) / 4u;
        if ((uint64_t)(end - p) < packed) {
            cursor->pos = (uint64_t)(p - data);
            return PGEN_ERR_TRUNCATED;
        }
        genotypes = p;
        p += packed;
    }

    uint64_t written = 0;
    for (uint64_t group = 0; group < group_ct; ++group) {
        const uint8_t *base_ptr = bases + group * id_bytes;
        uint64_t sample = 0;
        for (uint32_t b = 0; b < id_bytes; ++b) {
            sample |= (uint64_t)base_ptr[b] << (8u * b);
        }
        const uint64_t remaining = entry_ct - written;
        const uint64_t in_group = remaining < DIFFLIST_GROUP ? remaining : DIFFLIST_GROUP;
        for (uint64_t k = 0; k < in_group; ++k) {
            if (k != 0) {
                uint64_t delta;
                DIFFLIST_ULEB128(delta);
                sample += delta;
            }
            if (sample >= sample_ct) {
                cursor->pos = (uint64_t)(p - data);
                return PGEN_ERR_SAMPLE_OUT_OF_RANGE;
            }
            uint8_t category = 0;
            if (genotypes != NULL) {
                const uint64_t index = written + k;
                category = (uint8_t)((genotypes[index >> 2] >> ((index & 3u) * 2u)) & 0x03u);
            }
            set_category(row, sample, category);
        }
        written += in_group;
    }
    cursor->pos = (uint64_t)(p - data);
    return PGEN_OK;

#undef DIFFLIST_ULEB128
}

static int decode_record(const uint8_t *record, uint64_t record_len, uint8_t vrtype,
                         uint64_t sample_ct, uint32_t id_bytes,
                         const uint8_t *ld_base, int have_ld_base,
                         uint8_t *restrict row) {
    if ((vrtype & UNSUPPORTED_TRACKS) != 0u) {
        return PGEN_ERR_UNSUPPORTED_TRACK;
    }
    const uint64_t genovec_bytes = (sample_ct + 3u) / 4u;
    const unsigned form = vrtype & 0x07u;
    Cursor cursor = {record, record_len, 0};

    if (form == 0u) {
        if (record_len < genovec_bytes) {
            return PGEN_ERR_TRUNCATED;
        }
        memcpy(row, record, genovec_bytes);
        return PGEN_OK;
    }

    if (form == 4u || form == 6u || form == 7u) {
        /* Only the exceptions to one background category are stored. */
        const uint8_t background = (form == 4u) ? 0u : ((form == 6u) ? 2u : 3u);
        const uint8_t fill = (uint8_t)(background * 0x55u);  /* 0x00, 0xAA, 0xFF */
        memset(row, fill, (size_t)genovec_bytes);
        return apply_difflist(&cursor, row, sample_ct, id_bytes, 1);
    }

    if (form == 2u || form == 3u) {
        if (!have_ld_base) {
            return PGEN_ERR_NO_LD_BASE;
        }
        memcpy(row, ld_base, (size_t)genovec_bytes);
        const int status = apply_difflist(&cursor, row, sample_ct, id_bytes, 1);
        if (status != PGEN_OK) {
            return status;
        }
        if (form == 3u) {
            /* Categories 0 and 2 swap; 1 and 3 are untouched. Both are even
             * codes differing only in bit 1, so flip that bit wherever bit 0
             * is clear, a byte at a time. */
            for (uint64_t i = 0; i < genovec_bytes; ++i) {
                const uint8_t v = row[i];
                const uint8_t low_clear = (uint8_t)(~v & 0x55u);
                row[i] = (uint8_t)(v ^ (uint8_t)(low_clear << 1));
            }
        }
        return PGEN_OK;
    }

    if (form == 1u) {
        if (record_len < 1u) {
            return PGEN_ERR_TRUNCATED;
        }
        const uint8_t mode = record[0];
        if (mode == 0u) {
            return PGEN_ERR_BAD_MODE_BYTE;
        }
        const uint8_t low = (uint8_t)((mode - 1u) / 3u);
        const uint8_t high = (uint8_t)(mode - 3u * low);
        if (!(low < high && high <= 3u)) {
            return PGEN_ERR_BAD_MODE_BYTE;
        }
        const uint64_t bitarray_bytes = (sample_ct + 7u) / 8u;
        if (record_len < 1u + bitarray_bytes) {
            return PGEN_ERR_TRUNCATED;
        }
        const uint8_t *bits = record + 1;
        const uint8_t low_fill = (uint8_t)(low * 0x55u);
        /* This is the one form that touches every sample rather than a sparse
         * list, and on the real cohort it is 21% of records -- so the obvious
         * spelling, a per-sample bit test feeding a read-modify-write, ran the
         * full 22,250 samples for a fifth of the file and dominated the
         * decode. Eight samples at a time instead: spread the bit byte into
         * two-bit slots, then select between the two fill patterns with it.
         * Nothing is read back, and the write is two whole bytes.
         *
         * Byte order is written out explicitly rather than storing a uint16,
         * because the genovec's sample order is defined by the format and must
         * not follow the host's. */
        const uint16_t low_pair = (uint16_t)(low * 0x5555u);
        const uint16_t high_pair = (uint16_t)(high * 0x5555u);
        const uint64_t whole_bytes = sample_ct >> 3;
        for (uint64_t g = 0; g < whole_bytes; ++g) {
            const uint16_t select = spread_bits8(bits[g]);
            const uint16_t word = (uint16_t)((low_pair & (uint16_t)~select)
                                             | (high_pair & select));
            row[2u * g] = (uint8_t)word;
            row[2u * g + 1u] = (uint8_t)(word >> 8);
        }
        /* The tail is the samples past the last whole bit byte; two whole
         * genovec bytes per bit byte can never run past genovec_bytes. */
        memset(row + 2u * whole_bytes, low_fill,
               (size_t)(genovec_bytes - 2u * whole_bytes));
        for (uint64_t sample = whole_bytes << 3; sample < sample_ct; ++sample) {
            if ((bits[sample >> 3] >> (sample & 7u)) & 1u) {
                set_category(row, sample, high);
            }
        }
        cursor.pos = 1u + bitarray_bytes;
        return apply_difflist(&cursor, row, sample_ct, id_bytes, 1);
    }

    return PGEN_ERR_RESERVED_FORM;
}

/* --- public C ABI ------------------------------------------------------- */

int torchgwas_pgen_abi(void) { return TORCHGWAS_PGEN_ABI; }

/* Decode `count` consecutive records.
 *
 * records/offsets/lengths describe an in-memory span of the .pgen; offsets
 * are relative to `records`. out holds `count` rows of `out_stride` bytes.
 * ld_base is read for the first LD-compressed record and updated as
 * non-LD-compressed records are decoded, so a caller can resume a run.
 * have_ld_base says whether ld_base is already populated; on return it says
 * whether it is.
 *
 * failed_index receives the record that failed, so the caller can report the
 * variant rather than just the reason.
 */
int torchgwas_pgen_decode_range(const uint8_t *records,
                                const uint64_t *offsets,
                                const uint32_t *lengths,
                                const uint8_t *vrtypes,
                                uint64_t count,
                                uint64_t sample_ct,
                                uint32_t id_bytes,
                                uint8_t *out,
                                uint64_t out_stride,
                                uint8_t *ld_base,
                                int *have_ld_base,
                                uint64_t *failed_index) {
    const uint64_t genovec_bytes = (sample_ct + 3u) / 4u;
    if (out_stride < genovec_bytes || id_bytes == 0u || id_bytes > 4u) {
        if (failed_index) *failed_index = 0;
        return PGEN_ERR_TRUNCATED;
    }
    int have_base = have_ld_base ? *have_ld_base : 0;
    for (uint64_t i = 0; i < count; ++i) {
        uint8_t *row = out + i * out_stride;
        /* Padding beyond the logical samples is zeroed so a consumer reading
         * whole 64-byte words sees hom-ref rather than stale bytes. */
        if (out_stride > genovec_bytes) {
            memset(row + genovec_bytes, 0, (size_t)(out_stride - genovec_bytes));
        }
        const uint8_t vrtype = vrtypes[i];
        const int status = decode_record(records + offsets[i], lengths[i], vrtype,
                                         sample_ct, id_bytes, ld_base, have_base, row);
        if (status != PGEN_OK) {
            if (failed_index) *failed_index = i;
            if (have_ld_base) *have_ld_base = have_base;
            return status;
        }
        const unsigned form = vrtype & 0x07u;
        if (form != 2u && form != 3u) {
            memcpy(ld_base, row, (size_t)genovec_bytes);
            have_base = 1;
        }
    }
    if (have_ld_base) *have_ld_base = have_base;
    return PGEN_OK;
}

/* Expand packed rows straight to signed hard calls: 0, 1, 2, and `missing`
 * for code 3.
 *
 * This exists because doing the remap afterwards in numpy is the dominant cost
 * of a native read. `np.take` with a 256-entry table is a gather, one
 * random-access lookup per element, and it measured 468 ms of a 551 ms read --
 * 85% -- to move 133 MB, which is 0.28 GB/s. Emitting the final value here
 * costs nothing over emitting the category: the shift and mask are already
 * done, and a four-entry table lives in a register. It also removes the wide
 * intermediate array entirely, so the caller allocates one block instead of
 * two.
 *
 * Sample subsetting is not handled here; a caller wanting a subset still takes
 * the category path and gathers columns afterwards.
 */
int torchgwas_pgen_expand_hardcall(const uint8_t *packed, uint64_t count,
                                   uint64_t stride, uint64_t sample_ct,
                                   int8_t *out, int8_t missing) {
    const int8_t table[4] = {0, 1, 2, missing};
    const uint64_t whole = sample_ct >> 2;
    const uint64_t remainder = sample_ct & 3u;
    for (uint64_t i = 0; i < count; ++i) {
        const uint8_t *row = packed + i * stride;
        int8_t *dst = out + i * sample_ct;
        for (uint64_t byte = 0; byte < whole; ++byte) {
            const unsigned bits = row[byte];
            dst[0] = table[bits & 0x03u];
            dst[1] = table[(bits >> 2) & 0x03u];
            dst[2] = table[(bits >> 4) & 0x03u];
            dst[3] = table[(bits >> 6) & 0x03u];
            dst += 4;
        }
        if (remainder) {
            const unsigned bits = row[whole];
            for (uint64_t k = 0; k < remainder; ++k) {
                dst[k] = table[(bits >> (2u * k)) & 0x03u];
            }
        }
    }
    return PGEN_OK;
}

/* Expand packed rows to one byte per sample, for validation and for callers
 * that still want the wide form. */
int torchgwas_pgen_expand(const uint8_t *packed, uint64_t count, uint64_t stride,
                          uint64_t sample_ct, uint8_t *out) {
    /* Four samples per byte, unrolled.
     *
     * The scalar loop this replaces called get_category once per sample, and
     * each call re-derived a byte index and a shift from the sample number: for
     * a 3,000 x 22,250 block that is 66.75 million index computations to unpack
     * 16.7 million bytes. Reading each byte once and emitting its four codes
     * lets the compiler hold it in a register and vectorise the stores. The
     * tail covers a sample count that is not a multiple of four.
     */
    const uint64_t whole = sample_ct >> 2;
    const uint64_t remainder = sample_ct & 3u;
    for (uint64_t i = 0; i < count; ++i) {
        const uint8_t *row = packed + i * stride;
        uint8_t *dst = out + i * sample_ct;
        for (uint64_t byte = 0; byte < whole; ++byte) {
            const unsigned bits = row[byte];
            dst[0] = (uint8_t)(bits & 0x03u);
            dst[1] = (uint8_t)((bits >> 2) & 0x03u);
            dst[2] = (uint8_t)((bits >> 4) & 0x03u);
            dst[3] = (uint8_t)((bits >> 6) & 0x03u);
            dst += 4;
        }
        if (remainder) {
            const unsigned bits = row[whole];
            for (uint64_t k = 0; k < remainder; ++k) {
                dst[k] = (uint8_t)((bits >> (2u * k)) & 0x03u);
            }
        }
    }
    return PGEN_OK;
}
