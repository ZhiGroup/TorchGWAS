// CPU-only fused row-threshold predicate, compiled without fast-math.
// No ISA-specific intrinsics, CUDA code, shared state or allocations.
#include <cmath>
#include <cstddef>
#include <cstdint>

#ifndef TORCHGWAS_PREDICATE_BUILD_KEY
#error "Build key required"
#endif
#define EXPORT extern "C" __attribute__((visibility("default")))

EXPORT const char* torchgwas_predicate_build_key() {
    return TORCHGWAS_PREDICATE_BUILD_KEY;
}

EXPORT void torchgwas_predicate(const float* values, const float* critical,
                               std::uint8_t* mask, std::size_t rows,
                               std::size_t columns, std::ptrdiff_t row_stride) {
    for (std::size_t row = 0; row < rows; ++row) {
        const float limit = critical[static_cast<std::ptrdiff_t>(row) * row_stride];
        for (std::size_t column = 0; column < columns; ++column) {
            const float absolute = std::fabs(values[row * columns + column]);
            // Boolean bitwise AND avoids a data-dependent branch and permits
            // compiler vectorization without changing IEEE comparison rules.
            mask[row * columns + column] = std::isfinite(absolute) & (absolute >= limit);
        }
    }
}
