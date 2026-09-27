// CPU-only experimental predicate. No fast-math or ISA-specific intrinsics.
#include <cmath>
#include <cstddef>
#include <cstdint>
extern "C" void host_predicate(const float* values, const float* critical,
                                std::uint8_t* mask, std::size_t rows,
                                std::size_t columns) {
    for (std::size_t row = 0; row < rows; ++row) {
        const float limit = critical[row];
        for (std::size_t column = 0; column < columns; ++column) {
            const float value = std::fabs(values[row * columns + column]);
            mask[row * columns + column] = std::isfinite(value) && value >= limit;
        }
    }
}
