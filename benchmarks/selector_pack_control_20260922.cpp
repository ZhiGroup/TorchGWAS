// Experiment-only selected-result packing. No predicate or scientific change.
// Input flat indices are owned, sorted, nonnegative and in range; rows reuses
// that buffer only after all reads at the current index. No allocations/GIL.
#include <cstddef>
#include <cstdint>

extern "C" void selector_pack_control(
    std::int64_t* rows, std::int64_t* columns,
    const float* values, const float* beta, const float* df,
    float* out_values, float* out_beta, float* out_df,
    std::size_t count, std::int64_t width, std::int64_t start,
    std::ptrdiff_t df_stride) {
    for (std::size_t i = 0; i < count; ++i) {
        const std::int64_t flat = rows[i];
        const std::int64_t row = flat / width;
        columns[i] = flat % width;
        out_values[i] = values[flat];
        if (beta) out_beta[i] = beta[flat];
        out_df[i] = df[row * df_stride];
        rows[i] = row + start;
    }
}
