"""Achieved GEMM rate from the DEVICE SPECIFICATION and the shape.

This replaces a lookup table. `pipeline_model.MEASURED_GEMM_TFLOPS_BY_WIDTH`
holds four numbers -- 8.35, 8.71, 20.57, 20.28 TFLOPS at design widths 36, 60,
156, 540 -- and has three problems that no amount of re-measuring fixes:

  * it CLAMPS outside the measured span, so at K=1024 (width 1052) it returns
    the width-540 rate. That single clamp makes the model 1.81x wrong on the
    K ladder (32.79 s predicted against 18.09 s measured), which is the largest
    single error left anywhere in the calculator;
  * it is one machine. An A100 or a 2080 Ti needs the whole table re-taken
    before the model says anything at all about them;
  * its own note says the values were taken at host load 155 with a competitor
    benchmark running and are LOWER BOUNDS awaiting a quiet re-take.

None of that is necessary, because an achieved GEMM rate is not a mystery
constant. It is the device's arithmetic throughput, capped by what the chosen
tiling can feed, reduced by how badly the shape divides into tiles and waves.
Every one of those is either a published architectural parameter or arithmetic
on the shape.

THE DERIVATION

1. PEAK. `SMs * FP32 lanes per SM * 2 (fused multiply-add) * clock`. SM count
   and clock come from the device; lanes per SM is fixed by the compute
   capability (see `FP32_LANES_PER_SM`). For an H100 SXM: 132 * 128 * 2 *
   1.98e9 = 66.9 TFLOP/s, which is the figure NVIDIA publishes.

2. PER-TILE MEMORY CEILING. A TM x TW output tile reduces over the whole K
   dimension, loading TM*K and K*TW elements to produce 2*TM*TW*K flops, so its
   arithmetic intensity is

       2 * TM * TW / (4 * (TM + TW))   FLOP per byte

   independent of K. A 32x32 tile is 8 FLOP/byte and on 3.35 TB/s can never
   exceed 26.8 TFLOP/s no matter how many SMs are idle waiting for it; a
   128x128 tile is 32 FLOP/byte and is compute-bound instead. This is the term
   that stops the model over-predicting the wide-and-flat shapes, and note that
   the GLOBAL intensity does not capture it: the whole product at chunk 4,096
   and width 156 is 74.8 FLOP/byte, far above the machine balance of 20.0, so
   a naive roofline calls it compute-bound and is wrong.

3. TILE QUANTIZATION. The output is (chunk x width) covered by whole tiles, so
   the useful fraction is `chunk/(ceil(chunk/TM)*TM) * width/(ceil(width/TW)*TW)`.
   Design width 156 against a 128-wide tile wastes 39% of the second tile.

4. WAVE QUANTIZATION. Tiles are distributed over SMs; a partial last wave
   leaves SMs idle, so the useful fraction is `tiles/(ceil(tiles/SMs)*SMs)`.
   64 tiles on 132 SMs is 48% of the machine.

5. TILE CHOICE. cuBLAS picks a kernel per shape, and picks close to the best
   available, so the model takes the MAXIMUM over a catalog of standard tile
   shapes and split-K factors rather than assuming one. Split-K matters
   precisely for the shapes here: a 512-variant chunk makes too few tiles to
   fill 132 SMs, and splitting the reduction is how the hardware gets filled.

WHAT THIS IS AND IS NOT

This is an upper BOUND, and bounds are what rooflines give. Validated against a
measured `torch.mm` sweep at design width 156, samples 35,365, chunk 256 to
16,384 (`benchmarks/direct_gemm_chunk_curve.py`):

    chunk      256    512   1024   2048   4096   8192  16384
    measured 24.72  33.66  37.34  35.25  38.24  37.23  38.55
    bound    39.5   42.2   52.7   52.7   52.7   52.7   54.2
    ratio     1.60   1.25   1.41   1.50   1.38   1.42   1.41

Never violated (every ratio above 1.0, as a bound must be) and flat to within
1.27x across a 64x span of chunk. The realized fraction is the one number this
module does not derive: real SGEMM does not reach its roofline, and how close
it gets is a property of the device and its driver.

`realized_fraction` is therefore an explicit argument and defaults to None,
which returns the bound and says so. It is a MACHINE RATE to be measured once
per device with a single large square GEMM -- a shape this model never
predicts, so calibrating on it and validating on the sweep above are
independent. It is deliberately NOT set from the sweep it would then be
validated against; that would be fitting.
"""
from __future__ import annotations

from math import ceil

# FP32 lanes per SM, fixed by the architecture. These are published
# specifications, not measurements: a Hopper or Ada SM issues 128 FP32
# fused-multiply-adds per clock, a Volta/Turing/A100 SM issues 64.
#
# Keyed by compute capability major.minor. Checked against the vendor's own
# headline numbers: A100 (108 SM, 64 lanes, 1.41 GHz) -> 19.5 TFLOP/s; H100 SXM
# (132, 128, 1.98) -> 66.9; RTX 2080 Ti (68, 64, 1.545) -> 13.45. All three
# match the published figures, which is the check that the table is right.
FP32_LANES_PER_SM = {
    (7, 0): 64,    # Volta, V100
    (7, 2): 64,
    (7, 5): 64,    # Turing, 2080 Ti
    (8, 0): 64,    # Ampere GA100, A100
    (8, 6): 128,   # Ampere GA10x
    (8, 7): 128,
    (8, 9): 128,   # Ada
    (9, 0): 128,   # Hopper, H100
    (10, 0): 128,  # Blackwell
    (12, 0): 128,
}

# Standard cuBLAS SGEMM output tiles. Not exhaustive and not claimed to be --
# the point is that the library chooses among shapes of roughly this family,
# so maximising over them approximates its choice. Adding a shape can only
# raise the bound, never lower it.
TILE_CATALOG = ((32, 32), (64, 32), (32, 64), (64, 64), (128, 64),
                (64, 128), (128, 128), (256, 128), (128, 256))

# Split-K factors. A skinny GEMM makes too few output tiles to fill the SMs,
# and splitting the reduction is the standard remedy. Capped at 16 because
# beyond that the reduction over partials stops being negligible and this
# module does not model that cost -- so the catalog stops where the assumption
# it rests on stops.
SPLIT_K_FACTORS = (1, 2, 3, 4, 6, 8, 12, 16)

BYTES_PER_FLOAT32 = 4


class DeviceSpec:
    """The handful of published numbers a GEMM bound needs.

    Constructed from `torch.cuda.get_device_properties` where possible, but
    kept as a plain object so the model can be evaluated for a device that is
    not present -- planning a run on a card you do not have in front of you is
    most of what a calculator is for.
    """

    def __init__(self, name: str, multiprocessors: int, capability: tuple,
                 clock_hz: float, memory_bandwidth_bytes_per_second: float,
                 fp32_lanes_per_sm: int | None = None):
        if multiprocessors <= 0 or clock_hz <= 0:
            raise ValueError("multiprocessors and clock must be positive")
        if memory_bandwidth_bytes_per_second <= 0:
            raise ValueError("memory bandwidth must be positive")
        self.name = name
        self.multiprocessors = int(multiprocessors)
        self.capability = tuple(capability)
        self.clock_hz = float(clock_hz)
        self.memory_bandwidth = float(memory_bandwidth_bytes_per_second)
        if fp32_lanes_per_sm is None:
            if self.capability not in FP32_LANES_PER_SM:
                raise ValueError(
                    f"unknown compute capability {self.capability}; pass "
                    f"fp32_lanes_per_sm explicitly rather than guessing -- the "
                    f"value differs by a factor of two between architectures")
            fp32_lanes_per_sm = FP32_LANES_PER_SM[self.capability]
        self.fp32_lanes_per_sm = int(fp32_lanes_per_sm)

    @property
    def peak_flops_per_second(self) -> float:
        """SMs x lanes x 2 (FMA) x clock. The published headline number."""
        return (self.multiprocessors * self.fp32_lanes_per_sm
                * 2.0 * self.clock_hz)

    @property
    def machine_balance_flops_per_byte(self) -> float:
        """Peak divided by bandwidth: the intensity a kernel must beat."""
        return self.peak_flops_per_second / self.memory_bandwidth

    def __repr__(self) -> str:
        return (f"DeviceSpec({self.name!r}, {self.multiprocessors} SM, "
                f"{self.peak_flops_per_second / 1e12:.1f} TFLOP/s peak)")


# Published specifications for the cards this project runs on or plans for.
# Clock is the boost clock and bandwidth the published figure; both are vendor
# numbers, so a disagreement with measurement is informative rather than
# something to tune away.
KNOWN_DEVICES = {
    "H100": DeviceSpec("NVIDIA H100 80GB HBM3", 132, (9, 0), 1.98e9, 3.35e12),
    "A100": DeviceSpec("NVIDIA A100 80GB", 108, (8, 0), 1.41e9, 2.039e12),
    "2080Ti": DeviceSpec("NVIDIA GeForce RTX 2080 Ti", 68, (7, 5),
                         1.545e9, 616e9),
}


def tile_intensity_flops_per_byte(tile_rows: int, tile_columns: int) -> float:
    """Arithmetic intensity of one output tile, independent of the reduction.

    A TM x TW tile loads TM*K + K*TW elements and produces 2*TM*TW*K flops, so
    K cancels. This is the quantity a global roofline misses: the product as a
    whole can look compute-bound while every tile that computes it is starved.
    """
    if tile_rows <= 0 or tile_columns <= 0:
        raise ValueError("tile dimensions must be positive")
    return (2.0 * tile_rows * tile_columns
            / (BYTES_PER_FLOAT32 * (tile_rows + tile_columns)))


def gemm_bound_flops_per_second(device: DeviceSpec, rows: int,
                                columns: int) -> dict:
    """Upper bound on achieved FLOP/s for a (rows x K) @ (K x columns) product.

    K does not appear: every term above is independent of the reduction length,
    which is why this is a function of the two shape dimensions alone.

    Returns the bound and the tiling that achieves it, because the tiling is
    the explanation -- "width 156 wastes 39% of a 128-wide tile" is actionable
    where a bare number is not.
    """
    if rows <= 0 or columns <= 0:
        raise ValueError("rows and columns must be positive")
    peak = device.peak_flops_per_second
    best = 0.0
    best_tiling = None
    for tile_rows, tile_columns in TILE_CATALOG:
        intensity = tile_intensity_flops_per_byte(tile_rows, tile_columns)
        ceiling = min(peak, device.memory_bandwidth * intensity)
        covered_rows = ceil(rows / tile_rows) * tile_rows
        covered_columns = ceil(columns / tile_columns) * tile_columns
        tile_efficiency = (rows / covered_rows) * (columns / covered_columns)
        base_tiles = ceil(rows / tile_rows) * ceil(columns / tile_columns)
        for split in SPLIT_K_FACTORS:
            tiles = base_tiles * split
            waves = ceil(tiles / device.multiprocessors)
            wave_efficiency = tiles / (waves * device.multiprocessors)
            achieved = ceiling * tile_efficiency * wave_efficiency
            if achieved > best:
                best = achieved
                best_tiling = {
                    'tile': (tile_rows, tile_columns), 'split_k': split,
                    'tile_efficiency': tile_efficiency,
                    'wave_efficiency': wave_efficiency,
                    'memory_bound': device.memory_bandwidth * intensity < peak,
                }
    return {'bound_flops_per_second': best, 'peak_flops_per_second': peak,
            'fraction_of_peak': best / peak if peak else 0.0, **best_tiling}


def gemm_flops_per_second(device: DeviceSpec, rows: int, columns: int,
                          realized_fraction: float | None = None) -> float:
    """Predicted achieved FLOP/s, or the bound when the fraction is unknown.

    `realized_fraction` is the device's realized share of its own roofline --
    one measured number, taken once with a large square GEMM. Passing None
    returns the BOUND, which is honest but optimistic by about 1.4x on the
    H100; a caller feeding this into a time model should know which it got,
    so there is no default value pretending otherwise.
    """
    bound = gemm_bound_flops_per_second(device, rows, columns)
    if realized_fraction is None:
        return bound['bound_flops_per_second']
    if not 0.0 < realized_fraction <= 1.0:
        raise ValueError("realized_fraction must be in (0, 1]")
    return bound['bound_flops_per_second'] * realized_fraction


def device_spec_from_torch(index: int = 0) -> DeviceSpec:
    """Read the live device's specification.

    Bandwidth is computed from the memory clock and bus width when the build
    exposes them (`2 *` for double data rate), and falls back to the published
    figure for a known device. It raises rather than guessing for an unknown
    one: a fabricated bandwidth would silently move the memory ceiling that
    this whole module exists to apply.
    """
    import torch

    properties = torch.cuda.get_device_properties(index)
    capability = (properties.major, properties.minor)
    clock = getattr(properties, "clock_rate", None)
    bandwidth = None
    memory_clock = getattr(properties, "memory_clock_rate", None)
    bus_width = getattr(properties, "memory_bus_width", None)
    if memory_clock and bus_width:
        bandwidth = 2.0 * (memory_clock * 1e3) * (bus_width / 8.0)
    if bandwidth is None or not clock:
        for known in KNOWN_DEVICES.values():
            if known.name in properties.name or properties.name in known.name:
                bandwidth = bandwidth or known.memory_bandwidth
                clock = clock or known.clock_hz / 1e3
                break
    if bandwidth is None or not clock:
        raise ValueError(
            f"cannot determine clock/bandwidth for {properties.name!r}; add it "
            f"to KNOWN_DEVICES from the published specification rather than "
            f"letting the model invent one")
    return DeviceSpec(properties.name, properties.multi_processor_count,
                      capability, clock * 1e3, bandwidth)
