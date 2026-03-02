# LDPC Decoder Benchmark: Conventional vs RMAS1 vs RMAS2

This repository provides a compact C-based LDPC simulation to compare four decoders under the same channel and matrix setup:

- **Conventional** normalized/offset Min-Sum (flooding schedule)
- **RMAS1** residual-magnitude-aided grouped schedule (group-priority layered updates + damping)
- **RMAS2** residual-magnitude-aided check-priority schedule (per-check residual sorting + damping)
- **AS** adaptive-scaling Min-Sum decoder (iteration-dependent normalization factor)

The draft algorithm note is documented in:

- `docs/rmas_ldpc_decoder_algorithm_draft.pdf`

## Repository layout

- `include/ldpc.h`: public API, decoder enums, parameter types.
- `src/ldpc.c`: matrix generation, channel model, and decoding kernels (Conventional/RMAS1/RMAS2).
- `src/main.c`: benchmark runner and CSV writer.
- `scripts/plot_results.py`: result visualization utility (multi-panel BER/FER/iterations + summary).
- `results/performance.csv`: numerical benchmark output for the default/custom matrix profile.
- `results/matrix_performance.csv`: cross-matrix benchmark output for IEEE-like profiles.
- `results/performance.svg`: always-generated performance report figure.
- `results/matrix_performance.svg`: BER figure comparing IEEE-like matrix profiles.
- `results/performance.png`: optional PNG export (when matplotlib is available).

## Build

```bash
make
```

## Run

```bash
./ldpc_eval
```

Optional arguments:

```bash
./ldpc_eval <n> <m> <dv> <frames>
```

Defaults:

- `n=512`
- `m=64`
- `dv=6`
- `frames=300`

The program sweeps SNR from **4.5 dB to 7.0 dB** and writes BER/FER/average iterations to `results/performance.csv` and `results/matrix_performance.csv`, then automatically generates `results/performance.svg` and `results/matrix_performance.svg` (and optionally `results/performance.png`).

> Plot generation requires Python 3. PNG additionally requires `matplotlib`.

## Decoder notes

Shared tunables (`ldpc_decoder_params_t`):

- `max_iters`
- `alpha` (normalized factor)
- `beta` (offset)

RMAS-specific knobs:

- `group_size`
- `damping`

RMAS1 schedules check-node groups by group residual magnitude each iteration. RMAS2 sorts all checks by residual magnitude and refreshes neighboring variable beliefs immediately after each selected check update. AS linearly ramps the normalization factor from `alpha` to `alpha_final` across iterations.

## License

MIT (`LICENSE`).
