# LDPC Decoder Benchmark: Conventional vs RMAS1 vs RMAS2

This repository provides a compact C-based LDPC simulation to compare:

- **Conventional** normalized/offset Min-Sum (flooding schedule)
- **RMAS1** grouped residual-aided message passing (layered-in-groups + damping)
- **RMAS2** fine-grained residual-aided message passing (per-check refresh + damping)

The algorithmic draft is documented in:

- `docs/ldpc_decoder_draft.pdf`

## Repository layout

- `include/ldpc.h`: public API, decoder enums, parameter types.
- `src/ldpc.c`: matrix generation, channel model, and decoding kernels.
- `src/main.c`: benchmark runner and CSV writer.
- `scripts/plot_results.py`: result visualization utility.
- `results/performance.csv`: numerical benchmark output.
- `results/performance.svg`: always-generated BER/FER/iterations performance curves.
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

The program sweeps SNR from **4.5 dB to 7.0 dB** and writes:

- BER
- FER
- average iteration count

into `results/performance.csv`, then generates a summary figure in `results/performance.svg` (and optionally `results/performance.png`).

> Plot generation requires Python 3. PNG additionally requires `matplotlib`.

## Decoder notes

Shared tunables (`ldpc_decoder_params_t`):

- `max_iters`
- `alpha` (normalized factor)
- `beta` (offset)

RMAS-specific knobs:

- `group_size`
- `damping`

RMAS2 uses the same struct and applies a denser per-check APP refresh schedule intended to emulate a more aggressive residual propagation.

## License

MIT (`LICENSE`).
