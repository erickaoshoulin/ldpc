# LDPC Conventional vs Proposed Decoder (C Implementation)

This repository now contains a **modular and parameterized C implementation** to evaluate two LDPC decoding flows inspired by the paper draft in `draft.pdf`:

- **Conventional**: flooding-style normalized/offset min-sum.
- **Proposed**: grouped/layered variable-node-centric style update with damping.

> Note: The original draft PDF is hardware-focused. This implementation intentionally focuses on algorithmic behavior and software-side BER/FER evaluation only.

## Repository layout

- `include/ldpc.h`: public APIs and parameter structs.
- `src/ldpc.c`: LDPC matrix generation, AWGN model, decoding kernels.
- `src/main.c`: evaluation driver that compares conventional vs proposed decoder.
- `results/performance.csv`: generated performance table.

## Build

```bash
make
```

## Run evaluation

```bash
./ldpc_eval
```

Optional arguments:

```bash
./ldpc_eval <n> <m> <dv> <frames>
```

Default values:

- `n=512` code length
- `m=64` parity checks
- `dv=6` variable-node degree
- `frames=300` Monte-Carlo frames per SNR point

The executable sweeps SNR from **4.5 dB to 7.0 dB** and writes:

- BER (bit error rate)
- FER (frame error rate)
- average iterations

to `results/performance.csv`.

## Algorithm parameterization

Both decoders share these tunables in `ldpc_decoder_params_t`:

- `max_iters`
- `alpha` (normalized min-sum factor)
- `beta` (offset min-sum correction)

Proposed decoder adds:

- `group_size` (check-node group processing size)
- `damping` (message damping factor)

## Open-source

This project is now distributed under the MIT License (`LICENSE`).
