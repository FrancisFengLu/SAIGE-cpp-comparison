# LOCO null-model format (step 1 → step 2 contract)

Status: spec, 2026-08-10. Both sides must match this exactly.

## Background

LOCO does not re-fit the variance component. `theta`/tau is estimated once on
the full-genome GRM; then, holding tau fixed, only the fixed effects are
re-solved per chromosome (R: `Get_Coef_LOCO`, `SAIGE_fitGLMM_fast.R:185-231`).

So only the quantities derived from the fixed-effect solve vary by chromosome.
`mu` changes, therefore `V = mu(1-mu)` changes, therefore every product
containing `V` changes.

## Directory layout

```
<modelFile>/
  nullmodel.json          # gains loco fields, see below
  X.arma                  # chromosome-INVARIANT, written once
  y.arma                  # chromosome-INVARIANT, written once
  mu.arma                 # full-genome fit (unchanged, always written)
  res.arma
  V.arma
  offset.arma
  XV.arma
  XVX.arma
  XVX_inv.arma
  XVX_inv_XV.arma
  XXVX_inv.arma
  S_a.arma
  chr1/                   # only when LOCO ran
    mu.arma res.arma V.arma offset.arma
    XV.arma XVX.arma XVX_inv.arma XVX_inv_XV.arma XXVX_inv.arma S_a.arma
  chr2/
  ...
  chr22/
```

**Per-chromosome set** (10 files): `mu`, `res`, `V`, `offset`, `XV`, `XVX`,
`XVX_inv`, `XVX_inv_XV`, `XXVX_inv`, `S_a`.

**Never duplicated into `chr<j>/`**: `X`, `y`. At N=165k with 20 covariates,
duplicating `X` alone would add ~570 MB for no reason.

**Top-level copies are always written**, LOCO or not. They hold the
full-genome fit, which step 2 falls back to for non-autosomal chromosomes
(mirroring R's `readInGLMM.R:107-113`).

Rationale for directories over a stacked cube: step 2 is run one chromosome
per process, so each process reads only its own `chr<j>/` (~84 MB at UKB
scale) instead of all 22 (~1.8 GB). It also needs no new reader code — the
existing `.arma` loader is pointed at a subdirectory. R's `isLowMemLOCO` mode
uses the same per-chromosome-file structure.

## nullmodel.json additions

```json
{
  "loco": true,
  "loco_chroms": [1, 2, ..., 22]
}
```

- `loco` — `true` only if LOCO actually ran and the `chr<j>/` directories were
  written. Never set optimistically from the config flag.
- `loco_chroms` — the autosomes actually present, ascending. Absent or `[]`
  when `loco` is false. Step 2 must use this list rather than assuming 1..22,
  since the input may not cover every autosome.

## Step 2 behaviour

New config keys: `LOCO` (bool, default false until step 1 ships LOCO) and
`chrom` (string).

Load order: read top-level files, then if LOCO is active overwrite the
per-chromosome set from `chr<chrom>/`.

Guards, matching R:

| condition | behaviour | R reference |
|---|---|---|
| `LOCO: true`, model has `loco: false` | error, exit non-zero | `readInGLMM.R:79` |
| `LOCO: true`, `chrom` empty | error, exit non-zero | `readInGLMM.R:82` |
| `chrom` not in `loco_chroms` (e.g. X, Y, MT) | silent fallback to the full-genome fit | `readInGLMM.R:107-113` |
| `LOCO: false` | ignore `chr<j>/` entirely | — |

Marker filtering: when LOCO is on, single-variant tests must restrict markers
to `chrom` (R: `SAIGE_SPATest_Marker.R:44-50`). Region tests need no extra
logic — they inherit the swapped model through the same loader.

## Compatibility

A non-LOCO run must produce byte-identical output to the current build. The
`chr<j>/` directories are purely additive; a step 2 build that predates this
spec ignores them and still works.
