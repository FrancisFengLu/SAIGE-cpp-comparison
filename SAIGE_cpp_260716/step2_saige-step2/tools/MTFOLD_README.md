# mtFoldQuantProj acceptance scripts

The A/B harness for the `mtFoldQuantProj` config key (see
`saige_mt.hpp`, `MTContext::foldQuant`). Full write-up:
`SAIGE-work/optimization/torchgwas2/MT_FOLD_QUANT.md`.

| file | what it does |
|---|---|
| `mtfold_gen_cfg.py` | writes one step-2 config; `q:<n>` / `b:<n>` pick n quantitative / binary traits, `SGS=1` adds `outputFormat: sgs` |
| `mtfold_run_one.sh` | one run (`BIN=`, `RUNROOT=` overridable) |
| `mtfold_bench.sh` | the P=8 / P=32 wall-clock A/B on g1m |
| `mtfold_cmp_out.py` | text-output field comparison of two run directories |
| `mtfold_cmp_sgs.py` | **fp64** comparison via the `.sgs` binary output -- the text prints p to 7 digits, so -log10 p there cannot resolve better than 4.34e-7 |
| `mtfold_sgsread.py` | minimal python reader for a `.sgs` trait file |
| `mtfold_armaio.py` | reads an `arma_binary` matrix |
| `mtfold_residual.py` | fit / analytic / naive residual of the p x p fold map against the stored `XVX_inv_XV` |
| `mtfold_stat_ab.py` | fp64 numpy A/B of the fitted vs naive fold on simulated genotypes |

The model and data paths are the ones on the dev box
(`/opt/saige/logs/tg2_step2/...`, `/opt/saige/data/mid`); edit the constants at
the top of `mtfold_gen_cfg.py` / `mtfold_residual.py` elsewhere.
