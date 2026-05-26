# Usage example (extdata)

End-to-end example on the bundled SAIGE `extdata` test set
(`/media/leelabsg-storage0/UKBB_WORK/SAIGE_cpp/extdata`): 1000 samples, 128868
markers, quantitative trait. Step-1 fits the null model; Step-2 runs
single-variant association chained from that null. Verified to run end-to-end.

## 0. Build (see ../README.md for the full env)
```bash
ENV=/data/home/seokhojeong/.local/share/mamba/envs/saige-build
export CONDA_PREFIX=$ENV PATH=$ENV/bin:$PATH
export LD_LIBRARY_PATH=$ENV/lib:$ENV/lib/R/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$ENV/lib CPATH=$ENV/include PKG_CONFIG_PATH=$ENV/lib/pkgconfig CXX=g++
( cd ../step1_saige-null  && make clean && make )                                   # -> saige-null
( cd ../step2_saige-step2 && CXXFLAGS="-std=c++17 -O3 -DARMA_USE_SUPERLU=1 -DARMA_64BIT_WORD=1 -MMD -MP" make clean && make )  # -> saige-step2
mkdir -p output
```

## 1. Step-1 — null model (`step1_quant_extdata.yaml`)
```bash
../step1_saige-null/saige-null -c step1_quant_extdata.yaml
```
Dense-GRM fit, quantitative trait, covariates x1+x2. Converges in ~5 iterations
(tau ≈ [0.20, 0.49]; small run-to-run variation from the 30-vector Monte-Carlo
trace is expected). Produces in `output/`:
- `step1_quant_null/`              — null model directory (`nullmodel.json` + `*.arma`); this is Step-2's `modelFile`
- `step1_quant_vr.varianceRatio.txt` — variance ratio (single value); this is Step-2's `varianceRatioFile`

## 2. Step-2 — single-variant association (`step2_singleassoc_extdata.yaml`)
```bash
../step2_saige-step2/saige-step2 step2_singleassoc_extdata.yaml
```
Chains from the Step-1 null (`modelFile` = `output/step1_quant_null`). Tests all
128868 markers from the same PLINK set and writes per-marker results
(`CHR POS MarkerID … BETA SE Tstat var p.value N`) to
`output/step2_singleassoc.txt`. Peak RSS ~50 MB.

## Notes
- Paths in the YAMLs are absolute (this machine). Edit `out_prefix` / `outputFile`
  to redirect outputs.
- This example uses the **single-variant** path (dense-fit null → single variance
  ratio). For **SAIGE-GENE+ region/gene tests**, fit Step-1 with
  `use_sparse_grm_to_fit: true` (categorical variance ratio) and run Step-2 with a
  `groupFile` + `annotationList` + `maxMAFList` (see the `config_*` templates in
  `../step2_saige-step2/`). NOTE: GENE+ region tests have a large peak RSS
  (~N×N; ~240 GB at ~165k samples) and need a high-memory node.
