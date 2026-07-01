# UniDec ratio scripts

Two scripts for quantifying the abundance ratio of two co-eluting ~23.5 kDa
protein species (A ≈ 23,412 Da, B ≈ 23,658 Da, ~250 Da apart) from replicate
native-MS `.raw` files. Both drive the external
[UniDec](https://github.com/michaelmarty/UniDec) package (`pip install unidec`)
via its Python engine (`import unidec.engine`). They do **not** depend on the
`brainpy` package in this repository; they simply live here on this branch.

> **Platform note:** UniDec's Thermo `.raw` reader is effectively Windows-only.
> On other platforms, convert to `.mzML` (e.g. with `msconvert`).

## Requirements

```bash
pip install unidec numpy matplotlib scipy
```
`scipy` is required for the Gaussian peak-fit paths (`--integration-mode gaussian`
and the per-file quality metric).

## `mz_ratio_direct.py` — model-free ratio (no deconvolution)

Integrates the raw m/z signal directly, per charge state, in a safe zone
(default z = 15–22), and reports `Σ area_A / Σ area_B`.

```bash
python mz_ratio_direct.py --data /path/to/raw_folder --out results_direct
```

## `unidec_ratio_optimizer.py` — deconvolution + parameter sweep

Runs UniDec and sweeps deconvolution parameters. Two modes:

### Global (default) — one shared parameter set
Finds the single parameter set that minimises the cross-replicate **%CV** of
the A/B ratio.

```bash
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results --quick
```

### Per-file — optimise each file independently (`--per-file`)
Optimises parameters separately for every `.raw` file. Because cross-replicate
%CV is undefined once each file uses different parameters, files are scored by
a **ratio-independent peak-fit quality** metric:

```
quality = 0.5 * separation + 0.25 * (fitR2_A + fitR2_B)
```

- `separation` = `1 - valley/min(apexA, apexB)` — trough depth between the two
  fitted peaks. This penalises over-smoothing that merges A and B (a false win).
- `fitR2_A/B` = goodness of a local Gaussian+baseline fit to each peak.

The metric never looks at the ratio value, so per-file tuning cannot bias the
quantitation.

```bash
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results --per-file
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results --per-file --quick
```

## What the sweep tunes (v3)

The sweep focuses on the deconvolution smoothing/focusing priors the UniDec app
exposes, over higher ranges than before:

| App label            | Config attr | Swept values (full) |
|----------------------|-------------|---------------------|
| Peak FWHM            | `mzsig`     | 0.4, 0.6, 1.0 |
| Charge Smooth Width  | `zzsig`     | 1.0, 2.0, 4.0 |
| Point Smooth Width   | `psig`      | 0.0, 1.0, 2.0 |
| Beta                 | `beta`      | 0.0, 50.0, 100.0 |
| Subtract Curved width| `subbuff`   | 50, 100, 150 |

Fixed defaults (from `Params`): `massbins=0.1` Da, `mzbins=1.0` Th,
`smooth=2.0`, `subtype=1` (curved). Background subtraction now defaults to
**Subtract Curved** with a realistic width, so curved is actually tested.

## Outputs

**Global mode** → `sweep_results.csv`, `best_config.json`, `best_per_file.csv`,
`ratio_vs_cv.png`.

**Per-file mode** → `per_file_optimized.csv` (each file's own best params +
ratio + quality), `best_config_per_file.json`, `per_file_sweep_full.csv` (full
audit of every trial), `per_file_ratios.png`.

⚠ Always visually inspect the winning mass spectrum — confirm A and B are
cleanly separated. In per-file mode, ratios come from different parameter sets,
so cross-check against the global-mode result before trusting them.
