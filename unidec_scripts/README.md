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

## `unidec_ratio_optimizer.py` — optimise UniDec processing conditions

Runs UniDec and sweeps its processing parameters. The **objective is peak-fit
quality** — how cleanly A and B deconvolve into well-shaped, well-separated
peaks:

```
quality = 0.5 * separation + 0.25 * (fitR2_A + fitR2_B)
```

- `separation` = `1 - valley/min(apexA, apexB)` — trough depth between the two
  fitted peaks. Penalises over-smoothing that merges A and B (a false win).
- `fitR2_A/B` = goodness of a local Gaussian+baseline fit to each peak.

The metric never looks at the ratio value, so tuning the conditions **cannot
bias the quantitation**.

### Default — one shared set of conditions (maximise mean quality)
Finds the single parameter set with the highest mean quality across all files.
Use this to pick one best set of UniDec processing conditions.

```bash
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results --quick
```

### Per-file — optimise each file independently (`--per-file`)
Optimises parameters separately for every `.raw` file, each maximising that
file's own peak-fit quality.

```bash
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results --per-file
```

### By-concentration (`--by-concentration`) — CV within each concentration
For a concentration series with 2–3 replicates each. Files are grouped by a
concentration token in the filename, and **within each concentration** the
program finds the parameter set that minimises the cross-replicate %CV of the
A/B ratio. You get **one winning parameter set per concentration**.

The token is scientific notation with `p` as the decimal point, e.g.
`1p00e-4` → 1.0e-4 µg/mL, `2e-4` → 2e-4 µg/mL. Files sharing a token are
treated as replicates of the same concentration; files with no token are
reported and excluded. Override the pattern with `--conc-pattern` if your
naming differs.

```bash
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results --by-concentration
```

No ratio band is applied here — each concentration has its own (unknown) true
A/B ratio, so the only objective is minimising that concentration's replicate
CV.

### Legacy %CV objective (`--cv-mode`)
One shared set minimising cross-replicate %CV across *all* files (ignores
concentration). Retained but off by default.

```bash
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results --cv-mode
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

**Default quality mode** → `quality_sweep_results.csv` (all sets ranked by mean
quality), `best_config.json`, `best_per_file.csv` (per-file ratio + quality for
the winner), `quality_vs_ratio.png`.

**`--per-file`** → `per_file_optimized.csv` (each file's own best params + ratio
+ quality), `best_config_per_file.json`, `per_file_sweep_full.csv` (full audit
of every trial), `per_file_ratios.png`.

**`--by-concentration`** → `by_concentration_best.csv` (one winning set per
concentration), `best_config_by_concentration.json`,
`by_concentration_per_file.csv` (per-replicate ratios under each winner),
`by_concentration_sweep_full.csv` (full audit), `ratio_vs_concentration.png`.

**`--cv-mode`** → `sweep_results.csv`, `best_config.json`, `best_per_file.csv`,
`ratio_vs_cv.png`.

⚠ Always visually inspect the winning mass spectrum — confirm A and B are
cleanly separated. In `--per-file` mode, ratios come from different parameter
sets, so cross-check against the default (shared-conditions) result before
trusting them.
