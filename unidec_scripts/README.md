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

## `unidec_batch_apply.py` — apply fixed conditions (production step)

Once the optimiser has found reproducible conditions, use this to **apply** them
to data — no searching. It reads a `shared_conditions.csv` (or the full
`best_config_shared.json`) from an optimiser run and deconvolves every `.raw`
file with those exact parameters, reporting the test/IS ratio per file. This is
what you run on future single-sample-per-concentration data.

```bash
# from the tuned-conditions CSV (geometry filled from the CLI defaults):
python unidec_batch_apply.py --config results_.../shared_conditions.csv \
    --data /path/to/raw --out applied

# exact reproduction from the full JSON (geometry included, no CLI needed):
python unidec_batch_apply.py --config results_.../best_config_shared.json \
    --data /path/to/raw --out applied
```

`shared_conditions.csv` stores the *tuned* parameters only; the geometry
(centroids, m/z window, mass bounds, charge range) is not swept and is taken
from the same CLI defaults as the optimiser — pass the same
`--centroidA/--centroidB/--minmz/--maxmz/--zrange` if you changed them. The JSON
config stores the complete parameter set, so it needs no geometry flags.

Outputs: `batch_per_file.csv` (ratio + peak quality per file),
`batch_by_concentration.csv` (mean/SD/%CV per concentration, if filenames carry
tokens), `batch_summary.txt` (conditions used + calibration slope/intercept/R²),
and `ratio_vs_concentration.png`. It reuses the deconvolution/ratio/quality code
from `unidec_ratio_optimizer.py`, so results match the optimiser exactly.

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

### By-concentration (`--by-concentration`) — the calibration workflow
For a concentration series with 2–3 replicates each (test light chain 23412 Da
vs internal-standard light chain 23658 Da; ratio = test/IS). Files are grouped
by a concentration token in the filename, and the program reports **two answers
side by side**:

- **Per-concentration:** for each concentration, the parameter set minimising
  that concentration's cross-replicate %CV.
- **Shared:** one parameter set that minimises the *aggregated* within-
  concentration %CV across **all** concentrations at once.

The shared result is the more trustworthy basis for future single-sample runs:
with only 2–3 replicates, a per-concentration CV minimum is a fragile target
prone to overfitting, whereas the shared set is fit against every replicate of
every concentration. The per-concentration answer is kept as a **drift
diagnostic** — if the winning parameters barely change across concentrations,
the shared set is clearly right.

Both use an overfitting-resistant selection: a **peak-quality floor**
(`--quality-floor`, default 0.5) excludes degenerate sets where A and B merged,
and among sets within a **CV margin** of the best (`--cv-margin`, default 0.25 =
25%) it picks the *highest-quality* one rather than the razor-thin minimum.

The token is scientific notation with `p` as the decimal point, e.g.
`1p00e-4` → 1.0e-4 µg/mL. Files sharing a token are replicates; files with no
token are reported and excluded (override with `--conc-pattern`). No ratio band
is applied — each concentration has its own true ratio.

```bash
python unidec_ratio_optimizer.py --data /path/to/raw_folder --out results --by-concentration
```

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
| Peak FWHM            | `mzsig`     | 0.3, 0.4, 0.5 |
| Charge Smooth Width  | `zzsig`     | 1.0, 5.0, 10.0 |
| Point Smooth Width   | `psig`      | 0.0, 1.0, 2.0 |
| Beta                 | `beta`      | 0.0, 50.0, 100.0 |
| Gaussian Smoothing   | `smooth`    | 1.0, 2.0, 5.0 |
| Subtract Curved width| `subbuff`   | 50, 100, 150 |

Full grid = 3×3×3×3×3×3 = **729** parameter sets. `--quick` = 2×2×2×2×2×1 =
**32** sets.

Fixed defaults (from `Params`): `massbins=0.1` Da, `mzbins=1.0` Th,
`subtype=1` (curved), m/z window `700–2000`. Background subtraction defaults to
**Subtract Curved** with a realistic width, so curved is actually tested.

## Performance

The raw-spectrum load (`open_file`) is parameter-independent, so each `.raw` is
opened **once** and the loaded engine is reused across the whole grid; only
`process_data()`/`run_unidec()` — which depend on the swept parameters — re-run
per grid point. This removes redundant file I/O but does **not** make the sweep
"instant": the deconvolution itself is the swept quantity and runs once per grid
point (729 in the full grid, per file). Pass `--fresh-engine` to open a new
engine every call (slower; only for ruling out engine state carry-over).

> If you want a genuinely cache-once, integrate-in-memory sweep, that is the
> `mz_ratio_direct.py` model — it integrates the RAW spectrum directly, so the
> only per-parameter cost is cheap integration. It skips deconvolution entirely.

## Outputs

The output directory gets a `_YYYYmmdd_HHMMSS` timestamp appended automatically
(e.g. `results_20260702_182910`) so repeat runs never overwrite each other.
Pass `--no-timestamp` to disable that.

Every mode except `--per-file` writes a `cv_percent` column (a single file has
no cross-replicate CV, so `--per-file` has none). In default quality mode the
`cv_percent` is descriptive — reported, not the optimisation target.

**Default quality mode** → `quality_sweep_results.csv` (all sets ranked by mean
quality, with a `cv_percent` column), `best_config.json`, `best_per_file.csv`
(per-file ratio + quality for the winner), `quality_vs_ratio.png`.

**`--per-file`** → `per_file_optimized.csv` (each file's own best params + ratio
+ quality), `best_config_per_file.json`, `per_file_sweep_full.csv` (full audit
of every trial), `per_file_ratios.png`.

**`--by-concentration`** → `comparison.csv` (headline: per-concentration vs
shared CV/ratio, and whether params differ), `per_concentration_best.csv`,
`shared_conditions.csv` (the single shared set + calibration slope/intercept/R²),
`best_config_by_concentration.json`, `best_config_shared.json`,
`per_file_ratios.csv` (per-replicate ratios for both selections),
`by_concentration_sweep_full.csv` (full audit), `ratio_vs_concentration.png`
(per-concentration vs shared calibration curves with linear-fit R²).

**`--cv-mode`** → `sweep_results.csv`, `best_config.json`, `best_per_file.csv`,
`ratio_vs_cv.png`.

⚠ Always visually inspect the winning mass spectrum — confirm A and B are
cleanly separated. In `--per-file` mode, ratios come from different parameter
sets, so cross-check against the default (shared-conditions) result before
trusting them.
