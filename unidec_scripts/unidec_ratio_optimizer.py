#!/usr/bin/env python3
"""
unidec_ratio_optimizer.py  –  v3 (UniDec 8.2.1, confirmed API)
===============================================================
Sweeps UniDec deconvolution parameters across replicate .raw files to quantify
the ratio between two protein species (A ~23410 Da, B ~23660 Da, ~250 Da apart).

Optimisation objective
-----------------------
The objective is UniDec PROCESSING QUALITY: how cleanly the two species (A, B)
deconvolve into well-shaped, well-separated peaks. This is scored by a metric
that is deliberately INDEPENDENT of the ratio VALUE, so tuning the processing
conditions cannot bias the quantitation -- see peak_quality().

Three modes
-----------
  QUALITY, shared (default):  find ONE parameter set that maximises the MEAN
                      peak-fit quality across all files. Use this to pick a
                      single best set of UniDec processing conditions.

  QUALITY, per-file (--per-file):  optimise parameters SEPARATELY for each
                      .raw file, each maximising that file's own peak-fit
                      quality.

  %CV (--cv-mode):    the legacy objective -- find ONE parameter set that
                      minimises the cross-replicate %CV of the A/B ratio.
                      Retained but off by default.

What changed vs v2  (driven by real-data feedback)
--------------------------------------------------
  * subtype default is now 1 (SUBTRACT CURVED), with subbuff swept around the
    app's real curved-subtraction width (~100). v2 only paired curved
    subtraction with subbuff 0-20, so curved never had a fair chance even
    though it is visually the best background model for this data.
  * Gaussian smoothing (`smooth`) now defaults >1.0 and is applied by default,
    matching how the UniDec app is typically run on this data.
  * m/z linearisation bin (`mzbins`) default 1.0 Th (was 0.0 = raw spacing).
  * mass sampling (`massbins`) default 0.1 Da (was 1.0) for finer A/B mass axis.
  * The sweep now optimises the deconvolution SMOOTHING/FOCUSING priors the
    app calls Beta, Charge Smooth Width, and Point Smooth Width:
        beta   (softmax charge focusing / noise suppression)
        zzsig  (charge smooth width)
        psig   (point smooth width)
    over HIGHER ranges than v2 (v2 landed on values that were too low).

Confirmed against:  unidec 8.2.1, Python 3.11, Windows
Real API entry:     import unidec.engine as engine; u = engine.UniDec()

USAGE
-----
  # DEFAULT: one shared set of UniDec conditions, maximise mean quality:
  python unidec_ratio_optimizer.py --data "D:\\path\\to\\data" --out results

  # PER-FILE: optimise each replicate independently by peak-fit quality:
  python unidec_ratio_optimizer.py --data "D:\\..." --out results --per-file

  # LEGACY %CV objective:
  python unidec_ratio_optimizer.py --data "D:\\..." --out results --cv-mode

  # quick first-pass grid (16 combos), any mode:
  python unidec_ratio_optimizer.py --data "D:\\..." --out results --quick

OUTPUTS (default quality mode)
------------------------------
  quality_sweep_results.csv  all parameter sets ranked by mean quality
  best_config.json           winning parameters
  best_per_file.csv          per-file ratio + quality for the winner
  quality_vs_ratio.png       scatter of quality vs ratio

OUTPUTS (--per-file)
--------------------
  per_file_optimized.csv     each file's own best params + ratio + quality
  best_config_per_file.json  winning params per file
  per_file_sweep_full.csv    full audit: every (file, param set) tried
  per_file_ratios.png        per-file ratio bar chart

OUTPUTS (--cv-mode)
-------------------
  sweep_results.csv          all parameter sets ranked by %CV
  best_config.json / best_per_file.csv / ratio_vs_cv.png
"""

from __future__ import annotations
import argparse, csv, itertools, json, os, re, sys, traceback, warnings
from dataclasses import dataclass, asdict, field
import numpy as np

warnings.filterwarnings("ignore")   # suppress unidec's mpld3 / psims noise

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:
    _HAVE_MPL = False

# ---------------------------------------------------------------------------
# UniDec engine  (confirmed: unidec.engine.UniDec, v8.2.1)
# ---------------------------------------------------------------------------

def make_engine():
    """Return a fresh UniDec engine instance.

    UniDec's engine parses sys.argv on construction (it looks for -f/-o/-c
    flags for its own CLI). We temporarily replace sys.argv with just the
    script name while the engine is constructed, then restore it.
    """
    try:
        import unidec.engine as _eng
    except ImportError as e:
        raise SystemExit(
            "Cannot import unidec.engine. Make sure UniDec 8.x is installed:\n"
            "    pip install unidec\n"
            f"(error: {e})"
        )
    _saved_argv = sys.argv[:]
    sys.argv = sys.argv[:1]   # hide our flags from UniDec's arg parser
    try:
        eng = _eng.UniDec()
    finally:
        sys.argv = _saved_argv
    return eng


# ---------------------------------------------------------------------------
# Parameter container  (attribute names verified against dir(u.config))
# ---------------------------------------------------------------------------

@dataclass
class Params:
    # --- m/z input window (covers z=10..40 for ~23.5 kDa with headroom) ---
    minmz: float = 600.0
    maxmz: float = 2500.0

    # --- mass output grid: bracket both species tightly ---
    masslb: float = 22500.0
    massub: float = 24500.0
    massbins: float = 0.1       # Da/bin (was 1.0); 0.1 gives a finer A/B axis

    # --- charge range ---
    startz: int = 8
    endz: int = 45

    # --- native charge clipping (wide = effectively off) ---
    nativezlb: int = -1000
    nativezub: int = 1000

    # --- deconvolution kernel ---
    mzsig: float = 1.0          # FWHM of peak shape in m/z (key tuning param)
    psfun: int = 0              # 0=Gaussian, 1=Lorentzian, 2=split-Gaussian
    zzsig: float = 2.0          # CHARGE SMOOTH WIDTH (was 1.0; swept higher)
    psig: float = 1.0           # POINT SMOOTH WIDTH   (new; swept)
    beta: float = 50.0          # softmax charge focusing / noise suppression (new; swept)

    # --- convergence ---
    numit: int = 50

    # --- preprocessing ---
    smooth: float = 2.0         # Gaussian smoothing (was 0.0; app value usually >1.0)
    mzbins: float = 1.0         # m/z linearisation bin size Th (was 0.0 = raw spacing)
    linflag: int = 2            # interpolation mode for linearisation

    # --- background subtraction ---
    # subtype: 0=none, 1=CURVED, 2=linear. Curved with subbuff~100 matches the
    # app default and is the best visual background model for this data.
    subtype: int = 1            # was 2 (linear); now 1 (curved)
    subbuff: float = 100.0      # curved-subtraction width (was 0.0)

    # --- noise / normalisation ---
    intthresh: float = 0.0
    datanorm: int = 1
    adductmass: float = 1.007276

    # --- post-decon integration (not passed to UniDec) ---
    #   'window'   : wide fixed-window trapezoidal sum
    #   'gaussian' : narrow local Gaussian + linear-baseline fit per peak
    integration_mode: str = "window"

    # --- window-mode settings ---
    baseline_mode: str = "flat"     # 'none' | 'flat' | 'linear'
    bl_left_lo:  float = 22800.0
    bl_left_hi:  float = 23200.0
    bl_right_lo: float = 24000.0
    bl_right_hi: float = 24400.0
    int_lo_A:  float = 23200.0
    int_split: float = 23535.0
    int_hi_B:  float = 23900.0

    # --- gaussian-mode + quality-metric settings ---
    centroidA_expected: float = 23412.0
    centroidB_expected: float = 23658.0
    peak_half_window:   float = 20.0   # Da; narrow window around each centroid

    def key_str(self):
        base = (f"bins={self.massbins} mzsig={self.mzsig} zzsig={self.zzsig} "
                f"psig={self.psig} beta={self.beta} "
                f"sub={self.subtype}/{self.subbuff} smooth={self.smooth}")
        if self.integration_mode == "gaussian":
            return base + f" int=gaussian(hw={self.peak_half_window})"
        return base + f" int=window bl={self.baseline_mode}"


# ---------------------------------------------------------------------------
# Apply params to a UniDec config object
# ---------------------------------------------------------------------------

_WARNED_MISSING = set()


def _set_cfg(cfg, name, val):
    """Set a config attribute, warning once if the installed UniDec lacks it."""
    if hasattr(cfg, name):
        setattr(cfg, name, val)
    elif name not in _WARNED_MISSING:
        _WARNED_MISSING.add(name)
        print(f"  [WARN] u.config has no attribute '{name}' in this UniDec "
              f"version; skipping it.", file=sys.stderr)


def apply_params(cfg, p: Params):
    """Push every relevant Params field onto u.config."""
    cfg.minmz       = p.minmz
    cfg.maxmz       = p.maxmz
    cfg.masslb      = p.masslb
    cfg.massub      = p.massub
    cfg.massbins    = p.massbins
    cfg.startz      = p.startz
    cfg.endz        = p.endz
    cfg.nativezlb   = p.nativezlb
    cfg.nativezub   = p.nativezub
    cfg.mzsig       = p.mzsig
    cfg.psfun       = p.psfun
    cfg.zzsig       = p.zzsig
    # newer priors: guard in case of an older config object
    _set_cfg(cfg, "psig", p.psig)
    _set_cfg(cfg, "beta", p.beta)
    cfg.numit       = p.numit
    cfg.smooth      = p.smooth
    cfg.mzbins      = p.mzbins
    cfg.linflag     = p.linflag
    cfg.subtype     = p.subtype
    cfg.subbuff     = p.subbuff
    cfg.intthresh   = p.intthresh
    cfg.datanorm    = p.datanorm
    cfg.adductmass  = p.adductmass


# ---------------------------------------------------------------------------
# Deconvolution + integration
# ---------------------------------------------------------------------------

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


def baseline_correct(mass, inten, p: Params):
    """Estimate and subtract a local baseline from the mass spectrum."""
    if p.baseline_mode == "none":
        return inten, 0.0

    left_sel  = (mass >= p.bl_left_lo)  & (mass <= p.bl_left_hi)
    right_sel = (mass >= p.bl_right_lo) & (mass <= p.bl_right_hi)
    bl_left   = np.median(inten[left_sel])  if left_sel.any()  else 0.0
    bl_right  = np.median(inten[right_sel]) if right_sel.any() else 0.0

    if p.baseline_mode == "flat":
        bl = 0.5 * (bl_left + bl_right)
        corrected = np.maximum(inten - bl, 0.0)
        return corrected, bl

    bl_arr = np.interp(mass,
                       [0.5*(p.bl_left_lo+p.bl_left_hi),
                        0.5*(p.bl_right_lo+p.bl_right_hi)],
                       [bl_left, bl_right])
    corrected = np.maximum(inten - bl_arr, 0.0)
    return corrected, float(np.mean(bl_arr))


def integrate_region(mass, inten, lo, hi):
    sel = (mass >= lo) & (mass <= hi)
    if sel.sum() < 2:
        return 0.0
    return float(_trapz(inten[sel], mass[sel]))


# --- Gaussian local peak-fit -----------------------------------------------

def _gaussian(x, amp, center, sigma, base, slope):
    return amp * np.exp(-0.5 * ((x - center) / sigma) ** 2) + base + slope * (x - center)


def _fit_local_peak_full(mass, inten, expected_center, half_window, verbose=False):
    """
    Fit amp*Gaussian + linear baseline to a narrow window around expected_center.

    Returns (area, popt, success, r2, apex) where:
      area = analytic Gaussian area (amp*sigma*sqrt(2pi)), baseline EXCLUDED
      r2   = goodness of fit over the window (1 = perfect)
      apex = max raw intensity in the window (used for valley/separation)
    """
    from scipy.optimize import curve_fit

    sel = (mass >= expected_center - half_window) & (mass <= expected_center + half_window)
    x = mass[sel]
    y = inten[sel]
    if x.size < 8:
        return 0.0, None, False, float("nan"), float("nan")

    edge_n = max(2, x.size // 6)
    base0 = float(np.median(np.concatenate([y[:edge_n], y[-edge_n:]])))
    amp0 = max(float(y.max() - base0), 1e-6)
    sigma0 = half_window / 4.0

    p0 = [amp0, expected_center, sigma0, base0, 0.0]
    bounds = (
        [0.0, expected_center - half_window * 0.6, 0.5, -np.inf, -np.inf],
        [np.inf, expected_center + half_window * 0.6, half_window, np.inf, np.inf],
    )
    try:
        popt, _ = curve_fit(_gaussian, x, y, p0=p0, bounds=bounds, maxfev=5000)
    except Exception as exc:
        if verbose:
            print(f"    [fit_local_peak] failed at {expected_center}: {exc}")
        return 0.0, None, False, float("nan"), float("nan")

    amp, center, sigma, base, slope = popt
    area = float(amp * sigma * np.sqrt(2 * np.pi))
    yhat = _gaussian(x, *popt)
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    apex = float(y.max())
    return area, popt, True, r2, apex


def fit_local_peak(mass, inten, expected_center, half_window, verbose=False):
    """Back-compat wrapper: returns (area, popt, success)."""
    area, popt, ok, _r2, _apex = _fit_local_peak_full(
        mass, inten, expected_center, half_window, verbose)
    return area, popt, ok


def compute_ratio(mass, inten_raw, p: Params, verbose: bool = False):
    """Compute ratio A/B using whichever integration_mode is selected."""
    if p.integration_mode == "gaussian":
        aA, _, okA = fit_local_peak(mass, inten_raw, p.centroidA_expected,
                                    p.peak_half_window, verbose)
        aB, _, okB = fit_local_peak(mass, inten_raw, p.centroidB_expected,
                                    p.peak_half_window, verbose)
        if okA and okB and aB > 0:
            return aA / aB, aA, aB, np.nan
        if verbose:
            print("    [compute_ratio] gaussian fit failed, falling back to window")

    inten_bc, bl_val = baseline_correct(mass, inten_raw, p)
    aA = integrate_region(mass, inten_bc, p.int_lo_A,  p.int_split)
    aB = integrate_region(mass, inten_bc, p.int_split, p.int_hi_B)
    if aB <= 0:
        return np.nan, aA, aB, bl_val
    return aA / aB, aA, aB, bl_val


def peak_quality(mass, inten, p: Params, verbose: bool = False):
    """
    Ratio-INDEPENDENT measure of how cleanly A and B deconvolved.

    Combines, into a single score in [0, 1]:
      * fit quality : local Gaussian+baseline R^2 for A and B (are the peaks
                      real, well-shaped Gaussians?)
      * separation  : 1 - valley/min(apexA, apexB), the depth of the trough
                      between the two fitted centres (are A and B RESOLVED,
                      or smeared into one blob?).

    The separation term is what stops over-smoothing (huge zzsig/beta merging
    A and B) from scoring as a false win. This score never looks at the ratio
    VALUE, so choosing params by it cannot bias the quantitation.

    Returns dict(scoreQ, sep, r2A, r2B) or None if either peak fails to fit.
    """
    aA, poptA, okA, r2A, apexA = _fit_local_peak_full(
        mass, inten, p.centroidA_expected, p.peak_half_window, verbose)
    aB, poptB, okB, r2B, apexB = _fit_local_peak_full(
        mass, inten, p.centroidB_expected, p.peak_half_window, verbose)
    if not (okA and okB):
        return None

    cA, cB = poptA[1], poptB[1]
    lo, hi = sorted((cA, cB))
    seg = (mass >= lo) & (mass <= hi)
    if seg.sum() >= 3 and min(apexA, apexB) > 0:
        valley = float(np.min(inten[seg]))
        sep = float(np.clip(1.0 - valley / min(apexA, apexB), 0.0, 1.0))
    else:
        sep = 0.0

    r2A_c = float(np.clip(r2A, 0.0, 1.0)) if np.isfinite(r2A) else 0.0
    r2B_c = float(np.clip(r2B, 0.0, 1.0)) if np.isfinite(r2B) else 0.0
    score = 0.5 * sep + 0.25 * (r2A_c + r2B_c)
    return dict(scoreQ=score, sep=sep, r2A=r2A_c, r2B=r2B_c)


def deconvolve_one(path: str, p: Params, verbose: bool = False,
                   workdir: str | None = None):
    """
    Open one .raw file, apply params, run UniDec, return
    (mass, intensity, ratio, aA, aB, baseline, r_squared).
    """
    import io, contextlib, shutil, tempfile

    own_tmpdir = workdir is None
    if own_tmpdir:
        tmpdir = tempfile.mkdtemp(prefix="unidec_opt_")
    else:
        tmpdir = workdir
        os.makedirs(tmpdir, exist_ok=True)

    try:
        fname    = os.path.basename(path)
        tmp_path = os.path.join(tmpdir, fname)
        if not os.path.exists(tmp_path):
            shutil.copy2(path, tmp_path)

        u = make_engine()
        u.silent = not verbose

        def _run():
            u.open_file(tmp_path)
            apply_params(u.config, p)
            u.process_data()
            u.run_unidec()

        if not verbose:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                _run()
        else:
            _run()

        massdat = np.asarray(u.data.massdat)
        if massdat.ndim != 2 or massdat.shape[0] < 5:
            raise RuntimeError("UniDec returned an empty mass distribution")

        mass  = massdat[:, 0]
        inten = massdat[:, 1]
        r2    = float(u.config.error)   # confirmed: u.config.error in v8.2.1

        ratio, aA, aB, bl = compute_ratio(mass, inten, p, verbose)
        return mass, inten, ratio, aA, aB, bl, r2

    finally:
        if own_tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Parameter grid
# ---------------------------------------------------------------------------

def build_grid(args) -> list[Params]:
    """
    Construct the sweep grid.

    The sweep now focuses on the three deconvolution priors the app exposes as
    Beta, Charge Smooth Width, and Point Smooth Width, over ranges HIGHER than
    v2 (which landed on values that were too low), plus peak width and the
    curved-subtraction width. Background subtraction is fixed to CURVED
    (subtype=1) via Params default, with subbuff swept around the app's ~100
    default. massbins / mzbins / smooth are taken from the Params defaults
    (0.1 Da / 1.0 Th / 2.0) rather than swept.
    """
    if args.quick:
        grid = dict(
            mzsig   = [0.6, 1.0],
            zzsig   = [1.0, 2.0],     # charge smooth width
            psig    = [0.0, 1.0],     # point smooth width
            beta    = [0.0, 50.0],    # softmax charge focusing
            subbuff = [100.0],        # curved-subtraction width
        )
    else:
        grid = dict(
            mzsig   = [0.4, 0.6, 1.0],
            zzsig   = [1.0, 2.0, 4.0],        # charge smooth width (higher)
            psig    = [0.0, 1.0, 2.0],        # point smooth width
            beta    = [0.0, 50.0, 100.0],     # softmax charge focusing (higher)
            subbuff = [50.0, 100.0, 150.0],   # curved-subtraction width
        )

    keys  = list(grid.keys())
    combos = list(itertools.product(*[grid[k] for k in keys]))
    out = []
    for combo in combos:
        p = Params(
            minmz      = args.minmz,
            maxmz      = args.maxmz,
            masslb     = args.centroidA - args.mass_pad,
            massub     = args.centroidB + args.mass_pad,
            startz     = args.zrange[0],
            endz       = args.zrange[1],
            int_lo_A   = args.centroidA - 210,
            int_split  = 0.5 * (args.centroidA + args.centroidB),
            int_hi_B   = args.centroidB + 240,
            integration_mode    = args.integration_mode,
            centroidA_expected  = args.centroidA,
            centroidB_expected  = args.centroidB,
        )
        for k, v in zip(keys, combo):
            setattr(p, k, v)
        out.append(p)
    return out


# ---------------------------------------------------------------------------
# GLOBAL mode: one shared parameter set, minimise cross-replicate %CV
# ---------------------------------------------------------------------------

def evaluate(p: Params, files: list[str], args) -> dict | None:
    ratios, r2s, per_file = [], [], []

    for path in files:
        fname = os.path.basename(path)
        try:
            _, _, ratio, aA, aB, bl, r2 = deconvolve_one(path, p, args.verbose)
            ratios.append(ratio)
            r2s.append(r2)
            per_file.append((fname, ratio, aA, aB, bl, r2))
        except Exception as exc:
            per_file.append((fname, np.nan, np.nan, np.nan, np.nan, np.nan))
            if args.verbose:
                print(f"    [WARN] {fname}: {exc}")

    ratios = np.array(ratios, float)
    good   = ratios[np.isfinite(ratios) & (ratios > 0)]
    if good.size < max(2, int(0.5 * len(files))):
        return None

    mean_r = float(np.mean(good))
    std_r  = float(np.std(good, ddof=1))
    cv     = 100.0 * std_r / mean_r
    mean_r2 = float(np.nanmean(r2s))

    return dict(mean_ratio=mean_r, std_ratio=std_r, cv_percent=cv,
                n_good=int(good.size), mean_r2=mean_r2, per_file=per_file)


def optimise(files: list[str], args):
    grid = build_grid(args)
    n = len(grid)
    print(f"[GLOBAL] Sweeping {n} parameter sets × {len(files)} files "
          f"= {n*len(files)} deconvolutions")
    print(f"Ratio band accepted: [{args.ratio_band[0]}, {args.ratio_band[1]}]")
    print()

    rows = []
    for i, p in enumerate(grid, 1):
        res = evaluate(p, files, args)
        if res is None:
            continue
        rows.append((p, res))
        tag = "*IN-BAND*" if args.ratio_band[0] <= res["mean_ratio"] <= args.ratio_band[1] else ""
        if args.verbose or i % max(1, n//20) == 0:
            print(f"  [{i:>4}/{n}] ratio={res['mean_ratio']:.3f}  "
                  f"CV={res['cv_percent']:.2f}%  R²={res['mean_r2']:.3f}  "
                  f"{tag}  | {p.key_str()}")

    if not rows:
        raise SystemExit("All parameter sets failed. Check file paths and UniDec install.")

    lo, hi = args.ratio_band
    in_band = [(p, r) for p, r in rows if lo <= r["mean_ratio"] <= hi]

    if in_band:
        in_band.sort(key=lambda pr: pr[1]["cv_percent"])
        winner = in_band[0]
        basis  = f"lowest %CV with mean ratio in [{lo}, {hi}]"
    else:
        def score(pr):
            r = pr[1]
            drift = abs(np.log(max(r["mean_ratio"], 1e-6) / args.target_ratio))
            return r["cv_percent"] + args.lam * 100.0 * drift
        rows.sort(key=score)
        winner = rows[0]
        basis  = (f"no set in [{lo},{hi}]; minimised CV + "
                  f"{args.lam}×|log(ratio/{args.target_ratio})|")

    return rows, winner, basis


# ---------------------------------------------------------------------------
# PER-FILE mode: optimise each file independently by peak-fit quality
# ---------------------------------------------------------------------------

def evaluate_per_file(p: Params, path: str, args) -> dict | None:
    """Deconvolve ONE file with ONE param set; return ratio + quality, or None."""
    try:
        mass, inten, ratio, aA, aB, bl, r2 = deconvolve_one(path, p, args.verbose)
    except Exception as exc:
        if args.verbose:
            print(f"    [WARN] {os.path.basename(path)}: {exc}")
        return None
    if not (np.isfinite(ratio) and ratio > 0):
        return None
    q = peak_quality(mass, inten, p, args.verbose)
    if q is None:
        return None
    return dict(ratio=ratio, aA=aA, aB=aB, bl=bl, r2=r2,
                quality=q["scoreQ"], sep=q["sep"], r2A=q["r2A"], r2B=q["r2B"])


def optimise_per_file(files: list[str], args):
    """
    For each file, sweep the grid and keep the param set with the highest
    peak-fit quality. Returns (winners, full_records) where
      winners      = [(fname, Params, result_dict), ...]
      full_records = [(fname, Params, result_dict), ...]  (every valid trial)
    """
    grid = build_grid(args)
    n = len(grid)
    print(f"[PER-FILE] Sweeping {n} parameter sets on EACH of {len(files)} files "
          f"= {n*len(files)} deconvolutions")
    print("Objective: maximise per-file peak-fit quality (ratio-independent).")
    print()

    winners, full_records = [], []
    for path in files:
        fname = os.path.basename(path)
        best = None
        n_valid = 0
        for p in grid:
            res = evaluate_per_file(p, path, args)
            if res is None:
                continue
            n_valid += 1
            full_records.append((fname, p, res))
            if best is None or res["quality"] > best[1]["quality"]:
                best = (p, res)
        if best is None:
            print(f"  [WARN] {fname:<40} no valid parameter set")
            continue
        winners.append((fname, best[0], best[1]))
        print(f"  {fname:<40} ratio={best[1]['ratio']:.4f}  "
              f"quality={best[1]['quality']:.3f} "
              f"(sep={best[1]['sep']:.2f}, R²A={best[1]['r2A']:.2f}, "
              f"R²B={best[1]['r2B']:.2f})  n_valid={n_valid}\n"
              f"      -> {best[0].key_str()}")

    if not winners:
        raise SystemExit("No file produced a valid parameter set. "
                         "Check file paths and UniDec install.")
    return winners, full_records


# ---------------------------------------------------------------------------
# QUALITY mode (default): one shared parameter set, maximise mean peak quality
# ---------------------------------------------------------------------------

def evaluate_quality(p: Params, files: list[str], args) -> dict | None:
    """Aggregate peak-fit quality of ONE param set across all files."""
    per_file, quals = [], []
    for path in files:
        fname = os.path.basename(path)
        res = evaluate_per_file(p, path, args)
        per_file.append((fname, res))
        if res is not None:
            quals.append(res["quality"])

    if len(quals) < max(2, int(0.5 * len(files))):
        return None

    quals  = np.array(quals, float)
    ratios = np.array([r["ratio"] for _, r in per_file if r is not None], float)
    mean_ratio = float(np.mean(ratios)) if ratios.size else float("nan")
    std_ratio  = float(np.std(ratios, ddof=1)) if ratios.size > 1 else 0.0
    desc_cv    = 100.0 * std_ratio / mean_ratio if mean_ratio else float("nan")

    return dict(mean_quality=float(np.mean(quals)),
                min_quality=float(np.min(quals)),
                n_good=int(quals.size),
                mean_ratio=mean_ratio, std_ratio=std_ratio,
                desc_cv=desc_cv, per_file=per_file)


def optimise_global_quality(files: list[str], args):
    grid = build_grid(args)
    n = len(grid)
    print(f"[QUALITY] Sweeping {n} parameter sets × {len(files)} files "
          f"= {n*len(files)} deconvolutions")
    print("Objective: maximise MEAN peak-fit quality (one shared parameter set).")
    print()

    rows = []
    for i, p in enumerate(grid, 1):
        res = evaluate_quality(p, files, args)
        if res is None:
            continue
        rows.append((p, res))
        if args.verbose or i % max(1, n // 20) == 0:
            print(f"  [{i:>4}/{n}] Q={res['mean_quality']:.3f} "
                  f"(worst {res['min_quality']:.3f})  ratio={res['mean_ratio']:.3f}"
                  f"  | {p.key_str()}")

    if not rows:
        raise SystemExit("All parameter sets failed. Check file paths and UniDec install.")

    # Rank by mean quality; tie-break by the worst single file's quality.
    rows.sort(key=lambda pr: (pr[1]["mean_quality"], pr[1]["min_quality"]),
              reverse=True)
    winner = rows[0]
    basis  = "highest mean peak-fit quality across files (one shared set)"
    return rows, winner, basis


# ---------------------------------------------------------------------------
# BY-CONCENTRATION mode: group replicates by a filename token, minimise the
# cross-replicate %CV WITHIN each concentration group (one winner per group)
# ---------------------------------------------------------------------------

# Matches a scientific-notation concentration token such as "1p00e-4", "2e-4",
# "10p0e-5"  (mantissa, optional 'p' decimal, 'e', optional sign, exponent).
_CONC_RE = re.compile(r"\d+(?:p\d+)?e[-+]?\d+", re.IGNORECASE)


def parse_concentration(fname: str, pattern: re.Pattern | None = None):
    """Return the (lower-cased) concentration token in a filename, or None."""
    m = (pattern or _CONC_RE).search(fname)
    return m.group(0).lower() if m else None


def conc_to_float(token: str):
    """Convert a token like '1p00e-4' to the float 1e-4 (None if unparseable)."""
    try:
        return float(token.replace("p", "."))
    except (ValueError, AttributeError):
        return None


def group_by_concentration(files, pattern=None):
    """Group file paths by concentration token. Returns (groups, ungrouped)."""
    groups: dict[str, list[str]] = {}
    ungrouped: list[str] = []
    for path in files:
        key = parse_concentration(os.path.basename(path), pattern)
        if key is None:
            ungrouped.append(path)
        else:
            groups.setdefault(key, []).append(path)
    return groups, ungrouped


def optimise_by_concentration(files, args):
    """
    For each concentration group (>=2 replicates), sweep the grid and keep the
    parameter set with the lowest cross-replicate %CV of the A/B ratio.

    Returns winners = [(conc, group_files, Params, res, rows), ...] where
    `res` is the evaluate() summary for the winning set and `rows` is the full
    ranked sweep for that group (for the audit CSV).
    """
    pattern = re.compile(args.conc_pattern, re.IGNORECASE) if args.conc_pattern else _CONC_RE
    groups, ungrouped = group_by_concentration(files, pattern)

    if ungrouped:
        print("[WARN] no concentration token found in these files (excluded):")
        for p in ungrouped:
            print(f"        {os.path.basename(p)}")
        print()

    if not groups:
        raise SystemExit(
            "No concentration tokens found. Expected a filename part like "
            "'1p00e-4'. Override the pattern with --conc-pattern if needed.")

    grid = build_grid(args)
    print(f"[BY-CONCENTRATION] {len(groups)} concentration group(s), "
          f"{len(grid)} parameter sets each:")
    for conc in sorted(groups, key=lambda c: (conc_to_float(c) or 0.0, c)):
        print(f"        {conc:<12} {len(groups[conc])} replicate(s): "
              + ", ".join(os.path.basename(p) for p in groups[conc]))
    print()

    winners = []
    for conc in sorted(groups, key=lambda c: (conc_to_float(c) or 0.0, c)):
        gfiles = groups[conc]
        if len(gfiles) < 2:
            print(f"  [WARN] {conc}: only {len(gfiles)} replicate; "
                  f"cannot compute a cross-replicate CV. Skipping.")
            continue

        rows = []
        for p in grid:
            res = evaluate(p, gfiles, args)   # cross-replicate CV within group
            if res is not None:
                rows.append((p, res))
        if not rows:
            print(f"  [WARN] {conc}: all parameter sets failed.")
            continue

        # Objective within a concentration: lowest %CV. No ratio band here,
        # since each concentration has its own (unknown) true A/B ratio.
        rows.sort(key=lambda pr: pr[1]["cv_percent"])
        wp, wr = rows[0]
        winners.append((conc, gfiles, wp, wr, rows))
        print(f"  {conc:<12} best CV={wr['cv_percent']:.2f}%  "
              f"ratio={wr['mean_ratio']:.3f}  (n={wr['n_good']})  | {wp.key_str()}")

    if not winners:
        raise SystemExit("No concentration group could be optimised "
                         "(need >=2 replicates per group).")
    return winners


# ---------------------------------------------------------------------------
# Output writing — GLOBAL (%CV mode)
# ---------------------------------------------------------------------------

def write_outputs(rows, winner, basis, args):
    os.makedirs(args.out, exist_ok=True)
    wp, wr = winner

    sweep_path = os.path.join(args.out, "sweep_results.csv")
    with open(sweep_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank","cv_percent","mean_ratio","std_ratio","n_good",
                    "mean_r2","massbins","mzsig","zzsig","psig","beta",
                    "subtype","subbuff","smooth","integration_mode",
                    "baseline_mode","peak_half_window"])
        for rank, (p, r) in enumerate(
                sorted(rows, key=lambda pr: pr[1]["cv_percent"]), 1):
            w.writerow([rank, f"{r['cv_percent']:.3f}", f"{r['mean_ratio']:.4f}",
                        f"{r['std_ratio']:.4f}", r["n_good"],
                        f"{r['mean_r2']:.4f}", p.massbins, p.mzsig, p.zzsig,
                        p.psig, p.beta, p.subtype, p.subbuff, p.smooth,
                        p.integration_mode, p.baseline_mode, p.peak_half_window])

    best_json = os.path.join(args.out, "best_config.json")
    with open(best_json, "w") as fh:
        json.dump({
            "selection_basis": basis,
            "summary": {k: v for k, v in wr.items() if k != "per_file"},
            "params": asdict(wp)
        }, fh, indent=2)

    perfile_path = os.path.join(args.out, "best_per_file.csv")
    with open(perfile_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file","ratio_A_over_B","area_A","area_B","baseline","r2"])
        for name, ratio, aA, aB, bl, r2 in wr["per_file"]:
            fmt = lambda x: f"{x:.5g}" if np.isfinite(x) else "NA"
            w.writerow([name, fmt(ratio), fmt(aA), fmt(aB), fmt(bl), fmt(r2)])

    plot_path = None
    if _HAVE_MPL:
        try:
            cvs   = [r["cv_percent"]   for _, r in rows]
            ratios = [r["mean_ratio"]  for _, r in rows]
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.scatter(ratios, cvs, s=18, alpha=0.55, label="all param sets")
            ax.axvspan(args.ratio_band[0], args.ratio_band[1],
                       alpha=0.12, color="green", label="accepted ratio band")
            ax.scatter([wr["mean_ratio"]], [wr["cv_percent"]], s=160,
                       marker="*", color="crimson", zorder=6, label="winner")
            ax.axvline(args.target_ratio, ls="--", color="grey",
                       lw=1, label=f"target ratio {args.target_ratio}")
            ax.set_xlabel("Mean A/B ratio")
            ax.set_ylabel("Cross-replicate %CV")
            ax.set_title("UniDec parameter sweep: %CV vs ratio")
            ax.legend(fontsize=8)
            fig.tight_layout()
            plot_path = os.path.join(args.out, "ratio_vs_cv.png")
            fig.savefig(plot_path, dpi=140)
            plt.close(fig)
        except Exception:
            pass

    return sweep_path, best_json, perfile_path, plot_path


# ---------------------------------------------------------------------------
# Output writing — PER-FILE
# ---------------------------------------------------------------------------

_PARAM_COLS = ["massbins","mzsig","zzsig","psig","beta","subtype","subbuff",
               "smooth","integration_mode","baseline_mode","peak_half_window"]


def _param_row(p: Params):
    return [p.massbins, p.mzsig, p.zzsig, p.psig, p.beta, p.subtype,
            p.subbuff, p.smooth, p.integration_mode, p.baseline_mode,
            p.peak_half_window]


def write_outputs_per_file(winners, full_records, args):
    os.makedirs(args.out, exist_ok=True)

    # Per-file winners
    opt_path = os.path.join(args.out, "per_file_optimized.csv")
    with open(opt_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file","ratio_A_over_B","area_A","area_B","quality",
                    "separation","fitR2_A","fitR2_B","baseline"] + _PARAM_COLS)
        for fname, p, r in winners:
            fmt = lambda x: f"{x:.5g}" if np.isfinite(x) else "NA"
            w.writerow([fname, fmt(r["ratio"]), fmt(r["aA"]), fmt(r["aB"]),
                        f"{r['quality']:.4f}", f"{r['sep']:.4f}",
                        f"{r['r2A']:.4f}", f"{r['r2B']:.4f}", fmt(r["bl"])]
                       + _param_row(p))

    # Winning config per file
    best_json = os.path.join(args.out, "best_config_per_file.json")
    with open(best_json, "w") as fh:
        json.dump({
            fname: {
                "summary": {k: r[k] for k in
                            ("ratio","quality","sep","r2A","r2B","aA","aB")},
                "params": asdict(p),
            }
            for fname, p, r in winners
        }, fh, indent=2)

    # Full audit trail: every (file, param set) tried
    full_path = os.path.join(args.out, "per_file_sweep_full.csv")
    with open(full_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file","quality","ratio_A_over_B","separation",
                    "fitR2_A","fitR2_B"] + _PARAM_COLS)
        for fname, p, r in sorted(full_records,
                                  key=lambda t: (t[0], -t[2]["quality"])):
            w.writerow([fname, f"{r['quality']:.4f}", f"{r['ratio']:.5g}",
                        f"{r['sep']:.4f}", f"{r['r2A']:.4f}", f"{r['r2B']:.4f}"]
                       + _param_row(p))

    # Plot: per-file ratios (with descriptive mean line)
    plot_path = None
    ratios = np.array([r["ratio"] for _, _, r in winners], float)
    if _HAVE_MPL and ratios.size > 1:
        try:
            mean_r = float(np.mean(ratios))
            fnames = [f.replace(" Scan[1].raw", "") for f, _, _ in winners]
            fig, ax = plt.subplots(figsize=(max(7, len(fnames)*0.7), 4.5))
            ax.bar(range(len(ratios)), ratios, color="#185FA5", alpha=0.85)
            ax.axhline(mean_r, color="crimson", ls="--", lw=1.5,
                       label=f"mean {mean_r:.3f}")
            ax.set_xticks(range(len(fnames)))
            ax.set_xticklabels(fnames, rotation=40, ha="right", fontsize=8)
            ax.set_ylabel("A/B ratio (per-file optimised)")
            ax.set_title("Per-file optimised A/B ratio")
            ax.legend(fontsize=8)
            fig.tight_layout()
            plot_path = os.path.join(args.out, "per_file_ratios.png")
            fig.savefig(plot_path, dpi=140)
            plt.close(fig)
        except Exception:
            pass

    return opt_path, best_json, full_path, plot_path


# ---------------------------------------------------------------------------
# Output writing — QUALITY (default)
# ---------------------------------------------------------------------------

def write_outputs_quality(rows, winner, basis, args):
    os.makedirs(args.out, exist_ok=True)
    wp, wr = winner

    sweep_path = os.path.join(args.out, "quality_sweep_results.csv")
    with open(sweep_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank","mean_quality","min_quality","mean_ratio",
                    "descriptive_cv_percent","n_good"] + _PARAM_COLS)
        for rank, (p, r) in enumerate(rows, 1):
            w.writerow([rank, f"{r['mean_quality']:.4f}", f"{r['min_quality']:.4f}",
                        f"{r['mean_ratio']:.4f}", f"{r['desc_cv']:.3f}",
                        r["n_good"]] + _param_row(p))

    best_json = os.path.join(args.out, "best_config.json")
    with open(best_json, "w") as fh:
        json.dump({
            "selection_basis": basis,
            "summary": {k: v for k, v in wr.items() if k != "per_file"},
            "params": asdict(wp),
        }, fh, indent=2)

    perfile_path = os.path.join(args.out, "best_per_file.csv")
    with open(perfile_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file","ratio_A_over_B","area_A","area_B","quality",
                    "separation","fitR2_A","fitR2_B"])
        for fname, r in wr["per_file"]:
            if r is None:
                w.writerow([fname,"NA","NA","NA","NA","NA","NA","NA"])
                continue
            fmt = lambda x: f"{x:.5g}" if np.isfinite(x) else "NA"
            w.writerow([fname, fmt(r["ratio"]), fmt(r["aA"]), fmt(r["aB"]),
                        f"{r['quality']:.4f}", f"{r['sep']:.4f}",
                        f"{r['r2A']:.4f}", f"{r['r2B']:.4f}"])

    plot_path = None
    if _HAVE_MPL:
        try:
            qs = [r["mean_quality"] for _, r in rows]
            rr = [r["mean_ratio"]   for _, r in rows]
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.scatter(rr, qs, s=18, alpha=0.55, label="all param sets")
            ax.scatter([wr["mean_ratio"]], [wr["mean_quality"]], s=160,
                       marker="*", color="crimson", zorder=6, label="winner")
            ax.set_xlabel("Mean A/B ratio")
            ax.set_ylabel("Mean peak-fit quality")
            ax.set_title("UniDec parameter sweep: quality vs ratio")
            ax.legend(fontsize=8)
            fig.tight_layout()
            plot_path = os.path.join(args.out, "quality_vs_ratio.png")
            fig.savefig(plot_path, dpi=140)
            plt.close(fig)
        except Exception:
            pass

    return sweep_path, best_json, perfile_path, plot_path


# ---------------------------------------------------------------------------
# Output writing — BY-CONCENTRATION
# ---------------------------------------------------------------------------

def write_outputs_by_concentration(winners, args):
    os.makedirs(args.out, exist_ok=True)

    # One winning parameter set per concentration
    best_path = os.path.join(args.out, "by_concentration_best.csv")
    with open(best_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["concentration","conc_value","best_cv_percent","mean_ratio",
                    "std_ratio","n_good"] + _PARAM_COLS)
        for conc, gfiles, p, res, rows in winners:
            cv = conc_to_float(conc)
            w.writerow([conc, f"{cv:.6g}" if cv is not None else "NA",
                        f"{res['cv_percent']:.3f}", f"{res['mean_ratio']:.4f}",
                        f"{res['std_ratio']:.4f}", res["n_good"]] + _param_row(p))

    # Per-replicate ratios under each concentration's winner
    perfile_path = os.path.join(args.out, "by_concentration_per_file.csv")
    with open(perfile_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["concentration","file","ratio_A_over_B","area_A","area_B",
                    "baseline","r2"])
        for conc, gfiles, p, res, rows in winners:
            for name, ratio, aA, aB, bl, r2 in res["per_file"]:
                fmt = lambda x: f"{x:.5g}" if np.isfinite(x) else "NA"
                w.writerow([conc, name, fmt(ratio), fmt(aA), fmt(aB),
                            fmt(bl), fmt(r2)])

    # Winning config per concentration (JSON)
    best_json = os.path.join(args.out, "best_config_by_concentration.json")
    with open(best_json, "w") as fh:
        json.dump({
            conc: {
                "conc_value": conc_to_float(conc),
                "summary": {k: v for k, v in res.items() if k != "per_file"},
                "params": asdict(p),
            }
            for conc, gfiles, p, res, rows in winners
        }, fh, indent=2)

    # Full audit: every (concentration, param set) tried
    full_path = os.path.join(args.out, "by_concentration_sweep_full.csv")
    with open(full_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["concentration","cv_percent","mean_ratio","std_ratio",
                    "n_good"] + _PARAM_COLS)
        for conc, gfiles, p, res, rows in winners:
            for pp, rr in sorted(rows, key=lambda pr: pr[1]["cv_percent"]):
                w.writerow([conc, f"{rr['cv_percent']:.3f}",
                            f"{rr['mean_ratio']:.4f}", f"{rr['std_ratio']:.4f}",
                            rr["n_good"]] + _param_row(pp))

    # Plot: ratio vs concentration (the point of a concentration series)
    plot_path = None
    if _HAVE_MPL:
        try:
            pts = [(conc_to_float(conc), res["mean_ratio"], res["std_ratio"],
                    res["cv_percent"])
                   for conc, _, _, res, _ in winners
                   if conc_to_float(conc) is not None]
            pts.sort(key=lambda t: t[0])
            if pts:
                xs  = [t[0] for t in pts]
                ys  = [t[1] for t in pts]
                yer = [t[2] for t in pts]
                fig, ax = plt.subplots(figsize=(7, 5))
                ax.errorbar(xs, ys, yerr=yer, fmt="o-", color="#185FA5",
                            capsize=3, label="mean A/B ± SD")
                ax.set_xscale("log")
                ax.set_xlabel("Concentration (µg/mL)")
                ax.set_ylabel("A/B ratio (per-concentration optimum)")
                ax.set_title("A/B ratio vs concentration")
                for x, y, _, cv in pts:
                    ax.annotate(f"CV {cv:.1f}%", (x, y), fontsize=7,
                                textcoords="offset points", xytext=(4, 5))
                ax.legend(fontsize=8)
                fig.tight_layout()
                plot_path = os.path.join(args.out, "ratio_vs_concentration.png")
                fig.savefig(plot_path, dpi=140)
                plt.close(fig)
        except Exception:
            pass

    return best_path, best_json, perfile_path, full_path, plot_path


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------

def print_report(winner, basis):
    wp, wr = winner
    print("\n" + "="*70)
    print("WINNER (global — one shared parameter set)")
    print("="*70)
    print(f"  Basis          : {basis}")
    print(f"  Mean A/B ratio : {wr['mean_ratio']:.3f}")
    print(f"  %CV            : {wr['cv_percent']:.2f}%  (n={wr['n_good']})")
    print(f"  Mean R²        : {wr['mean_r2']:.4f}")
    print(f"\n  massbins       : {wp.massbins}")
    print(f"  mzsig          : {wp.mzsig}")
    print(f"  zzsig (charge) : {wp.zzsig}")
    print(f"  psig  (point)  : {wp.psig}")
    print(f"  beta           : {wp.beta}")
    print(f"  subtype        : {wp.subtype}  (0=none, 1=curved, 2=linear)")
    print(f"  subbuff        : {wp.subbuff}")
    print(f"  smooth         : {wp.smooth}")
    print(f"  integration    : {wp.integration_mode}")
    if wp.integration_mode == "gaussian":
        print(f"  peak_half_window: {wp.peak_half_window} Da")
    else:
        print(f"  baseline_mode  : {wp.baseline_mode}")
    print()
    print("  Per-replicate ratios:")
    for name, ratio, _, _, bl, _ in wr["per_file"]:
        s = f"{ratio:.3f}" if np.isfinite(ratio) else "FAILED"
        print(f"    {name:<38} {s}  (bl≈{bl:.4f})" if np.isfinite(bl) else
              f"    {name:<38} {s}")
    print()
    print("  ⚠  Always visually inspect the winning mass spectrum.")
    print("     A low CV achieved by over-smoothing A and B into one peak")
    print("     is a false win. Confirm both peaks are cleanly separated.")
    print("="*70)


def print_report_quality(winner, basis):
    wp, wr = winner
    print("\n" + "="*70)
    print("WINNER (quality — one shared parameter set)")
    print("="*70)
    print(f"  Basis          : {basis}")
    print(f"  Mean quality   : {wr['mean_quality']:.3f}  "
          f"(worst file {wr['min_quality']:.3f})")
    print(f"  Mean A/B ratio : {wr['mean_ratio']:.3f}")
    print(f"  Descriptive CV : {wr['desc_cv']:.2f}%  (reported only, NOT optimised)")
    print(f"  Files used     : {wr['n_good']}")
    print(f"\n  massbins       : {wp.massbins}")
    print(f"  mzsig          : {wp.mzsig}")
    print(f"  zzsig (charge) : {wp.zzsig}")
    print(f"  psig  (point)  : {wp.psig}")
    print(f"  beta           : {wp.beta}")
    print(f"  subtype        : {wp.subtype}  (0=none, 1=curved, 2=linear)")
    print(f"  subbuff        : {wp.subbuff}")
    print(f"  smooth         : {wp.smooth}")
    print(f"  integration    : {wp.integration_mode}")
    print()
    print("  Per-file under the winning conditions:")
    for fname, r in wr["per_file"]:
        if r is None:
            print(f"    {fname:<38} FAILED")
            continue
        print(f"    {fname:<38} ratio={r['ratio']:.3f}  Q={r['quality']:.3f} "
              f"(sep={r['sep']:.2f}, R²A={r['r2A']:.2f}, R²B={r['r2B']:.2f})")
    print()
    print("  ⚠  Quality rewards clean, well-separated A and B peaks, but always")
    print("     visually confirm the winning mass spectrum before trusting it.")
    print("="*70)


def print_report_per_file(winners):
    ratios = np.array([r["ratio"] for _, _, r in winners], float)
    mean_r = float(np.mean(ratios))
    std_r  = float(np.std(ratios, ddof=1)) if ratios.size > 1 else 0.0
    cv     = 100.0 * std_r / mean_r if mean_r else float("nan")

    print("\n" + "="*70)
    print("PER-FILE WINNERS (each file optimised independently)")
    print("="*70)
    print(f"  Files optimised : {len(winners)}")
    print(f"  Mean A/B ratio  : {mean_r:.3f}")
    print(f"  Descriptive CV  : {cv:.2f}%   (NOT the optimisation target —")
    print(f"                     each file used its own best parameters)")
    print()
    for fname, p, r in winners:
        print(f"  {fname:<38} ratio={r['ratio']:.3f}  Q={r['quality']:.3f}  "
              f"| zzsig={p.zzsig} psig={p.psig} beta={p.beta} "
              f"mzsig={p.mzsig} sub={p.subtype}/{p.subbuff}")
    print()
    print("  ⚠  Per-file tuning means the ratios come from DIFFERENT parameter")
    print("     sets. The quality metric is ratio-independent (fit R² + peak")
    print("     separation), so it will not bias the ratio, but differing")
    print("     smoothing between files can still introduce systematic offsets.")
    print("     Cross-check against the global-mode result before trusting it.")
    print("="*70)


def print_report_by_concentration(winners):
    print("\n" + "="*70)
    print("PER-CONCENTRATION WINNERS (minimise within-group cross-replicate %CV)")
    print("="*70)
    for conc, gfiles, p, res, rows in winners:
        print(f"\n  Concentration {conc}  (n={res['n_good']})")
        print(f"    Best %CV      : {res['cv_percent']:.2f}%")
        print(f"    Mean A/B ratio: {res['mean_ratio']:.3f}")
        print(f"    Params        : zzsig={p.zzsig} psig={p.psig} beta={p.beta} "
              f"mzsig={p.mzsig} sub={p.subtype}/{p.subbuff} smooth={p.smooth}")
        for name, ratio, _, _, bl, _ in res["per_file"]:
            s = f"{ratio:.3f}" if np.isfinite(ratio) else "FAILED"
            print(f"      {name:<40} {s}")
    print()
    print("  ⚠  Each concentration is optimised independently, so the winning")
    print("     parameters may differ between concentrations. Inspect the mass")
    print("     spectra and the ratio_vs_concentration.png trend before use.")
    print("="*70)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Optimise UniDec 8.x parameters for A/B ratio quantitation "
                    "of two close protein species, globally or per-file.")
    ap.add_argument("--data",  required=True,
                    help="Folder containing replicate .raw files")
    ap.add_argument("--out",   default="unidec_ratio_results")
    ap.add_argument("--per-file", action="store_true", dest="per_file",
                    help="Optimise parameters SEPARATELY for each file, each "
                         "maximising that file's own peak-fit quality "
                         "(default optimises one shared set by MEAN quality).")
    ap.add_argument("--cv-mode", action="store_true", dest="cv_mode",
                    help="Legacy objective: one shared parameter set that "
                         "minimises cross-replicate %%CV of the A/B ratio "
                         "(instead of the default quality objective).")
    ap.add_argument("--by-concentration", action="store_true",
                    dest="by_concentration",
                    help="Group replicates by a concentration token in the "
                         "filename (e.g. '1p00e-4') and, WITHIN each "
                         "concentration, find the parameter set that minimises "
                         "cross-replicate %%CV. One winner per concentration.")
    ap.add_argument("--conc-pattern", default=None, dest="conc_pattern",
                    help="Regex for the concentration token (default matches "
                         r"forms like '1p00e-4'/'2e-4').")
    ap.add_argument("--centroidA",   type=float, default=23412.0)
    ap.add_argument("--centroidB",   type=float, default=23658.0)
    ap.add_argument("--target-ratio",type=float, default=3.0, dest="target_ratio")
    ap.add_argument("--ratio-band",  type=float, nargs=2,
                    default=(1.8, 4.5), dest="ratio_band",
                    help="(global mode) accept only param sets whose mean "
                         "ratio is in this range")
    ap.add_argument("--integration-mode", choices=["window", "gaussian"],
                    default="window", dest="integration_mode",
                    help="How the mass spectrum is turned into a ratio.")
    ap.add_argument("--zrange",      type=int, nargs=2, default=(8, 45))
    ap.add_argument("--minmz",       type=float, default=600.0)
    ap.add_argument("--maxmz",       type=float, default=2500.0)
    ap.add_argument("--mass-pad",    type=float, default=910.0, dest="mass_pad",
                    help="Da of headroom beyond each centroid for masslb/massub")
    ap.add_argument("--lam",         type=float, default=1.0,
                    help="(global mode) penalty weight for ratio drift when "
                         "nothing lands in band")
    ap.add_argument("--quick",       action="store_true",
                    help="Run the small 16-combo grid (first pass)")
    ap.add_argument("--verbose",     action="store_true")
    return ap.parse_args(argv)


def list_files(data_dir):
    """List .raw / .mzML / .mzXML input spectra in data_dir."""
    exts = (".raw", ".mzml", ".mzxml")
    files = sorted(
        os.path.join(data_dir, f)
        for f in os.listdir(data_dir)
        if f.lower().endswith(exts) and os.path.isfile(os.path.join(data_dir, f))
    )
    if not files:
        raise SystemExit(
            f"No .raw/.mzML/.mzXML files found in {data_dir!r}\n"
            "If you have two-column text spectra, convert them to .mzML with "
            "msconvert."
        )
    return files


def main(argv=None):
    args = parse_args(argv)
    files = list_files(args.data)
    print(f"Found {len(files)} replicate spectra.")
    for f in files:
        print(f"  {os.path.basename(f)}")
    print()

    if args.by_concentration:
        winners = optimise_by_concentration(files, args)
        paths = write_outputs_by_concentration(winners, args)
        print_report_by_concentration(winners)
    elif args.cv_mode:
        rows, winner, basis = optimise(files, args)
        paths = write_outputs(rows, winner, basis, args)
        print_report(winner, basis)
    elif args.per_file:
        winners, full_records = optimise_per_file(files, args)
        paths = write_outputs_per_file(winners, full_records, args)
        print_report_per_file(winners)
    else:
        rows, winner, basis = optimise_global_quality(files, args)
        paths = write_outputs_quality(rows, winner, basis, args)
        print_report_quality(winner, basis)

    print(f"\nOutputs written to {args.out}\\")
    for pth in paths:
        if pth:
            print(f"  {os.path.basename(pth)}")


if __name__ == "__main__":
    main()
