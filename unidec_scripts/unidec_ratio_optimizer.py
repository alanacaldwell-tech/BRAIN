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
from datetime import datetime
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
    minmz: float = 700.0
    maxmz: float = 2000.0

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
    # NOTE: 23411.8 / 23658.7 are the MONOISOTOPIC masses; at ~23.5 kDa the
    # observed (deconvolved) centroid sits near the average mass, so we keep the
    # rounded values as window centers. peak_half_window (±20 Da) covers the gap.
    centroidA_expected: float = 23412.0   # test antibody light chain (Da)
    centroidB_expected: float = 23658.0   # internal-standard light chain (Da)
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


# Per-file engine cache. The RAW spectrum load (copy .raw + open_file) is
# parameter-INDEPENDENT, so we open each file exactly once and reuse the loaded
# engine across the whole parameter grid. Only process_data()/run_unidec() --
# which DO depend on the swept parameters -- are re-run per grid point.
_ENGINE_CACHE: dict[str, tuple] = {}   # path -> (engine, tmpdir)


def _open_engine(path: str, verbose: bool):
    """Copy .raw to an isolated temp dir and open it once. Returns (u, tmpdir)."""
    import io, contextlib, shutil, tempfile
    tmpdir = tempfile.mkdtemp(prefix="unidec_opt_")
    tmp_path = os.path.join(tmpdir, os.path.basename(path))
    shutil.copy2(path, tmp_path)
    u = make_engine()
    u.silent = not verbose
    if not verbose:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            u.open_file(tmp_path)
    else:
        u.open_file(tmp_path)
    return u, tmpdir


def cleanup_engines():
    """Tear down all cached engines and their temp dirs (call at end of run)."""
    import shutil
    for _u, tmpdir in _ENGINE_CACHE.values():
        shutil.rmtree(tmpdir, ignore_errors=True)
    _ENGINE_CACHE.clear()


def deconvolve_one(path: str, p: Params, verbose: bool = False,
                   reuse: bool = True):
    """
    Deconvolve one .raw file with parameter set `p`; return
    (mass, intensity, ratio, aA, aB, baseline, r_squared).

    If `reuse` (default), the file is opened once and the loaded engine is
    cached and reused for every subsequent parameter set on the same file --
    the raw load is parameter-independent, so this avoids re-reading the file
    729× per grid. Pass reuse=False (--fresh-engine) to open a throwaway engine
    per call, which exactly reproduces the previous behaviour.
    """
    import io, contextlib, shutil

    own_tmpdir = None
    if reuse:
        if path not in _ENGINE_CACHE:
            _ENGINE_CACHE[path] = _open_engine(path, verbose)
        u, _tmpdir = _ENGINE_CACHE[path]
    else:
        u, own_tmpdir = _open_engine(path, verbose)

    try:
        def _run():
            apply_params(u.config, p)
            u.process_data()    # re-derives processed data from cached rawdata
            u.run_unidec()      # the parameter-dependent deconvolution

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
        if own_tmpdir is not None:
            shutil.rmtree(own_tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Parameter grid
# ---------------------------------------------------------------------------

def build_grid(args) -> list[Params]:
    """
    Construct the sweep grid.

    The sweep covers peak width (mzsig), the deconvolution priors the app calls
    Charge Smooth Width (zzsig), Point Smooth Width (psig) and Beta, Gaussian
    smoothing (smooth), and the curved-subtraction width (subbuff). Background
    subtraction is fixed to CURVED (subtype=1) via Params default. massbins and
    mzbins are taken from the Params defaults (0.1 Da / 1.0 Th) rather than
    swept.
    """
    if args.quick:
        grid = dict(
            mzsig   = [0.3, 0.5],
            zzsig   = [1.0, 10.0],    # charge smooth width
            psig    = [0.0, 1.0],     # point smooth width
            beta    = [0.0, 50.0],    # softmax charge focusing
            smooth  = [1.0, 2.0],     # Gaussian smoothing
            subbuff = [100.0],        # curved-subtraction width
        )
    else:
        grid = dict(
            mzsig   = [0.3, 0.4, 0.5],
            zzsig   = [1.0, 5.0, 10.0],       # charge smooth width
            psig    = [0.0, 1.0, 2.0],        # point smooth width
            beta    = [0.0, 50.0, 100.0],     # softmax charge focusing
            smooth  = [1.0, 2.0, 5.0],        # Gaussian smoothing
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
            _, _, ratio, aA, aB, bl, r2 = deconvolve_one(
                path, p, args.verbose, reuse=not args.fresh_engine)
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
        mass, inten, ratio, aA, aB, bl, r2 = deconvolve_one(
            path, p, args.verbose, reuse=not args.fresh_engine)
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


def evaluate_group(p: Params, files: list[str], args) -> dict | None:
    """
    Evaluate ONE parameter set across a concentration group's replicates.

    Returns both the cross-replicate %CV of the A/B ratio AND the mean peak-fit
    quality (used as a floor so a degenerate low-CV set where A and B merged
    can't win). None if fewer than 2 replicates give a usable ratio.
    """
    ratios, quals, per_file = [], [], []
    for path in files:
        fname = os.path.basename(path)
        try:
            mass, inten, ratio, aA, aB, bl, r2 = deconvolve_one(
                path, p, args.verbose, reuse=not args.fresh_engine)
        except Exception as exc:
            per_file.append((fname, np.nan, np.nan, np.nan, np.nan, np.nan))
            if args.verbose:
                print(f"      [WARN] {fname}: {exc}")
            continue
        q = peak_quality(mass, inten, p, args.verbose)
        ratios.append(ratio)
        quals.append(q["scoreQ"] if q else np.nan)
        per_file.append((fname, ratio, aA, aB, bl, r2))

    ratios = np.array(ratios, float)
    good   = ratios[np.isfinite(ratios) & (ratios > 0)]
    if good.size < 2:
        return None

    mean_r = float(np.mean(good))
    std_r  = float(np.std(good, ddof=1))
    cv     = 100.0 * std_r / mean_r
    quals  = np.array(quals, float)
    has_q  = bool(np.isfinite(quals).any())
    return dict(cv_percent=cv, mean_ratio=mean_r, std_ratio=std_r,
                n_good=int(good.size),
                mean_quality=float(np.nanmean(quals)) if has_q else np.nan,
                min_quality=float(np.nanmin(quals)) if has_q else np.nan,
                per_file=per_file)


def _robust_pick(items, cv_of, q_of, args):
    """
    Overfitting-resistant selection from a list of candidates.

    1. Apply a peak-quality floor (drop candidates whose quality < args.quality_floor);
       if that empties the pool (or quality is unavailable), keep all candidates.
    2. Among candidates within args.cv_margin (fractional) of the best %CV,
       pick the one with the HIGHEST quality -- i.e. don't chase a razor-thin,
       possibly-overfit CV minimum; among near-best-CV sets take the cleanest,
       which is more likely to generalise to a future single sample.

    Returns the chosen item, or None if `items` is empty.
    """
    if not items:
        return None
    qs = [q_of(it) for it in items]
    have_q = any(np.isfinite(q) for q in qs)

    pool = items
    if have_q and args.quality_floor is not None:
        filt = [it for it in items
                if np.isfinite(q_of(it)) and q_of(it) >= args.quality_floor]
        if filt:
            pool = filt

    best_cv = min(cv_of(it) for it in pool)
    thresh  = best_cv * (1.0 + args.cv_margin)
    near    = [it for it in pool if cv_of(it) <= thresh + 1e-12]
    if have_q:
        return max(near, key=lambda it: q_of(it) if np.isfinite(q_of(it)) else -1.0)
    return min(near, key=cv_of)


def optimise_by_concentration(files, args):
    """
    Sweep the grid once per concentration group and report TWO answers side by
    side:

      * PER-CONCENTRATION: for each concentration, the (robustly selected)
        parameter set that minimises that concentration's cross-replicate %CV.
      * SHARED: one parameter set that works across ALL concentrations at once
        (minimises the aggregated within-concentration %CV). This is the more
        overfitting-resistant answer, since it is fit against every replicate
        of every concentration rather than the 2-3 replicates of one group.

    Both use _robust_pick (quality floor + CV-margin -> highest quality).
    Returns a result dict consumed by write_outputs / print_report.
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

    usable  = [c for c in sorted(groups, key=lambda c: (conc_to_float(c) or 0.0, c))
               if len(groups[c]) >= 2]
    skipped = [c for c in groups if len(groups[c]) < 2]
    for c in skipped:
        print(f"[WARN] {c}: only {len(groups[c])} replicate; need >=2 for a CV. Skipping.")
    if not usable:
        raise SystemExit("No concentration group has >=2 replicates.")

    grid = build_grid(args)
    print(f"[BY-CONCENTRATION] {len(usable)} usable concentration(s) × "
          f"{len(grid)} parameter sets:")
    for c in usable:
        print(f"        {c:<12} {len(groups[c])} replicate(s)")
    print()

    # matrix[conc][pid] = evaluate_group result (or None)
    matrix = {c: {} for c in usable}
    for c in usable:
        gfiles = groups[c]
        print(f"  sweeping concentration {c} ...")
        for pid, p in enumerate(grid):
            matrix[c][pid] = evaluate_group(p, gfiles, args)

    # --- per-concentration winners ---
    per_conc = {}
    for c in usable:
        rows = [(pid, grid[pid], res) for pid, res in matrix[c].items() if res]
        if not rows:
            print(f"  [WARN] {c}: all parameter sets failed.")
            continue
        pid, p, res = _robust_pick(
            rows, cv_of=lambda it: it[2]["cv_percent"],
            q_of=lambda it: it[2]["mean_quality"], args=args)
        per_conc[c] = (pid, p, res, rows)
        print(f"  [per-conc] {c:<12} CV={res['cv_percent']:.2f}%  "
              f"ratio={res['mean_ratio']:.3f}  Q={res['mean_quality']:.3f}")

    # --- shared winner: a pid present (non-None) in EVERY usable concentration ---
    shared_cands = []
    for pid, p in enumerate(grid):
        reslist = [matrix[c].get(pid) for c in usable]
        if any(r is None for r in reslist):
            continue
        cvs   = np.array([r["cv_percent"]  for r in reslist], float)
        quals = np.array([r["mean_quality"] for r in reslist], float)
        shared_cands.append(dict(
            pid=pid, p=p,
            agg_cv=float(np.mean(cvs)), worst_cv=float(np.max(cvs)),
            mean_quality=float(np.nanmean(quals)) if np.isfinite(quals).any() else np.nan,
            per_conc={c: matrix[c][pid] for c in usable}))
    if not shared_cands:
        print("[WARN] no single parameter set succeeded on every concentration; "
              "shared-conditions result unavailable.")
        shared = None
    else:
        shared = _robust_pick(
            shared_cands, cv_of=lambda it: it["agg_cv"],
            q_of=lambda it: it["mean_quality"], args=args)
        print(f"\n  [shared]  agg CV={shared['agg_cv']:.2f}%  "
              f"(worst {shared['worst_cv']:.2f}%)  Q={shared['mean_quality']:.3f}  "
              f"| {shared['p'].key_str()}")

    return dict(usable=usable, skipped=skipped, groups=groups, grid=grid,
                matrix=matrix, per_conc=per_conc, shared=shared,
                ungrouped=ungrouped,
                _qfloor=args.quality_floor, _cvmargin=args.cv_margin)


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
                    "cv_percent","n_good"] + _PARAM_COLS)
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

def _linfit_r2(x, y):
    """Linear fit y=slope*x+intercept; return (slope, intercept, r2) or Nones."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if x.size < 2 or np.ptp(x) == 0:
        return None, None, None
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), float(r2)


def write_outputs_by_concentration(result, args):
    os.makedirs(args.out, exist_ok=True)
    usable  = result["usable"]
    matrix  = result["matrix"]
    per_conc = result["per_conc"]         # conc -> (pid, p, res, rows)
    shared   = result["shared"]           # dict or None
    fmt = lambda x: f"{x:.5g}" if np.isfinite(x) else "NA"

    def sh_res(c):
        return shared["per_conc"][c] if shared else None

    # --- side-by-side comparison (the headline) ---
    cmp_path = os.path.join(args.out, "comparison.csv")
    with open(cmp_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["concentration","conc_value",
                    "perconc_cv","perconc_ratio","perconc_quality",
                    "shared_cv","shared_ratio","shared_quality",
                    "params_differ"])
        for c in usable:
            pc = per_conc.get(c)
            sr = sh_res(c)
            cv = conc_to_float(c)
            differ = ("NA" if (pc is None or shared is None)
                      else str(pc[1].key_str() != shared["p"].key_str()))
            row = [c, f"{cv:.6g}" if cv is not None else "NA"]
            row += ([fmt(pc[2]["cv_percent"]), fmt(pc[2]["mean_ratio"]),
                     fmt(pc[2]["mean_quality"])] if pc else ["NA","NA","NA"])
            row += ([fmt(sr["cv_percent"]), fmt(sr["mean_ratio"]),
                     fmt(sr["mean_quality"])] if sr else ["NA","NA","NA"])
            row += [differ]
            w.writerow(row)

    # --- per-concentration winning params ---
    perconc_path = os.path.join(args.out, "per_concentration_best.csv")
    with open(perconc_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["concentration","conc_value","cv_percent","mean_ratio",
                    "mean_quality","n_good"] + _PARAM_COLS)
        for c in usable:
            pc = per_conc.get(c)
            if not pc:
                continue
            _pid, p, res, _rows = pc
            cv = conc_to_float(c)
            w.writerow([c, f"{cv:.6g}" if cv is not None else "NA",
                        f"{res['cv_percent']:.3f}", f"{res['mean_ratio']:.4f}",
                        f"{res['mean_quality']:.4f}", res["n_good"]] + _param_row(p))

    # --- shared conditions (single set) + calibration linearity ---
    shared_path = os.path.join(args.out, "shared_conditions.csv")
    sh_slope = sh_int = sh_r2 = None
    if shared:
        xs = [conc_to_float(c) for c in usable if conc_to_float(c) is not None]
        ys = [shared["per_conc"][c]["mean_ratio"] for c in usable
              if conc_to_float(c) is not None]
        sh_slope, sh_int, sh_r2 = _linfit_r2(xs, ys)
        with open(shared_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["agg_cv_percent","worst_cv_percent","mean_quality",
                        "calibration_slope","calibration_intercept",
                        "calibration_r2"] + _PARAM_COLS)
            w.writerow([f"{shared['agg_cv']:.3f}", f"{shared['worst_cv']:.3f}",
                        f"{shared['mean_quality']:.4f}",
                        fmt(sh_slope) if sh_slope is not None else "NA",
                        fmt(sh_int) if sh_int is not None else "NA",
                        fmt(sh_r2) if sh_r2 is not None else "NA"]
                       + _param_row(shared["p"]))
    else:
        shared_path = None

    # --- per-replicate ratios (both selections) ---
    perfile_path = os.path.join(args.out, "per_file_ratios.csv")
    with open(perfile_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["selection","concentration","file","ratio_A_over_B",
                    "area_A","area_B","baseline","r2"])
        for c in usable:
            pc = per_conc.get(c)
            if pc:
                for name, ratio, aA, aB, bl, r2 in pc[2]["per_file"]:
                    w.writerow(["per_conc", c, name, fmt(ratio), fmt(aA),
                                fmt(aB), fmt(bl), fmt(r2)])
            if shared:
                for name, ratio, aA, aB, bl, r2 in shared["per_conc"][c]["per_file"]:
                    w.writerow(["shared", c, name, fmt(ratio), fmt(aA),
                                fmt(aB), fmt(bl), fmt(r2)])

    # --- JSON configs ---
    json_perconc = os.path.join(args.out, "best_config_by_concentration.json")
    with open(json_perconc, "w") as fh:
        json.dump({
            c: {"conc_value": conc_to_float(c),
                "summary": {k: v for k, v in per_conc[c][2].items() if k != "per_file"},
                "params": asdict(per_conc[c][1])}
            for c in usable if per_conc.get(c)
        }, fh, indent=2)

    json_shared = None
    if shared:
        json_shared = os.path.join(args.out, "best_config_shared.json")
        with open(json_shared, "w") as fh:
            json.dump({
                "agg_cv_percent": shared["agg_cv"],
                "worst_cv_percent": shared["worst_cv"],
                "mean_quality": shared["mean_quality"],
                "calibration": {"slope": sh_slope, "intercept": sh_int, "r2": sh_r2},
                "per_concentration": {
                    c: {k: v for k, v in shared["per_conc"][c].items() if k != "per_file"}
                    for c in usable},
                "params": asdict(shared["p"]),
            }, fh, indent=2)

    # --- full audit: every (concentration, param set) ---
    full_path = os.path.join(args.out, "by_concentration_sweep_full.csv")
    with open(full_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["concentration","cv_percent","mean_ratio","mean_quality",
                    "n_good"] + _PARAM_COLS)
        for c in usable:
            rows = [(pid, res) for pid, res in matrix[c].items() if res]
            for pid, res in sorted(rows, key=lambda t: t[1]["cv_percent"]):
                w.writerow([c, f"{res['cv_percent']:.3f}",
                            f"{res['mean_ratio']:.4f}",
                            f"{res['mean_quality']:.4f}", res["n_good"]]
                           + _param_row(result["grid"][pid]))

    # --- plot: per-conc vs shared calibration curves ---
    plot_path = None
    if _HAVE_MPL:
        try:
            fig, ax = plt.subplots(figsize=(7.5, 5))
            xs_pc = [conc_to_float(c) for c in usable if per_conc.get(c)
                     and conc_to_float(c) is not None]
            ys_pc = [per_conc[c][2]["mean_ratio"] for c in usable if per_conc.get(c)
                     and conc_to_float(c) is not None]
            er_pc = [per_conc[c][2]["std_ratio"] for c in usable if per_conc.get(c)
                     and conc_to_float(c) is not None]
            if xs_pc:
                ax.errorbar(xs_pc, ys_pc, yerr=er_pc, fmt="o", color="#1D9E75",
                            capsize=3, label="per-concentration optimum")
            if shared:
                xs_s = [conc_to_float(c) for c in usable if conc_to_float(c) is not None]
                ys_s = [shared["per_conc"][c]["mean_ratio"] for c in usable
                        if conc_to_float(c) is not None]
                er_s = [shared["per_conc"][c]["std_ratio"] for c in usable
                        if conc_to_float(c) is not None]
                order = np.argsort(xs_s)
                xs_s = list(np.array(xs_s)[order]); ys_s = list(np.array(ys_s)[order])
                er_s = list(np.array(er_s)[order])
                ax.errorbar(xs_s, ys_s, yerr=er_s, fmt="s-", color="#185FA5",
                            capsize=3, label="shared conditions")
                if sh_slope is not None:
                    xx = np.array([min(xs_s), max(xs_s)])
                    ax.plot(xx, sh_slope * xx + sh_int, "--", color="crimson",
                            lw=1, label=f"linear fit R²={sh_r2:.3f}")
            ax.set_xlabel("Concentration (µg/mL)")
            ax.set_ylabel("A/B ratio (test / internal standard)")
            ax.set_title("Calibration: per-concentration vs shared conditions")
            ax.legend(fontsize=8)
            fig.tight_layout()
            plot_path = os.path.join(args.out, "ratio_vs_concentration.png")
            fig.savefig(plot_path, dpi=140)
            plt.close(fig)
        except Exception:
            pass

    return (cmp_path, perconc_path, shared_path, perfile_path,
            json_perconc, json_shared, full_path, plot_path)


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


def print_report_by_concentration(result):
    usable   = result["usable"]
    per_conc = result["per_conc"]
    shared   = result["shared"]

    print("\n" + "="*70)
    print("BY-CONCENTRATION — per-concentration vs shared conditions")
    print("="*70)

    # Side-by-side table
    print(f"\n  {'conc':<10} {'per-conc CV':>12} {'ratio':>8}   "
          f"{'shared CV':>10} {'ratio':>8}   {'params differ':>13}")
    for c in usable:
        pc = per_conc.get(c)
        sr = shared["per_conc"][c] if shared else None
        pc_cv = f"{pc[2]['cv_percent']:.2f}%" if pc else "NA"
        pc_r  = f"{pc[2]['mean_ratio']:.3f}"   if pc else "NA"
        sh_cv = f"{sr['cv_percent']:.2f}%"     if sr else "NA"
        sh_r  = f"{sr['mean_ratio']:.3f}"      if sr else "NA"
        differ = ("NA" if (pc is None or shared is None)
                  else ("yes" if pc[1].key_str() != shared["p"].key_str() else "no"))
        print(f"  {c:<10} {pc_cv:>12} {pc_r:>8}   {sh_cv:>10} {sh_r:>8}   {differ:>13}")

    if shared:
        wp = shared["p"]
        print(f"\n  SHARED conditions (one set for all concentrations):")
        print(f"    agg CV {shared['agg_cv']:.2f}%  (worst {shared['worst_cv']:.2f}%)  "
              f"mean Q {shared['mean_quality']:.3f}")
        print(f"    zzsig={wp.zzsig} psig={wp.psig} beta={wp.beta} mzsig={wp.mzsig} "
              f"smooth={wp.smooth} sub={wp.subtype}/{wp.subbuff}")

    print()
    print("  Selection: quality floor + within CV-margin -> highest quality")
    print(f"    (--quality-floor {result.get('_qfloor','?')}, "
          f"--cv-margin {result.get('_cvmargin','?')})")
    print("  ⚠  With 2-3 replicates the per-concentration CV is a fragile target;")
    print("     the SHARED result is fit against every replicate and is the more")
    print("     trustworthy basis for future single-sample runs. If 'params")
    print("     differ' is 'no' across the board, the shared set is clearly right.")
    print("     Check ratio_vs_concentration.png for calibration linearity.")
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
    ap.add_argument("--out",   default="unidec_ratio_results",
                    help="Output directory base name. A _YYYYmmdd_HHMMSS "
                         "timestamp is appended unless --no-timestamp is set.")
    ap.add_argument("--no-timestamp", action="store_true", dest="no_timestamp",
                    help="Do NOT append a timestamp to --out (may overwrite "
                         "a previous run).")
    ap.add_argument("--fresh-engine", action="store_true", dest="fresh_engine",
                    help="Open a new UniDec engine for every parameter set "
                         "instead of loading each .raw once and reusing it "
                         "across the grid. Slower; use only to rule out engine "
                         "state carry-over.")
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
    ap.add_argument("--quality-floor", type=float, default=0.5, dest="quality_floor",
                    help="(by-concentration) minimum mean peak-fit quality a "
                         "parameter set must reach to be eligible; guards against "
                         "degenerate low-CV winners where A and B merged "
                         "(default 0.5). Set 0 to disable.")
    ap.add_argument("--cv-margin", type=float, default=0.25, dest="cv_margin",
                    help="(by-concentration) among sets within this fractional "
                         "margin of the best %%CV, pick the highest-quality one "
                         "rather than the razor-thin CV minimum (default 0.25 = "
                         "within 25%%). Improves generalisation.")
    ap.add_argument("--centroidA",   type=float, default=23412.0,
                    help="Test antibody light-chain mass (Da). Ratio is A/B.")
    ap.add_argument("--centroidB",   type=float, default=23658.0,
                    help="Internal-standard light-chain mass (Da).")
    ap.add_argument("--target-ratio",type=float, default=3.0, dest="target_ratio")
    ap.add_argument("--ratio-band",  type=float, nargs=2,
                    default=(1.8, 4.5), dest="ratio_band",
                    help="(global mode) accept only param sets whose mean "
                         "ratio is in this range")
    ap.add_argument("--integration-mode", choices=["window", "gaussian"],
                    default="window", dest="integration_mode",
                    help="How the mass spectrum is turned into a ratio.")
    ap.add_argument("--zrange",      type=int, nargs=2, default=(8, 45))
    ap.add_argument("--minmz",       type=float, default=700.0)
    ap.add_argument("--maxmz",       type=float, default=2000.0)
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

    # Timestamp the output directory so repeat runs don't overwrite each other.
    if not args.no_timestamp:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.out = f"{args.out.rstrip('/' + chr(92))}_{ts}"
    print(f"Output directory: {args.out}\n")

    files = list_files(args.data)
    print(f"Found {len(files)} replicate spectra.")
    for f in files:
        print(f"  {os.path.basename(f)}")
    print()

    if not args.fresh_engine:
        print("Engine reuse: each .raw is loaded once and reused across the "
              "grid (use --fresh-engine to disable).\n")

    try:
        if args.by_concentration:
            result = optimise_by_concentration(files, args)
            paths = write_outputs_by_concentration(result, args)
            print_report_by_concentration(result)
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
    finally:
        cleanup_engines()

    print(f"\nOutputs written to {args.out}\\")
    for pth in paths:
        if pth:
            print(f"  {os.path.basename(pth)}")


if __name__ == "__main__":
    main()
