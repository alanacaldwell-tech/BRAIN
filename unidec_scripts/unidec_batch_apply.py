#!/usr/bin/env python3
"""
unidec_batch_apply.py
=====================
Batch-process .raw files through UniDec using a FIXED set of conditions (no
optimisation). Point it at a `shared_conditions.csv` (or `best_config_shared.json`)
produced by unidec_ratio_optimizer.py --by-concentration, and it deconvolves
every file with those exact parameters and reports the test/IS ratio per file.

This is the production step: once the optimiser has found reproducible
conditions, use this to apply them to new data (including single samples per
concentration).

Reuses the deconvolution / ratio / quality / grouping code from
unidec_ratio_optimizer.py so the two stay in lockstep.

USAGE
-----
  # apply shared conditions to a folder of .raw files:
  python unidec_batch_apply.py --config results_.../shared_conditions.csv \\
      --data "D:\\path\\to\\raw" --out applied

  # exact reproduction from the full JSON config (includes geometry):
  python unidec_batch_apply.py --config results_.../best_config_shared.json \\
      --data "D:\\..." --out applied

CONFIG FILE
-----------
  shared_conditions.csv : contains the TUNED parameters (massbins, mzsig, zzsig,
      psig, beta, subtype, subbuff, smooth, integration_mode, baseline_mode,
      peak_half_window). Geometry NOT stored there (centroids, m/z window, mass
      bounds, charge range) is taken from the CLI defaults below -- pass the
      same --centroidA/--centroidB/--minmz/--maxmz/--zrange you optimised with
      if you changed them.
  best_config_shared.json : contains the COMPLETE parameter set (geometry
      included) -- exact reproduction, no CLI geometry needed.

OUTPUTS
-------
  batch_per_file.csv          ratio + peak quality for every file
  batch_by_concentration.csv  per-concentration mean/SD/%CV (if filenames carry
                              concentration tokens like '1p00e-4')
  batch_summary.txt           the conditions used + overall stats + calibration
  ratio_vs_concentration.png  calibration curve (if grouped)
"""

from __future__ import annotations
import argparse, csv, json, os, sys
from dataclasses import asdict

# --- import the optimiser module (same directory) for shared machinery -------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import unidec_ratio_optimizer as opt
except ImportError as e:
    raise SystemExit(f"Cannot import unidec_ratio_optimizer (must sit beside "
                     f"this script): {e}")

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:
    _HAVE_MPL = False


# ---------------------------------------------------------------------------
# Build a Params object from a config file
# ---------------------------------------------------------------------------

# The tuned parameter columns that shared_conditions.csv carries.
_TUNED = {
    "massbins": float, "mzsig": float, "zzsig": float, "psig": float,
    "beta": float, "subtype": int, "subbuff": float, "smooth": float,
    "integration_mode": str, "baseline_mode": str, "peak_half_window": float,
}


def _apply_geometry(p: "opt.Params", args):
    """Fill non-swept geometry the same way build_grid() did."""
    p.centroidA_expected = args.centroidA
    p.centroidB_expected = args.centroidB
    p.minmz  = args.minmz
    p.maxmz  = args.maxmz
    p.startz = args.zrange[0]
    p.endz   = args.zrange[1]
    p.masslb = args.centroidA - args.mass_pad
    p.massub = args.centroidB + args.mass_pad
    p.int_lo_A  = args.centroidA - 210
    p.int_split = 0.5 * (args.centroidA + args.centroidB)
    p.int_hi_B  = args.centroidB + 240
    return p


def params_from_config(path, args):
    """Return (Params, source_label). Accepts shared_conditions.csv or a JSON."""
    if path.lower().endswith(".json"):
        with open(path) as fh:
            data = json.load(fh)
        pdict = data.get("params", data)   # accept raw params dict too
        # Keep only real Params fields, in case of extra keys.
        fields = set(asdict(opt.Params()).keys())
        pdict = {k: v for k, v in pdict.items() if k in fields}
        return opt.Params(**pdict), f"{os.path.basename(path)} (full config)"

    # CSV: read the single data row, coerce the tuned columns, fill geometry.
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{path} has no data rows.")
    row = rows[0]
    missing = [k for k in _TUNED if k not in row]
    if missing:
        raise SystemExit(f"{path} is missing expected columns: {missing}. "
                         f"Is this a shared_conditions.csv?")
    p = opt.Params()
    for k, caster in _TUNED.items():
        try:
            setattr(p, k, caster(row[k]))
        except (ValueError, TypeError) as e:
            raise SystemExit(f"Bad value for {k!r} in {path}: {row[k]!r} ({e})")
    _apply_geometry(p, args)
    return p, f"{os.path.basename(path)} (tuned params + CLI geometry)"


# ---------------------------------------------------------------------------
# Batch run
# ---------------------------------------------------------------------------

def process_file(path, p, args):
    """Deconvolve one file with fixed params; return a result dict (or error)."""
    mass, inten, ratio, aA, aB, bl, r2 = opt.deconvolve_one(
        path, p, args.verbose, reuse=not args.fresh_engine)
    q = opt.peak_quality(mass, inten, p, args.verbose)
    return dict(
        ratio=ratio, aA=aA, aB=aB, bl=bl, r2=r2,
        quality=(q["scoreQ"] if q else np.nan),
        sep=(q["sep"] if q else np.nan),
        r2A=(q["r2A"] if q else np.nan),
        r2B=(q["r2B"] if q else np.nan),
    )


def run(files, p, args):
    results = []   # (fname, conc_token, res_or_None)
    for path in files:
        fname = os.path.basename(path)
        conc  = opt.parse_concentration(fname)
        try:
            res = process_file(path, p, args)
            results.append((fname, conc, res))
            print(f"  {fname:<44} ratio={res['ratio']:.4f}  Q={res['quality']:.3f}")
        except Exception as exc:
            results.append((fname, conc, None))
            print(f"  [FAIL] {fname}: {exc}")
            if args.verbose:
                import traceback; traceback.print_exc()
    return results


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_outputs(results, p, source, args):
    os.makedirs(args.out, exist_ok=True)
    fmt = lambda x: f"{x:.5g}" if (x is not None and np.isfinite(x)) else "NA"

    # Per-file
    perfile = os.path.join(args.out, "batch_per_file.csv")
    with open(perfile, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "concentration", "conc_value", "ratio_A_over_B",
                    "area_A", "area_B", "quality", "separation",
                    "fitR2_A", "fitR2_B", "baseline", "unidec_r2", "status"])
        for fname, conc, res in results:
            cv = opt.conc_to_float(conc) if conc else None
            if res is None:
                w.writerow([fname, conc or "NA", fmt(cv), "NA","NA","NA","NA",
                            "NA","NA","NA","NA","NA", "FAILED"])
                continue
            w.writerow([fname, conc or "NA", fmt(cv),
                        fmt(res["ratio"]), fmt(res["aA"]), fmt(res["aB"]),
                        fmt(res["quality"]), fmt(res["sep"]),
                        fmt(res["r2A"]), fmt(res["r2B"]),
                        fmt(res["bl"]), fmt(res["r2"]), "ok"])

    # Per-concentration (if any tokens present)
    grouped = {}
    for fname, conc, res in results:
        if conc and res is not None and np.isfinite(res["ratio"]) and res["ratio"] > 0:
            grouped.setdefault(conc, []).append(res)

    byconc = None
    cal_slope = cal_int = cal_r2 = None
    if grouped:
        byconc = os.path.join(args.out, "batch_by_concentration.csv")
        with open(byconc, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["concentration","conc_value","n","mean_ratio",
                        "std_ratio","cv_percent","mean_quality"])
            xs, ys = [], []
            for conc in sorted(grouped, key=lambda c: (opt.conc_to_float(c) or 0.0, c)):
                rs = grouped[conc]
                ratios = np.array([r["ratio"] for r in rs], float)
                quals  = np.array([r["quality"] for r in rs], float)
                mean_r = float(np.mean(ratios))
                std_r  = float(np.std(ratios, ddof=1)) if ratios.size > 1 else 0.0
                cvp    = 100.0 * std_r / mean_r if mean_r else float("nan")
                cv_val = opt.conc_to_float(conc)
                w.writerow([conc, f"{cv_val:.6g}" if cv_val is not None else "NA",
                            ratios.size, f"{mean_r:.4f}", f"{std_r:.4f}",
                            f"{cvp:.3f}", f"{np.nanmean(quals):.4f}"])
                if cv_val is not None:
                    xs.append(cv_val); ys.append(mean_r)
            cal_slope, cal_int, cal_r2 = opt._linfit_r2(xs, ys)

    # Summary
    summary = os.path.join(args.out, "batch_summary.txt")
    ok = [r for _, _, r in results if r is not None]
    with open(summary, "w") as fh:
        fh.write("UniDec batch apply — fixed-condition processing\n")
        fh.write("=" * 52 + "\n")
        fh.write(f"Conditions source : {source}\n")
        fh.write(f"Files processed   : {len(ok)}/{len(results)} succeeded\n\n")
        fh.write("Parameters applied:\n")
        for k in ("massbins","mzsig","zzsig","psig","beta","subtype","subbuff",
                  "smooth","integration_mode","baseline_mode","peak_half_window",
                  "minmz","maxmz","masslb","massub","startz","endz",
                  "centroidA_expected","centroidB_expected"):
            fh.write(f"  {k:<20} {getattr(p, k)}\n")
        if cal_r2 is not None:
            fh.write("\nCalibration (mean ratio vs concentration):\n")
            fh.write(f"  slope     {cal_slope:.6g}\n")
            fh.write(f"  intercept {cal_int:.6g}\n")
            fh.write(f"  R^2       {cal_r2:.5f}\n")

    # Plot
    plot = None
    if _HAVE_MPL and grouped:
        try:
            pts = []
            for conc in grouped:
                cv_val = opt.conc_to_float(conc)
                if cv_val is None:
                    continue
                ratios = np.array([r["ratio"] for r in grouped[conc]], float)
                pts.append((cv_val, float(np.mean(ratios)),
                            float(np.std(ratios, ddof=1)) if ratios.size > 1 else 0.0))
            pts.sort(key=lambda t: t[0])
            if pts:
                xs = [t[0] for t in pts]; ys = [t[1] for t in pts]; er = [t[2] for t in pts]
                fig, ax = plt.subplots(figsize=(7, 5))
                ax.errorbar(xs, ys, yerr=er, fmt="o-", color="#185FA5",
                            capsize=3, label="mean A/B ± SD")
                if cal_slope is not None:
                    xx = np.array([min(xs), max(xs)])
                    ax.plot(xx, cal_slope * xx + cal_int, "--", color="crimson",
                            lw=1, label=f"linear fit R²={cal_r2:.3f}")
                ax.set_xlabel("Concentration (µg/mL)")
                ax.set_ylabel("A/B ratio (test / internal standard)")
                ax.set_title("Batch-applied conditions: calibration")
                ax.legend(fontsize=8)
                fig.tight_layout()
                plot = os.path.join(args.out, "ratio_vs_concentration.png")
                fig.savefig(plot, dpi=140)
                plt.close(fig)
        except Exception:
            pass

    return perfile, byconc, summary, plot


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Batch-process .raw files through UniDec with FIXED "
                    "conditions from shared_conditions.csv (no optimisation).")
    ap.add_argument("--config", required=True,
                    help="shared_conditions.csv or best_config_shared.json")
    ap.add_argument("--data", required=True, help="Folder of .raw files")
    ap.add_argument("--out", default="unidec_batch_results",
                    help="Output directory base name (timestamp appended "
                         "unless --no-timestamp).")
    ap.add_argument("--no-timestamp", action="store_true", dest="no_timestamp")
    ap.add_argument("--fresh-engine", action="store_true", dest="fresh_engine",
                    help="Open a new engine per file instead of reusing one.")
    # Geometry (only used when --config is a CSV; must match what you optimised)
    ap.add_argument("--centroidA", type=float, default=23412.0,
                    help="Test light-chain mass (Da). Ratio is A/B.")
    ap.add_argument("--centroidB", type=float, default=23658.0,
                    help="Internal-standard light-chain mass (Da).")
    ap.add_argument("--minmz", type=float, default=700.0)
    ap.add_argument("--maxmz", type=float, default=2000.0)
    ap.add_argument("--zrange", type=int, nargs=2, default=(8, 45))
    ap.add_argument("--mass-pad", type=float, default=910.0, dest="mass_pad")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args(argv)


def main(argv=None):
    from datetime import datetime
    args = parse_args(argv)
    if not args.no_timestamp:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.out = f"{args.out.rstrip('/' + chr(92))}_{ts}"

    p, source = params_from_config(args.config, args)
    files = opt.list_files(args.data)

    print(f"Conditions: {source}")
    print(f"  massbins={p.massbins} mzsig={p.mzsig} zzsig={p.zzsig} psig={p.psig} "
          f"beta={p.beta} sub={p.subtype}/{p.subbuff} smooth={p.smooth} "
          f"int={p.integration_mode}")
    print(f"Output directory: {args.out}")
    print(f"Processing {len(files)} file(s)"
          + ("" if args.fresh_engine else "  (engine reuse on)") + "\n")

    try:
        results = run(files, p, args)
    finally:
        opt.cleanup_engines()

    perfile, byconc, summary, plot = write_outputs(results, p, source, args)

    ok = sum(1 for _, _, r in results if r is not None)
    print(f"\nDone: {ok}/{len(results)} files processed.")
    print(f"Outputs in {args.out}\\")
    for pth in (perfile, byconc, summary, plot):
        if pth:
            print(f"  {os.path.basename(pth)}")


if __name__ == "__main__":
    main()
