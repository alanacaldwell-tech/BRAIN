#!/usr/bin/env python3
"""
mz_ratio_direct.py
==================
Direct m/z-domain quantitation of two co-eluting protein species from
native-MS Orbitrap profile data (.raw files), without deconvolution.

Why no deconvolution
--------------------
UniDec collapses charge-state ladders into a zero-charge mass spectrum by
fitting a kernel model across all charge states simultaneously. For two
species with DIFFERENT charge envelopes (species A peaks at z~21, B at z~18),
the model redistributes intensity in a kernel-width-dependent way that
produces ratios of 1.4–8.2 depending on parameters — unstable and not
anchored to the raw signal.

Direct integration avoids this entirely: for each charge state z in a
user-specified "safe zone" (default z=15–22, where both envelopes have
strong, overlapping signal), we:
  1. Compute the expected m/z for each species:
       mz = (mass + z * adduct_mass) / z
  2. Integrate the raw profile intensity in a small window (±half_width m/z)
     around each expected position, after subtracting a local linear baseline
     estimated from flanking regions immediately outside the window.
  3. Sum integrated areas across all charge states:
       ratio = sum_z(area_A(z)) / sum_z(area_B(z))

This mirrors what CDMS does conceptually — accumulate counts (here: intensity)
across all charge states that contribute signal, without a model assumption
about the peak shape across z.

Usage
-----
    python mz_ratio_direct.py --data D:\\path\\to\\folder --out results_direct

    # custom masses / charge range / window:
    python mz_ratio_direct.py --data D:\\...\\test --out results_direct ^
        --massA 23412 --massB 23658 ^
        --zmin 15 --zmax 22 ^
        --half-width 1.5 ^
        --baseline-flanks 3.0

Outputs
-------
    direct_per_file.csv    per-replicate ratio + per-z breakdown
    direct_summary.txt     mean, std, %CV, winning parameters
    direct_zprofile.csv    mean area_A and area_B per z (diagnostic)
"""

import argparse, csv, os, sys, warnings
from dataclasses import dataclass
import numpy as np

warnings.filterwarnings("ignore")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:
    _HAVE_MPL = False

ADDUCT = 1.007276   # proton mass (Da)

# ---------------------------------------------------------------------------
# Engine wrapper
# ---------------------------------------------------------------------------

def make_engine():
    try:
        import unidec.engine as _eng
    except ImportError as e:
        raise SystemExit(f"Cannot import unidec.engine: {e}")
    saved = sys.argv[:]
    sys.argv = sys.argv[:1]
    try:
        return _eng.UniDec()
    finally:
        sys.argv = saved


def load_raw_mz(path):
    """Open a .raw file and return the raw m/z profile as (mz, intensity) arrays."""
    import io, contextlib
    u = make_engine()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        u.open_file(path)
    raw = np.asarray(u.data.rawdata)
    if raw.ndim != 2 or raw.shape[1] < 2:
        raise RuntimeError(f"Unexpected rawdata shape {raw.shape} for {path}")
    mz    = raw[:, 0]
    inten = raw[:, 1]
    order = np.argsort(mz)
    return mz[order], inten[order]


# ---------------------------------------------------------------------------
# Per-charge-state integration
# ---------------------------------------------------------------------------

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


def integrate_peak(mz, inten, center, half_width, baseline_flanks):
    """
    Integrate one isotope envelope at `center` m/z.

    Window:   [center - half_width, center + half_width]
    Baseline: median of two flanking regions immediately outside the window,
              each of width `baseline_flanks` m/z, linearly interpolated
              across the integration window and subtracted before summing.

    Returns (area, baseline_mean, n_points_in_window).
    """
    lo = center - half_width
    hi = center + half_width

    sel = (mz >= lo) & (mz <= hi)
    if sel.sum() < 3:
        return 0.0, 0.0, 0

    # Flanking baseline regions
    left_sel  = (mz >= lo - baseline_flanks) & (mz < lo)
    right_sel = (mz > hi) & (mz <= hi + baseline_flanks)
    bl_left  = np.median(inten[left_sel])  if left_sel.sum() > 0 else inten[sel][0]
    bl_right = np.median(inten[right_sel]) if right_sel.sum() > 0 else inten[sel][-1]

    # Linear baseline across the window
    x_win = mz[sel]
    bl_arr = np.interp(x_win, [lo, hi], [bl_left, bl_right])
    corrected = np.maximum(inten[sel] - bl_arr, 0.0)
    area = float(_trapz(corrected, x_win))
    return area, float(0.5 * (bl_left + bl_right)), int(sel.sum())


def quantitate_one_file(path, args):
    """
    Load raw m/z data and compute ratio A/B by direct per-z integration.

    Returns a dict with ratio, per-z breakdown, and diagnostics.
    """
    mz, inten = load_raw_mz(path)

    z_results = []
    for z in range(args.zmin, args.zmax + 1):
        mzA = (args.massA + z * ADDUCT) / z
        mzB = (args.massB + z * ADDUCT) / z

        # Skip charge states where expected peaks fall outside data range
        margin = args.half_width + args.baseline_flanks + 1.0
        if mzA < mz.min() + margin or mzA > mz.max() - margin:
            continue
        if mzB < mz.min() + margin or mzB > mz.max() - margin:
            continue

        areaA, blA, nA = integrate_peak(mz, inten, mzA, args.half_width, args.baseline_flanks)
        areaB, blB, nB = integrate_peak(mz, inten, mzB, args.half_width, args.baseline_flanks)

        z_results.append(dict(z=z, mzA=mzA, mzB=mzB,
                              areaA=areaA, areaB=areaB,
                              blA=blA, blB=blB, nA=nA, nB=nB))

    if not z_results:
        return None

    total_A = sum(r['areaA'] for r in z_results)
    total_B = sum(r['areaB'] for r in z_results)

    if total_B <= 0:
        return None

    return dict(
        file=os.path.basename(path),
        ratio=total_A / total_B,
        total_A=total_A,
        total_B=total_B,
        z_results=z_results,
        n_z=len(z_results),
    )


# ---------------------------------------------------------------------------
# Multi-file run
# ---------------------------------------------------------------------------

def run(files, args):
    results = []
    for path in files:
        fname = os.path.basename(path)
        try:
            res = quantitate_one_file(path, args)
            if res is None:
                print(f"  [WARN] {fname}: no usable charge states found")
                continue
            results.append(res)
            print(f"  {fname:<40}  ratio={res['ratio']:.4f}  "
                  f"(sumA={res['total_A']:.1f}  sumB={res['total_B']:.1f}  n_z={res['n_z']})")
        except Exception as exc:
            print(f"  [FAIL] {fname}: {exc}")
            if args.verbose:
                import traceback; traceback.print_exc()

    return results


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_outputs(results, args):
    os.makedirs(args.out, exist_ok=True)
    if not results:
        print("No results to write."); return

    ratios = np.array([r['ratio'] for r in results])
    mean_r = float(np.mean(ratios))
    std_r  = float(np.std(ratios, ddof=1)) if len(ratios) > 1 else 0.0
    cv     = 100.0 * std_r / mean_r if mean_r else float('nan')

    # Per-file CSV
    perfile = os.path.join(args.out, "direct_per_file.csv")
    with open(perfile, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "ratio_A_over_B", "total_area_A", "total_area_B", "n_z_used"])
        for r in results:
            w.writerow([r['file'], f"{r['ratio']:.5f}",
                        f"{r['total_A']:.3f}", f"{r['total_B']:.3f}", r['n_z']])

    # Per-z profile averaged across files (diagnostic)
    zprofile = os.path.join(args.out, "direct_zprofile.csv")
    all_zs = sorted({zr['z'] for r in results for zr in r['z_results']})
    with open(zprofile, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["z", "mean_areaA", "std_areaA", "mean_areaB", "std_areaB",
                    "mean_ratio_z", "n_files"])
        for z in all_zs:
            aAs = [zr['areaA'] for r in results for zr in r['z_results'] if zr['z']==z]
            aBs = [zr['areaB'] for r in results for zr in r['z_results'] if zr['z']==z]
            if not aAs: continue
            aAs, aBs = np.array(aAs), np.array(aBs)
            rz = aAs / np.where(aBs > 0, aBs, np.nan)
            w.writerow([z, f"{aAs.mean():.3f}", f"{aAs.std():.3f}",
                        f"{aBs.mean():.3f}", f"{aBs.std():.3f}",
                        f"{np.nanmean(rz):.4f}", len(aAs)])

    # Detailed per-z per-file CSV (for full audit trail)
    detail = os.path.join(args.out, "direct_per_z_detail.csv")
    with open(detail, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file","z","mzA","mzB","areaA","areaB","ratio_z",
                    "blA","blB","nA","nB"])
        for r in results:
            for zr in r['z_results']:
                rz = zr['areaA']/zr['areaB'] if zr['areaB']>0 else float('nan')
                w.writerow([r['file'], zr['z'],
                            f"{zr['mzA']:.4f}", f"{zr['mzB']:.4f}",
                            f"{zr['areaA']:.4f}", f"{zr['areaB']:.4f}",
                            f"{rz:.4f}" if np.isfinite(rz) else "NA",
                            f"{zr['blA']:.4f}", f"{zr['blB']:.4f}",
                            zr['nA'], zr['nB']])

    # Summary text
    summary = os.path.join(args.out, "direct_summary.txt")
    with open(summary, "w") as fh:
        fh.write("Direct m/z integration — quantitation summary\n")
        fh.write("=" * 52 + "\n")
        fh.write(f"Species A expected mass : {args.massA} Da\n")
        fh.write(f"Species B expected mass : {args.massB} Da\n")
        fh.write(f"Charge range used       : z={args.zmin}–{args.zmax}\n")
        fh.write(f"Integration half-width  : ±{args.half_width} m/z\n")
        fh.write(f"Baseline flank width    : {args.baseline_flanks} m/z\n")
        fh.write(f"N replicates            : {len(results)}\n\n")
        fh.write(f"Mean A/B ratio          : {mean_r:.4f}\n")
        fh.write(f"Std dev                 : {std_r:.4f}\n")
        fh.write(f"%CV                     : {cv:.2f}%\n\n")
        fh.write("Per-replicate ratios:\n")
        for r in results:
            fh.write(f"  {r['file']:<40}  {r['ratio']:.4f}\n")

    # Plot: per-file ratios + per-z profile
    if _HAVE_MPL and len(results) > 1:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

        # Left: per-file ratios
        fnames = [r['file'].replace(" Scan[1].raw","") for r in results]
        rat    = [r['ratio'] for r in results]
        ax1.bar(range(len(rat)), rat, color='#185FA5', alpha=0.8)
        ax1.axhline(mean_r, color='crimson', ls='--', lw=1.5, label=f"mean {mean_r:.3f}")
        ax1.set_xticks(range(len(fnames)))
        ax1.set_xticklabels(fnames, rotation=40, ha='right', fontsize=8)
        ax1.set_ylabel("A/B ratio")
        ax1.set_title(f"Per-replicate ratio  (CV={cv:.1f}%)")
        ax1.legend(fontsize=8)

        # Right: mean per-z areas (charge envelope shape)
        z_vals, mA_vals, mB_vals = [], [], []
        for z in all_zs:
            aAs = [zr['areaA'] for r in results for zr in r['z_results'] if zr['z']==z]
            aBs = [zr['areaB'] for r in results for zr in r['z_results'] if zr['z']==z]
            if aAs:
                z_vals.append(z)
                mA_vals.append(np.mean(aAs))
                mB_vals.append(np.mean(aBs))
        ax2.plot(z_vals, mA_vals, 'o-', color='#185FA5', label='Species A', ms=4)
        ax2.plot(z_vals, mB_vals, 's--', color='#1D9E75', label='Species B', ms=4)
        ax2.set_xlabel("charge state z")
        ax2.set_ylabel("mean integrated area")
        ax2.set_title("Charge envelope (mean across replicates)")
        ax2.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "direct_plots.png"), dpi=140)
        plt.close(fig)

    print(f"\nResults ({len(results)} files):")
    print(f"  Mean ratio : {mean_r:.4f}")
    print(f"  %CV        : {cv:.2f}%")
    print(f"\nOutputs in {args.out}\\")
    print(f"  direct_per_file.csv")
    print(f"  direct_per_z_detail.csv")
    print(f"  direct_zprofile.csv")
    print(f"  direct_summary.txt")
    if _HAVE_MPL:
        print(f"  direct_plots.png")
    print()
    print("Key diagnostics to check in direct_per_z_detail.csv:")
    print("  1. Per-z ratio should be consistent across z=15-22; spikes outside")
    print("     that range confirm the safe-zone choice is appropriate.")
    print("  2. areaB should follow a smooth bell curve across z; erratic values")
    print("     at specific z indicate interference (adducts, co-eluters).")
    print("  3. blA and blB should be low relative to areaA/areaB; high baselines")
    print("     indicate congested m/z space at that charge state.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def list_files(data_dir):
    exts = (".raw", ".mzml", ".mzxml")
    files = sorted(
        os.path.join(data_dir, f)
        for f in os.listdir(data_dir)
        if f.lower().endswith(exts) and os.path.isfile(os.path.join(data_dir, f))
    )
    if not files:
        raise SystemExit(f"No .raw/.mzML files found in {data_dir!r}")
    return files


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Direct m/z-domain A/B ratio quantitation — no deconvolution.")
    ap.add_argument("--data",   required=True,
                    help="Folder of replicate .raw files")
    ap.add_argument("--out",    default="results_direct")
    ap.add_argument("--massA",  type=float, default=23412.0,
                    help="Expected zero-charge mass of species A (Da)")
    ap.add_argument("--massB",  type=float, default=23658.0,
                    help="Expected zero-charge mass of species B (Da)")
    ap.add_argument("--zmin",   type=int,   default=15,
                    help="Lowest charge state to include (default 15)")
    ap.add_argument("--zmax",   type=int,   default=22,
                    help="Highest charge state to include (default 22)")
    ap.add_argument("--half-width", type=float, default=1.5, dest="half_width",
                    help="Integration half-width in m/z around each expected peak (default 1.5)")
    ap.add_argument("--baseline-flanks", type=float, default=3.0, dest="baseline_flanks",
                    help="Width of baseline flanking regions in m/z (default 3.0)")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    files = list_files(args.data)
    print(f"Found {len(files)} files.")
    print(f"Integrating z={args.zmin}–{args.zmax}, window ±{args.half_width} m/z, "
          f"baseline flanks {args.baseline_flanks} m/z\n")
    results = run(files, args)
    write_outputs(results, args)


if __name__ == "__main__":
    main()
