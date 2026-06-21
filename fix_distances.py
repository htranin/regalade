"""
fix_distances.py — REGALADEv2 distance re-prioritization + Simbad flagging.

D priority:  2 > 4 > 7 > 5 > 6 > 8 > keep existing
z priority:  4 > 7 > 8 > keep existing  (cats 2, 5, 6 excluded)

Usage:
    python fix_distances.py [--full]

Without --full, operates on the _sample.fits file.
"""

import argparse
import os
import subprocess
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from tqdm import tqdm
import astropy.units as u

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SAMPLE_IN  = "regalade_within_1500Mpc_sample.fits"
FULL_IN    = "regalade_within_1500Mpc.fits"
SAMPLE_OUT = "regalade_v2_sample.fits"
FULL_OUT   = "regalade_v2.fits"

# (cat_index, fits_file, ra_col, dec_col, D_col, D_err_col, z_col)
# z_col=None means derive z from D using the approximation formula.
CATALOGS = [
    (2, "2_HECATEv2.fits",          "RA",  "DEC", "Dist",  "e_Dist", None),
    (4, "4_desi_dr1_v2026.fits",    "ra",  "dec", "D",     None,     "z"),
    (7, "7_NEDLVS_20250602_z.fits", "ra",  "dec", "D",     "D_err",  None),
    (5, "3_cosmicflows_input.fits", "ra",  "dec", "D",     "Derr",   None),
    (6, "6_NEDLVS_20250602_D.fits", "ra",  "dec", "D",     "D_err",  None),
    (8, "8_NEDLVS_20250602_rest.fits","ra","dec", "D",     "D_err",  None),
]

D_PRIORITY = [2, 4, 7, 5, 6, 8]
Z_PRIORITY = [4, 7, 8]

LARGE_GALAXY_R1   = 60.0   # arcsec
RADIUS_LARGE      = 10.0   # arcsec
RADIUS_SMALL      = 5.0    # arcsec
OFFSET_WARN       = 5.0    # arcsec (flag in match_offset)

SIMBAD_MAX_D      = 200.0  # Mpc
SIMBAD_RADIUS     = 5.0    # arcsec
SIMBAD_Z_THRESH   = 0.10   # fractional term in discrepancy criterion
SIMBAD_Z_TOL      = 0.001  # absolute tolerance floor: flag if |z - z_simbad| > 0.1*z + SIMBAD_Z_TOL
STILTS_JAR        = "C:/Users/hgtra/OneDrive/Documents/stilts.jar"
SIMBAD_REF_IDX    = 16   # ref_z_in value meaning "z taken from Simbad"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def z_from_d(d):
    """Approximate z from comoving distance d (Mpc)."""
    return 2.3485e-4 * d + 1.69e-5 * d**1.348 - 1.336e-5 * d**1.403


def load_catalog(fname, ra_col, dec_col, D_col, D_err_col, z_col):
    """Return (SkyCoord, D_arr, D_err_arr, z_arr) for a catalog file."""
    with fits.open(fname) as hdul:
        data = hdul[1].data
        ra  = data[ra_col].astype(float)
        dec = data[dec_col].astype(float)
        D   = data[D_col].astype(float)
        D_err = data[D_err_col].astype(float) if D_err_col else np.full(len(D), np.nan)
        z   = data[z_col].astype(float) if z_col else None
    # Mask out bad (non-positive) distances
    valid = D > 0
    ra, dec, D, D_err = ra[valid], dec[valid], D[valid], D_err[valid]
    if z is not None:
        z = z[valid]
    sc = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    return sc, D, D_err, z


def match_galaxies(regal_sc, r1_arcsec, cat_sc):
    """
    Return (matched_idx, sep_arcsec, within_mask) where within_mask is True
    for galaxies whose nearest neighbor in cat_sc is within the allowed radius.
    Radius depends on R1: large galaxies (R1 > 60") use 10", others 5".
    """
    idx, sep2d, _ = regal_sc.match_to_catalog_sky(cat_sc)
    sep_arcsec = sep2d.arcsec
    radius = np.where(r1_arcsec > LARGE_GALAXY_R1, RADIUS_LARGE, RADIUS_SMALL)
    within = sep_arcsec <= radius
    return idx, sep_arcsec, within


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(use_full=False):
    in_file  = FULL_IN  if use_full else SAMPLE_IN
    out_file = FULL_OUT if use_full else SAMPLE_OUT

    print(f"Input:  {in_file}")
    print(f"Output: {out_file}")

    # --- Load REGALADE ---
    with fits.open(in_file) as hdul:
        hdr = hdul[1].header.copy()
        orig = hdul[1].data

    n = len(orig)
    print(f"REGALADE rows: {n:,}")

    regal_sc = SkyCoord(ra=orig["gal_ra"] * u.deg, dec=orig["gal_dec"] * u.deg)
    r1_arcsec = orig["R1"].astype(float)

    # Working copies of columns to modify
    D_orig      = orig["D"].astype(float).copy()   # keep original for magnitude/mass correction
    D_out       = orig["D"].astype(float).copy()
    D_input_out = orig["D_input"].astype(float).copy()
    D_err_out   = orig["D_err"].astype(float).copy()
    ref_D_out   = orig["ref_D_in"].astype(np.int32).copy()
    z_out       = orig["z"].astype(float).copy()
    ref_z_out   = orig["ref_D_in"].astype(np.int32).copy()  # default to same source as D
    match_offset_out = np.full(n, np.nan, dtype=float)

    # --- Load all catalogs and pre-compute matches ---
    cat_data = {}   # cat_index -> (sc, D, D_err, z, matched_idx, sep, within)
    for (cidx, fname, ra_col, dec_col, D_col, D_err_col, z_col) in CATALOGS:
        print(f"  Loading catalog {cidx}: {fname} ...", end=" ", flush=True)
        sc, D, D_err, z = load_catalog(fname, ra_col, dec_col, D_col, D_err_col, z_col)
        idx, sep, within = match_galaxies(regal_sc, r1_arcsec, sc)
        print(f"{within.sum():,} matches")
        cat_data[cidx] = (sc, D, D_err, z, idx, sep, within)

    # --- Apply D priority ---
    D_assigned = np.zeros(n, dtype=bool)
    for cidx in D_PRIORITY:
        sc, D, D_err, z, idx, sep, within = cat_data[cidx]
        # Only update galaxies not yet assigned
        update_mask = within & ~D_assigned
        if not update_mask.any():
            continue
        cat_idx = idx[update_mask]
        D_out[update_mask]       = D[cat_idx]
        D_input_out[update_mask] = D[cat_idx]
        D_err_out[update_mask]   = D_err[cat_idx]
        ref_D_out[update_mask]   = cidx
        match_offset_out[update_mask] = sep[update_mask]
        D_assigned[update_mask]  = True
        print(f"  D priority {cidx}: assigned {update_mask.sum():,} galaxies")

    # --- Catastrophic outlier guard: galaxy too small for its distance ---
    # If D < 100 / R1 (R1 in arcsec, D in Mpc), the implied physical size is
    # unrealistically small; fall back to D_tmean.
    D_tmean = orig["D_tmean"].astype(float)
    outlier_mask = (D_out > 0) & (r1_arcsec > 0) & (D_out < 100.0 / r1_arcsec)
    D_out[outlier_mask] = D_tmean[outlier_mask]
    print(f"  Outlier guard (D < 100/R1): {outlier_mask.sum():,} galaxies replaced with D_tmean")

    # --- Update absgmag and logM for all galaxies where D changed ---
    absgmag_out = orig["absgmag"].astype(float).copy()
    logM_out    = orig["logM"].astype(float).copy()
    changed = (D_out > 0) & (D_orig > 0) & (D_out != D_orig)
    log_ratio = np.zeros(n, dtype=float)
    log_ratio[changed] = np.log10(D_out[changed] / D_orig[changed])
    absgmag_out[changed] -= 5.0 * log_ratio[changed]
    logM_out[changed]    += 2.0 * log_ratio[changed]
    print(f"  absgmag/logM updated: {changed.sum():,} galaxies")

    # --- Apply z priority ---
    z_assigned = np.zeros(n, dtype=bool)
    for cidx in Z_PRIORITY:
        sc, D, D_err, z_cat, idx, sep, within = cat_data[cidx]
        update_mask = within & ~z_assigned
        if not update_mask.any():
            continue
        cat_idx = idx[update_mask]
        if z_cat is not None:
            z_out[update_mask] = z_cat[cat_idx]
        else:
            z_out[update_mask] = z_from_d(D[cat_idx])
        ref_z_out[update_mask] = cidx
        z_assigned[update_mask] = True
        print(f"  z priority {cidx}: assigned {update_mask.sum():,} galaxies")

    # Galaxies whose z source is a distance-only catalog have no valid redshift
    no_z_mask = np.isin(ref_z_out, [5, 6])
    z_out[no_z_mask] = np.nan
    print(f"  z set to NaN (cat 5/6, all distances): {no_z_mask.sum():,}")

    # --- Simbad matching via STILTS cdsskymatch ---
    print(f"\nQuerying Simbad for galaxies with D < {SIMBAD_MAX_D} Mpc ...")
    simbad_z_out    = np.full(n, np.nan, dtype=float)
    simbad_flag_out = np.zeros(n, dtype=bool)

    near_mask = D_out < SIMBAD_MAX_D
    near_idx  = np.where(near_mask)[0]
    print(f"  {len(near_idx):,} galaxies within {SIMBAD_MAX_D} Mpc")


    nearby_file  = "nearby.fits"
    matches_file = "simbad_matches.fits"

    # Write subset with a row index so we can merge back
    row_idx = near_idx.astype(np.int32)
    nearby_hdu = fits.BinTableHDU.from_columns([
        fits.Column(name="row_idx", format="J",  array=row_idx),
        fits.Column(name="gal_ra",  format="D",  array=orig["gal_ra"][near_idx].astype(float)),
        fits.Column(name="gal_dec", format="D",  array=orig["gal_dec"][near_idx].astype(float)),
    ])
    nearby_hdu.writeto(nearby_file,overwrite=True)

    cmd = [
        "java", "-jar", STILTS_JAR, "cdsskymatch",
        f"in={nearby_file}",
        "ra=gal_ra", "dec=gal_dec",
        f"radius={SIMBAD_RADIUS}",
        "find=all",
        "cdstable=SIMBAD",
        f"out={matches_file}",
    ]
    print(f"  Running: {' '.join(cmd)}")
    ret = subprocess.run(cmd, capture_output=True, text=True)
    if ret.returncode != 0:
        print(f"  stilts error:\n{ret.stderr}")
    elif not os.path.exists(matches_file):
        print("  No Simbad matches file produced.")
    else:
        with fits.open(matches_file) as hdul:
            m = hdul[1].data
            print(f"  Simbad returned {len(m):,} matched rows")
            print(f"  Columns: {hdul[1].columns.names}")
        with fits.open(matches_file) as hdul:
            m = hdul[1].data
            cols_lower = {c.lower(): c for c in m.dtype.names}
            # STILTS cdsskymatch returns Simbad's redshift as 'redshift'
            z_col = cols_lower.get("redshift") or cols_lower.get("rvz_redshift") or cols_lower.get("z_value")
            if z_col is None:
                print(f"  WARNING: no redshift column found in Simbad output. Columns: {list(m.dtype.names)}")
            else:
                # Group by row_idx, keep match with highest nbref
                best = {}  # row_idx -> row with max nbref that has a valid z
                for row in m:
                    gi      = int(row["row_idx"])
                    z_s     = float(row[z_col])
                    nbref   = int(row["nbref"]) if not np.ma.is_masked(row["nbref"]) else 0
                    if np.isnan(z_s):
                        continue
                    if gi not in best or nbref > best[gi][1]:
                        best[gi] = (z_s, nbref)

                # Pass 1: store Simbad z; replace z for distance-only catalogs (5, 6)
                n_simbad_replaced = 0
                for gi, (z_simbad, _) in best.items():
                    simbad_z_out[gi] = z_simbad
                    if ref_z_out[gi] in (5, 6) and z_simbad > 0:
                        z_out[gi]     = z_simbad
                        ref_z_out[gi] = SIMBAD_REF_IDX
                        n_simbad_replaced += 1
                print(f"  z replaced from Simbad (was cat 5/6): {n_simbad_replaced}")

                # Pass 2: flag remaining discrepancies
                for gi, (z_simbad, _) in best.items():
                    z_g = z_out[gi]
                    if z_simbad > 0 and not np.isnan(z_g):
                        if abs(z_g - z_simbad) > SIMBAD_Z_THRESH * z_g + SIMBAD_Z_TOL:
                            simbad_flag_out[gi] = True

    n_flagged = simbad_flag_out.sum()
    print(f"  Simbad flagged: {n_flagged} galaxies with |Δz - z_simbad| > {SIMBAD_Z_THRESH*100:.0f}%·z + {SIMBAD_Z_TOL}")

    # --- Build output FITS ---
    print("\nBuilding output FITS ...")

    # Columns that have been updated
    overrides = {
        "D":       (D_out.astype(np.float32),       "E"),
        "D_input": (D_input_out.astype(np.float32), "E"),
        "D_err":   (D_err_out.astype(np.float32),   "E"),
        "ref_D_in":(ref_D_out,                      "J"),
        "z":       (z_out.astype(np.float32),        "E"),
        "absgmag": (absgmag_out.astype(np.float32), "E"),
        "logM":    (logM_out.astype(np.float32),     "E"),
    }

    def fits_format(dtype, name, nbytes):
        if np.issubdtype(dtype, np.floating):
            return "E"
        if dtype == np.dtype("bool"):
            return "L"
        if dtype == np.dtype("int8") or dtype == np.dtype("uint8"):
            return "B"
        if dtype == np.dtype("int16"):
            return "I"
        if np.issubdtype(dtype, np.integer):
            return "J"
        if dtype.kind == "S":              # byte string
            return f"{dtype.itemsize}A"
        if dtype.kind == "U":              # unicode string
            return f"{dtype.itemsize // 4}A"
        raise TypeError(f"Unsupported dtype {dtype} for column '{name}'")

    # Pre-compute format for every original column once
    col_formats = {}
    for col_def in orig.columns:
        cname = col_def.name
        if cname in overrides:
            col_formats[cname] = overrides[cname][1]
        else:
            col_formats[cname] = fits_format(orig[cname].dtype, cname, orig[cname].dtype.itemsize)

    # Split into N_BINS RA bins, write each, then concatenate with STILTS tcat
    N_BINS   = 10
    ra_all   = orig["gal_ra"].astype(float)
    bin_edges = np.linspace(0.0, 360.0, N_BINS + 1)
    bin_files = []

    for b in tqdm(range(N_BINS)):
        mask = (ra_all >= bin_edges[b]) & (ra_all < bin_edges[b + 1])
        if b == N_BINS - 1:          # include RA=360 in last bin
            mask = (ra_all >= bin_edges[b])
        idx  = np.where(mask)[0]
        if len(idx) == 0:
            continue
        bin_file = out_file.replace(".fits", f"_bin{b:02d}.fits")
        bin_files.append(bin_file)
        print(f"  Writing bin {b}: RA [{bin_edges[b]:.1f}, {bin_edges[b+1]:.1f}) — {len(idx):,} rows -> {bin_file}")

        bin_cols = []
        for col_def in orig.columns:
            cname = col_def.name
            if cname in overrides:
                arr = overrides[cname][0][idx]
            else:
                arr = orig[cname][idx]
            bin_cols.append(fits.Column(name=cname, format=col_formats[cname], array=arr))
        bin_cols += [
            fits.Column(name="ref_z_in",              format="J", array=ref_z_out[idx]),
            fits.Column(name="match_offset",          format="E", array=match_offset_out[idx].astype(np.float32), unit="arcsec"),
            fits.Column(name="simbad_z",              format="E", array=simbad_z_out[idx].astype(np.float32)),
            fits.Column(name="f_simbad_zdiscrepancy", format="L", array=simbad_flag_out[idx]),
        ]
        hdu = fits.BinTableHDU.from_columns(fits.ColDefs(bin_cols), header=hdr)
        hdu.writeto(bin_file, overwrite=True)

    # Concatenate with STILTS tcat
    print(f"\nConcatenating {len(bin_files)} bins with STILTS tcat ...")
    cmd = (
        ["java", "-jar", STILTS_JAR, "tcat"]
        + [f"in={f}" for f in bin_files]
        + ["ifmt=fits", f"out={out_file}", "ofmt=fits"]
    )
    ret = subprocess.run(cmd, capture_output=True, text=True)
    if ret.returncode != 0:
        print(f"  stilts tcat error:\n{ret.stderr}")
    else:
        for f in bin_files:
            os.remove(f)
        print(f"Saved: {out_file}")

    # --- Summary ---
    print("\n=== Summary ===")
    print(f"D updated:  {D_assigned.sum():,} / {n:,} galaxies")
    print(f"z updated:  {z_assigned.sum():,} / {n:,} galaxies (+ {(ref_z_out == SIMBAD_REF_IDX).sum():,} replaced from Simbad)")
    ref_D_vals, ref_D_counts = np.unique(ref_D_out, return_counts=True)
    print("ref_D_in distribution (after):")
    for v, c in zip(ref_D_vals, ref_D_counts):
        print(f"  cat {v:2d}: {c:,}")
    ref_z_vals, ref_z_counts = np.unique(ref_z_out[ref_z_out >= 0], return_counts=True)
    print("ref_z_in distribution (newly assigned):")
    for v, c in zip(ref_z_vals, ref_z_counts):
        print(f"  cat {v:2d}: {c:,}")
    print(f"Simbad matches:  {(~np.isnan(simbad_z_out)).sum():,}")
    print(f"Simbad flagged:  {n_flagged:,}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true", help="Process full catalog instead of sample")
    args = parser.parse_args()
    main(use_full=args.full)
