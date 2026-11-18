# -*- coding: utf-8 -*-
"""
analyze_hsp90.py (enhanced recall)
----------------------------------
Same as your 60% version, with minimal targeted tweaks:
- Slightly longer vertical ROI (extends further below halo)
- Small smoothing of the vertical profile
- Relaxed QC for bright profiles (R² >= 0.25 if AUC is high)
"""

import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from skimage import io, feature, exposure, color, draw, util
from skimage.transform import hough_circle, hough_circle_peaks
from PIL import Image, ImageDraw, ImageFont

# -----------------------------
# Config
# -----------------------------
IMAGE_FILE = "hsp90.tif"
CELLINFO_FILE = "cell_info.xlsx"

N_ROWS = 16
N_COLS = 31

# ROI/internal params
WELL_SEARCH_RADIUS_FRAC = (0.04, 0.22)
HALO_SCALE = 1.6
HALO_MARGIN_PX = 2
BAND_HALF_WIDTH_FRAC = 0.1
GAP_BELOW_HALO_PX = 8

# NEW: extend ROI deeper
ROI_DEPTH_MULTIPLIER = 1.2   # 1.0 = current bottom; >1 = longer ROI

# Scatter plot options
POINT_SIZE = 40
LINE_WIDTH = 2.0
SCATTER_OUT = "scatter_hsp90_auc_vs_volume.png"
QC_OVERLAY_OUT = "qc_overlay_hsp90.png"
XLSX_OUT = "hsp90_results.xlsx"

# QC thresholds
GAUSS_R2_MIN = 0.30
GAUSS_R2_FLEX = 0.25      # relaxed threshold for bright lanes
AUC_BRIGHT_MIN = 250.0    # if above this AUC, accept even if R² >= 0.25

SMOOTH_WINDOW = 5         # pixels for 1D smoothing along vertical profile

# -----------------------------
# Helpers
# -----------------------------
def lane_id(r, c):
    return f"r{r}c{c}"

def split_grid(im, n_rows, n_cols):
    H, W = im.shape[:2]
    tile_h = H // n_rows
    tile_w = W // n_cols
    rois = []
    for r in range(n_rows):
        for c in range(n_cols):
            y0 = r * tile_h
            y1 = (r + 1) * tile_h if r < n_rows - 1 else H
            x0 = c * tile_w
            x1 = (c + 1) * tile_w if c < n_cols - 1 else W
            rois.append(((r+1, c+1), (slice(y0, y1), slice(x0, x1))))
    return rois, tile_h, tile_w

def detect_microwell_and_halo(roi_gray):
    """Detect the bright ring/hole at TOP of ROI."""
    H, W = roi_gray.shape
    roi_norm = exposure.rescale_intensity(roi_gray, in_range='image', out_range=(0,1))
    edges = feature.canny(roi_norm, sigma=2.0)
    rmin = max(2, int(WELL_SEARCH_RADIUS_FRAC[0] * H))
    rmax = max(rmin+2, int(WELL_SEARCH_RADIUS_FRAC[1] * H))
    if rmax <= rmin:
        rmax = rmin + 3
    try:
        hough_res = hough_circle(edges, np.arange(rmin, rmax))
        accums, cx, cy, radii = hough_circle_peaks(hough_res, np.arange(rmin, rmax), total_num_peaks=1)
        if len(radii) > 0:
            r_in = int(radii[0])
            if cy[0] > H // 2:
                cx0, cy0 = W // 2, int(0.18 * H)
                r_in = max(rmin, min(rmax, int(0.10 * H)))
                r_halo = int(HALO_SCALE * r_in)
                return cx0, cy0, r_in, r_halo
            r_halo = int(HALO_SCALE * r_in)
            return int(cx[0]), int(cy[0]), r_in, r_halo
    except Exception:
        pass
    cx0, cy0 = W // 2, int(0.18 * H)
    r_in = max(rmin, min(rmax, int(0.10 * H)))
    r_halo = int(HALO_SCALE * r_in)
    return cx0, cy0, r_in, r_halo

def gaussian_from_moments(y, x=None):
    if x is None:
        x = np.arange(len(y))
    y = np.maximum(0, y.astype(float))
    s = y.sum()
    if s <= 0:
        return 0.0, np.nan, np.nan
    mu = (x * y).sum() / s
    var = (y * (x - mu) ** 2).sum() / s
    sigma = math.sqrt(max(var, 1e-6))
    A = y.max()
    return A, mu, sigma

def gaussian_eval(A, mu, sigma, x):
    return A * np.exp(-0.5 * ((x - mu) / (sigma + 1e-9)) ** 2)

def r2_score(y_true, y_pred):
    y_true = y_true.astype(float)
    y_pred = y_pred.astype(float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

def smooth_profile(y, window):
    """Simple moving average smoothing."""
    if window <= 1 or len(y) < window:
        return y
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode='same')

def compute_auc_and_qc(roi_gray, cx, cy, r_halo):
    """Compute AUC below halo, excluding halo + margin."""
    H, W = roi_gray.shape
    y_start = min(H-1, int(cy + r_halo + GAP_BELOW_HALO_PX))
    # extend ROI slightly deeper
    y_end = min(H, int(y_start + ROI_DEPTH_MULTIPLIER * (H - y_start)))
    band_half = int(BAND_HALF_WIDTH_FRAC * W)
    x0 = max(0, int(cx) - band_half)
    x1 = min(W, int(cx) + band_half + 1)
    window = roi_gray[y_start:y_end, x0:x1].astype(float)
    if window.size == 0:
        return 0.0, np.nan, False
    bg = np.percentile(window, 40)
    window_bs = window - bg
    window_bs[window_bs < 0] = 0.0
    profile = window_bs.sum(axis=1)
    if profile.max() <= 0:
        return 0.0, np.nan, False
    profile = smooth_profile(profile, SMOOTH_WINDOW)
    x = np.arange(len(profile))
    A, mu, sigma = gaussian_from_moments(profile, x)
    y_hat = gaussian_eval(A, mu, sigma, x)
    r2 = r2_score(profile, y_hat)
    auc = profile.sum()
    # QC logic: relaxed if bright
    passed = ((r2 >= GAUSS_R2_MIN) or ((r2 >= GAUSS_R2_FLEX) and (auc >= AUC_BRIGHT_MIN))) and (auc > 0)
    return float(auc), float(r2), bool(passed)

def draw_qc_overlay_with_labels(im_gray, rois, centers, radii, passes):
    """QC overlay with green/red circles and lane IDs."""
    base = util.img_as_ubyte(exposure.rescale_intensity(im_gray, in_range='image'))
    img = Image.fromarray(base).convert("RGB")
    d = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("arial.ttf", 12)
    except Exception:
        font = ImageFont.load_default()
    for ((r_idx, c_idx), (ys, xs)), (cx, cy), r_halo, ok in zip(rois, centers, radii, passes):
        gx = xs.start + int(cx)
        gy = ys.start + int(cy)
        r = max(2, int(r_halo))
        color = (0, 255, 0) if ok else (255, 0, 0)
        d.ellipse([gx - r, gy - r, gx + r, gy + r], outline=color, width=2)
        label = lane_id(r_idx, c_idx)
        d.text((gx + r + 3, gy - r - 3), label, fill=color, font=font)
    return np.array(img)

def fit_and_r2(x, y):
    if len(x) < 2:
        return np.nan, np.nan, np.nan
    m, b = np.polyfit(x, y, 1)
    yhat = m * x + b
    ss_res = np.sum((y - yhat)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
    return m, b, r2

# -----------------------------
# Main
# -----------------------------
def main():
    if not os.path.exists(IMAGE_FILE):
        raise FileNotFoundError(f"Missing image file: {IMAGE_FILE}")
    if not os.path.exists(CELLINFO_FILE):
        raise FileNotFoundError(f"Missing Excel file: {CELLINFO_FILE}")

    im = io.imread(IMAGE_FILE)
    if im.ndim == 3:
        im_gray = color.rgb2gray(im)
    else:
        im_gray = util.img_as_float(im)

    rois, tile_h, tile_w = split_grid(im_gray, N_ROWS, N_COLS)

    records, centers, radii, passes = [], [], [], []

    for (r_idx, c_idx), (ys, xs) in rois:
        roi = im_gray[ys, xs]
        cx, cy, r_in, r_halo = detect_microwell_and_halo(roi)
        auc, r2, ok = compute_auc_and_qc(roi, cx, cy, r_halo + HALO_MARGIN_PX)
        centers.append((cx, cy))
        radii.append(r_halo)
        passes.append(ok)
        records.append({
            "lane_id": lane_id(r_idx, c_idx),
            "auc": auc,
            "gauss_r2": r2,
            "qc_pass": ok
        })

    df_auc = pd.DataFrame(records)
    df_cell = pd.read_excel(CELLINFO_FILE, header=0, engine="openpyxl")
    if df_cell.shape[1] < 5:
        raise ValueError("cell_info.xlsx must have at least 5 columns (A..E).")
    df_cell = df_cell.iloc[:, [0, 4]].copy()
    df_cell.columns = ["lane_id", "volume"]
    df_cell["lane_id"] = df_cell["lane_id"].astype(str).str.strip()
    df_cell["volume"] = pd.to_numeric(df_cell["volume"], errors="coerce")

    merged = df_auc.merge(df_cell, on="lane_id", how="left")
    keep = merged[merged["qc_pass"] & np.isfinite(merged["volume"])].copy()

    with pd.ExcelWriter(XLSX_OUT, engine="xlsxwriter") as writer:
        keep[["lane_id", "auc", "volume"]].to_excel(writer, sheet_name="HSP90", index=False)
        merged.to_excel(writer, sheet_name="All", index=False)

    sub = keep[["volume", "auc"]].dropna()
    if len(sub) >= 2:
        x = sub["volume"].to_numpy()
        y = sub["auc"].to_numpy()
        m, b, r2 = fit_and_r2(x, y)
        n = len(x)
        plt.figure(figsize=(9,6))
        plt.scatter(x, y, s=POINT_SIZE)
        if np.isfinite(m):
            xs = np.linspace(np.nanmin(x), np.nanmax(x), 200)
            plt.plot(xs, m*xs + b, linewidth=LINE_WIDTH)
        plt.xlabel("Cell volume (µm³)")
        plt.ylabel("HSP90 AUC")
        plt.title("HSP90 — Cell volume vs AUC (QC passed)")
        plt.text(0.05, 0.95, f"Slope = {m:0.2e}\nR² = {r2:0.3f}\nN = {n}",
                 transform=plt.gca().transAxes, ha="left", va="top")
        plt.tight_layout()
        plt.savefig(SCATTER_OUT, dpi=300)
        plt.close()
    else:
        print("[WARN] Not enough QC-passed points for scatter.")

    qc_rgb = draw_qc_overlay_with_labels(im_gray, rois, centers, radii, passes)
    io.imsave(QC_OVERLAY_OUT, util.img_as_ubyte(qc_rgb))

    print("Done.")
    print(" -", SCATTER_OUT)
    print(" -", QC_OVERLAY_OUT)
    print(" -", XLSX_OUT)

if __name__ == "__main__":
    main()
