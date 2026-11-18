# -*- coding: utf-8 -*-
"""
merge_and_plot_hsp90_quartiles.py
---------------------------------
Merge AUC/Volume with GeneScore and make:
  - All cells plot (color = volume, shared color scale)
  - Four extra plots by volume quartiles:
        Q1: 0–25% (smallest volumes)
        Q2: 25–50%
        Q3: 50–75%
        Q4: 75–100% (largest volumes)
"""

import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- Config ----------
AUC_VOL_XLSX   = "hsp90_results.xlsx"
GENESCORE_XLSX = "HSP90AA1_geneScores_proj_d7.xlsx"

OUT_MERGED_XLSX = "hsp90_auc_volume_genescore_merged.xlsx"
OUT_SCATTER_ALL = "scatter_auc_vs_genescore_all.png"
OUT_SCATTER_Q1  = "scatter_auc_vs_genescore_vol_q0_25.png"
OUT_SCATTER_Q2  = "scatter_auc_vs_genescore_vol_q25_50.png"
OUT_SCATTER_Q3  = "scatter_auc_vs_genescore_vol_q50_75.png"
OUT_SCATTER_Q4  = "scatter_auc_vs_genescore_vol_q75_100.png"

POINT_SIZE = 40
COLORMAP   = "viridis"

# ---------- Helpers ----------
lane_pat = re.compile(r"r\s*(\d+)[_\-]?\s*c\s*(\d+)", re.IGNORECASE)

def to_lane_id(s: str) -> str | None:
    if not isinstance(s, str):
        s = str(s)
    m = lane_pat.search(s)
    if not m:
        return None
    r, c = m.group(1), m.group(2)
    return f"r{int(r)}c{int(c)}"

def read_auc_vol_table(path: str) -> pd.DataFrame:
    xls = pd.ExcelFile(path)
    sheet = "HSP90" if "HSP90" in xls.sheet_names else ("All" if "All" in xls.sheet_names else xls.sheet_names[0])
    df = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")

    cols = {c.lower(): c for c in df.columns}
    lane_col = cols.get("lane_id") or "lane_id"
    auc_col  = cols.get("auc")     or "auc"
    vol_col  = cols.get("volume")  or "volume"

    if "qc_pass" in df.columns and sheet.lower() == "all":
        df = df[df["qc_pass"] == True]

    out = df[[lane_col, auc_col, vol_col]].copy()
    out.columns = ["lane_id", "auc", "volume"]
    out["lane_id"] = out["lane_id"].astype(str).str.strip().str.lower()
    out["auc"]     = pd.to_numeric(out["auc"], errors="coerce")
    out["volume"]  = pd.to_numeric(out["volume"], errors="coerce")
    return out.dropna(subset=["lane_id", "auc", "volume"])

def read_gene_scores(path: str) -> pd.DataFrame:
    df = pd.read_excel(path, engine="openpyxl")
    lc = {c.lower(): c for c in df.columns}
    barcode_col   = lc.get("cellbarcode") or lc.get("barcode") or df.columns[0]
    genescore_col = lc.get("genescore") or (df.columns[2] if df.shape[1] >= 3 else df.columns[-1])

    out = df[[barcode_col, genescore_col]].copy()
    out.columns = ["cellBarcode", "GeneScore"]
    out["lane_id"] = out["cellBarcode"].astype(str).map(to_lane_id)
    out = out.dropna(subset=["lane_id"])
    out["lane_id"]   = out["lane_id"].str.lower()
    out["GeneScore"] = pd.to_numeric(out["GeneScore"], errors="coerce")
    return out.dropna(subset=["GeneScore"])[["lane_id", "GeneScore"]]

def plot_scatter(df, out_png, title_prefix, xlim, ylim, vmin, vmax):
    """Generate a scatter plot for a given subset with shared color scale."""
    x = df["auc"].to_numpy(float)
    y = df["GeneScore"].to_numpy(float)
    z = df["volume"].to_numpy(float)

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(x, y, c=z, s=POINT_SIZE, cmap=COLORMAP,
                     vmin=vmin, vmax=vmax, edgecolor="none")
    cb = plt.colorbar(sc)
    cb.set_label("Cell volume (µm³)")

    plt.xlabel("AUC (HSP90)")
    plt.ylabel("GeneScore (HSP90AA1)")
    plt.title(f"{title_prefix}  N={len(df)}")

    # Label high AUC or high GeneScore
    for _, row in df.iterrows():
        if row["auc"] > 50 or row["GeneScore"] > 3:
            plt.text(row["auc"] + 1, row["GeneScore"], row["lane_id"],
                     fontsize=8, color="black", ha="left", va="center")

    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

# ---------- Main ----------
def main():
    if not os.path.exists(AUC_VOL_XLSX):
        raise FileNotFoundError(f"Missing file: {AUC_VOL_XLSX}")
    if not os.path.exists(GENESCORE_XLSX):
        raise FileNotFoundError(f"Missing file: {GENESCORE_XLSX}")

    df_aucvol = read_auc_vol_table(AUC_VOL_XLSX)
    df_gene   = read_gene_scores(GENESCORE_XLSX)

    merged = df_aucvol.merge(df_gene, on="lane_id", how="inner")
    merged = merged.replace([np.inf, -np.inf], np.nan).dropna(subset=["auc", "GeneScore", "volume"])
    merged = merged[merged["GeneScore"] > 0].copy()
    merged = merged[merged["auc"] <= 200].copy()

    # Save merged table
    with pd.ExcelWriter(OUT_MERGED_XLSX, engine="xlsxwriter") as writer:
        merged[["lane_id", "auc", "GeneScore", "volume"]].to_excel(writer, sheet_name="merged", index=False)

    # Global color and axis limits (shared across all plots)
    vmin, vmax = merged["volume"].min(), merged["volume"].max()
    xlim = (np.nanmin(merged["auc"]) * 0.95,  np.nanmax(merged["auc"]) * 1.05)
    ylim = (np.nanmin(merged["GeneScore"]) * 0.95, np.nanmax(merged["GeneScore"]) * 1.05)

    # ---- All cells ----
    plot_scatter(merged, OUT_SCATTER_ALL, "All Cells", xlim, ylim, vmin, vmax)

    # ---- Quartiles by volume ----
    q25, q50, q75 = merged["volume"].quantile([0.25, 0.50, 0.75]).to_list()

    # Q1: smallest 25% (v <= q25)
    df_q1 = merged[merged["volume"] <= q25]
    # Q2: 25–50% (q25 < v <= q50)
    df_q2 = merged[(merged["volume"] > q25) & (merged["volume"] <= q50)]
    # Q3: 50–75% (q50 < v <= q75)
    df_q3 = merged[(merged["volume"] > q50) & (merged["volume"] <= q75)]
    # Q4: largest 25% (v > q75)
    df_q4 = merged[merged["volume"] > q75]

    plot_scatter(df_q1, OUT_SCATTER_Q1, f"Volume Quartile Q1 (≤ {q25:.0f} µm³)", xlim, ylim, vmin, vmax)
    plot_scatter(df_q2, OUT_SCATTER_Q2, f"Volume Quartile Q2 ({q25:.0f}–{q50:.0f} µm³)", xlim, ylim, vmin, vmax)
    plot_scatter(df_q3, OUT_SCATTER_Q3, f"Volume Quartile Q3 ({q50:.0f}–{q75:.0f} µm³)", xlim, ylim, vmin, vmax)
    plot_scatter(df_q4, OUT_SCATTER_Q4, f"Volume Quartile Q4 (> {q75:.0f} µm³)", xlim, ylim, vmin, vmax)

    print(f"Merged rows: {len(merged)}")
    print("Quartile cutoffs (µm³): Q1≤{q1:.0f}, Q2≤{q2:.0f}, Q3≤{q3:.0f}".format(
        q1=q25, q2=q50, q3=q75))
    print("Saved:", OUT_MERGED_XLSX)
    print("Saved:", OUT_SCATTER_ALL)
    print("Saved:", OUT_SCATTER_Q1)
    print("Saved:", OUT_SCATTER_Q2)
    print("Saved:", OUT_SCATTER_Q3)
    print("Saved:", OUT_SCATTER_Q4)

if __name__ == "__main__":
    main()
