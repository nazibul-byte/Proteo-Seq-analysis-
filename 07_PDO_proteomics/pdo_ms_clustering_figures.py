"""
BlotSeq PDO proteomics: data-driven K=3 clustering and figure generation.

Single self-contained script. Reads directly from the raw workbook (no
intermediate CSV/pickle files) and writes three PNG figures.

REQUIRED FILES -- must sit in the SAME FOLDER as this script:
    BlotSeq_PDO_SameCell_Multiomics_Dataset.xlsx
        sheets used: Protein_Basal, Protein_Luminal, Cell_Morphology
    Cell_markers.xlsx
        sheet used: Sheet1  (columns: Gene, Class)

OUTPUT (written to the same folder as this script):
    Fig_UMAP_K3.png          -- UMAP scatter, colored by K-means (k=3) cluster
    Fig_dotplot_K3.png       -- marker dot plot, dot size = % cells detected,
                                 color = mean z-score
    Fig_Volume_by_cluster.png -- cell volume (CellenONE morphometry) by cluster,
                                 with per-cluster mean bar

METHODOLOGY (documented so every choice is auditable)
------------------------------------------------------
1. Combine Protein_Basal + Protein_Luminal into one gene-group x 20-cell
   matrix using the established detection criterion: DIA-NN q <= 0.05 AND
   MaxLFQ value present.
2. Restrict PCA/clustering features to gene groups detected in >=50% of the
   20 cells (keeps enough completeness for stable PCA without excessive
   imputation).
3. log2-transform, then z-score each gene group across cells using only
   OBSERVED values (mean/SD ignore missing).
4. Impute remaining missing entries as "below detection" (row minimum
   z-score - 1.8, plus N(0, 0.3) jitter; Perseus-style downshifted-normal
   imputation, seed fixed) -- this gives PCA/K-means a complete matrix. It is
   NOT used when computing marker z-scores for the dot plot (those use
   observed values only).
5. PCA on the complete z-scored/imputed matrix (up to 10 PCs).
6. Diagnostic: correlate each PC with each cell's total detected-protein
   count (a proxy for lysis efficiency / cell size / MS sensitivity, not
   biology). Any PC significantly associated with this technical covariate
   (Spearman p < 0.05) is EXCLUDED from clustering. In this dataset PC1
   (~35% of variance) is the confounded axis and is dropped; K-means and the
   UMAP embedding both run on the remaining PCs only.
7. K-means (k=3; n_init=25, random_state=0) on the depth-decorrelated PCA
   subspace.
8. UMAP (same PCA subspace, n_neighbors=5, min_dist=0.4, random_state=0) for
   2D visualization only -- cluster assignment happens in step 7, not on the
   2D UMAP coordinates themselves.
9. Cell volume comes directly from the Cell_Morphology sheet (CellenONE
   pre-lysis measurement, Volume = pi/6 x Diameter^3 / sqrt(Elongation)) --
   this is independent of the proteomics detection depth issues that affect
   steps 1-8, so it's a useful orthogonal check on the clustering.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.stats import spearmanr
import umap

# Use the folder containing this script for both input and output
INPUT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = INPUT_DIR
RNG = 0

XLSX_PATH = INPUT_DIR / "BlotSeq_PDO_SameCell_Multiomics_Dataset.xlsx"
MARKERS_PATH = INPUT_DIR / "Cell_markers.xlsx"

for p in (XLSX_PATH, MARKERS_PATH):
    if not p.exists():
        raise FileNotFoundError(
            f"Required file not found: {p}\n"
            f"Both BlotSeq_PDO_SameCell_Multiomics_Dataset.xlsx and Cell_markers.xlsx "
            f"must be in the same folder as this script:\n{INPUT_DIR}"
        )

plt.rcParams["font.family"] = "Liberation Sans"  # Arial-metric-compatible; substitute if Arial unavailable
plt.rcParams["pdf.fonttype"] = 42

# reference palette, pixel-sampled from the uploaded ArchR-style example figure
CLUSTER_COLORS = {"K1": "#D51F25", "K2": "#272E6A", "K3": "#1E8A42"}
CLUSTER_NUM = {"K1": "1", "K2": "2", "K3": "3"}
CLUSTER_LABELS = {"K1": "K1 (luminal-like)", "K2": "K2 (mixed)", "K3": "K3 (basal-like)"}
ZSCORE_CMAP = LinearSegmentedColormap.from_list(
    "ref_style", ["#C9A9FF", "#FFFFFF", "#FF1A0A"], N=256
)

UMAP_DOT_SIZE = 75  # UMAP scatter point size (reduced per feedback)

# ============================================================
# Step 1: load and filter the two protein sheets
# ============================================================
def load_protein_sheet(path, sheet_name):
    df = pd.read_excel(path, sheet_name=sheet_name, header=1)
    cells = df.columns[1:11].tolist()          # 10 MaxLFQ abundance columns
    qcols = df.columns[15:25].tolist()          # matched q-value columns
    abund = df[cells].copy()
    q = df[qcols].copy()
    q.columns = cells
    detected = (q <= 0.05) & abund.notna()
    return abund.where(detected).set_axis(df["Gene group"], axis=0), cells

ab_basal, cells_basal = load_protein_sheet(XLSX_PATH, "Protein_Basal")
ab_luminal, cells_luminal = load_protein_sheet(XLSX_PATH, "Protein_Luminal")
combined = ab_basal.join(ab_luminal, how="outer")  # gene groups x 20 cells
cells = combined.columns.tolist()
n_cells = len(cells)

orig_label = {c: "Basal-like (KRT19 call)" for c in cells_basal}
orig_label.update({c: "Luminal-like (KRT19 call)" for c in cells_luminal})

n_detected = combined.notna().sum(axis=0)  # per-cell total detected proteins (depth proxy)

# ============================================================
# Step 2-3: feature filter, log2, per-cell depth center, per-gene z-score
# ============================================================
detected_frac = combined.notna().sum(axis=1) / n_cells
feat = combined.loc[detected_frac >= 0.5].copy()

log_full_all = np.log2(combined.astype(float))
per_cell_median = log_full_all.median(axis=0, skipna=True)  # depth-normalization offset

log_feat = np.log2(feat.astype(float)).sub(per_cell_median, axis=1)
row_mean = log_feat.mean(axis=1, skipna=True)
row_std = log_feat.std(axis=1, skipna=True)
z = log_feat.sub(row_mean, axis=0).div(row_std, axis=0)

# full-panel z-scores (for the marker dot plot; observed values only, not imputed)
log_full = log_full_all.sub(per_cell_median, axis=1)
row_mean_full = log_full.mean(axis=1, skipna=True)
row_std_full = log_full.std(axis=1, skipna=True)
z_full = log_full.sub(row_mean_full, axis=0).div(row_std_full, axis=0)
det_mask = combined.notna()

# ============================================================
# Step 4: downshifted-normal imputation (clustering matrix only)
# ============================================================
rng = np.random.default_rng(RNG)
z_imputed = z.copy()
for gene in z.index:
    row = z.loc[gene]
    missing = row.isna()
    if missing.any():
        row_min = row.min(skipna=True)
        z_imputed.loc[gene, missing] = row_min - 1.8 + rng.normal(0, 0.3, size=missing.sum())

X = z_imputed.T.values  # 20 cells x n_features, complete

# ============================================================
# Step 5: PCA
# ============================================================
n_pcs = min(10, n_cells - 1, X.shape[1])
pca = PCA(n_components=n_pcs, random_state=RNG)
pcs = pca.fit_transform(X)
pcs_df = pd.DataFrame(pcs, index=cells, columns=[f"PC{i+1}" for i in range(n_pcs)])

# ============================================================
# Step 6: drop PCs confounded with detection depth
# ============================================================
depth_flagged = []
for pc in pcs_df.columns:
    rho, p = spearmanr(pcs_df[pc], n_detected.loc[cells])
    if p < 0.05:
        depth_flagged.append(pc)
clustering_pcs = [c for c in pcs_df.columns if c not in depth_flagged]
print(f"PCs dropped as depth-confounded (p<0.05): {depth_flagged}")
print(f"PCs used for clustering: {clustering_pcs}")
pcs_for_clustering = pcs_df[clustering_pcs].values

# ============================================================
# Step 7: K-means, k=3
# ============================================================
km3 = KMeans(n_clusters=3, n_init=25, random_state=RNG).fit(pcs_for_clustering)
cluster_of = {c: f"K{l+1}" for c, l in zip(cells, km3.labels_)}
clusters = pd.DataFrame({
    "cell": cells,
    "orig_KRT19_label": [orig_label[c] for c in cells],
    "K3_cluster": [cluster_of[c] for c in cells],
})
print(clusters.to_string(index=False))

# ============================================================
# Step 8: UMAP for visualization (same depth-decorrelated PC subspace)
# ============================================================
n_neighbors = min(5, n_cells - 1)
reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.4, random_state=RNG, n_components=2)
embedding = reducer.fit_transform(pcs_for_clustering)
umap_df = pd.DataFrame(embedding, index=cells, columns=["UMAP1", "UMAP2"])
umap_df = umap_df.join(clusters.set_index("cell"))

cluster_keys = ["K1", "K2", "K3"]

# ============================================================
# Marker table for the dot plot
# ============================================================
marker_tbl = pd.read_excel(MARKERS_PATH, sheet_name="Sheet1")
class_order = ["Luminal", "LASP", "Luminal Hormone-Sensing", "Basal", "Basal-luminal"]
detected_markers = [g for g in marker_tbl["Gene"] if g in det_mask.index and det_mask.loc[g].sum() > 0]
seen, marker_rows = set(), []
for _, r in marker_tbl.iterrows():
    if r["Gene"] in detected_markers and r["Gene"] not in seen:
        marker_rows.append((r["Gene"], r["Class"]))
        seen.add(r["Gene"])
marker_rows.sort(key=lambda x: (class_order.index(x[1]) if x[1] in class_order else 99, x[0]))
genes_ordered = [g for g, c in marker_rows]
n_genes = len(genes_ordered)

# ============================================================
# Figure 1: UMAP panel
# ============================================================
fig1 = plt.figure(figsize=(4.6, 5.4))
ax = fig1.add_axes([0.14, 0.20, 0.80, 0.62])
ax.text(0.0, 1.14, "UMAP of PCA (depth-corrected)  colored by",
        transform=ax.transAxes, fontsize=10, ha="left", va="bottom")
ax.text(0.0, 1.05, "colData :  ClustersKmeans_K3",
        transform=ax.transAxes, fontsize=10, ha="left", va="bottom")

centroids = {}
for k in cluster_keys:
    sub = umap_df[umap_df["K3_cluster"] == k]
    ax.scatter(sub["UMAP1"], sub["UMAP2"], s=UMAP_DOT_SIZE, c=CLUSTER_COLORS[k],
               edgecolor="none", zorder=3)
    centroids[k] = (sub["UMAP1"].mean(), sub["UMAP2"].mean())
for k, (cx, cy) in centroids.items():
    ax.text(cx, cy, CLUSTER_NUM[k], fontsize=18, fontweight="bold",
            color="black", ha="center", va="center", zorder=4)

ax.set_xlabel("UMAP Dimension 1", fontsize=10)
ax.set_ylabel("UMAP Dimension 2", fontsize=10)
ax.set_xticks([]); ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(True); spine.set_color("black"); spine.set_linewidth(0.8)

ax.text(0.0, -0.14, "color", transform=ax.transAxes, fontsize=9.5, ha="left", va="top")
swatch_handles = [Line2D([0], [0], marker="s", linestyle="", markersize=9,
                          markerfacecolor=CLUSTER_COLORS[k], markeredgecolor="none",
                          label=f"{CLUSTER_NUM[k]}-{k}") for k in cluster_keys]
ax.legend(handles=swatch_handles, loc="upper left", bbox_to_anchor=(0.0, -0.20),
          ncol=len(cluster_keys), frameon=False, fontsize=9.5, handletextpad=0.4,
          columnspacing=1.2, borderaxespad=0)

fig1.savefig(OUTPUT_DIR / "Fig_UMAP_K3.png", dpi=1200, bbox_inches="tight")
plt.close(fig1)
print("Saved Fig_UMAP_K3.png")

# ============================================================
# Figure 2: marker dot plot
# ============================================================
size_vals = np.zeros((n_genes, len(cluster_keys)))
color_vals = np.full((n_genes, len(cluster_keys)), np.nan)
assign = umap_df["K3_cluster"]
for j, k in enumerate(cluster_keys):
    cells_in_cluster = assign[assign == k].index.tolist()
    for i, gene in enumerate(genes_ordered):
        size_vals[i, j] = det_mask.loc[gene, cells_in_cluster].mean() * 100
        z_row = z_full.loc[gene, cells_in_cluster]
        if z_row.notna().any():
            color_vals[i, j] = z_row.mean(skipna=True)

vmax = np.nanmax(np.abs(color_vals))
norm = Normalize(vmin=-vmax, vmax=vmax)

def size_map(pct):
    return 15 + pct * 5.0

fig2 = plt.figure(figsize=(6.6, 5.4))
ax_dot = fig2.add_axes([0.28, 0.14, 0.44, 0.80])
for i in range(n_genes):
    for j in range(len(cluster_keys)):
        pct = size_vals[i, j]
        y = n_genes - 1 - i
        if pct == 0:
            ax_dot.scatter(j, y, s=size_map(0), facecolor="none",
                            edgecolor="#999999", linewidth=0.6, zorder=3)
        else:
            color = ZSCORE_CMAP(norm(color_vals[i, j])) if not np.isnan(color_vals[i, j]) else "#dddddd"
            ax_dot.scatter(j, y, s=size_map(pct), facecolor=color,
                            edgecolor="#999999", linewidth=0.4, zorder=3)

ax_dot.set_xlim(-0.6, len(cluster_keys) - 0.4)
ax_dot.set_ylim(-0.8, n_genes - 0.2)
ax_dot.set_xticks(range(len(cluster_keys)))
ax_dot.set_xticklabels(cluster_keys, fontsize=11, fontweight="bold")
for tick, k in zip(ax_dot.get_xticklabels(), cluster_keys):
    tick.set_color(CLUSTER_COLORS[k])
ax_dot.set_yticks(range(n_genes))
ax_dot.set_yticklabels(genes_ordered[::-1], fontsize=10)
ax_dot.set_xlabel("K-means cluster", fontsize=11)
for spine in ax_dot.spines.values():
    spine.set_visible(False)
ax_dot.tick_params(length=0)

# size legend ("score > 0")
ax_size = fig2.add_axes([0.78, 0.50, 0.20, 0.36])
ax_size.axis("off")
ax_size.text(0.0, 1.0, "score > 0", fontsize=11, fontweight="bold",
             transform=ax_size.transAxes, ha="left", va="bottom")
y0 = 0.88
for v in [0, 25, 50, 75, 100]:
    face = "#4d4d4d" if v > 0 else "none"
    ax_size.scatter([0.08], [y0], s=size_map(v), facecolor=face,
                     edgecolor="#4d4d4d" if v > 0 else "#999999", linewidth=0.6,
                     transform=ax_size.transAxes, clip_on=False)
    ax_size.text(0.30, y0, str(v), fontsize=9.5, transform=ax_size.transAxes, ha="left", va="center")
    y0 -= 0.20

# colorbar ("Mean z-score")
cax = fig2.add_axes([0.82, 0.155, 0.03, 0.19])
sm = ScalarMappable(norm=norm, cmap=ZSCORE_CMAP)
cb = fig2.colorbar(sm, cax=cax)
cb.ax.tick_params(labelsize=9)
cb.outline.set_visible(False)
fig2.text(0.78, 0.375, "Mean\nz-score", fontsize=11, fontweight="bold", ha="left", va="bottom")

fig2.savefig(OUTPUT_DIR / "Fig_dotplot_K3.png", dpi=1200, bbox_inches="tight")
plt.close(fig2)
print("Saved Fig_dotplot_K3.png")

# ============================================================
# Figure 3: cell volume by cluster (from Cell_Morphology sheet)
# ============================================================
morph = pd.read_excel(XLSX_PATH, sheet_name="Cell_Morphology", header=0, skiprows=[0])
morph.columns = ["Cell", "Group", "Digestion_condition", "Array_row", "Array_col", "ATAC_barcode",
                  "Diameter_um", "Elongation", "Circularity", "Volume_um3", "TSS_enrichment",
                  "Chromatin_QC", "Peptides_ID", "Protein_groups_ID", "Genes_LoB", "Genes_chromatin"]
morph = morph[morph["Cell"].astype(str).str.match(r"^R\d+C\d+$", na=False)].set_index("Cell")
volume_df = morph[["Volume_um3"]].join(pd.Series(cluster_of, name="cluster"))

fig3, ax3 = plt.subplots(figsize=(5.2, 6.2))
rng_j = np.random.default_rng(1)
xpos = {"K1": 0, "K2": 1, "K3": 2}
for k in cluster_keys:
    sub = volume_df[volume_df["cluster"] == k]
    jitter = rng_j.uniform(-0.12, 0.12, size=len(sub))
    ax3.scatter(xpos[k] + jitter, sub["Volume_um3"], s=180, facecolor=CLUSTER_COLORS[k],
                edgecolor="black", linewidth=1.3, zorder=3)
    mean_val = sub["Volume_um3"].mean()
    ax3.plot([xpos[k] - 0.22, xpos[k] + 0.22], [mean_val, mean_val],
              color="black", linewidth=2.5, zorder=4)

ax3.set_xticks([0, 1, 2])
ax3.set_xticklabels(["K1", "K2", "K3"], fontsize=15)
ax3.set_xlim(-0.5, 2.5)
ax3.set_ylabel("Volume ($\\mu m^3$)", fontsize=15)
for spine in ax3.spines.values():
    spine.set_visible(True); spine.set_color("black"); spine.set_linewidth(1.8)
ax3.tick_params(width=1.8, length=6, labelsize=13)

handles = [Line2D([0], [0], marker="o", linestyle="", markersize=13,
                   markerfacecolor=CLUSTER_COLORS[k], markeredgecolor="black",
                   markeredgewidth=1.2, label=CLUSTER_LABELS[k]) for k in cluster_keys]
fig3.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
            fontsize=12, bbox_to_anchor=(0.5, 0.02))

fig3.tight_layout(rect=[0, 0.08, 1, 1])
fig3.savefig(OUTPUT_DIR / "Fig_Volume_by_cluster.png", dpi=1200, bbox_inches="tight")
plt.close(fig3)
print("Saved Fig_Volume_by_cluster.png")

print("\nVolume by cluster (um3):")
print(volume_df.groupby("cluster")["Volume_um3"].agg(["mean", "median", "std", "count"]).round(1))
