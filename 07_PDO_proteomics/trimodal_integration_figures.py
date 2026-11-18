"""
trimodal_integration_figures.py

Adapts three previously-established trimodal figures (evidence composition,
chromatin-vs-protein, joint embedding) to the data-driven K1/K2/K3 proteome
clustering instead of the original binary KRT19-only Basal-like/Luminal-like
call. Single self-contained script, reads directly from the raw workbook.

CHROMATIN INTEGRATION RULE: any figure that uses ATAC GeneScores is restricted
to the 16 cells with TSS enrichment >= 0.9 (4 cells excluded: R16C24, R16C25,
R16C28, R4C31). The K1/K2/K3 cluster assignment itself is unaffected -- it is
derived from protein data alone, using all 20 cells, exactly as before.

INPUT (must sit in the same folder as this script):
    BlotSeq_PDO_SameCell_Multiomics_Dataset.xlsx
        sheets: Protein_Basal, Protein_Luminal, ATAC_Basal, ATAC_Luminal, Cell_Morphology
    Cell_markers.xlsx (not required for these three figures, kept for consistency
        with the rest of the pipeline if this script is extended later)

OUTPUT (written next to this script):
    Fig1_evidence_composition_K3.png
    Fig2_chromatin_vs_protein_K3.png
    Fig4_joint_embedding_K3.png

Figure logic (unchanged from the reference scripts except for the grouping
variable and the TSS filter):
  Fig1 -- per-cell stacked bar of gene-cell observations (protein only /
          chromatin only / both), over the "shared" gene set (the 474-gene
          ATAC panel, which per the workbook ReadMe was itself restricted to
          genes also seen at the protein level). Grouped by K1/K2/K3 instead
          of by original label.
  Fig2 -- chromatin GeneScore vs protein abundance for every gene-cell pair
          where BOTH were detected in the same physical cell. Reports raw
          Spearman rho and a rho computed after centering both variables
          within cell AND within gene (isolates the same-cell residual from
          between-gene and between-cell variation).
  Fig4 -- joint PCA: chromatin block (3 unit-variance PCs) + proteome block
          (3 unit-variance PCs) + morphometry block (log10 Volume, Diameter,
          z-scored) = explicit 3:3:2 weighting regardless of raw feature
          counts. Marker area = cell volume, color = K1/K2/K3.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy import stats
from scipy.stats import spearmanr

INPUT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = INPUT_DIR
RNG = 0

XLSX_PATH = INPUT_DIR / "BlotSeq_PDO_SameCell_Multiomics_Dataset.xlsx"
if not XLSX_PATH.exists():
    raise FileNotFoundError(f"Required file not found next to this script: {XLSX_PATH}")

# ---- palette: K1/K2/K3 (pixel-sampled earlier, established through this thread)
#      + the modality-stack colors pixel-sampled from the reference Fig1 legend
PAL = {
    "K1": "#D51F25", "K2": "#272E6A", "K3": "#1E8A42",
    "protein": "#CAA228", "chromatin": "#8FA9C4", "dual": "#3B3B3B",
    "grey": "#999999",
}
CLUSTER_LABELS = {"K1": "K1 (luminal-like)", "K2": "K2 (mixed)", "K3": "K3 (basal-like)"}

RC = {
    "font.family": "Liberation Sans",  # Arial-metric-compatible; substitute if unavailable
    "pdf.fonttype": 42,
    "axes.linewidth": 1.6,
    "xtick.major.width": 1.6, "ytick.major.width": 1.6,
    "xtick.major.size": 6, "ytick.major.size": 6,
    "font.size": 11,
}
mpl_rc = RC

def style_axes(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.6)
    ax.tick_params(width=1.6, length=6)

# ============================================================
# Step 1: rebuild K1/K2/K3 clustering from protein data (all 20 cells)
# ============================================================
def load_protein_sheet(path, sheet_name):
    df = pd.read_excel(path, sheet_name=sheet_name, header=1)
    cells = df.columns[1:11].tolist()
    qcols = df.columns[15:25].tolist()
    abund = df[cells].copy()
    q = df[qcols].copy()
    q.columns = cells
    detected = (q <= 0.05) & abund.notna()
    return abund.where(detected).set_axis(df["Gene group"], axis=0), cells

ab_basal, cells_basal = load_protein_sheet(XLSX_PATH, "Protein_Basal")
ab_luminal, cells_luminal = load_protein_sheet(XLSX_PATH, "Protein_Luminal")
combined = ab_basal.join(ab_luminal, how="outer")  # gene groups x 20 cells
cells_all = combined.columns.tolist()
n_cells_all = len(cells_all)
n_detected = combined.notna().sum(axis=0)

detected_frac = combined.notna().sum(axis=1) / n_cells_all
feat = combined.loc[detected_frac >= 0.5].copy()
log_full_all = np.log2(combined.astype(float))
per_cell_median = log_full_all.median(axis=0, skipna=True)
log_feat = np.log2(feat.astype(float)).sub(per_cell_median, axis=1)
row_mean = log_feat.mean(axis=1, skipna=True)
row_std = log_feat.std(axis=1, skipna=True)
z = log_feat.sub(row_mean, axis=0).div(row_std, axis=0)

rng = np.random.default_rng(RNG)
z_imputed = z.copy()
for gene in z.index:
    row = z.loc[gene]
    missing = row.isna()
    if missing.any():
        row_min = row.min(skipna=True)
        z_imputed.loc[gene, missing] = row_min - 1.8 + rng.normal(0, 0.3, size=missing.sum())
X = z_imputed.T.values

n_pcs = min(10, n_cells_all - 1, X.shape[1])
pca = PCA(n_components=n_pcs, random_state=RNG)
pcs = pca.fit_transform(X)
pcs_df = pd.DataFrame(pcs, index=cells_all, columns=[f"PC{i+1}" for i in range(n_pcs)])
depth_flagged = [pc for pc in pcs_df.columns if spearmanr(pcs_df[pc], n_detected.loc[cells_all])[1] < 0.05]
clustering_pcs = [c for c in pcs_df.columns if c not in depth_flagged]
km3 = KMeans(n_clusters=3, n_init=25, random_state=RNG).fit(pcs_df[clustering_pcs].values)
cluster_of = pd.Series({c: f"K{l+1}" for c, l in zip(cells_all, km3.labels_)}, name="cluster")
print("K1/K2/K3 (from protein data, all 20 cells):")
print(cluster_of.groupby(cluster_of).count())

# ============================================================
# Step 2: chromatin (ATAC) matrix, all 20 cells x 474-gene panel
# ============================================================
atac_b = pd.read_excel(XLSX_PATH, sheet_name="ATAC_Basal", header=1).set_index("Gene").iloc[:, :10]
atac_l = pd.read_excel(XLSX_PATH, sheet_name="ATAC_Luminal", header=1).set_index("Gene").iloc[:, :10]
atac = atac_b.join(atac_l, how="outer")  # 474 genes x 20 cells

# ============================================================
# Step 3: morphology / TSS -- determine the 16-cell TSS>=0.9 subset
# ============================================================
morph = pd.read_excel(XLSX_PATH, sheet_name="Cell_Morphology", header=0, skiprows=[0])
morph.columns = ["Cell", "Group", "Digestion_condition", "Array_row", "Array_col", "ATAC_barcode",
                  "Diameter_um", "Elongation", "Circularity", "Volume_um3", "TSS_enrichment",
                  "Chromatin_QC", "Peptides_ID", "Protein_groups_ID", "Genes_LoB", "Genes_chromatin"]
morph = morph[morph["Cell"].astype(str).str.match(r"^R\d+C\d+$", na=False)].set_index("Cell")

tss_pass = morph.index[morph["TSS_enrichment"] >= 0.9].tolist()
tss_excluded = [c for c in cells_all if c not in tss_pass]
print(f"\n{len(tss_pass)} cells pass TSS >= 0.9 (used for all chromatin-involving figures)")
print(f"{len(tss_excluded)} cells excluded: {tss_excluded}")

cluster16 = cluster_of.loc[tss_pass]
print("\nK1/K2/K3 composition after TSS filter:")
print(cluster16.groupby(cluster16).count())

# order cells within each cluster group for plotting (Fig1)
cluster_order = ["K1", "K2", "K3"]
cells_ordered = []
group_bounds = []
for k in cluster_order:
    members = sorted(cluster16[cluster16 == k].index.tolist())
    cells_ordered.extend(members)
    group_bounds.append(len(cells_ordered))

# ============================================================
# Step 4: resolve ATAC gene symbols against protein "Gene group" rows
# ============================================================
def resolve_gene_row(g, index):
    matches = [idx for idx in index if g in str(idx).split(";")]
    return matches[0] if matches else None

shared_genes = atac.index.tolist()  # 474-gene panel, already curated to be seen at protein level
gene_to_protrow = {g: resolve_gene_row(g, combined.index) for g in shared_genes}
resolved = {g: r for g, r in gene_to_protrow.items() if r is not None}
print(f"\n{len(resolved)}/{len(shared_genes)} ATAC panel genes resolved to a protein Gene-group row")

# detection matrix for the shared panel, restricted to resolvable genes
prot_shared = pd.DataFrame(index=list(resolved.keys()), columns=cells_all, dtype=float)
for g, r in resolved.items():
    prot_shared.loc[g] = combined.loc[r]
prot_det_shared = prot_shared.notna()          # protein detected (q<=0.05 & MaxLFQ present)
atac_acc_shared = (atac.loc[list(resolved.keys())] > 0)  # chromatin accessible

# ============================================================
# Figure 1: evidence composition, grouped by K1/K2/K3, 16-cell TSS subset
# ============================================================
det16 = prot_det_shared[cells_ordered].values
acc16 = atac_acc_shared[cells_ordered].values

n_prot_only = (~acc16 & det16).sum(axis=0)
n_atac_only = (acc16 & ~det16).sum(axis=0)
n_dual = (acc16 & det16).sum(axis=0)

print(f"\nFig1: {len(resolved)} shared genes x {len(cells_ordered)} cells")
print(f"  protein only   {n_prot_only.sum()}")
print(f"  chromatin only {n_atac_only.sum()}")
print(f"  both           {n_dual.sum()}")

plt.rcParams.update(mpl_rc)
fig1, ax1 = plt.subplots(figsize=(7.4, 3.6))
fig1.subplots_adjust(right=0.70, top=0.82)

x = np.arange(len(cells_ordered))
w = 0.78
ax1.bar(x, n_prot_only, w, color=PAL["protein"], ec="k", lw=0.5, zorder=3)
ax1.bar(x, n_atac_only, w, bottom=n_prot_only, color=PAL["chromatin"], ec="k", lw=0.5, zorder=3)
ax1.bar(x, n_dual, w, bottom=n_prot_only + n_atac_only, color=PAL["dual"], ec="k", lw=0.5, zorder=3)

for b in group_bounds[:-1]:
    ax1.axvline(b - 0.5, ls=":", c=PAL["grey"], lw=1.1, zorder=1)

ax1.set_xticks(x)
ax1.set_xticklabels(cells_ordered, rotation=90, fontsize=7)
for tick, cell in zip(ax1.get_xticklabels(), cells_ordered):
    tick.set_color(PAL[cluster16[cell]])
ax1.set_ylabel("Gene\u2013cell observations")
ax1.set_xlim(-0.7, len(cells_ordered) - 0.3)
ax1.set_ylim(0, (n_prot_only + n_atac_only + n_dual).max() * 1.08)
style_axes(ax1)

# cluster group headers, colored, above each block
prev = 0
for k, b in zip(cluster_order, group_bounds):
    cx = (prev + b - 1) / 2
    ax1.text(cx, 1.10, k, transform=ax1.get_xaxis_transform(), fontsize=13, fontweight="bold",
             color=PAL[k], ha="center", va="bottom")
    prev = b

handles = [plt.Rectangle((0, 0), 1, 1, fc=PAL["protein"], ec="k", lw=0.6),
           plt.Rectangle((0, 0), 1, 1, fc=PAL["chromatin"], ec="k", lw=0.6),
           plt.Rectangle((0, 0), 1, 1, fc=PAL["dual"], ec="k", lw=0.6),
           plt.Line2D([], [], marker="s", ls="", mfc="none", mec=PAL["K1"], mew=1.6, ms=8),
           plt.Line2D([], [], marker="s", ls="", mfc="none", mec=PAL["K2"], mew=1.6, ms=8),
           plt.Line2D([], [], marker="s", ls="", mfc="none", mec=PAL["K3"], mew=1.6, ms=8)]
labels = ["Proteome only", "Chromatin only", "Both modalities",
          "Cell label colour: K1", "Cell label colour: K2", "Cell label colour: K3"]
leg = ax1.legend(handles, labels, loc="center left", bbox_to_anchor=(1.03, 0.5),
                  frameon=False, fontsize=7.0, title="Key", title_fontsize=7.6,
                  handletextpad=0.9, labelspacing=1.2, borderaxespad=0)
leg.get_title().set_ha("left")

fig1.savefig(OUTPUT_DIR / "Fig1_evidence_composition_K3.png", dpi=1200, bbox_inches="tight")
plt.close(fig1)
print("Saved Fig1_evidence_composition_K3.png")

# ============================================================
# Figure 2: chromatin GeneScore vs protein abundance, dual-detected pairs
# ============================================================
rows = []
for c in tss_pass:
    for g in resolved:
        a = atac.loc[g, c]
        p = prot_shared.loc[g, c]
        if a > 0 and pd.notna(p):
            rows.append((c, cluster16[c], g, float(a), float(p)))
dual = pd.DataFrame(rows, columns=["cell", "cluster", "gene", "GeneScore", "abundance"])
dual["logGS"] = np.log10(dual.GeneScore)
dual["logAb"] = np.log10(dual.abundance)

rho_raw, p_raw = stats.spearmanr(dual.abundance, dual.GeneScore)

d = dual.copy()
d["gs_c"] = d.logGS - d.groupby("cell").logGS.transform("mean")
d["ab_c"] = d.logAb - d.groupby("cell").logAb.transform("mean")
d["gs_cg"] = d.gs_c - d.groupby("gene").gs_c.transform("mean")
d["ab_cg"] = d.ab_c - d.groupby("gene").ab_c.transform("mean")
rho_cg, p_cg = stats.spearmanr(d.gs_cg, d.ab_cg)

print(f"\nFig2: {len(dual)} dual observations, {dual.gene.nunique()} genes, {dual.cell.nunique()} cells")
print(f"  raw                      rho = {rho_raw:+.3f}, P = {p_raw:.4f}")
print(f"  centred within cell+gene rho = {rho_cg:+.3f}, P = {p_cg:.3f}")

fig2, ax2 = plt.subplots(figsize=(6.4, 4.0))
fig2.subplots_adjust(right=0.62)
col = dual["cluster"].map(PAL).values
ax2.scatter(dual.abundance, dual.GeneScore, s=18, c=col, ec="k", lw=0.3, alpha=0.85, zorder=3)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("Protein abundance in the same cell (MaxLFQ)")
ax2.set_ylabel("Chromatin gene activity score\nin the same nucleus")
style_axes(ax2)

handles = [Line2D([], [], marker="o", ls="", mfc=PAL[k], mec="k", mew=0.5, ms=7) for k in cluster_order]
labels = [CLUSTER_LABELS[k] for k in cluster_order]
leg = ax2.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.04, 1.0),
                  frameon=False, fontsize=7.2, title="Key", title_fontsize=7.6,
                  handletextpad=0.9, labelspacing=1.3, borderaxespad=0)
leg.get_title().set_ha("left")

stats_txt = (f"n = {len(dual)} gene\u2013cell pairs\n"
             f"{dual.gene.nunique()} genes, {dual.cell.nunique()} cells\n"
             f"raw \u03c1 = {rho_raw:+.2f}, P = {p_raw:.3f}\n"
             f"within cell+gene \u03c1 = {rho_cg:+.2f}, P = {p_cg:.2f}")
ax2.annotate(stats_txt, xy=(1.04, 0.0), xycoords="axes fraction",
             fontsize=6.6, ha="left", va="bottom", linespacing=1.5)

fig2.savefig(OUTPUT_DIR / "Fig2_chromatin_vs_protein_K3.png", dpi=1200, bbox_inches="tight")
plt.close(fig2)
print("Saved Fig2_chromatin_vs_protein_K3.png")

# ============================================================
# Figure 4: joint embedding, 16-cell TSS subset, K1/K2/K3 colored
# ============================================================
N_PC_MOLECULAR = 3
MIN_CELLS_PROT = 12
MIN_CELLS_ATAC = 5

def unit_pcs(Xb, k):
    Xz = (Xb - Xb.mean(axis=0)) / (Xb.std(axis=0) + 1e-9)
    S = PCA(n_components=k, random_state=RNG).fit_transform(Xz)
    return S / S.std(axis=0)

vol16 = morph.loc[tss_pass, "Volume_um3"].values
dia16 = morph.loc[tss_pass, "Diameter_um"].values

Pl = np.log10(combined[tss_pass].astype(float))
prot_feat = Pl.loc[Pl.notna().sum(axis=1) >= MIN_CELLS_PROT].T  # cells x genes
prot_feat = prot_feat.fillna(np.nanmin(prot_feat.values))

atac_feat = atac.loc[(atac[tss_pass] > 0).sum(axis=1) >= MIN_CELLS_ATAC, tss_pass].T  # cells x genes

morph_feat = np.column_stack([stats.zscore(np.log10(vol16)), stats.zscore(dia16)])

B_atac = unit_pcs(atac_feat.values, N_PC_MOLECULAR)
B_prot = unit_pcs(prot_feat.values, N_PC_MOLECULAR)
J = np.column_stack([B_atac, B_prot, morph_feat])
blocks = np.array(["chromatin"] * N_PC_MOLECULAR + ["proteome"] * N_PC_MOLECULAR +
                   ["morphometry"] * morph_feat.shape[1])

pj = PCA(n_components=4, random_state=RNG).fit(J)
S = pj.transform(J)
ev = 100 * pj.explained_variance_ratio_

print(f"\nFig4: {len(tss_pass)} cells; proteome features {prot_feat.shape[1]}, "
      f"chromatin features {atac_feat.shape[1]}")
print(f"joint PC variance explained: {np.round(ev, 1)}")
cluster_arr = cluster16.loc[tss_pass].values
for k in range(4):
    w = pj.components_[k] ** 2
    w = w / w.sum()
    row = {b: 100 * w[blocks == b].sum() for b in ["chromatin", "proteome", "morphometry"]}
    kw_stat, kw_p = stats.kruskal(*[S[cluster_arr == c, k] for c in cluster_order])
    print(f"  PC{k+1} ({ev[k]:.1f}% var): " + "  ".join(f"{b} {v:.0f}%" for b, v in row.items()) +
          f"   |  Kruskal-Wallis across K1/K2/K3, P = {kw_p:.4f}")

fig4, ax4 = plt.subplots(figsize=(5.8, 4.0))
fig4.subplots_adjust(right=0.64)
size = 30 + 260 * (vol16 - vol16.min()) / (vol16.max() - vol16.min())
col = pd.Series(cluster_arr).map(PAL).values

ax4.axhline(0, ls=":", c=PAL["grey"], lw=1.1, zorder=1)
ax4.axvline(0, ls=":", c=PAL["grey"], lw=1.1, zorder=1)
ax4.scatter(S[:, 0], S[:, 1], s=size, c=col, ec="k", lw=0.6, alpha=0.92, zorder=3)

ax4.set_xlabel(f"Joint PC1 ({ev[0]:.1f}%)")
ax4.set_ylabel(f"Joint PC2 ({ev[1]:.1f}%)")
lim = np.abs(S[:, :2]).max() * 1.18
ax4.set_xlim(-lim, lim)
ax4.set_ylim(-lim, lim)
style_axes(ax4)

handles = [Line2D([], [], marker="o", ls="", mfc=PAL[k], mec="k", mew=0.6, ms=8) for k in cluster_order] + [
    Line2D([], [], marker="o", ls="", mfc="0.6", mec="k", mew=0.6, ms=5),
    Line2D([], [], marker="o", ls="", mfc="0.6", mec="k", mew=0.6, ms=12)]
labels = [CLUSTER_LABELS[k] for k in cluster_order] + [
    f"Volume {vol16.min():,.0f} \u00b5m$^3$", f"Volume {vol16.max():,.0f} \u00b5m$^3$"]
leg = ax4.legend(handles, labels, loc="center left", bbox_to_anchor=(1.06, 0.5),
                  frameon=False, fontsize=7.2, title="Key", title_fontsize=7.6,
                  handletextpad=0.9, labelspacing=1.3, borderaxespad=0)
leg.get_title().set_ha("left")

pd.DataFrame(S, index=tss_pass, columns=[f"JointPC{k+1}" for k in range(4)]).round(4).to_csv(
    OUTPUT_DIR / "Fig4_joint_embedding_K3_scores.csv")
fig4.savefig(OUTPUT_DIR / "Fig4_joint_embedding_K3.png", dpi=1200, bbox_inches="tight")
plt.close(fig4)
print("Saved Fig4_joint_embedding_K3.png")
