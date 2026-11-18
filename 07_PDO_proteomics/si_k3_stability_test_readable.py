import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
from scipy.stats import spearmanr
from itertools import combinations
from pathlib import Path

DATA_PATH = str(Path(__file__).resolve().parent / 'BlotSeq_PDO_SameCell_Multiomics_Dataset.xlsx')
OUT_PNG = str(Path(__file__).resolve().parent / 'SI_Fig_K3_stability_readable.png')

CLUSTER_COLORS = {'K1': '#D51F25', 'K2': '#272E6A', 'K3': '#1E8A42'}
GREY = '#9E9E9E'
SEEDS = list(range(100))
KS = [2, 3, 4, 5, 6]

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['pdf.fonttype'] = 42


def load_protein_sheet(path, sheet_name):
    df = pd.read_excel(path, sheet_name=sheet_name, header=1)
    cells = df.columns[1:11].tolist()
    qcols = df.columns[15:25].tolist()
    abund = df[cells].copy()
    q = df[qcols].copy()
    q.columns = cells
    detected = (q <= 0.05) & abund.notna()
    return abund.where(detected).set_axis(df['Gene group'], axis=0), cells


def style_axis(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(2.0)
        spine.set_color('black')
    ax.tick_params(width=2.0, length=6, labelsize=11)


def build_preprocessed_space():
    ab_basal, cells_basal = load_protein_sheet(DATA_PATH, 'Protein_Basal')
    ab_luminal, cells_luminal = load_protein_sheet(DATA_PATH, 'Protein_Luminal')

    combined = ab_basal.join(ab_luminal, how='outer')
    cells = combined.columns.tolist()
    n_cells = len(cells)

    orig_label = np.array(['Basal-like' if c in cells_basal else 'Luminal-like' for c in cells])
    n_detected = combined.notna().sum(axis=0).values

    detected_frac = combined.notna().sum(axis=1) / n_cells
    feat = combined.loc[detected_frac >= 0.5].copy()

    log_full_all = np.log2(combined.astype(float))
    per_cell_median = np.nanmedian(log_full_all.values, axis=0)

    log_feat = np.log2(feat.astype(float)).values - per_cell_median
    row_mean = np.nanmean(log_feat, axis=1, keepdims=True)
    row_std = np.nanstd(log_feat, axis=1, ddof=1, keepdims=True)
    z = (log_feat - row_mean) / row_std

    row_min = np.nanmin(z, axis=1, keepdims=True)
    missing = np.isnan(z)
    n_pcs = min(10, n_cells - 1, z.shape[0])

    return cells, orig_label, n_detected, z, row_min, missing, n_pcs


def preprocess_one_seed(seed, cells, n_detected, z, row_min, missing, n_pcs):
    rng = np.random.default_rng(seed)
    noise = rng.normal(0, 0.3, size=z.shape)
    z_imputed = np.where(missing, row_min - 1.8 + noise, z)
    pcs = PCA(n_components=n_pcs, random_state=0).fit_transform(z_imputed.T)
    keep_idx = [i for i in range(n_pcs) if spearmanr(pcs[:, i], n_detected).pvalue >= 0.05]
    return pcs[:, keep_idx]


def semantic_map_for_k3(reference_labels, orig_label):
    labs = np.unique(reference_labels)
    lum_frac = {lab: np.mean(orig_label[reference_labels == lab] == 'Luminal-like') for lab in labs}
    k1_numeric = max(lum_frac, key=lum_frac.get)
    remaining = [lab for lab in labs if lab != k1_numeric]
    basal_frac = {lab: np.mean(orig_label[reference_labels == lab] == 'Basal-like') for lab in remaining}
    k3_numeric = max(remaining, key=lambda lab: (basal_frac[lab], -np.sum(reference_labels == lab)))
    k2_numeric = [lab for lab in labs if lab not in {k1_numeric, k3_numeric}][0]
    return {k1_numeric: 'K1', k2_numeric: 'K2', k3_numeric: 'K3'}


def main():
    cells, orig_label, n_detected, z, row_min, missing, n_pcs = build_preprocessed_space()
    n_cells = len(cells)

    Y_by_seed = {seed: preprocess_one_seed(seed, cells, n_detected, z, row_min, missing, n_pcs) for seed in SEEDS}

    Y0 = Y_by_seed[0]
    km0 = KMeans(n_clusters=3, n_init=25, random_state=0).fit(Y0)
    reference_labels = km0.labels_
    label_map = semantic_map_for_k3(reference_labels, orig_label)
    semantic_labels = np.array([label_map[x] for x in reference_labels])

    ordered_cells = []
    for k in ['K1', 'K2', 'K3']:
        ordered_cells.extend([c for c, lab in zip(cells, semantic_labels) if lab == k])
    ordered_idx = [cells.index(c) for c in ordered_cells]
    ordered_semantic = [semantic_labels[i] for i in ordered_idx]

    silhouette_summary = []
    ari_summary = []
    labels_by_k = {}

    for k in KS:
        all_labels = []
        sil_values = []
        for seed in SEEDS:
            Y = Y_by_seed[seed]
            km = KMeans(n_clusters=k, n_init=25, random_state=seed).fit(Y)
            all_labels.append(km.labels_)
            sil_values.append(silhouette_score(Y, km.labels_))
        labels_by_k[k] = all_labels
        silhouette_summary.append({
            'k': k,
            'mean': float(np.mean(sil_values)),
            'low': float(np.percentile(sil_values, 2.5)),
            'high': float(np.percentile(sil_values, 97.5))
        })
        ari_values = np.array([adjusted_rand_score(all_labels[i], all_labels[j]) for i, j in combinations(range(len(all_labels)), 2)])
        ari_summary.append({
            'k': k,
            'mean': float(np.mean(ari_values)),
            'low': float(np.percentile(ari_values, 2.5)),
            'high': float(np.percentile(ari_values, 97.5))
        })

    consensus = np.zeros((n_cells, n_cells), dtype=float)
    for labels in labels_by_k[3]:
        consensus += (labels[:, None] == labels[None, :]).astype(float)
    consensus /= len(labels_by_k[3])
    consensus = consensus[np.ix_(ordered_idx, ordered_idx)]

    fig = plt.figure(figsize=(15.2, 5.2))

    # Panel A
    ax1 = fig.add_axes([0.055, 0.20, 0.245, 0.68])
    ks = np.array([d['k'] for d in silhouette_summary])
    means = np.array([d['mean'] for d in silhouette_summary])
    lows = np.array([d['low'] for d in silhouette_summary])
    highs = np.array([d['high'] for d in silhouette_summary])
    ax1.plot(ks, means, color='black', linewidth=1.6, zorder=2)
    ax1.errorbar(ks, means, yerr=[means - lows, highs - means], fmt='none', ecolor=GREY, elinewidth=1.3, capsize=3, zorder=1)
    for x, y in zip(ks, means):
        ax1.scatter(x, y, s=95, facecolor=(CLUSTER_COLORS['K2'] if x == 3 else 'white'), edgecolor=(CLUSTER_COLORS['K2'] if x == 3 else 'black'), linewidth=1.7, zorder=3)
    ax1.axvline(3, color=GREY, linestyle=':', linewidth=1.4)
    ax1.set_xticks(ks)
    ax1.set_xlabel('Number of clusters (k)', fontsize=12)
    ax1.set_ylabel('Mean silhouette score', fontsize=12)
    ax1.set_title('Clustering compactness', fontsize=13, fontweight='bold', pad=12)
    style_axis(ax1)
    ax1.set_ylim(min(lows) - 0.012, max(highs) + 0.020)
    ax1.text(-0.15, 1.05, 'A', transform=ax1.transAxes, fontsize=18, fontweight='bold')
    ax1.text(3.05, means[list(ks).index(3)] + 0.035, 'Chosen\nk = 3', color=CLUSTER_COLORS['K2'], fontsize=10, fontweight='bold', ha='center', va='bottom')

    # Panel B
    ax2 = fig.add_axes([0.355, 0.20, 0.245, 0.68])
    ks = np.array([d['k'] for d in ari_summary])
    means = np.array([d['mean'] for d in ari_summary])
    lows = np.array([d['low'] for d in ari_summary])
    highs = np.array([d['high'] for d in ari_summary])
    ax2.plot(ks, means, color='black', linewidth=1.6, zorder=2)
    ax2.errorbar(ks, means, yerr=[means - lows, highs - means], fmt='none', ecolor=GREY, elinewidth=1.3, capsize=3, zorder=1)
    for x, y in zip(ks, means):
        ax2.scatter(x, y, s=95, facecolor=(CLUSTER_COLORS['K2'] if x == 3 else 'white'), edgecolor=(CLUSTER_COLORS['K2'] if x == 3 else 'black'), linewidth=1.7, zorder=3)
    ax2.axvline(3, color=GREY, linestyle=':', linewidth=1.4)
    ax2.set_xticks(ks)
    ax2.set_xlabel('Number of clusters (k)', fontsize=12)
    ax2.set_ylabel('Mean pairwise ARI', fontsize=12)
    ax2.set_title('Clustering stability across reruns', fontsize=13, fontweight='bold', pad=12)
    style_axis(ax2)
    ax2.set_ylim(0.18, 1.045)
    ax2.text(-0.15, 1.05, 'B', transform=ax2.transAxes, fontsize=18, fontweight='bold')
    ax2.text(3.0, 0.24, 'Stable\n3-state\nsolution', color=CLUSTER_COLORS['K2'], fontsize=10, fontweight='bold', ha='center', va='bottom')

    # Panel C
    ax3 = fig.add_axes([0.665, 0.20, 0.245, 0.68])
    cmap = LinearSegmentedColormap.from_list('consensus_blues', ['#FFFFFF', '#DCE6F4', CLUSTER_COLORS['K2']])
    im = ax3.imshow(consensus, vmin=0, vmax=1, cmap=cmap, origin='upper', interpolation='nearest', aspect='auto')
    ax3.set_title('k = 3 consensus co-clustering', fontsize=12.8, fontweight='bold', pad=16)
    ax3.set_xticks(range(n_cells))
    ax3.set_yticks(range(n_cells))
    ax3.set_xticklabels(ordered_cells, rotation=90, fontsize=7.5)
    ax3.set_yticklabels(ordered_cells, fontsize=7.5)
    for tick, sem in zip(ax3.get_xticklabels(), ordered_semantic):
        tick.set_color(CLUSTER_COLORS[sem])
    for tick, sem in zip(ax3.get_yticklabels(), ordered_semantic):
        tick.set_color(CLUSTER_COLORS[sem])
    cluster_sizes = [sum(semantic_labels == 'K1'), sum(semantic_labels == 'K2'), sum(semantic_labels == 'K3')]
    separators = [cluster_sizes[0] - 0.5, cluster_sizes[0] + cluster_sizes[1] - 0.5]
    for sep in separators:
        ax3.axhline(sep, color='black', linewidth=1.5)
        ax3.axvline(sep, color='black', linewidth=1.5)
    start = 0
    for k, size in zip(['K1', 'K2', 'K3'], cluster_sizes):
        mid = start + (size - 1) / 2
        ax3.text(mid, -2.2, f'{k}\n(n={size})', color=CLUSTER_COLORS[k], ha='center', va='bottom', fontsize=10, fontweight='bold', clip_on=False)
        start += size
    style_axis(ax3)
    ax3.tick_params(length=0)
    ax3.text(-0.15, 1.05, 'C', transform=ax3.transAxes, fontsize=18, fontweight='bold')

    cax = fig.add_axes([0.925, 0.26, 0.015, 0.52])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Co-clustering\nfrequency', fontsize=11)
    cb.ax.tick_params(labelsize=10, width=1.2, length=4)

    fig.text(0.49, 0.065, 'k = 3 was selected as a compromise between compactness and resolution: k = 2 was slightly more compact but collapsed structure, whereas k ≥ 4 showed lower reproducibility.', ha='center', va='center', fontsize=10.5)

    fig.savefig(OUT_PNG, dpi=600, bbox_inches='tight')
    print('Saved:', OUT_PNG)

if __name__ == '__main__':
    main()
