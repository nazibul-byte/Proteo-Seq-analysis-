"""
fig2_chromatin_protein.py

Chromatin gene activity score against protein abundance, one point per gene-cell
pair where both modalities returned a value from the same physical cell.

This is the measurement no existing single-cell platform produces: antibody-based
multiomic assays cannot supply an unbiased protein axis, and mass-spectrometry
single-cell proteomics consumes the nucleus before it can be sequenced.

Two correlation estimates are reported. The raw Spearman pools between-gene and
between-cell variation. Centring both variables within cell and within gene
isolates the same-cell residual, which is the quantity of interest.

Output: 1200 dpi PNG + SVG.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt

import blotseq_io as bio

# --------------------------------------------------------------------- paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR  = SCRIPT_DIR
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'figures')
STEM = 'Fig2_chromatin_vs_protein'

# --------------------------------------------------------------------- data
D = bio.build_matrices(INPUT_DIR)
cells, shared = D['cells'], D['shared']

rows = []
for c in cells:
    for g in shared:
        a = D['atac'].loc[c, g]
        p = D['ab'].loc[c, g]
        if a > 0 and pd.notna(p):
            rows.append((c, bio.LINEAGE[c], g, float(a), float(p)))
dual = pd.DataFrame(rows, columns=['cell', 'lineage', 'gene', 'GeneScore', 'abundance'])
dual['logGS'] = np.log10(dual.GeneScore)
dual['logAb'] = np.log10(dual.abundance)

rho_raw, p_raw = stats.spearmanr(dual.GeneScore, dual.logAb)

# double centring: remove each cell's mean and each gene's mean from both axes
d = dual.copy()
d['gs_c'] = d.logGS - d.groupby('cell').logGS.transform('mean')
d['ab_c'] = d.logAb - d.groupby('cell').logAb.transform('mean')
d['gs_cg'] = d.gs_c - d.groupby('gene').gs_c.transform('mean')
d['ab_cg'] = d.ab_c - d.groupby('gene').ab_c.transform('mean')
rho_cg, p_cg = stats.spearmanr(d.gs_cg, d.ab_cg)

print(f'{len(dual)} dual observations, {dual.gene.nunique()} genes, {dual.cell.nunique()} cells')
print(f'  raw                     rho = {rho_raw:+.3f}, P = {p_raw:.4f}')
print(f'  centred within cell+gene rho = {rho_cg:+.3f}, P = {p_cg:.3f}')
print(f'  protein span {dual.logAb.max()-dual.logAb.min():.2f} decades')

# --------------------------------------------------------------------- plot
mpl.rcParams.update(bio.RC)
fig, ax = plt.subplots(figsize=(6.2, 3.9))
fig.subplots_adjust(right=0.60)

col = np.where(dual.lineage.values == 'Luminal', bio.PAL['luminal'], bio.PAL['basal'])
ax.scatter(dual.abundance, dual.GeneScore, s=17, c=col,
           ec='k', lw=0.3, alpha=0.85, zorder=3)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Protein abundance in the same cell (MaxLFQ)')
ax.set_ylabel('Chromatin gene activity score\nin the same nucleus')
bio.style_axes(ax)

# ---- external key, right of the axes ---------------------------------------
handles = [plt.Line2D([], [], marker='o', ls='', mfc=bio.PAL['basal'], mec='k', mew=0.5, ms=7),
           plt.Line2D([], [], marker='o', ls='', mfc=bio.PAL['luminal'], mec='k', mew=0.5, ms=7)]
labels = ['Basal-like', 'Luminal-like']
leg = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.04, 1.0),
                frameon=False, fontsize=7.2, title='Key', title_fontsize=7.6,
                handletextpad=0.9, labelspacing=1.3, borderaxespad=0)
leg.get_title().set_ha('left')

stats_txt = (f'n = {len(dual)} gene\u2013cell pairs\n'
             f'{dual.gene.nunique()} genes, {dual.cell.nunique()} cells\n'
             f'raw \u03c1 = {rho_raw:+.2f}, P = {p_raw:.3f}\n'
             f'within cell+gene \u03c1 = {rho_cg:+.2f}, P = {p_cg:.2f}')
ax.annotate(stats_txt, xy=(1.04, 0.0), xycoords='axes fraction',
            fontsize=6.6, ha='left', va='bottom', linespacing=1.5)

os.makedirs(OUTPUT_DIR, exist_ok=True)
for ext in ('png', 'svg'):
    fig.savefig(os.path.join(OUTPUT_DIR, f'{STEM}.{ext}'), dpi=1200)
print('saved', STEM)
