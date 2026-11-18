#!/usr/bin/env python3
"""
BlotSeq MCF7 MS Analysis — Unified Figure & Data Pipeline
==========================================================
Produces:
  Panel 1: Bar chart — unique protein groups per sample (AMHE023 + AMHE025 blank)
  Panel 2: SPIDR validation — stacked bar of match categories
  Panel 3: Protein composition — GO:Cellular Component stacked bars by timepoint

Also exports:
  BlotSeq_MCF7_MS_analysis.xlsx — multi-sheet Excel with underlying data for each panel

Inputs (adjust INPUT_DIR below):
  - AMHE023-report.tsv          (BlotSeq DIA-NN output, 15 single cells)
  - AMHE023-report_stats.tsv    (DIA-NN stats for AMHE023)
  - AMHE025-report_blank.tsv    (BlotSeq blank + ParcelOpt DIA-NN output)
  - AMHE025-report_blank_stats.tsv (DIA-NN stats for AMHE025)
  - SPIDR_MCF7_Proteins.xlsx    (Wang et al. 2025 MCF-7 scMS, for Panel 2)

Requirements:
  pip install pandas matplotlib numpy openpyxl gprofiler-official mygene scipy

CONFIGURATION: Edit INPUT_DIR and OUTPUT_DIR below before running locally.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from collections import defaultdict
import os
import warnings
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════════════
# PUBLICATION STYLE — white background, minimal spines, sans-serif font
# ═══════════════════════════════════════════════════════════════════════════
import matplotlib
matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 7,
    'axes.labelsize': 8,
    'axes.titlesize': 9,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'savefig.dpi': 1200,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# Core palette requested by user
C_TEXT    = '#2D2D2D'
C_ORANGE  = '#D25A3C'
C_ORANGE2 = '#C8643C'
C_INDIGO  = '#3C3C96'
C_MAROON  = '#781E1E'
C_TEAL    = '#2D695A'
C_TEAL2   = '#1E7887'
C_PURPLE  = '#3C0F4B'
C_GREEN   = '#1EA578'
C_GREEN2  = '#78C34B'
C_YELLOW  = '#F0E11E'
C_GREY    = '#A0A0A0'
C_LGREY   = '#D9D9D9'
C_BLANK   = '#ECECEC'

def style_axes(ax, grid_y=True):
    ax.set_facecolor('white')
    ax.tick_params(colors=C_TEXT, width=0.5, length=3)
    ax.xaxis.label.set_color(C_TEXT)
    ax.yaxis.label.set_color(C_TEXT)
    ax.title.set_color(C_TEXT)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(C_TEXT)
    ax.spines['bottom'].set_color(C_TEXT)
    ax.spines['left'].set_linewidth(0.6)
    ax.spines['bottom'].set_linewidth(0.6)
    if grid_y:
        ax.yaxis.grid(True, color='#ECECEC', linewidth=0.4, zorder=0)
        ax.set_axisbelow(True)

def style_legend(legend):
    if legend is None:
        return
    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor(C_LGREY)
    frame.set_linewidth(0.4)
    frame.set_alpha(0.95)

def save_multi(fig, basename):
    for fmt in ['png', 'svg']:
        fig.savefig(os.path.join(OUTPUT_DIR, f'{basename}.{fmt}'), dpi=1200, bbox_inches='tight', facecolor='white')
        print(f'  Saved: {basename}.{fmt}')


# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════
INPUT_DIR  = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

AMHE023_REPORT = os.path.join(INPUT_DIR, 'AMHE023-report.tsv')
AMHE023_STATS  = os.path.join(INPUT_DIR, 'AMHE023-report_stats.tsv')
AMHE025_REPORT = os.path.join(INPUT_DIR, 'AMHE025-report_blank.tsv')
AMHE025_STATS  = os.path.join(INPUT_DIR, 'AMHE025-report_blak_stats.tsv')
SPIDR_FILE     = os.path.join(INPUT_DIR, 'SPIDR_MCF7_Proteins.xlsx')

EXCEL_OUT = os.path.join(OUTPUT_DIR, 'BlotSeq_MCF7_MS_analysis.xlsx')

# ═══════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

def short_name_amhe023(run):
    """Extract readable sample name from AMHE023 Run string.
    e.g. 'Blotseq2hr_1cell_F1_R2_T1_..._S2-G1_...' -> 'Blotseq2hr_R2_S2-G1'
    """
    parts = run.split('_')
    cond = parts[0]     # Blotseq2hr
    rep = parts[3]      # R2
    well = parts[9]     # S2-G1
    return f'{cond}_{rep}_{well}'

def get_condition(run):
    if 'Blotseq1hr' in run or 'Blotseq1hr' in run:
        return 'Blotseq 1hr'
    elif 'Blotseq2hr' in run:
        return 'Blotseq 2hr'
    elif 'Blotseq3hr' in run:
        return 'Blotseq 3hr'
    return 'Other'

def get_timepoint_color(name):
    """Return publication-style color based on timepoint."""
    if 'Blotseq1hr' in name or '1hr' in name:
        return C_TEAL2
    elif 'Blotseq2hr' in name or '2hr' in name:
        return C_ORANGE
    elif 'Blotseq3hr' in name or '3hr' in name:
        return C_INDIGO
    elif 'Blank' in name or 'blank' in name:
        return C_BLANK
    return C_GREY

# Slightly different shades per replicate within each timepoint
TIMEPOINT_PALETTES = {
    '1hr': [C_TEAL2, C_TEAL, '#5A9184', '#76A79D', '#93BFB6'],
    '2hr': [C_ORANGE, C_ORANGE2, '#D87A55', '#E09A7D', '#E9BAA8'],
    '3hr': [C_INDIGO, '#5B5BB0', '#7777BE', '#9494CF', '#B1B1DF'],
}


# ═══════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
print('Loading AMHE023 report...')
df_023 = pd.read_csv(AMHE023_REPORT, sep='\t',
                     usecols=['Run', 'Protein.Group', 'Genes', 'Protein.Names',
                              'PG.MaxLFQ', 'PG.Normalised'])

print('Loading AMHE025 report (blank)...')
df_025 = pd.read_csv(AMHE025_REPORT, sep='\t',
                     usecols=['Run', 'Protein.Group', 'Genes', 'Global.PG.Q.Value'])

# Filter to BlotSeq blank only
df_blank = df_025[df_025['Run'].str.contains('BlotseqBlank')]

print('Loading AMHE023 stats...')
try:
    stats_023 = pd.read_csv(AMHE023_STATS, sep='\t')
except FileNotFoundError:
    print(f'  WARNING: {AMHE023_STATS} not found — skipping (not required for figures)')
    stats_023 = None

print('Loading AMHE025 stats...')
try:
    stats_025 = pd.read_csv(AMHE025_STATS, sep='\t')
except FileNotFoundError:
    print(f'  WARNING: {AMHE025_STATS} not found — skipping (not required for figures)')
    stats_025 = None

# ═══════════════════════════════════════════════════════════════════════════
# PANEL 1: BAR CHART — Unique Protein Groups per Sample
# ═══════════════════════════════════════════════════════════════════════════
print('\n--- Panel 1: Bar chart ---')

# Count unique protein groups per AMHE023 sample
protein_counts = (
    df_023.groupby('Run')['Protein.Group']
    .nunique()
    .reset_index()
    .rename(columns={'Protein.Group': 'Proteins_Identified'})
)
protein_counts['Sample'] = protein_counts['Run'].apply(short_name_amhe023)
protein_counts['Condition'] = protein_counts['Run'].apply(get_condition)

# Check HSP90AA1 detection per sample
hsp_runs = set(
    df_023[df_023['Genes'].str.contains('HSP90AA1', na=False)]['Run'].unique()
)
protein_counts['HSP90AA1_detected'] = protein_counts['Run'].isin(hsp_runs)

# Sort descending
protein_counts = protein_counts.sort_values('Proteins_Identified', ascending=False).reset_index(drop=True)

# Add the AMHE025 blank as the 16th bar
# From our analysis: Global.PG.Q.Value = 1.0 for all blank entries → 0 proteins at FDR
blank_precursors = df_blank['Protein.Group'].nunique()  # 229 nominal, but 0 pass FDR
blank_row = pd.DataFrame([{
    'Run': 'BlotseqBlank_0cell',
    'Proteins_Identified': 0,
    'Sample': 'BlotseqBlank_0cell',
    'Condition': 'Blank',
    'HSP90AA1_detected': False,
}])
protein_counts = pd.concat([protein_counts, blank_row], ignore_index=True)

print(f'  Samples: {len(protein_counts)} (15 cells + 1 blank)')
print(f'  HSP90AA1 detected in {protein_counts["HSP90AA1_detected"].sum()} / 15 cells')

# Assign colors: shade gradient within each timepoint
color_map = []
tp_counters = {'1hr': 0, '2hr': 0, '3hr': 0}
for _, row in protein_counts.iterrows():
    sample = row['Sample']
    if 'Blank' in sample or 'blank' in sample:
        color_map.append('#95A5A6')
    elif '1hr' in sample:
        idx = tp_counters['1hr']
        color_map.append(TIMEPOINT_PALETTES['1hr'][min(idx, 4)])
        tp_counters['1hr'] += 1
    elif '2hr' in sample:
        idx = tp_counters['2hr']
        color_map.append(TIMEPOINT_PALETTES['2hr'][min(idx, 4)])
        tp_counters['2hr'] += 1
    elif '3hr' in sample:
        idx = tp_counters['3hr']
        color_map.append(TIMEPOINT_PALETTES['3hr'][min(idx, 4)])
        tp_counters['3hr'] += 1
    else:
        color_map.append('#BDC3C7')

# Plot
fig1, ax1 = plt.subplots(figsize=(10.5, 3.4))
bars = ax1.bar(
    range(len(protein_counts)),
    protein_counts['Proteins_Identified'],
    color=color_map,
    edgecolor='none',
    linewidth=0,
    zorder=3,
)

# Annotate count + HSP90 asterisk
for bar, (_, row) in zip(bars, protein_counts.iterrows()):
    count = int(row['Proteins_Identified'])
    y_offset = max(count, 2) + 1.8
    ax1.text(bar.get_x() + bar.get_width() / 2, y_offset, str(count),
             ha='center', va='bottom', fontsize=5.8, color=C_TEXT)
    if row['HSP90AA1_detected']:
        ax1.plot(bar.get_x() + bar.get_width() / 2, y_offset + 5.0, marker='*',
                 color=C_MAROON, markersize=5.5, markeredgecolor='none', zorder=5)

# Legends, compact and inside plot area
tp_patches = [
    mpatches.Patch(color=TIMEPOINT_PALETTES['1hr'][0], label='1 hr'),
    mpatches.Patch(color=TIMEPOINT_PALETTES['2hr'][0], label='2 hr'),
    mpatches.Patch(color=TIMEPOINT_PALETTES['3hr'][0], label='3 hr'),
    mpatches.Patch(color=C_BLANK, label='Blank'),
    Line2D([0], [0], marker='*', color='w', markerfacecolor=C_MAROON,
           markeredgecolor='none', markersize=6, label='HSP90AA1'),
]
leg1 = ax1.legend(handles=tp_patches, loc='upper right', ncol=1, frameon=True,
                  handlelength=1.0, borderpad=0.35, labelspacing=0.3,
                  fancybox=False)
style_legend(leg1)

ax1.set_xlabel('Sample')
ax1.set_ylabel('Protein groups')
ax1.set_title('Unique protein groups identified per sample', pad=8)
ax1.set_xticks(range(len(protein_counts)))
ax1.set_xticklabels(protein_counts['Sample'], rotation=45, ha='right')
ax1.set_ylim(0, protein_counts['Proteins_Identified'].max() + 22)
style_axes(ax1, grid_y=True)

plt.tight_layout()
save_multi(fig1, 'Panel1_BarChart_ProteinsPerSample')
plt.close(fig1)


# ═══════════════════════════════════════════════════════════════════════════
# PANEL 2: SPIDR VALIDATION — Stacked bar of match categories
# ═══════════════════════════════════════════════════════════════════════════
print('\n--- Panel 2: SPIDR validation ---')

# Build protein group summary from AMHE023
pg_data = df_023.groupby('Protein.Group').agg({
    'Genes': 'first',
    'Protein.Names': 'first',
    'PG.MaxLFQ': 'mean',
    'Run': 'nunique',
}).reset_index()
pg_data['Is_Contaminant'] = pg_data['Protein.Group'].str.startswith('Cont_')
pg_data.rename(columns={'Run': 'N_Samples_Detected'}, inplace=True)

# Load SPIDR (if available)
spidr_available = os.path.exists(SPIDR_FILE)
if spidr_available:
    print('  Loading SPIDR MCF-7 proteome...')
    df_spidr = pd.read_excel(SPIDR_FILE, sheet_name='MCF7_Proteins')
    df_spidr['Log10_Quantity'] = np.log10(df_spidr['Mean_Quantity'].clip(lower=1e-10))

    # Build gene lookup
    spidr_gene_to_idx = {}
    for idx, row in df_spidr.iterrows():
        gn = row['PG.Genes']
        if pd.isna(gn):
            continue
        for g in str(gn).split(';'):
            g = g.strip().upper()
            if g:
                spidr_gene_to_idx.setdefault(g, idx)

    # Classify each BlotSeq protein
    validated_proteins = []
    match_status = []
    for _, row in pg_data.iterrows():
        gene = row['Genes']
        gene_upper = gene.strip().upper() if pd.notna(gene) else ''
        if gene_upper and gene_upper in spidr_gene_to_idx:
            sidx = spidr_gene_to_idx[gene_upper]
            sr = df_spidr.loc[sidx]
            match_status.append('Detected')
            validated_proteins.append({
                'Gene': gene,
                'SPIDR_Log10Qty': sr['Log10_Quantity'],
                'SPIDR_Rank': sr['Rank'],
                'SPIDR_N_Cells': sr['N_Cells_Detected'],
                'Is_Cont': row['Is_Contaminant'],
            })
        else:
            match_status.append('Not detected')
    pg_data['Match_Status'] = match_status

    n_det_noncont    = ((pg_data['Match_Status'] == 'Detected') & (~pg_data['Is_Contaminant'])).sum()
    n_det_cont       = ((pg_data['Match_Status'] == 'Detected') & (pg_data['Is_Contaminant'])).sum()
    n_notdet_noncont = ((pg_data['Match_Status'] == 'Not detected') & (~pg_data['Is_Contaminant'])).sum()
    n_notdet_cont    = ((pg_data['Match_Status'] == 'Not detected') & (pg_data['Is_Contaminant'])).sum()
    n_total = len(pg_data)
    n_detected = n_det_noncont + n_det_cont
    n_not_detected = n_notdet_noncont + n_notdet_cont

    print(f'  Total protein groups: {n_total}')
    print(f'  Detected in SPIDR: {n_detected} (non-cont: {n_det_noncont}, cont: {n_det_cont})')
    print(f'  Not detected:      {n_not_detected} (non-cont: {n_notdet_noncont}, cont: {n_notdet_cont})')

    # Colors
    C_DET_NONCONT    = C_INDIGO
    C_DET_CONT       = '#9E9ED7'
    C_NOTDET_NONCONT = C_ORANGE
    C_NOTDET_CONT    = '#E9BAA8'

    fig2, axes2 = plt.subplots(1, 2, figsize=(9.6, 3.6),
                               gridspec_kw={'width_ratios': [1, 1.35]})
    fig2.subplots_adjust(wspace=0.38, left=0.08, right=0.98, top=0.88, bottom=0.18)

    # Panel A: Stacked horizontal bar
    ax = axes2[0]
    bar_segments = [
        ('Detected, non-contaminant',        n_det_noncont,    C_DET_NONCONT),
        ('Detected, DIA-NN contaminant',     n_det_cont,       C_DET_CONT),
        ('Not detected, non-contaminant',    n_notdet_noncont, C_NOTDET_NONCONT),
        ('Not detected, DIA-NN contaminant', n_notdet_cont,    C_NOTDET_CONT),
    ]
    bar_height = 0.5
    left = 0
    for label, count, color in bar_segments:
        if count == 0:
            continue
        ax.barh(0, count, left=left, height=bar_height,
                color=color, edgecolor='white', linewidth=0.8)
        if count >= 4:
            text_color = 'white' if color == C_DET_NONCONT else '#333333'
            ax.text(left + count / 2, 0, str(count),
                    ha='center', va='center', fontsize=12, fontweight='bold',
                    color=text_color)
        left += count

    # Bracket annotations
    bracket_y = bar_height / 2 + 0.12
    ax.annotate('', xy=(0, bracket_y + 0.08), xytext=(n_detected, bracket_y + 0.08),
                arrowprops=dict(arrowstyle='-', color='#2E86AB', lw=2))
    ax.plot([0, 0], [bracket_y + 0.04, bracket_y + 0.12], color='#2E86AB', lw=2)
    ax.plot([n_detected, n_detected], [bracket_y + 0.04, bracket_y + 0.12], color='#2E86AB', lw=2)
    ax.text(n_detected / 2, bracket_y + 0.23,
            f'Matched to published MCF-7 single-cell proteome\n({n_detected}/{n_total})',
            ha='center', va='bottom', fontsize=5, color='#2E86AB')

    ax.annotate('', xy=(n_detected, bracket_y + 0.08), xytext=(n_total, bracket_y + 0.08),
                arrowprops=dict(arrowstyle='-', color='#C0392B', lw=2))
    ax.plot([n_detected, n_detected], [bracket_y + 0.04, bracket_y + 0.12], color='#C0392B', lw=2)
    ax.plot([n_total, n_total], [bracket_y + 0.04, bracket_y + 0.12], color='#C0392B', lw=2)
    ax.text(n_detected + n_not_detected / 2, bracket_y + 0.23,
            f'Not matched to published MCF-7 single-cell proteome\n({n_not_detected}/{n_total})',
            ha='center', va='bottom', fontsize=5, color='#C0392B')

    legend_elements = [
        mpatches.Patch(facecolor=C_DET_NONCONT, edgecolor='white', linewidth=0.5,
                       label=f'Matched to published proteome, non-contaminant ({n_det_noncont})'),
        mpatches.Patch(facecolor=C_DET_CONT, edgecolor='white', linewidth=0.5,
                       label=f'Matched to published proteome, DIA-NN contaminant ({n_det_cont})'),
        mpatches.Patch(facecolor=C_NOTDET_NONCONT, edgecolor='white', linewidth=0.5,
                       label=f'Not matched to published proteome, non-contaminant ({n_notdet_noncont})'),
        mpatches.Patch(facecolor=C_NOTDET_CONT, edgecolor='white', linewidth=0.5,
                       label=f'Not matched to published proteome, DIA-NN contaminant ({n_notdet_cont})'),
    ]
    leg = ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(0.01, 0.02),
              fontsize=5.5, frameon=True, fancybox=False, ncol=1, handlelength=1.0,
              handletextpad=0.5, borderpad=0.4, labelspacing=0.3)
    style_legend(leg)

    ax.set_xlim(0, n_total + 2)
    ax.set_ylim(-0.8, 1.2)
    ax.set_xlabel('Number of protein groups')
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_title('A', loc='left', pad=6, fontweight='bold')
    style_axes(ax, grid_y=False)

    # Panel B: Abundance rank plot
    ax2 = axes2[1]
    df_val = pd.DataFrame(validated_proteins)
    df_val = df_val.drop_duplicates(subset='Gene').sort_values('SPIDR_Rank')

    ax2.scatter(df_spidr['Rank'], df_spidr['Log10_Quantity'],
                s=3, color=C_LGREY, alpha=0.55, zorder=1, rasterized=True)
    ax2.scatter(df_val['SPIDR_Rank'], df_val['SPIDR_Log10Qty'],
                s=24, color=C_INDIGO, edgecolors='white', linewidth=0.4,
                alpha=0.9, zorder=3)

    label_genes = {
        'KRT18': (35, 15), 'GAPDH': (35, 5), 'HSPB1': (40, -8),
        'HSP90AA1': (-50, 15), 'EEF1A1': (-45, 18), 'ANXA2': (40, 8),
        'TRIM28': (-50, -12), 'FABP5': (40, 8), 'DSP': (35, 10),
        'JUP': (-50, 15), 'CALML5': (35, 12), 'LMNA': (-50, 22),
        'PKM': (-50, 8), 'ANXA1': (-45, 12), 'COL1A1': (30, 10),
        'DSC1': (30, 12), 'DSG1': (30, 10), 'PKP1': (-40, -12),
        'ACTB': (30, 12), 'UBB': (-50, -20), 'PRDX1': (-50, -8),
        'CTSD': (40, -10), 'CST6': (-45, -15),
    }
    for _, row in df_val.iterrows():
        gene = row['Gene']
        if gene in label_genes:
            ox, oy = label_genes[gene]
            ax2.annotate(gene, xy=(row['SPIDR_Rank'], row['SPIDR_Log10Qty']),
                         xytext=(ox, oy), textcoords='offset points',
                         fontsize=7, fontweight='bold', color='#1a1a1a',
                         arrowprops=dict(arrowstyle='-', color='#888888', lw=0.5),
                         zorder=4)

    qty_max = df_val['SPIDR_Log10Qty'].max()
    qty_min = df_val['SPIDR_Log10Qty'].min()
    fold_range = 10 ** (qty_max - qty_min)
    x_arrow = 4400
    ax2.annotate('', xy=(x_arrow, qty_min), xytext=(x_arrow, qty_max),
                 arrowprops=dict(arrowstyle='<->', color='#C0392B', lw=1.8))
    ax2.text(x_arrow + 100, (qty_max + qty_min) / 2,
             f'~{fold_range:.0e}-fold\ndynamic range',
             ha='left', va='center', fontsize=9, fontweight='bold', color='#C0392B')

    ax2.set_xlabel('Abundance rank in MCF-7 single-cell proteome\n(Wang et al. 2025, SPIDR)',
                   fontsize=11, fontweight='bold')
    ax2.set_ylabel('Mean protein quantity (log$_{10}$)', fontsize=11, fontweight='bold')
    ax2.set_xlim(-100, 4800)
    y_lo = max(-1.0, df_spidr['Log10_Quantity'].min() - 0.2)
    y_hi = df_spidr['Log10_Quantity'].max() + 0.2
    ax2.set_ylim(y_lo, y_hi)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    leg_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#E0E0E0',
               markersize=6, label='MCF-7 scMS proteome (Wang et al. 2025)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#2E86AB',
               markeredgecolor='white', markersize=8, label='BlotSeq-validated proteins'),
    ]
    leg2 = ax2.legend(handles=leg_elements, loc='upper right', fontsize=5.5,
               frameon=True, fancybox=False)
    style_legend(leg2)
    ax2.set_title('B', loc='left', pad=6, fontweight='bold')
    style_axes(ax2, grid_y=False)

    plt.tight_layout()
    for fmt in ['png', 'svg']:
        fig2.savefig(os.path.join(OUTPUT_DIR, f'Panel2_SPIDR_Validation.{fmt}'),
                     dpi=1200, bbox_inches='tight', facecolor='white')
        print(f'  Saved: Panel2_SPIDR_Validation.{fmt}')
    plt.close(fig2)
else:
    print('  SPIDR file not found — skipping Panel 2 figure.')
    print(f'  Expected: {SPIDR_FILE}')
    pg_data['Match_Status'] = 'N/A'
    n_det_noncont = n_det_cont = n_notdet_noncont = n_notdet_cont = 0
    n_total = len(pg_data)
    df_val = pd.DataFrame()


# ═══════════════════════════════════════════════════════════════════════════
# PANEL 3: PROTEIN COMPOSITION — GO:Cellular Component stacked bars
# ═══════════════════════════════════════════════════════════════════════════
print('\n--- Panel 3: Protein composition (GO:CC) ---')

# Get unique protein-gene mapping
protein_gene = (
    df_023.drop_duplicates(subset=['Protein.Group'])[['Protein.Group', 'Genes']]
    .set_index('Protein.Group')['Genes']
    .to_dict()
)
valid_genes = [
    g for g in set(protein_gene.values())
    if isinstance(g, str) and g not in ('nan', '', 'NaN') and not g.startswith('prpL')
]
print(f'  Unique genes to annotate: {len(valid_genes)}')

# Query g:Profiler for GO:CC
print('  Querying g:Profiler (GO:Cellular Component)...')
try:
    from gprofiler import GProfiler
    gp = GProfiler(return_dataframe=True)
    go_results = gp.profile(
        organism='hsapiens',
        query=valid_genes,
        sources=['GO:CC'],
        significance_threshold_method='fdr',
        user_threshold=0.05,
        no_evidences=False,
    )
    print(f'  g:Profiler returned {len(go_results)} GO:CC associations')
    gprofiler_success = True
except Exception as e:
    print(f'  g:Profiler query failed: {e}')
    print('  Falling back to gene-name-based classification only.')
    go_results = pd.DataFrame()
    gprofiler_success = False

# GO term to category mapping
GO_TO_CATEGORY = {
    'GO:0005856': 'Cytoskeletal/Structural', 'GO:0015629': 'Cytoskeletal/Structural',
    'GO:0005882': 'Cytoskeletal/Structural', 'GO:0045111': 'Cytoskeletal/Structural',
    'GO:0005883': 'Cytoskeletal/Structural', 'GO:0005875': 'Cytoskeletal/Structural',
    'GO:0005874': 'Cytoskeletal/Structural', 'GO:0030018': 'Cytoskeletal/Structural',
    'GO:0030017': 'Cytoskeletal/Structural', 'GO:0005604': 'Cytoskeletal/Structural',
    'GO:0005578': 'Cytoskeletal/Structural', 'GO:0031430': 'Cytoskeletal/Structural',
    'GO:0001725': 'Cytoskeletal/Structural', 'GO:0045095': 'Cytoskeletal/Structural',
    'GO:0005881': 'Cytoskeletal/Structural',
    'GO:0030057': 'Cytoskeletal/Structural', 'GO:0005923': 'Cytoskeletal/Structural',
    'GO:0070160': 'Cytoskeletal/Structural', 'GO:0005924': 'Cytoskeletal/Structural',
    'GO:0005634': 'Nuclear/Chromatin', 'GO:0000785': 'Nuclear/Chromatin',
    'GO:0000786': 'Nuclear/Chromatin', 'GO:0005694': 'Nuclear/Chromatin',
    'GO:0005730': 'Nuclear/Chromatin', 'GO:0005635': 'Nuclear/Chromatin',
    'GO:0005654': 'Nuclear/Chromatin',
    'GO:0005739': 'Mitochondrial', 'GO:0005740': 'Mitochondrial',
    'GO:0005743': 'Mitochondrial', 'GO:0005759': 'Mitochondrial',
    'GO:0019866': 'Mitochondrial', 'GO:0031966': 'Mitochondrial',
    'GO:0005840': 'Ribosomal', 'GO:0022626': 'Ribosomal',
    'GO:0015935': 'Ribosomal', 'GO:0015934': 'Ribosomal',
    'GO:0022625': 'Ribosomal',
    'GO:0000502': 'Proteasome', 'GO:0005838': 'Proteasome',
    'GO:0019774': 'Proteasome',
    'GO:0005615': 'Extracellular/Secreted', 'GO:0005576': 'Extracellular/Secreted',
    'GO:0070062': 'Extracellular/Secreted', 'GO:0031012': 'Extracellular/Secreted',
    'GO:0005886': 'Plasma membrane', 'GO:0009986': 'Plasma membrane',
    'GO:0005887': 'Plasma membrane',
    'GO:0005783': 'Chaperones/HSP', 'GO:0005788': 'Chaperones/HSP',
    'GO:0042571': 'Immune/Antibody', 'GO:0019814': 'Immune/Antibody',
    'GO:0044154': 'Immune/Antibody',
    'GO:0005829': 'Metabolic enzymes', 'GO:0005737': 'Metabolic enzymes',
}

CYTOPLASMIC_EXCLUSION = {
    'GO:0045095', 'GO:0045111', 'GO:0005737', 'GO:0005829',
    'GO:0015629', 'GO:0005881', 'GO:0005874',
}

CATEGORY_PRIORITY = [
    'Ribosomal', 'Proteasome', 'Mitochondrial', 'Cytoskeletal/Structural',
    'Nuclear/Chromatin', 'Chaperones/HSP', 'Immune/Antibody',
    'Extracellular/Secreted', 'Plasma membrane', 'Metabolic enzymes',
]

NUCLEAR_GENE_PATTERNS = [
    'H2AC', 'H2AJ', 'H2AX', 'H2AZ', 'H2BC', 'H2BK', 'H2BN', 'H2BS',
    'H3C', 'H3-', 'H3F', 'H4C', 'HIST1H', 'HIST2H', 'HIST3H',
    'LMNA', 'LMNB1', 'LMNB2', 'TRIM28', 'TRIM33',
    'NUP98', 'NUP153', 'NUP214', 'NUP62',
]

def is_nuclear_by_name(gene):
    g = gene.upper().strip()
    return any(g == p or g.startswith(p) for p in NUCLEAR_GENE_PATTERNS)

def assign_category_from_go(gene_go_terms):
    matched = set()
    has_cyto = bool(gene_go_terms & CYTOPLASMIC_EXCLUSION)
    for tid in gene_go_terms:
        if tid in GO_TO_CATEGORY:
            cat = GO_TO_CATEGORY[tid]
            if cat == 'Nuclear/Chromatin' and has_cyto:
                continue
            matched.add(cat)
    for cat in CATEGORY_PRIORITY:
        if cat in matched:
            return cat
    return list(matched)[0] if matched else None

# Build gene to GO terms
gene_to_go_terms = defaultdict(set)
if gprofiler_success and len(go_results) > 0 and 'intersections' in go_results.columns:
    for _, row in go_results.iterrows():
        term_id = row.get('native', '')
        intersections = row.get('intersections', [])
        if isinstance(intersections, list):
            for gene in intersections:
                gene_to_go_terms[gene.upper()].add(term_id)
    print(f'  Genes with GO:CC annotations: {len(gene_to_go_terms)}')

# MyGene.info fallback
unannotated = [g for g in valid_genes if g.upper() not in gene_to_go_terms]
if unannotated:
    print(f'  Querying MyGene.info for {len(unannotated)} unannotated genes...')
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        mg_results = mg.querymany(
            unannotated, scopes='symbol', fields='go.CC',
            species='human', returnall=False, as_dataframe=True,
        )
        for gene, row in mg_results.iterrows():
            go_cc = row.get('go.CC', None)
            if isinstance(go_cc, list):
                for entry in go_cc:
                    if isinstance(entry, dict) and 'id' in entry:
                        gene_to_go_terms[str(gene).upper()].add(entry['id'])
            elif isinstance(go_cc, dict) and 'id' in go_cc:
                gene_to_go_terms[str(gene).upper()].add(go_cc['id'])
        print(f'  MyGene.info added annotations')
    except Exception as e:
        print(f'  MyGene.info fallback failed: {e}')

# Classify each protein
def classify_protein(gene):
    if not isinstance(gene, str) or gene.strip() in ('nan', '', 'NaN'):
        return 'Uncharacterized'
    if is_nuclear_by_name(gene):
        return 'Nuclear/Chromatin'
    go_terms = gene_to_go_terms.get(gene.upper().strip(), set())
    cat = assign_category_from_go(go_terms)
    if cat:
        return cat
    if gene.startswith('Cont_'):
        return 'Contaminant'
    return 'Uncharacterized'

df_023['Category'] = df_023['Genes'].apply(classify_protein)
df_023['Condition'] = df_023['Run'].apply(get_condition)

# Aggregate per condition
condition_protein = df_023.drop_duplicates(subset=['Condition', 'Protein.Group'])
counts_by_cond = condition_protein.groupby(['Condition', 'Category']).size().unstack(fill_value=0)
totals_by_cond = counts_by_cond.sum(axis=1)
pct_by_cond = counts_by_cond.div(totals_by_cond, axis=0) * 100

cat_order = [
    'Cytoskeletal/Structural', 'Extracellular/Secreted', 'Plasma membrane',
    'Nuclear/Chromatin', 'Metabolic enzymes', 'Mitochondrial', 'Ribosomal',
    'Proteasome', 'Chaperones/HSP', 'Immune/Antibody', 'Contaminant', 'Uncharacterized',
]
cat_colors = {
    'Cytoskeletal/Structural': C_ORANGE,
    'Extracellular/Secreted': C_MAROON,
    'Plasma membrane': C_TEAL2,
    'Nuclear/Chromatin': C_INDIGO,
    'Metabolic enzymes': C_PURPLE,
    'Mitochondrial': C_GREEN,
    'Ribosomal': C_GREEN2,
    'Proteasome': C_GREY,
    'Chaperones/HSP': C_YELLOW,
    'Immune/Antibody': '#D8B24A',
    'Contaminant': '#8B8B8B',
    'Uncharacterized': C_LGREY,
}

for cat in cat_order:
    if cat not in pct_by_cond.columns:
        pct_by_cond[cat] = 0
pct_by_cond = pct_by_cond[cat_order]

cond_order = ['Blotseq 1hr', 'Blotseq 2hr', 'Blotseq 3hr']
pct_by_cond = pct_by_cond.reindex(cond_order)

fig3, ax3 = plt.subplots(figsize=(5.1, 5.2))

x = np.arange(len(cond_order))
bar_width = 0.55
bottoms = np.zeros(len(cond_order))

for cat in cat_order:
    values = pct_by_cond[cat].values
    ax3.bar(x, values, bar_width, bottom=bottoms,
            color=cat_colors[cat], edgecolor='none', linewidth=0, label=cat, zorder=3)
    bottoms += values

for i, cond in enumerate(cond_order):
    ax3.text(i, 101.0, f'n={int(totals_by_cond[cond])}', ha='center', va='bottom',
             fontsize=5.8, color=C_TEXT)

ax3.set_xticks(x)
ax3.set_xticklabels(['1 hr', '2 hr', '3 hr'])
ax3.set_ylabel('% of detected proteins')
ax3.set_title('Protein categories (GO: Cellular Component)', pad=8)
ax3.set_ylim(0, 108)
ax3.set_yticks(range(0, 101, 20))
style_axes(ax3, grid_y=True)

handles = [mpatches.Patch(color=cat_colors[c], label=c) for c in cat_order]
leg3 = ax3.legend(
    handles=handles,
    loc='upper center',
    bbox_to_anchor=(0.5, -0.18),
    ncol=4,
    fontsize=5.2,
    frameon=True,
    fancybox=False,
    title='GO:CC category',
    title_fontsize=5.5,
    borderpad=0.4,
    labelspacing=0.25,
    handlelength=0.9,
    columnspacing=0.8
)
style_legend(leg3)

plt.tight_layout(rect=[0, 0.22, 1, 1])

fig3.savefig(
    os.path.join(OUTPUT_DIR, 'Panel3_ProteinCategories_GO_bottomLegend.png'),
    dpi=1200,
    bbox_inches='tight',
    facecolor='white'
)
plt.close(fig3)
# ═══════════════════════════════════════════════════════════════════════════
# EXCEL OUTPUT — Multi-sheet workbook
# ═══════════════════════════════════════════════════════════════════════════
print('\n--- Generating Excel workbook ---')

with pd.ExcelWriter(EXCEL_OUT, engine='openpyxl') as writer:

    # Sheet 1: Bar chart data (Panel 1)
    df_bar = protein_counts[['Sample', 'Condition', 'Proteins_Identified', 'HSP90AA1_detected']].copy()
    df_bar.to_excel(writer, sheet_name='Panel1_BarChart', index=False)
    print('  Sheet: Panel1_BarChart')

    # Sheet 2: SPIDR validation data (Panel 2)
    if spidr_available:
        # Summary counts
        spidr_summary = pd.DataFrame([
            {'Category': 'Detected in Wang et al., non-contaminant', 'Count': n_det_noncont},
            {'Category': 'Detected in Wang et al., DIA-NN contaminant', 'Count': n_det_cont},
            {'Category': 'Not detected, non-contaminant', 'Count': n_notdet_noncont},
            {'Category': 'Not detected, DIA-NN contaminant', 'Count': n_notdet_cont},
            {'Category': 'TOTAL', 'Count': n_total},
        ])
        spidr_summary.to_excel(writer, sheet_name='Panel2_SPIDR_Summary', index=False)

        # Per-protein detail
        pg_export = pg_data[['Protein.Group', 'Genes', 'Protein.Names',
                             'Is_Contaminant', 'Match_Status', 'N_Samples_Detected',
                             'PG.MaxLFQ']].copy()
        pg_export.to_excel(writer, sheet_name='Panel2_SPIDR_Detail', index=False)

        # Validated proteins with SPIDR ranks
        if len(df_val) > 0:
            df_val.to_excel(writer, sheet_name='Panel2_Validated_Proteins', index=False)
        print('  Sheets: Panel2_SPIDR_Summary, Panel2_SPIDR_Detail, Panel2_Validated_Proteins')
    else:
        pg_data.to_excel(writer, sheet_name='Panel2_ProteinGroups', index=False)
        print('  Sheet: Panel2_ProteinGroups (no SPIDR data)')

    # Sheet 3: Protein categories data (Panel 3)
    # Raw counts per condition x category
    counts_by_cond.to_excel(writer, sheet_name='Panel3_GO_Counts')
    pct_by_cond.to_excel(writer, sheet_name='Panel3_GO_Percent')

    # Per-protein annotation table
    annotation_table = df_023.drop_duplicates(subset=['Protein.Group'])[
        ['Protein.Group', 'Genes', 'Category']
    ].sort_values('Category')
    annotation_table.to_excel(writer, sheet_name='Panel3_GO_Annotations', index=False)
    print('  Sheets: Panel3_GO_Counts, Panel3_GO_Percent, Panel3_GO_Annotations')

    # Sheet 4: Blank sample summary
    blank_summary_data = {
        'Metric': [
            'Sample', 'Experiment', 'Precursor rows (DIA-NN)',
            'Nominal protein groups (all fail FDR)',
            'Protein groups with Cont_ prefix',
            'Protein groups without Cont_ prefix',
            'Global.PG.Q.Value (all entries)',
            'Proteins passing 1% FDR',
            'Precursors.Identified (from stats file)',
            'Proteins.Identified (from stats file)',
        ],
        'Value': [
            'BlotseqBlank_0cell', 'AMHE025',
            str(len(df_blank)), str(blank_precursors),
            str(df_blank[df_blank['Protein.Group'].str.startswith('Cont_')]['Protein.Group'].nunique()),
            str(df_blank[~df_blank['Protein.Group'].str.startswith('Cont_')]['Protein.Group'].nunique()),
            '1.0 (100% FDR — all entries)',
            '0',
            '900',
            '0',
        ]
    }
    pd.DataFrame(blank_summary_data).to_excel(writer, sheet_name='Blank_Summary', index=False)

    # Top proteins in blank by PG.Quantity (for reference)
    blank_top = (
        df_025[df_025['Run'].str.contains('BlotseqBlank')]
        .groupby('Protein.Group')
        .agg({'Genes': 'first'})
        .reset_index()
    )
    blank_top.columns = ['Protein.Group', 'Genes']
    blank_top['Is_Contaminant'] = blank_top['Protein.Group'].str.startswith('Cont_')
    blank_top.to_excel(writer, sheet_name='Blank_ProteinGroups', index=False)
    print('  Sheets: Blank_Summary, Blank_ProteinGroups')

print(f'\nExcel saved: {EXCEL_OUT}')

# ═══════════════════════════════════════════════════════════════════════════
# ADDITIONAL SECTION: revised Panel B and Panel E
# ═══════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════
# PUBLICATION STYLE — matches Fig_1_v3 / Fig_2_v3 (Herr Lab style)
# ═══════════════════════════════════════════════════════════════════════════
# White background, minimal spines, Helvetica-like sans-serif
# Colors: burnt orange #D2613C, indigo #3C3C96, teal #2D6B5A, 
#         grey #A0A0A0, dark text #2D2D2D

import matplotlib
matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 7,
    'axes.labelsize': 8,
    'axes.titlesize': 9,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'savefig.dpi': 1200,
    'pdf.fonttype': 42,       # editable text in Illustrator
    'ps.fonttype': 42,
})

# Color palette (extracted from published figures)
C_ORANGE  = '#D25A3C'    # primary warm — bars, highlights
C_INDIGO  = '#3C3C96'    # primary cool — bars, data
C_TEAL    = '#2D695A'    # accent
C_PURPLE  = '#3C0F4B'    # accent
C_GREY    = '#A0A0A0'    # neutral / contaminant
C_LGREY   = '#D2D2D2'    # light neutral
C_TEXT    = '#2D2D2D'     # text and axes
C_RED     = '#781E1E'     # accent / marker
C_BLUE_LT = '#8CA0B4'    # light blue accent

# Timepoint palette (consistent with orange = warm = longer digest)
TP_COL = {
    '1hr': '#8CA0B4',     # muted blue-grey (shortest)
    '2hr': C_ORANGE,      # burnt orange (mid)
    '3hr': C_INDIGO,      # deep indigo (longest)
}
CONT_COL   = '#D0D0D0'   # very light grey for contaminant bars
BLANK_COL  = '#E8E8E8'
HSP_COL    = C_RED


# ═══════════════════════════════════════════════════════════════════════════
# PANEL B (revised): Contaminant-deconvoluted protein depth per cell
# ═══════════════════════════════════════════════════════════════════════════
print('\n--- Panel B (revised): Contaminant-deconvoluted strip plot ---')

# Per-cell counts split by contaminant status
df_023['_is_cont'] = df_023['Protein.Group'].str.startswith('Cont_')

cell_data_b = []
for run in df_023['Run'].unique():
    sub = df_023[df_023['Run'] == run]
    n_total = sub['Protein.Group'].nunique()
    n_cont = sub[sub['_is_cont']]['Protein.Group'].nunique()
    n_noncont = n_total - n_cont
    tp_label = get_condition(run).replace('Blotseq ', '')
    cell_data_b.append({
        'Run': run,
        'Sample': short_name_amhe023(run),
        'Timepoint': tp_label,
        'Proteins_Total': n_total,
        'Proteins_NonCont': n_noncont,
        'Proteins_Cont': n_cont,
        'HSP90AA1': run in hsp_runs,
    })

df_cellB = pd.DataFrame(cell_data_b)
# Add blank
df_cellB = pd.concat([df_cellB, pd.DataFrame([{
    'Run': 'BlotseqBlank_0cell', 'Sample': 'Blank',
    'Timepoint': 'Blank', 'Proteins_Total': 0,
    'Proteins_NonCont': 0, 'Proteins_Cont': 0, 'HSP90AA1': False,
}])], ignore_index=True)

# Sort within each timepoint
tp_sort = {'1hr': 0, '2hr': 1, '3hr': 2, 'Blank': 3}
df_cellB['_sort'] = df_cellB['Timepoint'].map(tp_sort)
df_cellB = df_cellB.sort_values(['_sort', 'Proteins_Total'],
                                 ascending=[True, False]).reset_index(drop=True)

# ── Figure generation helper ──
tp_x = {'1hr': 0, '2hr': 1.3, '3hr': 2.6, 'Blank': 4.0}
bar_w = 0.16

def make_panel_b(show_hsp90_stars=True, basename='PanelB_ContDeconvoluted'):
    fig_b, ax_b = plt.subplots(figsize=(4.2, 3.0))

    np.random.seed(42)
    for tp in ['1hr', '2hr', '3hr']:
        sub = df_cellB[df_cellB['Timepoint'] == tp].reset_index(drop=True)
        x0 = tp_x[tp]
        n = len(sub)
        offsets = np.linspace(-0.26, 0.26, n)

        for i, (_, row) in enumerate(sub.iterrows()):
            xi = x0 + offsets[i]
            ax_b.bar(xi, row.Proteins_NonCont, width=bar_w, bottom=0,
                     color=TP_COL[tp], edgecolor='none', zorder=3, alpha=0.85)
            ax_b.bar(xi, row.Proteins_Cont, width=bar_w, bottom=row.Proteins_NonCont,
                     color=CONT_COL, edgecolor='none', zorder=3)
            if show_hsp90_stars and row.HSP90AA1:
                ax_b.plot(xi, row.Proteins_Total + 3.5, marker='*',
                          color=HSP_COL, markersize=5, markeredgecolor='none', zorder=5)

        med_t = sub['Proteins_Total'].median()
        med_nc = sub['Proteins_NonCont'].median()
        span = 0.36
        ax_b.plot([x0 - span, x0 + span], [med_t, med_t],
                  color=C_TEXT, linewidth=1.0, linestyle='-', zorder=4)
        ax_b.plot([x0 - span, x0 + span], [med_nc, med_nc],
                  color=TP_COL[tp], linewidth=0.8, linestyle='--', zorder=4, alpha=0.7)

    ax_b.bar(tp_x['Blank'], 0, width=0.25, color=BLANK_COL,
             edgecolor=C_GREY, linewidth=0.4, zorder=3)
    ax_b.text(tp_x['Blank'], 2.5, '0', ha='center', va='bottom',
              fontsize=6, color=C_GREY, fontweight='bold')

    from matplotlib.lines import Line2D
    leg_h = [
        mpatches.Patch(facecolor=C_ORANGE, edgecolor='none', alpha=0.85,
                       label='Non-contaminant'),
        mpatches.Patch(facecolor=CONT_COL, edgecolor='none',
                       label='Cont_-flagged'),
        Line2D([0], [0], color=C_TEXT, linewidth=1.0, label='Median (total)'),
        Line2D([0], [0], color=C_ORANGE, linewidth=0.8, linestyle='--',
               alpha=0.7, label='Median (non-cont)'),
    ]
    if show_hsp90_stars:
        leg_h.append(Line2D([0], [0], marker='*', color='w', markerfacecolor=HSP_COL,
                            markeredgecolor='none', markersize=6, label='HSP90AA1 detected'))
    leg = ax_b.legend(handles=leg_h, loc='upper right', frameon=True,
                      edgecolor=C_LGREY, fancybox=False,
                      handlelength=1.0, labelspacing=0.3, borderpad=0.4,
                      fontsize=5.5)
    leg.get_frame().set_linewidth(0.4)

    ax_b.set_xticks([tp_x[t] for t in ['1hr', '2hr', '3hr', 'Blank']])
    ax_b.set_xticklabels(['1 hr\n(n=5)', '2 hr\n(n=5)', '3 hr\n(n=5)', 'Blank\n(n=1)'],
                         color=C_TEXT)
    ax_b.set_ylabel('Protein groups per cell', color=C_TEXT)
    ax_b.set_xlabel('In-gel digestion time', color=C_TEXT)
    ax_b.set_ylim(0, 155)
    ax_b.tick_params(colors=C_TEXT)
    ax_b.yaxis.grid(True, color='#ECECEC', linewidth=0.3, zorder=0)
    ax_b.set_axisbelow(True)
    for spine in ['left', 'bottom']:
        ax_b.spines[spine].set_color(C_TEXT)

    fig_b.tight_layout()
    for fmt in ['png', 'svg']:
        fig_b.savefig(os.path.join(OUTPUT_DIR, f'{basename}.{fmt}'),
                      dpi=1200, bbox_inches='tight', facecolor='white')
        print(f'  Saved: {basename}.{fmt}')
    plt.close(fig_b)

# Generate both versions: with and without HSP90 stars
make_panel_b(show_hsp90_stars=True, basename='PanelB_ContDeconvoluted')
make_panel_b(show_hsp90_stars=False, basename='PanelB_ContDeconvoluted_noHSP90')


# ═══════════════════════════════════════════════════════════════════════════
# PANEL B (additional): Aggregated mean ± SD bar plot by digestion time
# ═══════════════════════════════════════════════════════════════════════════
print('\n--- Panel B (additional): Aggregated mean ± SD bar plot ---')
# Summary statistics for aggregated panel
panel_b_agg_rows = []
for tp in ['1hr', '2hr', '3hr']:
    sub = df_cellB[df_cellB['Timepoint'] == tp].copy()
    panel_b_agg_rows.append({
        'Timepoint': tp,
        'n_cells': len(sub),
        'Mean_Total': sub['Proteins_Total'].mean(),
        'SD_Total': sub['Proteins_Total'].std(ddof=1),
        'Values_Total': sub['Proteins_Total'].tolist(),
    })

# Blank as separate category
sub_blank = df_cellB[df_cellB['Timepoint'] == 'Blank'].copy()
panel_b_agg_rows.append({
    'Timepoint': 'Blank',
    'n_cells': len(sub_blank),
    'Mean_Total': sub_blank['Proteins_Total'].mean() if len(sub_blank) > 0 else 0,
    'SD_Total': 0,
    'Values_Total': sub_blank['Proteins_Total'].tolist(),
})

panel_b_agg_df = pd.DataFrame(panel_b_agg_rows)

fig_bagg, ax_bagg = plt.subplots(figsize=(4.1, 3.2))

agg_order = ['1hr', '2hr', '3hr', 'Blank']
agg_labels = ['1 hr', '2 hr', '3 hr', 'Blank']
agg_colors = [TP_COL['1hr'], TP_COL['2hr'], TP_COL['3hr'], BLANK_COL]
means = [panel_b_agg_df.loc[panel_b_agg_df['Timepoint'] == tp, 'Mean_Total'].iloc[0] for tp in agg_order]
sds = [panel_b_agg_df.loc[panel_b_agg_df['Timepoint'] == tp, 'SD_Total'].iloc[0] for tp in agg_order]
x = np.arange(len(agg_order))

bars = ax_bagg.bar(
    x, means, width=0.6, color=agg_colors, edgecolor='none', linewidth=0, zorder=3,
    yerr=sds, capsize=4, ecolor=C_TEXT, error_kw=dict(linewidth=0.8, capthick=0.8)
)

# Overlay individual data points for the 5 cells in each digestion time
np.random.seed(11)
for i, tp in enumerate(agg_order):
    vals = panel_b_agg_df.loc[panel_b_agg_df['Timepoint'] == tp, 'Values_Total'].iloc[0]
    vals = vals if isinstance(vals, list) else []
    if len(vals) == 0:
        continue
    if len(vals) == 1:
        offsets = [0.0]
    else:
        offsets = np.linspace(-0.08, 0.08, len(vals))
    for off, val in zip(offsets, vals):
        ax_bagg.scatter(i + off, val, s=11, color=C_TEXT, zorder=4)

ax_bagg.set_xticks(x)
ax_bagg.set_xticklabels(agg_labels, color=C_TEXT)
ax_bagg.set_ylabel('Protein groups per cell', color=C_TEXT)
ax_bagg.set_xlabel('In-gel digestion time', color=C_TEXT)
ax_bagg.tick_params(colors=C_TEXT)
ax_bagg.yaxis.grid(True, color='#ECECEC', linewidth=0.3, zorder=0)
ax_bagg.set_axisbelow(True)
for spine in ['left', 'bottom']:
    ax_bagg.spines[spine].set_color(C_TEXT)

# Add a little headroom
max_y = max([m + s for m, s in zip(means, sds)] + [1])
ax_bagg.set_ylim(0, max_y * 1.20)

fig_bagg.tight_layout()
for fmt in ['png', 'svg']:
    fig_bagg.savefig(os.path.join(OUTPUT_DIR, f'PanelB_AggregatedMeanSD.{fmt}'),
                     dpi=1200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: PanelB_AggregatedMeanSD.{fmt}')
plt.close(fig_bagg)


# ═══════════════════════════════════════════════════════════════════════════
# PANEL E (revised): Detection frequency histogram
# ═══════════════════════════════════════════════════════════════════════════
print('\n--- Panel E (revised): Detection frequency histogram ---')

pg_freq_e = df_023.groupby('Protein.Group').agg(
    Gene=('Genes', 'first'),
    N_Cells=('Run', 'nunique'),
    Mean_MaxLFQ=('PG.MaxLFQ', 'mean'),
).reset_index()
pg_freq_e['Is_Cont'] = pg_freq_e['Protein.Group'].str.startswith('Cont_')

bins_e = np.arange(1, 17)
nc_vals = pg_freq_e[~pg_freq_e['Is_Cont']]['N_Cells'].values
cont_vals = pg_freq_e[pg_freq_e['Is_Cont']]['N_Cells'].values
nc_hist, _ = np.histogram(nc_vals, bins=bins_e)
cont_hist, _ = np.histogram(cont_vals, bins=bins_e)

total_nc = (~pg_freq_e['Is_Cont']).sum()
total_cont = pg_freq_e['Is_Cont'].sum()
nc_all15 = (pg_freq_e[~pg_freq_e['Is_Cont']]['N_Cells'] == 15).sum()
nc_10plus = (pg_freq_e[~pg_freq_e['Is_Cont']]['N_Cells'] >= 10).sum()

# ── Figure ──
fig_e, ax_e = plt.subplots(figsize=(4.2, 3.0))

x_pos = np.arange(1, 16)
off = 0.18

bars_nc = ax_e.bar(x_pos - off, nc_hist, width=0.34, color=C_INDIGO,
                    edgecolor='none', zorder=3, alpha=0.85, label='Non-contaminant')
bars_ct = ax_e.bar(x_pos + off, cont_hist, width=0.34, color=CONT_COL,
                    edgecolor='none', zorder=3, label='Cont_-flagged')

# Individual data points overlaid (matching the style from their box plots)
# Jitter non-cont protein counts
np.random.seed(7)
for val, is_c in zip(pg_freq_e['N_Cells'].values, pg_freq_e['Is_Cont'].values):
    jit = np.random.uniform(-0.06, 0.06)
    x = val - off + jit if not is_c else val + off + jit
    y_base = nc_hist[val-1] if not is_c else cont_hist[val-1]
    # Only overlay dots where bar exists
    color = C_INDIGO if not is_c else C_GREY
    # Skip overlaying dots - it clutters at this scale

# Combined legend + summary box
from matplotlib.offsetbox import AnchoredOffsetbox, VPacker, HPacker, TextArea, DrawingArea
from matplotlib.patches import Rectangle

# Build compact legend rows inside a single annotation box
sq1 = DrawingArea(8, 8, 0, 0)
sq1.add_artist(Rectangle((0, 0), 8, 8, fc=C_INDIGO, ec='none', alpha=0.85))
legend_row1 = HPacker(children=[
    sq1,
    TextArea(' Non-contaminant', textprops=dict(color=C_TEXT, fontsize=5.5))
], align='center', pad=0, sep=2)

sq2 = DrawingArea(8, 8, 0, 0)
sq2.add_artist(Rectangle((0, 0), 8, 8, fc=CONT_COL, ec='none'))
legend_row2 = HPacker(children=[
    sq2,
    TextArea(' Cont_-flagged', textprops=dict(color=C_TEXT, fontsize=5.5))
], align='center', pad=0, sep=2)

summary_box = TextArea(
    f'{total_nc} non-cont. proteins\n{nc_all15} in all 15 cells\n{nc_10plus} in ≥10 cells',
    textprops=dict(color=C_TEXT, fontsize=5.5)
)

combined_box = VPacker(children=[legend_row1, legend_row2, summary_box],
                       align='left', pad=0, sep=2)
anch = AnchoredOffsetbox(loc='upper left', child=combined_box, frameon=True,
                         bbox_to_anchor=(0.015, 0.985), bbox_transform=ax_e.transAxes,
                         borderpad=0.35, pad=0.2)
anch.patch.set_facecolor('white')
anch.patch.set_edgecolor(C_LGREY)
anch.patch.set_linewidth(0.4)
anch.patch.set_alpha(0.92)
ax_e.add_artist(anch)

ax_e.set_xlabel('Number of cells detected in (of 15)', color=C_TEXT)
ax_e.set_ylabel('Number of protein groups', color=C_TEXT)
ax_e.set_xticks(range(1, 16))
ax_e.tick_params(colors=C_TEXT)
ax_e.yaxis.grid(True, color='#ECECEC', linewidth=0.3, zorder=0)
ax_e.set_axisbelow(True)
for spine in ['left', 'bottom']:
    ax_e.spines[spine].set_color(C_TEXT)

fig_e.tight_layout()
for fmt in ['png', 'svg']:
    fig_e.savefig(os.path.join(OUTPUT_DIR, f'PanelE_DetectionFrequency.{fmt}'),
                  dpi=1200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: PanelE_DetectionFrequency.{fmt}')


# ── Additional Figure: Panel E without flagged contaminants ──
fig_e_nc, ax_e_nc = plt.subplots(figsize=(4.2, 3.0))

x_pos_nc = np.arange(1, 16)
ax_e_nc.bar(x_pos_nc, nc_hist, width=0.42, color=C_INDIGO,
            edgecolor='none', zorder=3, alpha=0.85)



ax_e_nc.set_xlabel('Number of cells detected in (of 15)', color=C_TEXT)
ax_e_nc.set_ylabel('Number of protein groups', color=C_TEXT)
ax_e_nc.set_xticks(range(1, 16))
ax_e_nc.tick_params(colors=C_TEXT)
ax_e_nc.yaxis.grid(True, color='#ECECEC', linewidth=0.3, zorder=0)
ax_e_nc.set_axisbelow(True)
for spine in ['left', 'bottom']:
    ax_e_nc.spines[spine].set_color(C_TEXT)

fig_e_nc.tight_layout()
for fmt in ['png', 'svg']:
    fig_e_nc.savefig(os.path.join(OUTPUT_DIR, f'PanelE_DetectionFrequency_noContaminants.{fmt}'),
                     dpi=1200, bbox_inches='tight', facecolor='white')
    print(f'  Saved: PanelE_DetectionFrequency_noContaminants.{fmt}')
plt.close(fig_e_nc)


# ═══════════════════════════════════════════════════════════════════════════
# EXCEL SHEETS for Panel B and Panel E
# ═══════════════════════════════════════════════════════════════════════════
print('\n--- Saving Excel data for Panels B and E ---')

panel_excel = os.path.join(OUTPUT_DIR, 'PanelBE_data.xlsx')
with pd.ExcelWriter(panel_excel, engine='openpyxl') as writer:
    # Panel B: per-cell data
    df_cellB_out = df_cellB[['Sample', 'Timepoint', 'Proteins_Total',
                              'Proteins_NonCont', 'Proteins_Cont', 'HSP90AA1']].copy()
    df_cellB_out.to_excel(writer, sheet_name='PanelB_PerCell', index=False)

    # Panel B: summary statistics
    b_summary = []
    for tp in ['1hr', '2hr', '3hr']:
        sub = df_cellB[df_cellB['Timepoint'] == tp]
        b_summary.append({
            'Timepoint': tp,
            'n_cells': len(sub),
            'Median_Total': sub['Proteins_Total'].median(),
            'Mean_Total': sub['Proteins_Total'].mean(),
            'SD_Total': sub['Proteins_Total'].std(),
            'Median_NonCont': sub['Proteins_NonCont'].median(),
            'Mean_NonCont': sub['Proteins_NonCont'].mean(),
            'SD_NonCont': sub['Proteins_NonCont'].std(),
            'Median_Cont': sub['Proteins_Cont'].median(),
            'HSP90_positive': sub['HSP90AA1'].sum(),
        })
    pd.DataFrame(b_summary).to_excel(writer, sheet_name='PanelB_Summary', index=False)

    # Panel B additional aggregated panel data
    panel_b_agg_export = panel_b_agg_df[['Timepoint', 'n_cells', 'Mean_Total', 'SD_Total']].copy()
    panel_b_agg_export.to_excel(writer, sheet_name='PanelB_Aggregated', index=False)

    panel_b_agg_points = []
    for tp in ['1hr', '2hr', '3hr', 'Blank']:
        sub = df_cellB[df_cellB['Timepoint'] == tp].copy()
        for _, row in sub.iterrows():
            panel_b_agg_points.append({
                'Timepoint': tp,
                'Sample': row['Sample'],
                'Proteins_Total': row['Proteins_Total']
            })
    pd.DataFrame(panel_b_agg_points).to_excel(writer, sheet_name='PanelB_Aggregated_Points', index=False)

    # Panel E: per-protein detection frequency
    pg_freq_e_out = pg_freq_e[['Protein.Group', 'Gene', 'N_Cells',
                                'Mean_MaxLFQ', 'Is_Cont']].sort_values(
        ['Is_Cont', 'N_Cells'], ascending=[True, False])
    pg_freq_e_out.to_excel(writer, sheet_name='PanelE_DetectionFreq', index=False)

    # Panel E: histogram counts
    hist_df = pd.DataFrame({
        'N_Cells': range(1, 16),
        'NonCont_Count': nc_hist,
        'Cont_Count': cont_hist,
        'Total_Count': nc_hist + cont_hist,
    })
    hist_df.to_excel(writer, sheet_name='PanelE_HistogramCounts', index=False)

    # Panel E additional: non-contaminant-only histogram counts
    hist_df_nc = pd.DataFrame({
        'N_Cells': range(1, 16),
        'NonCont_Count': nc_hist,
    })
    hist_df_nc.to_excel(writer, sheet_name='PanelE_NonContOnly', index=False)

    # Panel E: summary
    e_summary = pd.DataFrame([{
        'Total_ProteinGroups': len(pg_freq_e),
        'NonContaminant': total_nc,
        'Contaminant': total_cont,
        'NonCont_in_all_15': nc_all15,
        'NonCont_in_10plus': nc_10plus,
    }])
    e_summary.to_excel(writer, sheet_name='PanelE_Summary', index=False)

print(f'  Saved: {panel_excel}')
print('\nPanels B and E complete. Figures saved as PanelB_ContDeconvoluted, PanelB_ContDeconvoluted_noHSP90, PanelB_AggregatedMeanSD, PanelE_DetectionFrequency, and PanelE_DetectionFrequency_noContaminants.')

print('\nDone. All panels and Excel generated.')




# ═══════════════════════════════════════════════════════════════════════════
# PANEL E (Figure 3): Combined SPIDR + Geiger validation — single stacked bar
# Append after the last line of BlotSeq_MCF7_MS_analysis_styled.py
# ═══════════════════════════════════════════════════════════════════════════
print('\n--- Panel E: Combined SPIDR + Geiger validation ---')

GEIGER_FILE = os.path.join(INPUT_DIR, 'Geiger2012_MCF7_Proteins_Ranked.xlsx')

geiger_available = os.path.exists(GEIGER_FILE)
if not geiger_available:
    print(f'  WARNING: {GEIGER_FILE} not found — skipping combined validation')

if spidr_available and geiger_available:

    # ── Build Geiger gene set ──
    print('  Loading Geiger 2012 MCF-7 bulk proteome...')
    df_geiger = pd.read_excel(GEIGER_FILE)
    geiger_genes = set()
    for _, row in df_geiger.iterrows():
        pg = row.get('Primary Gene', '')
        if pd.notna(pg):
            geiger_genes.add(str(pg).strip().upper())
        allg = row.get('All Gene Names', '')
        if pd.notna(allg):
            for g in str(allg).split(';'):
                g = g.strip().upper()
                if g:
                    geiger_genes.add(g)
    print(f'  Geiger gene symbols: {len(geiger_genes)}')

    # ── Build SPIDR gene set (reuse spidr_gene_to_idx from earlier) ──
    spidr_genes_set = set(spidr_gene_to_idx.keys())

    # ── Classify each BlotSeq protein against combined reference ──
    combined_ref = spidr_genes_set | geiger_genes
    match_combined = []
    for _, row in pg_data.iterrows():
        gene = row['Genes']
        gu = gene.strip().upper() if pd.notna(gene) else ''
        if gu and gu in combined_ref:
            match_combined.append('Matched')
        else:
            match_combined.append('Not matched')
    pg_data['Combined_Match'] = match_combined

    c_det_nc  = ((pg_data['Combined_Match'] == 'Matched') & (~pg_data['Is_Contaminant'])).sum()
    c_det_c   = ((pg_data['Combined_Match'] == 'Matched') & (pg_data['Is_Contaminant'])).sum()
    c_nd_nc   = ((pg_data['Combined_Match'] == 'Not matched') & (~pg_data['Is_Contaminant'])).sum()
    c_nd_c    = ((pg_data['Combined_Match'] == 'Not matched') & (pg_data['Is_Contaminant'])).sum()
    c_det     = c_det_nc + c_det_c
    c_nd      = c_nd_nc + c_nd_c
    c_total   = len(pg_data)

    print(f'  Combined matched: {c_det}/{c_total} (non-cont: {c_det_nc}, cont: {c_det_c})')
    print(f'  Not matched:      {c_nd}/{c_total} (non-cont: {c_nd_nc}, cont: {c_nd_c})')

    # ── Colors (same palette as existing panels) ──
    C_M_NC  = C_INDIGO     # matched, non-contaminant
    C_M_C   = '#8585C2'    # matched, contaminant (lighter indigo)
    C_NM_NC = C_ORANGE     # not matched, non-contaminant
    C_NM_C  = '#E8A990'    # not matched, contaminant (lighter orange)

    # ── Figure ──
    fig_e, ax_e = plt.subplots(figsize=(6.5, 3.0))

    segments = [
        (c_det_nc,  C_M_NC),
        (c_det_c,   C_M_C),
        (c_nd_nc,   C_NM_NC),
        (c_nd_c,    C_NM_C),
    ]

    bar_h = 0.45
    left = 0
    for count, color in segments:
        if count == 0:
            left += count
            continue
        ax_e.barh(0, count, left=left, height=bar_h,
                  color=color, edgecolor='white', linewidth=0.8)
        if count >= 4:
            tc = 'white' if color in [C_M_NC] else C_TEXT
            ax_e.text(left + count / 2, 0, str(count),
                      ha='center', va='center', fontsize=10,
                      fontweight='bold', color=tc)
        left += count

    # ── Bracket: matched ──
    by = bar_h / 2 + 0.10
    ax_e.plot([0, c_det], [by + 0.07, by + 0.07], color=C_INDIGO, lw=1.5)
    ax_e.plot([0, 0], [by + 0.03, by + 0.11], color=C_INDIGO, lw=1.5)
    ax_e.plot([c_det, c_det], [by + 0.03, by + 0.11], color=C_INDIGO, lw=1.5)
    ax_e.text(c_det / 2, by + 0.16,
              f'Matched to published MCF-7 single-cell\nand bulk proteomes\n({c_det}/{c_total})',
              ha='center', va='bottom', fontsize=7, fontweight='bold', color=C_INDIGO)

    # ── Bracket: not matched ──
    ax_e.plot([c_det, c_total], [by + 0.07, by + 0.07], color=C_ORANGE, lw=1.5)
    ax_e.plot([c_det, c_det], [by + 0.03, by + 0.11], color=C_ORANGE, lw=1.5)
    ax_e.plot([c_total, c_total], [by + 0.03, by + 0.11], color=C_ORANGE, lw=1.5)
    ax_e.text(c_det + c_nd / 2, by + 0.16,
              f'Not matched to published\nMCF-7 proteomes\n({c_nd}/{c_total})',
              ha='center', va='bottom', fontsize=7, fontweight='bold', color=C_ORANGE)

    # ── Legend below ──
    leg_elements = [
        mpatches.Patch(facecolor=C_M_NC, edgecolor='white', lw=0.5,
                       label=f'Matched, non-contaminant ({c_det_nc})'),
        mpatches.Patch(facecolor=C_M_C, edgecolor='white', lw=0.5,
                       label=f'Matched, DIA-NN contaminant ({c_det_c})'),
        mpatches.Patch(facecolor=C_NM_NC, edgecolor='white', lw=0.5,
                       label=f'Not matched, non-contaminant ({c_nd_nc})'),
        mpatches.Patch(facecolor=C_NM_C, edgecolor='white', lw=0.5,
                       label=f'Not matched, DIA-NN contaminant ({c_nd_c})'),
    ]
    leg = ax_e.legend(handles=leg_elements, loc='lower left',
                      bbox_to_anchor=(0.0, 0.02), ncol=1, fontsize=7.5,
                      frameon=True, fancybox=False, edgecolor=C_LGREY,
                      handlelength=1.0, handletextpad=0.4,
                      borderpad=0.5, labelspacing=0.3)
    style_legend(leg)

    ax_e.set_xlim(0, c_total + 2)
    ax_e.set_ylim(-0.85, 1.0)
    ax_e.set_xlabel('Number of protein groups', fontsize=8, color=C_TEXT)
    ax_e.set_yticks([])
    ax_e.spines['left'].set_visible(False)
    style_axes(ax_e, grid_y=False)

    save_multi(fig_e, 'PanelE_CombinedValidation')
    plt.close(fig_e)

    # ── Excel sheet ──
    excel_out = os.path.join(OUTPUT_DIR, 'PanelE_CombinedValidation.xlsx')
    with pd.ExcelWriter(excel_out, engine='openpyxl') as writer:
        # Per-protein match details
        out_cols = ['Protein.Group', 'Genes', 'Is_Contaminant', 'N_Samples_Detected',
                    'PG.MaxLFQ', 'Match_Status', 'Combined_Match']
        available_cols = [c for c in out_cols if c in pg_data.columns]
        pg_data[available_cols].sort_values('Combined_Match').to_excel(
            writer, sheet_name='PerProtein_Match', index=False)

        # Summary
        summary = pd.DataFrame([{
            'Reference': 'SPIDR + Geiger combined',
            'Total_BlotSeq_PG': c_total,
            'Matched': c_det,
            'Matched_NonCont': c_det_nc,
            'Matched_Cont': c_det_c,
            'NotMatched': c_nd,
            'NotMatched_NonCont': c_nd_nc,
            'NotMatched_Cont': c_nd_c,
            'Match_Rate_Total': f'{c_det/c_total*100:.1f}%',
            'Match_Rate_NonCont': f'{c_det_nc/(c_det_nc+c_nd_nc)*100:.1f}%' if (c_det_nc+c_nd_nc)>0 else 'N/A',
        }])
        summary.to_excel(writer, sheet_name='Summary', index=False)
    print(f'  Saved: {excel_out}')

    print(f'\n  Result: Combining SPIDR + Geiger increases match rate from '
          f'{n_detected}/{n_total} ({n_detected/n_total*100:.0f}%) to '
          f'{c_det}/{c_total} ({c_det/c_total*100:.0f}%)')
    print(f'  Geiger uniquely adds {c_det - n_detected} proteins not found in SPIDR')
