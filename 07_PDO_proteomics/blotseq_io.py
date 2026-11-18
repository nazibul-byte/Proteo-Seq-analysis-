"""
blotseq_io.py
Shared loaders for the BlotSeq PDO trimodal figures.

Reads the four source workbooks directly. Place this file in the same folder as
fig1_evidence_composition.py, fig2_chromatin_protein.py and
fig3_measurement_orthogonality.py, then set INPUT_DIR in each figure script.

Python 3.11, requires: pandas, numpy, scipy, matplotlib, openpyxl
"""
import os
import numpy as np
import pandas as pd
from openpyxl import load_workbook

# ----------------------------------------------------------------- constants
F_PROTEIN = 'BlotSeq_Organoid_LoB_Gene_Abundance_QValue.xlsx'
F_ATAC    = 'PDO_GeneScores_524matched_genes_All_cells.xlsx'
F_SIZE    = 'Cell_size.xlsx'
F_SUMMARY = 'BlotSeq_Organoid_Summary_Comparison.xlsx'

# chip map: microwell coordinate -> experimental group
CHIP = {
    'R2C31': 'Luminal-HV', 'R3C31': 'Luminal-HV', 'R9C31': 'Luminal-HV',
    'R11C1': 'Luminal-HV', 'R11C2': 'Luminal-HV',
    'R2C29': 'Luminal-LV', 'R4C31': 'Luminal-LV', 'R9C30': 'Luminal-LV',
    'R13C1': 'Luminal-LV', 'R13C2': 'Luminal-LV',
    'R16C29': 'Basal-LV', 'R14C30': 'Basal-LV', 'R16C26': 'Basal-LV',
    'R16C23': 'Basal-LV', 'R14C24': 'Basal-LV',
    'R16C28': 'Basal-HV', 'R14C28': 'Basal-HV', 'R16C25': 'Basal-HV',
    'R16C24': 'Basal-HV', 'R14C23': 'Basal-HV',
}
LINEAGE = {k: v.split('-')[0] for k, v in CHIP.items()}
DIGEST  = {k: v.split('-')[1] for k, v in CHIP.items()}

# cells excluded on TSS enrichment < 0.9
TSS_FAIL = ['R16C28', 'R16C25', 'R16C24', 'R4C31']
KEEP = [c for c in CHIP if c not in TSS_FAIL]

# display order: Basal-HV, Basal-LV, Luminal-HV, Luminal-LV
_GRP_ORDER = {'Basal-HV': 0, 'Basal-LV': 1, 'Luminal-HV': 2, 'Luminal-LV': 3}
KEEP = sorted(KEEP, key=lambda c: (_GRP_ORDER[CHIP[c]], c))

# ----------------------------------------------------------------- palette
PAL = {
    'basal':      '#9B7EBD',   # lineage colours, identical across all figures
    'luminal':    '#F2A93B',
    'protein':    '#C9A227',   # modality colours
    'chromatin':  '#8FA9C4',
    'dual':       '#3B3B3B',
    'grey':       '#7A7A7A',
    'div_lo':     '#2C5F8A',   # diverging map endpoints
    'div_hi':     '#B5322E',
}

# ----------------------------------------------------------------- rc params
RC = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8, 'axes.labelsize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 1.4,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.width': 1.4, 'ytick.major.width': 1.4,
    'xtick.major.size': 3.2, 'ytick.major.size': 3.2,
    'savefig.bbox': 'tight',
}


def barcode(cell):
    """R16C23 -> Device11#d11_r16_c23"""
    r, c = cell.upper().replace('R', '').split('C')
    return f'Device11#d11_r{int(r)}_c{int(c)}'


# ----------------------------------------------------------------- loaders
def load_protein(input_dir):
    """LoB-filtered MaxLFQ abundances.

    Sheet layout, one block of three columns per replicate:
      row 1  'Rep n'
      row 2  microwell coordinate (cell sheets only)
      row 3  TSS enrichment       (cell sheets only)
      row 4  Gene | Abundance | Q-value
      row 5+ data

    Returns (abundance, qvalue) DataFrames, genes x cells.
    """
    wb = load_workbook(os.path.join(input_dir, F_PROTEIN), read_only=True, data_only=True)
    rec = []
    for sheet in ['Basal-HV', 'Basal-LV', 'Luminal-LV', 'Luminal-HV']:
        rows = list(wb[sheet].iter_rows(values_only=True))
        n_blocks = sum(1 for c in rows[0] if c and str(c).startswith('Rep'))
        for b in range(n_blocks):
            c0 = 3 * b
            cell_id = rows[1][c0]
            for r in rows[4:]:
                if r[c0] is None:
                    continue
                rec.append((cell_id, r[c0], r[c0 + 1], r[c0 + 2]))
    wb.close()
    df = pd.DataFrame(rec, columns=['cell', 'gene', 'abundance', 'qvalue'])
    df['abundance'] = pd.to_numeric(df['abundance'], errors='coerce')
    df['qvalue'] = pd.to_numeric(df['qvalue'], errors='coerce')
    # DIA-NN reports gene groups such as 'ACTB;ACTG1' where peptides map to more
    # than one symbol. Expand to one row per member symbol so the protein layer can
    # be matched against single-symbol ATAC gene names.
    df['gene'] = df['gene'].astype(str).str.split(';')
    df = df.explode('gene')
    df['gene'] = df['gene'].str.strip()
    ab = df.pivot_table(index='gene', columns='cell', values='abundance', aggfunc='max')
    qv = df.pivot_table(index='gene', columns='cell', values='qvalue', aggfunc='min')
    return ab, qv


def load_genescores(input_dir):
    """ArchR gene activity scores, cells x genes.

    Excel had coerced three gene symbols into currency-formatted numbers; the
    original symbol survives in the number_format string, e.g. '[$MLF]0.00'
    with value 2 -> 'MLF2'. Repaired on load (MLF2, PGK1, VNN1).
    """
    wb = load_workbook(os.path.join(input_dir, F_ATAC), data_only=True)
    ws = wb['Sheet 1 - PDO_GeneScores_524vma']
    header = []
    for cell in next(ws.iter_rows(min_row=2, max_row=2)):
        v = cell.value
        if v is not None and not isinstance(v, str):
            fmt = cell.number_format
            if fmt.startswith('[$') and ']' in fmt:
                v = fmt[2:fmt.index(']')] + str(int(v))
        header.append(v)
    body = [list(r)[:len(header)] for r in ws.iter_rows(min_row=3, values_only=True) if r[0]]
    wb.close()
    gs = pd.DataFrame(body, columns=header).set_index('cellBarcode')
    return gs.apply(pd.to_numeric, errors='coerce')


def load_morphology(input_dir):
    """CellenONE pre-lysis morphometry, one row per microwell."""
    wb = load_workbook(os.path.join(input_dir, F_SIZE), read_only=True, data_only=True)
    rows = list(wb['Sheet1 (2)'].iter_rows(values_only=True))
    wb.close()
    v = pd.DataFrame(rows[1:], columns=list(rows[0]))
    v = v[['Lane ID', 'Diameter', 'Elongation', 'Circularity', 'Vol']].dropna(subset=['Lane ID'])
    v['cell'] = v['Lane ID'].str.upper()
    for c in ['Diameter', 'Elongation', 'Circularity', 'Vol']:
        v[c] = pd.to_numeric(v[c], errors='coerce')
    return v.set_index('cell')


def load_summary(input_dir):
    """Per-cell TSS enrichment and identification counts."""
    s = pd.read_excel(os.path.join(input_dir, F_SUMMARY), sheet_name='Summary')
    s = s[s['Chip (R#C#)'].astype(str).str.match(r'^R\d+C\d+$', na=False)]
    s = s.set_index('Chip (R#C#)')
    for c in ['TSS', 'Peptides detected', 'Protein groups detected',
              'Genes (Binary Exclusion)', 'Genes (LoB Exclusion)']:
        s[c] = pd.to_numeric(s[c], errors='coerce')
    return s


def build_matrices(input_dir, cells=None):
    """Convenience wrapper returning aligned matrices for a cell list.

    Returns dict with:
      cells   list of microwell coordinates, in display order
      ab      protein abundance, cells x genes (NaN = not retained)
      qv      protein q-value, cells x genes
      atac    gene activity score, cells x genes (0 = no signal)
      shared  genes present in both the ATAC panel and the protein gene list
      morph   morphometry, cells x [Diameter, Elongation, Circularity, Vol]
      summ    summary table, cells x [TSS, counts]
    """
    cells = list(KEEP if cells is None else cells)
    ab, qv = load_protein(input_dir)
    gs = load_genescores(input_dir)
    morph = load_morphology(input_dir)
    summ = load_summary(input_dir)
    atac = gs.loc[[barcode(c) for c in cells]].copy()
    atac.index = cells
    ab = ab.reindex(columns=cells).T
    qv = qv.reindex(columns=cells).T
    shared = [g for g in atac.columns if g in ab.columns]
    return dict(cells=cells, ab=ab, qv=qv, atac=atac, shared=shared,
                morph=morph.loc[cells], summ=summ.loc[cells])


def style_axes(ax):
    """Full box, inward ticks on all four sides."""
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(1.4)
    ax.tick_params(which='both', direction='in', top=True, right=True)
