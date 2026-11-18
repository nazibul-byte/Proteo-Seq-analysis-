"""
BlotSeq organoid MS: blank exclusion via Limit of Blank (LoB), with
per-gene abundance and q-value carried through to the output.

METHOD
------
Rather than excluding a gene if it is detected in a blank at ALL (too
aggressive for common background contaminants such as keratins, which
can be genuinely expressed at high levels in real samples while also
appearing at low trace levels in blanks), this script computes a
per-gene, per-volume-condition Limit of Blank (LoB):

    LoB(gene) = mean(blank replicates) + Z * SD(blank replicates)

A gene is only retained in a given replicate if its abundance in that
replicate EXCEEDS the LoB computed for that gene from the matching
blank condition (HV samples checked against Blank-HV; LV samples
checked against Blank-LV).

REFERENCES
----------
1. Armbruster, D.A. & Pry, T. "Limit of Blank, Limit of Detection and
   Limit of Quantitation." Clin Biochem Rev. 2008;29(Suppl 1):S49-52.
   PMID: 18852857.
   -> Source of the LoB = mean(blank) + 1.645*SD(blank) formalism
      (one-sided 95% confidence bound on blank-level signal).

2. Mule, M.P., Martins, A.J. & Tsang, J.S. "Normalizing and denoising
   protein expression data from droplet-based single cell profiling."
   Nat Commun 13, 2099 (2022). PMID: 35440536.
   -> Single-cell protein quantification precedent: uses "empty
      droplets" (no-cell controls, directly analogous to blank wells)
      to establish a background distribution, and calls a protein
      genuinely detected only when its signal is a defined number of
      standard deviations above that empty/blank background level.

Both are documented on the ReadMe sheet of the output workbook.

INPUT: report.tsv (raw DIA-NN long-format output for experiment AMHE027)
OUTPUT: BlotSeq_Organoid_LoB_Gene_Abundance_QValue.xlsx
"""

import re
import os
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

# ============================================================
# CONFIGURATION - edit these for your machine
# ============================================================
INPUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = 'report.tsv'
OUTPUT_FILE = 'BlotSeq_Organoid_LoB_Gene_Abundance_QValue.xlsx'

PG_QVALUE_THRESHOLD = 0.01     # Global.PG.Q.Value threshold for a raw detection call
LOB_Z_SCORE = 1.645            # 1.645 = one-sided 95% CI (Armbruster & Pry, 2008)
# ============================================================

CONDITIONS = ['BasalHV', 'BasalLV', 'LuminalLV', 'LuminalHV', 'BlankHV', 'BlankLV']
DISPLAY_NAME = {'BasalHV': 'Basal-HV', 'BasalLV': 'Basal-LV', 'LuminalLV': 'Luminal-LV',
                'LuminalHV': 'Luminal-HV', 'BlankHV': 'Blank-HV', 'BlankLV': 'Blank-LV'}
REP_COUNT = {'BasalHV': 5, 'BasalLV': 5, 'LuminalLV': 5, 'LuminalHV': 5, 'BlankHV': 3, 'BlankLV': 3}
VOLUME_OF = {'BasalHV': 'HV', 'LuminalHV': 'HV', 'BasalLV': 'LV', 'LuminalLV': 'LV'}


def parse_condition_replicate(file_name):
    """Extract Condition and Replicate from a DIA-NN File.Name string.
    e.g. '...BasalLV_2hrdigestionhig_F1_R3_T1_DIA...' -> ('BasalLV', 3)"""
    base = file_name.split('/')[-1]
    m = re.match(r'([A-Za-z]+)_2hrdigestionhig_F1_R(\d+)_T1_DIA', base)
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def load_data(filepath):
    print(f'Loading {filepath} ...')
    report = pd.read_csv(filepath, sep='\t')
    parsed = report['File.Name'].apply(lambda x: pd.Series(parse_condition_replicate(x), index=['Condition', 'Replicate']))
    report = pd.concat([report, parsed], axis=1)
    n_unparsed = report['Condition'].isna().sum()
    if n_unparsed:
        print(f'  WARNING: {n_unparsed} rows did not match the expected filename pattern and were dropped.')
    report = report.dropna(subset=['Condition'])

    # Raw detection filter: passes Global.PG.Q.Value threshold AND has a gene annotation.
    # (Rows with no gene annotation are excluded entirely, not kept as blank gene names -
    #  DIA-NN could not assign roughly 400 protein groups to a gene symbol in this dataset.)
    passing = report[(report['Global.PG.Q.Value'] <= PG_QVALUE_THRESHOLD) & (report['Genes'].notna())].copy()
    print(f'  {len(report)} total rows, {len(passing)} pass Global.PG.Q.Value <= {PG_QVALUE_THRESHOLD} with a gene annotation.')
    return passing


def build_matrices(passing):
    """Build Genes x (Condition,Replicate) matrices for abundance (Genes.MaxLFQ)
    and q-value (GG.Q.Value). Both are already confirmed 1:1 per (run, gene)."""
    dedup = passing.drop_duplicates(subset=['Condition', 'Replicate', 'Genes'])
    abundance = dedup.pivot_table(index='Genes', columns=['Condition', 'Replicate'],
                                   values='Genes.MaxLFQ', aggfunc='first')
    qvalue = dedup.pivot_table(index='Genes', columns=['Condition', 'Replicate'],
                                values='GG.Q.Value', aggfunc='first')
    return abundance, qvalue


def compute_lob(abundance, volume):
    """LoB(gene) = mean(blank replicates for this volume) + Z * SD(blank replicates).
    Missing values (gene not detected in that blank replicate) are treated as 0,
    consistent with treating non-detection as zero signal for that run."""
    blank_cond = 'BlankHV' if volume == 'HV' else 'BlankLV'
    n_reps = REP_COUNT[blank_cond]
    cols = [(blank_cond, r) for r in range(1, n_reps + 1) if (blank_cond, r) in abundance.columns]
    blank_vals = abundance[cols].fillna(0)
    lob = blank_vals.mean(axis=1) + LOB_Z_SCORE * blank_vals.std(axis=1, ddof=1)
    return lob.fillna(0)  # genes never seen anywhere have undefined SD; treat as LoB=0


def build_gene_lists(abundance, qvalue, lob_hv, lob_lv):
    """For each condition/replicate, return a sorted list of (gene, abundance, qvalue)
    for genes whose abundance exceeds the volume-matched LoB."""
    results = {}
    for cond in ['BasalHV', 'BasalLV', 'LuminalLV', 'LuminalHV']:
        lob = lob_hv if VOLUME_OF[cond] == 'HV' else lob_lv
        for rep in range(1, REP_COUNT[cond] + 1):
            col = (cond, rep)
            if col not in abundance.columns:
                results[(cond, rep)] = []
                continue
            ab_col = abundance[col].dropna()
            above = ab_col[ab_col > lob.reindex(ab_col.index).fillna(0)]
            rows = sorted(
                [(gene, above[gene], qvalue.loc[gene, col] if col in qvalue.columns and gene in qvalue.index else None)
                 for gene in above.index]
            )
            results[(cond, rep)] = rows
    return results


def build_raw_gene_lists(abundance, qvalue):
    """Raw (non-blank-excluded) gene/abundance/qvalue lists, for the Blank-HV/LV
    reference sheets and for comparison."""
    results = {}
    for cond in CONDITIONS:
        for rep in range(1, REP_COUNT[cond] + 1):
            col = (cond, rep)
            if col not in abundance.columns:
                results[(cond, rep)] = []
                continue
            ab_col = abundance[col].dropna()
            rows = sorted(
                [(gene, ab_col[gene], qvalue.loc[gene, col] if col in qvalue.columns and gene in qvalue.index else None)
                 for gene in ab_col.index]
            )
            results[(cond, rep)] = rows
    return results


# ---------------------- Workbook writing ----------------------

ARIAL = 'Arial'
HEADER_FONT = Font(name=ARIAL, bold=True, color='FFFFFF', size=10)
HEADER_FILL = PatternFill(start_color='305496', end_color='305496', fill_type='solid')
LOB_FILL = PatternFill(start_color='548235', end_color='548235', fill_type='solid')  # green = LoB-filtered
SUBHEADER_FONT = Font(name=ARIAL, bold=True, size=9, color='FFFFFF')
SUBHEADER_FILL = PatternFill(start_color='4472C4', end_color='4472C4', fill_type='solid')
BODY_FONT = Font(name=ARIAL, size=10)
NOTE_FONT = Font(name=ARIAL, size=9, italic=True, color='555555')
TITLE_FONT = Font(name=ARIAL, bold=True, size=14)
THIN = Side(style='thin', color='B7B7B7')
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)


def write_triplet_sheet(ws, groups, header_fill, note_text):
    col = 1
    max_len = 0
    for label, rows in groups:
        ws.cell(row=1, column=col, value=label)
        ws.cell(row=1, column=col).font = HEADER_FONT
        ws.cell(row=1, column=col).fill = header_fill
        ws.merge_cells(start_row=1, start_column=col, end_row=1, end_column=col + 2)
        ws.cell(row=1, column=col).alignment = Alignment(horizontal='center')
        for j, sub in enumerate(['Gene', 'Abundance', 'Q-value']):
            c = ws.cell(row=2, column=col + j, value=sub)
            c.font = SUBHEADER_FONT
            c.fill = SUBHEADER_FILL
            c.alignment = Alignment(horizontal='center')
            c.border = BORDER
        for i, (gene, ab, qv) in enumerate(rows, start=3):
            ws.cell(row=i, column=col, value=gene).font = BODY_FONT
            ws.cell(row=i, column=col + 1, value=round(float(ab), 1) if ab is not None else None).font = BODY_FONT
            ws.cell(row=i, column=col + 2, value=round(float(qv), 6) if qv is not None else None).font = BODY_FONT
        max_len = max(max_len, len(rows))
        ws.column_dimensions[get_column_letter(col)].width = 14
        ws.column_dimensions[get_column_letter(col + 1)].width = 13
        ws.column_dimensions[get_column_letter(col + 2)].width = 11
        col += 3
    ws.cell(row=max_len + 4, column=1, value=note_text).font = NOTE_FONT


def write_readme(ws, lob_z):
    lines = [
        ("BlotSeq organoid MS - blank exclusion via Limit of Blank (LoB)", TITLE_FONT),
        ("", BODY_FONT),
        ("Method:", SUBHEADER_FONT),
        (f"For each gene, LoB = mean(blank replicates) + {lob_z} x SD(blank replicates), computed", BODY_FONT),
        ("separately for HV blanks and LV blanks (volume-matched). A gene is retained in a given", BODY_FONT),
        ("replicate only if its abundance (Genes.MaxLFQ) exceeds the LoB computed from the", BODY_FONT),
        ("matching-volume blank condition. This replaces binary 'exclude if detected in blank at all'", BODY_FONT),
        ("filtering, which incorrectly discards real signal for proteins (e.g. keratins) that are", BODY_FONT),
        ("common low-level background contaminants but genuinely, substantially expressed in real samples.", BODY_FONT),
        ("", BODY_FONT),
        ("References:", SUBHEADER_FONT),
        ("1. Armbruster, D.A. & Pry, T. Limit of Blank, Limit of Detection and Limit of Quantitation.", BODY_FONT),
        ("   Clin Biochem Rev. 2008;29(Suppl 1):S49-52. PMID: 18852857.", BODY_FONT),
        ("   (Source of the LoB = mean(blank) + 1.645*SD(blank) formula.)", BODY_FONT),
        ("2. Mule, M.P., Martins, A.J. & Tsang, J.S. Normalizing and denoising protein expression", BODY_FONT),
        ("   data from droplet-based single cell profiling. Nat Commun 13, 2099 (2022). PMID: 35440536.", BODY_FONT),
        ("   (Single-cell protein quantification precedent: empty droplets / no-cell controls used as", BODY_FONT),
        ("   the background distribution, with detection called at a defined number of SD above it.)", BODY_FONT),
        ("", BODY_FONT),
        ("Sheets:", SUBHEADER_FONT),
        ("- Basal-HV, Basal-LV, Luminal-LV, Luminal-HV: LoB-filtered gene / abundance / q-value,", BODY_FONT),
        ("  5 replicates each (green header).", BODY_FONT),
        ("- Blank-HV, Blank-LV: raw (unfiltered) blank gene / abundance / q-value, shown for", BODY_FONT),
        ("  transparency, since these are what the LoB values were computed from.", BODY_FONT),
        ("- Abundance = Genes.MaxLFQ (DIA-NN native gene-level quantity).", BODY_FONT),
        ("- Q-value = GG.Q.Value (DIA-NN native gene-group q-value); may occasionally exceed 0.01", BODY_FONT),
        ("  since raw detection used Global.PG.Q.Value (protein-group level) as the initial filter.", BODY_FONT),
    ]
    for i, (text, font) in enumerate(lines, start=1):
        c = ws.cell(row=i, column=1, value=text)
        c.font = font
    ws.column_dimensions['A'].width = 100


def main():
    filepath = f'{INPUT_DIR}/{INPUT_FILE}'
    passing = load_data(filepath)
    abundance, qvalue = build_matrices(passing)

    lob_hv = compute_lob(abundance, 'HV')
    lob_lv = compute_lob(abundance, 'LV')

    lob_filtered = build_gene_lists(abundance, qvalue, lob_hv, lob_lv)
    raw_lists = build_raw_gene_lists(abundance, qvalue)

    for cond in ['BasalHV', 'BasalLV', 'LuminalLV', 'LuminalHV']:
        counts = [len(lob_filtered[(cond, r)]) for r in range(1, REP_COUNT[cond] + 1)]
        print(f'{DISPLAY_NAME[cond]} (LoB-filtered): {counts}')

    wb = Workbook()
    wb.remove(wb.active)

    ws_readme = wb.create_sheet('ReadMe')
    write_readme(ws_readme, LOB_Z_SCORE)

    note = ("Gene retained only if Genes.MaxLFQ exceeds this gene's Limit of Blank "
            "(mean(blank) + %.3f x SD(blank), volume-matched). See ReadMe sheet for method and references." % LOB_Z_SCORE)
    for cond in ['BasalHV', 'BasalLV', 'LuminalLV', 'LuminalHV']:
        ws = wb.create_sheet(DISPLAY_NAME[cond])
        groups = [(f'Rep {r}', lob_filtered[(cond, r)]) for r in range(1, REP_COUNT[cond] + 1)]
        write_triplet_sheet(ws, groups, LOB_FILL, note)

    blank_note = "Raw (unfiltered) gene / abundance / q-value for this blank condition. Used to compute the LoB thresholds applied to the condition sheets."
    for cond in ['BlankHV', 'BlankLV']:
        ws = wb.create_sheet(DISPLAY_NAME[cond])
        groups = [(f'Rep {r}', raw_lists[(cond, r)]) for r in range(1, REP_COUNT[cond] + 1)]
        write_triplet_sheet(ws, groups, HEADER_FILL, blank_note)

    outpath = f'{OUTPUT_DIR}/{OUTPUT_FILE}'
    wb.save(outpath)
    print(f'\nSaved: {outpath}')


if __name__ == '__main__':
    main()
