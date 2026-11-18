"""Reproducible exploratory Proteo-Seq variance analysis.

Run from this folder: python variance_analysis.py
Inputs are read only. No missing protein abundances are imputed.
Outputs: CSV audit tables, publication-quality charts, run manifest, worked examples.
The revised baseline uses diameter and elongation and v1 baseline inclusion rules. A separately labelled matched-gene
analysis adds GeneScore and recomputes ALL allocations on the SAME eligible cells.
These are descriptive associations, not causal or biological/technical noise fractions.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import math
import platform
import sys
import io

import numpy as np
import pandas as pd
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLOCKS = ["Morphology", "Processing condition G", "Identification count"]
COLORS = {"Morphology": "#4C91AF", "Processing condition G": "#D59A58",
          "Identification count": "#A69B80", "GeneScore": "#9273AC",
          "Unexplained": "#D8DDE1"}


def read_block(wb, sheet, ncell, q=False):
    rows = list(wb[sheet].values)
    cells = [str(x) for x in rows[1][1:1+ncell]]
    body = [r for r in rows[2:] if r[0] is not None]
    genes = [str(r[0]) for r in body]
    if len(genes) != len(set(genes)):
        raise ValueError(f"Duplicate exact row keys in {sheet}")
    a = pd.DataFrame([r[1:1+ncell] for r in body], index=genes,
                     columns=cells).apply(pd.to_numeric, errors="coerce")
    a = a.where(a > 0) if sheet.startswith("Protein") else a
    if q:
        qv = pd.DataFrame([r[15:25] for r in body], index=genes,
                          columns=cells).apply(pd.to_numeric, errors="coerce")
        a = a.where(qv <= 0.05)
    return a


def load_data(path, system):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    if system == "MCF7":
        groups, ncell = ["1hr", "2hr", "3hr"], 5
    else:
        groups, ncell = ["Basal", "Luminal"], 10
    p = pd.concat([read_block(wb, "Protein_"+g, ncell, system=="PDO")
                   for g in groups], axis=1)
    a = pd.concat([read_block(wb, "ATAC_"+g, ncell) for g in groups], axis=1)
    rows = list(wb["Cell_Morphology"].values)
    m = pd.DataFrame([r for r in rows[2:] if isinstance(r[0], str)
                      and r[0].startswith("R") and "C" in r[0]], columns=rows[1])
    m = m.rename(columns={"Cell (R#C#)": "cell", "Volume (um3)": "volume", "Diameter (um)": "diameter",
         "Elongation": "elongation", "Digestion condition": "processing",
         "Group": "recorded_group", "Protein groups (this workbook)": "count",
         "Protein groups identified": "count", "TSS enrichment": "TSS"})
    m = m.set_index("cell")
    for col in ["diameter", "volume", "elongation", "count"]:
        m[col] = pd.to_numeric(m[col], errors="raise")
        if not np.all(np.isfinite(m[col]) & (m[col] > 0)):
            raise ValueError(f"Invalid {col} in {system}")
    if system == "PDO":
        m["TSS"] = pd.to_numeric(m["TSS"], errors="raise")
    if not (set(p.columns) == set(a.columns) == set(m.index)):
        raise ValueError(f"Cell identifiers do not align in {system}")
    wb.close()
    return p.loc[:, m.index], a.loc[:, m.index], m


def make_blocks(m, system):
    b = {"Morphology": np.column_stack([np.log2(m.diameter), np.log2(m.elongation)])}
    if system == "MCF7":
        b["Processing condition G"] = np.column_stack([(m.processing == "2hr").astype(float),
                                                 (m.processing == "3hr").astype(float)])
    else:
        b["Processing condition G"] = (m.processing == "Low volume").astype(float).to_numpy()[:,None]
    b["Identification count"] = np.log2(m["count"].to_numpy())[:, None]
    return b


def design(blocks, subset, n):
    return np.column_stack([np.ones(n)] + [blocks[k] for k in subset])


def fit(y, x):
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    predicted = x @ beta
    residual = y - predicted
    sse = float(residual @ residual)
    sst = float(np.sum((y-y.mean())**2))
    rank = int(np.linalg.matrix_rank(x))
    r2 = 1-sse/sst
    adjusted = 1-(1-r2)*(len(y)-1)/(len(y)-rank) if len(y)>rank else np.nan
    return dict(beta=beta, predicted=predicted, residual=residual, SSE=sse, SST=sst,
                R2=r2, adjusted_R2=adjusted, rank=rank)


def allocation(y, blocks):
    """Exact Shapley/LMG allocation over predictor BLOCKS, not individual columns."""
    names = list(blocks)
    fits = {}
    for n in range(len(names)+1):
        for subset in itertools.combinations(names, n):
            f = fit(y, design(blocks, subset, len(y)))
            fits[frozenset(subset)] = f
    gains = {name: [] for name in names}
    order_rows = []
    for order in itertools.permutations(names):
        prev = frozenset()
        record = {"order": " -> ".join(order)}
        for name in order:
            nxt = prev | {name}
            gain = fits[nxt]["R2"]-fits[prev]["R2"]
            if gain < -1e-8:
                raise AssertionError("A nested least-squares fit lost R2")
            gains[name].append(gain)
            record[name] = gain
            prev = nxt
        order_rows.append(record)
    phi = {name: float(np.mean(v)) for name, v in gains.items()}
    full = fits[frozenset(names)]
    phi["Unexplained"] = 1-full["R2"]
    if not np.isclose(sum(phi.values()), 1, atol=1e-8):
        raise AssertionError("Allocations do not sum to one")
    return phi, fits, order_rows


def loo(y, x):
    """OLS leave-one-cell-out predictions; reject folds with lost model rank."""
    rank = np.linalg.matrix_rank(x)
    pred = np.empty(len(y))
    for i in range(len(y)):
        keep = np.arange(len(y)) != i
        if np.linalg.matrix_rank(x[keep]) < rank:
            return None
        pred[i] = x[i] @ np.linalg.lstsq(x[keep], y[keep], rcond=None)[0]
    press = np.sum((y-pred)**2)
    sst = np.sum((y-y.mean())**2)
    return float(1-press/sst), pred


def eligible_gene_names(p, a):
    upper = pd.Index([g.upper() for g in p.index])
    ambiguous = set(upper[upper.duplicated(keep=False)])
    return [g for g in p.index if ";" not in g and g==g.upper()
            and g.upper() not in ambiguous and g in a.index]


def save_figure(fig, stem):
    for extension in ["png"]:
        buffer = io.BytesIO()
        fig.savefig(buffer, format=extension, dpi=1200, bbox_inches="tight", facecolor="white")
        Path(str(stem)+"."+extension).write_bytes(buffer.getvalue())
    plt.close(fig)


def bar_chart(summary, components, title, subtitle, stem):
    fig, ax = plt.subplots(figsize=(7.6, 5.5))
    systems = ["MCF7", "PDO"]
    bottom = np.zeros(2)
    for comp in components:
        heights = np.array([summary.loc[s, comp] for s in systems])*100
        bars = ax.bar([0,1], heights, bottom=bottom, width=.55,
                      color=COLORS[comp], edgecolor="white", linewidth=1.1, label=comp)
        for b, h, low in zip(bars, heights, bottom):
            if h >= 4:
                ax.text(b.get_x()+b.get_width()/2, low+h/2, f"{h:.1f}%",
                        ha="center", va="center", fontsize=11)
        bottom += heights
    labels = [f"{s}\n{int(summary.loc[s,'n_cells'])} cells; "
              f"{int(summary.loc[s,'n_proteins'])} protein rows" for s in systems]
    ax.set_xticks([0,1], labels)
    ax.set_ylim(0,100)
    ax.set_ylabel("Mean allocated log-protein variance (%)")
    ax.set_title(title, loc="left", fontsize=14, weight="bold", pad=40)
    ax.text(0, 1.035, subtitle, transform=ax.transAxes, fontsize=9, va="bottom")
    ax.spines[["top","right"]].set_visible(False)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01,1), frameon=False, fontsize=9)
    fig.subplots_adjust(left=.11,right=.70,bottom=.20,top=.79)
    fig.text(.11,.025,"Descriptive in-sample allocation; equal weight per protein.\n"
             "Fractions are not causal or calibrated biological/technical noise components.",fontsize=9)
    save_figure(fig, stem)


def main():
    home = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mcf7", type=Path, default=home/"inputs"/"MCF7_SameCell_Multiomics_Dataset(1).xlsx")
    ap.add_argument("--pdo", type=Path, default=home.parent/"07_PDO_proteomics"/"BlotSeq_PDO_SameCell_Multiomics_Dataset.xlsx")
    ap.add_argument("--out", type=Path, default=home/"results")
    ap.add_argument("--mcf7-min-cells", type=int, default=9)
    ap.add_argument("--pdo-min-cells", type=int, default=11)
    ap.add_argument("--tss", type=float, default=.9, help="Historical PDO exploratory cutoff; not high-quality ATAC certification")
    ap.add_argument("--permutations", type=int, default=399)
    ap.add_argument("--seed", type=int, default=20261025)
    args = ap.parse_args()
    if not 7<=args.mcf7_min_cells<=15 or not 7<=args.pdo_min_cells<=20 or args.permutations<1:
        ap.error("minimum cells must be 7-15 (MCF7) or 7-20 (PDO); permutations must be positive")
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    data = {s:load_data(path,s) for s,path in [("MCF7",args.mcf7),("PDO",args.pdo)]}
    records, subrecords, orderrecords, summaries, diagnostics = [], [], [], [], []
    matched_records, matched_subsets, matched_orders, cv_records, skipped = [],[],[],[],[]
    matched_summaries, nullrecords, examples = [], [], {}
    for system,(p,a,m) in data.items():
        print(f"Analyzing {system}: {len(m)} indexed cells", flush=True)
        p.to_csv(out/f"{system}_filtered_protein_abundance.csv", index_label="protein")
        a.to_csv(out/f"{system}_GeneScores.csv", index_label="gene")
        m.to_csv(out/f"{system}_cell_metadata.csv")
        b_all = make_blocks(m,system)
        eligible = p.index[p.notna().sum(axis=1)>=(args.mcf7_min_cells if system=="MCF7" else args.pdo_min_cells)]
        local_records, null_cache = [], []
        for g in eligible:
            obs = p.loc[g].notna().to_numpy()
            y = np.log2(p.loc[g].to_numpy(dtype=float)[obs])
            if np.var(y)==0:
                skipped.append(dict(system=system,protein=g,analysis="baseline",reason="constant protein"));continue
            b = {k:v[obs] for k,v in b_all.items()}
            phi,fits,orders = allocation(y,b)
            full = fits[frozenset(b)]
            row = dict(system=system,protein=g,n_obs=len(y),rank=full["rank"],
                       SST=full["SST"],SSE=full["SSE"],R2=full["R2"],
                       adjusted_R2=full["adjusted_R2"],**phi)
            basecv=loo(y,design(b,list(b),len(y)))
            row["LOO_R2"]=basecv[0] if basecv else np.nan
            local_records.append(row);records.append(row)
            null_cache.append((obs,y))
            for ss,f in fits.items():
                subrecords.append(dict(system=system,protein=g,subset=" + ".join(k for k in b if k in ss) or "Intercept",R2=f["R2"],SSE=f["SSE"],SST=f["SST"],rank=f["rank"]))
            orderrecords += [dict(system=system,protein=g,**r) for r in orders]
            if g=="DSP":
                cells=m.index[obs]
                table=m.loc[cells,["diameter","elongation","count","processing"]].copy()
                table["abundance"]=p.loc[g,cells]
                table["log2_abundance"]=y
                table["predicted"]=full["predicted"]
                table["residual"]=full["residual"]
                table["residual_squared"]=full["residual"]**2
                table["centered_y_squared"]=(y-y.mean())**2
                table.to_csv(out/f"worked_example_{system}_DSP_cells.csv")
                examples[system]=dict(protein=g,n_obs=len(y),mean_log2=float(y.mean()),SST=full["SST"],SSE=full["SSE"],R2=full["R2"],allocations=phi,orders=orders)
        tab=pd.DataFrame(local_records)
        summary=dict(system=system,n_cells=len(m),n_proteins=len(tab))
        summary.update({k:float(tab[k].mean()) for k in BLOCKS+["Unexplained","R2","adjusted_R2"]})
        summaries.append(summary)
        # Joint permutation of predictor rows: preserves relationships among predictors.
        # Each protein retains its original response and observed-cell mask.
        null=[]
        for z in range(args.permutations):
            perm=rng.permutation(len(m));bp={k:v[perm] for k,v in b_all.items()}
            val=np.mean([fit(y,design({k:v[obs] for k,v in bp.items()},BLOCKS,len(y)))["R2"] for obs,y in null_cache])
            null.append(val);nullrecords.append(dict(system=system,draw=z+1,mean_R2=val))
        diagnostics.append(dict(system=system,mean_R2=summary["R2"],mean_adjusted_R2=summary["adjusted_R2"],
            permuted_mean_R2=float(np.mean(null)),permuted_p95_R2=float(np.quantile(null,.95)),
            median_baseline_LOO_R2=float(tab.LOO_R2.median()),n_baseline_LOO=int(tab.LOO_R2.notna().sum())))
        # GeneScore extension: same observed cells in baseline and extended fit.
        cells=m.index[m.TSS>=args.tss] if system=="PDO" else m.index
        mm=m.loc[cells];bb=make_blocks(mm,system)
        local_matched=[]
        for g in eligible_gene_names(p,a):
            pv=p.loc[g,cells].to_numpy(dtype=float); av=a.loc[g,cells].to_numpy(dtype=float)
            ok=np.isfinite(pv)&(pv>0)&np.isfinite(av)&(av>=0)
            reason=None
            if ok.sum()<10:reason="fewer than 10 paired observations"
            if reason:
                skipped.append(dict(system=system,protein=g,analysis="matched",reason=reason));continue
            y=np.log2(pv[ok]);z=np.log2(1+av[ok])
            b={k:v[ok] for k,v in bb.items()}
            b["GeneScore"]=z[:,None]
            x0=design(b,BLOCKS,len(y));x1=design(b,list(b),len(y))
            if np.std(z)==0 or np.std(y)==0:reason="constant GeneScore or protein"
            elif len(y)-np.linalg.matrix_rank(x1)<2:reason="fewer than 2 residual degrees of freedom"
            if reason:
                skipped.append(dict(system=system,protein=g,analysis="matched",reason=reason));continue
            cv0,cv1=loo(y,x0),loo(y,x1)
            if cv0 is None or cv1 is None:
                skipped.append(dict(system=system,protein=g,analysis="matched",reason="leave-one-out training rank loss"));continue
            f0,f1=fit(y,x0),fit(y,x1)
            phi,fits,orders=allocation(y,b)
            row=dict(system=system,protein=g,n_obs=len(y),rank=f1["rank"],
                baseline_R2=f0["R2"],extended_R2=f1["R2"],
                incremental_GeneScore_R2=f1["R2"]-f0["R2"],**phi)
            matched_records.append(row);local_matched.append(row)
            matched_subsets += [dict(system=system,protein=g,subset=" + ".join(k for k in b if k in ss) or "Intercept",R2=f["R2"],SSE=f["SSE"],SST=f["SST"],rank=f["rank"]) for ss,f in fits.items()]
            matched_orders += [dict(system=system,protein=g,**r) for r in orders]
            cv_records.append(dict(system=system,protein=g,n_obs=len(y),baseline_LOO_R2=cv0[0],extended_LOO_R2=cv1[0],delta_LOO_R2=cv1[0]-cv0[0]))
        mt=pd.DataFrame(local_matched)
        ms=dict(system=system,n_cells=len(cells),n_proteins=len(mt))
        ms.update({k:float(mt[k].mean()) for k in BLOCKS+["GeneScore","Unexplained","baseline_R2","extended_R2","incremental_GeneScore_R2"]})
        matched_summaries.append(ms)
    tables={"baseline_per_protein":records,"baseline_subset_fits":subrecords,
        "baseline_entry_orders":orderrecords,"baseline_summary":summaries,
        "fit_diagnostics":diagnostics,"permuted_design_R2":nullrecords,
        "matched_per_gene":matched_records,"matched_subset_fits":matched_subsets,
        "matched_entry_orders":matched_orders,"matched_summary":matched_summaries,
        "chromatin_prediction":cv_records,"excluded_models":skipped}
    for name,rows in tables.items():pd.DataFrame(rows).to_csv(out/(name+".csv"),index=False)
    cv=pd.DataFrame(cv_records)
    cv_summary=cv.groupby("system").agg(n_genes=("protein","size"),median_baseline_LOO_R2=("baseline_LOO_R2","median"),median_extended_LOO_R2=("extended_LOO_R2","median"),median_delta_LOO_R2=("delta_LOO_R2","median"))
    cv_summary.to_csv(out/"chromatin_prediction_summary.csv")
    # Same-row baseline sensitivity: keeps panels aligned across systems without selecting 10 proteins.
    bt=pd.DataFrame(records)
    shared=set(bt.loc[bt.system=="MCF7","protein"]) & set(bt.loc[bt.system=="PDO","protein"])
    bt[bt.protein.isin(shared)].groupby("system")[BLOCKS+["Unexplained"]].mean().to_csv(out/"shared_eligible_panel_summary.csv")
    (out/"shared_eligible_panel.txt").write_text("\n".join(sorted(shared)),encoding="utf-8")
    plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"font.size":10,"svg.fonttype":"none","pdf.fonttype":42})
    bar_chart(pd.DataFrame(summaries).set_index("system"),BLOCKS+["Unexplained"],
              "Protein variance: morphology and recorded covariates",
              f"MCF7: >= {args.mcf7_min_cells}/15 cells; PDO: >= {args.pdo_min_cells}/20 cells; observed entries only",
              out/"01_baseline_variance")
    bar_chart(pd.DataFrame(matched_summaries).set_index("system"),BLOCKS+["GeneScore","Unexplained"],
              "Exploratory variance allocation with matched chromatin",
              "Separate matched-gene panel; identical cells per nested fit; 24 entry orders",
              out/"02_matched_chromatin_variance")
    fig,axes=plt.subplots(1,2,figsize=(9.0,4.8))
    dt=pd.DataFrame(diagnostics).set_index("system")
    xx=np.arange(2);width=.23
    for j,(col,label,color) in enumerate([("mean_R2","Training fit","#4C91AF"),("permuted_mean_R2","Shuffled-design mean","#B8BEC4"),("mean_adjusted_R2","Rank-adjusted fit","#D59A58")]):
        vals=dt.loc[["MCF7","PDO"],col].to_numpy()*100
        bars=axes[0].bar(xx+(j-1)*width,vals,width,label=label,color=color)
        axes[0].bar_label(bars,fmt="%.1f",fontsize=8,padding=2)
    axes[0].set_xticks(xx,["MCF7","PDO"]);axes[0].set_ylabel("Mean R² (%)")
    axes[0].set_title("A  How much fit arises from model flexibility?",loc="left",fontsize=10)
    axes[0].set_ylim(0,85)
    axes[0].legend(frameon=False,fontsize=8,loc="upper right")
    colors=["#4C91AF","#9273AC"]
    vals=cv_summary.loc[["MCF7","PDO"],"median_delta_LOO_R2"].to_numpy()
    bars=axes[1].bar(xx,vals,.55,color=colors)
    axes[1].bar_label(bars,fmt="%.3f",padding=3,fontsize=10)
    axes[1].axhline(0,color="black",lw=.8);axes[1].set_xticks(xx,["MCF7","PDO"])
    axes[1].set_ylabel("Median change in held-out R²")
    axes[1].set_title("B  Does GeneScore improve held-out prediction?",loc="left",fontsize=10)
    axes[1].set_ylim(min(vals)*1.35,max(.15,max(vals)*1.2))
    for ax in axes:ax.spines[["top","right"]].set_visible(False)
    fig.subplots_adjust(wspace=.45,top=.86,bottom=.23)
    fig.text(.12,.04,"Negative held-out change means worse median prediction after adding GeneScore.\n"
             "Small cohorts and unstable baseline models limit biological interpretation.",fontsize=9)
    save_figure(fig,out/"03_fit_and_prediction_checks")
    from sensitivity_analysis import run_sensitivity
    run_sensitivity(data["PDO"], out)
    manifest={"python":sys.version,"platform":platform.platform(),"numpy":np.__version__,"pandas":pd.__version__,
        "openpyxl":openpyxl.__version__,"matplotlib":matplotlib.__version__,
        "parameters":{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
        "inputs":{s:{"name":path.name,"sha256":hashlib.sha256(path.read_bytes()).hexdigest()} for s,path in [("MCF7",args.mcf7),("PDO",args.pdo)]},
        "important_notes":["Source gene-group rows retained case-sensitively; not all rows are distinct endogenous human proteins.","PDO q<=0.05 additional filtering; MCF7 lacks equivalent per-entry q columns.","No missing protein imputation; zero GeneScores retained; GeneScore transformation log2(1+C).", "MCF7 minimum remains 9/15; PDO minimum is 11/20. Baseline retains v1 inclusion rules; rank and holdout limitations are diagnostics. K19 is excluded from all PDO predictors. Matched panel retains separate >=10 paired-cell eligibility.","PDO matched analysis uses historical exploratory TSS>=0.9; manuscript QC thresholds require reconciliation.","GeneScore is local accessibility, not genome-wide chromatin or transcription rate.","Diameter and elongation are a combined block; count is a coverage proxy, not measured recovery efficiency.","No causal partition or calibrated intrinsic/extrinsic noise estimate."]}
    (out/"run_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    (out/"worked_examples.json").write_text(json.dumps(examples,indent=2),encoding="utf-8")
    print("\nBASELINE",pd.DataFrame(summaries).to_string(index=False),sep="\n")
    print("\nMATCHED",pd.DataFrame(matched_summaries).to_string(index=False),sep="\n")
    print("\nPREDICTION",cv_summary.to_string(),sep="\n")
    print("\nDIAGNOSTICS",pd.DataFrame(diagnostics).to_string(index=False),sep="\n")
    print(f"\nSaved results to {out.resolve()}")


if __name__ == "__main__":
    main()
