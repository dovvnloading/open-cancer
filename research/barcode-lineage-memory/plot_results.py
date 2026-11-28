#!/usr/bin/env python3
"""Original scientific figure from freshly recomputed aggregate results."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", type=Path, default=ROOT/".compute/barcode-lineage-memory")
    parser.add_argument("--out", type=Path, default=ROOT/"results/barcode-lineage-memory/run-20261128")
    args = parser.parse_args()
    hct = json.loads((args.reports/"HCTbc.json").read_text())
    exp = json.loads((args.reports/"HCTbc_exploratory.json").read_text())
    original = [f["gain"] for f in hct["primary"]["folds"]]
    trimmed = [next(v["lower_gain"] for v in f["fixed_predictor_target_only_bounds"] if v["epsilon"] == .01)
               for f in hct["primary"]["folds"]]
    donor = next(a for a in exp["ablations"] if a["mode"] == "donor_feature_only" and a["floor"] == 1e-5)
    late = next(a for a in exp["ablations"] if a["mode"] == "late_only" and a["floor"] == 1e-5)
    masses = [100*(1-f["donor_feature_mass_retained"]) for f in donor["folds"]]
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10.5,
                         "axes.spines.top":False, "axes.spines.right":False,
                         "text.color":"#182638", "axes.labelcolor":"#263549", "svg.fonttype":"none"})
    fig, axes = plt.subplots(1,2,figsize=(12.8,5.3),gridspec_kw={"width_ratios":[1,1.1]})
    fig.suptitle("Rare early predictors carry the gain in late-lineage forecasts", x=.06, y=.98,
                 ha="left", fontsize=17, fontweight="bold")
    fig.text(.06,.905,"HCT116 · four held-out treated cultures · recorded primary verdict: FAIL · reconstructed and rerun", fontsize=11,color="#546475")
    ax = axes[0]
    x = np.arange(4)
    ax.bar(x-.16,original,width=.30,color="#2865A0",label="Original forecast gain")
    ax.bar(x+.16,trimmed,width=.30,color="#45A199",label="Worst 1% late-read deletion")
    ax.set_xticks(x,["DT1","DT2","DT3","DT4"])
    ax.set_ylim(0,1.04)
    ax.set_ylabel("Incremental log score (nats per read)")
    ax.set_title("A  Late-read checks retain positive gain",loc="left",fontsize=12,pad=16,fontweight="bold")
    ax.legend(loc="upper left",frameon=False,fontsize=9)
    ax.text(.02,-.22,"Predictions fixed; this check challenges only late outcomes.",transform=ax.transAxes,
            fontsize=9.3,color="#546475",va="top")
    ax=axes[1]
    values=[hct["primary"]["mean_gain"],late["mean_gain"],donor["mean_gain"]]
    ax.bar(range(3),values,color=["#2865A0","#45A199","#C4514A"],width=.58)
    ax.set_xticks(range(3),["Original","Rare late reads\nremoved","Rare donor predictors\nremoved"])
    ax.set_ylim(0,.41)
    ax.set_ylabel("Mean incremental log score (nats per read)")
    ax.set_title("B  An early-predictor check eliminates it",loc="left",fontsize=12,pad=16,fontweight="bold")
    for i,value in enumerate(values):
        ax.text(i,value+.013,f"{value:.3f}",ha="center",fontsize=12,fontweight="bold",color="#C4514A" if i==2 else "#182638")
    ax.text(.02,-.22,f"Exploratory refit: {min(masses):.2f}–{max(masses):.2f}% of donor read mass removed.\nOwn early data, late outcomes, baseline and dictionary unchanged.",
            transform=ax.transAxes,fontsize=9.3,color="#546475",va="top")
    for ax in axes:
        ax.yaxis.grid(True,color="#E8EDF2",zorder=0)
        ax.set_axisbelow(True)
    fig.text(.06,.025,"Panel B removes entries below 0.001% frequency per barcode. Fragility is established; contamination causation is not.",fontsize=9.5,color="#546475")
    fig.subplots_adjust(left=.07,right=.98,top=.79,bottom=.27,wspace=.28)
    args.out.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.out/"forecast_robustness.png",dpi=180,facecolor="white")
    fig.savefig(args.out/"forecast_robustness.svg",facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
