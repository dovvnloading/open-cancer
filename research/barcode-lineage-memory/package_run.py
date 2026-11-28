#!/usr/bin/env python3
"""Package original aggregate calculations explicitly, with source/code checksums."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
REPORTS=("HCTbc.json","SW6bc.json","HCTbc_exploratory.json","SW6bc_exploratory.json","independent_verification.json","predictor_counterexamples.json")
PROTOCOL_COMMIT="d6de5c46b4199483509b1ea88cf35aa19639918b"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,default=ROOT/".compute/barcode-lineage-memory")
    parser.add_argument("--out",type=Path,default=ROOT/"results/barcode-lineage-memory/run-20261128")
    args=parser.parse_args()
    args.input=args.input.expanduser().resolve()
    args.out=args.out.expanduser().resolve()
    args.out.mkdir(parents=True,exist_ok=True)
    if not (args.out/"reproduce.py").is_file():
        raise ValueError("Run-specific reproduction entry point required before packaging")
    for name in REPORTS:
        shutil.copyfile(args.input/name,args.out/name)
    hct=json.loads((args.out/"HCTbc.json").read_text())
    sw=json.loads((args.out/"SW6bc.json").read_text())
    exploration=json.loads((args.out/"HCTbc_exploratory.json").read_text())
    verification=json.loads((args.out/"independent_verification.json").read_text())
    for check in verification["checks"]:
        suffix="_exploratory" if check.get("exploratory_A2") else ""
        if not check["matched"] or digest(args.out/f"{check['cell_line']}{suffix}.json") != check["report_sha256"]:
            raise ValueError("Independent replay must match the exact reports being packaged")
    donor=next(v for v in exploration["ablations"] if v["mode"]=="donor_feature_only" and v["floor"]==1e-5)
    late=next(v for v in exploration["ablations"] if v["mode"]=="late_only" and v["floor"]==1e-5)
    removed=[1-f["donor_feature_mass_retained"] for f in donor["folds"]]
    temporal=exploration["temporal_diagnostic"]
    halves=[f for h in temporal["home_results"] for f in h.get("receiver_half_results",[]) if f["available"]]
    if (hct["decision"]["verdict"]!="FAIL" or len(hct["primary"]["folds"])!=4
            or not all(f["gain"]>0 for f in hct["primary"]["folds"])
            or len(donor["folds"])!=4 or any(f["gain"]!=0 for f in donor["folds"])
            or late["mean_gain"]<=0):
        raise ValueError("Run-specific narrative requires the documented failed gate and donor-only sensitivity pattern")
    observed=[f["gain"] for f in hct["primary"]["folds"]]+[
        hct["primary"]["mean_gain"],hct["decision"]["mean_fixed_predictor_one_percent_lower_gain"],
        donor["mean_gain"],late["mean_gain"],sw["primary"]["mean_gain"],temporal["equal_home_mean_gain"]]+removed
    earlier=[.14424247849217992,.04254064448634634,.9093840407316872,.14649162119995474,
             .31066469622754206,.3032144246585811,0.,.3114309239165016,0.,1.4803276373883425,
             .00457817568539598,.003465069460411585,.0033411302679190458,.0019795896037803073]
    comparison_error=max(abs(a-b) for a,b in zip(observed,earlier,strict=True))
    summary={
        "review_status":"unreviewed","reconstruction_with_known_targets":True,
        "earlier_reported_history_recovered":False,
        "earlier_reported_headline_values_reproduced":comparison_error<1e-10,
        "maximum_absolute_difference_from_earlier_reported_headline_values":comparison_error,
        "primary_hct_verdict":hct["decision"]["verdict"],
        "hct_mean_gain_nats_per_read":hct["primary"]["mean_gain"],
        "hct_primary_fold_gains":[f["gain"] for f in hct["primary"]["folds"]],
        "hct_target_only_one_percent_mean_lower_gain":hct["decision"]["mean_fixed_predictor_one_percent_lower_gain"],
        "exploratory_donor_only_frequency_floor":1e-5,
        "exploratory_donor_only_mean_gain":donor["mean_gain"],
        "exploratory_donor_removed_mass_fraction_range":[min(removed),max(removed)],
        "exploratory_late_only_same_floor_mean_gain":late["mean_gain"],
        "sw_mean_incremental_gain":sw["primary"]["mean_gain"],
        "exploratory_temporal_conditional_mean_gain":temporal["equal_home_mean_gain"],
        "exploratory_temporal_positive_home_count":sum(h.get("mean_gain",0)>0 for h in temporal["home_results"]),
        "exploratory_temporal_home_count":temporal["available_home_cultures"],
        "exploratory_temporal_single_occupied_halves":sum(f["receiver_test_occupied_barcodes"]==1 for f in halves),
        "exploratory_temporal_evaluated_halves_not_independent":len(halves),
        "independent_checks":len(verification["checks"]),
        "independent_maximum_absolute_numeric_difference":max(c["maximum_absolute_numeric_difference"] for c in verification["checks"]),
        "interpretation":"The specified incremental predictor fails its recorded robustness gate and is sensitive to a small tail of donor measurements. Technical spillover is plausible, not uniquely established. No new biological mechanism or clinical result is claimed.",
    }
    write_json(args.out/"summary.json",summary)
    files=sorted(HERE.glob("*.py"))+sorted(HERE.glob("*.md"))+sorted((HERE/"data").glob("*.json"))+[HERE/"requirements.txt"]
    write_json(args.out/"code_fingerprints.json",{
        "reconstruction_design_commit":PROTOCOL_COMMIT,
        "original_local_commits_not_recovered":True,
        "files":[{"path":str(p.relative_to(ROOT)),"sha256":digest(p)} for p in files]})
    report=(HERE/"README.md").read_text()
    def relocate(match):
        label,target=match.groups()
        if "://" in target or target.startswith("#"):
            return match.group(0)
        path,separator,anchor=target.partition("#")
        moved=Path(os.path.relpath((HERE/path).resolve(),args.out)).as_posix()
        return f"[{label}]({moved}{separator}{anchor})"
    (args.out/"research_report.md").write_text(re.sub(r"\[([^\]]+)\]\(([^)]+)\)",relocate,report))
    source="whiting-barcode-2025"
    article="https://doi.org/10.1038/s41467-025-59479-7"
    evidence=[{"source_id":source,"record_url":article,"relation":"context"}]
    findings=[
        {"id":"primary-gate-failed","statement":f"Reconstructed HCT116 nested forecasts improve all four held-out cultures by mean {hct['primary']['mean_gain']:.10f} nats/read, but the recorded descriptive robustness gate fails its required abundance-floor refit.",
         "evidence":evidence,"review_status":"unreviewed","limitations":"Newly calculated known-result reconstruction. Four cultures share a parent; reads are not independent cells. Source filtering and missing full parental/untreated controls constrain interpretation."},
        {"id":"rare-donor-predictor-sensitivity","statement":f"Exploratory donor-only filtering below frequency 1e-5 removes {100*min(removed):.8f}–{100*max(removed):.8f}% of early donor read mass and eliminates every incremental gain while preserving own early inputs, later outcomes, baseline predictions and dictionary.",
         "evidence":evidence,"review_status":"unreviewed","limitations":"The sensitivity of this specified fitted predictor does not prove that the removed reads are artifacts or that shared biological information is absent. This is an exploratory follow-up, not a new blinded primary test."},
        {"id":"outcome-check-does-not-protect-predictors","statement":f"The fixed-predictor adverse 1% outcome-mass deletion lower mean remains {hct['decision']['mean_fixed_predictor_one_percent_lower_gain']:+.10f} nats/read. An explicit synthetic construction with a sufficiently small predictor probability floor shows why arbitrarily small donor contamination can evade every target-only deletion check below one.",
         "evidence":evidence,"review_status":"unreviewed","limitations":"The construction is probability algebra, not an empirical contamination estimate. A fixed-floor continuity bound states its limit. Sparse temporal resemblance and first-record flowcell support do not uniquely identify spillover causation or reject the source mechanistic analysis."},
    ]
    descriptions={"summary.json":"Fresh reconstruction headline calculations and explicit limits",
                  "code_fingerprints.json":"Fingerprints of reconstructed code, design and source metadata",
                  "reproduce.py":"Run-specific authenticated source replay and artifact integrity entry point",
                  "research_report.md":"Substantive original report of the reconstructed computational investigation",
                  "forecast_robustness.png":"Original scientific chart of predictor and outcome sensitivities",
                  "forecast_robustness.svg":"Scalable original chart of the freshly calculated comparisons",
                  "predictor_counterexamples.json":"Explicitly synthetic probability examples and limiting calculations",
                  "independent_verification.json":"Separate probability-space and count-weighted empirical replay"}
    names=list(REPORTS)+["summary.json","code_fingerprints.json","research_report.md","forecast_robustness.png","forecast_robustness.svg","reproduce.py"]
    manifest={
        "schema_version":1,"task_id":"barcode-lineage-memory","run_id":"run-20261128",
        "submitted_by":{"type":"agent","name":"Open Cancer computational reconstruction","model_or_tool":"GPT-6 Astra Pro; separate Python mathematical implementations"},
        "created_at":datetime.now(timezone.utc).isoformat(),
        "inputs":[{"source_id":source,"accessed_on":"2026-11-28","url":article,
                   "query_or_filter":"Publisher Supplementary Data 2 HCTbc and SW6bc paired treated-culture barcode count matrices; original downloads untracked.",
                   "source_release":"Workbook SHA256 5c61f93a15134b795e6655ecab18f2d30e7d52daed33168ca0a60406303ee100"}],
        "execution":{"code_path":"reproduce.py","command":"python results/barcode-lineage-memory/run-20261128/reproduce.py --fetch",
                     "environment":"Python3.12; NumPy>=1.24,<3; OpenPyXL>=3.1,<4 for source preparation; Matplotlib optional for figures"},
        "artifacts":[{"path":name,"sha256":digest(args.out/name),"description":descriptions.get(name,"Original aggregate per-culture calculations and all declared sensitivity outcomes")} for name in names],
        "findings":findings,
        "limitations":"Unreviewed nonclinical reconstruction with prior reported targets known. The previous local history and logs were unavailable; all evidence in this bundle was recomputed. No source tables, barcode identities, raw reads or patient records are redistributed. No biological mechanism, historical priority, exact hopping rate, causal source error or treatment benefit is established. Negative and unavailable results remain explicit.",
    }
    write_json(args.out/"manifest.json",manifest)
    print(json.dumps({"artifacts":len(names),"verdict":summary["primary_hct_verdict"],"prior_headlines_reproduced":summary["earlier_reported_headline_values_reproduced"]}))


if __name__ == "__main__":
    main()
