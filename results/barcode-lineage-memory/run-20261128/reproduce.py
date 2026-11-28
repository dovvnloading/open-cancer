#!/usr/bin/env python3
"""Replay the reconstructed aggregate results from authenticated source data."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESEARCH = ROOT/"research/barcode-lineage-memory"
REPORTS = ("HCTbc.json","SW6bc.json","HCTbc_exploratory.json","SW6bc_exploratory.json","predictor_counterexamples.json")


def compare(actual, expected, path="root"):
    if isinstance(expected,dict):
        for key,value in expected.items():
            if key not in ("execution","cache_sha256"):
                compare(actual[key],value,path+"."+key)
    elif isinstance(expected,list):
        if len(actual) != len(expected):
            raise AssertionError(path+": changed list length")
        for i,(a,b) in enumerate(zip(actual,expected)):
            compare(a,b,f"{path}[{i}]")
    elif isinstance(expected,float):
        if not math.isclose(actual,expected,rel_tol=1e-9,abs_tol=1e-10):
            raise AssertionError(f"{path}: {actual} != {expected}")
    elif actual != expected:
        raise AssertionError(f"{path}: {actual!r} != {expected!r}")


def check_artifacts():
    manifest=json.loads((HERE/"manifest.json").read_text())
    for item in manifest["artifacts"]:
        path=(HERE/item["path"]).resolve()
        if not path.is_relative_to(HERE) or hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise AssertionError(f"Aggregate artifact fingerprint mismatch: {item['path']}")
    pins=json.loads((HERE/"code_fingerprints.json").read_text())
    for item in pins["files"]:
        path=(ROOT/item["path"]).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise AssertionError(f"Analysis fingerprint mismatch: {item['path']}")
    print(f"Matched {len(manifest['artifacts'])} aggregate artifacts and {len(pins['files'])} analysis/source fingerprints.",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir",type=Path,default=ROOT/".compute/barcode-source")
    parser.add_argument("--out",type=Path,default=ROOT/".compute/barcode-lineage-memory-replay")
    parser.add_argument("--fetch",action="store_true")
    parser.add_argument("--check-artifacts",action="store_true")
    args=parser.parse_args()
    args.cache_dir=args.cache_dir.expanduser().resolve()
    args.out=args.out.expanduser().resolve()
    check_artifacts()
    if args.check_artifacts:
        return
    if args.out.is_relative_to(ROOT) and not args.out.is_relative_to(ROOT/".compute"):
        raise ValueError("Replay output inside the checkout must remain in its untracked .compute directory")
    args.out.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1")
    def run(name,*arguments):
        subprocess.run([sys.executable,"-B",str(RESEARCH/name),*map(str,arguments)],cwd=ROOT,env=env,check=True)
    if args.fetch:
        run("fetch_source.py","--cache-dir",args.cache_dir)
    run("analysis.py","--cache-dir",args.cache_dir,"--out",args.out)
    run("exploratory.py","--cache-dir",args.cache_dir,"--out",args.out)
    run("predictor_counterexample.py","--output",args.out/"predictor_counterexamples.json")
    for name in REPORTS:
        compare(json.loads((args.out/name).read_text()),json.loads((HERE/name).read_text()),name)
    run("verifier.py","--cache-dir",args.cache_dir,"--reports-dir",args.out,
        "--exploratory-hct","--out",args.out/"independent_verification.json")
    print("Reproduced four empirical reports and synthetic examples; independent mathematical replay passed.")


if __name__ == "__main__":
    main()
