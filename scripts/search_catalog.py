#!/usr/bin/env python3
"""Offline, deterministic dataset discovery for humans and agents."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def search_datasets(records, *, source=None, disease=None, text=None):
    """Exact source/ontology filters and substring free-text matching."""
    result = []
    for obj in records:
        if source and obj.get("source_id") != source:
            continue
        dis = obj.get("disease", {})
        if disease and dis.get("ontology_id", "").casefold() != disease.casefold():
            continue
        if text:
            searchable = " ".join([obj.get("id", ""), obj.get("title", ""),
                                   obj.get("description", ""), dis.get("name", ""),
                                   obj.get("cohort_definition", "")]).casefold()
            if text.casefold() not in searchable:
                continue
        result.append(obj)
    return sorted(result, key=lambda row: row["id"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", help="Exact source ID, such as civic or gdc")
    parser.add_argument("--disease", help="Exact ontology ID, e.g. DOID:3910")
    parser.add_argument("--text", help="Case-insensitive substring")
    parser.add_argument("--format", choices=["json", "tsv"], default="tsv")
    args = parser.parse_args()
    records = [json.loads(path.read_text(encoding="utf-8"))
               for path in sorted((ROOT / "catalog/datasets").glob("*.json"))]
    matches = search_datasets(records, source=args.source, disease=args.disease,
                              text=args.text)
    if args.format == "json":
        print(json.dumps({"count": len(matches), "datasets": matches}, indent=2,
                         ensure_ascii=False))
    else:
        print("ID\tSOURCE\tACCESS\tTITLE\tLANDING_URL")
        for obj in matches:
            print("\t".join([obj["id"], obj["source_id"], obj["access"]["level"],
                             obj["title"], obj["landing_url"]]))


if __name__ == "__main__":
    main()
