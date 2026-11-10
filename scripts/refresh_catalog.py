#!/usr/bin/env python3
"""Discover public TCGA project/study metadata and CIViC nightly export headers."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys
import tempfile
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "catalog" / "datasets"
GDC = "https://api.gdc.cancer.gov/projects"
CBIO = "https://www.cbioportal.org/api/studies"
CIVIC = "https://civicdb.org/downloads/nightly/"
GDC_FIELDS = "project_id,name,program.name,summary.data_categories.data_category"
PANCAN = "_tcga_pan_can_atlas_2018"
EDITORIAL = {
    "gdc": "gdc-tcga-luad-project",
    "cbioportal": "cbioportal-luad-tcga-pancan-atlas-2018",
    "civic": "civic-nightly-clinical-evidence",
}
# Filename stems, catalog names, and identifying TSV fields in each release.
EXPORTS = {
    "ClinicalEvidenceSummaries": ("clinical-evidence", "clinical evidence", "evidence_id"),
    "MolecularProfileSummaries": ("molecular-profiles", "molecular profile", "molecular_profile_id"),
    "AssertionSummaries": ("assertions", "assertion", "assertion_id"),
}


MAX_JSON_BYTES = 16 * 1024 * 1024


def get_json(url):
    req = Request(url, headers={"Accept": "application/json", "User-Agent": "open-cancer-catalog/1"})
    with urlopen(req, timeout=45) as response:
        payload = response.read(MAX_JSON_BYTES + 1)
    if len(payload) > MAX_JSON_BYTES:
        raise ValueError(f"oversized metadata response at {url}")
    return json.loads(payload)


def cbio_total():
    req = Request(f"{CBIO}?projection=META", headers={"Accept": "application/json",
                                                       "User-Agent": "open-cancer-catalog/1"})
    with urlopen(req, timeout=45) as response:
        value = response.headers.get("total-count")
    if not value or not value.isdigit():
        raise ValueError("cBioPortal did not report its study count")
    return int(value)


def civic_header(url):
    req = Request(url, headers={"User-Agent": "open-cancer-catalog/1"})
    with urlopen(req, timeout=45) as response:
        header = response.readline(8192)
    if not header.endswith(b"\n"):
        raise ValueError(f"missing or oversized CIViC TSV header at {url}")
    return header.decode("utf-8-sig").strip().split("\t")


def base_record(source, date):
    # Reuse the curated source access, licensing and documentation fields.
    template = json.loads((DATASETS / f"{EDITORIAL[source]}.json").read_text(encoding="utf-8"))
    return {k: template[k] for k in ("source_id", "documentation_url", "access", "reuse", "release")} | {
        "verified_on": date,
    }


def gdc_record(project, date, release):
    accession = project["project_id"]
    if not re.fullmatch(r"TCGA-[A-Z0-9]+", accession) or project["program"]["name"] != "TCGA":
        raise ValueError(f"unexpected TCGA project: {accession!r}")
    disease = project["name"]
    categories = {item["data_category"] for item in project["summary"]["data_categories"]}
    modality_categories = {
        "genomics": {"Simple Nucleotide Variation", "Copy Number Variation", "DNA Methylation",
                     "Somatic Structural Variation", "Structural Variation"},
        "clinical": {"Clinical"},
        "transcriptomics": {"Transcriptome Profiling"},
        "proteomics": {"Proteome Profiling"},
        "imaging": {"Imaging"},
    }
    modalities = [k for k, values in modality_categories.items() if categories & values]
    if not disease or not modalities:
        raise ValueError(f"incomplete project metadata for {accession}")
    query = urlencode({"fields": GDC_FIELDS})
    return base_record("gdc", date) | {
        "id": f"gdc-{accession.lower()}-project",
        "accession": {"namespace": "NCI GDC project ID", "value": accession},
        "title": f"{disease} ({accession}) in NCI Genomic Data Commons",
        "description": f"GDC project metadata for {disease} in The Cancer Genome Atlas. "
                       "The project summary identifies the available data categories.",
        "landing_url": f"https://portal.gdc.cancer.gov/projects/{accession}",
        "retrieval_url": f"{GDC}/{accession}?{query}",
        "release": {"label": f"{release} at retrieval; live GDC project metadata", "mutable": True},
        "disease": {"name": disease},
        "modalities": modalities,
        "unit_of_observation": "One GDC project",
        "cohort_definition": f"All data indexed under GDC project {accession}. "
                             "File access and assay coverage vary within the project.",
    }


def gdc_records(accession, date):
    if accession:
        if not re.fullmatch(r"TCGA-[A-Z0-9]+", accession):
            raise ValueError("expected a TCGA project accession such as TCGA-BRCA")
        query = urlencode({"fields": GDC_FIELDS})
        payload = get_json(f"{GDC}/{accession}?{query}")
        projects = [payload["data"]]
        if projects[0]["project_id"] != accession:
            raise ValueError("GDC returned a different project accession")
    else:
        filters = {"op": "in", "content": {"field": "program.name", "value": ["TCGA"]}}
        query = urlencode({"filters": json.dumps(filters, separators=(",", ":")),
                           "fields": GDC_FIELDS, "size": 100, "format": "json"})
        payload = get_json(f"{GDC}?{query}")
        data = payload["data"]
        projects = data["hits"]
        if len(projects) != data["pagination"]["total"]:
            raise ValueError("incomplete GDC TCGA project list")
    if payload.get("warnings"):
        raise ValueError(f"GDC returned query warnings: {payload['warnings']}")
    status = get_json("https://api.gdc.cancer.gov/status")
    release = status["data_release"]
    if status["status"] != "OK" or not release:
        raise ValueError("GDC did not report a current release")
    return [gdc_record(project, date, release) for project in projects]


def cbio_record(study, date):
    accession = study["studyId"]
    if not re.fullmatch(r"[a-z0-9]+" + PANCAN, accession) or study["publicStudy"] is not True:
        raise ValueError(f"expected a public PanCancer Atlas study, received {accession!r}")
    fields = {
        "genomics": ["sequencedSampleCount", "cnaSampleCount"],
        "transcriptomics": ["mrnaRnaSeqSampleCount", "mrnaRnaSeqV2SampleCount", "mrnaMicroarraySampleCount"],
        "proteomics": ["rppaSampleCount", "massSpectrometrySampleCount"],
    }
    modalities = [k for k, counts in fields.items() if any(study[count] > 0 for count in counts)]
    if not modalities or not study["cancerType"]["name"] or not study["name"]:
        raise ValueError(f"incomplete study metadata for {accession}")
    record = base_record("cbioportal", date) | {
        "id": "cbioportal-" + accession.removesuffix(PANCAN) + "-tcga-pancan-atlas-2018",
        "accession": {"namespace": "cBioPortal study ID", "value": accession},
        "title": study["name"] + " in cBioPortal",
        "description": f"Public study metadata for {study['name']}. Modality tags reflect "
                       "assays with a positive profiled-sample count in the study summary.",
        "landing_url": f"https://www.cbioportal.org/study/summary?id={accession}",
        "retrieval_url": f"{CBIO}/{accession}?projection=DETAILED",
        "disease": {"name": study["cancerType"]["name"]},
        "modalities": modalities,
        "unit_of_observation": "One cBioPortal study",
        "cohort_definition": f"Samples included by cBioPortal in {accession}. "
                             "Each molecular profile has its own set of assayed samples.",
        "release": {"label": f"TCGA PanCancer Atlas study, portal import {study['importDate']}",
                    "mutable": True},
    }
    if study.get("pmid"):
        pmids = list(dict.fromkeys(x.strip() for x in study["pmid"].split(",")))
        if not all(re.fullmatch(r"[1-9][0-9]*", pmid) for pmid in pmids):
            raise ValueError(f"invalid publication identifiers for {accession}")
        if len(pmids) == 1:
            record["citation_url"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmids[0]}/"
        else:
            query = urlencode({"term": " OR ".join(f"{pmid}[uid]" for pmid in pmids)})
            record["citation_url"] = f"https://pubmed.ncbi.nlm.nih.gov/?{query}"
    return record


def cbio_records(accession, date):
    if accession:
        if not re.fullmatch(r"[a-z0-9]+" + PANCAN, accession):
            raise ValueError(f"expected a cBioPortal study ID ending in {PANCAN}")
        study = get_json(f"{CBIO}/{accession}?projection=DETAILED")
        if study["studyId"] != accession:
            raise ValueError("cBioPortal returned a different study accession")
        return [cbio_record(study, date)]
    # Fetch study summaries with aggregate assay counts, never sample records.
    # One bounded page also avoids duplicates or omissions across a changing sort.
    limit = 2000
    studies = get_json(f"{CBIO}?projection=DETAILED&pageSize={limit}")
    if not isinstance(studies, list) or len(studies) >= limit or len(studies) != cbio_total():
        raise ValueError("incomplete cBioPortal study list")
    if len({s["studyId"] for s in studies}) != len(studies):
        raise ValueError("duplicate cBioPortal study in upstream collection")
    return [cbio_record(s, date) for s in studies if s["studyId"].endswith(PANCAN)]


def civic_records(accession, date):
    records = []
    for stem, (slug, unit, id_field) in EXPORTS.items():
        filename = f"nightly-{stem}.tsv"
        if accession and accession != filename:
            continue
        url = CIVIC + filename
        header = civic_header(url)
        if not {id_field, "last_review_date", "is_flagged"}.issubset(header):
            raise ValueError(f"unexpected columns for {filename}")
        records.append(base_record("civic", date) | {
            "id": f"civic-nightly-{slug}",
            "accession": {"namespace": "CIViC export filename", "value": filename},
            "title": f"CIViC nightly {unit} summaries",
            "description": f"Nightly CIViC {unit} table. The export header identifies "
                           f"records with {id_field} and includes flags and review dates.",
            "landing_url": "https://civicdb.org/releases",
            "retrieval_url": url,
            "disease": {"name": "Multiple cancer types (pan-cancer)"},
            "modalities": ["variant-evidence"],
            "unit_of_observation": f"One CIViC {unit} per row",
            "cohort_definition": f"All {unit} records included in the named nightly CIViC "
                                 "export, across diseases and molecular profiles.",
        })
    return records


CONNECTORS = {"gdc": gdc_records, "cbioportal": cbio_records, "civic": civic_records}


def discover(source=None, accession=None):
    if accession and not source:
        raise ValueError("select --source when looking up an accession")
    records = {}
    date = datetime.now(timezone.utc).date().isoformat()
    for name in ([source] if source else CONNECTORS):
        batch = CONNECTORS[name](accession, date)
        if not batch:
            raise ValueError(f"no matching {name} accessions")
        for record in batch:
            key = record["id"]
            if key in records:
                raise ValueError(f"duplicate upstream accession for {key}")
            records[key] = record
    return dict(sorted(records.items()))


def save_new_directory(output, records):
    if output.exists():
        raise ValueError(f"choose a new output directory: {output}")
    # Make the complete batch visible only after every record has been written.
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        stage = Path(tmp) / "records"
        stage.mkdir()
        for key, record in records.items():
            (stage / f"{key}.json").write_text(json.dumps(record, ensure_ascii=False, indent=2)
                                               + "\n", encoding="utf-8")
        stage.rename(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=CONNECTORS)
    parser.add_argument("--accession", help="Exact upstream project ID, study ID, or export filename")
    parser.add_argument("--output", type=Path, help="New directory for reviewable dataset records")
    args = parser.parse_args()
    try:
        records = discover(args.source, args.accession)
        if args.output:
            save_new_directory(args.output, records)
            print(f"Wrote {len(records)} upstream metadata records to {args.output}", file=sys.stderr)
        else:
            print(json.dumps(list(records.values()), indent=2, ensure_ascii=False))
    except (ValueError, OSError, TypeError, KeyError) as exc:
        parser.exit(1, f"catalog refresh failed: {exc}\n")


if __name__ == "__main__":
    main()
