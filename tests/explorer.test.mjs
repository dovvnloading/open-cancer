import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import test from "node:test";
import {filterDatasets, pageFor, permalink, pilotMeasures, reviewsFor, scientificReviewSummary, searchText} from "../explorer/catalog.mjs";

function json(path) {
  return JSON.parse(readFileSync(new URL("../" + path, import.meta.url), "utf8"));
}

const gdc = json("catalog/datasets/gdc-tcga-luad-egfr-ssm-occurrences.json");
const civic = json("catalog/datasets/civic-nightly-clinical-evidence.json");
const rows = [gdc, civic];
const searchable = new Map(rows.map(row => [row.id, searchText(row, json("catalog/sources/" + row.source_id + ".json").name)]));
const runDirectory = "results/egfr-luad-evidence-coverage/pilot-20261110/";

test("accession, ontology and source-name tokens combine with exact filters", () => {
  assert.deepEqual(filterDatasets(rows, searchable, {q: "doid:3910 ENSG00000146648 nci", source: "gdc", modality: "genomics", access: "mixed"}), [gdc]);
  assert.deepEqual(filterDatasets(rows, searchable, {q: "egfr", access: "open"}), []);
  assert.deepEqual(filterDatasets(rows, searchable, {q: "", source: "unknown-source"}), []);
  assert.deepEqual(filterDatasets(rows, searchable, {q: "", access: "open", modality: "variant-evidence"}), [civic]);
});

test("permalinks preserve catalog filters and remove competing views", () => {
  const params = new URLSearchParams({q: "EGFR & LUAD", source: "gdc", page: "2", view: "pilot"});
  const detail = new URLSearchParams(permalink(params, "dataset", gdc.id));
  assert.equal(detail.get("q"), "EGFR & LUAD");
  assert.equal(detail.get("dataset"), gdc.id);
  assert.equal(detail.has("view"), false);
  const back = new URLSearchParams(permalink(detail));
  assert.equal(back.has("dataset"), false);
  assert.equal(back.get("source"), "gdc");
  assert.equal(back.get("page"), "2");
  assert.equal(params.get("view"), "pilot");
});

test("pagination handles stale, malformed and out-of-range page numbers", () => {
  for (const value of [null, "", "bad", "-1", "2.5", "Infinity"]) assert.equal(pageFor(value, 3), 1);
  assert.equal(pageFor("2", 3), 2);
  assert.equal(pageFor("150", 3), 3);
});

test("occurrence rows and distinct cases have separate units and missing counts differ from zero", () => {
  const coverage = json(runDirectory + "coverage.json");
  coverage.gdc.egfr_occurrences.somatic_mutation_occurrence_rows = null;
  coverage.gdc.egfr_occurrences.distinct_case_ids_with_recorded_occurrence = null;
  coverage.cbioportal.egfr_mutation_call_rows = 0;
  const measures = pilotMeasures(coverage);
  assert.deepEqual(measures.map(row => row[1]), ["Unavailable", "Unavailable", 0, coverage.civic.egfr_accepted_evidence_exact_luad]);
  assert.match(measures[0][2], /occurrence rows/i);
  assert.match(measures[1][2], /distinct case IDs/i);
});

test("review display includes adverse decisions and excludes examples and stale hashes", () => {
  const pilot = {task_id: "test-task", run_id: "test-run", manifest_sha256: "a".repeat(64)};
  const disputed = {task_id: "test-task", run_id: "test-run", reviewed_manifest_sha256: pilot.manifest_sha256, synthetic: false, decision: "disputed"};
  assert.deepEqual(reviewsFor(pilot, [disputed, {...disputed, synthetic: true}, {...disputed, run_id: "other"}, {...disputed, reviewed_manifest_sha256: "b".repeat(64)}]), [disputed]);
});

// The pilot crosswalk is a curated link, not an inference that all source rows
// match EGFR or even describe LUAD in pan-cancer resources.
test("EGFR LUAD search discovers all four specifically linked source records", () => {
  const ids = [
    "gdc-tcga-luad-egfr-ssm-occurrences", "gdc-tcga-luad-project",
    "cbioportal-luad-tcga-pancan-atlas-2018", "civic-nightly-clinical-evidence"
  ];
  const records = ids.map(id => json("catalog/datasets/" + id + ".json"));
  const unrelated = json("catalog/datasets/gdc-tcga-kirc-project.json");
  const all = [...records, unrelated];
  const lookup = new Map(all.map(r => [r.id, searchText(r, r.source_id)]));
  assert.deepEqual(
    new Set(filterDatasets(all, lookup, {q:"EGFR LUAD"}).map(r=>r.id)),
    new Set(ids)
  );
  assert.equal(records[3].research_context[0].relationship, "upstream-source");
  assert.match(records[3].research_context[0].note, /pan-cancer/i);
});

test("critical scientific review never signals independent verification", () => {
  const base = {scope:"scientific-interpretation", reviewer:{type:"human"},
    checks:Object.fromEntries(["source_provenance","cohort_and_unit","biological_interpretation","uncertainty"]
      .map(k=>[k,{status:"passed"}]))};
  assert.match(scientificReviewSummary([]).heading, /pending/i);
  assert.match(scientificReviewSummary([{...base,decision:"needs-revision"}]).heading, /revisions requested/i);
  assert.match(scientificReviewSummary([{...base,decision:"disputed"}]).heading, /disputed/i);
  assert.match(scientificReviewSummary([{...base,decision:"scientific-assessment-completed"}]).heading, /verification pending/i);
  assert.match(scientificReviewSummary([{...base,reviewer:{type:"agent"},decision:"scientific-assessment-completed"}]).heading, /pending/i);
  assert.match(scientificReviewSummary([{...base,decision:"scientific-assessment-completed",
    checks: {biological_interpretation:{status:"not-assessed"}}}]).heading, /pending/i);
  assert.match(scientificReviewSummary([{...base,decision:"scientific-assessment-completed"},
    {...base,decision:"disputed"}]).heading, /disputed/i);
});
