"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const {encodeLedger, decodeLedger} = require("../scripts/catalog_ledger.cjs");

test("large, public metadata observations round-trip under the GitHub issue body limit", () => {
  const records = Object.fromEntries(Array.from({length: 125}, (_, i) => [
    "cbioportal-example-" + i,
    {
      title: "Open public TCGA study metadata record " + i + " — cohort",
      accession: {value: "tcga_example_" + i},
      citation_url: "https://pubmed.ncbi.nlm.nih.gov/?term=" +
        Array.from({length: 20}, (_, j) => "29625048%5Buid%5D+OR+" + j).join(""),
      retrieval_url: "https://www.cbioportal.org/api/studies/" + i + "?projection=DETAILED",
      release: {label: "Mutable source import, public metadata only"},
      modalities: ["genomics", "transcriptomics", "proteomics"],
    },
  ]));
  const report = {
    checked_at: "2026-11-10T20:00:00+00:00",
    sources: {cbioportal: {ok: true, last_healthy: "2026-11-10", observations: records}},
    alerts: [],
  };
  const encoded = encodeLedger(report);
  assert.ok(encoded.length < 50_000, "Issue marker should leave space for alert summaries");
  assert.deepEqual(decodeLedger("## Source health\n" + encoded), report);
});

test("no ledger and legacy plain-JSON ledger decode safely", () => {
  assert.equal(decodeLedger(null), null);
  assert.equal(decodeLedger("No previous ledger"), null);
  assert.deepEqual(decodeLedger("<!-- catalog-health-json\n{\"sources\":{\"gdc\":{}}}\n-->"),
                   {sources: {gdc: {}}});
});

test("malformed compressed ledger fails rather than silently treating it as healthy", () => {
  assert.throws(() => decodeLedger("<!-- catalog-health-gzip-base64\nYWJj\n-->"),
                /incorrect header check|unknown compression method|unexpected end|invalid/i);
});
