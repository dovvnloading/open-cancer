// Pure catalog operations shared by the Explorer and its regression tests.
export function permalink(params, key, value) {
  const query = new URLSearchParams(params);
  query.delete("dataset");
  query.delete("view");
  if (key) query.set(key, value);
  return "?" + query.toString();
}

export function searchText(row, sourceName) {
  return [row.id, row.title, row.description, row.accession.namespace, row.accession.value,
    row.disease.name, row.disease.ontology_id, row.cohort_definition,
    row.unit_of_observation, row.source_id, sourceName, ...row.modalities].join(" ").toLowerCase();
}

export function filterDatasets(rows, searchable, {q, source, modality, access}) {
  const words = q.toLowerCase().trim().split(/\s+/).filter(Boolean);
  return rows.filter(row => (!source || row.source_id === source)
    && (!access || row.access.level === access)
    && (!modality || row.modalities.includes(modality))
    && words.every(word => searchable.get(row.id).includes(word)));
}

export function pageFor(value, pages) {
  const page = Number(value);
  return Number.isSafeInteger(page) ? Math.min(pages, Math.max(1, page)) : 1;
}

export function pilotMeasures(coverage) {
  const gdc = coverage.gdc.egfr_occurrences;
  const values = [
    ["GDC", gdc.somatic_mutation_occurrence_rows, "Simple somatic mutation occurrence rows in TCGA-LUAD"],
    ["GDC", gdc.distinct_case_ids_with_recorded_occurrence, "Distinct case IDs among those occurrences"],
    ["cBioPortal", coverage.cbioportal.egfr_mutation_call_rows, "EGFR mutation-call rows in the PanCancer Atlas LUAD study"],
    ["CIViC", coverage.civic.egfr_accepted_evidence_exact_luad, "Accepted evidence items matching EGFR and the exact LUAD filter"]
  ];
  return values.map(([name, value, unit]) => [name, value ?? "Unavailable", unit]);
}

export function reviewsFor(pilot, reviews) {
  return reviews.filter(review => !review.synthetic && review.task_id === pilot.task_id
    && review.run_id === pilot.run_id && review.reviewed_manifest_sha256 === pilot.manifest_sha256);
}
