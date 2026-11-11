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
    row.unit_of_observation, row.source_id, sourceName, ...row.modalities,
    ...(row.research_context || []).flatMap(context =>
      [context.gene_symbol, context.disease_label, context.task_id])
  ].join(" ").toLowerCase();
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

/** Publication of an assessment is not verification of expertise or independence. */
export function scientificReviewSummary(reviews) {
  const expert = reviews.filter(r => r.scope === "scientific-interpretation" && r.reviewer?.type === "human");
  const disputed = expert.some(r => r.decision === "disputed");
  const revision = expert.some(r => r.decision === "needs-revision");
  const assessed = expert.some(r => r.decision === "scientific-assessment-completed"
    && ["source_provenance", "cohort_and_unit", "biological_interpretation", "uncertainty"]
      .every(key => ["passed", "concern"].includes(r.checks?.[key]?.status)));
  if (disputed) return {
    heading: "Scientific interpretation disputed",
    explanation: "An expert review disputes the analysis. Independent scientific validity is not established; inspect the conflicting review and its evidence."
  };
  if (revision) return {
    heading: "Scientific revisions requested",
    explanation: "An expert reviewer requested corrections. Independent scientific validity remains unestablished; inspect the review and outstanding work."
  };
  if (assessed) return {
    heading: "Scientific assessment submitted — verification pending",
    explanation: "A human scientific assessment is on record, but Open Cancer does not independently verify reviewer credentials or independence in this API. This is not scientific endorsement or clinical validation."
  };
  return {
    heading: "Independent scientific interpretation pending",
    explanation: "No independently credential-verified scientific interpretation is established. A computational methods reproduction does not change this status."
  };
}
