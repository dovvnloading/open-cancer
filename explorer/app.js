import {filterDatasets, pageFor, permalink, pilotMeasures, reviewsFor, searchText} from "./catalog.mjs";

const main = document.querySelector("main");
const params = new URLSearchParams(location.search);
const project = "https://github.com/dovvnloading/open-cancer";
const pilotTask = "egfr-luad-evidence-coverage";
const pilotRun = "pilot-20261110";
const pilotDatasets = [
  "gdc-tcga-luad-project", "gdc-tcga-luad-egfr-ssm-occurrences",
  "cbioportal-luad-tcga-pancan-atlas-2018", "civic-nightly-clinical-evidence"
];
const pageSize = 12;
const portalNames = {gdc: "GDC", civic: "CIViC", cbioportal: "cBioPortal"};
let index;
let sources;

function el(tag, content, className) {
  const node = document.createElement(tag);
  if (content !== undefined) node.textContent = content;
  if (className) node.className = className;
  return node;
}

function link(label, href, className) {
  const node = el("a", label, className);
  const url = new URL(href, location.href);
  if (url.protocol === "https:" || (url.protocol === "http:" && url.origin === location.origin)) {
    node.href = url.href;
  }
  return node;
}

function group(tag, className, ...children) {
  const node = el(tag, undefined, className);
  node.append(...children);
  return node;
}

function list(items, className = "link-list") {
  return group("ul", className, ...items.map(item => group("li", "", item)));
}

function titleCase(value) {
  const text = value.replaceAll("-", " ");
  return text.charAt(0).toUpperCase() + text.slice(1);
}

function sourceName(id) {
  return portalNames[id] || sources.get(id)?.name || id;
}

function coveragePanel() {
  const counts = new Map();
  for (const row of index.datasets) counts.set(row.source_id, (counts.get(row.source_id) || 0) + 1);
  const entries = [...counts].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
  const heading = group("div", "coverage-title", el("strong", index.datasets.length, "catalog-number"), "indexed datasets");
  const rows = el("ol");
  for (const [id, count] of entries) {
    const bar = el("i");
    bar.style.width = (100 * count / index.datasets.length) + "%";
    const track = group("span", "track", bar);
    track.setAttribute("aria-hidden", "true");
    rows.append(group("li", "", link(sourceName(id), "?source=" + encodeURIComponent(id)), track, el("span", count, "count")));
  }
  const panel = group("div", "coverage", heading, rows);
  panel.setAttribute("aria-label", "Dataset counts by source");
  return panel;
}

function catalogHero() {
  const content = group("div", "", el("h1", "Open Cancer Explorer"),
    el("p", "Find public cancer research datasets by study, disease or accession. Inspect the cohort, access conditions and original source for every record."));
  return group("div", "hero", group("div", "wrap hero-inner", content, coveragePanel()));
}

function tableHead(labels) {
  const row = el("tr");
  for (const label of labels) {
    const cell = el("th", label);
    cell.scope = "col";
    row.append(cell);
  }
  return group("thead", "", row);
}

function selectField(name, label, values, defaultLabel) {
  const control = el("select");
  control.name = name;
  control.add(new Option(defaultLabel, ""));
  for (const [value, text] of values) control.add(new Option(text, value));
  const initial = params.get(name) || "";
  if (initial && !values.some(([value]) => value === initial)) control.add(new Option(initial, initial));
  control.value = initial;
  return group("label", "", label, control);
}

function catalog() {
  document.querySelector("#catalog-nav").setAttribute("aria-current", "page");
  const rows = [...index.datasets].sort((a, b) => a.title.localeCompare(b.title, "en") || a.id.localeCompare(b.id));
  const searchable = new Map(rows.map(row => [row.id, searchText(row, sources.get(row.source_id)?.name)]));
  const sourceIds = [...new Set(rows.map(row => row.source_id))].sort();
  const modalities = [...new Set(rows.flatMap(row => row.modalities))].sort();
  const levels = [...new Set(rows.map(row => row.access.level))].sort();
  const search = el("input");
  search.type = "search";
  search.name = "q";
  search.placeholder = "e.g. EGFR, TCGA-LUAD, DOID:3910";
  search.value = params.get("q") || "";
  const form = group("form", "filters", group("label", "search-field", "Search datasets", search),
    selectField("source", "Source", sourceIds.map(id => [id, sourceName(id)]), "All sources"),
    selectField("modality", "Modality", modalities.map(m => [m, titleCase(m)]), "All modalities"),
    selectField("access", "Access", levels.map(a => [a, titleCase(a)]), "All access levels"));
  form.setAttribute("role", "search");
  const count = el("span");
  count.setAttribute("role", "status");
  const clear = el("button", "Clear filters", "text-button");
  clear.type = "button";
  const table = group("table", "catalog-table", tableHead(["Dataset / accession", "Source", "Disease", "Access"]), el("tbody"));
  table.setAttribute("aria-label", "Matching datasets");
  const empty = group("div", "empty", el("h3", "No datasets match"), el("p", "Try another accession or remove a filter."));
  const pager = el("div", undefined, "pager");
  const summary = group("div", "result-summary", count, clear);
  main.replaceChildren(catalogHero(), group("section", "wrap catalog-body",
    el("h2", "Dataset catalog", "section-head"), form, summary, table, empty, pager));

  function render() {
    const filters = Object.fromEntries(new FormData(form));
    const matches = filterDatasets(rows, searchable, filters);
    const pages = Math.max(1, Math.ceil(matches.length / pageSize));
    const page = pageFor(params.get("page"), pages);
    if (page > 1) params.set("page", page);
    else params.delete("page");
    history.replaceState(null, "", location.pathname + (params.size ? "?" + params : ""));
    const start = (page - 1) * pageSize;
    const end = Math.min(start + pageSize, matches.length);
    count.textContent = matches.length > pageSize
      ? `Showing ${start + 1}–${end} of ${matches.length} datasets`
      : `${matches.length} dataset${matches.length === 1 ? "" : "s"}`;
    clear.hidden = !Object.values(filters).some(Boolean);
    table.hidden = !matches.length;
    empty.hidden = !!matches.length;
    const body = table.tBodies[0];
    body.replaceChildren();
    for (const row of matches.slice(start, end)) {
      const title = group("td", "", link(row.title, permalink(params, "dataset", row.id)), el("code", row.accession.value, "accession"));
      body.append(group("tr", "", title, el("td", sourceName(row.source_id)), el("td", row.disease.name),
        group("td", "", el("span", titleCase(row.access.level), "access " + row.access.level))));
    }
    pager.hidden = pages === 1;
    const changePage = (label, target) => {
      const button = el("button", label, "button");
      button.disabled = target < 1 || target > pages;
      button.onclick = () => {
        params.set("page", target);
        render();
        summary.tabIndex = -1;
        summary.focus();
      };
      return button;
    };
    pager.replaceChildren(changePage("Previous", page - 1), el("span", `Page ${page} of ${pages}`), changePage("Next", page + 1));
  }
  function update() {
    for (const [key, value] of new FormData(form)) {
      if (value) params.set(key, value);
      else params.delete(key);
    }
    params.delete("page");
    render();
  }
  form.addEventListener("input", update);
  form.addEventListener("submit", event => { event.preventDefault(); update(); });
  clear.onclick = () => {
    for (const field of form.elements) field.value = "";
    update();
    search.focus();
  };
  render();
}

function pageHead(title, description, meta) {
  document.title = title + " | Open Cancer Explorer";
  const inner = group("div", "wrap", link("← Dataset catalog", permalink(params), "back"),
    el("h1", title), el("p", description));
  if (meta) inner.append(el("div", meta, "meta"));
  return group("div", "page-head", inner);
}

function section(title, ...content) {
  return group("section", "", el("h2", title), ...content);
}

function facts(rows) {
  return group("dl", "fact-list", ...rows.map(([title, value]) =>
    group("div", "", el("dt", title), el("dd", value))));
}

function dataset(id) {
  document.querySelector("#catalog-nav").setAttribute("aria-current", "page");
  const row = index.datasets.find(record => record.id === id);
  if (!row) {
    main.replaceChildren(pageHead("Dataset unavailable", "This catalog has no dataset with this identifier.", id));
    return;
  }
  const disease = [row.disease.name, row.disease.ontology_id].filter(Boolean).join(" · ");
  const accessFacts = [
    ["Access level", titleCase(row.access.level)], ["Access conditions", row.access.notes],
    ["Reuse status", row.reuse.license_id || titleCase(row.reuse.status)], ["Reuse conditions", row.reuse.notes],
    ["Release", row.release.label], ["Release type", row.release.mutable ? "Mutable" : "Fixed"],
    ["Metadata verified", row.verified_on]
  ];
  const summary = group("div", "",
    section("Dataset description", facts([["Disease", disease], ["Modality", row.modalities.map(titleCase).join(", ")],
      ["Accession namespace", row.accession.namespace], ["Cohort definition", row.cohort_definition],
      ["Unit of observation", row.unit_of_observation]])),
    section("Access and provenance", facts(accessFacts)));
  const links = [link("Upstream dataset", row.landing_url), link("Source documentation", row.documentation_url),
    link("Retrieval endpoint", row.retrieval_url), link("Reuse terms", row.reuse.terms_url), link("Catalog record", row.record_url)];
  if (row.citation_url) links.push(link("Study citation", row.citation_url));
  const aside = group("aside", "side", el("h2", "Source links"), list(links));
  if (pilotDatasets.includes(id)) {
    aside.append(section("Related research", link("EGFR lung adenocarcinoma pilot", "?view=pilot")));
  }
  main.replaceChildren(pageHead(row.title, row.description, sourceName(row.source_id) + "  /  " + row.accession.value),
    group("div", "wrap page-body", summary, aside));
}

function renderPilot() {
  document.querySelector("#pilot-nav").setAttribute("aria-current", "page");
  const question = index.questions.find(task => task.id === pilotTask);
  const pilot = index.pilots.find(run => run.task_id === pilotTask && run.run_id === pilotRun);
  if (!question || !pilot || !pilot.coverage) {
    main.replaceChildren(pageHead("Pilot snapshot unavailable", "Consult the research question and its version history in the repository."),
      group("div", "wrap loading", link("Open research question", project + "/issues/3")));
    return;
  }
  const coverage = pilot.coverage;
  const body = el("tbody");
  for (const [name, value, unit] of pilotMeasures(coverage)) {
    body.append(group("tr", "", el("td", name), el("td", value), el("td", unit)));
  }
  const measures = group("table", "measure-table", tableHead(["Source", "Count", "Unit and selection"]), body);
  measures.setAttribute("aria-label", "Pilot snapshot source counts");
  const reviews = reviewsFor(pilot, index.reviews);
  const reviewInfo = reviews.length
    ? list(reviews.map(r => link(`${titleCase(r.scope)} · ${titleCase(r.decision)} · ${r.reviewer.name} (${r.reviewer.type})`, r.record_url)))
    : el("p", "No version-bound review has been submitted for this snapshot.");
  const observed = section("Recorded counts",
    el("p", "Each source counts a different kind of record. The cohorts and units do not support pooling these counts or estimating prevalence or treatment effects.", "lead-note"),
    measures);
  const reviewState = section("Review record", reviewInfo);
  if (!reviews.some(r => r.scope === "scientific-interpretation")) {
    reviewState.append(el("p", "Independent scientific interpretation is pending."));
  }
  reviewState.append(link("Scientific review protocol", project + "/blob/main/REVIEW_PROTOCOL.md"));
  const content = group("div", "", observed, section("Interpretation limits", list(coverage.limitations, "notes")), reviewState);
  const sourceLinks = pilot.artifacts.map(a => link(a.path, a.url));
  const aside = group("aside", "side", el("h2", "Reproduce the pilot"),
    el("p", "The manifest records the inputs, query filters, execution command and artifact checksums."),
    list([link("Run manifest", pilot.manifest_url), link("Analysis code", pilot.code_url),
      ...sourceLinks, link("Research question", question.record_url), link("Follow the research discussion", project + "/issues/3")]));
  const datasets = pilotDatasets.map(id => index.datasets.find(row => row.id === id)).filter(Boolean);
  if (datasets.length) aside.append(section("Dataset records",
    list(datasets.map(row => link(row.title, permalink(params, "dataset", row.id))))));
  const header = pageHead("EGFR in lung adenocarcinoma", question.question,
    [pilot.run_id, pilot.created_at.slice(0, 10), coverage.gene.ensembl_id, coverage.disease.doid].join("  /  "));
  header.firstChild.append(el("span", "Snapshot status · " + titleCase(coverage.status), "status"));
  main.replaceChildren(header, group("div", "wrap page-body", content, aside));
}

async function start() {
  try {
    const response = await fetch("./api/v1/index.json");
    if (!response.ok) throw Error("Catalog unavailable");
    index = await response.json();
    if (index.version !== 1 || !["datasets", "sources", "questions", "pilots", "reviews"].every(key => Array.isArray(index[key]))) {
      throw Error("Unrecognized index");
    }
    sources = new Map(index.sources.map(row => [row.id, row]));
    if (params.has("dataset")) dataset(params.get("dataset"));
    else if (params.get("view") === "pilot") renderPilot();
    else catalog();
  } catch {
    main.replaceChildren(group("div", "wrap loading", el("h1", "Explorer unavailable"),
      el("p", "The catalog could not be loaded. Reload the page or browse the source records."),
      link("Browse records on GitHub", project + "/tree/main/catalog/datasets")));
  } finally {
    main.setAttribute("aria-busy", "false");
  }
}

start();
