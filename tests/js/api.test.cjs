const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");
const { test } = require("node:test");
const path = require("node:path");

const source = fs.readFileSync(path.join(__dirname, "../../src/phytclust/gui/static/js/api.js"), "utf8")
  .replace(/^import .*$/gm, "").replace(/export /g, "");

function setup() {
  const elements = new Map();
  const element = (id, value = "") => {
    const node = { value, checked: false, hidden: true, disabled: false, innerHTML: "Run", children: [], replaceChildren() { this.children = []; this.innerHTML = ""; }, appendChild(child) { this.children.push(child); } };
    elements.set(id, node);
    return node;
  };
  const newickEl = element("newick", "(a:1,b:1);");
  element("extra-k", "2");
  const button = element("btn-run");
  const plot = element("optimalk_plot");
  const marker = element("results-stale-indicator");
  const mode = { value: "k" };
  const flags = { failParse: false, failTree: false, failLayouts: false };
  const messages = [];
  const state = { isRunning: false, latestApiData: null, latestOptimalKData: null, latestRunId: null, lastRunSignature: null, runHistory: [], CURRENT_CLUSTERS: {}, CLUSTER_COLORS: [], SEARCH_TERM: "" };
  const data = { newick: "(a:1,b:1);", mode: "k", run_id: "new", k_values: [2], clusters: [{ a: 0, b: 1 }], scores: null };
  const calls = { requests: 0, tree: 0, scores: 0 };
  const context = vm.createContext({
    state, newickEl, extraOutgroupEl: element("outgroup"), extraRootTaxonEl: element("root"), extraResolutionEl: element("resolution"),
    document: { getElementById: (id) => elements.get(id) || null, querySelector: () => ({ dataset: { mode: mode.value } }), createElement: () => ({}) },
    showStatus: (message, type) => messages.push({ message, type }), showToast: (message, type) => messages.push({ message, type }),
    readCheckParam: (id) => elements.get(id)?.checked || false,
    readSelectParam: (id) => elements.get(id)?.value || "",
    estimateLeafCount: () => 2,
    parseNewick: () => { if (flags.failParse) throw new Error("Bad returned tree"); return { name: "root" }; },
    accumulateBranchLength: (tree) => { tree._bl = 0; },
    generateClusterColors: (count) => Array.from({ length: count }, () => "blue"),
    computeLayouts: () => { if (flags.failLayouts) throw new Error("Layout failed"); state.HIER_CART = {}; state.HIER_CIRC = {}; },
    drawTree: () => { calls.tree++; if (flags.failTree) { flags.failTree = false; throw new Error("Tree drawing failed"); } },
    clearTree() {}, populateClusterSelector() {}, updateClusterEditorAvailability() {}, populateCompareSelectors() {},
    drawMiniScores() {}, drawOptimalK: () => { calls.scores++; },
    performance: { now: () => 0 }, console: { warn() {}, error() {} },
    fetch: async () => { calls.requests++; return { ok: true, status: 200, text: async () => JSON.stringify(data) }; },
  });
  vm.runInContext(source, context);
  return { context, state, data, calls, flags, elements, element, messages, mode, button, plot, marker };
}

function previousResult(fixture) {
  const old = { run_id: "old", newick: "old tree", clusters: [{ old: 0 }] };
  fixture.state.latestApiData = old;
  fixture.state.latestRunId = "old";
  fixture.state.NEWICK_RAW_TREE = { name: "old" };
  fixture.state.CURRENT_CLUSTERS = old.clusters[0];
  return old;
}

test("validation details become readable errors and preserve the previous run", async () => {
  const fixture = setup();
  const old = previousResult(fixture);
  fixture.context.fetch = async () => ({ ok: false, status: 422, text: async () => JSON.stringify({ detail: [{ loc: ["body", "top_n"], msg: "Must be positive" }] }) });
  await fixture.context.runPhytClust();
  assert.equal(fixture.state.latestApiData, old);
  assert.match(fixture.messages.at(-1).message, /top_n: Must be positive/);
  assert.equal(fixture.state.isRunning, false);
  assert.equal(fixture.button.innerHTML, "Run");
});

for (const response of ["<html>error</html>", "null", "{}", '{"newick":"tree","clusters":[{"a":"bad"}]}']) {
  test(`unreadable or incomplete response preserves previous result: ${response}`, async () => {
    const fixture = setup();
    const old = previousResult(fixture);
    fixture.context.fetch = async () => ({ ok: true, status: 200, text: async () => response });
    await fixture.context.runPhytClust();
    assert.equal(fixture.state.latestApiData, old);
    assert.equal(fixture.state.latestRunId, "old");
    assert.equal(fixture.state.runHistory.length, 0);
    assert.equal(fixture.button.disabled, false);
  });
}

test("fixed-k run clears previous score data", async () => {
  const fixture = setup();
  fixture.state.latestOptimalKData = { scores: [2, 4] };
  await fixture.context.runPhytClust();
  assert.equal(fixture.state.latestOptimalKData, null);
  assert.equal(fixture.state.latestRunId, "new");
  assert.equal(fixture.plot.children[0].textContent, "Scores are not calculated for a fixed cluster count.");
  assert.equal(fixture.marker.hidden, true);
});

test("disabled settings are ignored and fractional counts are rejected", () => {
  const fixture = setup();
  fixture.element("extra-topn", "invalid");
  fixture.element("extra-support-weight", "invalid");
  fixture.element("extra-prominence-weight", "invalid");
  assert.equal(fixture.context.buildRunRequest().k, 2);
  fixture.mode.value = "global";
  fixture.element("extra-topn", "1.5");
  assert.throws(() => fixture.context.buildRunRequest(), /Peak count must be a positive whole number/);
  fixture.element("extra-topn", "3");
  fixture.element("extra-ranking-mode", "raw");
  const request = fixture.context.buildRunRequest();
  assert.equal(request.prominence_weight, undefined);
  assert.equal(request.support_weight, undefined);
});

test("invalid active settings do not submit a request", async () => {
  const fixture = setup();
  fixture.element("extra-k", "1.5");
  await fixture.context.runPhytClust();
  assert.equal(fixture.calls.requests, 0);
  assert.match(fixture.messages[0].message, /Cluster count/);
});

test("outgroup validation is left to the server", async () => {
  const fixture = setup();
  fixture.elements.get("outgroup").value = "internal clade";
  await fixture.context.runPhytClust();
  assert.equal(fixture.calls.requests, 1);
  assert.equal(fixture.state.latestRunId, "new");
});

test("layout and tree failures restore result state", async () => {
  for (const failure of ["failLayouts", "failTree", "failParse"]) {
    const fixture = setup();
    const old = previousResult(fixture);
    fixture.flags[failure] = true;
    await fixture.context.runPhytClust();
    assert.equal(fixture.state.latestApiData, old);
    assert.equal(fixture.state.CURRENT_CLUSTERS, old.clusters[0]);
    assert.equal(fixture.state.runHistory.length, 0);
    assert.equal(fixture.state.isRunning, false);
  }
});

test("edits during a request leave the new result marked stale", async () => {
  const fixture = setup();
  fixture.context.fetch = async () => {
    fixture.elements.get("extra-k").value = "3";
    return { ok: true, status: 200, text: async () => JSON.stringify(fixture.data) };
  };
  await fixture.context.runPhytClust();
  assert.equal(fixture.marker.hidden, false);
});

test("large cluster maps avoid spreading values into function arguments", async () => {
  const fixture = setup();
  fixture.data.clusters = [Object.fromEntries(Array.from({ length: 150000 }, (_, index) => [`leaf${index}`, index % 3]))];
  await fixture.context.runPhytClust();
  assert.equal(fixture.state.latestRunId, "new");
  assert.equal(fixture.state.CLUSTER_COLORS.length, 3);
});

test("repeated clicks issue only one request", async () => {
  const fixture = setup();
  let finish;
  fixture.context.fetch = () => { fixture.calls.requests++; return new Promise((resolve) => { finish = resolve; }); };
  const first = fixture.context.runPhytClust();
  await fixture.context.runPhytClust();
  assert.equal(fixture.calls.requests, 1);
  finish({ ok: true, status: 200, text: async () => JSON.stringify(fixture.data) });
  await first;
  assert.equal(fixture.button.disabled, false);
});

test("no peaks is a successful empty selection and clears old scores", async () => {
  const fixture = setup();
  fixture.mode.value = "global";
  fixture.data.mode = "global";
  fixture.data.clusters = [];
  fixture.data.k_values = [];
  fixture.data.scores = [];
  await fixture.context.runPhytClust();
  assert.equal(fixture.state.latestOptimalKData, null);
  assert.equal(fixture.state.latestRunId, "new");
  assert.ok(fixture.messages.some(({ message, type }) => message.startsWith("No partitions selected") && type === "info"));
});
