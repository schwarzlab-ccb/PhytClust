const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const { test } = require('node:test');
const path = require('node:path');

const source = fs.readFileSync(path.join(__dirname, '../../src/phytclust/gui/static/js/api.js'), 'utf8')
  .replace(/^import .*$/gm, '').replace(/export /g, '');

function setup() {
  const elements = new Map();
  const element = (id, value = '') => {
    const node = { value, checked: false, hidden: true, disabled: false, innerHTML: 'Run', children: [], replaceChildren() { this.children = []; this.innerHTML = ''; }, appendChild(child) { this.children.push(child); } };
    elements.set(id, node);
    return node;
  };
  const newickEl = element('newick', '(a:1,b:1);');
  element('k', '2');
  const button = element('btn-run');
  const plot = element('scores-plot');
  const marker = element('results-stale-indicator');
  const mode = { value: 'k' };
  const flags = { failParse: false, failTree: false, failLayouts: false };
  const messages = [];
  const logEntries = [];
  const state = { isRunning: false, lastRunData: null, lastScoresData: null, lastRunId: null, lastRunSignature: null, runHistory: [], currentClusters: {}, clusterColours: [], searchTerm: '' };
  const data = { newick: '(a:1,b:1);', mode: 'k', run_id: 'new', k_values: [2], clusters: [{ a: 0, b: 1 }], scores: null };
  const calls = { requests: 0, tree: 0, scores: 0, miniScores: [] };
  const context = vm.createContext({
    state, newickEl, outgroupEl: element('outgroup'), rootTaxonEl: element('root'),
    document: { getElementById: (id) => elements.get(id) || null, querySelector: () => ({ dataset: { mode: mode.value } }), createElement: () => ({}) },
    showStatus: (message, type) => messages.push({ message, type }), showToast: (message, type) => messages.push({ message, type }),
    logStatus: (message, type) => logEntries.push({ message, type }),
    setCurrentClusters(map) {
      state.runClusters = map;
      state.currentClusters = map;
      let maximum = -1;
      for (const value of Object.values(map)) maximum = Math.max(maximum, value);
      state.clusterColours = Array.from({ length: maximum + 1 }, () => 'blue');
    },
    resetClusterSelector() { calls.selectorReset = (calls.selectorReset || 0) + 1; },
    drawComparison() { calls.comparison = (calls.comparison || 0) + 1; },
    d3: { zoomIdentity: 'identity' },
    readCheckParam: (id) => elements.get(id)?.checked || false,
    readSelectParam: (id) => elements.get(id)?.value || '',
    estimateLeafCount: () => 2,
    parseNewick: () => { if (flags.failParse) throw new Error('Bad returned tree'); return { name: 'root' }; },
    accumulateBranchLength: (tree) => { tree._bl = 0; },
    computeLayouts: () => { if (flags.failLayouts) throw new Error('Layout failed'); state.cartesianHierarchy = {}; state.circularHierarchy = {}; },
    drawTree: () => { calls.tree++; if (flags.failTree) { flags.failTree = false; throw new Error('Tree drawing failed'); } },
    clearTree() {}, populateClusterSelector() {}, updateClusterEditorAvailability() {}, populateCompareSelectors() {},
    drawMiniScores(data) { calls.miniScores.push(data); }, drawScores: () => { calls.scores++; },
    performance: { now: () => 0 }, console: { warn() {}, error() {} },
    fetch: async () => { calls.requests++; return { ok: true, status: 200, text: async () => JSON.stringify(data) }; },
  });
  vm.runInContext(source, context);
  return { context, state, data, calls, flags, elements, element, messages, logEntries, mode, button, plot, marker };
}

function previousResult(fixture) {
  const old = { run_id: 'old', newick: 'old tree', clusters: [{ old: 0 }] };
  fixture.state.lastRunData = old;
  fixture.state.lastRunId = 'old';
  fixture.state.newickTree = { name: 'old' };
  fixture.state.currentClusters = old.clusters[0];
  return old;
}

test('validation details become readable errors and preserve the previous run', async () => {
  const fixture = setup();
  const old = previousResult(fixture);
  fixture.context.fetch = async () => ({ ok: false, status: 422, text: async () => JSON.stringify({ detail: [{ loc: ['body', 'top_n'], msg: 'Must be positive' }] }) });
  await fixture.context.runPhytClust();
  assert.equal(fixture.state.lastRunData, old);
  assert.match(fixture.messages.at(-1).message, /top_n: Must be positive/);
  assert.equal(fixture.state.isRunning, false);
  assert.equal(fixture.button.innerHTML, 'Run');
});

for (const response of ['<html>error</html>', 'null', '{}', '{"newick":"tree","clusters":[{"a":"bad"}]}']) {
  test(`unreadable or incomplete response preserves previous result: ${response}`, async () => {
    const fixture = setup();
    const old = previousResult(fixture);
    fixture.context.fetch = async () => ({ ok: true, status: 200, text: async () => response });
    await fixture.context.runPhytClust();
    assert.equal(fixture.state.lastRunData, old);
    assert.equal(fixture.state.lastRunId, 'old');
    assert.equal(fixture.state.runHistory.length, 0);
    assert.equal(fixture.button.disabled, false);
  });
}

test('fixed-k run clears previous score data', async () => {
  const fixture = setup();
  fixture.state.lastScoresData = { scores: [2, 4] };
  await fixture.context.runPhytClust();
  assert.equal(fixture.state.lastScoresData, null);
  assert.equal(fixture.state.lastRunId, 'new');
  assert.equal(fixture.plot.children[0].textContent, 'Scores are not calculated for a fixed cluster count.');
  assert.equal(fixture.marker.hidden, true);
});

test('disabled settings are ignored and fractional counts are rejected', () => {
  const fixture = setup();
  fixture.element('top-n', 'invalid');
  fixture.element('support-weight', 'invalid');
  fixture.element('prominence-weight', 'invalid');
  assert.equal(fixture.context.buildRunRequest().k, 2);
  fixture.mode.value = 'global';
  fixture.element('top-n', '1.5');
  assert.throws(() => fixture.context.buildRunRequest(), /Peak count must be a positive whole number/);
  fixture.element('top-n', '3');
  fixture.element('ranking-mode', 'raw');
  const request = fixture.context.buildRunRequest();
  assert.equal(request.prominence_weight, undefined);
  assert.equal(request.support_weight, undefined);
});

test('invalid active settings do not submit a request', async () => {
  const fixture = setup();
  fixture.element('k', '1.5');
  await fixture.context.runPhytClust();
  assert.equal(fixture.calls.requests, 0);
  assert.match(fixture.messages[0].message, /Cluster count/);
});

test('outgroup validation is left to the server', async () => {
  const fixture = setup();
  fixture.elements.get('outgroup').value = 'internal clade';
  await fixture.context.runPhytClust();
  assert.equal(fixture.calls.requests, 1);
  assert.equal(fixture.state.lastRunId, 'new');
});

test('layout and tree failures restore result state', async () => {
  for (const failure of ['failLayouts', 'failTree', 'failParse']) {
    const fixture = setup();
    const old = previousResult(fixture);
    fixture.flags[failure] = true;
    await fixture.context.runPhytClust();
    assert.equal(fixture.state.lastRunData, old);
    assert.equal(fixture.state.currentClusters, old.clusters[0]);
    assert.equal(fixture.state.runHistory.length, 0);
    assert.equal(fixture.state.isRunning, false);
  }
});

test('edits during a request leave the new result marked stale', async () => {
  const fixture = setup();
  fixture.context.fetch = async () => {
    fixture.elements.get('k').value = '3';
    return { ok: true, status: 200, text: async () => JSON.stringify(fixture.data) };
  };
  await fixture.context.runPhytClust();
  assert.equal(fixture.marker.hidden, false);
});

test('large cluster maps avoid spreading values into function arguments', async () => {
  const fixture = setup();
  fixture.data.clusters = [Object.fromEntries(Array.from({ length: 150000 }, (_, index) => [`leaf${index}`, index % 3]))];
  await fixture.context.runPhytClust();
  assert.equal(fixture.state.lastRunId, 'new');
  assert.equal(fixture.state.clusterColours.length, 3);
});

test('repeated clicks issue only one request', async () => {
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

test('no peaks is a successful empty selection and clears old scores', async () => {
  const fixture = setup();
  fixture.mode.value = 'global';
  fixture.data.mode = 'global';
  fixture.data.clusters = [];
  fixture.data.k_values = [];
  fixture.data.scores = [];
  await fixture.context.runPhytClust();
  assert.equal(fixture.state.lastScoresData, null);
  assert.equal(fixture.state.lastRunId, 'new');
  assert.ok(fixture.messages.some(({ message, type }) => message.startsWith('No partitions selected') && type === 'info'));
});

test('server notes become a warning notification and warning log entries', async () => {
  const fixture = setup();
  fixture.data.notes = ["Node name 'A' is duplicated; renaming to 'A_1'.", '1 of 6 branches have a negative length; they are set to 0.0.'];
  await fixture.context.runPhytClust();
  assert.equal(fixture.logEntries.map(({ type }) => type).join(), 'warning,warning');
  assert.equal(fixture.logEntries[0].message, fixture.data.notes[0]);
  const warning = fixture.messages.find(({ type }) => type === 'warning');
  assert.match(warning.message, /renaming to 'A_1'/);
  assert.match(warning.message, /\+1 more/);
});

test('runs without notes show no warning', async () => {
  const fixture = setup();
  fixture.data.notes = [];
  await fixture.context.runPhytClust();
  assert.equal(fixture.logEntries.length, 0);
  assert.ok(!fixture.messages.some(({ type }) => type === 'warning'));
});

test('a multi-resolution request answered by a global search is labelled as such', async () => {
  const fixture = setup();
  fixture.mode.value = 'resolution';
  fixture.data.mode = 'global';
  fixture.data.notes = [];
  await fixture.context.runPhytClust();
  assert.equal(fixture.state.runHistory[0].mode, 'global');
  assert.match(fixture.state.runHistory[0].label, /^resolution → global/);
  assert.match(fixture.logEntries[0].message, /global search was run instead/);
  // The library's own explanation is not repeated.
  const notes = fixture.context.runNotes({ mode: 'global', notes: ['Resolution mode needs at least 50 scored cluster counts.'] }, 'resolution');
  assert.equal(notes.length, 1);
  assert.equal(fixture.context.runNotes({ mode: 'resolution', notes: [' ', 3] }, 'resolution').length, 0);
});

test('clearing results removes every result view', async () => {
  const fixture = setup();
  await fixture.context.runPhytClust();
  const leafCount = fixture.element('leaf-count-label', '2 leaves');
  leafCount.textContent = '2 leaves';
  fixture.plot.children = [{}];
  fixture.state.compareConfigs = [{ id: 1 }];
  fixture.context.clearResults();
  for (const name of ['newickTree', 'cartesianHierarchy', 'circularHierarchy', 'lastRunData', 'lastRunId', 'lastRunSignature', 'lastScoresData']) {
    assert.equal(fixture.state[name], null, name);
  }
  assert.equal(JSON.stringify(fixture.state.currentClusters), '{}');
  assert.equal(JSON.stringify(fixture.state.clusterColours), '[]');
  assert.equal(JSON.stringify(fixture.state.compareConfigs), '[]');
  assert.equal(fixture.calls.selectorReset, 1);
  assert.equal(fixture.calls.miniScores.at(-1), null);
  assert.equal(fixture.calls.comparison, 1);
  assert.equal(fixture.plot.children.length, 0);
  assert.equal(leafCount.textContent, '');
  assert.equal(fixture.marker.hidden, true);
});
