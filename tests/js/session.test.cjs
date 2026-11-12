const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { test } = require('node:test');
const root = path.join(__dirname, '../../src/phytclust/gui/static/js');
const read = (file) => fs.readFileSync(path.join(root, file), 'utf8').replace(/^import .*;\n/gm, '').replace(/export /g, '');
function setup() {
  const inputs = new Map();
  const documentEvents = [];
  const buttons = ['global', 'k', 'resolution'].map(mode => ({ dataset: { mode }, click() { active = this; } }));
  let active = buttons[0];
  const context = vm.createContext({ structuredClone, TextEncoder, d3: { zoomIdentity: {} }, Event: class { constructor(type) { this.type = type; } }, document: {
    getElementById: id => inputs.get(id) || null,
    querySelector: () => active,
    querySelectorAll: () => buttons,
    dispatchEvent: event => documentEvents.push(event.type),
  } });
  vm.runInContext(read('state.js') + read('session.js') + '\nglobalThis.state = state;', context);
  function input(id, value = '') {
    const control = { value, checked: false, options: [], events: [], dispatchEvent(event) { this.events.push(event.type); } };
    inputs.set(id, control);
    return control;
  }
  return { app: context, input, documentEvents };
}

test('numeric collection rejects partial numbers and nonfinite values', () => {
  const { app, input } = setup();
  const k = input('k', '3oops');
  input('max-k-limit', 'Infinity');
  assert.equal(app.collectSession().params.k, null);
  assert.equal(app.collectSession().params.max_k_limit, null);
  k.value = '3';
  assert.equal(app.collectSession().params.k, 3);
});

test('invalid documents do not mutate controls', () => {
  const { app, input } = setup();
  const k = input('k', '7');
  assert.equal(app.applySession({ params: { k: 3 } }).applied, 0);
  assert.equal(k.value, '7');
});

test('invalid checkbox, numeric, mode, and select values are skipped', () => {
  const { app, input } = setup();
  const k = input('k', '7');
  const checkbox = input('use-support');
  const select = input('ranking-mode', 'raw');
  select.options = [{ value: 'raw' }, { value: 'adjusted' }];
  const result = app.applySession({ version: 2, params: { mode: '"]broken', k: '3', support: { use: 'false' }, peak: { ranking_mode: 'unknown' } } });
  assert.equal(result.errors.length, 4);
  assert.equal(result.applied, 0);
  assert.equal(k.value, '7');
  assert.equal(checkbox.checked, false);
  assert.equal(select.value, 'raw');
});

test('valid settings restore and notify listeners; missing controls are not counted', () => {
  const { app, input } = setup();
  const k = input('k');
  const tree = input('newick-input');
  const result = app.applySession({ version: 1, params: { mode: 'k', k: 3, top_n: 2 }, newick: ' (a,b); ' });
  assert.equal(result.applied, 2);
  assert.equal(result.hadNewick, true);
  assert.equal(k.value, '3');
  assert.deepEqual(k.events, ['input', 'change']);
  assert.equal(tree.value, '(a,b);');
});

test('appearance failures are atomic and cannot alter inherited properties', () => {
  const { app } = setup();
  const original = app.state.render;
  for (const appearance of [{ labels: { fontSize: 12, show: 'false' } }, { labels: null }, { layout: 'bad' }, JSON.parse('{"__proto__":{"polluted":true}}')]) {
    assert.equal(app.applySession({ version: 2, params: {}, appearance }).errors.length, 1);
    assert.equal(app.state.render, original);
    assert.equal(app.state.render.labels.fontSize, 9);
  }
  assert.equal(app.applySession({ version: 2, params: {}, appearance: { labels: { fontSize: 12 } } }).errors.length, 0);
  assert.equal(app.state.render.labels.fontSize, 12);
});

test('saved selection is an actual cluster count, preserved until the next run', () => {
  const { app, input } = setup();
  input('cluster-select', '1');
  app.state.lastRunData = { k_values: [3, 8] };
  assert.equal(app.collectSession().view.selected_k, 8);
  app.applySession({ version: 2, params: {}, view: { selected_k: 8, cluster_mode: 'all' } });
  assert.equal(app.state.pendingSessionView.selected_k, 8);
  app.state.lastRunData = null;
  assert.equal(app.collectSession().view.selected_k, 8);
  // Version 1 stored an index under selected_k, so it cannot safely be reused.
  app.applySession({ version: 1, params: {}, view: { selected_k: 1 } });
  assert.equal(app.state.pendingSessionView, null);
});

test('size warning uses UTF-8 byte length', () => {
  const { app } = setup();
  assert.equal(app.sessionSizeWarning({}), null);
  assert.equal(app.sessionSizeWarning({ newick: 'a' }), null);
  assert.match(app.sessionSizeWarning({ newick: 'é'.repeat(300000) }), /586 KB/);
});

test('fresh results restore a saved cluster count from the all-partitions view', () => {
  const { app, input } = setup();
  const select = input('cluster-select');
  select.appendChild = () => {};
  const controls = input('cluster-select-controls');
  controls.classList = { add() {}, remove() {} };
  app.document.createElement = () => ({});
  app.generateClusterColours = () => ['red'];
  app.drawTree = () => {};
  vm.runInContext(read('utils.js') + read('views/cluster_editor.js') + '\nupdateClusterEditorAvailability = () => {};', app);
  app.applySession({ version: 2, params: {}, view: { selected_k: 4, cluster_mode: 'all' } });
  const data = { clusters: [{ a: 0 }], k_values: [1], all_clusters: [{ a: 0 }, { a: 3 }], all_ks: [1, 4] };
  app.state.lastRunData = data;
  app.populateClusterSelector(data);
  assert.equal(select.value, '1');
  assert.equal(app.state.runClusters, data.all_clusters[1]);
  assert.equal(app.state.pendingSessionView, null);
});

test('saving an imported selection before rerunning preserves it over old results', () => {
  const { app, input } = setup();
  input('cluster-select', '0');
  app.state.lastRunData = { k_values: [3, 8] };
  app.applySession({ version: 2, params: {}, view: { selected_k: 8, cluster_mode: 'all' } });
  assert.equal(app.collectSession().view.selected_k, 8);
  assert.equal(app.collectSession().view.cluster_mode, 'all');
});

test('appearance import synchronizes toolbar controls and requests a redraw', () => {
  const { app, input, documentEvents } = setup();
  const layout = input('layout-mode', 'rectangular');
  const colourMode = input('colour-mode', 'bars');
  const result = app.applySession({ version: 2, params: {}, appearance: { layout: 'circular', clusters: { colorMode: 'boxes' } } });
  assert.equal(result.valid, true);
  assert.equal(layout.value, 'circular');
  assert.equal(colourMode.value, 'boxes');
  assert.equal(app.state.render.layout, layout.value);
  assert.deepEqual(documentEvents, ['phytclust:appearance-restored']);
});

test('invalid documents are distinguished from valid sessions with skipped fields', () => {
  const { app, input } = setup();
  input('k', '7');
  assert.equal(app.applySession({}).valid, false);
  const partial = app.applySession({ version: 2, params: { k: 'invalid' } });
  assert.equal(partial.valid, true);
  assert.equal(partial.errors.length, 1);
});

test('every toolbar colour mode survives a session round trip', () => {
  for (const mode of ['bars', 'boxes', 'nodes']) {
    const saved = setup();
    saved.app.state.render.clusters.colorMode = mode;
    saved.app.state.render.layout = 'circular';
    const session = JSON.parse(JSON.stringify(saved.app.collectSession({ includeNewick: false })));
    const restored = setup();
    const colourMode = restored.input('colour-mode', 'bars');
    // This fixture has no sidebar controls, so restore only the saved appearance.
    const result = restored.app.applySession({ version: session.version, params: {}, appearance: session.appearance });
    assert.equal(result.errors.length, 0, mode + ': ' + result.errors.join('; '));
    assert.equal(restored.app.state.render.clusters.colorMode, mode);
    assert.equal(restored.app.state.render.layout, 'circular');
    assert.equal(colourMode.value, mode);
  }
  const { app } = setup();
  assert.equal(app.applySession({ version: 2, params: {}, appearance: { clusters: { colorMode: 'stripes' } } }).errors.length, 1);
});

test('branch colour round-trips between theme default and a custom colour', () => {
  const { app } = setup();
  const version = 2;

  app.state.render.branches.color = '#aabbcc';
  const themed = app.applySession({ version, params: {}, appearance: { branches: { color: null } } });
  assert.equal(themed.errors.length, 0, themed.errors.join('; '));
  assert.equal(app.state.render.branches.color, null);

  const custom = app.applySession({ version, params: {}, appearance: { branches: { color: '#123456' } } });
  assert.equal(custom.errors.length, 0, custom.errors.join('; '));
  assert.equal(app.state.render.branches.color, '#123456');

  const rejected = app.applySession({ version, params: {}, appearance: { branches: { color: 'not-a-colour' } } });
  assert.equal(rejected.errors.length, 1);
  assert.ok(rejected.errors.join(' ').includes('branches.color'));
});
