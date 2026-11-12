const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { test } = require('node:test');
const root = path.join(__dirname, '../../src/phytclust/gui/static/js');
const read = (file) => fs.readFileSync(path.join(root, file), 'utf8').replace(/^import .*;\n/gm, '').replace(/export /g, '');

function element(extra = {}) {
  return {
    value: '', checked: false, textContent: '', innerHTML: '', style: {}, children: [],
    classList: { removed: [], remove(name) { this.removed.push(name); }, add() {} },
    appendChild(child) { this.children.push(child); },
    setAttribute(name, value) { this[name] = value; },
    getAttribute(name) { return this[name]; },
    ...extra,
  };
}

function setup() {
  const elements = new Map();
  const context = vm.createContext({
    structuredClone, d3: { zoomIdentity: {} },
    document: {
      getElementById: (id) => elements.get(id) || null,
      createElement: () => element(),
    },
    generateClusterColours: (count) => Array.from({ length: count }, (_, index) => 'c' + index),
    drawTree() {}, showToast() {},
  });
  vm.runInContext(read('state.js') + read('utils.js') + read('views/cluster_editor.js') + read('ui/export.js') + '\nglobalThis.state = state;', context);
  return { app: context, elements, add: (id, extra) => { const node = element(extra); elements.set(id, node); return node; } };
}

const RUN = { a: 0, b: 0, c: 1, d: 2, e: 2, f: 2, g: 3 };

test('small clusters are outliers: singletons by default, else below the threshold', () => {
  const { app } = setup();
  assert.equal(JSON.stringify(app.markOutlierClusters(RUN)), JSON.stringify({ a: 0, b: 0, c: -1, d: 2, e: 2, f: 2, g: -1 }));
  assert.equal(JSON.stringify(app.markOutlierClusters(RUN, { threshold: 3 })), JSON.stringify({ a: -1, b: -1, c: -1, d: 2, e: 2, f: 2, g: -1 }));
  assert.equal(JSON.stringify(app.markOutlierClusters(RUN, { enabled: false })), JSON.stringify(RUN));
  // Ids that are already negative (uploaded files) are left alone.
  assert.equal(app.markOutlierClusters({ x: -1, y: -1, z: 4 }).x, -1);
});

test('the drawn map and the exported TSV mark the same leaves', () => {
  for (const [threshold, checked] of [[null, true], [3, true], [3, false]]) {
    const { app, add } = setup();
    add('mark-outliers', { checked });
    app.state.lastRunData = { outlier_size_threshold: threshold };
    app.setCurrentClusters(RUN);
    const rows = app.currentTSV().text.trim().split('\n').slice(1).map((row) => row.split('\t'));
    assert.equal(rows.length, 7);
    for (const [name, clusterId] of rows) assert.equal(Number(clusterId), app.state.currentClusters[name], `${name} at ${threshold}/${checked}`);
    assert.equal(app.currentTSV().k, 4);
  }
});

test('the run threshold is used, not an edited sidebar value', () => {
  const { app, add } = setup();
  add('mark-outliers', { checked: true });
  add('outlier-threshold', { value: '5' });
  app.state.lastRunData = { outlier_size_threshold: null };
  app.setCurrentClusters(RUN);
  assert.equal(app.state.currentClusters.a, 0);
});

test('switching "Mark outliers" off restores run ids without changing colours', () => {
  const { app, add } = setup();
  const mark = add('mark-outliers', { checked: true });
  app.state.lastRunData = { outlier_size_threshold: null };
  app.setCurrentClusters(RUN);
  const colours = app.state.clusterColours;
  assert.equal(colours.length, 4);
  assert.equal(app.state.currentClusters.c, -1);
  mark.checked = false;
  app.refreshOutlierMarks();
  assert.equal(app.state.currentClusters.c, 1);
  assert.equal(app.state.clusterColours, colours);
  assert.equal(app.state.runClusters, RUN);
});

test('the cluster editor lists run outliers and "Include outliers" hides them', () => {
  const { app, add } = setup();
  add('mark-outliers', { checked: true });
  const body = add('cluster-editor-body');
  const include = add('cluster-filter-show-outliers', { checked: true });
  app.state.lastRunData = { outlier_size_threshold: null };
  app.setCurrentClusters(RUN);
  assert.equal(app.getCurrentClusterIds().join(), '-1,0,2');
  app.renderClusterEditorRows();
  assert.equal(body.children.length, 3);
  include.checked = false;
  body.children = [];
  app.renderClusterEditorRows();
  assert.equal(body.children.length, 2);
});

test('large maps do not spread values into function arguments', () => {
  const { app } = setup();
  const map = Object.fromEntries(Array.from({ length: 150000 }, (_, index) => [`leaf${index}`, index % 3]));
  app.setCurrentClusters(map);
  assert.equal(app.state.clusterColours.length, 3);
});

test('resetting the selector hides the result controls', () => {
  const { app, add } = setup();
  const controls = add('cluster-select-controls');
  const select = add('cluster-select', { innerHTML: '<option>k = 3</option>' });
  const toggle = add('cluster-view-toggle');
  const edit = add('btn-edit-clusters');
  toggle.style.display = 'inline-flex';
  app.state.clusterViewMode = 'all';
  app.setCurrentClusters({});
  app.resetClusterSelector();
  assert.deepEqual(controls.classList.removed, ['visible']);
  assert.equal(select.innerHTML, '');
  assert.equal(toggle.style.display, 'none');
  assert.equal(edit.style.display, 'none');
  assert.equal(app.state.clusterViewMode, 'peaks');
});
