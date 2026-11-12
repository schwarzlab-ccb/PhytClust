const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { test } = require('node:test');
const root = path.join(__dirname, '../../src/phytclust/gui/static/js');
const read = (file) => fs.readFileSync(path.join(root, file), 'utf8').replace(/^import .*;\n/gm, '').replace(/export /g, '');

function setup(leaves) {
  const context = vm.createContext({
    structuredClone, d3: { zoomIdentity: {} },
    document: { getElementById: () => null },
    showToast() {}, escapeHtmlAttr: String, drawComparisonBarsInto() {},
    generateClusterColours: () => [], getThemeColours: () => ({}), withAlpha: (colour, alpha) => `${colour}@${alpha}`,
  });
  vm.runInContext(read('state.js') + read('utils.js') + read('views/compare.js') + '\nglobalThis.state = state;', context);
  context.state.newickTree = { children: leaves.map((name) => ({ name })) };
  return context;
}

const LEAVES = ['a', 'b', 'c', 'd', 'e', 'f'];

test('leaves missing from a comparison file stay unassigned and are counted', () => {
  const app = setup(LEAVES);
  const parsed = app.parseClusterTableText('leaf,cluster\na,5\nb,5\n');
  const matched = app.matchClustersToLeaves(parsed.map);
  assert.equal(JSON.stringify(matched.clusters), JSON.stringify({ a: 5, b: 5 }));
  assert.equal(matched.missing, 4);
  assert.equal(matched.total, 6);
  assert.equal(JSON.stringify(app.normalizeClustersForLeaves(parsed.map)), JSON.stringify({ a: 5, b: 5 }));
  assert.equal(app.matchClustersToLeaves({ x: 1 }), null);
});

test('the agreement matrix uses only leaves assigned in both maps', () => {
  const app = setup(LEAVES);
  const run = { a: 0, b: 0, c: 1, d: 1, e: 2, f: 2 };
  const file = app.normalizeClustersForLeaves({ a: 5, b: 5 });
  const table = app.buildContingency(run, file);
  assert.equal(table.rowIds.join(), '0');
  assert.equal(table.colIds.join(), '5');
  assert.equal(table.cells[0][5], 2);
  const reverse = app.buildContingency(file, run);
  assert.equal(reverse.rowIds.join(), '5');
  assert.equal(reverse.colIds.join(), '0');
});

test('the comparison TSV leaves unassigned leaves blank', () => {
  const app = setup(LEAVES);
  const tsv = app.buildComparisonTSV([{ title: 'file', clusters: { a: 5, b: 5 } }]);
  assert.equal(tsv.split('\n')[3], 'c\t');
});

test('the current comparisons resolve k bars from the run and keep file bars as loaded', () => {
  const app = setup(LEAVES);
  app.state.lastRunData = { all_ks: [2, 3], all_clusters: [{ a: 0, b: 0, c: 1, d: 1, e: 1, f: 1 }, { a: 0, b: 0, c: 1, d: 1, e: 2, f: 2 }] };
  app.state.compareConfigs = [
    { id: 1, title: 'k=3', source: 'k', k: 3 },
    { id: 2, title: '', source: 'k', k: 2 },
    { id: 3, title: '', source: 'file', fileName: 'mine.csv', fileClusters: { a: 5, b: 5 } },
    { id: 4, title: 'missing', source: 'k', k: 9 },
    { id: 5, title: 'empty', source: 'file', fileClusters: null },
  ];
  const comparisons = app.buildCurrentComparisons();
  assert.equal(comparisons.map((c) => c.title).join(), 'k=3,k=2,mine.csv');
  assert.equal(comparisons[0].clusters.e, 2);
  assert.equal(comparisons[1].clusters.e, 1);
  assert.deepEqual(JSON.parse(JSON.stringify(comparisons[2].clusters)), { a: 5, b: 5 });
  assert.equal(app.buildComparisonTSV(comparisons).split('\n')[0], 'leaf\tk=3\tk=2\tmine.csv');
});
