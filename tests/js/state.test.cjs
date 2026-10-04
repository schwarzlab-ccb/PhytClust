const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { test } = require('node:test');
const source = fs.readFileSync(path.join(__dirname, '../../src/phytclust/gui/static/js/state.js'), 'utf8').replace(/export /g, '');
function loadState() {
  const context = vm.createContext({ d3: { zoomIdentity: {} }, structuredClone });
  vm.runInContext(source + '\nglobalThis.state = state; globalThis.selection = SELECTED_CLUSTER_IDS;', context);
  return context;
}

test('display controls read and update existing settings', () => {
  const app = loadState();
  assert.equal(app.getRenderOption('labels.fontSize'), 9);
  assert.equal(app.setRenderOption('labels.fontSize', 12), 12);
  assert.equal(app.state.render.labels.fontSize, 12);
});

test('invalid paths cannot change inherited properties', () => {
  const app = loadState();
  for (const name of ['', 'labels.missing', 'missing.size', '__proto__.polluted', 'labels.constructor.prototype.polluted']) {
    assert.throws(() => app.setRenderOption(name, true), /display setting/);
    assert.throws(() => app.getRenderOption(name), /display setting/);
  }
});

test('partial overrides leave the original settings untouched', () => {
  const app = loadState();
  const copy = app.renderOptionsWithOverrides({ labels: { fontSize: 15 }, branches: { color: '#fff' } });
  assert.equal(copy.labels.fontSize, 15);
  assert.equal(copy.labels.show, true);
  assert.equal(app.state.render.labels.fontSize, 9);
  copy.nodes.leafRadius = 100;
  assert.equal(app.state.render.nodes.leafRadius, 3);
});

test('malformed overrides fail before changing live settings', () => {
  const app = loadState();
  const original = app.state.render;
  for (const overrides of [null, [], { labels: null }, { labels: 2 }, { missing: 2 }, { labels: { fontSize: {} } }, JSON.parse('{"__proto__":{"polluted":true}}')]) {
    assert.throws(() => app.renderOptionsWithOverrides(overrides), /Display|display/);
    assert.equal(app.state.render, original);
  }
});

test('export takes a snapshot with overrides and redraws the original settings', () => {
  const app = loadState();
  const original = app.state.render;
  const redrawSizes = [];
  const snapshot = app.withRenderOverrides({ labels: { fontSize: 20 } }, () => redrawSizes.push(app.state.render.labels.fontSize), () => app.state.render.labels.fontSize);
  assert.equal(snapshot, 20);
  assert.deepEqual(redrawSizes, [20, 9]);
  assert.equal(app.state.render, original);
});

test('failed snapshots and initial redraws restore the original settings', () => {
  for (const failDuringRedraw of [false, true]) {
    const app = loadState();
    const original = app.state.render;
    let redrawCount = 0;
    assert.throws(() => app.withRenderOverrides({}, () => {
      if (++redrawCount === 1 && failDuringRedraw) throw new Error('redraw failed');
    }, () => { throw new Error('snapshot failed'); }), /failed/);
    assert.equal(app.state.render, original);
    assert.equal(redrawCount, 2);
  }
});

test('box adjustments are independent even for inherited property names', () => {
  const app = loadState();
  const first = app.getBoxAdjust('__proto__');
  assert.equal(first.dx, 0);
  first.dx = 5;
  assert.equal(app.getBoxAdjust('__proto__').dx, 5);
  assert.equal(app.getBoxAdjust('constructor').dx, 0);
  assert.equal(app.getBoxAdjust(1), app.getBoxAdjust('1'));
});

test('cluster focus accepts numeric IDs from display controls', () => {
  const app = loadState();
  assert.equal(app.hasClusterFocus(null), true);
  assert.equal(app.hasClusterFocus(2), true);
  app.selection.add(2);
  assert.equal(app.hasClusterFocus('2'), true);
  assert.equal(app.hasClusterFocus(3), false);
  assert.equal(app.hasClusterFocus(null), false);
});
