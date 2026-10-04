const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { test } = require('node:test');
const d3 = require('../../src/phytclust/gui/static/vendor/d3.v7.min.js');
const source = fs.readFileSync(path.join(__dirname, '../../src/phytclust/gui/static/js/colors.js'), 'utf8').replace(/^import .*;\n/gm, '').replace(/export /g, '');
function setup(storage = { getItem: () => null, setItem() {} }) {
  const attributes = new Map();
  const warnings = [];
  const context = vm.createContext({ d3, state: { render: { branches: { color: null } } }, console: { warn: message => warnings.push(message) }, localStorage: storage,
    document: { documentElement: { setAttribute: (key, value) => attributes.set(key, value), getAttribute: key => attributes.get(key) } },
    getComputedStyle: () => ({ getPropertyValue: () => '' }),
  });
  vm.runInContext(source + '\nglobalThis.baseColors = BASE_COLORS;', context);
  return { app: context, warnings, attributes };
}

test('palette order and expanded shades match Python defaults', () => {
  const { app } = setup();
  const colors = app.generateClusterColors(33);
  assert.equal(colors.length, 33);
  assert.deepEqual(Array.from(colors.slice(0, 16)), Array.from(app.baseColors));
  const shade = d3.color(colors[16]);
  const red = d3.color(app.baseColors[0]);
  assert.equal(shade.r, Math.round(red.r + (255 - red.r) * 0.14));
  assert.equal(shade.opacity, 0.88);
  assert.equal(d3.color(colors[32]).opacity, 0.76);
});

test('opacity supports hex, short hex, RGB, RGBA, HSL, and named colours', () => {
  const { app } = setup();
  for (const color of ['#ff0000', '#f00', 'rgb(255, 0, 0)', 'rgba(255, 0, 0, 0.5)', 'hsl(0, 100%, 50%)', 'red']) {
    assert.equal(app.withAlpha(color, 0.25), 'rgba(255, 0, 0, 0.25)');
  }
  assert.equal(app.withAlpha('transparent', 0), 'rgba(0, 0, 0, 0)');
});

test('lightening retains an existing opacity', () => {
  const { app } = setup();
  const color = d3.color(app.adjustLight('rgba(0, 0, 0, 0.4)', 0.5));
  assert.equal(color.r, 128);
  assert.equal(color.opacity, 0.4);
  assert.equal(d3.color(app.adjustLight('#fff', 0)).r, 255);
});

test('invalid colours and fractions fail clearly', () => {
  const { app } = setup();
  assert.throws(() => app.withAlpha('bad', 0.5), /Invalid colour/);
  for (const value of [-1, 2, NaN, Infinity, '0.5']) {
    assert.throws(() => app.withAlpha('#fff', value), /Opacity/);
    assert.throws(() => app.adjustLight('#fff', value), /Lightening fraction/);
  }
});

test('empty and malformed palette counts fail without looping', () => {
  const { app } = setup();
  for (const count of [-1, 1.5, Infinity, NaN, Number.MAX_SAFE_INTEGER + 1]) assert.throws(() => app.generateClusterColors(count), /Colour count/);
  app.baseColors.length = 0;
  assert.equal(app.generateClusterColors(0).length, 0);
  assert.throws(() => app.generateClusterColors(1), /palette is empty/);
});

test('repeat warning begins after seven palette rounds', () => {
  const { app, warnings } = setup();
  app.generateClusterColors(112);
  assert.equal(warnings.length, 0);
  const colors = app.generateClusterColors(129);
  assert.equal(warnings.length, 1);
  assert.equal(colors[112], colors[128]);
});

test('expanded custom colours retain their original opacity', () => {
  const { app } = setup();
  app.baseColors.splice(0, app.baseColors.length, 'rgba(10, 20, 30, 0.5)');
  assert.equal(d3.color(app.generateClusterColors(2)[1]).opacity, 0.44);
});

test('shuffling leaves the source palette unchanged', () => {
  const { app } = setup();
  const original = ['a', 'b', 'c'];
  const shuffled = app.shuffle(original);
  assert.deepEqual(original, ['a', 'b', 'c']);
  assert.deepEqual(Array.from(shuffled).sort(), original);
  assert.notEqual(shuffled, original);
});

test('theme starts light and honors a saved selection', () => {
  const { app, attributes } = setup();
  app.initTheme();
  assert.equal(attributes.get('data-theme'), 'light');
  const saved = setup({ getItem: () => 'dark', setItem() {} });
  saved.app.initTheme();
  assert.equal(saved.attributes.get('data-theme'), 'dark');
});

test('theme can load and toggle when browser storage is unavailable', () => {
  const { app, attributes, warnings } = setup({ getItem() { throw Error('blocked'); }, setItem() { throw Error('blocked'); } });
  app.initTheme();
  assert.equal(attributes.get('data-theme'), 'light');
  assert.equal(app.toggleTheme(), 'dark');
  assert.equal(attributes.get('data-theme'), 'dark');
  assert.equal(warnings.length, 2);
});

test('theme colours keep defaults and respect custom branch colour', () => {
  const { app } = setup();
  assert.equal(app.getThemeColors().branch, '#000000');
  app.state.render.branches.color = '#f00';
  assert.equal(app.getThemeColors().branch, '#f00');
});
