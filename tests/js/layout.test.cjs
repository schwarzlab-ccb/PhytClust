const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { test } = require('node:test');
const source = fs.readFileSync(path.join(__dirname, '../../src/phytclust/gui/static/js/tree/layout.js'), 'utf8').replace(/^import .*;\n/gm, '').replace(/export /g, '');
const app = vm.createContext({});
vm.runInContext(source, app);

test('branch accumulation preserves distances and treats missing lengths as zero', () => {
  const tree = { length: 100, children: [{ length: 3, children: [{ length: 4 }, {}] }, { length: 2 }] };
  app.accumulateBranchLength(tree, 5);
  assert.equal(tree._bl, 5);
  assert.equal(tree.children[0]._bl, 8);
  assert.equal(tree.children[0].children[0]._bl, 12);
  assert.equal(tree.children[0].children[1]._bl, 8);
  assert.equal(tree.children[1]._bl, 7);
});

test('a deep tree accumulates every edge without overflowing the call stack', () => {
  const tree = {};
  let last = tree;
  for (let index = 0; index < 25000; index++) {
    const child = { length: 1 };
    last.children = [child];
    last = child;
  }
  app.accumulateBranchLength(tree);
  assert.equal(tree._bl, 0);
  assert.equal(last._bl, 25000);
});
