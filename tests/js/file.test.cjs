const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { test } = require('node:test');
const root = path.join(__dirname, '../../src/phytclust/gui/static/js');
const read = name => fs.readFileSync(path.join(root, name), 'utf8').replace(/^import .*;\n/gm, '').replace(/export /g, '');
function setup(text) {
  const statuses = [], toasts = [], events = [];
  const newick = { value: '', dispatchEvent(event) { events.push(event.type); } };
  const app = vm.createContext({ structuredClone, d3: { zoomIdentity: {} },
    Event: class { constructor(type) { this.type = type; } },
    document: { getElementById: () => null, querySelectorAll: () => [] }, newickEl: newick,
    showStatus: (message, type) => statuses.push({ message, type }),
    showToast: (message, type) => toasts.push({ message, type }),
    FileReader: class { readAsText() { this.onload({ target: { result: text } }); } },
  });
  vm.runInContext(read('state.js') + read('session.js') + read('ui/file.js'), app);
  return { app, statuses, toasts, newick, events };
}

test('valid JSON that is not a session is rejected without success status', () => {
  const { app, statuses } = setup('{}');
  app.loadFile({ name: 'session.json' });
  assert.equal(statuses.at(-1).type, 'danger');
  assert.equal(statuses.some(item => item.type === 'success'), false);
});

test('malformed JSON reports an error', () => {
  const { app, statuses } = setup('{bad');
  app.loadFile({ name: 'session.json' });
  assert.equal(statuses.at(-1).type, 'danger');
});

test('a valid session reports success', () => {
  const { app, statuses } = setup('{"version":2,"params":{}}');
  app.loadFile({ name: 'session.json' });
  assert.equal(statuses.at(-1).type, 'success');
});

test('partial restoration reports a warning instead of success', () => {
  const { app, statuses, toasts } = setup('{"version":2,"params":{"k":"bad"}}');
  app.loadFile({ name: 'session.json' });
  assert.equal(statuses.at(-1).type, 'warning');
  assert.equal(toasts[0].type, 'warning');
});

test('ordinary Newick files still load and notify listeners', () => {
  const { app, statuses, newick, events } = setup('(a:1,b:2);');
  app.loadFile({ name: 'tree.nwk' });
  assert.equal(newick.value, '(a:1,b:2);');
  assert.deepEqual(events, ['change']);
  assert.equal(statuses.at(-1).type, 'success');
});
