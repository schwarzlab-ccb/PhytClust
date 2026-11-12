const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { test } = require('node:test');
const source = fs.readFileSync(path.join(__dirname, '../../src/phytclust/gui/static/js/ui/status.js'), 'utf8').replace(/^import .*;\n/gm, '').replace(/export /g, '');
function setup() {
  const toasts = [], events = [], classes = [];
  const statusEl = { textContent: '' };
  const context = vm.createContext({ statusEl, statusDot: { classList: { add: name => classes.push(name) } },
    showToast: (message, type) => toasts.push({ message, type }),
    window: { dispatchEvent: event => events.push(event.type) },
    Event: class { constructor(type) { this.type = type; } },
  });
  vm.runInContext(source, context);
  return { app: context, toasts, events, classes, statusEl };
}

test('warnings show a warning notification and warning status indicator', () => {
  const { app, toasts, classes, statusEl } = setup();
  app.showStatus('Some settings were skipped.', 'warning');
  assert.equal(statusEl.textContent, 'Some settings were skipped.');
  assert.deepEqual(classes, ['warning']);
  assert.equal(toasts[0].type, 'warning');
  assert.equal(app.getStatusLog()[0].type, 'warning');
});

test('errors retain error severity and notify log listeners', () => {
  const { app, toasts, classes, events } = setup();
  app.showStatus('Invalid session.', 'danger');
  assert.deepEqual(classes, ['error']);
  assert.equal(toasts[0].type, 'danger');
  assert.deepEqual(events, ['phytclust:log-update']);
});

test('log-only entries reach the activity log without a notification', () => {
  const { app, toasts, events, statusEl } = setup();
  statusEl.textContent = 'k = 3';
  app.logStatus("Node name 'A' is duplicated; renaming to 'A_1'.", 'warning');
  assert.equal(statusEl.textContent, 'k = 3');
  assert.equal(toasts.length, 0);
  assert.equal(app.getStatusLog()[0].type, 'warning');
  assert.deepEqual(events, ['phytclust:log-update']);
});

test('the activity log is bounded and can be cleared', () => {
  const { app } = setup();
  for (let index = 0; index < 55; index++) app.showStatus(String(index), 'info');
  assert.equal(app.getStatusLog().length, 50);
  assert.equal(app.getStatusLog()[0].message, '5');
  app.clearStatusLog();
  assert.equal(app.getStatusLog().length, 0);
});
