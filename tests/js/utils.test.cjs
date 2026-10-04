const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const { test } = require("node:test");

const controls = new Map();
const context = vm.createContext({ document: { getElementById: (id) => controls.get(id) || null } });
const source = fs.readFileSync(path.join(__dirname, "../../src/phytclust/gui/static/js/utils.js"), "utf8").replace(/export /g, "");
vm.runInContext(source, context);
const plain = (value) => JSON.parse(JSON.stringify(value));

test("ordinary trees keep names, lengths, supports, and child order", () => {
  const tree = context.parseNewick("((a:1,b:2)95:3,c:4)Root:0;");
  assert.deepEqual(plain(tree), { name: "Root", length: 0, children: [{ name: "95", length: 3, children: [{ name: "a", length: 1 }, { name: "b", length: 2 }] }, { name: "c", length: 4 }] });
  assert.deepEqual(plain(context.collectLeafNamesData(tree)), ["a", "b", "c"]);
});

test("quoted punctuation and doubled quotes stay inside labels", () => {
  const text = "('A,(B): C':1,'O''Brien':2,\"D,E\":3);";
  assert.deepEqual(plain(context.collectLeafNamesData(context.parseNewick(text))), ["A,(B): C", "O'Brien", "D,E"]);
  assert.equal(context.estimateLeafCount(text), 3);
});

test("comments do not become labels, lengths, or extra leaves", () => {
  const text = "[&R]('a':1[notes, here],b[outer[inner, comment]]:2)95[comment];";
  const tree = context.parseNewick(text);
  assert.deepEqual(plain(context.collectLeafNamesData(tree)), ["a", "b"]);
  assert.equal(tree.children[0].length, 1);
  assert.equal(tree.name, "95");
  assert.equal(context.estimateLeafCount(text), 2);
});

test("single leaves and optional semicolons are accepted", () => {
  assert.equal(context.parseNewick("a:1e-3").length, 0.001);
  assert.equal(context.estimateLeafCount("a:1;"), 1);
  assert.equal(context.parseNewick("(a:.5,b:+2E1)").children[1].length, 20);
});

for (const text of ["", ";", "(a,b;", "a,b;", "a);", "(a:1:2,b);", "(a:1x,b);", "(a:Infinity,b);", "(a:,b);", "a:", "('a,b);", "(a[comment,b);", "(a,b);(c,d);", "a'other';"]) {
  test(`invalid tree reports an error: ${text}`, () => {
    assert.throws(() => context.parseNewick(text), /Newick|branch length|tree|name/i);
  });
}

test("empty input has zero estimated leaves", () => {
  assert.equal(context.estimateLeafCount(""), 0);
  assert.equal(context.estimateLeafCount("  ; [comment] "), 0);
  assert.equal(context.estimateLeafCount("('unclosed"), 0);
});

test("deep trees parse and collect leaves without recursion", () => {
  const depth = 25000;
  const tree = context.parseNewick("(".repeat(depth) + "a:1" + ")".repeat(depth) + ";");
  const names = ["existing"];
  assert.equal(context.collectLeafNamesData(tree, names), names);
  assert.deepEqual(names, ["existing", "a"]);
});

test("collapsed nodes still expose all leaves when requested", () => {
  const tree = context.parseNewick("(a,b);");
  tree._collapsed = true;
  assert.equal(context.getVisibleChildren(tree), null);
  assert.equal(context.getAllChildren(tree), tree.children);
  assert.deepEqual(plain(context.collectLeafNamesData(tree)), ["a", "b"]);
});

test("HTML attributes retain zero and escape both quote styles", () => {
  assert.equal(context.escapeHtmlAttr(0), "0");
  assert.equal(context.escapeHtmlAttr(false), "false");
  assert.equal(context.escapeHtmlAttr(null), "");
  assert.equal(context.escapeHtmlAttr(`A&\"'<B>`), "A&amp;&quot;&#39;&lt;B&gt;");
});

test("reset restores the supported default options", () => {
  for (const id of ["extra-k", "extra-outlier", "extra-resolution", "extra-exclude-k2", "extra-ranking-mode", "extra-outlier-ratio-mode", "extra-polytomy-mode"]) controls.set(id, { value: "changed", checked: true });
  context.resetExtraParams();
  assert.equal(controls.get("extra-k").value, "");
  assert.equal(controls.get("extra-outlier").checked, true);
  assert.equal(controls.get("extra-exclude-k2").checked, false);
  assert.equal(controls.get("extra-ranking-mode").value, "adjusted");
  assert.equal(controls.get("extra-outlier-ratio-mode").value, "exp");
  assert.equal(controls.get("extra-polytomy-mode").value, "soft");
});
