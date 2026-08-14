/* ============================================================
   PhytClust – utils.js
   Pure utility functions — no DOM, no state dependencies.
   ============================================================ */

export function escapeRegExp(str) {
  return str.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

export function isOutgroupInNewick(newick, outgroup) {
  if (!outgroup) return true;
  var name = escapeRegExp(outgroup);
  return new RegExp("[\\(,]\\s*'?" + name + "'?\\s*(?:[:,\\)])").test(newick);
}

export function estimateLeafCount(newick) {
  var inQuote = null,
    commas = 0;
  for (var i = 0; i < newick.length; i++) {
    var ch = newick[i];
    if (inQuote) {
      if (ch === inQuote) inQuote = null;
      continue;
    }
    if (ch === '"' || ch === "'") {
      inQuote = ch;
      continue;
    }
    if (ch === ",") commas++;
  }
  return Math.max(0, commas + 1);
}

export function resetExtraParams() {
  [
    "extra-k",
    "extra-outgroup",
    "extra-root-taxon",
    "extra-topn",
    "extra-bins",
    "extra-maxk",
    "extra-maxklimit",
    "extra-prominence-weight",
    "extra-min-cluster-size",
    "extra-min-prominence",
    "extra-outlier-threshold",
    "extra-min-support",
    "extra-support-weight",
  ].forEach(function (id) {
    var el = document.getElementById(id);
    if (el) el.value = "";
  });
  ["extra-outlier", "extra-optimize-polytomies"].forEach(function (id) {
    var el = document.getElementById(id);
    if (el) el.checked = true;
  });
  [
    "extra-resolution",
    "extra-compute-all",
    "extra-outlier-prefer-fewer",
    "extra-no-split-zero",
    "extra-use-support",
    "extra-relative-prom",
  ].forEach(function (id) {
    var el = document.getElementById(id);
    if (el) el.checked = false;
  });
  var rankMode = document.getElementById("extra-ranking-mode");
  if (rankMode) rankMode.value = "adjusted";
  var ratioMode = document.getElementById("extra-outlier-ratio-mode");
  if (ratioMode) ratioMode.value = "exp";
}

export function parseNewick(newick) {
  var tokens = newick.split(/\s*(;|\(|\)|,|:)\s*/);
  var stack = [];
  var current = {};
  for (var i = 0; i < tokens.length; i++) {
    var token = tokens[i];
    switch (token) {
      case "(":
        var child = {};
        current.children = [child];
        stack.push(current);
        current = child;
        break;
      case ",":
        var sibling = {};
        stack[stack.length - 1].children.push(sibling);
        current = sibling;
        break;
      case ")":
        current = stack.pop();
        break;
      case ":":
        var length = parseFloat(tokens[++i]);
        current.length = isNaN(length) ? 0 : length;
        break;
      default:
        if (token && token !== ";")
          current.name = token.replace(/^(['"])(.*)\1$/, "$2");
    }
  }
  return current;
}

export function getVisibleChildren(node) {
  if (node && node._collapsed) return null;
  return node && node.children && node.children.length ? node.children : null;
}

export function getAllChildren(node) {
  return node && node.children && node.children.length ? node.children : null;
}

export function readIntParam(id) {
  var el = document.getElementById(id);
  if (!el) return null;
  var raw = (el.value || "").trim();
  if (!raw) return null;
  var v = parseInt(raw, 10);
  return isNaN(v) ? null : v;
}

export function readFloatParam(id) {
  var el = document.getElementById(id);
  if (!el) return null;
  var raw = (el.value || "").trim();
  if (!raw) return null;
  var v = parseFloat(raw);
  return isNaN(v) ? null : v;
}

export function readCheckParam(id) {
  var el = document.getElementById(id);
  return el ? el.checked : false;
}

export function readSelectParam(id) {
  var el = document.getElementById(id);
  return el ? el.value : null;
}

export function escapeHtmlAttr(value) {
  return String(value || "")
    .replace(/&/g, "&amp;")
    .replace(/"/g, "&quot;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

export function collectLeafNamesData(node, out) {
  if (!node) return;
  const kids = node.children && node.children.length ? node.children : null;
  if (!kids) {
    if (node.name) out.push(node.name);
    return;
  }
  for (const c of kids) collectLeafNamesData(c, out);
}
