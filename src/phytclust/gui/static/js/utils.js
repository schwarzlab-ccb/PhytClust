/* Parse trees and read or reset GUI controls. */

export function resetParams() {
  [
    "k",
    "outgroup",
    "root-taxon",
    "top-n",
    "bins",
    "max-k",
    "max-k-limit",
    "prominence-weight",
    "min-cluster-size",
    "min-prominence",
    "outlier-threshold",
    "min-support",
    "support-weight",
  ].forEach(function (id) {
    const element = document.getElementById(id);
    if (element) element.value = "";
  });
  ["mark-outliers"].forEach(function (id) {
    const element = document.getElementById(id);
    if (element) element.checked = true;
  });
  [
    "compute-all",
    "outlier-prefer-fewer",
    "no-split-zero",
    "use-support",
    "relative-prominence",
    "exclude-k2",
  ].forEach(function (id) {
    const element = document.getElementById(id);
    if (element) element.checked = false;
  });
  const rankingMode = document.getElementById("ranking-mode");
  if (rankingMode) rankingMode.value = "adjusted";
  const penaltyMode = document.getElementById("outlier-ratio-mode");
  if (penaltyMode) penaltyMode.value = "exp";
  const polytomyMode = document.getElementById("polytomy-mode");
  if (polytomyMode) polytomyMode.value = "soft";
}

// Read labels without splitting punctuation inside quotes or comments.
function* newickTokens(text) {
  let position = 0;
  while (position < text.length) {
    const character = text[position];
    if (/\s/.test(character)) {
      position++;
    } else if (character === "[") {
      let depth = 1;
      position++;
      while (position < text.length && depth) {
        if (text[position] === "[") depth++;
        if (text[position] === "]") depth--;
        position++;
      }
      if (depth) throw new Error("The Newick comment is not closed.");
    } else if ("(),:;".includes(character)) {
      yield { type: character };
      position++;
    } else if (character === "'" || character === '"') {
      const quote = character;
      let label = "";
      let closed = false;
      position++;
      while (position < text.length) {
        const next = text[position++];
        if (next !== quote) {
          label += next;
        } else if (text[position] === quote) {
          label += quote;
          position++;
        } else {
          closed = true;
          break;
        }
      }
      if (!closed) throw new Error("The quoted Newick name is not closed.");
      yield { type: "label", value: label, quoted: true };
    } else {
      const start = position;
      while (position < text.length && !"(),:;[]'\"".includes(text[position])) position++;
      if (position === start) throw new Error("Unexpected character in Newick tree.");
      const label = text.slice(start, position).trim();
      if (label) yield { type: "label", value: label, quoted: false };
    }
  }
}

const DECIMAL_NUMBER = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i;

export function parseNewick(newick) {
  const root = {};
  let current = root;
  const ancestors = [];
  let canStartGroup = true;
  let canName = true;
  let canSetLength = true;
  let expectingLength = false;
  let finished = false;
  let hasTree = false;
  for (const token of newickTokens(newick)) {
    if (finished) throw new Error("Only one Newick tree can be loaded at a time.");
    if (expectingLength) {
      if (token.type !== "label" || token.quoted || !DECIMAL_NUMBER.test(token.value) || !Number.isFinite(Number(token.value))) {
        throw new Error("A branch length must be a finite number.");
      }
      current.length = Number(token.value);
      expectingLength = false;
      canSetLength = false;
      continue;
    }
    switch (token.type) {
      case "(": {
        if (!canStartGroup) throw new Error("Unexpected opening parenthesis in Newick tree.");
        const child = {};
        current.children = [child];
        ancestors.push(current);
        current = child;
        canStartGroup = canName = canSetLength = true;
        hasTree = true;
        break;
      }
      case ",": {
        if (!ancestors.length) throw new Error("A Newick comma must separate child nodes.");
        const sibling = {};
        ancestors[ancestors.length - 1].children.push(sibling);
        current = sibling;
        canStartGroup = canName = canSetLength = true;
        break;
      }
      case ")":
        if (!ancestors.length) throw new Error("Unexpected closing parenthesis in Newick tree.");
        current = ancestors.pop();
        canStartGroup = false;
        canName = canSetLength = true;
        break;
      case ":":
        if (!canSetLength) throw new Error("A node has more than one branch length.");
        expectingLength = true;
        canStartGroup = canName = false;
        hasTree = true;
        break;
      case ";":
        if (ancestors.length) throw new Error("The Newick parentheses are not balanced.");
        finished = true;
        break;
      case "label":
        if (!canName) throw new Error("Unexpected name after a node or branch length.");
        current.name = token.value;
        canStartGroup = canName = false;
        hasTree = true;
        break;
    }
  }
  if (!hasTree) throw new Error("The Newick tree is empty.");
  if (expectingLength) throw new Error("A branch length is missing after ':'.");
  if (ancestors.length) throw new Error("The Newick parentheses are not balanced.");
  return root;
}

export function estimateLeafCount(newick) {
  let commas = 0;
  let hasTree = false;
  try {
    for (const token of newickTokens(newick)) {
      if (token.type === ",") commas++;
      if (token.type !== ";") hasTree = true;
    }
  } catch {
    return 0;
  }
  return hasTree ? commas + 1 : 0;
}

export function getVisibleChildren(node) {
  return node?._collapsed ? null : getAllChildren(node);
}

export function getAllChildren(node) {
  return node?.children?.length ? node.children : null;
}

export function readCheckParam(id) {
  return document.getElementById(id)?.checked || false;
}

export function readSelectParam(id) {
  return document.getElementById(id)?.value ?? null;
}

export function escapeHtmlAttr(value) {
  return String(value ?? "")
    .replace(/&/g, "&amp;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

/**
  * Copy of a leaf-to-cluster map with outlier leaves set to -1: clusters
  * smaller than `threshold`, or singletons when there is none. Negative ids
  * are kept.
 */
export function markOutlierClusters(clusterMap, { enabled = true, threshold = null } = {}) {
  const entries = Object.entries(clusterMap || {});
  const marked = {};
  if (!enabled) {
    for (const [name, clusterId] of entries) marked[name] = clusterId;
    return marked;
  }
  const minimumSize = Number.isSafeInteger(threshold) && threshold >= 1 ? threshold : 2;
  const sizes = new Map();
  for (const [, clusterId] of entries) sizes.set(clusterId, (sizes.get(clusterId) || 0) + 1);
  for (const [name, clusterId] of entries) {
    marked[name] = Number(clusterId) >= 0 && sizes.get(clusterId) < minimumSize ? -1 : clusterId;
  }
  return marked;
}

/** Boundaries halfway between neighbouring sorted positions, extended one half-step past each end. */
export function midpointBoundaries(positions) {
  const count = positions.length;
  const boundaries = new Array(count + 1);
  for (let index = 1; index < count; index++) boundaries[index] = (positions[index - 1] + positions[index]) / 2;
  boundaries[0] = positions[0] - (boundaries[1] - positions[0]);
  boundaries[count] = positions[count - 1] + (positions[count - 1] - boundaries[count - 1]);
  return boundaries;
}

/** Runs of neighbouring items with the same cluster id, as {start, end, cid} with `end` exclusive; unassigned items are skipped. */
export function clusterRuns(count, clusterIdAt) {
  const runs = [];
  if (!count) return runs;
  let runStart = 0;
  let runCid = clusterIdAt(0);
  for (let index = 1; index <= count; index++) {
    const cid = index < count ? clusterIdAt(index) : Symbol("END");
    if (cid === runCid) continue;
    if (runCid != null) runs.push({ start: runStart, end: index, cid: runCid });
    runStart = index;
    runCid = index < count ? cid : null;
  }
  return runs;
}

export function collectLeafNamesData(node, names = []) {
  const pending = node ? [node] : [];
  while (pending.length) {
    const current = pending.pop();
    const children = getAllChildren(current);
    if (children) {
      for (let index = children.length - 1; index >= 0; index--) pending.push(children[index]);
    } else if (current.name !== undefined && current.name !== null && current.name !== "") {
      names.push(String(current.name));
    }
  }
  return names;
}
