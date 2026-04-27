/* ============================================================
   PhytClust – state.js
   All mutable application state in one place.
   ============================================================ */

export const state = {
  // --- Data state ---------------------------------------------------------
  HIER_CART: null,
  HIER_CIRC: null,
  CLUSTER_COLORS: [],
  CURRENT_CLUSTERS: {},
  NEWICK_RAW_TREE: null,

  isRunning: false,
  latestOptimalKData: null,
  latestApiData: null,
  latestRunId: null,
  LAST_TREE_SVG: null,
  LAST_TREE_ZOOM: null,
  LAST_ZOOM_LAYER: null,
  LAST_TREE_TRANSFORM: d3.zoomIdentity,
  LAST_COMPARE_SVG: null,
  LAST_COMPARE_ZOOM: null,
  LAST_COMPARE_LAYER: null,
  SEARCH_TERM: "",

  CLUSTER_VIEW_MODE: "peaks",
  COMPARE_CONFIGS: [],
  COMPARE_CONFIG_NEXT_ID: 1,
  COMPARE_HIGHLIGHT_CHANGES: false,
  COMPARE_REFERENCE_INDEX: 0,
  COMPARE_CHANGED_LEAVES: new Set(),

  CTX_TARGET_DATA: null,
  runHistory: [],

  // --- Render options (single source of truth for "what the tree looks like") ---
  // Read/write via state.render.* directly, or via setRenderOption /
  // getRenderOption for dotted-path access.
  render: {
    layout: "rectangular", // "rectangular" | "circular"
    branches: {
      width: 1.2,
      color: null,       // null = use cluster colour
      colorByClusters: false,
    },
    labels: {
      show: true,
      internalShow: false,
      fontSize: 9,
      axisFontSize: 10,
    },
    nodes: {
      leafRadius: 3.0,
      internalRadius: 1.8,
    },
    clusters: {
      colorMode: "bars", // "bars" | "boxes" | ...
      labelFontSize: 8,
      showBoxLabels: true,
      showOutlierBoxes: true,
      boxAlpha: 0.12,
      boxPadV: 8,
      boxPadH: 6,
      boxCornerRadius: 6,
    },
    scale: {
      width: 1.0,
      height: 1.0,
    },
    compare: {
      barWidth: 16,
      barGap: 6,
      showColumnTitles: true,
    },
  },
};

/**
 * Set a render option by dotted path, e.g. setRenderOption("labels.fontSize", 12).
 * Does not trigger a redraw — callers handle that. Returns the value written.
 */
export function setRenderOption(path, value) {
  const parts = path.split(".");
  let cur = state.render;
  for (let i = 0; i < parts.length - 1; i++) cur = cur[parts[i]];
  cur[parts[parts.length - 1]] = value;
  return value;
}

/**
 * Read a render option by dotted path, e.g. getRenderOption("labels.fontSize").
 */
export function getRenderOption(path) {
  const parts = path.split(".");
  let cur = state.render;
  for (const k of parts) cur = cur[k];
  return cur;
}

/**
 * Return a deep copy of state.render with the given overrides merged in.
 * Useful for "publication preset" exports that want a tweaked snapshot
 * without mutating the live render state.
 */
export function renderOptionsWithOverrides(overrides = {}) {
  const clone = JSON.parse(JSON.stringify(state.render));
  const merge = (dst, src) => {
    for (const [k, v] of Object.entries(src)) {
      if (v && typeof v === "object" && !Array.isArray(v)) {
        if (!dst[k] || typeof dst[k] !== "object") dst[k] = {};
        merge(dst[k], v);
      } else {
        dst[k] = v;
      }
    }
  };
  merge(clone, overrides);
  return clone;
}

/**
 * Temporarily swap state.render for a copy with `overrides` merged in,
 * invoke `fn`, then restore the original. Useful for export presets:
 * the live view is briefly redrawn under the preset, snapshotted, and
 * restored. Returns whatever `fn` returns.
 *
 * `redraw` is called once after applying overrides and once after restoring.
 */
export function withRenderOverrides(overrides, redraw, fn) {
  const original = state.render;
  state.render = renderOptionsWithOverrides(overrides);
  try {
    if (typeof redraw === "function") redraw();
    return fn();
  } finally {
    state.render = original;
    if (typeof redraw === "function") redraw();
  }
}

export const LABEL_PAD = 6;

export const EXAMPLE_NEWICK =
  "(((A:5, B:3)C1:6, (C:3, D:7)D1:4)A13:22, (((E:7, F:13)E12:5, G:6)B23:10, H:60):35):0;";

export const SELECTED_CLUSTER_IDS = new Set();
export const BOX_LABEL_MAP = {}; // cid -> custom name, editable via context menu
export const BOX_ADJUST_MAP = {}; // cid -> {dx, dy, padX, padY} in px

// Per-node customizations: keyed by data node reference
export const NODE_CUSTOM = new WeakMap();

export function hasClusterFocus(cid) {
  if (cid == null) return SELECTED_CLUSTER_IDS.size === 0;
  return (
    SELECTED_CLUSTER_IDS.size === 0 || SELECTED_CLUSTER_IDS.has(Number(cid))
  );
}

export function getBoxAdjust(cid) {
  const key = String(cid);
  if (!BOX_ADJUST_MAP[key]) {
    BOX_ADJUST_MAP[key] = { dx: 0, dy: 0, padX: 0, padY: 0 };
  }
  return BOX_ADJUST_MAP[key];
}
