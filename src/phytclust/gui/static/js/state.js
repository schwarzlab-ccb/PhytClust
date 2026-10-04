// Shared tree data, display settings, and selection state.

export const state = {
  // Tree data and current results.
  HIER_CART: null,
  HIER_CIRC: null,
  CLUSTER_COLORS: [],
  CURRENT_CLUSTERS: {},
  NEWICK_RAW_TREE: null,

  isRunning: false,
  latestOptimalKData: null,
  latestApiData: null,
  lastRunSignature: null,
  latestRunId: null,
  pendingSessionView: null,
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

  // Display settings.
  render: {
    layout: "rectangular", // rectangular, cladogram, or circular
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
      colorMode: "bars", // bars or boxes
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

// Find an existing display setting without following inherited properties.
function renderOptionLocation(path) {
  if (typeof path !== "string" || !path) {
    throw new Error("A display setting path is required.");
  }
  const names = path.split(".");
  let settings = state.render;
  for (const name of names.slice(0, -1)) {
    if (!settings || typeof settings !== "object" ||
        !Object.hasOwn(settings, name)) {
      throw new Error(`Unknown display setting: ${path}`);
    }
    settings = settings[name];
  }
  const name = names[names.length - 1];
  if (!settings || typeof settings !== "object" ||
      !Object.hasOwn(settings, name)) {
    throw new Error(`Unknown display setting: ${path}`);
  }
  return { settings, name };
}

/** Set an existing display setting. The caller redraws the view. */
export function setRenderOption(path, value) {
  const { settings, name } = renderOptionLocation(path);
  settings[name] = value;
  return value;
}

/** Read a display setting such as "labels.fontSize". */
export function getRenderOption(path) {
  const { settings, name } = renderOptionLocation(path);
  return settings[name];
}

/** Copy the display settings and apply partial overrides. */
export function renderOptionsWithOverrides(overrides = {}) {
  const copy = structuredClone(state.render);
  function applyOverrides(settings, changes, parentPath = "") {
    if (!changes || typeof changes !== "object" || Array.isArray(changes)) {
      throw new Error(`Display overrides must be an object: ${parentPath || "render"}`);
    }
    for (const [name, value] of Object.entries(changes)) {
      const path = parentPath ? `${parentPath}.${name}` : name;
      if (!Object.hasOwn(settings, name)) {
        throw new Error(`Unknown display setting: ${path}`);
      }
      if (settings[name] !== null && typeof settings[name] === "object") {
        applyOverrides(settings[name], value, path);
      } else {
        if (value !== null && typeof value === "object") {
          throw new Error(`Display setting requires a single value: ${path}`);
        }
        settings[name] = value;
      }
    }
  }
  applyOverrides(copy, overrides);
  return copy;
}

/** Redraw with temporary settings, take a synchronous snapshot, then restore. */
export function withRenderOverrides(overrides, redraw, takeSnapshot) {
  const originalSettings = state.render;
  state.render = renderOptionsWithOverrides(overrides);
  try {
    if (typeof redraw === "function") redraw();
    return takeSnapshot();
  } finally {
    state.render = originalSettings;
    if (typeof redraw === "function") redraw();
  }
}

export const LABEL_PAD = 6;

export const EXAMPLE_NEWICK =
  "(((A:5, B:3)C1:6, (C:3, D:7)D1:4)A13:22, (((E:7, F:13)E12:5, G:6)B23:10, H:60):35):0;";

export const SELECTED_CLUSTER_IDS = new Set();
export const BOX_LABEL_MAP = Object.create(null); // Cluster ID to custom name.
export const BOX_ADJUST_MAP = Object.create(null); // Cluster ID to box offsets in pixels.

// Custom settings for each tree node.
export const NODE_CUSTOM = new WeakMap();

export function hasClusterFocus(clusterId) {
  if (clusterId == null) return SELECTED_CLUSTER_IDS.size === 0;
  return (
    SELECTED_CLUSTER_IDS.size === 0 || SELECTED_CLUSTER_IDS.has(Number(clusterId))
  );
}

export function getBoxAdjust(clusterId) {
  const key = String(clusterId);
  if (!BOX_ADJUST_MAP[key]) {
    BOX_ADJUST_MAP[key] = { dx: 0, dy: 0, padX: 0, padY: 0 };
  }
  return BOX_ADJUST_MAP[key];
}
