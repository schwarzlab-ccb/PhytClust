// Shared tree data, display settings, and selection state.

export const state = {
  // Tree data and current results.
  cartesianHierarchy: null,
  circularHierarchy: null,
  clusterColours: [],
  currentClusters: {}, // Displayed map: outlier leaves are -1.
  runClusters: {}, // The run's map, before outlier marking.
  newickTree: null,

  isRunning: false,
  lastScoresData: null,
  lastRunData: null,
  lastRunSignature: null,
  lastRunId: null,
  pendingSessionView: null,
  lastTreeSvg: null,
  lastTreeZoom: null,
  lastZoomLayer: null,
  lastTreeTransform: d3.zoomIdentity,
  lastCompareSvg: null,
  lastCompareZoom: null,
  lastCompareLayer: null,
  searchTerm: "",

  clusterViewMode: "peaks",
  compareConfigs: [],
  compareConfigNextId: 1,
  compareHighlightChanges: false,
  compareReferenceIndex: 0,
  compareChangedLeaves: new Set(),

  contextTargetData: null,
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
      colorMode: "bars", // bars, boxes, or nodes
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

export const selectedClusterIds = new Set();
export const boxLabelMap = Object.create(null); // Cluster ID to custom name.
export const boxAdjustMap = Object.create(null); // Cluster ID to box offsets in pixels.

// Custom settings for each tree node.
export const nodeCustom = new WeakMap();

export function hasClusterFocus(clusterId) {
  if (clusterId == null) return selectedClusterIds.size === 0;
  return (
    selectedClusterIds.size === 0 || selectedClusterIds.has(Number(clusterId))
  );
}

export function getBoxAdjust(clusterId) {
  const key = String(clusterId);
  if (!boxAdjustMap[key]) {
    boxAdjustMap[key] = { dx: 0, dy: 0, padX: 0, padY: 0 };
  }
  return boxAdjustMap[key];
}
