/* PhytClust – tree/layout.js
   Branch-length accumulation and the cartesian / circular layout passes. */

import { state } from "../state.js";
import { treeHost } from "../dom.js";
import { getAllChildren } from "../utils.js";

/** Accumulate branch lengths without depending on JavaScript's call-stack limit. */
export function accumulateBranchLength(node, length = 0) {
  const pending = [[node, length]];
  while (pending.length) {
    const [current, distance] = pending.pop();
    current._bl = distance;
    for (const child of current.children || []) {
      pending.push([child, distance + (child.length || 0)]);
    }
  }
}

// Docked score preview: panel width (260) plus its offset from the edge (16).
const SCORE_PANEL_RESERVE = 276;

/** Width to keep free on the right of the tree while the score preview is docked there. */
export function scorePanelReserve() {
  const data = state.lastRunData;
  if (!data || !Array.isArray(data.scores) || !data.scores.length) return 0;
  const panel = document.getElementById("mini-scores-panel");
  if (panel && panel.style.position === "fixed") return 0;
  return SCORE_PANEL_RESERVE;
}

export function computeLayouts() {
  if (!state.newickTree) return;
  state.cartesianHierarchy = d3.hierarchy(state.newickTree, getAllChildren);
  state.circularHierarchy = d3.hierarchy(state.newickTree, getAllChildren);

  var width = treeHost.clientWidth || 800;
  var height = treeHost.clientHeight || 500;
  var margin = { top: 20, right: 80 + scorePanelReserve(), bottom: 20, left: 80 };
  var innerW = width - margin.left - margin.right;
  var innerH = height - margin.top - margin.bottom;
  var radius = Math.min(innerW, innerH) / 2;

  var maxBl = d3.max(state.cartesianHierarchy.descendants(), (d) => d.data._bl || 0) || 1;
  var blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);
  var blToR = d3.scaleLinear().domain([0, maxBl]).range([0, radius]);

  d3.cluster().size([innerH, 1])(state.cartesianHierarchy);
  // Cladogram: ignore branch lengths — place each node at its topological depth
  // so all leaves align at the right edge. d.y from d3.cluster is in [0,1].
  if ((state.render.layout || "rectangular") === "cladogram") {
    state.cartesianHierarchy.each((d) => {
      d._x = d.x;
      d._y = d.y * innerW;
    });
  } else {
    state.cartesianHierarchy.each((d) => {
      d._x = d.x;
      d._y = blToX(d.data._bl || 0);
    });
  }

  d3.cluster().size([2 * Math.PI, 1])(state.circularHierarchy);
  state.circularHierarchy.each((d) => {
    d._angle = d.x;
    d._radius = blToR(d.data._bl || 0);
  });
}
