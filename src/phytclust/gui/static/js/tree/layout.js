/* PhytClust – tree/layout.js
   Branch-length accumulation and the cartesian / circular layout passes. */

import { state } from "../state.js";
import { treeHost } from "../dom.js";
import { getAllChildren } from "../utils.js";

export function accumulateBranchLength(node, length) {
  length = length || 0;
  node._bl = length;
  if (node.children) {
    for (const c of node.children)
      accumulateBranchLength(c, length + (c.length || 0));
  }
}

export function computeLayouts() {
  if (!state.NEWICK_RAW_TREE) return;
  state.HIER_CART = d3.hierarchy(state.NEWICK_RAW_TREE, getAllChildren);
  state.HIER_CIRC = d3.hierarchy(state.NEWICK_RAW_TREE, getAllChildren);

  var width = treeHost.clientWidth || 800;
  var height = treeHost.clientHeight || 500;
  var margin = { top: 20, right: 80, bottom: 20, left: 80 };
  var innerW = width - margin.left - margin.right;
  var innerH = height - margin.top - margin.bottom;
  var radius = Math.min(innerW, innerH) / 2;

  var maxBl = d3.max(state.HIER_CART.descendants(), (d) => d.data._bl || 0) || 1;
  var blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);
  var blToR = d3.scaleLinear().domain([0, maxBl]).range([0, radius]);

  d3.cluster().size([innerH, 1])(state.HIER_CART);
  // Cladogram: ignore branch lengths — place each node at its topological depth
  // so all leaves align at the right edge. d.y from d3.cluster is in [0,1].
  if ((state.render.layout || "rectangular") === "cladogram") {
    state.HIER_CART.each((d) => {
      d._x = d.x;
      d._y = d.y * innerW;
    });
  } else {
    state.HIER_CART.each((d) => {
      d._x = d.x;
      d._y = blToX(d.data._bl || 0);
    });
  }

  d3.cluster().size([2 * Math.PI, 1])(state.HIER_CIRC);
  state.HIER_CIRC.each((d) => {
    d._angle = d.x;
    d._radius = blToR(d.data._bl || 0);
  });
}
