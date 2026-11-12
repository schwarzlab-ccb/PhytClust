/* Per-node metadata (custom overrides, display name, radius, cluster
   identity) and the node context menu. */

import { state, nodeCustom } from "../state.js";
import { collectLeafNamesData } from "../utils.js";

export function getNodeCustom(dataNode) {
  if (!nodeCustom.has(dataNode))
    nodeCustom.set(dataNode, {
      highlighted: false,
      radiusScale: 1,
      renamedTo: null,
      lockOriginalSize: false,
    });
  return nodeCustom.get(dataNode);
}

export function nodeDisplayName(dataNode) {
  const c = nodeCustom.has(dataNode) ? nodeCustom.get(dataNode) : null;
  return c && c.renamedTo != null ? c.renamedTo : dataNode.name || "";
}

export function countLeafDescendantsData(node) {
  if (!node) return 0;
  if (!node.children || !node.children.length) return 1;
  var total = 0;
  for (var i = 0; i < node.children.length; i++) {
    total += countLeafDescendantsData(node.children[i]);
  }
  return total;
}

export function customNodeRadius(d) {
  let base =
    d.data && d.data.children && d.data.children.length
      ? state.render.nodes.internalRadius
      : state.render.nodes.leafRadius;
  const c = nodeCustom.has(d.data) ? nodeCustom.get(d.data) : null;

  if (c && c.lockOriginalSize) {
    return base;
  }

  if (
    d.data &&
    d.data._collapsed &&
    d.data.children &&
    d.data.children.length
  ) {
    var leafCount = countLeafDescendantsData(d.data);
    if (leafCount > 1) {
      var autoCollapsed =
        state.render.nodes.internalRadius + state.render.nodes.leafRadius * 0.55 * Math.log2(leafCount);
      base = Math.max(base, Math.min(14, autoCollapsed));
    }
  }

  return base * (c ? c.radiusScale : 1);
}

export function representativeClusterIdFromData(dataNode) {
  if (!dataNode || !state.currentClusters) return null;
  const leaves = [];
  collectLeafNamesData(dataNode, leaves);
  if (!leaves.length) return null;
  let cid = null;
  for (const nm of leaves) {
    const v = state.currentClusters[nm];
    if (v == null) return null;
    if (cid == null) cid = v;
    else if (cid !== v) return null;
  }
  return cid;
}

export function showNodeContextMenu(event, d) {
  event.preventDefault();
  event.stopPropagation();
  state.contextTargetData = d.data;
  const menu = document.getElementById("node-context-menu");
  if (!menu) return;
  menu.style.display = "block";
  menu.style.left = event.clientX + "px";
  menu.style.top = event.clientY + "px";
  var ctxSlider = document.getElementById("ctx-size-slider");
  var ctxSliderVal = document.getElementById("ctx-size-value");
  var c = getNodeCustom(d.data);
  if (ctxSlider) ctxSlider.value = c.radiusScale;
  if (ctxSliderVal) ctxSliderVal.textContent = c.radiusScale.toFixed(1);

  var isLeaf = !(d.data.children && d.data.children.length);
  var hasClusters =
    state.currentClusters && Object.keys(state.currentClusters).length > 0;
  var isBoxMode = state.render.clusters.colorMode === "boxes";

  var copyBtn = document.getElementById("ctx-copy-subtree");
  var originalSizeBtn = document.getElementById("ctx-original-size");
  var clusterBtn = document.getElementById("ctx-select-cluster");
  var renameClusterBtn = document.getElementById("ctx-rename-cluster");
  var adjustBoxBtn = document.getElementById("ctx-adjust-box");

  if (copyBtn) copyBtn.style.display = isLeaf ? "none" : "";
  if (originalSizeBtn) {
    originalSizeBtn.style.display = isLeaf ? "none" : "";
    originalSizeBtn.textContent = c.lockOriginalSize
      ? "Use auto collapsed size"
      : "Keep original size";
  }
  if (clusterBtn) clusterBtn.style.display = hasClusters ? "" : "none";
  if (renameClusterBtn)
    renameClusterBtn.style.display = isBoxMode && hasClusters ? "" : "none";
  if (adjustBoxBtn)
    adjustBoxBtn.style.display = isBoxMode && hasClusters ? "" : "none";
}

export function hideNodeContextMenu() {
  const menu = document.getElementById("node-context-menu");
  if (menu) menu.style.display = "none";
  state.contextTargetData = null;
}

export function isSearchMatch(dataNode) {
  if (!state.searchTerm) return false;
  const name = nodeDisplayName(dataNode).toLowerCase();
  return name.includes(state.searchTerm.toLowerCase());
}
