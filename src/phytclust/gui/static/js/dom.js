/* ============================================================
   PhytClust – dom.js
   Cached DOM element references.
   ============================================================ */

export const newickEl = document.getElementById("newick-input");
export const resultEl = document.getElementById("result");
export const statusEl = document.getElementById("status-message");
export const statusDot = document.getElementById("status-dot");
export const treeHost = document.getElementById("tree_display");

export const extraKEl = document.getElementById("extra-k");
export const extraOutgroupEl = document.getElementById("extra-outgroup");
export const extraRootTaxonEl = document.getElementById("extra-root-taxon");
export const extraTopNEl = document.getElementById("extra-topn");
export const extraResolutionEl = document.getElementById("extra-resolution");
export const extraBinsEl = document.getElementById("extra-bins");
export const extraMaxKEl = document.getElementById("extra-maxk");
export const extraMaxKLimitEl = document.getElementById("extra-maxklimit");
export const extraPromWeightEl = document.getElementById("extra-prominence-weight");
export const extraMinClusterEl = document.getElementById("extra-min-cluster-size");
export const extraOutlierEl = document.getElementById("extra-outlier");

// D3 tooltip
export const d3Tooltip = d3.select("body").append("div").attr("class", "d3-tooltip");
