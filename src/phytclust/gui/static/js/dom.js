// Cached DOM element references.

export const newickEl = document.getElementById("newick-input");
export const statusEl = document.getElementById("status-message");
export const statusDot = document.getElementById("status-dot");
export const treeHost = document.getElementById("tree-display");

export const outgroupEl = document.getElementById("outgroup");
export const rootTaxonEl = document.getElementById("root-taxon");

// D3 tooltip
export const d3Tooltip = d3.select("body").append("div").attr("class", "d3-tooltip");
