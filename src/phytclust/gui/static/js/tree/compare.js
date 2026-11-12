/* Standalone renderer for the Compare view: one tree with cluster-comparison
   bars, independent of the main drawTree pipeline. */

import { state, LABEL_PAD } from "../state.js";
import { getVisibleChildren, midpointBoundaries, clusterRuns } from "../utils.js";
import { generateClusterColours, getThemeColours } from "../colors.js";
import { fitToViewport } from "./draw.js";

export function drawComparisonBarsInto(hostEl, comparisons) {
  hostEl.innerHTML = "";
  if (!state.newickTree) return;

  const tc = getThemeColours();
  const hier = d3.hierarchy(state.newickTree, getVisibleChildren);
  const width = hostEl.clientWidth || 840;
  const height = hostEl.clientHeight || 500;
  const compareOpts = state.render.compare || {};
  const colW = Number(compareOpts.barWidth) || 16;
  const gap = Number(compareOpts.barGap) != null ? Number(compareOpts.barGap) : 6;
  const showColTitles = compareOpts.showColumnTitles !== false;
  const nCols = Math.max(0, (comparisons || []).length);
  const margin = {
    top: 28,
    right: Math.max(170, nCols * (colW + gap) + 46),
    bottom: 64,
    left: 60,
  };
  const innerW = width - margin.left - margin.right;
  const innerH = height - margin.top - margin.bottom;

  const maxBl = d3.max(hier.descendants(), (d) => d.data._bl || 0) || 1;
  const blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);

  d3.cluster().size([innerH, 1])(hier);
  const useCladogram = (state.render.layout || "rectangular") === "cladogram";
  hier.each((d) => {
    d._x = d.x;
    d._y = useCladogram ? d.y * innerW : blToX(d.data._bl || 0);
  });

  const leaves = hier
    .descendants()
    .filter((d) => !(d.children && d.children.length))
    .sort((a, b) => a._x - b._x);
  const leafLabels = leaves.map((d) => (d.data.name || "").length);
  const maxLabelChars = d3.max(leafLabels) || 0;
  const estLabelWidth = state.render.labels.show
    ? maxLabelChars * state.render.labels.fontSize * 0.6 + LABEL_PAD + state.render.nodes.leafRadius
    : 20;

  const svg = d3
    .select(hostEl)
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("background", "transparent");
  const zoomLayer = svg.append("g");
  const g = zoomLayer
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);
  const zoom = d3
    .zoom()
    .scaleExtent([0.3, 8])
    .on("zoom", (event) => {
      zoomLayer.attr("transform", event.transform);
    });
  svg.call(zoom);

  state.lastCompareSvg = svg;
  state.lastCompareZoom = zoom;
  state.lastCompareLayer = zoomLayer;

  g.append("g")
    .selectAll(".tree-link")
    .data(hier.links())
    .enter()
    .append("path")
    .attr("class", "tree-link")
    .attr("fill", "none")
    .attr("stroke", tc.branch)
    .attr("stroke-opacity", 1)
    .attr("stroke-width", state.render.branches.width)
    .attr(
      "d",
      (d) =>
        "M" +
        d.source._y +
        "," +
        d.source._x +
        "V" +
        d.target._x +
        "H" +
        d.target._y,
    );

  const node = g
    .append("g")
    .selectAll(".tree-node")
    .data(hier.descendants())
    .enter()
    .append("g")
    .attr("class", "tree-node")
    .attr("transform", (d) => `translate(${d._y},${d._x})`);

  node
    .append("circle")
    .attr("r", (d) =>
      d.children && d.children.length ? state.render.nodes.internalRadius : state.render.nodes.leafRadius,
    )
    .attr("fill", tc.label)
    .attr("stroke", "none");

  const highlightChanged = !!state.compareHighlightChanges;
  const changedSet = state.compareChangedLeaves || new Set();
  const changedColour = "#f97316";

  node
    .filter((d) => !(d.children && d.children.length))
    .append("text")
    .attr("dy", 3)
    .attr("x", state.render.nodes.leafRadius + LABEL_PAD)
    .style("text-anchor", "start")
    .style("font-size", state.render.labels.fontSize + "px")
    .attr("fill", (d) =>
      highlightChanged && changedSet.has(d.data.name) ? changedColour : tc.label,
    )
    .attr("font-weight", (d) =>
      highlightChanged && changedSet.has(d.data.name) ? "700" : "normal",
    )
    .text((d) => (state.render.labels.show ? d.data.name || "" : ""));

  const barX =
    (d3.max(hier.descendants(), (d) => d._y) || 0) + estLabelWidth + 8;
  const ys = leaves.map((d) => d._x);
  if (leaves.length > 1) {
    const boundaries = midpointBoundaries(ys);

    function drawBarColumn(map, x, colIdx) {
      if (!map) return;
      const maxCid = Math.max(-1, ...Object.values(map).map(Number));
      const colours = maxCid >= 0 ? generateClusterColours(maxCid + 1) : [];
      for (const run of clusterRuns(leaves.length, (i) => map[leaves[i].data.name])) {
        const nCid = Number(run.cid);
        const isOutlier = nCid < 0;
        const fillColour = isOutlier
          ? tc.internal
          : colours.length
            ? colours[nCid % colours.length]
            : tc.label;
        g.append("rect")
          .attr("x", x)
          .attr("y", boundaries[run.start])
          .attr("width", colW)
          .attr("height", Math.max(1, boundaries[run.end] - boundaries[run.start]))
          .attr("fill", fillColour)
          .attr("opacity", 0.76)
          .attr("pointer-events", "none");
      }

      if (showColTitles) {
        const title =
          comparisons[colIdx] && comparisons[colIdx].title
            ? comparisons[colIdx].title
            : "Bar " + (colIdx + 1);
        const shortTitle = title.length > 14 ? title.slice(0, 12) + ".." : title;
        g.append("text")
          .attr("x", x + colW / 2)
          .attr("y", -8)
          .attr("text-anchor", "middle")
          .attr("fill", tc.internal)
          .attr("font-size", "10px")
          .text(shortTitle);
      }
    }

    for (let i = 0; i < (comparisons || []).length; i++) {
      drawBarColumn(comparisons[i].clusters, barX + i * (colW + gap), i);
    }
  }

  // Axis only for non-cladogram layouts
  if (!useCladogram) {
    const axisYPos = innerH + 20;
    const axisG = g
      .append("g")
      .attr("class", "branch-length-axis")
      .attr("transform", "translate(0, " + axisYPos + ")")
      .call(d3.axisBottom(blToX).ticks(5));
    axisG
      .selectAll("text")
      .attr("fill", tc.internal)
      .style("font-size", state.render.labels.axisFontSize + "px");
    axisG.selectAll("line").attr("stroke", tc.branch);
    axisG.selectAll("path").attr("stroke", tc.branch);
    axisG
      .append("text")
      .attr("x", innerW / 2)
      .attr("y", state.render.labels.axisFontSize + 20)
      .attr("text-anchor", "middle")
      .attr("font-size", state.render.labels.axisFontSize)
      .attr("fill", tc.internal)
      .text("Branch length");
  }
}

export function fitCompare() {
  fitToViewport(state.lastCompareSvg, state.lastCompareZoom, state.lastCompareLayer);
}
