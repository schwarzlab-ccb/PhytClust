/* PhytClust – tree/compare.js
   Standalone renderers for the Compare view: a tree drawn into an arbitrary
   host element, and the cluster-comparison bars. Independent of the main
   drawTree pipeline. */

import { state, LABEL_PAD } from "../state.js";
import { getVisibleChildren } from "../utils.js";
import { generateClusterColors, getThemeColors } from "../colors.js";

export function drawTreeInto(hostEl, clusters) {
  hostEl.innerHTML = "";
  if (!state.NEWICK_RAW_TREE) return;

  const tc = getThemeColors();
  const hier = d3.hierarchy(state.NEWICK_RAW_TREE, getVisibleChildren);
  const width = hostEl.clientWidth || 400;
  const height = hostEl.clientHeight || 400;
  const margin = { top: 16, right: 60, bottom: 16, left: 60 };
  const innerW = width - margin.left - margin.right;
  const innerH = height - margin.top - margin.bottom;

  const maxBl = d3.max(hier.descendants(), (d) => d.data._bl || 0) || 1;
  const blToX = d3.scaleLinear().domain([0, maxBl]).range([0, innerW]);

  d3.cluster().size([innerH, 1])(hier);
  hier.each((d) => {
    d._x = d.x;
    d._y = blToX(d.data._bl || 0);
  });

  const nClusters = clusters ? Math.max(0, ...Object.values(clusters)) + 1 : 0;
  const colors = nClusters > 0 ? generateClusterColors(nClusters) : [];

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
  svg.call(
    d3
      .zoom()
      .scaleExtent([0.3, 8])
      .on("zoom", (event) => {
        zoomLayer.attr("transform", event.transform);
      }),
  );

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

  const allNodes = hier.descendants();
  const node = g
    .append("g")
    .selectAll(".tree-node")
    .data(allNodes)
    .enter()
    .append("g")
    .attr("class", "tree-node")
    .attr("transform", (d) => `translate(${d._y},${d._x})`);

  node
    .append("circle")
    .attr("r", (d) => (d.children && d.children.length ? 1.5 : 2.5))
    .attr("fill", function (d) {
      const name = d.data && d.data.name;
      if (name && !(d.children && d.children.length) && clusters) {
        const cid = clusters[name];
        if (cid != null && colors.length) return colors[cid % colors.length];
      }
      return tc.label;
    });

  node
    .filter((d) => !(d.children && d.children.length))
    .append("text")
    .attr("dy", 3)
    .attr("x", 6)
    .style("text-anchor", "start")
    .style("font-size", "8px")
    .attr("fill", tc.label)
    .text((d) => d.data.name || "");

  // Side bars
  if (clusters && colors.length) {
    const leaves = allNodes
      .filter((d) => !(d.children && d.children.length))
      .sort((a, b) => a._x - b._x);
    if (leaves.length > 1) {
      const estLabel =
        d3.max(leaves, (d) => (d.data.name || "").length) * 5 + 12;
      const barX = (d3.max(allNodes, (d) => d._y) || 0) + estLabel;
      const ys = leaves.map((d) => d._x);
      const boundaries = new Array(leaves.length + 1);
      for (let i = 1; i < leaves.length; i++)
        boundaries[i] = (ys[i - 1] + ys[i]) / 2;
      boundaries[0] = ys[0] - (boundaries[1] - ys[0]);
      boundaries[leaves.length] =
        ys[leaves.length - 1] +
        (ys[leaves.length - 1] - boundaries[leaves.length - 1]);

      let runStart = 0,
        runCid = clusters[leaves[0].data.name];
      for (let i = 1; i <= leaves.length; i++) {
        const cid =
          i < leaves.length ? clusters[leaves[i].data.name] : Symbol("END");
        if (cid !== runCid) {
          if (runCid != null) {
            g.append("rect")
              .attr("x", barX)
              .attr("y", boundaries[runStart])
              .attr("width", 14)
              .attr("height", Math.max(1, boundaries[i] - boundaries[runStart]))
              .attr("fill", colors[runCid % colors.length])
              .attr("opacity", 0.7);
          }
          runStart = i;
          runCid = i < leaves.length ? cid : null;
        }
      }
    }
  }
}

/* ─────────────────────────────────────
   Draw one tree with side-by-side bars
   ───────────────────────────────────── */
export function drawComparisonBarsInto(hostEl, comparisons) {
  hostEl.innerHTML = "";
  if (!state.NEWICK_RAW_TREE) return;

  const tc = getThemeColors();
  const hier = d3.hierarchy(state.NEWICK_RAW_TREE, getVisibleChildren);
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

  state.LAST_COMPARE_SVG = svg;
  state.LAST_COMPARE_ZOOM = zoom;
  state.LAST_COMPARE_LAYER = zoomLayer;

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

  const highlightChanged = !!state.COMPARE_HIGHLIGHT_CHANGES;
  const changedSet = state.COMPARE_CHANGED_LEAVES || new Set();
  const changedColor = "#f97316";

  node
    .filter((d) => !(d.children && d.children.length))
    .append("text")
    .attr("dy", 3)
    .attr("x", state.render.nodes.leafRadius + LABEL_PAD)
    .style("text-anchor", "start")
    .style("font-size", state.render.labels.fontSize + "px")
    .attr("fill", (d) =>
      highlightChanged && changedSet.has(d.data.name) ? changedColor : tc.label,
    )
    .attr("font-weight", (d) =>
      highlightChanged && changedSet.has(d.data.name) ? "700" : "normal",
    )
    .text((d) => (state.render.labels.show ? d.data.name || "" : ""));

  const barX =
    (d3.max(hier.descendants(), (d) => d._y) || 0) + estLabelWidth + 8;
  const ys = leaves.map((d) => d._x);
  if (leaves.length > 1) {
    const boundaries = new Array(leaves.length + 1);
    for (let i = 1; i < leaves.length; i++)
      boundaries[i] = (ys[i - 1] + ys[i]) / 2;
    boundaries[0] = ys[0] - (boundaries[1] - ys[0]);
    boundaries[leaves.length] =
      ys[leaves.length - 1] +
      (ys[leaves.length - 1] - boundaries[leaves.length - 1]);

    function drawBarColumn(map, x, colIdx) {
      if (!map) return;
      const maxCid = Math.max(-1, ...Object.values(map).map(Number));
      const colors = maxCid >= 0 ? generateClusterColors(maxCid + 1) : [];
      let runStart = 0;
      let runCid = map[leaves[0].data.name];
      for (let i = 1; i <= leaves.length; i++) {
        const cid =
          i < leaves.length ? map[leaves[i].data.name] : Symbol("END");
        if (cid !== runCid) {
          if (runCid != null) {
            const nCid = Number(runCid);
            const isOutlier = nCid < 0;
            const fillColor = isOutlier
              ? tc.internal
              : colors.length
                ? colors[nCid % colors.length]
                : tc.label;
            g.append("rect")
              .attr("x", x)
              .attr("y", boundaries[runStart])
              .attr("width", colW)
              .attr("height", Math.max(1, boundaries[i] - boundaries[runStart]))
              .attr("fill", fillColor)
              .attr("opacity", 0.76)
              .attr("pointer-events", "none");
          }
          runStart = i;
          runCid = i < leaves.length ? cid : null;
        }
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

/* ─────────────────────────────────────
   Fit compare view to viewport
   ───────────────────────────────────── */
export function fitCompare() {
  const svg = state.LAST_COMPARE_SVG;
  const zoom = state.LAST_COMPARE_ZOOM;
  const layer = state.LAST_COMPARE_LAYER;
  if (!svg || !zoom || !layer) return;

  let bbox;
  try { bbox = layer.node().getBBox(); } catch { return; }
  if (!bbox || !isFinite(bbox.width) || bbox.width <= 0 || bbox.height <= 0) return;

  const svgNode = svg.node();
  const vw = (svgNode && svgNode.clientWidth) || 800;
  const vh = (svgNode && svgNode.clientHeight) || 600;
  const pad = 28;

  const scale = Math.min(
    (vw - pad * 2) / bbox.width,
    (vh - pad * 2) / bbox.height,
    4,
  );
  const tx = pad + (vw - pad * 2 - bbox.width * scale) / 2 - bbox.x * scale;
  const ty = pad + (vh - pad * 2 - bbox.height * scale) / 2 - bbox.y * scale;

  svg.transition().duration(280).call(zoom.transform, d3.zoomIdentity.translate(tx, ty).scale(scale));
}

