// Mini scores panel (overlay on tree page) + switchToK.

import { state } from "../state.js";
import { d3Tooltip } from "../dom.js";
import { getThemeColours } from "../colors.js";
import { showToast } from "../ui/toast.js";
import { populateClusterSelector, switchCluster } from "./cluster_editor.js";
import { scoreDataPoints, hasBinRanges, drawBinBands } from "./scores.js";

export function drawMiniScores(data) {
  var panel = document.getElementById("mini-scores-panel");
  var plotEl = document.getElementById("mini-scores-plot");
  if (!panel || !plotEl) return;

  var scores = data && data.scores ? data.scores : [];
  var peaks = data && data.peaks ? data.peaks : [];
  if (!Array.isArray(scores) || scores.length === 0) {
    panel.classList.remove("visible");
    return;
  }

  panel.classList.add("visible");
  plotEl.innerHTML = "";
  var tc = getThemeColours();

  var dataPoints = scoreDataPoints(scores);
  if (!dataPoints.length) return;

  var width = 252,
    height = 124;
  var margin = { top: 8, right: 8, bottom: 20, left: 36 };
  var innerW = width - margin.left - margin.right;
  var innerH = height - margin.top - margin.bottom;

  var svg = d3
    .select(plotEl)
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("background", "transparent");
  var g = svg
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  var xScale = d3
    .scaleLinear()
    .domain([d3.min(dataPoints, (d) => d.k), d3.max(dataPoints, (d) => d.k)])
    .range([0, innerW]);
  var yScale = d3
    .scaleLinear()
    .domain([
      d3.min(dataPoints, (d) => d.score),
      d3.max(dataPoints, (d) => d.score),
    ])
    .nice()
    .range([innerH, 0]);

  // Resolution-mode bin bands
  if (hasBinRanges(data)) {
    var ks = dataPoints.map((d) => d.k);
    drawBinBands(g, data.bin_ranges, {
      kMin: ks[0],
      kMax: ks[ks.length - 1],
      xLeft: xScale,
      xRight: xScale,
      innerWidth: innerW,
      innerHeight: innerH,
      dash: "2,2",
      strokeColour: tc.muted,
    });
  }

  // Axes (minimal)
  g.append("g")
    .attr("transform", "translate(0," + innerH + ")")
    .call(d3.axisBottom(xScale).ticks(4).tickFormat(d3.format("d")))
    .selectAll("text")
    .attr("fill", tc.internal)
    .style("font-size", "9px");
  g.append("g")
    .call(d3.axisLeft(yScale).ticks(3))
    .selectAll("text")
    .attr("fill", tc.internal)
    .style("font-size", "9px");
  g.selectAll(".domain").attr("stroke", tc.branch);
  g.selectAll(".tick line").attr("stroke", tc.branch);

  // Line
  var line = d3
    .line()
    .x((d) => xScale(d.k))
    .y((d) => yScale(d.score))
    .curve(d3.curveMonotoneX);
  g.append("path")
    .datum(dataPoints)
    .attr("fill", "none")
    .attr("stroke", tc.accent)
    .attr("stroke-width", 1.5)
    .attr("d", line);

  // Clickable points
  g.selectAll(".mini-point")
    .data(dataPoints)
    .enter()
    .append("circle")
    .attr("cx", (d) => xScale(d.k))
    .attr("cy", (d) => yScale(d.score))
    .attr("r", (d) => (peaks.indexOf(d.k) >= 0 ? 4 : 2.5))
    .attr("fill", (d) => (peaks.indexOf(d.k) >= 0 ? "#ef4444" : tc.accent))
    .attr("stroke", tc.bg)
    .attr("stroke-width", 1)
    .style("cursor", "pointer")
    .on("click", function (event, d) {
      switchToK(d.k);
    })
    .on("mouseover", function (event, d) {
      d3Tooltip
        .style("opacity", 1)
        .html("k=" + d.k + " score=" + d.score.toFixed(3));
    })
    .on("mousemove", function (event) {
      d3Tooltip
        .style("left", event.pageX + 10 + "px")
        .style("top", event.pageY + 10 + "px");
    })
    .on("mouseout", function () {
      d3Tooltip.style("opacity", 0);
    });
}

export function switchToK(k) {
  if (!state.lastRunData) return;
  if (state.lastRunData.all_clusters && state.lastRunData.all_ks) {
    var idx = state.lastRunData.all_ks.indexOf(k);
    if (idx >= 0) {
      if (state.clusterViewMode !== "all") {
        state.clusterViewMode = "all";
        populateClusterSelector(state.lastRunData);
        var toggleEl = document.getElementById("cluster-view-toggle");
        if (toggleEl) toggleEl.textContent = "Show peaks only";
      }
      var selectEl = document.getElementById("cluster-select");
      if (selectEl) selectEl.selectedIndex = idx;
      switchCluster(idx);
      return;
    }
  }
  var peakKs = state.lastRunData.k_values || state.lastRunData.ks;
  if (peakKs && state.lastRunData.clusters) {
    var idx2 = peakKs.indexOf(k);
    if (idx2 >= 0) {
      if (state.clusterViewMode !== "peaks") {
        state.clusterViewMode = "peaks";
        populateClusterSelector(state.lastRunData);
        var toggleEl2 = document.getElementById("cluster-view-toggle");
        if (toggleEl2) toggleEl2.textContent = "Show all k";
      }
      var selectEl2 = document.getElementById("cluster-select");
      if (selectEl2) selectEl2.selectedIndex = idx2;
      switchCluster(idx2);
      return;
    }
  }
  showToast(
    "k=" + k + " not cached. Enable 'Compute all k' to access any k.",
    "info",
    3000,
  );
}
