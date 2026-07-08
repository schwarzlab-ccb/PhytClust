/* ============================================================
   PhytClust – views/optimal_k.js
   Optimal k score plot.
   ============================================================ */

import { state } from "../state.js";
import { d3Tooltip } from "../dom.js";
import { getThemeColors, BASE_COLORS, withAlpha } from "../colors.js";

export function drawOptimalK(data) {
  if (!data) data = state.latestOptimalKData;
  var plotEl = document.getElementById("optimalk_plot");
  if (!plotEl) return;

  var scores = data && data.scores ? data.scores : [];
  var peaks = data && data.peaks ? data.peaks : [];
  const tc = getThemeColors();

  if (!Array.isArray(scores) || scores.length === 0) {
    plotEl.innerHTML =
      '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);">No scores available to plot.</div>';
    window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = null;
    return;
  }

  plotEl.innerHTML = "";
  var dataPoints = [];
  for (let i = 0; i < scores.length - 1; i++)
    dataPoints.push({ k: i + 2, score: scores[i + 1] });

  // With only a single score (k=1) there is nothing to plot at k>=2; guard so
  // the log/linear scale domains don't collapse to NaN and produce broken axes.
  if (dataPoints.length === 0) {
    plotEl.innerHTML =
      '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);">Not enough scores to plot (need k ≥ 2).</div>';
    window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = null;
    return;
  }

  var peakPoints = peaks
    .map((k) => {
      var idx = k - 2;
      return idx >= 0 && idx < dataPoints.length
        ? { k, score: dataPoints[idx].score }
        : null;
    })
    .filter((d) => d);

  var widthScaleEl = document.getElementById("optk-width-scale");
  var widthScale = widthScaleEl ? parseFloat(widthScaleEl.value || "1") : 1;
  if (!isFinite(widthScale) || widthScale < 1) widthScale = 1;

  // Container scrolls horizontally when the SVG is wider than the viewport.
  // Use `scroll` (not `auto`) so the user sees the scrollbar appear.
  plotEl.style.overflowX = widthScale > 1 ? "scroll" : "hidden";
  plotEl.style.overflowY = "hidden";

  // clientWidth must be measured BEFORE the SVG is appended, otherwise an
  // already-wide SVG from the previous render distorts the next baseline.
  var baseWidth = plotEl.clientWidth || 700;
  var width = Math.max(baseWidth, Math.round(baseWidth * widthScale));
  var height = plotEl.clientHeight || 420;
  var margin = { top: 24, right: 24, bottom: 48, left: 58 };
  var innerWidth = width - margin.left - margin.right;
  var innerHeight = height - margin.top - margin.bottom;

  var svg = d3
    .select(plotEl)
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("background", "transparent");
  var g = svg
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  var axisModeEl = document.getElementById("axis-mode");
  var defaultMode = scores.length > 50 ? "log" : "normal";
  var mode = defaultMode;
  if (axisModeEl) {
    if (!axisModeEl.__initialized) {
      axisModeEl.value = defaultMode;
      axisModeEl.__initialized = true;
    }
    if (axisModeEl.__userOverride) mode = axisModeEl.value || defaultMode;
    else axisModeEl.value = defaultMode;
  }

  var xScale, xAxis;
  if (mode === "log") {
    xScale = d3
      .scaleLog()
      .domain([2, d3.max(dataPoints, (d) => d.k)])
      .range([0, innerWidth]);
    xAxis = d3
      .axisBottom(xScale)
      .tickValues(dataPoints.map((d) => d.k))
      .tickFormat(d3.format("d"));
  } else {
    xScale = d3
      .scaleBand()
      .domain(dataPoints.map((d) => d.k))
      .range([0, innerWidth])
      .padding(0.2);
    var nPts = dataPoints.length;
    var dtick = Math.max(1, Math.ceil(nPts / Math.max(10, 10 * widthScale)));
    xAxis = d3
      .axisBottom(xScale)
      .tickValues(dataPoints.map((d) => d.k).filter((d, i) => i % dtick === 0))
      .tickFormat(d3.format("d"));
  }

  const yScale = d3
    .scaleLinear()
    .domain([
      d3.min(dataPoints, (d) => d.score),
      d3.max(dataPoints, (d) => d.score),
    ])
    .nice()
    .range([innerHeight, 0]);

  // Resolution-mode bin bands (behind grid + line)
  const isResolution =
    data &&
    data.mode === "resolution" &&
    Array.isArray(data.bin_ranges) &&
    data.bin_ranges.length > 0;
  let binLabels = null;
  if (isResolution) {
    const ks = dataPoints.map((d) => d.k);
    const kMin = ks[0];
    const kMax = ks[ks.length - 1];
    const clamp = (k) => Math.max(kMin, Math.min(kMax, k));
    const xLeft = (k) => xScale(clamp(k));
    const xRight = (k) =>
      mode === "log"
        ? xScale(clamp(k))
        : xScale(clamp(k)) + xScale.bandwidth();
    binLabels = [];
    data.bin_ranges.forEach(([lo, hi], i) => {
      const rawL = xLeft(lo);
      const rawR = xRight(hi);
      if (!isFinite(rawL) || !isFinite(rawR)) return;
      const x0 = Math.max(0, Math.min(rawL, rawR));
      const x1 = Math.min(innerWidth, Math.max(rawL, rawR));
      if (x1 <= x0) return;
      const colour = BASE_COLORS[i % BASE_COLORS.length];
      g.append("rect")
        .attr("class", "bin-band")
        .attr("x", x0)
        .attr("y", 0)
        .attr("width", x1 - x0)
        .attr("height", innerHeight)
        .attr("fill", withAlpha(colour, 0.1))
        .attr("stroke", "none");
      if (i < data.bin_ranges.length - 1) {
        g.append("line")
          .attr("class", "bin-boundary")
          .attr("x1", x1)
          .attr("x2", x1)
          .attr("y1", 0)
          .attr("y2", innerHeight)
          .attr("stroke", tc.muted)
          .attr("stroke-width", 1)
          .attr("stroke-dasharray", "3,3")
          .attr("opacity", 0.7);
      }
      binLabels.push({
        x: (x0 + x1) / 2,
        label: `CL${i + 1}`,
        range: `k = ${lo}–${hi}`,
        colour,
      });
    });
  }

  // Grid
  g.append("g")
    .attr("class", "grid")
    .attr("transform", `translate(0,${innerHeight})`)
    .call(d3.axisBottom(xScale).tickSize(-innerHeight).tickFormat(""))
    .selectAll("line")
    .attr("stroke", tc.border)
    .attr("stroke-dasharray", "2,2");
  g.append("g")
    .attr("class", "grid")
    .call(d3.axisLeft(yScale).tickSize(-innerWidth).tickFormat(""))
    .selectAll("line")
    .attr("stroke", tc.border)
    .attr("stroke-dasharray", "2,2");
  g.selectAll(".grid .domain").remove();

  // Axes
  var xAxisG = g
    .append("g")
    .attr("transform", `translate(0,${innerHeight})`)
    .call(xAxis);
  var yAxisG = g.append("g").call(d3.axisLeft(yScale).ticks(6));
  [xAxisG, yAxisG].forEach((ag) => {
    ag.selectAll("text").attr("fill", tc.internal).style("font-size", "11px");
    ag.selectAll("line").attr("stroke", tc.branch);
    ag.selectAll("path").attr("stroke", tc.branch);
  });

  if (binLabels && binLabels.length > 0) {
    const labelG = g
      .selectAll(".bin-label")
      .data(binLabels)
      .enter()
      .append("g")
      .attr("class", "bin-label")
      .attr("transform", (d) => `translate(${d.x}, 4)`);
    labelG
      .append("text")
      .attr("text-anchor", "middle")
      .attr("dominant-baseline", "hanging")
      .attr("font-size", 11)
      .attr("font-weight", 600)
      .attr("fill", (d) => d.colour)
      .text((d) => d.label)
      .append("title")
      .text((d) => `${d.label}: ${d.range}`);
  }

  g.append("text")
    .attr("x", innerWidth / 2)
    .attr("y", innerHeight + 38)
    .attr("text-anchor", "middle")
    .attr("font-size", 12)
    .attr("fill", tc.internal)
    .text("k (number of clusters)");
  g.append("text")
    .attr("transform", "rotate(-90)")
    .attr("x", -innerHeight / 2)
    .attr("y", -42)
    .attr("text-anchor", "middle")
    .attr("font-size", 12)
    .attr("fill", tc.internal)
    .text("CalBow Score");

  // Line + points
  const lineColor = tc.accent;
  const line = d3
    .line()
    .x((d) =>
      mode === "log" ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2,
    )
    .y((d) => yScale(d.score))
    .curve(d3.curveMonotoneX);
  g.append("path")
    .datum(dataPoints)
    .attr("fill", "none")
    .attr("stroke", lineColor)
    .attr("stroke-width", 2)
    .attr("d", line);

  g.selectAll(".score-point")
    .data(dataPoints)
    .enter()
    .append("circle")
    .attr("class", "score-point")
    .attr("cx", (d) =>
      mode === "log" ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2,
    )
    .attr("cy", (d) => yScale(d.score))
    .attr("r", 3.5)
    .attr("fill", lineColor)
    .attr("stroke", tc.bg)
    .attr("stroke-width", 1.5)
    .on("mouseover", (event, d) =>
      d3Tooltip
        .style("opacity", 1)
        .html(`k = ${d.k}<br/>score = ${d.score.toFixed(4)}`),
    )
    .on("mousemove", (event) =>
      d3Tooltip
        .style("left", event.pageX + 12 + "px")
        .style("top", event.pageY + 12 + "px"),
    )
    .on("mouseout", () => d3Tooltip.style("opacity", 0));

  // Peaks
  g.selectAll(".peak-point")
    .data(peakPoints)
    .enter()
    .append("path")
    .attr("class", "peak-point")
    .attr(
      "transform",
      (d) =>
        `translate(${mode === "log" ? xScale(d.k) : xScale(d.k) + xScale.bandwidth() / 2},${yScale(d.score)})`,
    )
    .attr("d", d3.symbol().type(d3.symbolDiamond).size(100))
    .attr("fill", "#ef4444")
    .attr("stroke", tc.bg)
    .attr("stroke-width", 1.5)
    .on("mouseover", (event, d) =>
      d3Tooltip
        .style("opacity", 1)
        .html(
          `<strong>Optimal k = ${d.k}</strong><br/>score = ${d.score.toFixed(4)}`,
        ),
    )
    .on("mousemove", (event) =>
      d3Tooltip
        .style("left", event.pageX + 12 + "px")
        .style("top", event.pageY + 12 + "px"),
    )
    .on("mouseout", () => d3Tooltip.style("opacity", 0));

  // Legend
  const legend = g
    .append("g")
    .attr("transform", `translate(${innerWidth - 130}, 8)`);
  legend
    .append("rect")
    .attr("x", -8)
    .attr("y", -8)
    .attr("width", 140)
    .attr("height", 44)
    .attr("fill", tc.bg)
    .attr("rx", 6)
    .attr("stroke", tc.border)
    .attr("opacity", 0.9);
  legend
    .append("line")
    .attr("x1", 0)
    .attr("y1", 6)
    .attr("x2", 18)
    .attr("y2", 6)
    .attr("stroke", lineColor)
    .attr("stroke-width", 2);
  legend
    .append("circle")
    .attr("cx", 9)
    .attr("cy", 6)
    .attr("r", 3)
    .attr("fill", lineColor);
  legend
    .append("text")
    .attr("x", 24)
    .attr("y", 10)
    .attr("font-size", 11)
    .attr("fill", tc.internal)
    .text("Score");
  legend
    .append("path")
    .attr("transform", "translate(9,26)")
    .attr("d", d3.symbol().type(d3.symbolDiamond).size(80))
    .attr("fill", "#ef4444");
  legend
    .append("text")
    .attr("x", 24)
    .attr("y", 30)
    .attr("font-size", 11)
    .attr("fill", tc.internal)
    .text("Optimal k");

  if (axisModeEl && !axisModeEl.__wired) {
    axisModeEl.__wired = true;
    axisModeEl.addEventListener("change", function () {
      axisModeEl.__userOverride = true;
      // Call with no argument so we re-read state.latestOptimalKData: the
      // handler is wired once and would otherwise close over the FIRST render's
      // `data`, redrawing a stale (previous-tree) plot after a new run.
      drawOptimalK();
    });
  }


  const okSvgNode = document.querySelector("#optimalk_plot svg");
  window.PHYTCLUST_ORIGINAL_OPTIMALK_SVG = okSvgNode
    ? new XMLSerializer().serializeToString(okSvgNode)
    : null;
}
