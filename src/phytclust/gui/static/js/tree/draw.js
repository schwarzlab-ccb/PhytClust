/* The main tree renderer: drawTree (rectangular + circular), fitTree, and the
   collapse / zoom / tooltip helpers it owns. Layout, per-node metadata and the
   Compare renderers live in tree/layout.js, tree/nodes.js and tree/compare.js. */

import { state, LABEL_PAD, selectedClusterIds, boxLabelMap, nodeCustom, hasClusterFocus, getBoxAdjust } from "../state.js";
import { treeHost } from "../dom.js";
import { escapeHtmlAttr, midpointBoundaries, clusterRuns } from "../utils.js";
import { getThemeColours } from "../colors.js";
import { customNodeRadius, isSearchMatch, nodeDisplayName, countLeafDescendantsData, representativeClusterIdFromData, showNodeContextMenu } from "./nodes.js";
import { computeLayouts, scorePanelReserve } from "./layout.js";

const HIGHLIGHT_COLOR = "#f59e0b";

function setEmptyStateVisible(visible) {
  // Inline display, so no stylesheet rule overrides it.
  var es = document.getElementById("tree-empty-state");
  if (es) es.style.display = visible ? "flex" : "none";
}

function positionTooltip(tip, cx, cy) {
  var tw = tip.offsetWidth || 200;
  var th = tip.offsetHeight || 80;
  var pad = 14;
  var x = cx + pad;
  var y = cy + pad;
  if (x + tw > window.innerWidth - 8) x = cx - tw - pad;
  if (y + th > window.innerHeight - 8) y = cy - th - pad;
  tip.style.left = x + "px";
  tip.style.top = y + "px";
}

export function clearTree() {
  treeHost.innerHTML = "";
  // Show the empty state again when the tree is cleared; drawTree hides it when it draws.
  setEmptyStateVisible(!state.newickTree);
  var tip = document.getElementById("node-tooltip");
  if (tip) tip.setAttribute("hidden", "");
}

export function clearAllCollapsedFlags(node) {
  if (!node) return;
  if (node._collapsed) delete node._collapsed;
  if (node.children && node.children.length) {
    for (const c of node.children) clearAllCollapsedFlags(c);
  }
}

function applyTreeZoomBounds(svg, zoom, zoomLayer, previousTransform) {
  if (!svg || !zoom || !zoomLayer) return;
  let bbox = null;
  try {
    bbox = zoomLayer.node().getBBox();
  } catch {
    bbox = null;
  }
  // Keeps at least half a viewport of the tree in view when panning.
  const svgNode = svg.node();
  const vw = (svgNode && svgNode.clientWidth) || 800;
  const vh = (svgNode && svgNode.clientHeight) || 600;
  const padX = Math.max(vw * 0.5, 200);
  const padY = Math.max(vh * 0.5, 200);
  if (
    bbox &&
    isFinite(bbox.x) &&
    isFinite(bbox.y) &&
    isFinite(bbox.width) &&
    isFinite(bbox.height) &&
    bbox.width > 0 &&
    bbox.height > 0
  ) {
    zoom.translateExtent([
      [bbox.x - padX, bbox.y - padY],
      [bbox.x + bbox.width + padX, bbox.y + bbox.height + padY],
    ]);
  }

  svg.call(zoom);
  svg.call(zoom.transform, previousTransform || d3.zoomIdentity);
}

function clearDescendantCollapsedFlags(node) {
  if (!node || !node.children || !node.children.length) return;
  for (const c of node.children) {
    if (c._collapsed) delete c._collapsed;
    clearDescendantCollapsedFlags(c);
  }
}

function toggleNodeCollapsedState(dataNode) {
  if (!dataNode || !dataNode.children || !dataNode.children.length) return;
  if (dataNode._collapsed) {
    delete dataNode._collapsed;
    return;
  }
  dataNode._collapsed = true;
  clearDescendantCollapsedFlags(dataNode);
}

function findMRCA(nodes) {
  if (!nodes.length) return null;
  if (nodes.length === 1) return nodes[0];
  function ancestors(n) {
    const a = [];
    let cur = n;
    while (cur) {
      a.push(cur);
      cur = cur.parent;
    }
    return a;
  }
  let common = ancestors(nodes[0]);
  for (let i = 1; i < nodes.length; i++) {
    const anc = new Set(ancestors(nodes[i]));
    common = common.filter((n) => anc.has(n));
  }
  return common[0];
}

export function drawTree() {
  clearTree();
  if (!state.newickTree) {
    setEmptyStateVisible(true);
    return;
  }
  setEmptyStateVisible(false);
  computeLayouts();

  const layoutMode = state.render.layout || "rectangular";
  const colourMode = state.render.clusters.colorMode || "bars";
  const tc = getThemeColours();

  const container = d3.select("#tree-display");
  const width = treeHost.clientWidth || 800;
  const height = treeHost.clientHeight || 500;
  const margin = { top: 20, right: 80, bottom: 20, left: 80 };

  const svg = container
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("background", "transparent")
    .style("cursor", "grab");

  // Catches drags on empty canvas, which would otherwise not pan.
  svg
    .append("rect")
    .attr("class", "tree-zoom-catcher")
    .attr("width", "100%")
    .attr("height", "100%")
    .attr("fill", "transparent")
    .attr("pointer-events", "all");

  const zoomLayer = svg.append("g");
  const g = zoomLayer
    .append("g")
    .attr(
      "transform",
      layoutMode === "circular"
        ? `translate(${width / 2},${height / 2})`
        : `translate(${margin.left},${margin.top})`,
    );

  const zoom = d3
    .zoom()
    .scaleExtent([0.3, 8])
    .on("start", (event) => {
      // Only flip the cursor for actual drag-pans, not wheel zooms.
      if (event.sourceEvent && event.sourceEvent.type === "mousedown") {
        svg.style("cursor", "grabbing");
      }
    })
    .on("zoom", (event) => {
      zoomLayer.attr("transform", event.transform);
      state.lastTreeTransform = event.transform;
    })
    .on("end", () => {
      svg.style("cursor", "grab");
    });
  const previousTransform = state.lastTreeTransform || d3.zoomIdentity;
  state.lastTreeSvg = svg;
  state.lastTreeZoom = zoom;
  state.lastZoomLayer = zoomLayer;

  function hasOriginalChildren(d) {
    return !!(
      d &&
      d.data &&
      Array.isArray(d.data.children) &&
      d.data.children.length > 0
    );
  }
  function isNodeVisible(d) {
    if (!d) return false;
    var cur = d.parent;
    while (cur) {
      if (cur.data && cur.data._collapsed) return false;
      cur = cur.parent;
    }
    return true;
  }
  function hasVisibleChildren(d) {
    if (!d || !d.children || !d.children.length) return false;
    if (d.data && d.data._collapsed) return false;
    for (var i = 0; i < d.children.length; i++) {
      if (isNodeVisible(d.children[i])) return true;
    }
    return false;
  }
  function isVisibleLeaf(d) {
    return !hasVisibleChildren(d);
  }
  function isHighlighted(d) {
    if (isSearchMatch(d.data)) return true;
    return nodeCustom.has(d.data) && nodeCustom.get(d.data).highlighted;
  }

  function isCollapsedInternal(d) {
    return !!(d.data && d.data._collapsed && d.data.children && d.data.children.length);
  }
  function wedgeHalfH(d) {
    var n = countLeafDescendantsData(d.data);
    return Math.max(8, Math.min(n * state.render.nodes.leafRadius * 1.2, 80)) * state.render.scale.height;
  }
  function wedgeWidth(d) {
    var n = countLeafDescendantsData(d.data);
    return Math.max(24, Math.min(n * state.render.nodes.leafRadius * 1.5, 100)) * state.render.scale.width;
  }
  function wedgeFill(d) {
    var rep = representativeClusterIdFromData(d.data);
    if (rep != null && Number(rep) >= 0 && state.clusterColours.length)
      return state.clusterColours[Number(rep) % state.clusterColours.length];
    return tc.internal;
  }
  function wedgePath(d) {
    var h = wedgeHalfH(d), w = wedgeWidth(d);
    return "M0,0 L" + w + "," + (-h) + " L" + w + "," + h + " Z";
  }

  function clusterColour(cid) {
    if (cid == null || !state.clusterColours.length) return null;
    if (Number(cid) < 0) return tc.internal;
    return state.clusterColours[cid % state.clusterColours.length];
  }

  function nodeFill(d) {
    const name = d.data && d.data.name;
    if (colourMode === "bars") return tc.label;
    if (name && !hasOriginalChildren(d)) {
      const cid = state.currentClusters ? state.currentClusters[name] : null;
      return clusterColour(cid) || tc.label;
    }
    const rep = representativeClusterIdFromData(d.data);
    return clusterColour(rep) || tc.label;
  }

  function nodeClusterId(d) {
    const name = d.data && d.data.name;
    if (name && !hasOriginalChildren(d)) {
      const cid = state.currentClusters ? state.currentClusters[name] : null;
      return cid == null ? null : Number(cid);
    }
    const rep = representativeClusterIdFromData(d.data);
    return rep == null ? null : Number(rep);
  }

  function clusterVisible(d) {
    if (selectedClusterIds.size === 0) return true;
    const cid = nodeClusterId(d);
    return cid != null && selectedClusterIds.has(cid);
  }

  // Box grouping key. Outliers share id -1, but each original small cluster
  // gets its own box so scattered outliers are not wrapped in one large box.
  function boxGroupKey(d, cid) {
    if (Number(cid) >= 0) return String(cid);
    const name = d.data && d.data.name;
    const rawCid = name != null && state.runClusters ? state.runClusters[name] : null;
    return rawCid != null && Number(rawCid) >= 0 ? "outlier:" + rawCid : String(cid);
  }

  // Branch colouring: colour edges internal to a cluster subtree.
  // The MRCA entry branch (parent outside cluster → MRCA) keeps the default colour.
  function linkStroke(link) {
    if (!state.render.branches.colorByClusters || !state.clusterColours.length)
      return state.render.branches.color || tc.branch;
    const repTarget = representativeClusterIdFromData(link.target.data);
    if (repTarget == null || Number(repTarget) < 0) return state.render.branches.color || tc.branch;
    const repSource = representativeClusterIdFromData(link.source.data);
    if (repSource !== repTarget) return state.render.branches.color || tc.branch;
    return state.clusterColours[Number(repTarget) % state.clusterColours.length];
  }

  function linkOpacity(link) {
    if (!state.render.branches.colorByClusters || selectedClusterIds.size === 0) return 1;
    const rep = representativeClusterIdFromData(link.target.data);
    if (rep == null) return 0.18;
    return selectedClusterIds.has(Number(rep)) ? 1 : 0.18;
  }

  // Wire hover tooltip on a node <g> selection. Shares local closures so it
  // can call nodeClusterId and nodeDisplayName directly.
  function wireTooltip(nodeSelection) {
    nodeSelection
      .on("mouseenter.tooltip", function (event, d) {
        var tip = document.getElementById("node-tooltip");
        if (!tip) return;
        var name = nodeDisplayName(d.data);
        var isLeafNode = !(d.data && Array.isArray(d.data.children) && d.data.children.length);
        var cid = nodeClusterId(d);
        var bl = d.data && d.data.length != null ? Number(d.data.length) : null;
        var lines = [];
        lines.push(
          '<div class="tt-name">' +
            (name
              ? escapeHtmlAttr(name)
              : isLeafNode
              ? "<em>unnamed leaf</em>"
              : "<em>internal node</em>") +
            "</div>",
        );
        if (cid !== null) {
          var cLabel = cid < 0 ? "Outlier" : "C" + cid;
          var sz = 0;
          if (state.currentClusters) {
            var cidStr = String(cid);
            var vals = Object.values(state.currentClusters);
            for (var vi = 0; vi < vals.length; vi++) {
              if (String(vals[vi]) === cidStr) sz++;
            }
          }
          lines.push(
            '<div class="tt-row"><span class="tt-key">Cluster</span><span>' +
              cLabel +
              (sz > 0 ? " &middot; " + sz + " leaves" : "") +
              "</span></div>",
          );
        }
        if (bl !== null) {
          lines.push(
            '<div class="tt-row"><span class="tt-key">Branch</span><span>' +
              bl.toFixed(5) +
              "</span></div>",
          );
        }
        tip.innerHTML = lines.join("");
        tip.removeAttribute("hidden");
        positionTooltip(tip, event.clientX, event.clientY);
      })
      .on("mousemove.tooltip", function (event) {
        var tip = document.getElementById("node-tooltip");
        if (tip && !tip.hasAttribute("hidden"))
          positionTooltip(tip, event.clientX, event.clientY);
      })
      .on("mouseleave.tooltip", function () {
        var tip = document.getElementById("node-tooltip");
        if (tip) tip.setAttribute("hidden", "");
      });
  }

  // ── CIRCULAR LAYOUT ──
  if (layoutMode === "circular") {
    const root = state.circularHierarchy;
    const visibleNodes = root.descendants().filter((d) => isNodeVisible(d));
    const visibleLinks = root
      .links()
      .filter((l) => isNodeVisible(l.source) && isNodeVisible(l.target));
    const maxRadiusPx = Math.min(width, height) / 2 - 30;
    const maxR = d3.max(root.descendants(), (d) => d._radius || 0) || 1;
    const rScale = d3.scaleLinear().domain([0, maxR]).range([0, maxRadiusPx]);

    function polarToXY(angle, radius) {
      const a = angle - Math.PI / 2;
      const rr = rScale(radius);
      return [
        rr * Math.cos(a) * state.render.scale.width,
        rr * Math.sin(a) * state.render.scale.height,
      ];
    }
    function nodeXY(d) {
      return polarToXY(d._angle || 0, d._radius || 0);
    }

    function radialElbowPath(link) {
      const s = link.source,
        t = link.target;
      const sa = s._angle || 0,
        sr = s._radius || 0;
      const ta = t._angle || 0,
        tr = t._radius || 0;
      const p0 = polarToXY(sa, sr),
        p1 = polarToXY(ta, sr),
        p2 = polarToXY(ta, tr);
      const delta = ta - sa,
        sweep = delta >= 0 ? 1 : 0;
      const largeArc = Math.abs(delta) > Math.PI ? 1 : 0;
      const arcR = rScale(sr);
      return `M${p0[0]},${p0[1]}A${arcR},${arcR} 0 ${largeArc},${sweep} ${p1[0]},${p1[1]}L${p2[0]},${p2[1]}`;
    }

    g.append("g")
      .selectAll(".tree-link")
      .data(visibleLinks)
      .enter()
      .append("path")
      .attr("class", "tree-link")
      .attr("fill", "none")
      .attr("stroke", linkStroke)
      .attr("stroke-opacity", linkOpacity)
      .attr("stroke-width", state.render.branches.width)
      .attr("d", radialElbowPath);

    const node = g
      .append("g")
      .selectAll(".tree-node")
      .data(visibleNodes)
      .enter()
      .append("g")
      .attr(
        "class",
        (d) =>
          "tree-node" + (isSearchMatch(d.data) ? " node-search-match" : ""),
      )
      .attr("transform", (d) => {
        const p = nodeXY(d);
        return `translate(${p[0]},${p[1]})`;
      });

    node
      .append("circle")
      .attr("r", (d) => isCollapsedInternal(d) ? 0 : customNodeRadius(d))
      .attr("fill", nodeFill)
      .attr("stroke", (d) => (isHighlighted(d) ? HIGHLIGHT_COLOR : "none"))
      .attr("stroke-width", (d) => (isHighlighted(d) ? 2.5 : 0))
      .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.22))
      .style("cursor", "pointer")
      .on("click", function (event, d) {
        if (!hasOriginalChildren(d)) return;
        toggleNodeCollapsedState(d.data);
        drawTree();
      })
      .on("contextmenu", showNodeContextMenu);

    node
      .filter(isCollapsedInternal)
      .append("path")
      .attr("class", "collapse-wedge")
      .attr("d", wedgePath)
      .attr("fill", wedgeFill)
      .attr("opacity", (d) => clusterVisible(d) ? 0.72 : 0.22)
      .attr("transform", (d) => "rotate(" + (((d._angle || 0) * 180) / Math.PI - 90) + ")")
      .style("cursor", "pointer")
      .on("click", function (event, d) {
        toggleNodeCollapsedState(d.data);
        drawTree();
      })
      .on("contextmenu", showNodeContextMenu);

    node
      .filter(isCollapsedInternal)
      .append("text")
      .attr("dy", "0.32em")
      .attr("font-size", Math.max(8, state.render.labels.fontSize - 2) + "px")
      .attr("font-weight", "700")
      .attr("fill", "white")
      .attr("pointer-events", "none")
      .attr("transform", (d) => {
        var deg = ((d._angle || 0) * 180) / Math.PI - 90;
        return "rotate(" + deg + ") translate(" + (wedgeWidth(d) * 0.5) + ",0)";
      })
      .style("text-anchor", "middle")
      .attr("opacity", (d) => clusterVisible(d) ? 0.95 : 0.35)
      .text((d) => String(countLeafDescendantsData(d.data)));

    node
      .filter((d) => isVisibleLeaf(d))
      .append("text")
      .attr("dy", "0.32em")
      .attr("font-size", state.render.labels.fontSize + "px")
      .attr("fill", (d) => (isHighlighted(d) ? HIGHLIGHT_COLOR : tc.label))
      .attr("font-weight", (d) => (isHighlighted(d) ? "bold" : "normal"))
      .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
      .attr("transform", function (d) {
        const a = d._angle || 0;
        let deg = (a * 180) / Math.PI - 90;
        const left = deg > 90 || deg < -90;
        return `rotate(${deg}) translate(${customNodeRadius(d) + LABEL_PAD},0) rotate(${left ? 180 : 0})`;
      })
      .style("text-anchor", (d) => {
        const deg = ((d._angle || 0) * 180) / Math.PI - 90;
        return deg > 90 || deg < -90 ? "end" : "start";
      })
      .text((d) => (state.render.labels.show ? nodeDisplayName(d.data) : ""));

    node
      .filter((d) => !isVisibleLeaf(d) && !isCollapsedInternal(d))
      .append("text")
      .attr("dy", "0.32em")
      .attr("font-size", state.render.labels.fontSize - 1 + "px")
      .attr("fill", (d) => (isHighlighted(d) ? HIGHLIGHT_COLOR : tc.internal))
      .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
      .attr("transform", (d) => {
        const a = d._angle || 0;
        let deg = (a * 180) / Math.PI - 90;
        return `rotate(${deg}) translate(0,${-(customNodeRadius(d) + 5)})`;
      })
      .style("text-anchor", "middle")
      .text((d) => (state.render.labels.internalShow ? nodeDisplayName(d.data) : ""));

    // Circular cluster visualization
    if (state.currentClusters && state.clusterColours.length) {
      const visLeaves = root
        .descendants()
        .filter((d) => isNodeVisible(d))
        .filter((d) => isVisibleLeaf(d))
        .sort((a, b) => (a._angle || 0) - (b._angle || 0));

      function visLeafCid(d) {
        const nm = d.data && d.data.name;
        const isOrigLeaf = nm && !(d.data.children && d.data.children.length);
        if (isOrigLeaf) return state.currentClusters[nm];
        return representativeClusterIdFromData(d.data);
      }

      if (colourMode === "bars" && visLeaves.length >= 2) {
        const boundaries = midpointBoundaries(visLeaves.map((d) => d._angle || 0));

        const leafLabelChars = visLeaves.reduce(
          (mx, d) => Math.max(mx, (d.data.name || "").length),
          0,
        );
        const estRadialLabelWidth = state.render.labels.show
          ? leafLabelChars * state.render.labels.fontSize * 0.55 +
            LABEL_PAD +
            state.render.nodes.leafRadius
          : 10;
        const ringInner = maxRadiusPx + estRadialLabelWidth + 6;
        const ringOuter = ringInner + 14;
        const arc = d3.arc().innerRadius(ringInner).outerRadius(ringOuter);

        for (const run of clusterRuns(visLeaves.length, (i) => visLeafCid(visLeaves[i]))) {
          var arcColour = clusterColour(run.cid) || tc.internal;
          var isOutlierArc = Number(run.cid) < 0;
          var arcEl = g
            .append("path")
            .attr(
              "d",
              arc.startAngle(boundaries[run.start]).endAngle(boundaries[run.end]),
            )
            .attr("fill", isOutlierArc ? "none" : arcColour)
            .attr("stroke", isOutlierArc ? arcColour : "none")
            .attr("stroke-width", isOutlierArc ? 1.5 : 0)
            .attr("opacity", hasClusterFocus(run.cid) ? 0.75 : 0.18);
          if (isOutlierArc) arcEl.attr("stroke-dasharray", "4,2");
        }
        g.selectAll(".tree-link, .tree-node").raise();
      }

      // Circular boxes: wedge-shaped arcs from MRCA radius to leaf tips
      if (colourMode === "boxes" && visLeaves.length >= 2) {
        var circClusterGroups = {};
        visLeaves.forEach(function (d) {
          var cid = visLeafCid(d);
          if (cid == null) return;
          if (!state.render.clusters.showOutlierBoxes && Number(cid) < 0) return;
          var key = boxGroupKey(d, cid);
          if (!circClusterGroups[key]) circClusterGroups[key] = { cid: cid, leaves: [] };
          circClusterGroups[key].leaves.push(d);
        });
        var boxesG = g.append("g").attr("class", "cluster-boxes");
        for (var circGroup of Object.values(circClusterGroups)) {
          var cid = circGroup.cid;
          var cLeaves = circGroup.leaves;
          var mrca = findMRCA(cLeaves);
          if (!mrca) continue;
          var leafAngles = cLeaves
            .map(function (d) {
              return d._angle || 0;
            })
            .sort(function (a, b) {
              return a - b;
            });
          var mrcaR = mrca._radius || 0;
          var maxLeafR =
            d3.max(cLeaves, function (d) {
              return d._radius || 0;
            }) || 0;

          var angPad = state.render.clusters.boxPadV * 0.01;
          var minAngle = leafAngles[0] - angPad;
          var maxAngle = leafAngles[leafAngles.length - 1] + angPad;
          if (cLeaves.length === 1) {
            var spread = 0.03;
            minAngle -= spread;
            maxAngle += spread;
          }

          var innerR = rScale(mrcaR) - state.render.clusters.boxPadH;
          var outerR =
            rScale(maxLeafR) + state.render.clusters.boxPadH + (state.render.labels.show ? 40 : 5);
          if (innerR < 0) innerR = 0;

          var isOutlierWedge = Number(cid) < 0;
          var wedgeColour = isOutlierWedge
            ? tc.internal
            : clusterColour(cid) || tc.internal;
          var wedgeArc = d3
            .arc()
            .innerRadius(innerR * state.render.scale.width)
            .outerRadius(outerR * state.render.scale.width)
            .startAngle(minAngle)
            .endAngle(maxAngle)
            .cornerRadius(state.render.clusters.boxCornerRadius);

          var wedge = boxesG
            .append("path")
            .attr("d", wedgeArc)
            .attr("fill", isOutlierWedge ? "none" : wedgeColour)
            .attr("opacity", hasClusterFocus(cid) ? state.render.clusters.boxAlpha : 0.02)
            .attr("stroke", wedgeColour)
            .attr("stroke-width", 1.5)
            .attr(
              "stroke-opacity",
              hasClusterFocus(cid) ? Math.min(state.render.clusters.boxAlpha * 3.5, 0.8) : 0.08,
            );
          if (isOutlierWedge) wedge.attr("stroke-dasharray", "5,3");

          if (state.render.clusters.showBoxLabels) {
            var midAngle = (minAngle + maxAngle) / 2;
            var labelR = outerR * state.render.scale.width + 4;
            var lx = labelR * Math.cos(midAngle - Math.PI / 2);
            var ly = labelR * Math.sin(midAngle - Math.PI / 2);
            var deg = (midAngle * 180) / Math.PI - 90;
            var flip = deg > 90 || deg < -90;
            boxesG
              .append("text")
              .attr("x", lx)
              .attr("y", ly)
              .attr("dy", "0.35em")
              .attr("font-size", state.render.clusters.labelFontSize + "px")
              .attr("font-weight", "600")
              .attr("fill", wedgeColour)
              .attr("font-style", isOutlierWedge ? "italic" : "normal")
              .attr("opacity", hasClusterFocus(cid) ? 1 : 0.25)
              .attr(
                "transform",
                "rotate(" +
                  (flip ? deg + 180 : deg) +
                  "," +
                  lx +
                  "," +
                  ly +
                  ")",
              )
              .style("text-anchor", flip ? "end" : "start")
              .text(
                boxLabelMap[cid] || (isOutlierWedge ? "outlier" : "C" + cid),
              );
          }
        }
        boxesG.attr("pointer-events", "none");
        boxesG.lower();
      }
    }
    wireTooltip(node);
    applyTreeZoomBounds(svg, zoom, zoomLayer, previousTransform);
    return;
  }

  // ── CARTESIAN (rectangular / cladogram) ──
  const root = state.cartesianHierarchy;
  const allNodes = root.descendants().filter((d) => isNodeVisible(d));
  const allLinks = root
    .links()
    .filter((l) => isNodeVisible(l.source) && isNodeVisible(l.target));
  const maxYCart = d3.max(allNodes, (d) => d._y || 0) || 0;
  const leafNames = allNodes
    .filter((d) => !d.children || !d.children.length)
    .map((d) => (d.data.name || "").length);
  const maxLabelChars = d3.max(leafNames) || 0;
  const estLabelWidth = state.render.labels.show
    ? maxLabelChars * state.render.labels.fontSize * 0.6 + LABEL_PAD + state.render.nodes.leafRadius
    : 20;
  const labelColumnX = maxYCart * state.render.scale.width + estLabelWidth + 12;

  g.append("g")
    .selectAll(".tree-link")
    .data(allLinks)
    .enter()
    .append("path")
    .attr("class", "tree-link")
    .attr("fill", "none")
    .attr("stroke", linkStroke)
    .attr("stroke-opacity", linkOpacity)
    .attr("stroke-width", state.render.branches.width)
    .attr(
      "d",
      (d) =>
        "M" +
        d.source._y * state.render.scale.width +
        "," +
        d.source._x * state.render.scale.height +
        "V" +
        d.target._x * state.render.scale.height +
        "H" +
        d.target._y * state.render.scale.width,
    );

  const node = g
    .append("g")
    .selectAll(".tree-node")
    .data(allNodes)
    .enter()
    .append("g")
    .attr(
      "class",
      (d) => "tree-node" + (isSearchMatch(d.data) ? " node-search-match" : ""),
    )
    .attr(
      "transform",
      (d) =>
        `translate(${(d._y || 0) * state.render.scale.width},${(d._x || 0) * state.render.scale.height})`,
    );

  node
    .append("circle")
    .attr("r", (d) => isCollapsedInternal(d) ? 0 : customNodeRadius(d))
    .attr("fill", nodeFill)
    .attr("stroke", (d) => (isHighlighted(d) ? HIGHLIGHT_COLOR : "none"))
    .attr("stroke-width", (d) => (isHighlighted(d) ? 2.5 : 0))
    .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.22))
    .style("cursor", "pointer")
    .on("click", function (event, d) {
      if (!hasOriginalChildren(d)) return;
      toggleNodeCollapsedState(d.data);
      drawTree();
    })
    .on("contextmenu", showNodeContextMenu);

  // Collapsed-node wedge (filled triangle pointing toward leaves)
  node
    .filter(isCollapsedInternal)
    .append("path")
    .attr("class", "collapse-wedge")
    .attr("d", wedgePath)
    .attr("fill", wedgeFill)
    .attr("opacity", (d) => clusterVisible(d) ? 0.72 : 0.22)
    .style("cursor", "pointer")
    .on("click", function (event, d) {
      toggleNodeCollapsedState(d.data);
      drawTree();
    })
    .on("contextmenu", showNodeContextMenu);

  node
    .filter(isCollapsedInternal)
    .append("text")
    .attr("dy", "0.32em")
    .attr("x", (d) => wedgeWidth(d) * 0.5)
    .style("text-anchor", "middle")
    .style("font-size", Math.max(8, state.render.labels.fontSize - 2) + "px")
    .attr("font-weight", "700")
    .attr("fill", "white")
    .attr("pointer-events", "none")
    .attr("opacity", (d) => clusterVisible(d) ? 0.95 : 0.35)
    .text((d) => String(countLeafDescendantsData(d.data)));

  // Leaf labels
  node
    .filter((d) => !(d.children && d.children.length))
    .append("text")
    .attr("dy", 3)
    .attr("x", (d) => customNodeRadius(d) + LABEL_PAD)
    .style("text-anchor", "start")
    .style("font-size", state.render.labels.fontSize + "px")
    .attr("fill", (d) => (isHighlighted(d) ? HIGHLIGHT_COLOR : tc.label))
    .attr("font-weight", (d) => (isHighlighted(d) ? "bold" : "normal"))
    .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
    .text((d) => (state.render.labels.show ? nodeDisplayName(d.data) : ""));

  // Internal node labels (not shown on collapsed nodes — wedge communicates that)
  node
    .filter((d) => !!(d.children && d.children.length) && !isCollapsedInternal(d))
    .append("text")
    .attr("dy", (d) => customNodeRadius(d) + 12)
    .attr("x", (d) => customNodeRadius(d) + 3)
    .style("text-anchor", "start")
    .style("font-size", state.render.labels.fontSize - 1 + "px")
    .style("font-style", "italic")
    .attr("fill", (d) => (isHighlighted(d) ? HIGHLIGHT_COLOR : tc.internal))
    .attr("font-weight", (d) => (isHighlighted(d) ? "bold" : "normal"))
    .attr("opacity", (d) => (clusterVisible(d) ? 1 : 0.3))
    .text((d) => (state.render.labels.internalShow ? nodeDisplayName(d.data) : ""));

  // Cluster side bars
  if (colourMode === "bars") {
    const leafNodes = allNodes
      .filter((d) => isVisibleLeaf(d))
      .sort((a, b) => (a._x || 0) - (b._x || 0));
    const barX = labelColumnX,
      barW = 20;
    const ys = leafNodes.map((d) => (d._x || 0) * state.render.scale.height);

    if (leafNodes.length === 1) {
      const cid = nodeClusterId(leafNodes[0]);
      if (cid != null) {
        var barColour = clusterColour(cid) || tc.internal;
        var barRect = g
          .append("rect")
          .attr("x", barX)
          .attr("y", ys[0] - 6)
          .attr("width", barW)
          .attr("height", 12)
          .attr("fill", Number(cid) < 0 ? "none" : barColour)
          .attr("stroke", Number(cid) < 0 ? barColour : "none")
          .attr("stroke-width", Number(cid) < 0 ? 1.5 : 0)
          .attr("opacity", 0.75);
        if (Number(cid) < 0) barRect.attr("stroke-dasharray", "4,2");
      }
    } else if (leafNodes.length > 1) {
      const boundaries = midpointBoundaries(ys);
      const barsG = g.append("g").attr("class", "cluster-bars");
      for (const run of clusterRuns(leafNodes.length, (i) => nodeClusterId(leafNodes[i]))) {
        const y0 = boundaries[run.start],
          y1 = boundaries[run.end];
        var runColour = clusterColour(run.cid) || tc.internal;
        var isOutlierBar = Number(run.cid) < 0;
        var barEl = barsG
          .append("rect")
          .attr("x", barX)
          .attr("y", y0)
          .attr("width", barW)
          .attr("height", Math.max(1, y1 - y0))
          .attr("fill", isOutlierBar ? "none" : runColour)
          .attr("stroke", isOutlierBar ? runColour : "none")
          .attr("stroke-width", isOutlierBar ? 1.5 : 0)
          .attr("opacity", hasClusterFocus(run.cid) ? 0.75 : 0.18);
        if (isOutlierBar) barEl.attr("stroke-dasharray", "4,2");
      }
      barsG.attr("pointer-events", "none");
    }
  }

  // Cluster boxes (MRCA mode)
  if (colourMode === "boxes" && state.currentClusters && state.clusterColours.length) {
    const clusterGroups = {};
    allNodes.forEach((d) => {
      if (!isVisibleLeaf(d)) return;
      const cid = nodeClusterId(d);
      if (cid == null) return;
      if (!state.render.clusters.showOutlierBoxes && Number(cid) < 0) return;
      const key = boxGroupKey(d, cid);
      if (!clusterGroups[key]) clusterGroups[key] = { cid, leaves: [] };
      clusterGroups[key].leaves.push(d);
    });
    var estLabelPx = state.render.labels.show
      ? maxLabelChars * state.render.labels.fontSize * 0.6 +
        LABEL_PAD +
        state.render.nodes.leafRadius +
        16
      : 16;
    var boxPadV = state.render.clusters.boxPadV;
    var boxPadL = state.render.clusters.boxPadH;
    var boxesG = g.append("g").attr("class", "cluster-boxes");
    for (var group of Object.values(clusterGroups)) {
      var cid = group.cid;
      var leaves = group.leaves;
      var mrca = findMRCA(leaves);
      if (!mrca) continue;
      var ys = leaves.map(function (d) {
        return (d._x || 0) * state.render.scale.height;
      });
      var minY = Math.min.apply(null, ys) - boxPadV;
      var maxY = Math.max.apply(null, ys) + boxPadV;
      var mrcaX = (mrca._y || 0) * state.render.scale.width - boxPadL;
      var maxLeafX = Math.max.apply(
        null,
        leaves.map(function (d) {
          return (d._y || 0) * state.render.scale.width;
        }),
      );
      var boxRight = maxLeafX + estLabelPx;
      var adj = getBoxAdjust(cid);
      minY += adj.dy - adj.padY;
      maxY += adj.dy + adj.padY;
      mrcaX += adj.dx - adj.padX;
      boxRight += adj.dx + adj.padX;
      var boxW = boxRight - mrcaX;
      var isOutlier = Number(cid) < 0;
      var colour = isOutlier
        ? tc.internal
        : state.clusterColours[cid % state.clusterColours.length];
      var rect = boxesG
        .append("rect")
        .attr("x", mrcaX)
        .attr("y", minY)
        .attr("width", Math.max(boxW, 8))
        .attr("height", Math.max(maxY - minY, 4))
        .attr("rx", state.render.clusters.boxCornerRadius)
        .attr("ry", state.render.clusters.boxCornerRadius)
        .attr("fill", isOutlier ? "none" : colour)
        .attr("opacity", hasClusterFocus(cid) ? state.render.clusters.boxAlpha : 0.02)
        .attr("stroke", colour)
        .attr("stroke-width", 1.5)
        .attr(
          "stroke-opacity",
          hasClusterFocus(cid) ? Math.min(state.render.clusters.boxAlpha * 3.5, 0.8) : 0.08,
        );
      if (isOutlier) {
        rect.attr("stroke-dasharray", "5,3");
      }
      if (state.render.clusters.showBoxLabels) {
        boxesG
          .append("text")
          .attr("x", boxRight + 4)
          .attr("y", (minY + maxY) / 2)
          .attr("dy", "0.35em")
          .attr("font-size", state.render.clusters.labelFontSize + "px")
          .attr("font-weight", "600")
          .attr("fill", colour)
          .attr("font-style", isOutlier ? "italic" : "normal")
          .attr("opacity", hasClusterFocus(cid) ? 1 : 0.25)
          .text(boxLabelMap[cid] || (isOutlier ? "outlier" : "C" + cid));
      }
    }
    boxesG.attr("pointer-events", "none");
    boxesG.lower();
  }

  // Branch-length axis — not meaningful for cladogram (equal-depth layout)
  if (layoutMode !== "cladogram") {
    var axisMaxX = d3.max(allNodes, (d) => d._x || 0) || 0;
    var axisMaxY = d3.max(allNodes, (d) => d._y || 0) || 0;
    var axisMaxBl = d3.max(allNodes, (d) => d.data._bl || 0) || 1;
    var axisYPos = axisMaxX * state.render.scale.height + 20;
    var blScale = d3
      .scaleLinear()
      .domain([0, axisMaxBl])
      .range([0, axisMaxY * state.render.scale.width]);

    var axisG = g
      .append("g")
      .attr("class", "branch-length-axis")
      .attr("transform", "translate(0, " + axisYPos + ")")
      .call(d3.axisBottom(blScale).ticks(5));
    axisG
      .selectAll("text")
      .attr("fill", tc.internal)
      .style("font-size", state.render.labels.axisFontSize + "px");
    axisG.selectAll("line").attr("stroke", tc.branch);
    axisG.selectAll("path").attr("stroke", tc.branch);
    axisG
      .append("text")
      .attr("x", (axisMaxY * state.render.scale.width) / 2)
      .attr("y", state.render.labels.axisFontSize + 20)
      .attr("text-anchor", "middle")
      .attr("font-size", state.render.labels.axisFontSize)
      .attr("fill", tc.internal)
      .text("Branch length");

    svg.attr("height", Math.max(height, axisYPos + 50 + margin.bottom));
  }

  wireTooltip(node);
  applyTreeZoomBounds(svg, zoom, zoomLayer, previousTransform);
}

/** Zoom `svg` so the content of `layer` fills the viewport, leaving `insetRight` pixels free on the right; returns the transform applied, or null. */
export function fitToViewport(svg, zoom, layer, insetRight = 0) {
  if (!svg || !zoom || !layer) return null;

  let bbox;
  try { bbox = layer.node().getBBox(); } catch { return null; }
  if (!bbox || !isFinite(bbox.width) || bbox.width <= 0 || bbox.height <= 0) return null;

  const svgNode = svg.node();
  const vw = ((svgNode && svgNode.clientWidth) || 800) - insetRight;
  const vh = (svgNode && svgNode.clientHeight) || 600;
  const pad = 28;

  const scale = Math.min(
    (vw - pad * 2) / bbox.width,
    (vh - pad * 2) / bbox.height,
    4,
  );
  const tx = pad + (vw - pad * 2 - bbox.width * scale) / 2 - bbox.x * scale;
  const ty = pad + (vh - pad * 2 - bbox.height * scale) / 2 - bbox.y * scale;

  const t = d3.zoomIdentity.translate(tx, ty).scale(scale);
  svg.transition().duration(280).call(zoom.transform, t);
  return t;
}

export function fitTree() {
  const layout = state.render.layout || "rectangular";
  const inset = layout === "circular" ? 0 : scorePanelReserve();
  const t = fitToViewport(state.lastTreeSvg, state.lastTreeZoom, state.lastZoomLayer, inset);
  if (t) state.lastTreeTransform = t;
}
