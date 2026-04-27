/* ============================================================
   PhytClust – views/compare.js
   Comparison mode: multi-bar tree comparisons.
   ============================================================ */

import { state } from "../state.js";
import { showToast } from "../ui/toast.js";
import { escapeHtmlAttr } from "../utils.js";
import { drawComparisonBarsInto } from "../tree/draw.js";
import { generateClusterColors, getThemeColors } from "../colors.js";

function listTreeLeaves() {
  const out = [];
  function walk(node) {
    if (!node) return;
    if (node.children && node.children.length) {
      for (const c of node.children) walk(c);
      return;
    }
    if (node.name) out.push(node.name);
  }
  walk(state.NEWICK_RAW_TREE);
  return out;
}

export function normalizeClustersForLeaves(rawMap) {
  if (!rawMap) return null;
  const leaves = listTreeLeaves();
  if (!leaves.length) return null;
  const leafSet = new Set(leaves);
  const out = {};
  let found = 0;
  Object.keys(rawMap).forEach((name) => {
    if (!leafSet.has(name)) return;
    const cid = parseInt(rawMap[name], 10);
    if (isNaN(cid)) return;
    out[name] = cid;
    found += 1;
  });
  if (!found) return null;
  for (const name of leaves) {
    if (out[name] == null) out[name] = 0;
  }
  return out;
}

export function getAvailableCompareKs(data) {
  if (!data) return [];
  if (data.all_ks && data.all_ks.length) return data.all_ks.slice();
  return (data.k_values || data.ks || []).slice();
}

export function getClustersForK(k) {
  if (!state.latestApiData) return null;
  if (state.latestApiData.all_clusters && state.latestApiData.all_ks) {
    const idx = state.latestApiData.all_ks.indexOf(k);
    if (idx >= 0) return state.latestApiData.all_clusters[idx];
  }
  const peakKs = state.latestApiData.k_values || state.latestApiData.ks;
  if (peakKs && state.latestApiData.clusters) {
    const idx = peakKs.indexOf(k);
    if (idx >= 0) return state.latestApiData.clusters[idx];
  }
  return null;
}

export function parseClusterTableText(rawText) {
  const lines = String(rawText || "")
    .replace(/\r\n/g, "\n")
    .split("\n")
    .map((s) => s.trim())
    .filter(Boolean);
  if (!lines.length) return { map: null, invalidRows: 0 };

  const delimiter = lines[0].indexOf("\t") >= 0 ? "\t" : ",";
  const map = {};
  let invalidRows = 0;
  let parsedRows = 0;

  lines.forEach((line, idx) => {
    let cols = line.split(delimiter).map((s) => s.trim());
    if (cols.length < 2 && delimiter === ",") {
      cols = line.split(/\s+/).map((s) => s.trim());
    }
    if (cols.length < 2) {
      invalidRows += 1;
      return;
    }

    let leaf = cols[0];
    let cidRaw = cols[1];

    const firstNum = parseInt(cols[0], 10);
    const secondNum = parseInt(cols[1], 10);
    if (!isNaN(firstNum) && isNaN(secondNum)) {
      leaf = cols[1];
      cidRaw = cols[0];
    }

    const cid = parseInt(cidRaw, 10);
    if (isNaN(cid)) {
      if (idx === 0) return;
      invalidRows += 1;
      return;
    }
    map[leaf] = cid;
    parsedRows += 1;
  });

  return { map: parsedRows ? map : null, invalidRows };
}

export function makeCompareConfig() {
  const ks = getAvailableCompareKs(state.latestApiData);
  const k = ks.length ? Number(ks[0]) : null;
  const id = state.COMPARE_CONFIG_NEXT_ID++;
  return {
    id,
    title: "Bar " + id,
    source: "k",
    k,
    fileName: "",
    fileClusters: null,
  };
}

export function ensureCompareConfigs(data) {
  const ks = getAvailableCompareKs(data);
  if (!state.COMPARE_CONFIGS.length) {
    const cfg = makeCompareConfig();
    cfg.title = ks.length ? "k=" + ks[0] : "Bar 1";
    cfg.k = ks.length ? Number(ks[0]) : null;
    state.COMPARE_CONFIGS = [cfg];
  }

  state.COMPARE_CONFIGS.forEach((cfg) => {
    if (
      cfg.source === "k" &&
      (cfg.k == null || ks.indexOf(Number(cfg.k)) < 0)
    ) {
      cfg.k = ks.length ? Number(ks[0]) : null;
    }
    if (cfg.source === "file" && cfg.fileClusters) {
      cfg.fileClusters = normalizeClustersForLeaves(cfg.fileClusters);
    }
  });
}

export function buildContingency(mapA, mapB) {
  const rowSet = new Set();
  const colSet = new Set();
  Object.values(mapA).forEach((v) => rowSet.add(v));
  Object.values(mapB).forEach((v) => colSet.add(v));
  const rowIds = Array.from(rowSet).sort((a, b) => a - b);
  const colIds = Array.from(colSet).sort((a, b) => a - b);

  const cells = {};
  rowIds.forEach((r) => {
    cells[r] = {};
    colIds.forEach((c) => { cells[r][c] = 0; });
  });
  const rowTotals = {};
  const colTotals = {};
  rowIds.forEach((r) => { rowTotals[r] = 0; });
  colIds.forEach((c) => { colTotals[c] = 0; });

  Object.keys(mapA).forEach((leaf) => {
    const r = mapA[leaf];
    const c = mapB[leaf];
    if (r == null || c == null) return;
    if (cells[r] == null) return;
    cells[r][c] = (cells[r][c] || 0) + 1;
    rowTotals[r] = (rowTotals[r] || 0) + 1;
    colTotals[c] = (colTotals[c] || 0) + 1;
  });

  return { rowIds, colIds, cells, rowTotals, colTotals };
}

function clusterLabel(cid) {
  return Number(cid) < 0 ? "†" : "C" + cid;
}

function colorWithAlpha(color, alpha) {
  if (!color) return "rgba(128,128,128," + alpha.toFixed(3) + ")";
  if (color.startsWith("rgba")) {
    return color.replace(/,\s*[\d.]+\)$/, "," + alpha.toFixed(3) + ")");
  }
  if (color.startsWith("rgb")) {
    return color.replace("rgb(", "rgba(").replace(")", "," + alpha.toFixed(3) + ")");
  }
  const r = parseInt(color.slice(1, 3), 16);
  const g = parseInt(color.slice(3, 5), 16);
  const b = parseInt(color.slice(5, 7), 16);
  return "rgba(" + r + "," + g + "," + b + "," + alpha.toFixed(3) + ")";
}

function renderAgreementPair(container, titleText, mapA, mapB, titleA, titleB) {
  const ct = buildContingency(mapA, mapB);
  if (!ct.rowIds.length || !ct.colIds.length) return;

  const maxCidA = Math.max(-1, ...ct.rowIds.map(Number));
  const colorsA = maxCidA >= 0 ? generateClusterColors(maxCidA + 1) : [];
  const tc = getThemeColors();

  const wrap = document.createElement("div");
  wrap.className = "agree-pair";

  const title = document.createElement("div");
  title.className = "agree-pair-title";
  title.textContent = titleA + "  →  " + titleB;
  wrap.appendChild(title);

  const table = document.createElement("table");
  table.className = "agree-table";

  // Header row
  const thead = table.createTHead();
  const hrow = thead.insertRow();
  const cornerTh = document.createElement("th");
  hrow.appendChild(cornerTh);
  ct.colIds.forEach((c) => {
    const th = document.createElement("th");
    th.className = "agree-col-label";
    const nonZeroRows = ct.rowIds.filter((r) => (ct.cells[r][c] || 0) > 0).length;
    th.innerHTML = clusterLabel(c) + (nonZeroRows >= 2
      ? ' <span class="badge-merge">←merge</span>'
      : "");
    hrow.appendChild(th);
  });

  // Body rows
  const tbody = table.createTBody();
  let splitCount = 0;
  let mergeCount = 0;

  ct.rowIds.forEach((r) => {
    const row = tbody.insertRow();
    const labelTd = document.createElement("td");
    labelTd.className = "agree-row-label";
    const nonZeroCols = ct.colIds.filter((c) => (ct.cells[r][c] || 0) > 0).length;
    const isSplit = nonZeroCols >= 2;
    if (isSplit) splitCount++;
    labelTd.innerHTML = clusterLabel(r) + (isSplit
      ? ' <span class="badge-split">split→' + nonZeroCols + "</span>"
      : "");
    row.appendChild(labelTd);

    ct.colIds.forEach((c) => {
      const count = ct.cells[r][c] || 0;
      const td = document.createElement("td");
      td.className = "agree-cell" + (count === 0 ? " agree-cell-zero" : "");
      const nCidA = Number(r);
      const isOutlier = nCidA < 0;
      if (count > 0 && ct.rowTotals[r] > 0) {
        const alpha = count / ct.rowTotals[r];
        const baseColor = isOutlier ? tc.internal : (colorsA.length ? colorsA[nCidA % colorsA.length] : "#888");
        td.style.background = colorWithAlpha(baseColor, alpha * 0.8 + 0.05);
      }
      const isExact = ct.rowTotals[r] === count && ct.colTotals[c] === count && count > 0;
      td.textContent = count > 0 ? String(count) : "";
      if (isExact) {
        const mark = document.createElement("span");
        mark.className = "agree-exact";
        mark.textContent = "≡";
        td.appendChild(mark);
      }
      row.appendChild(td);
    });
  });

  // Count merges from colTotals perspective
  ct.colIds.forEach((c) => {
    const nonZeroRows = ct.rowIds.filter((r) => (ct.cells[r][c] || 0) > 0).length;
    if (nonZeroRows >= 2) mergeCount++;
  });

  table.appendChild(tbody);
  wrap.appendChild(table);

  const summary = document.createElement("div");
  summary.className = "agree-summary";
  if (splitCount === 0 && mergeCount === 0) {
    summary.textContent = "All clusters map cleanly.";
  } else {
    const parts = [];
    if (splitCount > 0) parts.push(splitCount + " cluster" + (splitCount === 1 ? "" : "s") + " split");
    if (mergeCount > 0) parts.push(mergeCount + " cluster" + (mergeCount === 1 ? "" : "s") + " merged");
    summary.textContent = parts.join(", ") + ".";
  }
  wrap.appendChild(summary);
  container.appendChild(wrap);
}

export function renderAgreementSection(comparisons) {
  const host = document.getElementById("compare-agreement");
  if (!host) return;
  host.innerHTML = "";

  if (!comparisons || comparisons.length < 2) {
    const msg = document.createElement("div");
    msg.style.cssText = "font-size:12px;color:var(--pc-text-muted);padding:8px 0;";
    msg.textContent = "Add at least two comparison bars to see the agreement matrix.";
    host.appendChild(msg);
    return;
  }

  const refIdx = Math.min(
    Math.max(0, state.COMPARE_REFERENCE_INDEX || 0),
    comparisons.length - 1,
  );
  const ref = comparisons[refIdx];
  if (!ref || !ref.clusters) return;

  for (let i = 0; i < comparisons.length; i++) {
    if (i === refIdx) continue;
    const other = comparisons[i];
    if (!other.clusters) continue;
    renderAgreementPair(host, "", ref.clusters, other.clusters, ref.title, other.title);
  }
}

export function buildComparisonTSV(comparisons) {
  const leaves = listTreeLeaves();
  if (!leaves.length || !comparisons || !comparisons.length) return "";
  const header = ["leaf"].concat(comparisons.map((c) => c.title)).join("\t");
  const rows = leaves.map((leaf) => {
    const cols = [leaf].concat(
      comparisons.map((c) => {
        const cid = c.clusters ? c.clusters[leaf] : null;
        return cid == null ? "" : String(cid);
      }),
    );
    return cols.join("\t");
  });
  return header + "\n" + rows.join("\n");
}

function computeChangedLeaves(comparisons) {
  const changed = new Set();
  if (comparisons.length < 2) return changed;
  const leaves = listTreeLeaves();
  for (const leaf of leaves) {
    let first = null;
    let differs = false;
    for (const c of comparisons) {
      const cid = c.clusters ? c.clusters[leaf] : null;
      if (cid == null) continue;
      if (first == null) {
        first = cid;
      } else if (cid !== first) {
        differs = true;
        break;
      }
    }
    if (differs) changed.add(leaf);
  }
  return changed;
}

function updateCompareReferenceSelect(comparisons) {
  const wrap = document.getElementById("compare-reference-group");
  const sel = document.getElementById("compare-reference");
  if (!wrap || !sel) return;
  if (comparisons.length < 2) {
    wrap.style.display = "none";
    return;
  }
  wrap.style.display = "";
  sel.innerHTML = "";
  comparisons.forEach((c, i) => {
    const opt = document.createElement("option");
    opt.value = String(i);
    opt.textContent = c.title;
    sel.appendChild(opt);
  });
  if (state.COMPARE_REFERENCE_INDEX >= comparisons.length) {
    state.COMPARE_REFERENCE_INDEX = 0;
  }
  sel.value = String(state.COMPARE_REFERENCE_INDEX);
}

function updateCompareSummary(comparisons, changed) {
  const el = document.getElementById("compare-summary");
  if (!el) return;
  if (comparisons.length < 2) {
    el.textContent = "";
    return;
  }
  const total = listTreeLeaves().length;
  const same = total - changed.size;
  el.textContent = same + " stable · " + changed.size + " changed";
}

function updateCompareHeader(comparisons) {
  const leftLabel = document.getElementById("compare-left-label");
  if (!leftLabel) return;
  if (!comparisons.length) {
    leftLabel.textContent =
      "Add comparison bars from k results or uploaded files.";
    return;
  }
  leftLabel.textContent = comparisons.map((c) => c.title).join(" | ");
}

export function renderCompareConfigList() {
  const host = document.getElementById("compare-config-list");
  if (!host) return;
  if (!state.latestApiData) {
    host.innerHTML = "";
    return;
  }

  const ks = getAvailableCompareKs(state.latestApiData);
  host.innerHTML = "";

  state.COMPARE_CONFIGS.forEach((cfg) => {
    const card = document.createElement("div");
    card.className = "compare-config-card";
    card.dataset.cfgId = String(cfg.id);

    const titleGroup = document.createElement("div");
    titleGroup.className = "form-group";
    titleGroup.innerHTML =
      '<label>Title</label><input class="form-input compare-title" type="text" value="' +
      escapeHtmlAttr(cfg.title || "") +
      '" />';

    const sourceGroup = document.createElement("div");
    sourceGroup.className = "form-group";
    sourceGroup.innerHTML =
      '<label>Source</label><select class="form-input compare-source"><option value="k">Run result (k)</option><option value="file">CSV/TSV file</option></select>';

    const removeWrap = document.createElement("div");
    removeWrap.style.alignSelf = "end";
    const removeBtn = document.createElement("button");
    removeBtn.className = "btn-ghost";
    removeBtn.type = "button";
    removeBtn.textContent = "Remove";
    removeWrap.appendChild(removeBtn);

    const kGroup = document.createElement("div");
    kGroup.className = "form-group full-row compare-k-wrap";
    const kLabel = document.createElement("label");
    kLabel.textContent = "k value";
    const kSelect = document.createElement("select");
    kSelect.className = "form-input compare-k";
    ks.forEach((k) => {
      const opt = document.createElement("option");
      opt.value = String(k);
      opt.textContent = "k = " + k;
      kSelect.appendChild(opt);
    });
    if (cfg.k != null) kSelect.value = String(cfg.k);
    kGroup.appendChild(kLabel);
    kGroup.appendChild(kSelect);

    const fileGroup = document.createElement("div");
    fileGroup.className = "form-group full-row compare-file-wrap";
    const fLabel = document.createElement("label");
    fLabel.textContent = "Cluster file";
    const fInput = document.createElement("input");
    fInput.type = "file";
    fInput.className = "form-input compare-file";
    fInput.accept = ".csv,.tsv,.txt";
    const fHint = document.createElement("div");
    fHint.className = "hint";
    fHint.textContent = cfg.fileName
      ? "Loaded: " + cfg.fileName
      : "Two-column file: leaf_name, cluster_id";
    fileGroup.appendChild(fLabel);
    fileGroup.appendChild(fInput);
    fileGroup.appendChild(fHint);

    card.appendChild(titleGroup);
    card.appendChild(sourceGroup);
    card.appendChild(removeWrap);
    card.appendChild(kGroup);
    card.appendChild(fileGroup);
    host.appendChild(card);

    const titleInput = card.querySelector(".compare-title");
    const sourceSel = card.querySelector(".compare-source");
    sourceSel.value = cfg.source;

    function updateVisibility() {
      const useFile = sourceSel.value === "file";
      kGroup.style.display = useFile ? "none" : "";
      fileGroup.style.display = useFile ? "" : "none";
    }
    updateVisibility();

    titleInput.addEventListener("input", function () {
      cfg.title = titleInput.value.trim() || "Bar " + cfg.id;
      drawComparison();
    });

    sourceSel.addEventListener("change", function () {
      cfg.source = sourceSel.value;
      updateVisibility();
      drawComparison();
    });

    kSelect.addEventListener("change", function () {
      cfg.k = parseInt(kSelect.value, 10);
      if (!titleInput.value.trim() || /^k=\d+$/.test(titleInput.value.trim())) {
        cfg.title = "k=" + cfg.k;
        titleInput.value = cfg.title;
      }
      drawComparison();
    });

    fInput.addEventListener("change", function () {
      const file = fInput.files && fInput.files[0] ? fInput.files[0] : null;
      if (!file) return;
      const reader = new FileReader();
      reader.onload = function (e) {
        const txt = e && e.target ? e.target.result : "";
        const parsed = parseClusterTableText(txt);
        const normalized = normalizeClustersForLeaves(parsed.map);
        if (!normalized) {
          showToast(
            "Could not parse a valid leaf/cluster mapping from file.",
            "danger",
            3200,
          );
          return;
        }
        cfg.fileClusters = normalized;
        cfg.fileName = file.name;
        if (
          !titleInput.value.trim() ||
          /^Bar \d+$/.test(titleInput.value.trim())
        ) {
          cfg.title = file.name.replace(/\.[^.]+$/, "");
          titleInput.value = cfg.title;
        }
        if (parsed.invalidRows > 0) {
          showToast(
            "Loaded with " + parsed.invalidRows + " skipped row(s).",
            "info",
            2500,
          );
        }
        renderCompareConfigList();
        drawComparison();
      };
      reader.readAsText(file);
    });

    removeBtn.addEventListener("click", function () {
      state.COMPARE_CONFIGS = state.COMPARE_CONFIGS.filter((x) => x.id !== cfg.id);
      if (!state.COMPARE_CONFIGS.length) state.COMPARE_CONFIGS = [makeCompareConfig()];
      renderCompareConfigList();
      drawComparison();
    });
  });
}

export function populateCompareSelectors(data) {
  ensureCompareConfigs(data);
  renderCompareConfigList();
}

export function drawComparison() {
  const leftHost = document.getElementById("compare-tree-left");
  const agreeHost = document.getElementById("compare-agreement");
  const headerEl = document.getElementById("compare-left-label");
  const summaryEl = document.getElementById("compare-summary");

  if (!state.latestApiData || !state.NEWICK_RAW_TREE) {
    if (leftHost) {
      leftHost.innerHTML =
        '<div class="compare-empty">' +
        '<div class="compare-empty-title">No tree loaded</div>' +
        '<div class="compare-empty-text">Load a tree and click <strong>Run PhytClust</strong> first. ' +
        'Turn on <em>Compute all k</em> to compare across multiple k values.</div>' +
        "</div>";
    }
    if (agreeHost) agreeHost.innerHTML = "";
    if (headerEl) headerEl.textContent = "Add comparison bars from k results or uploaded files.";
    if (summaryEl) summaryEl.textContent = "";
    const refWrap = document.getElementById("compare-reference-group");
    if (refWrap) refWrap.style.display = "none";
    return;
  }

  if (!leftHost) return;

  const comparisons = [];
  state.COMPARE_CONFIGS.forEach((cfg) => {
    if (cfg.source === "k") {
      const raw = getClustersForK(Number(cfg.k));
      const clusters = normalizeClustersForLeaves(raw);
      if (clusters) {
        comparisons.push({
          title: cfg.title || (cfg.k != null ? "k=" + cfg.k : "k"),
          clusters,
        });
      }
      return;
    }
    if (cfg.source === "file" && cfg.fileClusters) {
      comparisons.push({
        title: cfg.title || cfg.fileName || "file",
        clusters: cfg.fileClusters,
      });
    }
  });

  state.COMPARE_CHANGED_LEAVES = computeChangedLeaves(comparisons);
  updateCompareHeader(comparisons);
  updateCompareReferenceSelect(comparisons);
  updateCompareSummary(comparisons, state.COMPARE_CHANGED_LEAVES);
  if (!comparisons.length) {
    leftHost.innerHTML =
      '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);font-size:14px;">No comparison bars configured yet.</div>';
    renderAgreementSection(comparisons);
    return;
  }
  drawComparisonBarsInto(leftHost, comparisons);
  renderAgreementSection(comparisons);
}
