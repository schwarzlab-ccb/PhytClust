/* ============================================================
   PhytClust – views/cluster_editor.js
   Cluster selector (peak/all switching) and bulk label editor.
   ============================================================ */

import { state, BOX_LABEL_MAP } from "../state.js";
import { showToast } from "../ui/toast.js";
import { generateClusterColors } from "../colors.js";
import { drawTree } from "../tree/draw.js";

export function getCurrentClusterIds() {
  if (!state.CURRENT_CLUSTERS) return [];
  const ids = new Set();
  Object.keys(state.CURRENT_CLUSTERS).forEach((leaf) => {
    const cid = Number(state.CURRENT_CLUSTERS[leaf]);
    if (isNaN(cid)) return;
    ids.add(cid);
  });
  return Array.from(ids).sort((a, b) => a - b);
}

function defaultClusterLabel(cid) {
  return Number(cid) < 0 ? "outlier" : "C" + cid;
}

export function updateClusterEditorAvailability() {
  const btn = document.getElementById("btn-edit-clusters");
  if (!btn) return;
  const hasAnyClusters = getCurrentClusterIds().length > 0;
  btn.style.display = hasAnyClusters ? "inline-flex" : "none";
}

export function renderClusterEditorRows() {
  const body = document.getElementById("cluster-editor-body");
  const filterInput = document.getElementById("cluster-filter-text");
  const showOutliersCb = document.getElementById(
    "cluster-filter-show-outliers",
  );
  if (!body) return;
  const ids = getCurrentClusterIds();
  const filter = (filterInput && filterInput.value ? filterInput.value : "")
    .toLowerCase()
    .trim();
  const includeOutliers = showOutliersCb ? showOutliersCb.checked : true;

  const visibleIds = ids.filter((cid) => {
    if (!includeOutliers && Number(cid) < 0) return false;
    if (!filter) return true;
    const label = (BOX_LABEL_MAP[cid] || defaultClusterLabel(cid))
      .toLowerCase()
      .trim();
    return String(cid).includes(filter) || label.includes(filter);
  });

  body.innerHTML = "";

  if (!visibleIds.length) {
    const tr = document.createElement("tr");
    const td = document.createElement("td");
    td.colSpan = 2;
    td.className = "cluster-editor-empty";
    td.textContent = ids.length
      ? "No clusters match this filter."
      : "No clusters available. Run clustering first.";
    tr.appendChild(td);
    body.appendChild(tr);
    return;
  }

  visibleIds.forEach((cid) => {
    const tr = document.createElement("tr");

    const tdId = document.createElement("td");
    tdId.textContent = String(cid);

    const tdLabel = document.createElement("td");
    const input = document.createElement("input");
    input.type = "text";
    input.setAttribute("data-cid", String(cid));
    input.value = BOX_LABEL_MAP[cid] || defaultClusterLabel(cid);
    tdLabel.appendChild(input);

    tr.appendChild(tdId);
    tr.appendChild(tdLabel);
    body.appendChild(tr);
  });
}

export function openClusterEditor() {
  const modal = document.getElementById("clusterEditorModal");
  const bg = document.getElementById("clusterEditorBackdrop");
  const labelSizeInput = document.getElementById("cluster-label-size-input");
  const filterInput = document.getElementById("cluster-filter-text");
  const showOutliersCb = document.getElementById(
    "cluster-filter-show-outliers",
  );
  if (filterInput) filterInput.value = "";
  if (showOutliersCb) showOutliersCb.checked = true;
  renderClusterEditorRows();
  if (labelSizeInput) {
    labelSizeInput.value = String(Math.round(state.render.clusters.labelFontSize));
  }
  if (modal) modal.classList.add("show");
  if (bg) bg.classList.add("show");
}

export function closeClusterEditor() {
  const modal = document.getElementById("clusterEditorModal");
  const bg = document.getElementById("clusterEditorBackdrop");
  if (modal) modal.classList.remove("show");
  if (bg) bg.classList.remove("show");
}

export function saveClusterEditor() {
  const body = document.getElementById("cluster-editor-body");
  const labelSizeInput = document.getElementById("cluster-label-size-input");
  if (!body) return;

  body.querySelectorAll("input[data-cid]").forEach((input) => {
    const cid = Number(input.getAttribute("data-cid"));
    const val = (input.value || "").trim();
    const fallback = defaultClusterLabel(cid);
    if (!val || val === fallback) delete BOX_LABEL_MAP[cid];
    else BOX_LABEL_MAP[cid] = val;
  });

  if (labelSizeInput) {
    const v = parseFloat(labelSizeInput.value);
    if (!isNaN(v)) {
      state.render.clusters.labelFontSize = Math.min(36, Math.max(6, v));
    }
  }

  drawTree();
  closeClusterEditor();
  showToast("Cluster labels updated.", "success", 1800);
}

export function populateClusterSelector(data) {
  const controls = document.getElementById("cluster-select-controls");
  const selectEl = document.getElementById("cluster-select");
  const labelEl = document.getElementById("cluster-select-label");
  const toggleEl = document.getElementById("cluster-view-toggle");
  if (!controls || !selectEl) return;

  const hasPeaks = (data.clusters || []).length > 1;
  const hasAll = !!(data.all_clusters && data.all_clusters.length > 1);

  if (toggleEl) toggleEl.style.display = hasAll ? "inline-flex" : "none";

  const peakKs = data.k_values || data.ks;

  const hasAny =
    (data.clusters && data.clusters.length > 0) ||
    (data.all_clusters && data.all_clusters.length > 0);
  if (!hasAny) {
    controls.classList.remove("visible");
    if (toggleEl) toggleEl.style.display = "none";
    updateClusterEditorAvailability();
    return;
  }

  if (hasAll && state.CLUSTER_VIEW_MODE === "all") {
    _fillClusterSelect(
      selectEl,
      labelEl,
      data.all_clusters,
      data.all_ks,
      false,
      data.all_alphas,
    );
  } else {
    if (data.clusters && data.clusters.length) {
      state.CLUSTER_VIEW_MODE = "peaks";
      _fillClusterSelect(
        selectEl,
        labelEl,
        data.clusters,
        peakKs,
        true,
        data.alphas,
      );
    } else {
      state.CLUSTER_VIEW_MODE = "all";
      _fillClusterSelect(
        selectEl,
        labelEl,
        data.all_clusters || [],
        data.all_ks || [],
        false,
        data.all_alphas,
      );
    }
  }
  controls.classList.add("visible");
  updateClusterEditorAvailability();
}

function _formatAlpha(v) {
  if (v == null || !isFinite(v)) return null;
  const abs = Math.abs(v);
  if (abs === 0) return "0";
  if (abs >= 100 || abs < 0.01) return v.toExponential(2);
  return v.toFixed(3);
}

function _fillClusterSelect(selectEl, labelEl, clusters, ks, showRank, alphas) {
  selectEl.innerHTML = "";
  for (let i = 0; i < clusters.length; i++) {
    const opt = document.createElement("option");
    opt.value = i;
    const kLabel = ks && ks[i] != null ? ks[i] : i + 1;
    const alphaStr =
      alphas && alphas[i] != null ? _formatAlpha(alphas[i]) : null;
    const alphaPart = alphaStr != null ? `, α = ${alphaStr}` : "";
    opt.textContent = showRank
      ? `k = ${kLabel} (rank ${i + 1}${alphaPart})`
      : `k = ${kLabel}${alphaPart}`;
    selectEl.appendChild(opt);
  }
  if (labelEl) labelEl.textContent = `1 / ${clusters.length}`;
}

export function switchCluster(idx) {
  if (!state.latestApiData) return;
  let clusters, ks;
  if (state.CLUSTER_VIEW_MODE === "all" && state.latestApiData.all_clusters) {
    clusters = state.latestApiData.all_clusters;
    ks = state.latestApiData.all_ks;
  } else {
    clusters = state.latestApiData.clusters || [];
    ks = state.latestApiData.k_values || state.latestApiData.ks || [];
  }
  if (idx < 0 || idx >= clusters.length) return;
  state.CURRENT_CLUSTERS = clusters[idx];
  if (Object.keys(state.CURRENT_CLUSTERS).length > 0) {
    state.CLUSTER_COLORS = generateClusterColors(
      Math.max(...Object.values(state.CURRENT_CLUSTERS)) + 1,
    );
  } else {
    state.CLUSTER_COLORS = [];
  }
  updateClusterEditorAvailability();
  drawTree();
  const labelEl = document.getElementById("cluster-select-label");
  if (labelEl) labelEl.textContent = `${idx + 1} / ${clusters.length}`;
}

export function cycleCluster(delta) {
  var selectEl = document.getElementById("cluster-select");
  if (!selectEl || !selectEl.options.length) return;
  var cur = selectEl.selectedIndex;
  var next = cur + delta;
  if (next < 0) next = selectEl.options.length - 1;
  if (next >= selectEl.options.length) next = 0;
  selectEl.selectedIndex = next;
  switchCluster(parseInt(selectEl.value, 10));
}

export function toggleClusterViewMode() {
  if (!state.latestApiData) return;
  state.CLUSTER_VIEW_MODE = state.CLUSTER_VIEW_MODE === "peaks" ? "all" : "peaks";
  populateClusterSelector(state.latestApiData);
  switchCluster(0);
  const toggleEl = document.getElementById("cluster-view-toggle");
  if (toggleEl)
    toggleEl.textContent =
      state.CLUSTER_VIEW_MODE === "all" ? "Show peaks only" : "Show all k";
}
