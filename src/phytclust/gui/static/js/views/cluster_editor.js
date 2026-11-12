/* ============================================================
   PhytClust – views/cluster_editor.js
   Cluster selector (peak/all switching) and bulk label editor.
   ============================================================ */

import { state, boxLabelMap } from "../state.js";
import { showToast } from "../ui/toast.js";
import { generateClusterColours } from "../colors.js";
import { drawTree } from "../tree/draw.js";
import { markOutlierClusters } from "../utils.js";

/**
 * Outlier rule for run results: the "Mark outliers" switch, and the size
 * threshold the displayed run used (singletons when it had none).
 */
export function outlierMarkingSettings() {
  const markEl = document.getElementById("mark-outliers");
  const threshold = state.lastRunData ? state.lastRunData.outlier_size_threshold : null;
  return {
    enabled: !markEl || !!markEl.checked,
    threshold: Number.isSafeInteger(threshold) ? threshold : null,
  };
}

/** Recompute the displayed cluster map (outliers as -1) from the run's map. */
export function refreshOutlierMarks() {
  state.currentClusters = markOutlierClusters(
    state.runClusters || {},
    outlierMarkingSettings(),
  );
}

/**
 * Show a run's leaf-to-cluster map. Colours follow the run's cluster ids, so
 * they do not shift when outlier marking is switched on or off.
 */
export function setCurrentClusters(clusterMap) {
  state.runClusters = clusterMap || {};
  refreshOutlierMarks();
  // A loop, not Math.max(...ids): large trees exceed the argument limit.
  let maximumId = -1;
  for (const clusterId of Object.values(state.runClusters)) {
    if (Number(clusterId) > maximumId) maximumId = Number(clusterId);
  }
  state.clusterColours = maximumId >= 0 ? generateClusterColours(maximumId + 1) : [];
}

/** Hide the result selector and cluster editor button. */
export function resetClusterSelector() {
  const controls = document.getElementById("cluster-select-controls");
  const selectEl = document.getElementById("cluster-select");
  const labelEl = document.getElementById("cluster-select-label");
  const toggleEl = document.getElementById("cluster-view-toggle");
  if (controls) controls.classList.remove("visible");
  if (selectEl) selectEl.innerHTML = "";
  if (labelEl) labelEl.textContent = "";
  if (toggleEl) {
    toggleEl.style.display = "none";
    toggleEl.textContent = "Show all k";
  }
  state.clusterViewMode = "peaks";
  updateClusterEditorAvailability();
}

export function getCurrentClusterIds() {
  if (!state.currentClusters) return [];
  const ids = new Set();
  Object.keys(state.currentClusters).forEach((leaf) => {
    const cid = Number(state.currentClusters[leaf]);
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
    const label = (boxLabelMap[cid] || defaultClusterLabel(cid))
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
    input.value = boxLabelMap[cid] || defaultClusterLabel(cid);
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
    if (!val || val === fallback) delete boxLabelMap[cid];
    else boxLabelMap[cid] = val;
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

  const savedView = state.pendingSessionView;
  if (savedView) state.clusterViewMode = savedView.cluster_mode;
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

  if (hasAll && state.clusterViewMode === "all") {
    fillClusterSelect(
      selectEl,
      labelEl,
      data.all_clusters,
      data.all_ks,
      false,
      data.all_alphas,
    );
  } else {
    if (data.clusters && data.clusters.length) {
      state.clusterViewMode = "peaks";
      fillClusterSelect(
        selectEl,
        labelEl,
        data.clusters,
        peakKs,
        true,
        data.alphas,
      );
    } else {
      state.clusterViewMode = "all";
      fillClusterSelect(
        selectEl,
        labelEl,
        data.all_clusters || [],
        data.all_ks || [],
        false,
        data.all_alphas,
      );
    }
  }
  if (savedView) {
    const counts = state.clusterViewMode === "all" ? data.all_ks : data.k_values || data.ks;
    const index = counts?.indexOf(savedView.selected_k) ?? -1;
    state.pendingSessionView = null;
    if (index >= 0) {
      selectEl.value = String(index);
      switchCluster(index);
    }
  }
  controls.classList.add("visible");
  updateClusterEditorAvailability();
}

function formatAlpha(v) {
  if (v == null || !isFinite(v)) return null;
  const abs = Math.abs(v);
  if (abs === 0) return "0";
  if (abs >= 100 || abs < 0.01) return v.toExponential(2);
  return v.toFixed(3);
}

function fillClusterSelect(selectEl, labelEl, clusters, ks, showRank, alphas) {
  selectEl.innerHTML = "";
  for (let i = 0; i < clusters.length; i++) {
    const opt = document.createElement("option");
    opt.value = i;
    const kLabel = ks && ks[i] != null ? ks[i] : i + 1;
    const alphaStr =
      alphas && alphas[i] != null ? formatAlpha(alphas[i]) : null;
    const alphaPart = alphaStr != null ? `, α = ${alphaStr}` : "";
    opt.textContent = showRank
      ? `k = ${kLabel} (rank ${i + 1}${alphaPart})`
      : `k = ${kLabel}${alphaPart}`;
    selectEl.appendChild(opt);
  }
  if (labelEl) labelEl.textContent = `1 / ${clusters.length}`;
}

export function switchCluster(idx) {
  if (!state.lastRunData) return;
  let clusters, ks;
  if (state.clusterViewMode === "all" && state.lastRunData.all_clusters) {
    clusters = state.lastRunData.all_clusters;
    ks = state.lastRunData.all_ks;
  } else {
    clusters = state.lastRunData.clusters || [];
    ks = state.lastRunData.k_values || state.lastRunData.ks || [];
  }
  if (idx < 0 || idx >= clusters.length) return;
  setCurrentClusters(clusters[idx]);
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
  if (!state.lastRunData) return;
  state.clusterViewMode = state.clusterViewMode === "peaks" ? "all" : "peaks";
  populateClusterSelector(state.lastRunData);
  switchCluster(0);
  const toggleEl = document.getElementById("cluster-view-toggle");
  if (toggleEl)
    toggleEl.textContent =
      state.clusterViewMode === "all" ? "Show peaks only" : "Show all k";
}
