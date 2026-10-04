/* Send clustering requests and update the result views. */

import { state } from "./state.js";
import { newickEl, extraOutgroupEl, extraRootTaxonEl, extraResolutionEl } from "./dom.js";
import { showToast } from "./ui/toast.js";
import { showStatus } from "./ui/status.js";
import { estimateLeafCount, parseNewick, readCheckParam, readSelectParam } from "./utils.js";
import { generateClusterColors } from "./colors.js";
import { drawTree, clearTree } from "./tree/draw.js";
import { accumulateBranchLength, computeLayouts } from "./tree/layout.js";
import { populateClusterSelector, updateClusterEditorAvailability } from "./views/cluster_editor.js";
import { populateCompareSelectors } from "./views/compare.js";
import { drawOptimalK } from "./views/optimal_k.js";
import { drawMiniScores } from "./views/scores_panel.js";

function responseError(detail, status) {
  if (typeof detail === "string" && detail.trim()) return detail;
  if (Array.isArray(detail)) {
    const messages = detail.map((issue) => {
      const field = Array.isArray(issue?.loc) ? issue.loc.filter((part) => part !== "body").join(".") : "";
      return typeof issue?.msg === "string" ? (field ? `${field}: ${issue.msg}` : issue.msg) : "";
    }).filter(Boolean);
    if (messages.length) return messages.join("; ");
  }
  return `Request failed (${status}).`;
}

async function postJson(url, payload) {
  let response;
  try {
    response = await fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
  } catch {
    throw new Error("Could not reach the server. Check that the GUI server is running.");
  }
  const text = await response.text();
  let data;
  try {
    data = JSON.parse(text);
  } catch {
    throw new Error(response.ok ? "The server returned an unreadable result." : `Request failed (${response.status}).`);
  }
  if (!response.ok) throw new Error(responseError(data?.detail, response.status));
  return data;
}

function currentMode() {
  return document.querySelector("#mode-selector .mode-btn.active")?.dataset.mode || "global";
}

function readNumber(id, label, integer = false) {
  const text = (document.getElementById(id)?.value || "").trim();
  if (!text) return null;
  const value = Number(text);
  if (!Number.isFinite(value) || (integer && (!Number.isSafeInteger(value) || value < 1))) {
    throw new Error(`${label} must be ${integer ? "a positive whole number" : "a number"}.`);
  }
  return value;
}

export function buildRunRequest() {
  const mode = currentMode();
  const payload = { newick: (newickEl.value || "").trim(), mode };
  const addNumber = (key, id, label, integer = false) => {
    const value = readNumber(id, label, integer);
    if (value !== null) payload[key] = value;
  };
  const outgroup = extraOutgroupEl?.value.trim();
  const rootTaxon = extraRootTaxonEl?.value.trim();
  if (outgroup) payload.outgroup = outgroup;
  if (rootTaxon) payload.root_taxon = rootTaxon;
  if (mode === "k") {
    payload.k = readNumber("extra-k", "Cluster count", true);
  } else {
    addNumber("top_n", "extra-topn", "Peak count", true);
    addNumber("max_k", "extra-maxk", "Maximum cluster count", true);
    addNumber("max_k_limit", "extra-maxklimit", "Maximum cluster fraction");
    const rankingMode = readSelectParam("extra-ranking-mode") || "adjusted";
    payload.ranking_mode = rankingMode;
    if (rankingMode === "adjusted") addNumber("prominence_weight", "extra-prominence-weight", "Prominence weight");
    addNumber("min_prominence", "extra-min-prominence", "Minimum prominence");
    payload.use_relative_prominence = readCheckParam("extra-relative-prom");
    payload.exclude_k2 = readCheckParam("extra-exclude-k2");
  }
  if (mode === "resolution") {
    addNumber("num_bins", "extra-bins", "Resolution class count", true);
    payload.by_resolution = true;
  }
  if (mode === "global" && readCheckParam("extra-compute-all")) payload.compute_all_clusters = true;
  addNumber("min_cluster_size", "extra-min-cluster-size", "Minimum cluster size", true);
  if (readCheckParam("extra-use-support")) {
    payload.use_branch_support = true;
    addNumber("min_support", "extra-min-support", "Minimum support");
    addNumber("support_weight", "extra-support-weight", "Support weight");
  }
  addNumber("outlier_size_threshold", "extra-outlier-threshold", "Outlier size threshold", true);
  if (readCheckParam("extra-outlier-prefer-fewer")) payload.outlier_prefer_fewer = true;
  const ratioMode = readSelectParam("extra-outlier-ratio-mode");
  if (ratioMode && ratioMode !== "exp") payload.outlier_ratio_mode = ratioMode;
  const polytomyMode = readSelectParam("extra-polytomy-mode");
  if (polytomyMode) payload.polytomy_mode = polytomyMode;
  if (readCheckParam("extra-no-split-zero")) payload.no_split_zero_length = true;
  return payload;
}

export function refreshResultsStale() {
  const marker = document.getElementById("results-stale-indicator");
  if (!marker) return;
  try {
    marker.hidden = !state.latestApiData || state.lastRunSignature === JSON.stringify(buildRunRequest());
  } catch {
    marker.hidden = !state.latestApiData;
  }
}

function prepareResult(data) {
  if (!data || typeof data !== "object" || typeof data.newick !== "string" || !data.newick.trim() || !Array.isArray(data.clusters)) {
    throw new Error("The server returned an incomplete result.");
  }
  for (const clusterMap of data.clusters) {
    if (!clusterMap || typeof clusterMap !== "object" || Array.isArray(clusterMap) || Object.values(clusterMap).some((value) => !Number.isSafeInteger(value) || value < -1)) {
      throw new Error("The server returned invalid cluster assignments.");
    }
  }
  const tree = parseNewick(data.newick);
  if (!tree || typeof tree !== "object") throw new Error("The returned tree could not be read.");
  accumulateBranchLength(tree);
  const clusterMap = data.clusters[0] || {};
  let maximumClusterId = -1;
  for (const clusterId of Object.values(clusterMap)) maximumClusterId = Math.max(maximumClusterId, clusterId);
  return { tree, clusterMap, colors: generateClusterColors(maximumClusterId + 1) };
}

function drawScoreViews(data) {
  try {
    drawMiniScores(data);
  } catch (error) {
    console.warn("Could not draw the score preview:", error);
    showToast("The score preview could not be drawn.", "warning", 3000);
  }
  const plotHost = document.getElementById("optimalk_plot");
  state.latestOptimalKData = Array.isArray(data.scores) && data.scores.length ? data : null;
  if (plotHost) plotHost.replaceChildren();
  if (!state.latestOptimalKData) {
    if (plotHost) {
      const message = document.createElement("div");
      message.className = "score-empty";
      message.textContent = data.mode === "k" ? "Scores are not calculated for a fixed cluster count." : "No scores are available for this run.";
      plotHost.appendChild(message);
    }
    return;
  }
  try {
    drawOptimalK(data);
  } catch (error) {
    console.warn("Could not draw the score plot:", error);
    showToast("The score plot could not be drawn.", "warning", 3000);
  }
}

export async function runPhytClust() {
  if (state.isRunning) return;
  let payload;
  try {
    payload = buildRunRequest();
    if (!payload.newick) throw new Error("Upload or paste a Newick tree.");
    if (payload.mode === "k" && payload.k === null) throw new Error("Enter a cluster count.");
  } catch (error) {
    showToast(error.message, "danger", 4000);
    return;
  }
  const runButton = document.getElementById("btn-run");
  const originalButtonContent = runButton?.innerHTML;
  const resultFields = ["NEWICK_RAW_TREE", "CURRENT_CLUSTERS", "CLUSTER_COLORS", "HIER_CART", "HIER_CIRC", "latestApiData", "latestRunId", "lastRunSignature", "latestOptimalKData", "CLUSTER_VIEW_MODE", "runHistory"];
  const previousState = Object.fromEntries(resultFields.map((name) => [name, state[name]]));
  let updatingViews = false;
  state.isRunning = true;
  if (extraResolutionEl) extraResolutionEl.checked = payload.mode === "resolution";
  showStatus("Running PhytClust...", "info");
  if (runButton) {
    runButton.disabled = true;
    runButton.innerHTML = '<span class="spinner"></span> Running...';
  }
  try {
    const startedAt = performance.now();
    const data = await postJson("/api/run", payload);
    const prepared = prepareResult(data);
    const elapsedSeconds = (performance.now() - startedAt) / 1000;
    updatingViews = true;
    state.NEWICK_RAW_TREE = prepared.tree;
    state.CURRENT_CLUSTERS = prepared.clusterMap;
    state.CLUSTER_COLORS = prepared.colors;
    computeLayouts();
    state.latestApiData = data;
    state.latestRunId = data.run_id || null;
    populateClusterSelector(data);
    updateClusterEditorAvailability();
    drawTree();
    state.lastRunSignature = JSON.stringify(payload);
    const leafCount = Object.keys(prepared.clusterMap).length;
    const clusterCount = new Set(Object.values(prepared.clusterMap)).size;
    const selectedCounts = data.k_values || data.ks || [];
    showStatus(leafCount ? `k = ${clusterCount}${selectedCounts.length > 1 ? ` (rank 1 of ${selectedCounts.length})` : ""} · ${leafCount} leaves · ${elapsedSeconds.toFixed(2)}s` : `No partitions selected · ${elapsedSeconds.toFixed(2)}s`, leafCount ? "success" : "info");
    state.runHistory = [{
      ts: Date.now(),
      mode: payload.mode,
      label: `${payload.mode === "k" ? `k=${payload.k}` : payload.mode} · ${elapsedSeconds.toFixed(1)}s`,
      nLeaves: estimateLeafCount(data.newick),
    }, ...state.runHistory].slice(0, 8);
    const leafCountLabel = document.getElementById("leaf-count-label");
    if (leafCountLabel) leafCountLabel.textContent = `${estimateLeafCount(data.newick)} leaves`;
    try {
      populateCompareSelectors(data);
    } catch (error) {
      console.warn("Could not update comparison choices:", error);
      showToast("Comparison choices could not be updated.", "warning", 3000);
    }
    drawScoreViews(data);
    refreshResultsStale();
  } catch (error) {
    if (updatingViews) {
      Object.assign(state, previousState, { isRunning: true });
      try {
        if (previousState.latestApiData) {
          populateClusterSelector(previousState.latestApiData);
          updateClusterEditorAvailability();
        }
        if (previousState.NEWICK_RAW_TREE) drawTree();
        else clearTree();
      } catch (restoreError) {
        console.warn("Could not redraw the previous tree:", restoreError);
      }
    }
    console.error(error);
    showStatus(`Error: ${error.message}${previousState.latestApiData ? " (previous result retained)" : ""}`, "danger");
    refreshResultsStale();
  } finally {
    state.isRunning = false;
    if (runButton) {
      runButton.disabled = false;
      runButton.innerHTML = originalButtonContent;
    }
  }
}
