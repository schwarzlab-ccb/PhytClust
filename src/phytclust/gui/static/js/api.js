/* Send clustering requests and update the result views. */

import { state } from "./state.js";
import { newickEl, outgroupEl, rootTaxonEl } from "./dom.js";
import { showToast } from "./ui/toast.js";
import { showStatus, logStatus } from "./ui/status.js";
import { estimateLeafCount, parseNewick, readCheckParam, readSelectParam } from "./utils.js";
import { drawTree, clearTree } from "./tree/draw.js";
import { accumulateBranchLength, computeLayouts } from "./tree/layout.js";
import { populateClusterSelector, updateClusterEditorAvailability, setCurrentClusters, resetClusterSelector } from "./views/cluster_editor.js";
import { populateCompareSelectors, drawComparison } from "./views/compare.js";
import { drawScores } from "./views/scores.js";
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

export async function postJson(url, payload) {
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
  const outgroup = outgroupEl?.value.trim();
  const rootTaxon = rootTaxonEl?.value.trim();
  if (outgroup) payload.outgroup = outgroup;
  if (rootTaxon) payload.root_taxon = rootTaxon;
  if (mode === "k") {
    payload.k = readNumber("k", "Cluster count", true);
  } else {
    addNumber("top_n", "top-n", "Peak count", true);
    addNumber("max_k", "max-k", "Maximum cluster count", true);
    addNumber("max_k_limit", "max-k-limit", "Maximum cluster fraction");
    const rankingMode = readSelectParam("ranking-mode") || "adjusted";
    payload.ranking_mode = rankingMode;
    if (rankingMode === "adjusted") addNumber("prominence_weight", "prominence-weight", "Prominence weight");
    addNumber("min_prominence", "min-prominence", "Minimum prominence");
    payload.use_relative_prominence = readCheckParam("relative-prominence");
    payload.exclude_k2 = readCheckParam("exclude-k2");
  }
  if (mode === "resolution") {
    addNumber("num_bins", "bins", "Resolution class count", true);
    payload.by_resolution = true;
  }
  if (mode === "global" && readCheckParam("compute-all")) payload.compute_all_clusters = true;
  addNumber("min_cluster_size", "min-cluster-size", "Minimum cluster size", true);
  if (readCheckParam("use-support")) {
    payload.use_branch_support = true;
    addNumber("min_support", "min-support", "Minimum support");
    addNumber("support_weight", "support-weight", "Support weight");
  }
  addNumber("outlier_size_threshold", "outlier-threshold", "Outlier size threshold", true);
  if (readCheckParam("outlier-prefer-fewer")) payload.outlier_prefer_fewer = true;
  const ratioMode = readSelectParam("outlier-ratio-mode");
  if (ratioMode && ratioMode !== "exp") payload.outlier_ratio_mode = ratioMode;
  const polytomyMode = readSelectParam("polytomy-mode");
  if (polytomyMode) payload.polytomy_mode = polytomyMode;
  if (readCheckParam("no-split-zero")) payload.no_split_zero_length = true;
  return payload;
}

export function refreshResultsStale() {
  const marker = document.getElementById("results-stale-indicator");
  if (!marker) return;
  try {
    marker.hidden = !state.lastRunData || state.lastRunSignature === JSON.stringify(buildRunRequest());
  } catch {
    marker.hidden = !state.lastRunData;
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
  return { tree, clusterMap: data.clusters[0] || {} };
}

/**
 * Messages to show after a run: input repairs the server reported, and a
 * note when a multi-resolution request was answered with a global search.
 */
export function runNotes(data, requestedMode) {
  const notes = Array.isArray(data?.notes)
    ? data.notes.filter((note) => typeof note === "string" && note.trim()).map((note) => note.trim())
    : [];
  if (requestedMode === "resolution" && data?.mode && data.mode !== "resolution" &&
      !notes.some((note) => /resolution/i.test(note))) {
    notes.push("This tree is too small for a multi-resolution search; a global search was run instead.");
  }
  return notes;
}

function reportRunNotes(notes) {
  if (!notes.length) return;
  for (const note of notes) logStatus(note, "warning");
  const more = notes.length > 1 ? ` (+${notes.length - 1} more in the activity log)` : "";
  showToast(`Input adjusted: ${notes[0]}${more}`, "warning", 8000);
}

/** Forget the current result and remove every view of it. */
export function clearResults() {
  state.newickTree = null;
  state.cartesianHierarchy = null;
  state.circularHierarchy = null;
  state.lastRunData = null;
  state.lastRunId = null;
  state.lastRunSignature = null;
  state.lastScoresData = null;
  state.pendingSessionView = null;
  state.compareConfigs = [];
  state.compareChangedLeaves = new Set();
  state.lastTreeTransform = d3.zoomIdentity;
  setCurrentClusters({});
  resetClusterSelector();
  drawMiniScores(null);
  const plotHost = document.getElementById("scores-plot");
  if (plotHost) plotHost.replaceChildren();
  const leafCountLabel = document.getElementById("leaf-count-label");
  if (leafCountLabel) leafCountLabel.textContent = "";
  clearTree();
  drawComparison();
  refreshResultsStale();
}

function drawScoreViews(data) {
  try {
    drawMiniScores(data);
  } catch (error) {
    console.warn("Could not draw the score preview:", error);
    showToast("The score preview could not be drawn.", "warning", 3000);
  }
  const plotHost = document.getElementById("scores-plot");
  state.lastScoresData = Array.isArray(data.scores) && data.scores.length ? data : null;
  if (plotHost) plotHost.replaceChildren();
  if (!state.lastScoresData) {
    if (plotHost) {
      const message = document.createElement("div");
      message.className = "score-empty";
      message.textContent = data.mode === "k" ? "Scores are not calculated for a fixed cluster count." : "No scores are available for this run.";
      plotHost.appendChild(message);
    }
    return;
  }
  try {
    drawScores(data);
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
  const resultFields = ["newickTree", "currentClusters", "runClusters", "clusterColours", "cartesianHierarchy", "circularHierarchy", "lastRunData", "lastRunId", "lastRunSignature", "lastScoresData", "clusterViewMode", "pendingSessionView", "runHistory"];
  const previousState = Object.fromEntries(resultFields.map((name) => [name, state[name]]));
  let updatingViews = false;
  state.isRunning = true;
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
    state.newickTree = prepared.tree;
    state.lastRunData = data;
    setCurrentClusters(prepared.clusterMap);
    computeLayouts();
    state.lastRunId = data.run_id || null;
    populateClusterSelector(data);
    updateClusterEditorAvailability();
    drawTree();
    state.lastRunSignature = JSON.stringify(payload);
    const leafCount = Object.keys(prepared.clusterMap).length;
    const clusterCount = new Set(Object.values(prepared.clusterMap)).size;
    const selectedCounts = data.k_values || data.ks || [];
    showStatus(leafCount ? `k = ${clusterCount}${selectedCounts.length > 1 ? ` (rank 1 of ${selectedCounts.length})` : ""} · ${leafCount} leaves · ${elapsedSeconds.toFixed(2)}s` : `No partitions selected · ${elapsedSeconds.toFixed(2)}s`, leafCount ? "success" : "info");
    const resultMode = typeof data.mode === "string" && data.mode ? data.mode : payload.mode;
    const modeLabel = payload.mode === "k" ? `k=${payload.k}` : resultMode !== payload.mode ? `${payload.mode} → ${resultMode}` : resultMode;
    state.runHistory = [{
      ts: Date.now(),
      mode: resultMode,
      label: `${modeLabel} · ${elapsedSeconds.toFixed(1)}s`,
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
    reportRunNotes(runNotes(data, payload.mode));
    refreshResultsStale();
  } catch (error) {
    if (updatingViews) {
      Object.assign(state, previousState, { isRunning: true });
      try {
        if (previousState.lastRunData) {
          populateClusterSelector(previousState.lastRunData);
          updateClusterEditorAvailability();
        }
        if (previousState.newickTree) drawTree();
        else clearTree();
      } catch (restoreError) {
        console.warn("Could not redraw the previous tree:", restoreError);
      }
    }
    console.error(error);
    showStatus(`Error: ${error.message}${previousState.lastRunData ? " (previous result retained)" : ""}`, "danger");
    refreshResultsStale();
  } finally {
    state.isRunning = false;
    if (runButton) {
      runButton.disabled = false;
      runButton.innerHTML = originalButtonContent;
    }
  }
}
