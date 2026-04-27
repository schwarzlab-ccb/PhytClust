/* ============================================================
   PhytClust – api.js
   API call + runPhytClust orchestration.
   ============================================================ */

import { state } from "./state.js";
import { newickEl, extraOutgroupEl, extraRootTaxonEl, extraResolutionEl } from "./dom.js";
import { showToast } from "./ui/toast.js";
import { showStatus } from "./ui/status.js";
import { estimateLeafCount, isOutgroupInNewick, parseNewick, readIntParam, readFloatParam, readCheckParam, readSelectParam } from "./utils.js";
import { generateClusterColors } from "./colors.js";
import { accumulateBranchLength, computeLayouts, drawTree, clearTree } from "./tree/draw.js";
import { populateClusterSelector, updateClusterEditorAvailability } from "./views/cluster_editor.js";
import { populateCompareSelectors } from "./views/compare.js";
import { drawOptimalK } from "./views/optimal_k.js";
import { drawMiniScores } from "./views/scores_panel.js";

async function apiPostJson(url, payload) {
  const res = await fetch(url, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  const text = await res.text();
  let data = null;
  try {
    data = text ? JSON.parse(text) : null;
  } catch {
    /* ignore */
  }
  if (!res.ok) {
    throw new Error(
      data && data.detail ? data.detail : `Request failed (${res.status})`,
    );
  }
  return data;
}

function getCurrentMode() {
  const activeBtn = document.querySelector("#mode-selector .mode-btn.active");
  return activeBtn ? activeBtn.dataset.mode : "global";
}

export async function runPhytClust() {
  if (state.isRunning) return;
  var newickText = (newickEl.value || "").trim();
  if (!newickText) {
    showToast("Please upload or paste a Newick tree.", "danger", 4000);
    return;
  }

  var numSamples = estimateLeafCount(newickText);
  var mode = getCurrentMode();

  var kVal = readIntParam("extra-k");
  var outgroupVal =
    (extraOutgroupEl ? (extraOutgroupEl.value || "").trim() : "") || null;
  var rootTaxonVal =
    (extraRootTaxonEl ? (extraRootTaxonEl.value || "").trim() : "") || null;
  var topNVal = readIntParam("extra-topn");
  var binsVal = readIntParam("extra-bins");
  var maxKVal = readIntParam("extra-maxk");
  var maxKLimitVal = readFloatParam("extra-maxklimit");
  var lambdaVal = readFloatParam("extra-lambda");
  var minClusterVal = readIntParam("extra-min-cluster-size");

  if (outgroupVal && !isOutgroupInNewick(newickText, outgroupVal)) {
    showStatus("Outgroup not found in Newick.", "danger");
    return;
  }

  if (extraResolutionEl) extraResolutionEl.checked = mode === "resolution";

  var payload = { newick: newickText, mode: mode };
  if (mode === "k") {
    if (kVal === null) {
      showToast("Please enter a value for k.", "danger");
      return;
    }
    payload.k = kVal;
  }
  if (outgroupVal) payload.outgroup = outgroupVal;
  if (rootTaxonVal) payload.root_taxon = rootTaxonVal;
  if (topNVal !== null) payload.top_n = topNVal;
  if (binsVal !== null) payload.num_bins = binsVal;
  if (maxKVal !== null) payload.max_k = maxKVal;
  if (maxKLimitVal !== null) payload.max_k_limit = maxKLimitVal;
  if (lambdaVal !== null) payload.lambda_weight = lambdaVal;
  if (minClusterVal !== null) payload.min_cluster_size = minClusterVal;
  if (mode === "resolution") payload.by_resolution = true;

  if (readCheckParam("extra-compute-all")) payload.compute_all_clusters = true;
  if (readCheckParam("extra-use-support")) payload.use_branch_support = true;
  var minSup = readFloatParam("extra-min-support");
  if (minSup !== null) payload.min_support = minSup;
  var supW = readFloatParam("extra-support-weight");
  if (supW !== null) payload.support_weight = supW;

  var outlierThresh = readIntParam("extra-outlier-threshold");
  if (outlierThresh !== null) payload.outlier_size_threshold = outlierThresh;
  if (readCheckParam("extra-outlier-prefer-fewer"))
    payload.outlier_prefer_fewer = true;
  var ratioMode = readSelectParam("extra-outlier-ratio-mode");
  if (ratioMode && ratioMode !== "exp") payload.outlier_ratio_mode = ratioMode;

  payload.optimize_polytomies = readCheckParam("extra-optimize-polytomies");
  if (readCheckParam("extra-no-split-zero"))
    payload.no_split_zero_length = true;

  var rankMode = readSelectParam("extra-ranking-mode");
  if (rankMode) payload.ranking_mode = rankMode;
  var minProm = readFloatParam("extra-min-prominence");
  if (minProm !== null) payload.min_prominence = minProm;
  if (readCheckParam("extra-relative-prom"))
    payload.use_relative_prominence = true;

  showStatus("Running PhytClust...", "info");
  clearTree();

  var runBtn = document.getElementById("btn-run");
  try {
    state.isRunning = true;
    if (runBtn) {
      runBtn.disabled = true;
      runBtn.innerHTML = '<span class="spinner"></span> Running...';
    }

    const t0 = performance.now();
    const data = await apiPostJson("/api/run", payload);
    const dt = (performance.now() - t0) / 1000;

    state.latestApiData = data;
    state.latestRunId = data.run_id || null;

    state.runHistory.unshift({
      ts: Date.now(),
      mode,
      label:
        mode === "k"
          ? "k=" + (data.k != null ? data.k : kVal) + " · " + dt.toFixed(1) + "s"
          : mode + " · " + dt.toFixed(1) + "s",
      nLeaves: numSamples,
    });
    if (state.runHistory.length > 8) state.runHistory.pop();

    showStatus(`Finished in ${dt.toFixed(2)}s`, "success");

    var lcLabel = document.getElementById("leaf-count-label");
    if (lcLabel && data.newick)
      lcLabel.textContent = estimateLeafCount(data.newick) + " leaves";

    if (data.newick) {
      state.NEWICK_RAW_TREE = parseNewick(data.newick);
      accumulateBranchLength(state.NEWICK_RAW_TREE);
      computeLayouts();
      populateClusterSelector(data);
      try {
        populateCompareSelectors(data);
      } catch (cmpErr) {
        console.warn("Compare panel init error:", cmpErr);
        showToast("Compare panel failed to initialize.", "danger", 2500);
      }

      const clusterMap =
        data.clusters && data.clusters.length > 0 ? data.clusters[0] : {};
      state.CURRENT_CLUSTERS = clusterMap;
      if (Object.keys(clusterMap).length > 0) {
        state.CLUSTER_COLORS = generateClusterColors(
          Math.max(...Object.values(clusterMap)) + 1,
        );
      } else {
        state.CLUSTER_COLORS = [];
        showStatus("No clusters found.", "danger");
      }
      updateClusterEditorAvailability();
      drawTree();
    } else {
      clearTree();
      showStatus("No Newick tree returned by API.", "danger");
    }

    try {
      drawMiniScores(data);
    } catch (e) {
      console.warn("Mini scores error:", e);
    }

    try {
      if (data && data.scores) {
        var plotHost = document.getElementById("optimalk_plot");
        if (plotHost) plotHost.innerHTML = "";
        state.latestOptimalKData = data;
        drawOptimalK(state.latestOptimalKData);
      } else {
        var plotHost2 = document.getElementById("optimalk_plot");
        if (plotHost2)
          plotHost2.innerHTML =
            '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:var(--pc-text-muted);font-size:14px;">Score plot unavailable in fixed k mode.</div>';
      }
    } catch (plotErr) {
      console.warn("Plot rendering error:", plotErr);
    }
  } catch (e) {
    console.error(e);
    showStatus("Error: " + e.message, "danger");
    state.latestOptimalKData = null;
    state.latestApiData = null;
  } finally {
    state.isRunning = false;
    if (runBtn) {
      runBtn.disabled = false;
      runBtn.innerHTML =
        '<svg width="14" height="14" viewBox="0 0 16 16" fill="currentColor"><polygon points="3,1 13,8 3,15"/></svg> Run PhytClust';
    }
  }
}
