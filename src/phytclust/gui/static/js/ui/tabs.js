// Switch between the tree, score, and comparison views.

import { state } from "../state.js";
import { drawScores } from "../views/scores.js";
import { drawComparison } from "../views/compare.js";

export function switchTab(tabName) {
  document
    .querySelectorAll(".tab-btn")
    .forEach((btn) =>
      btn.classList.toggle("active", btn.dataset.tab === tabName),
    );
  document
    .querySelectorAll(".viz-panel")
    .forEach((panel) => panel.classList.toggle("active", panel.id === tabName));

  var treeToolbar = document.getElementById("tree-toolbar");
  var scoresToolbar = document.getElementById("scores-toolbar");
  var compareToolbar = document.getElementById("compare-toolbar");
  var compareConfigList = document.getElementById("compare-config-list");
  if (treeToolbar)
    treeToolbar.style.display = tabName === "viewer" ? "" : "none";
  if (scoresToolbar)
    scoresToolbar.style.display = tabName === "viewer-scores" ? "" : "none";
  if (compareToolbar)
    compareToolbar.style.display = tabName === "compare" ? "" : "none";
  if (compareConfigList)
    compareConfigList.classList.toggle("visible", tabName === "compare");

  if (tabName === "viewer-scores" && state.lastScoresData)
    drawScores(state.lastScoresData);
  if (tabName === "compare") drawComparison();
}
