// Entry point: wire everything once the DOM is ready.

import { state, EXAMPLE_NEWICK, selectedClusterIds, boxLabelMap, nodeCustom, getBoxAdjust } from "./state.js";
import { newickEl } from "./dom.js";
import { resetParams, collectLeafNamesData, escapeHtmlAttr } from "./utils.js";
import { baseColours, shuffle, generateClusterColours, initTheme, toggleTheme } from "./colors.js";
import { showToast } from "./ui/toast.js";
import { showStatus, getStatusLog, clearStatusLog } from "./ui/status.js";
import { switchTab } from "./ui/tabs.js";
import { wireCollapse } from "./ui/collapsibles.js";
import { handleFileSelect, loadFile, reportSessionResult } from "./ui/file.js";
import { exportSvgFromEl, exportPngFromEl, exportTSV, saveToServer, copySvgToClipboard, copyRasterToClipboard } from "./ui/export.js";
import { wireAppearanceTrigger } from "./ui/appearance_panel.js";
import { collectSession, applySession, sessionSizeWarning } from "./session.js";
import { runPhytClust, refreshResultsStale, clearResults } from "./api.js";
import { drawTree, fitTree, clearAllCollapsedFlags } from "./tree/draw.js";
import { hideNodeContextMenu, getNodeCustom, nodeDisplayName } from "./tree/nodes.js";
import { fitCompare } from "./tree/compare.js";
import { drawComparison, makeCompareConfig, renderCompareConfigList, buildComparisonTSV, buildCurrentComparisons } from "./views/compare.js";
import { drawScores } from "./views/scores.js";
import {
  getCurrentClusterIds,
  updateClusterEditorAvailability,
  openClusterEditor,
  closeClusterEditor,
  saveClusterEditor,
  renderClusterEditorRows,
  switchCluster,
  refreshOutlierMarks,
  cycleCluster,
  toggleClusterViewMode,
} from "./views/cluster_editor.js";

// Initialize theme
initTheme();

// Wire context menu dismiss
document.addEventListener("click", hideNodeContextMenu);
document.addEventListener("contextmenu", function (e) {
  if (!e.target.closest(".tree-node")) hideNodeContextMenu();
});

document.addEventListener("DOMContentLoaded", function () {
  showStatus("Ready", "success");

  // ── Appearance panel ──
  // Build appearance controls for labels, branches, nodes, and cluster boxes.
  // Spec lives in ui/appearance_panel.js; controls flow through state.render
  // via setRenderOption + a redraw on every change.
  wireAppearanceTrigger({
    buttonId: ["appearance-toggle", "appearance-toggle-compare"],
    panelId: "appearance-panel",
    onChange: () => {
      drawTree();
      if (document.getElementById("compare").classList.contains("active"))
        drawComparison();
    },
  });

  document.addEventListener("phytclust:appearance-restored", () => {
    drawTree();
    if (document.getElementById("compare").classList.contains("active")) drawComparison();
  });

  // ── Sidebar toggle ──
  var sidebarToggle = document.getElementById("sidebar-toggle");
  var sidebar = document.getElementById("app-sidebar");
  if (sidebarToggle && sidebar) {
    sidebarToggle.addEventListener("click", function () {
      sidebar.classList.toggle("collapsed");
    });
  }

  // ── Theme toggle ──
  var themeToggle = document.getElementById("theme-toggle");
  if (themeToggle) {
    themeToggle.addEventListener("click", function () {
      toggleTheme();
      // Tree colours that resolve from CSS variables (branch, label, etc.)
      // need a redraw to pick up the new theme.
      try { drawTree(); } catch { /* tree may not be rendered yet */ }
    });
  }

  // ── File input ──
  var fileInput = document.getElementById("file-input");
  if (fileInput) fileInput.addEventListener("change", handleFileSelect);

  // ── Clear loaded tree (returns to empty state) ──
  var btnClearTree = document.getElementById("btn-clear-tree");
  if (btnClearTree) {
    btnClearTree.addEventListener("click", function () {
      if (newickEl) newickEl.value = "";
      if (fileInput) fileInput.value = "";
      var fileNameEl = document.getElementById("file-name");
      if (fileNameEl) fileNameEl.style.display = "none";
      selectedClusterIds.clear();
      clearResults();
      showStatus("Tree cleared", "info");
    });
  }

  // ── Newick paste toggle ──
  var newickToggle = document.getElementById("newick-toggle");
  var newickAreaEl = document.getElementById("newick-area");
  function setNewickAreaOpen(open) {
    if (!newickAreaEl || !newickToggle) return;
    newickAreaEl.hidden = !open;
    newickToggle.setAttribute("aria-expanded", String(open));
  }
  if (newickToggle) {
    newickToggle.addEventListener("click", function () {
      var opening = newickAreaEl && newickAreaEl.hidden;
      setNewickAreaOpen(opening);
      if (opening && newickEl) newickEl.focus();
    });
    // Auto-expand if newick already has content (page restored, session applied)
    if (newickEl && newickEl.value.trim()) setNewickAreaOpen(true);
  }

  // ── Drag and drop ──
  var uploadZone = document.getElementById("upload-zone");
  if (uploadZone) {
    uploadZone.addEventListener("dragover", function (e) {
      e.preventDefault();
      e.stopPropagation();
      uploadZone.classList.add("drag-over");
    });
    uploadZone.addEventListener("dragleave", function (e) {
      e.preventDefault();
      e.stopPropagation();
      uploadZone.classList.remove("drag-over");
    });
    uploadZone.addEventListener("drop", function (e) {
      e.preventDefault();
      e.stopPropagation();
      uploadZone.classList.remove("drag-over");
      if (e.dataTransfer.files && e.dataTransfer.files.length)
        loadFile(e.dataTransfer.files[0]);
    });
  }

  // ── Run button + keyboard shortcut ──
  var btnRun = document.getElementById("btn-run");
  if (btnRun)
    btnRun.addEventListener("click", function (e) {
      e.preventDefault();
      runPhytClust();
    });
  document.addEventListener("keydown", function (e) {
    if ((e.ctrlKey || e.metaKey) && e.key === "Enter") {
      e.preventDefault();
      runPhytClust();
    }
  });

  // ── Keyboard shortcuts modal ──
  // ? opens, Esc / backdrop / X close. Ignored while typing in an input so the
  // search box and Newick area still accept ? and /.
  var shortcutsModal = document.getElementById("shortcutsModal");
  var shortcutsBackdrop = document.getElementById("shortcutsBackdrop");
  function openShortcuts() {
    if (shortcutsModal) shortcutsModal.classList.add("show");
    if (shortcutsBackdrop) shortcutsBackdrop.classList.add("show");
  }
  function closeShortcuts() {
    if (shortcutsModal) shortcutsModal.classList.remove("show");
    if (shortcutsBackdrop) shortcutsBackdrop.classList.remove("show");
  }
  var btnShortcutsClose = document.getElementById("btn-shortcuts-close");
  if (btnShortcutsClose) btnShortcutsClose.addEventListener("click", closeShortcuts);
  if (shortcutsBackdrop) shortcutsBackdrop.addEventListener("click", closeShortcuts);

  // Single-key shortcuts; skipped while a text input has focus.
  document.addEventListener("keydown", function (e) {
    var t = e.target;
    if (
      t && (t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.tagName === "SELECT" || t.isContentEditable)
    ) return;
    if (e.ctrlKey || e.metaKey || e.altKey) return;

    if (e.key === "?") {
      e.preventDefault();
      shortcutsModal && shortcutsModal.classList.contains("show")
        ? closeShortcuts()
        : openShortcuts();
      return;
    }
    if (e.key === "/") {
      e.preventDefault();
      var s = document.getElementById("node-search");
      if (s) { s.focus(); s.select(); }
      return;
    }
    if (e.key === "t" || e.key === "T") {
      e.preventDefault();
      toggleTheme();
      try { drawTree(); } catch { /* tree may not be rendered yet */ }
      return;
    }
    if (e.key === "f" || e.key === "F") {
      e.preventDefault();
      var compareActive = document.getElementById("compare").classList.contains("active");
      compareActive ? fitCompare() : fitTree();
      return;
    }
    // Compare-tab-only shortcuts
    if (document.getElementById("compare").classList.contains("active")) {
      if (e.key === "d" || e.key === "D") {
        e.preventDefault();
        var cb = document.getElementById("compare-highlight-changes");
        if (cb) { cb.checked = !cb.checked; cb.dispatchEvent(new Event("change")); }
        return;
      }
      if (e.key === "," || e.key === ".") {
        e.preventDefault();
        var sel = document.getElementById("compare-reference");
        if (sel && sel.options.length > 1) {
          var dir = e.key === "." ? 1 : -1;
          var n = sel.options.length;
          sel.selectedIndex = (sel.selectedIndex + dir + n) % n;
          sel.dispatchEvent(new Event("change"));
        }
        return;
      }
    }
    if (e.key === "Escape") {
      // Close the shortcuts modal first if it's open; existing Esc handlers
      // deal with cluster focus / other modals.
      if (shortcutsModal && shortcutsModal.classList.contains("show")) {
        closeShortcuts();
      }
    }
  });

  // ── Empty state actions ──
  // The empty state is hidden once the viewer panel gains `.has-tree`,
  // which happens whenever drawTree successfully renders (see api.js).
  var btnEmptyExample = document.getElementById("btn-empty-example");
  if (btnEmptyExample) {
    btnEmptyExample.addEventListener("click", function () {
      newickEl.value = EXAMPLE_NEWICK;
      setNewickAreaOpen(true);
      newickEl.dispatchEvent(new Event("change", { bubbles: true }));
      runPhytClust();
    });
  }
  var btnEmptyHelp = document.getElementById("btn-empty-help");
  if (btnEmptyHelp) {
    btnEmptyHelp.addEventListener("click", function () {
      var helpBtn = document.getElementById("btn-help");
      if (helpBtn) helpBtn.click();
    });
  }

  // ── Reset button ──
  var resetBtn = document.getElementById("btn-reset-params");
  if (resetBtn)
    resetBtn.addEventListener("click", function (e) {
      e.preventDefault();
      resetParams();
      refreshResultsStale();
      showToast("Parameters reset to defaults.", "info", 2000);
    });

  // ── Outlier marking ──
  // "Mark outliers" changes how the current result is drawn and exported,
  // not the clustering itself, so redraw without asking for a new run.
  var markOutliersEl = document.getElementById("mark-outliers");
  if (markOutliersEl) {
    markOutliersEl.addEventListener("change", function () {
      refreshOutlierMarks();
      drawTree();
      renderClusterEditorRows();
    });
  }

  // ── Mode selector ──
  document.querySelectorAll("#mode-selector .mode-btn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      document
        .querySelectorAll("#mode-selector .mode-btn")
        .forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      var mode = btn.dataset.mode;
      var paramK = document.getElementById("param-k");
      var paramTopN = document.getElementById("param-top-n");
      var paramBins = document.getElementById("param-bins");
      if (paramK) paramK.style.display = mode === "k" ? "" : "none";
      if (paramTopN) paramTopN.style.display = mode === "k" ? "none" : "";
      if (paramBins)
        paramBins.style.display = mode === "resolution" ? "" : "none";
      refreshResultsStale();
    });
  });

  // ── Collapsible sections ──
  wireCollapse("refine-opts-toggle", "refine-opts-body");
  wireCollapse("expert-opts-toggle", "expert-opts-body");

  // ── Tabs ──
  document.querySelectorAll(".tab-btn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      switchTab(btn.dataset.tab);
      if (btn.dataset.tab === "session") refreshSessionPanel();
    });
  });
  switchTab("viewer");

  // ── Session JSON tab ──
  // A portable session document: params, appearance and optionally the Newick.
  function getSessionOpts() {
    const cb = document.getElementById("session-include-newick");
    return { includeNewick: cb ? cb.checked : true };
  }
  function refreshSessionPanel() {
    const sessionEl = document.getElementById("session-json");
    if (!sessionEl) return;
    const session = collectSession(getSessionOpts());
    sessionEl.textContent = JSON.stringify(session, null, 2);
    const warn = sessionSizeWarning(session);
    if (warn) console.info("[phytclust session] " + warn);
  }
  document.querySelectorAll(
    ".sidebar input, .sidebar select, .sidebar textarea, #session-include-newick",
  ).forEach(function (el) {
    for (const eventName of ["input", "change"]) {
      el.addEventListener(eventName, function () {
        refreshResultsStale();
        if (
          eventName === "change" &&
          document.getElementById("session").classList.contains("active")
        ) {
          refreshSessionPanel();
        }
      });
    }
  });

  var btnSessionDl = document.getElementById("btn-session-download");
  if (btnSessionDl) {
    btnSessionDl.addEventListener("click", function () {
      const session = collectSession(getSessionOpts());
      const blob = new Blob([JSON.stringify(session, null, 2)], {
        type: "application/json",
      });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "phytclust_session.json";
      document.body.appendChild(a);
      a.click();
      setTimeout(function () {
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
      }, 100);
      showToast("Session downloaded", "success", 1500);
    });
  }

  var btnSessionCopy = document.getElementById("btn-session-copy");
  if (btnSessionCopy) {
    btnSessionCopy.addEventListener("click", function () {
      const text = JSON.stringify(collectSession(getSessionOpts()), null, 2);
      navigator.clipboard.writeText(text).then(
        function () { showToast("Session copied", "success", 1500); },
        function () { showToast("Copy failed", "danger"); },
      );
    });
  }

  var sessionLoadInput = document.getElementById("session-load-input");
  if (sessionLoadInput) {
    sessionLoadInput.addEventListener("change", function () {
      const file = sessionLoadInput.files && sessionLoadInput.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = function (e) {
        try {
          const parsed = JSON.parse(e.target.result || "");
          const { valid, applied, errors, hadNewick } = applySession(parsed);
          reportSessionResult({ valid, applied, errors, hadNewick });
          if (hadNewick) setNewickAreaOpen(true);
          refreshSessionPanel();
        } catch (err) {
          showToast("Could not parse session JSON: " + err.message, "danger");
        } finally {
          sessionLoadInput.value = "";
        }
      };
      reader.readAsText(file);
    });
  }

  // ── Help sidebar ──
  var btnHelp = document.getElementById("btn-help");
  var helpSidebar = document.getElementById("help-sidebar");
  if (btnHelp && helpSidebar)
    btnHelp.addEventListener("click", function (e) {
      e.stopPropagation();
      helpSidebar.classList.toggle("open");
    });
  var btnHelpClose = document.getElementById("btn-help-close");
  if (btnHelpClose && helpSidebar)
    btnHelpClose.addEventListener("click", function () {
      helpSidebar.classList.remove("open");
    });
  if (helpSidebar) {
    document.addEventListener("click", function (e) {
      if (
        helpSidebar.classList.contains("open") &&
        !helpSidebar.contains(e.target) &&
        e.target !== btnHelp
      ) {
        helpSidebar.classList.remove("open");
      }
    });
    helpSidebar.addEventListener("click", function (e) {
      e.stopPropagation();
    });
  }

  // ── Activity log panel ──
  // Reverse-chronological list of showStatus messages (max 50).
  var btnActivityLog = document.getElementById("btn-activity-log");
  var activityLogPanel = document.getElementById("activity-log-panel");
  var btnActivityLogClear = document.getElementById("btn-activity-log-clear");

  function renderActivityLog() {
    var body = document.getElementById("activity-log-body");
    if (!body) return;
    var entries = getStatusLog();
    var runs = state.runHistory || [];
    var html = "";

    if (runs.length) {
      html += '<div class="log-section-title">Recent runs</div>';
      html += runs
        .map(function (r) {
          var t = new Date(r.ts);
          var ts = t.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
          return (
            '<div class="log-entry log-entry-run">' +
            '<span class="log-ts">' + ts + "</span>" +
            '<span class="log-msg">' + escapeHtmlAttr(r.label) +
            (r.nLeaves ? " · " + r.nLeaves + " leaves" : "") +
            "</span></div>"
          );
        })
        .join("");
      html += '<div class="log-section-title">Activity</div>';
    }

    if (!entries.length && !runs.length) {
      body.innerHTML = '<div class="log-entry-empty">No activity yet.</div>';
      return;
    }
    if (entries.length) {
      html += entries
        .slice()
        .reverse()
        .map(function (e) {
          var t = new Date(e.ts);
          var ts = t.toLocaleTimeString([], {
            hour: "2-digit", minute: "2-digit", second: "2-digit",
          });
          return (
            '<div class="log-entry log-entry-' + e.type + '">' +
            '<span class="log-ts">' + ts + "</span>" +
            '<span class="log-msg">' + escapeHtmlAttr(e.message) + "</span></div>"
          );
        })
        .join("");
    }
    body.innerHTML = html;
  }

  if (btnActivityLog && activityLogPanel) {
    btnActivityLog.addEventListener("click", function (e) {
      e.stopPropagation();
      var isOpen = !activityLogPanel.hidden;
      activityLogPanel.hidden = isOpen;
      if (!isOpen) {
        renderActivityLog();
        btnActivityLog.classList.remove("has-new");
      }
    });
    document.addEventListener("click", function (e) {
      if (
        activityLogPanel &&
        !activityLogPanel.hidden &&
        !activityLogPanel.contains(e.target) &&
        e.target !== btnActivityLog &&
        !btnActivityLog.contains(e.target)
      ) {
        activityLogPanel.hidden = true;
      }
    });
  }
  if (btnActivityLogClear) {
    btnActivityLogClear.addEventListener("click", function () {
      clearStatusLog();
      renderActivityLog();
    });
  }
  // Keep the panel live: re-render whenever showStatus pushes a new entry.
  window.addEventListener("phytclust:log-update", function () {
    if (btnActivityLog) btnActivityLog.classList.add("has-new");
    if (activityLogPanel && !activityLogPanel.hidden) renderActivityLog();
  });

  // ── About modal ──
  var btnAbout = document.getElementById("btn-about");
  var aboutModal = document.getElementById("aboutModal");
  var aboutBg = document.getElementById("aboutBackdrop");
  function openAbout() {
    if (aboutModal) aboutModal.classList.add("show");
    if (aboutBg) aboutBg.classList.add("show");
  }
  function closeAbout() {
    if (aboutModal) aboutModal.classList.remove("show");
    if (aboutBg) aboutBg.classList.remove("show");
  }
  if (btnAbout) btnAbout.addEventListener("click", openAbout);
  var btnAboutClose1 = document.getElementById("btn-about-close");
  var btnAboutClose2 = document.getElementById("btn-about-close2");
  if (btnAboutClose1) btnAboutClose1.addEventListener("click", closeAbout);
  if (btnAboutClose2) btnAboutClose2.addEventListener("click", closeAbout);
  if (aboutBg) aboutBg.addEventListener("click", closeAbout);

  // ── Help load example ──
  var helpExampleBtn = document.getElementById("btn-help-load-example");
  if (helpExampleBtn) {
    helpExampleBtn.addEventListener("click", function (e) {
      e.preventDefault();
      newickEl.value = EXAMPLE_NEWICK;
      setNewickAreaOpen(true);
      if (helpSidebar) helpSidebar.classList.remove("open");
      runPhytClust();
    });
  }

  // ── Save dropdown ──
  var btnSave = document.getElementById("btn-save");
  var saveDropdown = document.getElementById("save-dropdown");
  if (btnSave && saveDropdown) {
    btnSave.addEventListener("click", function (e) {
      e.stopPropagation();
      saveDropdown.classList.toggle("show");
    });
    document.addEventListener("click", function () {
      saveDropdown.classList.remove("show");
    });
    saveDropdown.addEventListener("click", function (e) {
      e.stopPropagation();
    });
    saveDropdown.querySelectorAll(".dropdown-item").forEach(function (item) {
      item.addEventListener("click", function () {
        saveDropdown.classList.remove("show");
        var type = item.dataset.saveType;
        if (type === "tsv") exportTSV();
        if (type === "tree_plot")
          exportSvgFromEl("#tree-display", "phytclust_tree.svg");
        if (type === "tree_png") {
          var dpiStr = prompt(
            "Enter DPI for publication-ready PNG (default 300):",
            "300",
          );
          if (dpiStr === null) return;
          var dpi = parseInt(dpiStr, 10);
          if (isNaN(dpi) || dpi < 72) dpi = 300;
          exportPngFromEl("#tree-display", "phytclust_tree.png", dpi);
        }
        if (type === "k_plot")
          exportSvgFromEl("#scores-plot", "phytclust_scores.svg");
        if (type === "save_server") saveToServer();
      });
    });
  }

  // ── Copy dropdown (shared between Viewer and Compare) ──
  // Wires a button -> dropdown pair to the SVG/PNG/JPG export functions
  // targeting `targetSelector`. Reads DPI from `dpiInputId` and the
  // "Publication style" preset from `publicationCheckboxId`.
  function wireCopyDropdown({ btnId, dropdownId, targetSelector, dpiInputId, publicationCheckboxId, redraw }) {
    var btn = document.getElementById(btnId);
    var dd = document.getElementById(dropdownId);
    if (!btn || !dd) return;
    btn.addEventListener("click", function (e) {
      e.stopPropagation();
      dd.classList.toggle("show");
    });
    document.addEventListener("click", function () {
      dd.classList.remove("show");
    });
    dd.addEventListener("click", function (e) { e.stopPropagation(); });
    dd.querySelectorAll(".dropdown-item").forEach(function (item) {
      item.addEventListener("click", function () {
        dd.classList.remove("show");
        var fmt = item.dataset.copyFmt;
        var dpiEl = document.getElementById(dpiInputId);
        var dpi = dpiEl ? parseInt(dpiEl.value, 10) : 300;
        if (isNaN(dpi) || dpi < 72) dpi = 300;
        var pubEl = document.getElementById(publicationCheckboxId);
        var opts = { usePreset: !!(pubEl && pubEl.checked), redraw: redraw };
        if (fmt === "svg") {
          copySvgToClipboard(targetSelector, opts);
        } else if (fmt === "png") {
          copyRasterToClipboard(targetSelector, dpi, "image/png", opts);
        } else if (fmt === "jpg") {
          copyRasterToClipboard(targetSelector, dpi, "image/jpeg", opts);
        }
      });
    });
  }

  wireCopyDropdown({
    btnId: "btn-copy-tree",
    dropdownId: "copy-dropdown",
    targetSelector: "#tree-display",
    dpiInputId: "copy-dpi",
    publicationCheckboxId: "copy-publication",
    redraw: drawTree,
  });
  wireCopyDropdown({
    btnId: "btn-copy-compare",
    dropdownId: "copy-dropdown-compare",
    targetSelector: "#compare-tree-left",
    dpiInputId: "copy-dpi-compare",
    publicationCheckboxId: "copy-publication-compare",
    redraw: drawComparison,
  });
  var btnCopyScores = document.getElementById("btn-copy-scores");
  if (btnCopyScores)
    btnCopyScores.addEventListener("click", function () {
      copySvgToClipboard("#scores-plot");
    });

  var scoresWidthScaleEl = document.getElementById("scores-width-scale");
  if (scoresWidthScaleEl) {
    scoresWidthScaleEl.addEventListener("change", function () {
      drawScores(state.lastScoresData);
    });
  }

  // ── Cluster selector ──
  var clusterSelect = document.getElementById("cluster-select");
  if (clusterSelect)
    clusterSelect.addEventListener("change", function () {
      switchCluster(parseInt(this.value, 10));
    });
  var clusterToggle = document.getElementById("cluster-view-toggle");
  if (clusterToggle)
    clusterToggle.addEventListener("click", toggleClusterViewMode);

  // ── Cluster bulk editor ──
  var editClustersBtn = document.getElementById("btn-edit-clusters");
  var clusterEditorBg = document.getElementById("clusterEditorBackdrop");
  var clusterEditorClose = document.getElementById("btn-cluster-editor-close");
  var clusterEditorCancel = document.getElementById(
    "btn-cluster-editor-cancel",
  );
  var clusterEditorSave = document.getElementById("btn-cluster-editor-save");
  if (editClustersBtn)
    editClustersBtn.addEventListener("click", function () {
      if (!getCurrentClusterIds().length) {
        showToast("No clusters available yet.", "info", 1600);
        return;
      }
      openClusterEditor();
    });
  if (clusterEditorClose)
    clusterEditorClose.addEventListener("click", closeClusterEditor);
  if (clusterEditorCancel)
    clusterEditorCancel.addEventListener("click", closeClusterEditor);
  if (clusterEditorSave)
    clusterEditorSave.addEventListener("click", saveClusterEditor);
  if (clusterEditorBg)
    clusterEditorBg.addEventListener("click", closeClusterEditor);
  var clusterFilterInput = document.getElementById("cluster-filter-text");
  var clusterFilterOutliersCb = document.getElementById(
    "cluster-filter-show-outliers",
  );
  if (clusterFilterInput)
    clusterFilterInput.addEventListener("input", renderClusterEditorRows);
  if (clusterFilterOutliersCb)
    clusterFilterOutliersCb.addEventListener("change", renderClusterEditorRows);

  // ── Cluster prev/next buttons ──
  var clusterPrev = document.getElementById("cluster-prev");
  var clusterNext = document.getElementById("cluster-next");
  if (clusterPrev)
    clusterPrev.addEventListener("click", function () {
      cycleCluster(-1);
    });
  if (clusterNext)
    clusterNext.addEventListener("click", function () {
      cycleCluster(1);
    });

  // ── Arrow key cycling through clusters ──
  document.addEventListener("keydown", function (e) {
    if (
      e.target.tagName === "INPUT" ||
      e.target.tagName === "TEXTAREA" ||
      e.target.tagName === "SELECT"
    )
      return;

    var isArrow =
      e.key === "ArrowLeft" ||
      e.key === "ArrowRight" ||
      e.key === "ArrowUp" ||
      e.key === "ArrowDown";

    if (isArrow && state.render.clusters.colorMode === "boxes" && selectedClusterIds.size > 0) {
      e.preventDefault();
      var step = 4;
      var resize = !!e.shiftKey;
      selectedClusterIds.forEach(function (cid) {
        var adj = getBoxAdjust(cid);
        if (resize) {
          if (e.key === "ArrowLeft") adj.padX = Math.max(0, adj.padX - step);
          if (e.key === "ArrowRight") adj.padX += step;
          if (e.key === "ArrowUp") adj.padY += step;
          if (e.key === "ArrowDown") adj.padY = Math.max(0, adj.padY - step);
        } else {
          if (e.key === "ArrowLeft") adj.dx -= step;
          if (e.key === "ArrowRight") adj.dx += step;
          if (e.key === "ArrowUp") adj.dy -= step;
          if (e.key === "ArrowDown") adj.dy += step;
        }
      });
      drawTree();
      return;
    }

    if (e.key === "Escape" && selectedClusterIds.size > 0) {
      selectedClusterIds.clear();
      drawTree();
      showToast("Cluster focus cleared", "info", 1500);
      return;
    }

    if (e.key === "ArrowLeft") {
      e.preventDefault();
      cycleCluster(-1);
    }
    if (e.key === "ArrowRight") {
      e.preventDefault();
      cycleCluster(1);
    }
  });

  // ── Compare add bar ──
  var compareAddBtn = document.getElementById("compare-add-bar");
  if (compareAddBtn) {
    compareAddBtn.addEventListener("click", function () {
      state.compareConfigs.push(makeCompareConfig());
      renderCompareConfigList();
      drawComparison();
    });
  }

  // ── Fit compare ──
  var btnFitCompare = document.getElementById("btn-fit-compare");
  if (btnFitCompare) btnFitCompare.addEventListener("click", function () { fitCompare(); });

  // ── Compare: highlight changed leaves ──
  var compareHighlightCb = document.getElementById("compare-highlight-changes");
  if (compareHighlightCb) {
    compareHighlightCb.checked = !!state.compareHighlightChanges;
    compareHighlightCb.addEventListener("change", function () {
      state.compareHighlightChanges = compareHighlightCb.checked;
      drawComparison();
    });
  }

  // ── Compare: reference bar selector ──
  var compareReferenceSel = document.getElementById("compare-reference");
  if (compareReferenceSel) {
    compareReferenceSel.addEventListener("change", function () {
      state.compareReferenceIndex = parseInt(compareReferenceSel.value, 10) || 0;
      drawComparison();
    });
  }

  // ── Copy TSV ──
  var btnCopyTsv = document.getElementById("btn-copy-tsv-compare");
  if (btnCopyTsv) {
    btnCopyTsv.addEventListener("click", function () {
      var tsv = buildComparisonTSV(buildCurrentComparisons());
      if (!tsv) { showToast("No comparison data to export.", "info", 2000); return; }
      navigator.clipboard.writeText(tsv).then(
        function () { showToast("Copied TSV to clipboard", "success", 1800); },
        function () { showToast("Copy failed", "danger"); },
      );
    });
  }

  // ── Download TSV ──
  var btnDownloadTsv = document.getElementById("btn-download-tsv-compare");
  if (btnDownloadTsv) {
    btnDownloadTsv.addEventListener("click", function () {
      var tsv = buildComparisonTSV(buildCurrentComparisons());
      if (!tsv) { showToast("No comparison data to export.", "info", 2000); return; }
      var blob = new Blob([tsv], { type: "text/tab-separated-values" });
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a");
      a.href = url;
      a.download = "comparison.tsv";
      document.body.appendChild(a);
      a.click();
      setTimeout(function () { document.body.removeChild(a); URL.revokeObjectURL(url); }, 100);
      showToast("Downloaded comparison.tsv", "success", 1800);
    });
  }

  // ── Layout mode ──
  var layoutSelect = document.getElementById("layout-mode");
  if (layoutSelect)
    layoutSelect.addEventListener("change", function () {
      state.render.layout = this.value;
      drawTree();
      if (document.getElementById("compare").classList.contains("active"))
        drawComparison();
    });

  // ── Colour mode ──
  var colourModeSelect = document.getElementById("colour-mode");
  if (colourModeSelect)
    colourModeSelect.addEventListener("change", function () {
      state.render.clusters.colorMode = this.value || "bars";
      drawTree();
    });

  // Branch width / colour live in the Appearance panel.

  // ── Palette selector ──
  const DEFAULT_PALETTE = baseColours.slice();
  let customPalette = null;
  function normaliseHexColour(s) {
    const v = (s || "").trim();
    if (/^#[0-9a-fA-F]{6}$/.test(v)) return v;
    if (/^[0-9a-fA-F]{6}$/.test(v)) return "#" + v;
    return null;
  }
  function updatePalettePreview(colours) {
    const box = document.getElementById("palette-preview");
    if (box)
      box.innerHTML = colours
        .map((c) => `<span title="${c}" style="background:${c};"></span>`)
        .join("");
  }

  function applyPalette(type) {
    const palettes = {
      default: DEFAULT_PALETTE,
      pastel: [
        "#ffb3ba",
        "#ffdfba",
        "#ffffba",
        "#baffc9",
        "#bae1ff",
        "#d7baff",
        "#ffcce6",
        "#c2f0c2",
      ],
      vivid: [
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#984ea3",
        "#ff7f00",
        "#ffff33",
        "#a65628",
        "#f781bf",
      ],
      dark: [
        "#4e79a7",
        "#59a14f",
        "#e15759",
        "#b07aa1",
        "#edc948",
        "#76b7b2",
        "#ff9da7",
        "#9c755f",
      ],
    };
    if (type === "custom" && customPalette && customPalette.length)
      baseColours.splice(0, baseColours.length, ...customPalette);
    else if (palettes[type])
      baseColours.splice(0, baseColours.length, ...palettes[type]);
    updatePalettePreview(baseColours);
    if (Object.keys(state.runClusters || {}).length > 0) {
      state.clusterColours = generateClusterColours(
        Math.max(...Object.values(state.runClusters).map(Number)) + 1,
      );
      drawTree();
    }
  }

  var paletteSelect = document.getElementById("palette-select");
  if (paletteSelect) {
    paletteSelect.addEventListener("change", function () {
      if (this.value === "custom") {
        var initial =
          customPalette && customPalette.length
            ? customPalette.join(",")
            : baseColours.join(",");
        var raw = prompt("Enter hex colours separated by commas:", initial);
        if (raw == null) {
          this.value = "default";
          applyPalette("default");
          return;
        }
        var parsed = raw.split(",").map(normaliseHexColour).filter(Boolean);
        if (parsed.length < 2) {
          alert("Please provide at least 2 valid hex colours.");
          this.value = "default";
          applyPalette("default");
          return;
        }
        customPalette = parsed;
        applyPalette("custom");
        return;
      }
      applyPalette(this.value);
    });
  }
  updatePalettePreview(baseColours);

  // ── Shuffle colours button ──
  var btnShuffle = document.getElementById("btn-shuffle-colours");
  if (btnShuffle) {
    btnShuffle.addEventListener("click", function () {
      var shuffled = shuffle(baseColours.slice());
      baseColours.splice(0, baseColours.length);
      for (var i = 0; i < shuffled.length; i++) baseColours.push(shuffled[i]);
      updatePalettePreview(baseColours);
      if (Object.keys(state.runClusters || {}).length > 0) {
        state.clusterColours = generateClusterColours(
          Math.max.apply(null, Object.values(state.runClusters).map(Number)) + 1,
        );
        drawTree();
      }
    });
  }

  // Label toggles, box options, node/label sizes, and scale sliders all
  // live in the Appearance panel and write through state.render directly.

  // ── Fit to screen ──
  var btnFitTree = document.getElementById("btn-fit-tree");
  if (btnFitTree) btnFitTree.addEventListener("click", function () { fitTree(); });

  // ── Reset view ──
  var btnResetView = document.getElementById("btn-reset-view");
  if (btnResetView) {
    btnResetView.addEventListener("click", function (e) {
      e.preventDefault();
      if (state.lastTreeSvg && state.lastTreeZoom)
        state.lastTreeSvg.transition()
          .duration(150)
          .call(state.lastTreeZoom.transform, d3.zoomIdentity);
      state.render.scale.width = 1.0;
      state.render.scale.height = 1.0;
      clearAllCollapsedFlags(state.newickTree);
      drawTree();
    });
  }

  // ── Search ──
  // After every keystroke (debounced) we redraw to refresh the highlight,
  // count matches by walking the tree leaves, and update the counter badge.
  var searchInput = document.getElementById("node-search");
  var searchCounter = document.getElementById("search-counter");
  function updateSearchCounter() {
    if (!searchCounter) return;
    var term = (state.searchTerm || "").toLowerCase();
    if (!term) {
      searchCounter.hidden = true;
      searchCounter.classList.remove("no-matches");
      return;
    }
    var n = 0;
    if (state.newickTree) {
      var leaves = [];
      collectLeafNamesData(state.newickTree, leaves);
      for (var i = 0; i < leaves.length; i++) {
        if (String(leaves[i]).toLowerCase().includes(term)) n++;
      }
    }
    searchCounter.hidden = false;
    searchCounter.textContent = n + (n === 1 ? " match" : " matches");
    searchCounter.classList.toggle("no-matches", n === 0);
  }
  if (searchInput) {
    let searchTimer = null;
    searchInput.addEventListener("input", function () {
      clearTimeout(searchTimer);
      searchTimer = setTimeout(function () {
        state.searchTerm = searchInput.value.trim();
        if (state.newickTree) drawTree();
        updateSearchCounter();
      }, 200);
    });
  }

  // ── Context menu buttons ──
  (function initContextMenuButtons() {
    function btn(id, fn) {
      var el = document.getElementById(id);
      if (el)
        el.addEventListener("click", function () {
          fn();
          hideNodeContextMenu();
        });
    }
    btn("ctx-rename", function () {
      if (!state.contextTargetData) return;
      var c = getNodeCustom(state.contextTargetData);
      var val = prompt("Rename node:", nodeDisplayName(state.contextTargetData));
      if (val != null) {
        c.renamedTo = val;
        drawTree();
      }
    });
    btn("ctx-highlight", function () {
      if (!state.contextTargetData) return;
      var c = getNodeCustom(state.contextTargetData);
      c.highlighted = !c.highlighted;
      drawTree();
    });
    var ctxSlider = document.getElementById("ctx-size-slider");
    var ctxSliderVal = document.getElementById("ctx-size-value");
    if (ctxSlider) {
      ctxSlider.addEventListener("input", function () {
        if (!state.contextTargetData) return;
        var c = getNodeCustom(state.contextTargetData);
        c.radiusScale = parseFloat(ctxSlider.value);
        if (ctxSliderVal) ctxSliderVal.textContent = c.radiusScale.toFixed(1);
        drawTree();
      });
      ctxSlider.addEventListener("click", function (e) {
        e.stopPropagation();
      });
    }
    btn("ctx-copy-subtree", function () {
      if (!state.contextTargetData) return;
      var names = [];
      collectLeafNamesData(state.contextTargetData, names);
      navigator.clipboard.writeText(names.join("\n")).then(
        () => showToast(names.length + " leaf names copied", "success", 2000),
        () => showToast("Copy failed", "danger"),
      );
    });
    btn("ctx-original-size", function () {
      if (!state.contextTargetData) return;
      var c = getNodeCustom(state.contextTargetData);
      c.lockOriginalSize = !c.lockOriginalSize;
      if (c.lockOriginalSize) c.radiusScale = 1;
      showToast(
        c.lockOriginalSize
          ? "Node locked to original size."
          : "Node returned to auto collapsed sizing.",
        "info",
        1600,
      );
      drawTree();
    });
    // Cluster id of the context-menu target (its first leaf for internal
    // nodes), or null after a toast when there is none.
    function ctxTargetClusterId() {
      if (!state.contextTargetData || !state.currentClusters) return null;
      var name = state.contextTargetData.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(state.contextTargetData, leaves);
        name = leaves[0];
      }
      var cid = name ? state.currentClusters[name] : null;
      if (cid == null) showToast("No cluster for this node", "info");
      return cid;
    }
    btn("ctx-select-cluster", function () {
      var cid = ctxTargetClusterId();
      if (cid == null) return;
      var nCid = Number(cid);
      if (selectedClusterIds.has(nCid)) {
        selectedClusterIds.delete(nCid);
        if (selectedClusterIds.size === 0) {
          showToast("Cluster focus cleared", "info", 1800);
        } else {
          showToast(
            "Focused clusters: " +
              Array.from(selectedClusterIds)
                .sort((a, b) => a - b)
                .join(", "),
            "info",
            1800,
          );
        }
      } else {
        selectedClusterIds.add(nCid);
        showToast(
          "Focused clusters: " +
            Array.from(selectedClusterIds)
              .sort((a, b) => a - b)
              .join(", "),
          "success",
          1800,
        );
      }
      drawTree();
    });
    btn("ctx-rename-cluster", function () {
      var cid = ctxTargetClusterId();
      if (cid == null) return;
      var current = boxLabelMap[cid] || "C" + cid;
      var val = prompt("Label for cluster " + cid + ":", current);
      if (val != null) {
        boxLabelMap[cid] = val;
        drawTree();
      }
    });
    btn("ctx-adjust-box", function () {
      var cid = ctxTargetClusterId();
      if (cid == null) return;
      var adj = getBoxAdjust(cid);
      var padStr = prompt(
        "Box padding for cluster " + cid + " (padX, padY in px):",
        adj.padX + ", " + adj.padY,
      );
      if (padStr != null) {
        var parts = padStr.split(",").map(function (s) {
          return parseFloat(s.trim());
        });
        if (parts.length >= 1 && !isNaN(parts[0])) adj.padX = parts[0];
        if (parts.length >= 2 && !isNaN(parts[1])) adj.padY = parts[1];
        drawTree();
      }
    });
    btn("ctx-reset", function () {
      if (!state.contextTargetData) return;
      nodeCustom.delete(state.contextTargetData);
      drawTree();
    });
  })();

  // ── Mini scores panel toggle + drag ──
  var miniPanel = document.getElementById("mini-scores-panel");
  var miniScoresToggle = document.getElementById("mini-scores-toggle");
  var miniScoresBody = document.getElementById("mini-scores-body");
  if (miniScoresToggle && miniScoresBody) {
    miniScoresToggle.addEventListener("click", function (e) {
      e.stopPropagation();
      miniScoresBody.classList.toggle("collapsed");
      miniScoresToggle.innerHTML = miniScoresBody.classList.contains(
        "collapsed",
      )
        ? "&#x25B2;"
        : "&#x25BC;";
    });
  }
  if (miniPanel) {
    var miniHeader = miniPanel.querySelector(".mini-scores-header");
    if (miniHeader) {
      var dragState = {
        dragging: false,
        startX: 0,
        startY: 0,
        origLeft: 0,
        origTop: 0,
      };
      miniHeader.addEventListener("mousedown", function (e) {
        if (e.target === miniScoresToggle) return;
        dragState.dragging = true;
        var rect = miniPanel.getBoundingClientRect();
        dragState.startX = e.clientX;
        dragState.startY = e.clientY;
        dragState.origLeft = rect.left;
        dragState.origTop = rect.top;
        e.preventDefault();
      });
      document.addEventListener("mousemove", function (e) {
        if (!dragState.dragging) return;
        var dx = e.clientX - dragState.startX;
        var dy = e.clientY - dragState.startY;
        miniPanel.style.position = "fixed";
        miniPanel.style.left = dragState.origLeft + dx + "px";
        miniPanel.style.top = dragState.origTop + dy + "px";
        miniPanel.style.right = "auto";
        miniPanel.style.bottom = "auto";
      });
      document.addEventListener("mouseup", function () {
        dragState.dragging = false;
      });
    }
  }

  updateClusterEditorAvailability();
});
