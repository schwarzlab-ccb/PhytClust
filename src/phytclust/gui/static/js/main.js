/* ============================================================
   PhytClust – main.js
   Entry point: DOM ready — wire everything.
   ============================================================ */

import { state, EXAMPLE_NEWICK, SELECTED_CLUSTER_IDS, BOX_LABEL_MAP, NODE_CUSTOM, getBoxAdjust } from "./state.js";
import { newickEl, extraResolutionEl } from "./dom.js";
import { resetExtraParams, collectLeafNamesData } from "./utils.js";
import { BASE_COLORS, shuffle, generateClusterColors, initTheme, toggleTheme } from "./colors.js";
import { showToast } from "./ui/toast.js";
import { showStatus, getStatusLog, clearStatusLog } from "./ui/status.js";
import { switchTab } from "./ui/tabs.js";
import { wireCollapse } from "./ui/collapsibles.js";
import { handleFileSelect, loadFile } from "./ui/file.js";
import { exportSvgFromEl, exportPngFromEl, exportTSV, saveToServer, copySvgToClipboard, copyRasterToClipboard } from "./ui/export.js";
import { wireAppearanceTrigger } from "./ui/appearance_panel.js";
import { collectSession, applySession, sessionSizeWarning } from "./session.js";
import { runPhytClust } from "./api.js";
import { drawTree, fitTree, fitCompare, clearAllCollapsedFlags, hideNodeContextMenu, showNodeContextMenu, getNodeCustom, nodeDisplayName, clearTree } from "./tree/draw.js";
import { drawComparison, makeCompareConfig, renderCompareConfigList, buildComparisonTSV, getClustersForK, normalizeClustersForLeaves } from "./views/compare.js";
import { drawOptimalK } from "./views/optimal_k.js";
import {
  getCurrentClusterIds,
  updateClusterEditorAvailability,
  openClusterEditor,
  closeClusterEditor,
  saveClusterEditor,
  renderClusterEditorRows,
  populateClusterSelector,
  switchCluster,
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

/* ─────────────────────────────────────
   DOM Ready – Wire Everything
   ───────────────────────────────────── */
document.addEventListener("DOMContentLoaded", function () {
  var appLoading = document.getElementById("app-loading");
  if (appLoading) appLoading.style.display = "none";
  showStatus("Ready", "success");

  // ── Appearance panel ──
  // Owns every visual option (labels, branches, nodes, cluster boxes, scale).
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
      state.NEWICK_RAW_TREE = null;
      state.latestApiData = null;
      state.CURRENT_CLUSTERS = {};
      state.CLUSTER_COLORS = [];
      if (newickEl) newickEl.value = "";
      if (fileInput) fileInput.value = "";
      var fileNameEl = document.getElementById("file-name");
      if (fileNameEl) fileNameEl.style.display = "none";
      clearTree();
      drawComparison();
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
  // Reuses the existing modal-custom skeleton; opened by `?`, closed by Esc /
  // backdrop click / X button. We deliberately ignore the keystroke when the
  // user is typing in an input/textarea so the search box and Newick paste
  // area still accept literal "?" or "/" characters.
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

  // Global single-key shortcuts. Skipped entirely when focus is in any kind
  // of text input so we don't hijack the user's typing.
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
  var resetBtn = document.getElementById("btn-extra-reset");
  if (resetBtn)
    resetBtn.addEventListener("click", function (e) {
      e.preventDefault();
      resetExtraParams();
      showToast("Parameters reset to defaults.", "info", 2000);
    });

  // ── Mode selector ──
  document.querySelectorAll("#mode-selector .mode-btn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      document
        .querySelectorAll("#mode-selector .mode-btn")
        .forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      var mode = btn.dataset.mode;
      var paramK = document.getElementById("param-k");
      var paramTopN = document.getElementById("param-topn");
      var paramBins = document.getElementById("param-bins");
      if (paramK) paramK.style.display = mode === "k" ? "" : "none";
      if (paramTopN) paramTopN.style.display = mode === "k" ? "none" : "";
      if (paramBins)
        paramBins.style.display = mode === "resolution" ? "" : "none";
      if (extraResolutionEl) extraResolutionEl.checked = mode === "resolution";
    });
  });

  // ── Collapsible sections ──
  wireCollapse("refine-opts-toggle", "refine-opts-body");
  wireCollapse("expert-opts-toggle", "expert-opts-body");

  // ── Tabs ──
  document.querySelectorAll(".tab-btn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      switchTab(btn.dataset.tab);
      if (btn.dataset.tab === "debug") refreshSessionPanel();
    });
  });
  switchTab("viewer");

  // ── Session JSON tab ──
  // The "JSON" tab shows a portable session document — sidebar params,
  // appearance state, optionally the Newick. It does NOT show the API
  // result dump (which used to live here); developers can use the browser
  // console for raw payloads.
  function getSessionOpts() {
    const cb = document.getElementById("session-include-newick");
    return { includeNewick: cb ? cb.checked : true };
  }
  function refreshSessionPanel() {
    const resultPane = document.getElementById("result");
    if (!resultPane) return;
    const session = collectSession(getSessionOpts());
    resultPane.textContent = JSON.stringify(session, null, 2);
    const warn = sessionSizeWarning(session);
    if (warn) console.info("[phytclust session] " + warn);
  }
  // Refresh whenever the user changes a sidebar input or the include-tree toggle.
  document.querySelectorAll(
    ".sidebar input, .sidebar select, .sidebar textarea, #session-include-newick",
  ).forEach(function (el) {
    el.addEventListener("change", function () {
      if (document.getElementById("debug").classList.contains("active")) {
        refreshSessionPanel();
      }
    });
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
          const { applied, errors, hadNewick } = applySession(parsed);
          if (errors.length) errors.forEach((m) => showToast(m, "info", 4000));
          if (hadNewick) setNewickAreaOpen(true);
          showToast(
            `Restored ${applied} setting${applied === 1 ? "" : "s"}` +
              (hadNewick ? " + tree" : ""),
            "success",
            2500,
          );
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
  // Opened by the clock button in the status bar. Shows a reverse-chronological
  // list of all showStatus messages (max 50). Auto-refreshes while open.
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
            '<span class="log-msg">' + r.label +
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
            '<span class="log-msg">' + e.message + "</span></div>"
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
          exportSvgFromEl("#tree_display", "phytclust_tree.svg");
        if (type === "tree_png") {
          var dpiStr = prompt(
            "Enter DPI for publication-ready PNG (default 300):",
            "300",
          );
          if (dpiStr === null) return;
          var dpi = parseInt(dpiStr, 10);
          if (isNaN(dpi) || dpi < 72) dpi = 300;
          exportPngFromEl("#tree_display", "phytclust_tree.png", dpi);
        }
        if (type === "k_plot")
          exportSvgFromEl("#optimalk_plot", "phytclust_scores.svg");
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
    targetSelector: "#tree_display",
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
  var btnCopyMaxk = document.getElementById("btn-copy-maxk");
  if (btnCopyMaxk)
    btnCopyMaxk.addEventListener("click", function () {
      copySvgToClipboard("#optimalk_plot");
    });

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

    if (isArrow && state.render.clusters.colorMode === "boxes" && SELECTED_CLUSTER_IDS.size > 0) {
      e.preventDefault();
      var step = 4;
      var resize = !!e.shiftKey;
      SELECTED_CLUSTER_IDS.forEach(function (cid) {
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

    if (e.key === "Escape" && SELECTED_CLUSTER_IDS.size > 0) {
      SELECTED_CLUSTER_IDS.clear();
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
      state.COMPARE_CONFIGS.push(makeCompareConfig());
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
    compareHighlightCb.checked = !!state.COMPARE_HIGHLIGHT_CHANGES;
    compareHighlightCb.addEventListener("change", function () {
      state.COMPARE_HIGHLIGHT_CHANGES = compareHighlightCb.checked;
      drawComparison();
    });
  }

  // ── Compare: reference bar selector ──
  var compareReferenceSel = document.getElementById("compare-reference");
  if (compareReferenceSel) {
    compareReferenceSel.addEventListener("change", function () {
      state.COMPARE_REFERENCE_INDEX = parseInt(compareReferenceSel.value, 10) || 0;
      drawComparison();
    });
  }

  function buildCurrentComparisons() {
    var comparisons = [];
    (state.COMPARE_CONFIGS || []).forEach(function (cfg) {
      var clusters = null;
      if (cfg.source === "k") {
        clusters = normalizeClustersForLeaves(getClustersForK(Number(cfg.k)));
      } else if (cfg.source === "file") {
        clusters = cfg.fileClusters || null;
      }
      if (clusters) comparisons.push({ title: cfg.title || "Bar " + cfg.id, clusters: clusters });
    });
    return comparisons;
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

  // ── Color mode ──
  var colorModeSelect = document.getElementById("color-mode");
  if (colorModeSelect)
    colorModeSelect.addEventListener("change", function () {
      state.render.clusters.colorMode = this.value || "bars";
      drawTree();
    });

  // Branch width / colour live in the Appearance panel.

  // ── Palette selector ──
  let CUSTOM_PALETTE = null;
  function normalizeHexColor(s) {
    const v = (s || "").trim();
    if (/^#[0-9a-fA-F]{6}$/.test(v)) return v;
    if (/^[0-9a-fA-F]{6}$/.test(v)) return "#" + v;
    return null;
  }
  function updatePalettePreview(colors) {
    const box = document.getElementById("palette-preview");
    if (box)
      box.innerHTML = colors
        .map((c) => `<span title="${c}" style="background:${c};"></span>`)
        .join("");
  }

  function applyPalette(type) {
    const palettes = {
      default: [
        "#b84b4b",
        "#849060",
        "#3d7c74",
        "#6e3f8a",
        "#ceb94b",
        "#3f648a",
        "#3f408a",
        "#da63aa",
      ],
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
    if (type === "custom" && CUSTOM_PALETTE && CUSTOM_PALETTE.length)
      BASE_COLORS.splice(0, BASE_COLORS.length, ...CUSTOM_PALETTE);
    else if (palettes[type])
      BASE_COLORS.splice(0, BASE_COLORS.length, ...palettes[type]);
    updatePalettePreview(BASE_COLORS);
    if (Object.keys(state.CURRENT_CLUSTERS || {}).length > 0) {
      state.CLUSTER_COLORS = generateClusterColors(
        Math.max(...Object.values(state.CURRENT_CLUSTERS)) + 1,
      );
      drawTree();
    }
  }

  var paletteSelect = document.getElementById("palette-select");
  if (paletteSelect) {
    paletteSelect.addEventListener("change", function () {
      if (this.value === "custom") {
        var initial =
          CUSTOM_PALETTE && CUSTOM_PALETTE.length
            ? CUSTOM_PALETTE.join(",")
            : BASE_COLORS.join(",");
        var raw = prompt("Enter hex colors separated by commas:", initial);
        if (raw == null) {
          this.value = "default";
          applyPalette("default");
          return;
        }
        var parsed = raw.split(",").map(normalizeHexColor).filter(Boolean);
        if (parsed.length < 2) {
          alert("Please provide at least 2 valid hex colors.");
          this.value = "default";
          applyPalette("default");
          return;
        }
        CUSTOM_PALETTE = parsed;
        applyPalette("custom");
        return;
      }
      applyPalette(this.value);
    });
  }
  updatePalettePreview(BASE_COLORS);

  // ── Shuffle colors button ──
  var btnShuffle = document.getElementById("btn-shuffle-colors");
  if (btnShuffle) {
    btnShuffle.addEventListener("click", function () {
      var shuffled = shuffle(BASE_COLORS.slice());
      BASE_COLORS.splice(0, BASE_COLORS.length);
      for (var i = 0; i < shuffled.length; i++) BASE_COLORS.push(shuffled[i]);
      updatePalettePreview(BASE_COLORS);
      if (Object.keys(state.CURRENT_CLUSTERS || {}).length > 0) {
        state.CLUSTER_COLORS = generateClusterColors(
          Math.max.apply(null, Object.values(state.CURRENT_CLUSTERS).map(Number)) + 1,
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
      if (state.LAST_TREE_SVG && state.LAST_TREE_ZOOM)
        state.LAST_TREE_SVG.transition()
          .duration(150)
          .call(state.LAST_TREE_ZOOM.transform, d3.zoomIdentity);
      state.render.scale.width = 1.0;
      state.render.scale.height = 1.0;
      clearAllCollapsedFlags(state.NEWICK_RAW_TREE);
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
    var term = (state.SEARCH_TERM || "").toLowerCase();
    if (!term) {
      searchCounter.hidden = true;
      searchCounter.classList.remove("no-matches");
      return;
    }
    var n = 0;
    if (state.NEWICK_RAW_TREE) {
      var leaves = [];
      collectLeafNamesData(state.NEWICK_RAW_TREE, leaves);
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
        state.SEARCH_TERM = searchInput.value.trim();
        if (state.NEWICK_RAW_TREE) drawTree();
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
      if (!state.CTX_TARGET_DATA) return;
      var c = getNodeCustom(state.CTX_TARGET_DATA);
      var val = prompt("Rename node:", nodeDisplayName(state.CTX_TARGET_DATA));
      if (val != null) {
        c.renamedTo = val;
        drawTree();
      }
    });
    btn("ctx-highlight", function () {
      if (!state.CTX_TARGET_DATA) return;
      var c = getNodeCustom(state.CTX_TARGET_DATA);
      c.highlighted = !c.highlighted;
      drawTree();
    });
    var ctxSlider = document.getElementById("ctx-size-slider");
    var ctxSliderVal = document.getElementById("ctx-size-value");
    if (ctxSlider) {
      ctxSlider.addEventListener("input", function () {
        if (!state.CTX_TARGET_DATA) return;
        var c = getNodeCustom(state.CTX_TARGET_DATA);
        c.radiusScale = parseFloat(ctxSlider.value);
        if (ctxSliderVal) ctxSliderVal.textContent = c.radiusScale.toFixed(1);
        drawTree();
      });
      ctxSlider.addEventListener("click", function (e) {
        e.stopPropagation();
      });
    }
    btn("ctx-copy-subtree", function () {
      if (!state.CTX_TARGET_DATA) return;
      var names = [];
      collectLeafNamesData(state.CTX_TARGET_DATA, names);
      navigator.clipboard.writeText(names.join("\n")).then(
        () => showToast(names.length + " leaf names copied", "success", 2000),
        () => showToast("Copy failed", "danger"),
      );
    });
    btn("ctx-original-size", function () {
      if (!state.CTX_TARGET_DATA) return;
      var c = getNodeCustom(state.CTX_TARGET_DATA);
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
    btn("ctx-select-cluster", function () {
      if (!state.CTX_TARGET_DATA || !state.CURRENT_CLUSTERS) return;
      var name = state.CTX_TARGET_DATA.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(state.CTX_TARGET_DATA, leaves);
        name = leaves[0];
      }
      var cid = name ? state.CURRENT_CLUSTERS[name] : null;
      if (cid == null) {
        showToast("No cluster for this node", "info");
        return;
      }
      var nCid = Number(cid);
      if (SELECTED_CLUSTER_IDS.has(nCid)) {
        SELECTED_CLUSTER_IDS.delete(nCid);
        if (SELECTED_CLUSTER_IDS.size === 0) {
          showToast("Cluster focus cleared", "info", 1800);
        } else {
          showToast(
            "Focused clusters: " +
              Array.from(SELECTED_CLUSTER_IDS)
                .sort((a, b) => a - b)
                .join(", "),
            "info",
            1800,
          );
        }
      } else {
        SELECTED_CLUSTER_IDS.add(nCid);
        showToast(
          "Focused clusters: " +
            Array.from(SELECTED_CLUSTER_IDS)
              .sort((a, b) => a - b)
              .join(", "),
          "success",
          1800,
        );
      }
      drawTree();
    });
    btn("ctx-rename-cluster", function () {
      if (!state.CTX_TARGET_DATA || !state.CURRENT_CLUSTERS) return;
      var name = state.CTX_TARGET_DATA.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(state.CTX_TARGET_DATA, leaves);
        name = leaves[0];
      }
      var cid = name ? state.CURRENT_CLUSTERS[name] : null;
      if (cid == null) {
        showToast("No cluster for this node", "info");
        return;
      }
      var current = BOX_LABEL_MAP[cid] || "C" + cid;
      var val = prompt("Label for cluster " + cid + ":", current);
      if (val != null) {
        BOX_LABEL_MAP[cid] = val;
        drawTree();
      }
    });
    btn("ctx-adjust-box", function () {
      if (!state.CTX_TARGET_DATA || !state.CURRENT_CLUSTERS) return;
      var name = state.CTX_TARGET_DATA.name;
      if (!name) {
        var leaves = [];
        collectLeafNamesData(state.CTX_TARGET_DATA, leaves);
        name = leaves[0];
      }
      var cid = name ? state.CURRENT_CLUSTERS[name] : null;
      if (cid == null) {
        showToast("No cluster for this node", "info");
        return;
      }
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
      if (!state.CTX_TARGET_DATA) return;
      NODE_CUSTOM.delete(state.CTX_TARGET_DATA);
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
