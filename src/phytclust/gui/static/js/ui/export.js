/* ============================================================
   PhytClust – ui/export.js
   Download, SVG/PNG export, clipboard helpers.
   ============================================================ */

import { state, withRenderOverrides } from "../state.js";
import { drawTree } from "../tree/draw.js";
import { showToast } from "./toast.js";

/**
 * Render-option overrides applied when "Publication style" is requested.
 * Bigger fonts, thicker strokes, fully-opaque cluster boxes — tuned for
 * inclusion in print figures rather than on-screen browsing.
 */
export const PUBLICATION_PRESET = {
  branches: { width: 2.0 },
  labels: { fontSize: 12, axisFontSize: 12 },
  nodes: { leafRadius: 4.0 },
  clusters: { boxAlpha: 0.25, labelFontSize: 11 },
};

/**
 * Run `fn` with the live tree temporarily redrawn under PUBLICATION_PRESET
 * (or no overrides if `usePreset` is false). The view flashes briefly while
 * the snapshot is taken, then the original render state is restored.
 *
 * `redraw` defaults to `drawTree`; callers exporting from the Compare tab
 * pass `drawComparison` so the preset actually applies to that view before
 * the snapshot.
 */
function withExportPreset(usePreset, fn, redraw = drawTree) {
  if (!usePreset) return fn();
  return withRenderOverrides(PUBLICATION_PRESET, redraw, fn);
}

export function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => {
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }, 100);
}

export function getTightSvgSnapshot(svgNode, paddingPx) {
  if (!svgNode) return null;
  const padding = Math.max(0, Number(paddingPx || 10));

  let zoomLayer = null;
  const firstChild = svgNode.firstElementChild;
  if (firstChild && firstChild.tagName.toLowerCase() === "g") {
    zoomLayer = firstChild;
  }

  let target = zoomLayer || svgNode.querySelector("g") || svgNode;

  let bbox = null;
  try {
    bbox = target.getBBox();
  } catch {
    bbox = null;
  }

  const fallbackW =
    parseFloat(svgNode.getAttribute("width")) || svgNode.clientWidth || 800;
  const fallbackH =
    parseFloat(svgNode.getAttribute("height")) || svgNode.clientHeight || 600;

  const x = bbox && isFinite(bbox.x) ? bbox.x - padding : 0;
  const y = bbox && isFinite(bbox.y) ? bbox.y - padding : 0;
  const w =
    bbox && isFinite(bbox.width) && bbox.width > 0
      ? bbox.width + padding * 2
      : fallbackW;
  const h =
    bbox && isFinite(bbox.height) && bbox.height > 0
      ? bbox.height + padding * 2
      : fallbackH;

  const clone = svgNode.cloneNode(true);
  if (!clone.getAttribute("xmlns")) {
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  }

  const cloneFirstChild = clone.firstElementChild;
  if (cloneFirstChild && cloneFirstChild.tagName.toLowerCase() === "g") {
    cloneFirstChild.removeAttribute("transform");
  }

  clone.setAttribute("viewBox", `${x} ${y} ${w} ${h}`);
  clone.setAttribute("width", String(w));
  clone.setAttribute("height", String(h));

  return {
    svgData: new XMLSerializer().serializeToString(clone),
    width: w,
    height: h,
  };
}

export function exportSvgFromEl(selector, filename, opts = {}) {
  withExportPreset(opts.usePreset, () => {
    const svgNode = document.querySelector(selector + " svg");
    if (!svgNode) {
      showToast("No visualization to export.", "danger");
      return;
    }
    const snap = getTightSvgSnapshot(svgNode, 10);
    if (!snap) {
      showToast("SVG export failed.", "danger");
      return;
    }
    downloadBlob(
      new Blob([snap.svgData], {
        type: "image/svg+xml",
      }),
      filename,
    );
    showToast("Downloaded " + filename, "success", 2000);
  }, opts.redraw);
}

export function exportPngFromEl(selector, filename, dpi, opts = {}) {
  let snap = null;
  let svgNode = null;
  withExportPreset(opts.usePreset, () => {
    svgNode = document.querySelector(selector + " svg");
    if (!svgNode) return;
    snap = getTightSvgSnapshot(svgNode, 10);
  }, opts.redraw);
  if (!svgNode) {
    showToast("No visualization to export.", "danger");
    return;
  }
  if (!snap) {
    showToast("PNG export failed.", "danger");
    return;
  }
  dpi = dpi || 300;
  const scale = dpi / 96;
  const svgData = snap.svgData;
  const svgW = snap.width;
  const svgH = snap.height;
  const canvas = document.createElement("canvas");
  canvas.width = Math.ceil(svgW * scale);
  canvas.height = Math.ceil(svgH * scale);
  const ctx = canvas.getContext("2d");
  const bgColor =
    getComputedStyle(document.documentElement)
      .getPropertyValue("--pc-bg")
      .trim() || "#ffffff";
  ctx.fillStyle = bgColor;
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.scale(scale, scale);
  const img = new Image();
  const blob = new Blob([svgData], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  img.onload = function () {
    ctx.drawImage(img, 0, 0, svgW, svgH);
    URL.revokeObjectURL(url);
    canvas.toBlob(function (pngBlob) {
      if (!pngBlob) {
        showToast("PNG export failed.", "danger");
        return;
      }
      downloadBlob(pngBlob, filename);
      showToast(
        "Downloaded " + filename + " at " + dpi + " DPI",
        "success",
        3000,
      );
    }, "image/png");
  };
  img.onerror = function () {
    URL.revokeObjectURL(url);
    showToast("PNG rendering failed.", "danger");
  };
  img.src = url;
}

export function copyRasterToClipboard(selector, dpi, mimeType, opts = {}) {
  var snap = null;
  var svgNode = null;
  withExportPreset(opts && opts.usePreset, function () {
    svgNode = document.querySelector(selector + " svg");
    if (!svgNode) return;
    snap = getTightSvgSnapshot(svgNode, 10);
  }, opts && opts.redraw);
  if (!svgNode) {
    showToast("Nothing to copy.", "danger");
    return;
  }
  if (!snap) {
    showToast("Copy failed.", "danger");
    return;
  }
  dpi = dpi || 300;
  mimeType = mimeType || "image/png";
  var fmtLabel = mimeType === "image/jpeg" ? "JPG" : "PNG";
  var scale = dpi / 96;
  var svgData = snap.svgData;
  var svgW = snap.width;
  var svgH = snap.height;
  var canvas = document.createElement("canvas");
  canvas.width = Math.ceil(svgW * scale);
  canvas.height = Math.ceil(svgH * scale);
  var ctx = canvas.getContext("2d");
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.scale(scale, scale);
  var img = new Image();
  var blob = new Blob([svgData], { type: "image/svg+xml;charset=utf-8" });
  var url = URL.createObjectURL(blob);
  img.onload = function () {
    ctx.drawImage(img, 0, 0, svgW, svgH);
    URL.revokeObjectURL(url);
    if (mimeType === "image/jpeg") {
      canvas.toBlob(
        function (jpgBlob) {
          if (!jpgBlob) {
            showToast("JPG render failed.", "danger");
            return;
          }
          downloadBlob(jpgBlob, "phytclust_tree.jpg");
          showToast("Downloaded JPG at " + dpi + " DPI", "success", 2000);
        },
        "image/jpeg",
        0.95,
      );
      return;
    }
    canvas.toBlob(function (pngBlob) {
      if (!pngBlob) {
        showToast(fmtLabel + " copy failed.", "danger");
        return;
      }
      try {
        var item = new ClipboardItem({ "image/png": pngBlob });
        navigator.clipboard.write([item]).then(
          function () {
            showToast(
              fmtLabel + " copied to clipboard (" + dpi + " DPI)",
              "success",
              2000,
            );
          },
          function () {
            showToast(
              "Clipboard write failed — try the Export menu instead.",
              "danger",
            );
          },
        );
      } catch (e) {
        showToast(
          "Clipboard API not supported — try the Export menu instead.",
          "danger",
        );
      }
    }, "image/png");
  };
  img.onerror = function () {
    URL.revokeObjectURL(url);
    showToast(fmtLabel + " rendering failed.", "danger");
  };
  img.src = url;
}

export function copySvgToClipboard(selector, opts = {}) {
  let snap = null;
  let svgNode = null;
  withExportPreset(opts.usePreset, () => {
    svgNode = document.querySelector(selector + " svg");
    if (!svgNode) return;
    snap = getTightSvgSnapshot(svgNode, 10);
  }, opts.redraw);
  if (!svgNode) {
    showToast("Nothing to copy.", "danger");
    return;
  }
  if (!snap) {
    showToast("Copy failed", "danger");
    return;
  }
  navigator.clipboard.writeText(snap.svgData).then(
    () => showToast("SVG copied to clipboard", "success", 2000),
    () => showToast("Copy failed", "danger"),
  );
}

export async function exportTSV() {
  try {
    const res = await fetch("/api/export_tsv", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ top_n: 1, outlier: true, run_id: state.latestRunId }),
    });
    if (!res.ok) {
      const err = await res.json().catch(() => ({}));
      throw new Error(err.detail || "Export failed");
    }
    downloadBlob(
      new Blob([await res.text()], { type: "text/tab-separated-values" }),
      "phytclust_results.tsv",
    );
    showToast("Downloaded TSV", "success", 2000);
  } catch (e) {
    showToast("TSV export failed: " + e.message, "danger");
  }
}

export async function saveToServer() {
  var el = document.getElementById("output-dir");
  const dir = el && el.value ? el.value : "results";
  try {
    const res = await fetch("/api/save", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        results_dir: dir,
        top_n: 1,
        outlier: true,
      }),
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
    showToast("Saved to " + dir, "success", 3000);
  } catch (e) {
    showToast("Save failed: " + e.message, "danger");
  }
}
