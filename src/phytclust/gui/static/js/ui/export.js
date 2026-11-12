// Download figures, assignments, and copied images.

import { state, withRenderOverrides } from "../state.js";
import { drawTree } from "../tree/draw.js";
import { showToast } from "./toast.js";
import { markOutlierClusters } from "../utils.js";
import { outlierMarkingSettings } from "../views/cluster_editor.js";
import { postJson } from "../api.js";

/** Larger labels, thicker branches, and a stronger cluster-box fill for exports. */
const PUBLICATION_PRESET = {
  branches: { width: 2.0 },
  labels: { fontSize: 12, axisFontSize: 12 },
  nodes: { leafRadius: 4.0 },
  clusters: { boxAlpha: 0.25, labelFontSize: 11 },
};

/**
  * Run `fn` with the tree redrawn under PUBLICATION_PRESET (none if
  * `usePreset` is false), then restore. `redraw` defaults to drawTree;
  * Compare passes drawComparison.
 */
function withExportPreset(usePreset, fn, redraw = drawTree) {
  if (!usePreset) return fn();
  return withRenderOverrides(PUBLICATION_PRESET, redraw, fn);
}

function downloadBlob(blob, filename) {
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

function getTightSvgSnapshot(svgNode, paddingPx) {
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

/**
 * Snapshot the SVG under `selector` and paint it onto a canvas at `dpi`.
 * `onCanvas(canvas)` runs once the image has rendered; `onError(message)`
 * reports why no canvas could be produced.
 */
function rasterize(selector, dpi, opts, { noSvgMessage, failMessage, renderFailMessage, onCanvas }) {
  let snap = null;
  let svgNode = null;
  withExportPreset(opts && opts.usePreset, () => {
    svgNode = document.querySelector(selector + " svg");
    if (!svgNode) return;
    snap = getTightSvgSnapshot(svgNode, 10);
  }, opts && opts.redraw);
  if (!svgNode) {
    showToast(noSvgMessage, "danger");
    return;
  }
  if (!snap) {
    showToast(failMessage, "danger");
    return;
  }
  const scale = dpi / 96;
  const svgW = snap.width;
  const svgH = snap.height;
  const canvas = document.createElement("canvas");
  canvas.width = Math.ceil(svgW * scale);
  canvas.height = Math.ceil(svgH * scale);
  const ctx = canvas.getContext("2d");
  // Same background as the page, so theme-coloured branches and labels stay visible.
  ctx.fillStyle =
    getComputedStyle(document.documentElement)
      .getPropertyValue("--pc-bg")
      .trim() || "#ffffff";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.scale(scale, scale);
  const img = new Image();
  const blob = new Blob([snap.svgData], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  img.onload = function () {
    ctx.drawImage(img, 0, 0, svgW, svgH);
    URL.revokeObjectURL(url);
    onCanvas(canvas);
  };
  img.onerror = function () {
    URL.revokeObjectURL(url);
    showToast(renderFailMessage, "danger");
  };
  img.src = url;
}

export function exportPngFromEl(selector, filename, dpi, opts = {}) {
  dpi = dpi || 300;
  rasterize(selector, dpi, opts, {
    noSvgMessage: "No visualization to export.",
    failMessage: "PNG export failed.",
    renderFailMessage: "PNG rendering failed.",
    onCanvas(canvas) {
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
    },
  });
}

export function copyRasterToClipboard(selector, dpi, mimeType, opts = {}) {
  dpi = dpi || 300;
  mimeType = mimeType || "image/png";
  var fmtLabel = mimeType === "image/jpeg" ? "JPG" : "PNG";
  rasterize(selector, dpi, opts, {
    noSvgMessage: "Nothing to copy.",
    failMessage: "Copy failed.",
    renderFailMessage: fmtLabel + " rendering failed.",
    onCanvas(canvas) {
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
    },
  });
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

/**
  * The displayed assignments as TSV, with the same outlier marking as the tree.
 */
export function currentTSV() {
  const raw = state.runClusters || {};
  const names = Object.keys(raw);
  if (!names.length) {
    return null;
  }
  const displayed = markOutlierClusters(raw, outlierMarkingSettings());
  const k = new Set(Object.values(raw)).size;
  const rows = names.map((name) => name + "\t" + displayed[name]);
  return {
    text: "Node Name\tclusters_k" + k + "\n" + rows.join("\n") + "\n",
    filename: `phytclust_k${k}.tsv`,
    k,
  };
}

export function exportTSV() {
  const current = currentTSV();
  if (!current) {
    showToast("Nothing to export yet. Run PhytClust first.", "danger");
    return;
  }
  downloadBlob(
    new Blob([current.text], { type: "text/tab-separated-values" }),
    current.filename,
  );
  showToast(`Downloaded k=${current.k} TSV`, "success", 2000);
}

export async function saveToServer() {
  const current = currentTSV();
  if (!current) {
    showToast("Nothing to save yet. Run PhytClust first.", "danger");
    return;
  }
  var el = document.getElementById("output-dir");
  const dir = el && el.value ? el.value : "results";
  try {
    const data = await postJson("/api/save", {
      results_dir: dir,
      filename: current.filename,
      tsv: current.text,
    });
    const savedDir = data && typeof data.results_dir === "string" ? data.results_dir : dir;
    const savedName = data && typeof data.filename === "string" ? data.filename : current.filename;
    showToast("Saved " + savedName + " to " + savedDir, "success", 4000);
  } catch (e) {
    showToast("Save failed: " + e.message, "danger");
  }
}
