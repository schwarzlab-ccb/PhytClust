/* ============================================================
   PhytClust – ui/file.js
   File handling: load tree (Newick) or restore session (JSON).
   ============================================================ */

import { newickEl } from "../dom.js";
import { showStatus } from "./status.js";
import { showToast } from "./toast.js";
import { applySession } from "../session.js";

export function handleFileSelect(evt) {
  var file = evt.target.files ? evt.target.files[0] : null;
  if (!file) return;
  loadFile(file);
}

function isJsonFile(file, text) {
  if (file && /\.json$/i.test(file.name || "")) return true;
  if (file && file.type === "application/json") return true;
  // Fall back to content-sniffing — a session JSON will start with {
  return typeof text === "string" && text.trim().startsWith("{");
}

export function loadFile(file) {
  var fileNameEl = document.getElementById("file-name");
  var fileBadge = document.getElementById("file-badge");
  if (fileNameEl && fileBadge) {
    fileBadge.textContent = file.name;
    fileNameEl.style.display = "flex";
  }
  showStatus("Loading...", "info");
  var reader = new FileReader();
  reader.onload = function (e) {
    var text = (e.target.result || "").trim();
    if (isJsonFile(file, text)) {
      try {
        var parsed = JSON.parse(text);
        var res = applySession(parsed);
        if (res.errors && res.errors.length) {
          res.errors.forEach((m) => showToast(m, "info", 4000));
        }
        showStatus(
          "Session restored (" + res.applied + " settings" +
            (res.hadNewick ? " + tree" : "") + ")",
          "success",
        );
      } catch (err) {
        showStatus("Not a valid session JSON: " + err.message, "danger");
      }
      return;
    }
    newickEl.value = text;
    newickEl.dispatchEvent(new Event("change", { bubbles: true }));
    showStatus("Tree loaded", "success");
  };
  reader.readAsText(file);
}
