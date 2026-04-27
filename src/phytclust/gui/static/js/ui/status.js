/* ============================================================
   PhytClust – ui/status.js
   Status bar updates + persistent activity log.
   ============================================================ */

import { statusEl, statusDot } from "../dom.js";
import { showToast } from "./toast.js";

const LOG = [];
const LOG_MAX = 50;

export function getStatusLog() {
  return LOG.slice();
}

export function clearStatusLog() {
  LOG.length = 0;
}

export function showStatus(message, type) {
  type = type || "info";
  LOG.push({ ts: Date.now(), message, type });
  if (LOG.length > LOG_MAX) LOG.shift();
  window.dispatchEvent(new Event("phytclust:log-update"));

  if (statusEl) statusEl.textContent = message;
  if (statusDot) {
    statusDot.className = "status-dot";
    if (type === "info") statusDot.classList.add("running");
    if (type === "success") statusDot.classList.add("ready");
    if (type === "danger") statusDot.classList.add("error");
  }
  if (type === "danger") showToast(message, "danger", 6000);
  else if (type === "success") showToast(message, "success", 3000);
}
