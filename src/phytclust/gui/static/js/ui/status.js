// Status messages and the current session’s activity log.

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

/** Add an activity-log entry without changing the status bar or showing a toast. */
export function logStatus(message, type) {
  LOG.push({ ts: Date.now(), message, type: type || "info" });
  if (LOG.length > LOG_MAX) LOG.shift();
  window.dispatchEvent(new Event("phytclust:log-update"));
}

export function showStatus(message, type) {
  type = type || "info";
  logStatus(message, type);

  if (statusEl) statusEl.textContent = message;
  if (statusDot) {
    statusDot.className = "status-dot";
    if (type === "info") statusDot.classList.add("running");
    if (type === "success") statusDot.classList.add("ready");
    if (type === "danger") statusDot.classList.add("error");
    if (type === "warning") statusDot.classList.add("warning");
  }
  if (type === "danger") showToast(message, "danger", 6000);
  else if (type === "warning") showToast(message, "warning", 4000);
}
