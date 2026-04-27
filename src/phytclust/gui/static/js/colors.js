/* ============================================================
   PhytClust – colors.js
   Color palette, generation, and theme helpers.
   ============================================================ */

import { state } from "./state.js";

export const BASE_COLORS = [
  "#b84b4b",
  "#849060",
  "#3d7c74",
  "#6e3f8a",
  "#ceb94b",
  "#3f648a",
  "#3f408a",
  "#da63aa",
  "#c06f2e",
  "#2f6f93",
  "#4f8f4a",
  "#ad5c7a",
  "#7a5d3b",
  "#2f8a85",
  "#8f4b7f",
  "#5b6bb3",
];

export function shuffle(arr) {
  let a = arr.slice();
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
}

export function adjustLight(hex, factor) {
  const num = parseInt(hex.slice(1), 16);
  let r = (num >> 16) + factor * 255;
  let g = ((num >> 8) & 0xff) + factor * 255;
  let b = (num & 0xff) + factor * 255;
  return `rgb(${Math.min(255, Math.max(0, r)) | 0}, ${Math.min(255, Math.max(0, g)) | 0}, ${Math.min(255, Math.max(0, b)) | 0})`;
}

export function withAlpha(color, alpha) {
  if (color.startsWith("rgb"))
    return color.replace("rgb", "rgba").replace(")", `, ${alpha})`);
  const num = parseInt(color.slice(1), 16);
  return `rgba(${num >> 16}, ${(num >> 8) & 0xff}, ${num & 0xff}, ${alpha})`;
}

export function generateClusterColors(nClusters) {
  let palette = BASE_COLORS.slice();
  if (nClusters <= palette.length) return palette.slice(0, nClusters);
  let colors = [];
  const minAlpha = 0.58;
  const alphaStep = 0.12;
  const lightStep = 0.14;
  const repeats = Math.ceil(nClusters / palette.length);
  for (let r = 0; r < repeats; r++) {
    const factor = r * lightStep;
    const alpha = Math.max(1 - r * alphaStep, minAlpha);
    palette.forEach((hex) => {
      const adjusted = adjustLight(hex, factor);
      colors.push(alpha < 1 ? withAlpha(adjusted, alpha) : adjusted);
    });
  }
  return colors.slice(0, nClusters);
}

export function getThemeColors() {
  const style = getComputedStyle(document.documentElement);
  const themeBranch =
    style.getPropertyValue("--pc-tree-branch").trim() || "#000000";
  return {
    branch: state.render.branches.color || themeBranch,
    label: style.getPropertyValue("--pc-tree-label").trim() || "#334155",
    internal: style.getPropertyValue("--pc-tree-internal").trim() || "#64748b",
    text: style.getPropertyValue("--pc-text").trim() || "#241b3d",
    muted: style.getPropertyValue("--pc-text-muted").trim() || "#9a91b2",
    bg: style.getPropertyValue("--pc-bg").trim() || "#f7f5fa",
    border: style.getPropertyValue("--pc-border").trim() || "#e4ddf0",
    accent: style.getPropertyValue("--pc-accent").trim() || "#6b61ac",
  };
}

const THEME_KEY = "phytclust-theme";

function preferredTheme() {
  // Only honour an explicit user choice. We deliberately ignore the OS
  // prefers-color-scheme so a user on a dark-mode system still sees the
  // intended light look until they opt in via the theme toggle.
  const stored = localStorage.getItem(THEME_KEY);
  if (stored === "light" || stored === "dark") return stored;
  return "light";
}

export function initTheme() {
  document.documentElement.setAttribute("data-theme", preferredTheme());
}

export function toggleTheme() {
  const current = document.documentElement.getAttribute("data-theme") || "light";
  const next = current === "dark" ? "light" : "dark";
  document.documentElement.setAttribute("data-theme", next);
  localStorage.setItem(THEME_KEY, next);
  return next;
}
