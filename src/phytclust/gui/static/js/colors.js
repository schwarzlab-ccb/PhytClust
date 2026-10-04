/* ============================================================
   PhytClust – colors.js
   Color palette, generation, and theme helpers.
   ============================================================ */

import { state } from "./state.js";

// Cluster colours share their order and expansion with viz/palette.py.
export const BASE_COLORS = [
  "#b84b4b", // red
  "#4f8f4a", // green
  "#da63aa", // pink
  "#ceb94b", // gold
  "#3f408a", // indigo
  "#c06f2e", // amber brown
  "#5b6bb3", // slate blue
  "#849060", // olive
  "#6e3f8a", // purple
  "#2f8a85", // cyan teal
  "#8f4b7f", // magenta plum
  "#3d7c74", // teal
  "#ad5c7a", // rose
  "#2f6f93", // ocean blue
  "#7a5d3b", // earth
  "#3f648a", // steel blue
];

export function shuffle(arr) {
  let a = arr.slice();
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
}

export function adjustLight(hex, fraction) {
  const value = parseInt(hex.slice(1), 16);
  const blend = (channel) => channel + (255 - channel) * fraction;
  return `rgb(${blend(value >> 16)}, ${blend((value >> 8) & 0xff)}, ${blend(value & 0xff)})`;
}

export function withAlpha(color, alpha) {
  if (color.startsWith("rgb"))
    return color.replace("rgb", "rgba").replace(")", `, ${alpha})`);
  const num = parseInt(color.slice(1), 16);
  return `rgba(${num >> 16}, ${(num >> 8) & 0xff}, ${num & 0xff}, ${alpha})`;
}

export function generateClusterColors(nClusters) {
  if (!Number.isInteger(nClusters) || nClusters < 0) {
    throw new RangeError("Colour count must be an integer zero or greater.");
  }
  let palette = BASE_COLORS.slice();
  if (nClusters > palette.length * 7) {
    console.warn(`Colours repeat beyond ${palette.length * 7} entries. Use cluster labels for larger partitions.`);
  }
  if (nClusters <= palette.length) return palette.slice(0, nClusters);
  let colors = [];
  const minAlpha = 0.58;
  const alphaStep = 0.12;
  const lightStep = 0.14;
  const repeats = Math.ceil(nClusters / palette.length);
  for (let r = 0; r < repeats; r++) {
    const factor = Math.min(r * lightStep, 0.84);
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
