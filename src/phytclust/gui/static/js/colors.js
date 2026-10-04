// Cluster palettes and display themes.

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

/** Return a shuffled copy of the colours. */
export function shuffle(colors) {
  const shuffled = colors.slice();
  for (let index = shuffled.length - 1; index > 0; index--) {
    const otherIndex = Math.floor(Math.random() * (index + 1));
    [shuffled[index], shuffled[otherIndex]] = [shuffled[otherIndex], shuffled[index]];
  }
  return shuffled;
}

function readColor(color) {
  const parsed = typeof color === "string" ? d3.color(color) : null;
  if (!parsed) throw new Error(`Invalid colour: ${color}`);
  const rgb = parsed.rgb();
  // D3 represents transparent channels as NaN.
  if (rgb.opacity === 0) {
    for (const channel of ["r", "g", "b"]) {
      if (Number.isNaN(rgb[channel])) rgb[channel] = 0;
    }
  }
  if (![rgb.r, rgb.g, rgb.b, rgb.opacity].every(Number.isFinite)) {
    throw new Error(`Invalid colour: ${color}`);
  }
  return rgb;
}

function checkFraction(value, name) {
  if (!Number.isFinite(value) || value < 0 || value > 1) {
    throw new RangeError(`${name} must be a number from 0 to 1.`);
  }
}

function lightenedColor(rgb, fraction, opacity = rgb.opacity) {
  const channels = [rgb.r, rgb.g, rgb.b].map(channel => channel + (255 - channel) * fraction);
  return d3.rgb(...channels, opacity).formatRgb();
}

/** Blend a CSS colour toward white, preserving its opacity. */
export function adjustLight(color, fraction) {
  checkFraction(fraction, "Lightening fraction");
  return lightenedColor(readColor(color), fraction);
}

/** Set a CSS colour's opacity without changing its RGB channels. */
export function withAlpha(color, alpha) {
  checkFraction(alpha, "Opacity");
  const rgb = readColor(color);
  return `rgba(${rgb.r}, ${rgb.g}, ${rgb.b}, ${alpha})`;
}

/** Expand the base palette with progressively lighter colours. */
export function generateClusterColors(clusterCount) {
  if (!Number.isSafeInteger(clusterCount) || clusterCount < 0) {
    throw new RangeError("Colour count must be a safe integer zero or greater.");
  }
  if (clusterCount === 0) return [];
  if (!BASE_COLORS.length) throw new Error("The cluster palette is empty.");
  const baseColors = BASE_COLORS.map(readColor);
  if (clusterCount > baseColors.length * 7) {
    console.warn(`Colours repeat beyond ${baseColors.length * 7} entries. Use cluster labels for larger partitions.`);
  }
  const colors = new Array(clusterCount);
  for (let index = 0; index < clusterCount; index++) {
    const baseIndex = index % baseColors.length;
    const roundIndex = Math.floor(index / baseColors.length);
    if (roundIndex === 0) {
      colors[index] = BASE_COLORS[baseIndex];
      continue;
    }
    const fraction = Math.min(roundIndex * 0.14, 0.84);
    const opacity = baseColors[baseIndex].opacity * Math.max(1 - roundIndex * 0.12, 0.58);
    colors[index] = lightenedColor(baseColors[baseIndex], fraction, opacity);
  }
  return colors;
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
  // Start in light mode unless the user saved a different theme.
  try {
    const storedTheme = localStorage.getItem(THEME_KEY);
    if (storedTheme === "light" || storedTheme === "dark") return storedTheme;
  } catch {
    console.warn("Could not read the saved theme preference.");
  }
  return "light";
}

export function initTheme() {
  document.documentElement.setAttribute("data-theme", preferredTheme());
}

export function toggleTheme() {
  const current = document.documentElement.getAttribute("data-theme") || "light";
  const next = current === "dark" ? "light" : "dark";
  document.documentElement.setAttribute("data-theme", next);
  try {
    localStorage.setItem(THEME_KEY, next);
  } catch {
    console.warn("Could not save the theme preference.");
  }
  return next;
}
