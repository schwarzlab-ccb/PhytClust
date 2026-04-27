/* ============================================================
   PhytClust – ui/appearance_panel.js
   A single declarative panel for all visual options.
   Reads/writes state.render via setRenderOption + getRenderOption,
   triggers a redraw on every change.
   ============================================================ */

import { state, setRenderOption, getRenderOption } from "../state.js";

/**
 * Declarative spec of every appearance control. To add a new option:
 *   1. Add the field to state.render in state.js.
 *   2. Add a row here.
 * The builder handles DOM, events, and initial value sync.
 */
const APPEARANCE_SPEC = [
  {
    section: "Labels",
    controls: [
      { kind: "checkbox", path: "labels.show", label: "Show leaf names" },
      { kind: "checkbox", path: "labels.internalShow", label: "Show internal names" },
      { kind: "number", path: "labels.fontSize", label: "Label size", min: 6, max: 36, step: 1 },
      { kind: "number", path: "labels.axisFontSize", label: "Axis size", min: 6, max: 18, step: 1 },
    ],
  },
  {
    section: "Branches",
    controls: [
      { kind: "number", path: "branches.width", label: "Branch width", min: 0.5, max: 6, step: 0.1 },
      { kind: "color", path: "branches.color", label: "Branch color", nullable: true },
      { kind: "checkbox", path: "branches.colorByClusters", label: "Colour by cluster" },
    ],
  },
  {
    section: "Nodes",
    controls: [
      { kind: "number", path: "nodes.leafRadius", label: "Leaf node size", min: 1, max: 10, step: 0.5 },
    ],
  },
  {
    section: "Cluster boxes",
    controls: [
      { kind: "checkbox", path: "clusters.showBoxLabels", label: "Show box labels" },
      { kind: "checkbox", path: "clusters.showOutlierBoxes", label: "Show outlier boxes" },
      { kind: "range", path: "clusters.boxAlpha", label: "Fill opacity", min: 0.02, max: 0.5, step: 0.02 },
      { kind: "range", path: "clusters.boxPadV", label: "Vertical padding", min: 0, max: 30, step: 1 },
      { kind: "range", path: "clusters.boxPadH", label: "Horizontal padding", min: 0, max: 30, step: 1 },
      { kind: "range", path: "clusters.boxCornerRadius", label: "Corner radius", min: 0, max: 20, step: 1 },
      { kind: "number", path: "clusters.labelFontSize", label: "Box label size", min: 6, max: 36, step: 1 },
    ],
  },
  {
    section: "Tree scale",
    controls: [
      { kind: "range", path: "scale.width", label: "Width", min: 0.1, max: 8.0, step: 0.1 },
      { kind: "range", path: "scale.height", label: "Height", min: 0.1, max: 8.0, step: 0.1 },
    ],
  },
  {
    section: "Comparison bars",
    controls: [
      { kind: "range", path: "compare.barWidth", label: "Bar width", min: 4, max: 60, step: 1 },
      { kind: "range", path: "compare.barGap", label: "Bar gap", min: 0, max: 24, step: 1 },
      { kind: "checkbox", path: "compare.showColumnTitles", label: "Show bar titles" },
    ],
  },
];

function buildControl(ctrl, onChange) {
  const row = document.createElement("div");
  row.className = "ap-row ap-row-" + ctrl.kind;

  const label = document.createElement("label");
  label.className = "ap-label";
  label.textContent = ctrl.label;

  let input;
  const current = getRenderOption(ctrl.path);

  switch (ctrl.kind) {
    case "checkbox": {
      input = document.createElement("input");
      input.type = "checkbox";
      input.checked = !!current;
      label.prepend(input);
      input.addEventListener("change", () => {
        setRenderOption(ctrl.path, input.checked);
        onChange();
      });
      row.appendChild(label);
      break;
    }
    case "number": {
      input = document.createElement("input");
      input.type = "number";
      if (ctrl.min != null) input.min = ctrl.min;
      if (ctrl.max != null) input.max = ctrl.max;
      if (ctrl.step != null) input.step = ctrl.step;
      input.value = current;
      input.addEventListener("input", () => {
        const v = parseFloat(input.value);
        if (Number.isFinite(v)) {
          setRenderOption(ctrl.path, v);
          onChange();
        }
      });
      row.appendChild(label);
      row.appendChild(input);
      break;
    }
    case "range": {
      input = document.createElement("input");
      input.type = "range";
      input.min = ctrl.min;
      input.max = ctrl.max;
      input.step = ctrl.step;
      input.value = current;
      const valueEcho = document.createElement("span");
      valueEcho.className = "ap-range-value";
      valueEcho.textContent = current;
      input.addEventListener("input", () => {
        const v = parseFloat(input.value);
        if (Number.isFinite(v)) {
          setRenderOption(ctrl.path, v);
          valueEcho.textContent = v;
          onChange();
        }
      });
      row.appendChild(label);
      row.appendChild(input);
      row.appendChild(valueEcho);
      break;
    }
    case "color": {
      input = document.createElement("input");
      input.type = "color";
      input.value = current || "#000000";
      input.addEventListener("input", () => {
        setRenderOption(ctrl.path, input.value);
        onChange();
      });
      row.appendChild(label);
      row.appendChild(input);
      if (ctrl.nullable) {
        const reset = document.createElement("button");
        reset.type = "button";
        reset.className = "ap-reset";
        reset.textContent = "Reset";
        reset.title = "Use cluster colour";
        reset.addEventListener("click", () => {
          setRenderOption(ctrl.path, null);
          onChange();
        });
        row.appendChild(reset);
      }
      break;
    }
    default:
      console.warn("Unknown control kind:", ctrl.kind);
  }

  return row;
}

/**
 * Build the panel DOM into `container` and wire all controls.
 * `onChange` is invoked after every option write — typically `drawTree`.
 */
export function mountAppearancePanel(container, onChange) {
  if (!container) return;
  container.innerHTML = "";
  container.classList.add("appearance-panel");

  for (const block of APPEARANCE_SPEC) {
    const section = document.createElement("section");
    section.className = "ap-section";

    const title = document.createElement("div");
    title.className = "ap-section-title";
    title.textContent = block.section;
    section.appendChild(title);

    for (const ctrl of block.controls) {
      section.appendChild(buildControl(ctrl, onChange));
    }
    container.appendChild(section);
  }
}

/**
 * Wire one or more trigger buttons to the same shared panel. Clicking any
 * trigger toggles visibility; clicking outside closes the panel.
 *
 * @param {Object} opts
 * @param {string|string[]} opts.buttonId  — id, or array of ids, of trigger buttons
 * @param {string} opts.panelId
 * @param {Function} opts.onChange         — invoked on every option write
 */
export function wireAppearanceTrigger({ buttonId, panelId, onChange }) {
  const panel = document.getElementById(panelId);
  if (!panel) return;
  const ids = Array.isArray(buttonId) ? buttonId : [buttonId];
  const buttons = ids.map((id) => document.getElementById(id)).filter(Boolean);
  if (!buttons.length) return;

  const setExpanded = (val) => {
    buttons.forEach((b) => b.setAttribute("aria-expanded", String(val)));
  };
  const close = () => {
    panel.classList.remove("open");
    setExpanded(false);
  };
  const open = () => {
    // Re-mount on every open so the controls reflect any external state
    // changes (e.g. Reset View, presets, programmatic toggles) since the
    // panel was last shown.
    mountAppearancePanel(panel, onChange);
    panel.classList.add("open");
    setExpanded(true);
  };

  buttons.forEach((btn) => {
    btn.addEventListener("click", (e) => {
      e.stopPropagation();
      panel.classList.contains("open") ? close() : open();
    });
  });

  document.addEventListener("click", (e) => {
    if (panel.contains(e.target)) return;
    if (buttons.some((b) => b === e.target || b.contains(e.target))) return;
    close();
  });

  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") close();
  });
}
