// Toggle collapsible sections.

export function wireCollapse(toggleId, bodyId) {
  var toggle = document.getElementById(toggleId);
  var body = document.getElementById(bodyId);
  if (!toggle || !body) return;
  toggle.addEventListener("click", function () {
    var isOpen = body.classList.contains("open");
    body.classList.toggle("open", !isOpen);
    // The chevron turns via CSS on [aria-expanded="true"] (sidebar.css).
    toggle.setAttribute("aria-expanded", String(!isOpen));
  });
}
