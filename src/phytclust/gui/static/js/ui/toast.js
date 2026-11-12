// Show temporary notifications.

export function showToast(message, type, duration) {
  type = type || "info";
  duration = duration || 4000;
  const container = document.getElementById("toast-container");
  if (!container) return;
  const toast = document.createElement("div");
  toast.className = "toast " + type;
  toast.textContent = message;
  container.appendChild(toast);
  setTimeout(function () {
    toast.classList.add("fade-out");
    setTimeout(function () {
      toast.remove();
    }, 300);
  }, duration);
}
