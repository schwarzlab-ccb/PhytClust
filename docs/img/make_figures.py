"""Regenerate the cluster-plot figures shown on the reference page.

Run from the repository root:  python docs/img/make_figures.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from phytclust import PhytClust  # noqa: E402
from phytclust.viz.cluster import plot_clusters  # noqa: E402

OUT = Path(__file__).parent
FIGURES = {
    "01_baseline": {},
    "02_bars": {"show_cluster_bars": True},
    "03_boxes": {"show_cluster_boxes": True},
    "04_branches_coloured": {"colour_branches_by_cluster": True},
    "05_cladogram": {"layout": "cladogram", "show_cluster_bars": True},
    "06_publication_style": {
        "layout": "cladogram",
        "colour_branches_by_cluster": True,
        "show_cluster_boxes": True,
    },
    "07_custom_palette": {"palette": ["#264653", "#2a9d8f", "#e9c46a"]},
}

pc = PhytClust("examples/sample_tree.nwk")
for name, options in FIGURES.items():
    plot_clusters(
        pc,
        k=3,
        show_terminal_labels=True,
        save=True,
        results_dir=str(OUT),
        filename=f"{name}.png",
        dpi=100,
        **options,
    )
