"""Shared test setup.

Selects a non-interactive matplotlib backend before any test imports pyplot.
Without this the suite blocks on plt.show() wherever a display is available.
"""

import matplotlib

matplotlib.use("Agg")
