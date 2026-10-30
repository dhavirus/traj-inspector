"""Draw the 256 px application icon (a trace with a state threshold) with matplotlib: python make_icon.py OUT.png"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

fig = plt.figure(figsize=(2.56, 2.56), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_axis_off()
fig.patch.set_facecolor("#0b3d5c")
x = np.linspace(0, 1, 200)
ax.plot(0.08 + 0.84 * x, 0.5 + 0.28 * np.sign(np.sin(14 * x)) * 0.8 + 0.04 * np.sin(90 * x), color="#f2f2f2", lw=3)
ax.axhline(0.5, color="#e69f00", lw=2.5, xmin=0.08, xmax=0.92)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
fig.savefig(sys.argv[1], dpi=100)
