# encoding: utf-8
"""
Make a topofile
"""

import numpy as np
import matplotlib.pyplot as plt
from clawpack.geoclaw import topotools


def make_3region_slope_topo(
    h1, h2, s,
    xlower, xupper,
    ylower, yupper,
    nx, ny,
    x0,
    outfile="simple_slope.tt3",
    plot=True,
):
    """
    Build a 2D topography with 3 regions:
    uniform h1 / slope / uniform h2.
    The h1 region is always on the left side.
    The slope always extends rightward from x0 to x1, and
    can either deepen or shoal depending on h2 - h1.
    The profile is uniform in y (Z depends only on x).
    """

    if s <= 0:
        raise ValueError("Assumes s > 0.")

    # If h1 == h2, the topography is uniform
    if np.isclose(h1, h2):
        x1 = x0
    else:
        # Slope end point measured to the right of x0
        x1 = x0 + abs(h2 - h1) / s

    # Grid
    x = np.linspace(xlower, xupper, nx)
    y = np.linspace(ylower, yupper, ny)

    # Water depth
    h = np.empty_like(x)

    if np.isclose(x0, x1):
        h[:] = h1
    else:
        left_flat = x <= x0
        slope = (x > x0) & (x < x1)
        right_flat = x >= x1

        h[left_flat] = h1
        h[right_flat] = h2

        # Linear interpolation: (x0, h1) -> (x1, h2)
        h[slope] = h1 + (h2 - h1) * (x[slope] - x0) / (x1 - x0)

    Z1d = -h
    Z = np.tile(Z1d, (ny, 1))

    # --- topo object ---
    topo = topotools.Topography()
    topo.x = x
    topo.y = y
    topo.Z = Z

    topo.write(outfile, topo_type=3)

    print(f"Topo written to: {outfile}")
    trend = "deepening" if h2 > h1 else "shoaling"
    print(f"Slope region: x0 = {x0:.1f} m -> x1 = {x1:.1f} m ({trend} to the right)")

    # --- Plot ---
    if plot:
        figfile = outfile.replace(".tt3", ".png")
        plot_topography(x, y, Z, x0, x1, figfile)

    return topo


def plot_topography(x, y, Z, x0, x1, figfile):
    """
    plot 2D map and cross-shore profile of the topography.
    """

    fig = plt.figure(figsize=(12, 5))

    # --- 2D map ---
    ax1 = fig.add_subplot(1, 2, 1)
    im = ax1.pcolormesh(x, y, Z, shading="auto")
    plt.colorbar(im, ax=ax1, label="Elevation (m)")
    ax1.set_title("Topography (2D)")
    ax1.set_xlabel("x (m)")
    ax1.set_ylabel("y (m)")
    ax1.axvline(x0, color="r", linestyle="--")
    ax1.axvline(x1, color="r", linestyle="--")

    # --- Cross section ---
    ax2 = fig.add_subplot(1, 2, 2)
    mid = len(y) // 2
    ax2.plot(x, Z[mid, :])
    ax2.set_title("Cross-shore profile")
    ax2.set_xlabel("x (m)")
    ax2.set_ylabel("Elevation (m)")
    ax2.axvline(x0, color="r", linestyle="--")
    ax2.axvline(x1, color="r", linestyle="--")
    ax2.grid()

    plt.tight_layout()
    plt.savefig(figfile, dpi=200, bbox_inches="tight")
    #plt.show()


if __name__ == "__main__":

    h1 = 200.0
    h2 = 2000.0
    s = 0.01

    xlower, xupper = 0.0, 400e3
    ylower, yupper = 0.0, 1000e3
    nx, ny = 2001, 11
    x0 = 150e3

    make_3region_slope_topo(
        h1, h2, s,
        xlower, xupper,
        ylower, yupper,
        nx, ny,
        x0,
        outfile="simple_slope.tt3",
        plot=True,
    )
