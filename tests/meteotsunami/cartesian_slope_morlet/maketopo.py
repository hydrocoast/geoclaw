# encoding: utf-8
"""
Make topofiles for the Cartesian slope case using topotype=1.
"""

import numpy as np
import matplotlib.pyplot as plt


def write_topo1_rectangle(x1, x2, y1, y2, z_nw, z_ne, z_sw, z_se, outfile):
    """
    Write a topotype=1 file using the four rectangle corners.
    """

    corners = [
        (x1, y2, z_nw),  # NW
        (x2, y2, z_ne),  # NE
        (x1, y1, z_sw),  # SW
        (x2, y1, z_se),  # SE
    ]

    with open(outfile, "w", encoding="ascii") as fp:
        for x, y, z in corners:
            fp.write(f"{x:.12e} {y:.12e} {z:.12e}\n")


def make_3region_slope_topo(
    h1,
    h2,
    s,
    xlower,
    xupper,
    ylower,
    yupper,
    x0,
    outfiles=None,
    plot=True,
):
    """
    Build a 2D topography with 3 regions:
    uniform h1 / slope / uniform h2.

    The h1 region is always on the left side.
    The slope always extends rightward from x0 to x1, and
    can either deepen or shoal depending on h2 - h1.
    The profile is uniform in y (Z depends only on x).

    The topography is written as three adjacent topotype=1 files:
    left flat, slope, and right flat.
    """

    if s <= 0:
        raise ValueError("Assumes s > 0.")

    if outfiles is None:
        outfiles = [
            "region1.tt1",
            "region2.tt1",
            "region3.tt1",
        ]

    if len(outfiles) != 3:
        raise ValueError("Expected exactly three output files.")

    if np.isclose(h1, h2):
        x1 = x0
    else:
        x1 = x0 + abs(h2 - h1) / s

    if x0 < xlower or x0 > xupper:
        raise ValueError("x0 must lie inside the domain.")
    if x1 < x0 or x1 > xupper:
        raise ValueError("The slope end point x1 must satisfy x0 <= x1 <= xupper.")

    z1 = -h1
    z2 = -h2

    write_topo1_rectangle(
        xlower, x0, ylower, yupper,
        z1, z1, z1, z1,
        outfiles[0],
    )
    write_topo1_rectangle(
        x0, x1, ylower, yupper,
        z1, z2, z1, z2,
        outfiles[1],
    )
    write_topo1_rectangle(
        x1, xupper, ylower, yupper,
        z2, z2, z2, z2,
        outfiles[2],
    )

    print("Topofiles written:")
    for outfile in outfiles:
        print(f"  {outfile}")

    trend = "deepening" if h2 > h1 else "shoaling"
    print(f"Slope region: x0 = {x0:.1f} m -> x1 = {x1:.1f} m ({trend} to the right)")



    if plot:
        nx, ny = 2001, 11
        x = np.linspace(xlower, xupper, nx)
        y = np.linspace(ylower, yupper, ny)

        h = np.empty_like(x)
        if np.isclose(x0, x1):
            h[:] = h1
        else:
            left_flat = x <= x0
            slope = (x > x0) & (x < x1)
            right_flat = x >= x1

            h[left_flat] = h1
            h[right_flat] = h2
            h[slope] = h1 + (h2 - h1) * (x[slope] - x0) / (x1 - x0)

        z = -np.tile(h, (ny, 1))
        plot_topography(x, y, z, x0, x1, "topo_slope.png")


def plot_topography(x, y, z, x0, x1, figfile):
    """
    Plot 2D map and cross-shore profile of the topography.
    """

    fig = plt.figure(figsize=(12, 5))

    ax1 = fig.add_subplot(1, 2, 1)
    im = ax1.pcolormesh(x, y, z, shading="auto")
    plt.colorbar(im, ax=ax1, label="Elevation (m)")
    ax1.set_title("Topography (2D)")
    ax1.set_xlabel("x (m)")
    ax1.set_ylabel("y (m)")
    ax1.axvline(x0, color="r", linestyle="--")
    ax1.axvline(x1, color="r", linestyle="--")

    ax2 = fig.add_subplot(1, 2, 2)
    mid = len(y) // 2
    ax2.plot(x, z[mid, :])
    ax2.set_title("Cross-shore profile")
    ax2.set_xlabel("x (m)")
    ax2.set_ylabel("Elevation (m)")
    ax2.axvline(x0, color="r", linestyle="--")
    ax2.axvline(x1, color="r", linestyle="--")
    ax2.grid()

    plt.tight_layout()
    plt.savefig(figfile, dpi=200, bbox_inches="tight")


if __name__ == "__main__":

    h1 = 200.0
    h2 = 2000.0
    s = 0.01

    xlower, xupper = 0.0, 400e3
    ylower, yupper = 0.0, 1000e3
    x0 = 150e3

    make_3region_slope_topo(
        h1,
        h2,
        s,
        xlower,
        xupper,
        ylower,
        yupper,
        x0,
        outfiles=[
            "region1.tt1",
            "region2.tt1",
            "region3.tt1",
        ],
        plot=True,
    )
