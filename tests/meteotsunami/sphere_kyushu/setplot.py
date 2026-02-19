"""
Set up figures for the spherical Kyushu meteotsunami test case.
"""

from __future__ import absolute_import


def setplot(plotdata=None):
    """
    Specify what is to be plotted at each frame.
    """

    from clawpack.visclaw import geoplot

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    plotdata.clearfigures()

    # Domain used in test_sphere_kyushu/setrun.py
    xlimits = [122.0, 132.0]
    ylimits = [25.0, 35.0]

    drytol = 1.0e-3
    ambient_pressure = 101.3e3
    pressure_index = 6  # aux(7) in Fortran -> aux[6] in Python

    def beforeframe(current_data):
        current_data.user["drytol"] = drytol

    plotdata.beforeframe = beforeframe

    def title_with_time(cd, title):
        from pylab import title
        title(f"t = {cd.t/60.0:.1f} min")

    # ------------------------------------------------------------------
    # Figure 0: Surface elevation + topography
    # ------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name="surface", figno=0)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.scaled = True
    plotaxes.afteraxes = lambda cd: title_with_time(cd, "Surface elevation")

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var = geoplot.surface
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.2
    plotitem.pcolor_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0, 0, 0]
    plotitem.patchedges_show = 0

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 3000.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0, 0, 0]
    plotitem.patchedges_show = 0

    # ------------------------------------------------------------------
    # Figure 1: Pressure anomaly field (aux pressure - ambient)
    # ------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name="pressure_anomaly", figno=1)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.scaled = True
    plotaxes.afteraxes = lambda cd: title_with_time(cd, "Atmospheric pressure anomaly")

    def pressure_anomaly(cd):
        return cd.aux[pressure_index, :, :] - ambient_pressure

    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var = pressure_anomaly
    plotitem.pcolor_cmap = "bwr"
    plotitem.pcolor_cmin = -200.0
    plotitem.pcolor_cmax = 200.0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0, 0, 0]
    plotitem.patchedges_show = 0

    # ------------------------------------------------------------------
    # Figure 300: gauge time series (surface)
    # ------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name="gauges", figno=300, type="each_gauge")
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = "auto"
    plotaxes.ylimits = "auto"
    plotaxes.title = "Gauge surface elevation"

    def gauge_afteraxes(cd):
        from pylab import grid, title, xlabel, ylabel
        xlabel("time (s)")
        ylabel("surface elevation (m)")
        title(f"Gauge {cd.gaugeno}")
        grid(True)

    plotaxes.afteraxes = gauge_afteraxes

    plotitem = plotaxes.new_plotitem(plot_type="1d_plot")
    # gauges: q[0,:]=h, q[3,:]=eta in GeoClaw gauge output convention
    plotitem.plot_var = 3
    plotitem.plotstyle = "b-"
    plotitem.kwargs = {"linewidth": 1.5}

    # ------------------------------------------------------------------
    # Print / html options
    # ------------------------------------------------------------------
    plotdata.printfigs = True
    plotdata.print_format = "png"
    plotdata.print_framenos = "all"
    plotdata.print_gaugenos = "all"
    plotdata.print_fignos = "all"
    plotdata.html = True
    plotdata.html_homelink = "../README.html"
    plotdata.latex = False
    plotdata.parallel = True

    return plotdata

