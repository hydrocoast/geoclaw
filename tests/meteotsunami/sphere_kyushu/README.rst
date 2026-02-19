.. _geoclaw_examples_meteotsunami_sphere_kyushu:

Meteotsunami Test Around Kyushu (Spherical Coordinates)
========================================================

This directory contains a spherical-coordinate GeoClaw test for pressure-forced
meteotsunami waves around Kyushu, Japan.

Main characteristics
--------------------

1. Coordinate system: spherical longitude-latitude,
   ``geo_data.coordinate_system = 2``.
2. Domain: longitude ``[122, 132]``, latitude ``[25, 35]``.
3. Base grid: ``200 x 200`` cells.
4. Topography: GEBCO ASCII file
   ``gebco_2025_n35_s25_w122_e132.asc`` read from ``$CLAW/geoclaw/scratch``.
5. Atmospheric forcing: custom planewave storm type
   ``storm_specification_type = 10`` (see ``param_airpressure.txt``).
6. Wind forcing is disabled; pressure forcing is enabled.

Build and run
-------------

Run everything (data, simulation, plots)::

    make all

Common individual targets::

    make .data
    make .output
    make .plots

Input files
-----------

1. ``setrun.py``: domain, numerics, AMR, gauges, forcing switches.
2. ``param_airpressure.txt``: plane wave pressure parameters.
3. ``setplot.py``: plotting setup for surface elevation, pressure anomaly, and gauges.

Planewave parameter file format
-------------------------------

``param_airpressure.txt`` uses the following order:

1. ``wave_shape``: ``SINE`` or ``HAT``
2. ``wave_amplitude`` [Pa]
3. ``wave_wavelength`` [m]
4. ``wave_count`` [-]
5. ``wave_speed`` [m/s]
6. ``wave_origin`` [lon, lat]
7. ``theta_deg`` [deg]
8. ``wave_cross_width`` [m]

Current case values
-------------------

1. ``wave_shape = SINE``
2. ``wave_amplitude = 200 Pa``
3. ``wave_wavelength = 50 km``
4. ``wave_count = 2.0``
5. ``wave_speed = 40 m/s``
6. ``wave_origin = (124.0, 28.0)``
7. ``theta_deg = 15``
8. ``wave_cross_width = 300 km``

Output and gauges
-----------------

1. Simulation duration: ``6 hours``.
2. Output times: every ``600 s`` (10 minutes), from ``0`` to ``21600 s``.
3. Gauges:

   - Gauge 1: ``(128.0, 30.0)``
   - Gauge 2: ``(129.0, 31.0)``
   - Gauge 3: ``(130.0, 32.0)``

Plotting
--------

``setplot.py`` is tailored to this case and includes:

1. Surface elevation map.
2. Pressure anomaly map (``aux(7) - ambient_pressure``) using ``bwr`` colormap.
3. Gauge time series of surface elevation.
