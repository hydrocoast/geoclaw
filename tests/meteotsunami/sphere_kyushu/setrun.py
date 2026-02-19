# encoding: utf-8
"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import datetime
import shutil
import gzip
from pathlib import Path

import numpy as np

import clawpack.clawutil as clawutil
from clawpack.geoclaw import fgout_tools
from clawpack.geoclaw import fgmax_tools
from clawpack.geoclaw import topotools
from clawpack.geoclaw.data import ForceDry
from clawpack.geoclaw.surge.storm import Storm

# Time Conversions
def days2seconds(days):
    return days * 60.0**2 * 24.0

def seconds2days(seconds):
    return seconds / (60.0**2 * 24.0)

# Scratch directory for storing topo and dtopo files:
scratch_dir = (Path(os.environ["CLAW"]) / "geoclaw" / "scratch").resolve()

# ------------------------------
def setrun(claw_pkg='geoclaw'):

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = clawutil.data.ClawRunData(claw_pkg, num_dim)

    # ------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    # ------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # Set single grid parameters first.
    # See below for AMR parameters.

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = 122.0    # west longitude
    clawdata.upper[0] = 132.0    # east longitude

    clawdata.lower[1] =  25.0    # south latitude
    clawdata.upper[1] =  35.0    # north latitude

    # Number of grid cells:
    clawdata.num_cells[0] = 200
    clawdata.num_cells[1] = 200

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    # First three are from shallow GeoClaw, fourth is friction and last 3 are
    # storm fields
    clawdata.num_aux = 3 + 1 + 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    # -------------
    # Initial time:
    # -------------
    clawdata.t0 = 0.0

    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data

    # -------------
    # Output times:
    # --------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 2

    if clawdata.output_style == 1:
        # Output nout frames at equally spaced times up to tfinal:
        # clawdata.tfinal = days2seconds(date2days('2008091400'))
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        #clawdata.output_times = [0.5, 1.0]
        clawdata.output_times = [i*600.0 for i in range(0,37)] # every 10 min, 0 to 360 min
        clawdata.tfinal = 3600.0*6.0

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = True

    clawdata.output_format = 'ascii'      # 'ascii' or 'binary'
    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'all'
    clawdata.output_aux_onlyonce = False    # output aux arrays only at t0

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 1.0e0

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 1000000

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'
    #      ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov'
    #      ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'
    #      ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'

    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0 

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif clawdata.checkpt_style == 2:
        # Specify a list of checkpoint times.
        clawdata.checkpt_times = [0.1, 0.15]

    elif clawdata.checkpt_style == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 2

    # List of refinement ratios at each level (length at least mxnest-1)
    amrdata.refinement_ratios_x = [3, 4]
    amrdata.refinement_ratios_y = [3, 4]
    amrdata.refinement_ratios_t = [3, 4]

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center', 'capacity', 'yleft', 'center', 'center',
                        'center', 'center']

    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0

    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting

    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # == setregions.data values ==
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    gauges = rundata.gaugedata.gauges
    gauges.append([1, 128.0, 30.0, rundata.clawdata.t0, rundata.clawdata.tfinal])
    gauges.append([2, 129.0, 31.0, rundata.clawdata.t0, rundata.clawdata.tfinal])
    gauges.append([3, 130.0, 32.0, rundata.clawdata.t0, rundata.clawdata.tfinal])

    # ------------------------------------------------------------------
    # GeoClaw specific parameters:
    # ------------------------------------------------------------------
    rundata = setgeo(rundata)

    return rundata
    # end of function setrun
    # ----------------------


# -------------------
def setgeo(rundata):
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3
    geo_data.rho = 1025.0
    geo_data.rho_air = 1.15
    geo_data.ambient_pressure = 101.3e3

    # == Forcing Options
    geo_data.coriolis_forcing = True
    geo_data.friction_forcing = True
    geo_data.friction_depth = 1e10

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-2

    # Refinement Criteria
    refine_data = rundata.refinement_data
    refine_data.wave_tolerance = 5.0e-2
    refine_data.speed_tolerance = [1.0, 2.0, 3.0, 4.0]
    refine_data.variable_dt_refinement_ratios = True

    # == settopo.data values ==
    topo_data = rundata.topo_data
    topo_data.topofiles = []
    # for topography, append lines of the form
    #   [topotype, fname]
    topofile = "gebco_2025_n35_s25_w122_e132.asc"
    topo_data.topofiles.append([3,os.path.join(scratch_dir,topofile)])

    # == fgout grids ==
    # new style as of v5.9.0 (old rundata.fixed_grid_data is deprecated)
    # set rundata.fgout_data.fgout_grids to be a 
    # list of objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:
    #rundata.fgout_data.fgout_grids = []

    # Fixed grid output
    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially
    ### fgout 1
    fgout = fgout_tools.FGoutGrid()
    fgout.fgno = 1
    fgout.output_format = 'binary'
    fgout.nx = rundata.clawdata.num_cells[0]
    fgout.ny = rundata.clawdata.num_cells[1]
    fgout.x1 = rundata.clawdata.lower[0]
    fgout.x2 = rundata.clawdata.upper[0]
    fgout.y1 = rundata.clawdata.lower[1]
    fgout.y2 = rundata.clawdata.upper[1]
    fgout.tstart = rundata.clawdata.t0
    fgout.tend = rundata.clawdata.tfinal
    fgout.nout = int((fgout.tend - fgout.tstart)/3600.0) * 2 + 1 # every 30 min

    # ============================
    # == fgmax.data values =======
    # ============================

    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 2  # Save depth and speed
    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Entire
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2  # uniform rectangular x-y grid
    fg.x1 = rundata.clawdata.lower[0]
    fg.x2 = rundata.clawdata.upper[0]
    fg.y1 = rundata.clawdata.lower[1]
    fg.y2 = rundata.clawdata.upper[1]
    fg.dx = 1.0/60.0e0
    fg.tstart_max = rundata.clawdata.t0     # when to start monitoring max values
    fg.tend_max = rundata.clawdata.tfinal   # when to stop monitoring max values
    fg.dt_check = 5.         # target time (sec) increment between updating max values
    fg.min_level_check = 1    # which levels to monitor max on
    fg.arrival_tol = 5.e-2    # tolerance for flagging arrival
    fg.interp_method = 0      # 0 ==> pw const in cells, recommended
    fgmax_grids.append(fg)    # written to fgmax_grids.data

    # ================
    #  Set Surge Data
    # ================
    data = rundata.surge_data

    # Source term controls - These are currently not respected
    data.wind_forcing = False
    data.drag_law = 1
    data.pressure_forcing = True

    data.display_landfall_time = False

    # AMR parameters, m/s and m respectively
    data.wind_refine = [999.99e3,999.99e4]
    data.R_refine = [999.99e3,999.99e4]

    # Storm parameters - Parameterized storm (Holland 1980)
    data.storm_specification_type = 10
    #data.storm_file = (Path() / 'pres.data').resolve()
    data.storm_file = os.path.join(os.getcwd(), 'param_airpressure.txt')

    # =======================
    #  Set Variable Friction
    # =======================
    data = rundata.friction_data

    # Variable friction
    data.variable_friction = False

    # Region based friction
    # Entire domain
    data.friction_regions.append([rundata.clawdata.lower,
                                  rundata.clawdata.upper,
                                  [np.inf, 0.0, -np.inf],
                                  [0.030, 0.022]])

    return rundata
    # end of function setgeo
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
        rundata = setrun()

    rundata.write()
