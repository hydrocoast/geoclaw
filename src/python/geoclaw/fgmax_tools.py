r"""
fgmax_tools module: $CLAW/geoclaw/src/python/geoclaw/fgmax_tools.py

Tools to specify an fgmax grid for keeping track of maximum flow depth, etc.
and to read in the fgmax output after doing a GeoClaw run.

"""

from __future__ import absolute_import
from __future__ import print_function
import os
from numpy import sqrt, ma
import numpy
from six.moves import range


class FGmaxGrid(object):

    """
    New class introduced in 5.2.1 to keep store information both about the
    fgmax input data and the output generated by a GeoClaw run.
    """

    def __init__(self):

        super(FGmaxGrid, self).__init__()

        # GeoClaw input values:
        self.id = ''  # identifier, optional
        self.point_style = None
        self.npts = None
        self.nx = None
        self.ny = None
        self.n12 = None
        self.n23 = None
        self.tstart_max =  0.
        self.tend_max = 1.e10   # when to stop monitoring max values
        self.dt_check = 10.     # target time (sec) increment between updating
                                # max values
        self.min_level_check = None    # which levels to monitor max on
        self.interp_method = 0    # 0 for pw const, 1 for bilinear
        self.arrival_tol = 1.e-2       # tolerance for flagging arrival
        #self.input_file_name = 'fgmax.txt'  # file for GeoClaw input data
        self.fgno = None  # FG number
        self.xy_fname = None   # optional file name for separate list of points
                             # when point_style==0, distinct from header file
        self.write_xy_fname = False # controls whether xy_fname is created
                                    # by self.write_input_data, or only header

        # Other possible GeoClaw inputs:
        self.x = None
        self.y = None
        self.X = None
        self.Y = None
        self.Z = None  # for topo DEM values if available
        self.fgmax_point = None  # for point_style==4
        self.force_dry_init = None  # =1 if wet, =0 if dry
        self.dx = None
        self.dy = None

        # possible output values that may be available after run:
        self.outdir = '_output'    # where to find GeoClaw output fgmax*.txt
        self.level = None
        self.X = None
        self.Y = None
        self.Z = None
        self.dx = None
        self.dy = None
        self.B = None
        self.h = None
        self.h_time = None
        self.s = None
        self.s_time = None
        self.hs = None
        self.hs_time = None
        self.hss = None
        self.hss_time = None
        self.hmin = None
        self.hmin_time = None
        self.arrival_time = None

        # possible derived quantities of interest:
        self.dz = None
        self.B0 = None
        self.eta = None
        self.h_onshore = None
        self.label = ''  # text for legend



    def read_fgmax_grids_data(self, fgno, data_file='fgmax_grids.data'):
        """
        Read input info for fgmax grid number fgno from the data file
        fgmax_grids.data, which should have been created by setrun.py.
        This file now contains info about all fgmax grids.
        """
        with open(data_file) as filep:
            lines = filep.readlines()
        #print('+++ opened fgmax_new.data with %i lines' % len(lines))
        fgmax_input = None
        for lineno,line in enumerate(lines):
            if 'fgno' in line:
                if int(line.split()[0]) == fgno:
                    fgmax_input = lines[lineno+1:]
                    #print('Found line %i: %s' % (lineno,line))
                    break

        if fgmax_input is None:
            raise ValueError('fgmax grid fgno = %i not found in %s' \
                             % (fgno, data_file))

        self.fgno = fgno
        self.tstart_max = float(fgmax_input[0].split()[0])
        self.tend_max = float(fgmax_input[1].split()[0])
        self.dt_check = float(fgmax_input[2].split()[0])
        self.min_level_check = int(fgmax_input[3].split()[0])
        self.arrival_tol = float(fgmax_input[4].split()[0])
        self.interp_method = int(fgmax_input[5].split()[0])
        self.point_style = point_style = int(fgmax_input[6].split()[0])
        print('Reading input for fgno=%i, point_style = %i ' \
                % (self.fgno, self.point_style))
        if point_style == 0:
            self.npts = npts = int(fgmax_input[7].split()[0])
            if npts == 0:
                self.xy_fname = fgmax_input[8][1:-2]  # strip quotes
                xy = numpy.loadtxt(self.xy_fname, skiprows=1)
                self.X = xy[:,0]
                self.Y = xy[:,1]
                if xy.shape[1] > 2:
                    self.Z = xy[:,2]  # in case DEM values also stored in input file
                else:
                    self.Z = None
                self.npts = npts = len(self.X)
                print('Read %i x,y points from \n    %s' % (npts, self.xy_fname))
            else:
                X = []; Y = []
                for k in range(8,8+npts):
                    xk = float(fgmax_input[k].split()[0])
                    yk = float(fgmax_input[k].split()[1])
                    X.append(xk)
                    Y.append(yk)
                self.X = numpy.array(X)
                self.Y = numpy.array(Y)
        elif point_style == 1:
            self.npts = npts = int(fgmax_input[7].split()[0])
            self.x1 = float(fgmax_input[8].split()[0])
            self.y1 = float(fgmax_input[8].split()[1])
            self.x2 = float(fgmax_input[9].split()[0])
            self.y2 = float(fgmax_input[9].split()[1])
        elif point_style == 2:
            self.nx = nx = int(fgmax_input[7].split()[0])
            self.ny = ny = int(fgmax_input[7].split()[1])
            self.x1 = float(fgmax_input[8].split()[0])
            self.y1 = float(fgmax_input[8].split()[1])
            self.x2 = float(fgmax_input[9].split()[0])
            self.y2 = float(fgmax_input[9].split()[1])
        elif point_style == 3:
            self.n12 = n12 = int(fgmax_input[7].split()[0])
            self.n23 = n23 = int(fgmax_input[7].split()[1])
            self.x1 = float(fgmax_input[8].split()[0])
            self.y1 = float(fgmax_input[8].split()[1])
            self.x2 = float(fgmax_input[9].split()[0])
            self.y2 = float(fgmax_input[9].split()[1])
            self.x3 = float(fgmax_input[10].split()[0])
            self.y3 = float(fgmax_input[10].split()[1])
            self.x4 = float(fgmax_input[11].split()[0])
            self.y4 = float(fgmax_input[11].split()[1])
        elif point_style == 4:
            self.xy_fname = fgmax_input[7][1:-2]  # strip quotes
            ## Need to read in topotype 3 file and set self.npts
            # xy = numpy.loadtxt(self.xy_fname, skiprows=1)
            # self.X = xy[:,0]
            # self.Y = xy[:,1]
            # if xy.shape[1] > 2:
            #     self.Z = xy[:,2]  # in case DEM values also stored in input file
            # else:
            #     self.Z = None
            # self.npts = npts = len(self.X)
            # print('Read %i x,y points from \n    %s' % (npts, self.xy_fname))


    def write_to_fgmax_data(self, fid):
        """
        Write the fgmax grid data to the file specified by `fid`, normally
        the `fgmax_grids.data` file that is read in by the GeoClaw Fortran code.
        """

        print("\n---------------------------------------------- ")
        point_style = self.point_style
        if point_style not in [0,1,2,3,4]:
            raise NotImplementedError("make_fgmax not implemented for point_style %i" \
                % point_style)

        # write header, independent of point_style:
        #fid = open(self.input_file_name,'w')
        fid.write("\n")
        fid.write("%i                           # fgno\n" % self.fgno)
        fid.write("%16.10e            # tstart_max\n"  % self.tstart_max)
        fid.write("%16.10e            # tend_max\n"  % self.tend_max)
        fid.write("%16.10e            # dt_check\n" % self.dt_check)
        fid.write("%i %s              # min_level_check\n" \
                            % (self.min_level_check,12*" "))

        fid.write("%16.10e            # arrival_tol\n" % self.arrival_tol)
        fid.write("%i %s              # interp_method\n" \
                            % (self.interp_method,12*" "))
        fid.write("%i %s              # point_style\n" \
                            % (self.point_style,12*" "))

        print('fgmax grid %i has point_style = %i' % (self.fgno, point_style))

        if point_style == 0:
            if self.xy_fname is not None:
                fid.write("0         # npts==0 ==> points in this file:\n")
                fid.write("'%s'\n" % self.xy_fname)
                print("points should be in file:")
                print("   %s" % self.xy_fname)
                if self.write_xy_fname:
                    if self.Z is not None:
                        xydata = numpy.vstack([self.X,self.Y,self.Z]).T
                    else:
                        xydata = numpy.vstack([self.X,self.Y])
                    numpy.savetxt(self.xy_fname, xydata,
                                  header='%8i' % len(self.X),
                                  comments='', fmt='%24.14e')
            else:
                #print("+++ expecting xy_fname")
                # list of points
                npts = self.npts

                print("unstructured grid of %s points" % npts)

                fid.write("%i                 # npts\n" % (npts))
                for k in range(npts):
                    fid.write("%22.12f   %22.12f \n" % (self.X[k],self.Y[k]))
            #fid.close()


        elif point_style==1:
            # 1d transect of points
            x1,x2 = self.x1, self.x2
            y1,y2 = self.y1, self.y2
            # require self.npts to be set and don't set dx since
            # ambiguous for general transect in lat-long.
            if self.npts is None:
                raise ValueError('With point_style==1 must set set npts')
            else:
                npts = self.npts
                if self.dx is not None:
                    print("*** Warning: With point_style==1 cannot set dx")


            print("1d fixed grid with %s points" % npts)


            fid.write("%i                 # npts\n" % (npts))
            fid.write("%g   %g            # x1, y1\n" % (x1,y1))
            fid.write("%g   %g            # x2, y2\n" % (x2,y2))
            #fid.close()


            #print("Created file ", self.input_file_name)
            print("   specifying fixed grid with %i points equally spaced from " \
                    % npts)
            print("   (%g,%g)  to  (%g,%g)" % (x1,y1,x2,y2))


        if point_style == 2:
            # 2d grid of points
            x1,x2 = self.x1, self.x2
            y1,y2 = self.y1, self.y2
            if self.nx is None:
                dx = self.dx
                nx = int(round((x2-x1)/dx)) + 1
                if abs((nx-1)*dx + x1 - x2) > 1e-6:
                    print("Warning: abs((nx-1)*dx + x1 - x2) = ", \
                          abs((nx-1)*dx + x1 - x2))
                    print("         old x2: %22.16e" % x2)
                    x2 = x1 + dx*(nx-1)
                    print("         resetting x2 to %22.16e" % x2)
            else:
                nx = self.nx
                dx = (x2-x1)/(nx+1.)
                if self.dx is not None:
                    print("*** Warning: dx specified over-ridden by: ",dx)

            if self.ny is None:
                dy = self.dy
                if dy is None:
                    dy = dx
                ny = int(round((y2-y1)/dy)) + 1
                if abs((ny-1)*dy + y1 - y2) > 1e-6:
                    print("Warning: abs((ny-1)*dy + y1 - y2) = ", \
                          abs((ny-1)*dy + y1 - y2))
                    print("         old y2: %22.16e" % y2)
                    y2 = y1 + dy*(ny-1)
                    print("         resetting y2 to %22.16e" % y2)
            else:
                ny = self.ny
                dy = (y2-y1)/(ny+1.)
                if self.dy is not None:
                    print("*** Warning: dy specified over-ridden by: ",dy)


            npts = nx*ny


            fid.write("%i  %i %s          # nx,ny\n" \
                                % (nx,ny,10*" "))
            fid.write("%16.10e   %20.10e            # x1, y1\n" % (x1,y1))
            fid.write("%16.10e   %20.10e            # x2, y2\n" % (x2,y2))
            #fid.close()


            #print("Created file ", self.input_file_name)
            print("   specifying fixed grid with shape %i by %i, with  %i points" \
                    % (nx,ny,npts))
            print("   lower left  = (%15.10f,%15.10f)" % (x1,y1))
            print("   upper right = (%15.10f,%15.10f)" % (x2,y2))
            print("   dx = %15.10e,  dy = %15.10e" % (dx,dy))


        elif point_style==3:
            # arbitrary quadrilateral
            x1,x2 = self.x1, self.x2
            y1,y2 = self.y1, self.y2
            x3,x4 = self.x3, self.x4
            y3,y4 = self.y3, self.y4
            if self.n12 is None:
                raise NotImplementedError("Need to set n12 and n23")
            else:
                npts = self.n12 * self.n23

            fid.write("%i  %i %s          # self.n12,self.n23\n" \
                                % (self.n12,self.n23,10*" "))
            fid.write("%16.10e   %20.10e            # x1, y1\n" % (x1,y1))
            fid.write("%16.10e   %20.10e            # x2, y2\n" % (x2,y2))
            fid.write("%16.10e   %20.10e            # x3, y3\n" % (x3,y3))
            fid.write("%16.10e   %20.10e            # x4, y4\n" % (x4,y4))
            #fid.close()


            #print("Created file ", self.input_file_name)
            print("   specifying fixed grid as a quadrilateral")
            print("       %i by %i, with  %i points" \
                    % (self.n12,self.n23,npts))
            print("   corner 1 = (%15.10f,%15.10f)" % (x1,y1))
            print("   corner 2 = (%15.10f,%15.10f)" % (x2,y2))
            print("   corner 3 = (%15.10f,%15.10f)" % (x3,y3))
            print("   corner 4 = (%15.10f,%15.10f)" % (x4,y4))


        elif point_style == 4:
            if self.xy_fname is not None:
                fid.write("'%s'\n" % self.xy_fname)
                print("points should be in file:")
                print("   %s" % self.xy_fname)
            else:
                raise ValueError('for point_style==4, require xy_fname')

    def read_output(self, fgno=None, outdir=None, verbose=True, 
                    indexing='ij'):
        r"""
        Read the GeoClaw results on the fgmax grid numbered *fgno*.
        
        indexing='ij' gives backward compatibility.
           X[i,j],Y[i,j] corresponds to point x[i],y[j]

        Alternatively, can set indexing=='xy' so that X,Y and other
        arrays have same layout as topo arrays:
           X[j,i],Y[j,i] corresponds to point x[i],y[j]
        This is useful if you want to save the fgmax results in same format as 
        topofiles, using topotools.Topography.write().
           
        """

        if indexing == 'xy':
            reshape_order = 'C'
        elif indexing == 'ij':
            reshape_order = 'F'
        else:
            raise InputError("*** indexing must by 'xy' or 'ij'")
            
        if self.point_style is None:
            raise InputError("*** point_style is not set, need to read input?")
        point_style = self.point_style

        if fgno is not None:
            self.fgno = fgno
        if outdir is not None:
            self.outdir = outdir

        # Require new style in v5.7.0, e.g. fgmax0001.txt etc.
        fname = os.path.join(self.outdir, 'fgmax%s.txt' \
                % str(self.fgno).zfill(4))

        if not os.path.isfile(fname):
            raise IOError("File not found: %s" % fname)

        print("Reading %s ..." % fname)
        d = numpy.loadtxt(fname)

        if point_style == 4:
            self.npts = d.shape[0]
            print('point_style == 4, found %i points ' % self.npts)


        # new format in v5.7.0, includes column for B = topo from aux array
        cols_expected = [7,9,15]

        ncols = d.shape[1]

        if ncols not in cols_expected:
            raise IOError("*** Unexpected number of columns %s in file %s" \
                    % (ncols, fname))

        ind_s = None
        ind_hs = None
        ind_hss = None
        ind_hmin = None
        ind_h_time = None
        ind_s_time = None
        ind_hs_time = None
        ind_hss_time = None
        ind_hmin_time = None

        ind_x = 0
        ind_y = 1
        ind_level = 2
        ind_B = 3  # added in new fname style
        ind_h = 4
        if ncols == 7:
            ind_h_time = 5
            ind_arrival_time = 6
        elif ncols == 9:
            ind_s = 5
            ind_h_time = 6
            ind_s_time = 7
            ind_arrival_time = 8
        elif ncols == 15:
            ind_s = 5
            ind_hs = 6
            ind_hss = 7
            ind_hmin = 8
            ind_h_time = 9
            ind_s_time = 10
            ind_hs_time = 11
            ind_hss_time = 12
            ind_hmin_time = 13
            ind_arrival_time = 14


        if point_style in [0,1,4]:
            fg_shape = (self.npts,)
        elif point_style == 2:
            if indexing == 'xy':
                fg_shape = (self.ny,self.nx)
            else:
                fg_shape = (self.nx,self.ny)

        elif point_style == 3:
            if indexing == 'xy':
                fg_shape = (self.n23,self.n12)
            else:
                fg_shape = (self.n12,self.n23)
        else:
            raise NotImplementedError("Not implemented for point_style %s" \
                % point_style)

        X = numpy.reshape(d[:,0],fg_shape,order=reshape_order)
        Y = numpy.reshape(d[:,1],fg_shape,order=reshape_order)
        y0 = 0.5*(Y.min() + Y.max())   # mid-latitude for scaling plots
        h = numpy.reshape(d[:,ind_h],fg_shape,order=reshape_order)

        # AMR level used for each fgmax value:
        level = numpy.reshape(d[:,ind_level].astype('int'),fg_shape,
                              order=reshape_order)

        # Set B = topo array
        B = numpy.reshape(d[:,ind_B],fg_shape,order=reshape_order)

        mask = (h < -1e50)  # points that were never set
        B = ma.masked_where(mask, B)
        h = ma.masked_where(mask, h)

        def set_q_time(ind_q, ind_q_time):
            q = numpy.reshape(d[:,ind_q],fg_shape,order=reshape_order)
            q = ma.masked_where(mask,q)
            q_time = numpy.reshape(d[:,ind_q_time],fg_shape,order=reshape_order)
            q_time = ma.masked_where(mask, q_time)
            return q, q_time

        self.h, self.h_time = set_q_time(ind_h, ind_h_time)
        if ind_s:
            self.s, self.s_time = set_q_time(ind_s, ind_s_time)
        if ind_hs:
            self.hs, self.hs_time = set_q_time(ind_hs, ind_hs_time)
            self.hss, self.hss_time = set_q_time(ind_hss, ind_hss_time)
            self.hmin, self.hmin_time = set_q_time(ind_hmin, ind_hmin_time)

        # last column is arrival times:
        arrival_time = numpy.reshape(d[:,ind_arrival_time],
                                     fg_shape,order=reshape_order)
        arrival_time = ma.masked_where(arrival_time < -1e50, arrival_time)
        arrival_time = ma.masked_where(mask, arrival_time)
        self.arrival_time = arrival_time

        self.level = level
        self.X = X
        self.Y = Y
        self.B = B
        self.h = h

        # do not set these, leave for user to do as desired:
        if 0:
            self.B0 = B  ## SHOULD MODIFY BY dz!

            if self.force_dry_init is not None:
                self.h_onshore = ma.masked_where(self.force_dry_init==0, self.h)
            else:
                self.h_onshore = ma.masked_where(self.B0 < 0., self.h)

        if point_style==4:
            #print('Returning lists, convert to masked arrays based on input grid')
            #print('   using to_arrays() function')
            try:
                self.ps4_to_arrays(verbose=verbose)
            except:
                print('*** Problem converting from 1d lists to 2d arrays,\n' \
                      + '    Trying to map onto grid specified by:\n    ', \
                      self.xy_fname)
                raise

        if self.X.ndim==2:
            if indexing=='xy':
                self.x = self.X[0,:]
                self.y = self.Y[:,0]
            else:
                # for indexing=='ij' this fixes bug in v5.9.0 version:
                self.x = self.X[:,0]
                self.y = self.Y[0,:]
        else:
                self.x = self.X
                self.y = self.Y


    def bounding_box(self):
        """
        Return the bounding box of the grid as a list [x1,x2,y1,y2]
        """
        x1 = self.X.min()
        x2 = self.X.max()
        y1 = self.Y.min()
        y2 = self.Y.max()
        return [x1,x2,y1,y2]

    def ps4_to_arrays(self, verbose=True):
        """
        for point_style==4, convert lists of fgmax values into masked arrays
        based on the topo_style==3 file self.xy_fname that was used to specify
        the fgmax points in the GeoClaw run.
        """

        from numpy import ma
        assert self.point_style==4, '*** Requires point_style==4'

        if self.X.ndim==2 or self.Y.ndim==2:
            print('*** X and Y already 2d, not converting')
            return

        x_1d = self.X
        y_1d = self.Y

        if verbose:
            print('Will map fgmax points onto masked arrays defined by file:')
            print('     %s' % self.xy_fname)

        from clawpack.geoclaw import topotools
        pts_chosen = topotools.Topography(path=self.xy_fname, topo_type=3)
        X = pts_chosen.X
        Y = pts_chosen.Y
        mask = numpy.logical_not(pts_chosen.Z)
        x1 = X.min()
        y1 = Y.min()

        # possible arrays from GeoClaw output to convert:
        zarrays = ['level','B','h','h_time','s','s_time','hs','hs_time',\
                   'hss','hss_time','hmin','hmin_time','arrival_time']

        dx = X[0,1] - X[0,0]
        dy = Y[1,0] - Y[0,0]
        if verbose:
            print('Deduced dx = %g, dy = %g'  % (dx,dy))

        for attr in zarrays:
            z_1d = getattr(self, attr, None)
            if z_1d is None:
                if verbose: print('not converting attribute %s == None' % attr)
            else:
                Z = ma.masked_array(data=numpy.empty(X.shape), mask=True)
                for k in range(len(x_1d)):
                    i = int(round((x_1d[k]-x1)/dx))
                    j = int(round((y_1d[k]-y1)/dy))
                    Z[j,i] = z_1d[k]
                if 0:
                    if not numpy.alltrue(mask == Z.mask):
                        print('*** converting to arrays gave unexpected mask for')
                        print('    Z array =  %s' % attr)
                setattr(self, attr, Z)
                if verbose: print('converted %s to 2d array' % attr)

        self.X = X
        self.Y = Y



    def interp_dz(self,dtopo_path,dtopo_type):
        """
        Compute approximate values of deformation dz on X,Y grid using
        a specified dtopo file.
        Also calculates B0 = B - dz, attempting to recover the pre-event
        topography from the GeoClaw run topography stored in B.
        """
        from clawpack.geoclaw import dtopotools
        from scipy.interpolate import RegularGridInterpolator

        dtopo = dtopotools.DTopography(dtopo_path, dtopo_type=dtopo_type)
        x1d = dtopo.X[0,:]
        y1d = dtopo.Y[:,0]
        dtopo_func = RegularGridInterpolator((x1d,y1d), dtopo.dZ[-1,:,:].T,
                        method='linear', bounds_error=False, fill_value=0.)
        dz = dtopo_func(list(zip(numpy.ravel(self.X), numpy.ravel(self.Y))))
        self.dz = numpy.reshape(dz, self.X.shape)
        print('Over fgmax extent, min(dz) = %.2f m, max(dz) = %.2f m' \
             % (dz.min(), dz.max()))

        # leave this for user to do if desired:
        # self.B0 = self.B - self.dz


def adjust_fgmax_1d(x1_desired, x2_desired, x1_domain, dx):
    """
    Adjust the upper and lower limits of a grid so that equally spaced
    grid points with spacing `dx` lie exactly at cell centers, so that
    no interpolation is needed for fgmax values.  Note that parameter
    names refer to `x` limits, but works equally well for `y` values.

    :Input:
     - x1_desired, x2_desired: approximate desired limits of fgmax grid
     - x1_domain:  lower edge of computational domain
     - dx: Mesh spacing on fine grid that fgmax grid should conform to
    :Output:
     - x1_new, x2_new: limits to set so (x2-x1) is integer multiple
       of dx and points are at cell centers of computational grid
     - npoints: number of points to specify, so that
       `linspace(x1_new, x2_new, npoints) gives points with spacing `dx`.
    """

    i1 = numpy.floor((x1_desired-x1_domain - 0.5*dx)/dx)
    x1_new = x1_domain + (i1 + 0.5)*dx
    i2 = numpy.floor((x2_desired-x1_domain + 0.5*dx)/dx)
    x2_new = x1_domain + (i2 + 0.5)*dx
    npoints = int(i2 - i1) + 1
    return x1_new, x2_new, npoints


def adjust_fgmax_grid(x1_desired, x2_desired, x1_domain, dx,
                      y1_desired, y2_desired, y1_domain, dy=None, verbose=True):


    if dy == None:
        dy = dx

    x1_new, x2_new, nx = adjust_fgmax_1d(x1_desired, x2_desired, x1_domain, dx)
    y1_new, y2_new, ny = adjust_fgmax_1d(y1_desired, y2_desired, y1_domain, dy)

    if verbose:
        print("x:")
        print("  moved %17.12f to %17.12f by %g" % (x1_desired, x1_new, abs(x1_desired-x1_new)))
        print("  moved %17.12f to %17.12f by %g" % (x2_desired, x2_new, abs(x2_desired-x2_new)))
        print("y:")
        print("  moved %17.12f to %17.12f by %g" % (y1_desired, y1_new, abs(y1_desired-y1_new)))
        print("  moved %17.12f to %17.12f by %g" % (y2_desired, y2_new, abs(y2_desired-y2_new)))
        #print "  "
        #print "fg.nx = %g" % nx
        #print "fg.ny = %g" % ny
        #print "fg.x1 = %17.12f" % x1_new
        #print "fg.x2 = %17.12f" % x2_new
        #print "fg.y1 = %17.12f" % y1_new
        #print "fg.y2 = %17.12f" % y2_new
    return x1_new, x2_new, nx, y1_new, y2_new, ny

