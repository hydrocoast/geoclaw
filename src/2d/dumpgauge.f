      subroutine dumpgauge(q,aux,xlow,ylow,nvar,mitot,mjtot,naux,mptr)

      use amr_module
      use geoclaw_module
      use gauges_module

      implicit double precision (a-h,o-z)

      integer bsearch
      dimension q(nvar,mitot,mjtot), var(maxvar)
      dimension aux(naux,mitot,mjtot)
      dimension eta(num_layers)
      dimension h(num_layers,4)

c  # see if this grid contains any gauges so data can be output
c  # may turn out this should be sorted, but for now do linear search
c
c  # array is sorted according to indices in mbestorder array
c  # so do binary search to find start. Could have many same source grids
c
      if (mgauges.eq.0) then
         return
         endif


      istart = bsearch(mptr) 
      if (istart .lt. 1) return      !this grid not used

c     # this stuff the same for all gauges on this grid
      tgrid = rnode(timemult,mptr)
      level = node(nestlevel,mptr)
      hx    =  hxposs(level)
      hy    =  hyposs(level)

      do 10 ii = istart, mgauges
        i = mbestorder(ii)   ! gauge number
        if (mptr .ne. mbestsrc(i)) go to 99  ! all done
        if (tgrid.lt.t1gauge(i) .or. tgrid.gt.t2gauge(i)) then
c          # don't output at this time for gauge i
           return
           endif
c
c
c ## prepare to do linear interp at gauge location to get vars
c ## should fancier (limited) interp be done??
c
        iindex =  int(.5 + (xgauge(i)-xlow)/hx)
        jindex =  int(.5 + (ygauge(i)-ylow)/hy)
        if ((iindex .lt. nghost .or. iindex .gt. mitot-nghost) .or.
     .      (jindex .lt. nghost .or. jindex .gt. mjtot-nghost))
     .    write(*,*)"ERROR in output of Gauge Data "
        xcent  = xlow + (iindex-.5)*hx
        ycent  = ylow + (jindex-.5)*hy
        xoff   = (xgauge(i)-xcent)/hx
        yoff   = (ygauge(i)-ycent)/hy
        if (xoff < 0. .or. xoff > 1 .or. yoff < 0. .or. yoff > 1.) then
            print *, " BIG PROBLEM in DUMPGAUGE", i
        endif

c ## Modified by RJL 12/31/09 to interpolate only where all four cells are
c ## wet, otherwise just take this cell value:

c Check for dry cells by comparing h to drytol2, which should be smaller
c than drytolerance to avoid oscillations since when h < drytolerance the
c velocities are zeroed out which can then lead to increase in h again.

        drytol2 = 0.1d0 * dry_tolerance(1)

          do m=1,num_layers
              layer_index = 3*(m-1)
              h(m,1) = q(layer_index+1,iindex,jindex) / rho(m)
              h(m,2) = q(layer_index+1,iindex+1,jindex) / rho(m)
              h(m,3) = q(layer_index+1,iindex,jindex+1) / rho(m)
              h(m,4) = q(layer_index+1,iindex+1,jindex+1) / rho(m)
              
              if ((h(m,1) < drytol2) .or.
     &            (h(m,2) < drytol2) .or.
     &            (h(m,3) < drytol2) .or.
     &            (h(m,4) < drytol2)) then
                  ! One of the cells is dry, so just use value from grid cell
                  ! that contains gauge rather than interpolating
                  
                  icell = int(1.d0 + (xgauge(i) - xlow) / hx)
                  jcell = int(1.d0 + (ygauge(i) - ylow) / hy)
                  do ivar=1,3
                      var(ivar + layer_index) = 
     &                       q(ivar + layer_index,icell,jcell) / rho(m)
                  enddo
                  if (m == num_layers) then
                      ! This is the bottom layer and we should figure out the
                      ! topography
                      topo = aux(1,icell,jcell)
                  endif
              else
                  ! Linear interpolation between four cells
                  do ivar=1,3
                      var(layer_index + ivar) = (1.d0 - xoff) * 
     &                   (1.d0 - yoff)
     &                 * q(layer_index + ivar,iindex,jindex) / rho(m)
     &                 + xoff*(1.d0 - yoff) 
     &                 * q(layer_index + ivar,iindex+1,jindex) / rho(m)
     &                 + (1.d0 - xoff) * yoff 
     &                 * q(layer_index + ivar,iindex,jindex+1) / rho(m)
     &                 + xoff * yoff 
     &                 * q(layer_index + ivar,iindex+1,jindex+1)/rho(m)
                  enddo
                  if (m == num_layers) then
                      topo = (1.d0 - xoff) * (1.d0 - yoff) 
     &                        * aux(1,iindex,jindex) 
     &                      + xoff * (1.d0 - yoff) 
     &                        * aux(1,iindex+1,jindex) 
     &                      + (1.d0 - xoff) * yoff 
     &                        * aux(1,iindex,jindex+1) 
     &                      + xoff * yoff 
     &                        * aux(1,iindex+1,jindex+1)
                  endif
              endif
          enddo
              
          ! Extract surfaces
          eta(num_layers) = var(3*num_layers-2) + topo
          do k=num_layers-1,1,-1
              eta(k) = var(3*k-2) + eta(k+1)
          enddo
              
          write(OUTGAUGEUNIT,100) igauge(i),level,tgrid, 
     &              (var(j),j=1,3*num_layers),(eta(j),j=1,num_layers)
  10  enddo
      
 100  format(2i5,15e15.7)
 
  99  return
 
      end subroutine dumpgauge
c
c --------------------------------------------------------------------
c
      subroutine setbestsrc()
c
c  ## called every time grids change, to set the best source grid
c  ## to find gauge data
c
c  ## lbase is grid level that didn't change but since fine
c  ## grid may have disappeared, still have to look starting
c  ## at coarsest level 1.
c
      use amr_module
      use gauges_module
      implicit double precision (a-h,o-z)

c
c ##  set source grid for each loc from coarsest level to finest.
c ##  that way finest src grid left and old ones overwritten
c ##  this code uses fact that grids do not overlap

c # for debugging, initialize sources to 0 then check that all set
      do i = 1, mgauges
         mbestsrc(i) = 0
      end do

 
      do 20 lev = 1, lfine  
          mptr = lstart(lev)
 5        do 10 i = 1, mgauges
            if ((xgauge(i) .ge. rnode(cornxlo,mptr)) .and.    
     .          (xgauge(i) .le. rnode(cornxhi,mptr)) .and.    
     .          (ygauge(i) .ge. rnode(cornylo,mptr)) .and.  
     .          (ygauge(i) .le. rnode(cornyhi,mptr)) )
     .      mbestsrc(i) = mptr
 10       continue

          mptr = node(levelptr, mptr)
          if (mptr .ne. 0) go to 5
 20   continue


      do i = 1, mgauges
        if (mbestsrc(i) .eq. 0) 
     .      write(6,*)"ERROR in setting grid src for gauge data",i
      end do

c
c     sort the source arrays for easy testing during integration
      call qsorti(mbestorder,mgauges,mbestsrc)

      return
      end
c
c ------------------------------------------------------------------------
c
      integer function bsearch(mptr)

      use gauges_module
      implicit double precision (a-h,o-z)

      bsearch = -1           ! signal if not found

      indexlo = 1
      indexhi = mgauges

 5    if (indexhi .lt. indexlo) go to 99
      mid = (indexlo + indexhi)/2

      if (mptr .gt. mbestsrc(mbestorder(mid))) then
          indexlo = mid+1
          go to 5
      else if (mptr .lt. mbestsrc(mbestorder(mid))) then
          indexhi = mid-1
          go to 5
      else    ! found the grid. find its first use in the array
        istart = mid


 10      if (istart .gt. 1) then
            if (mbestsrc(mbestorder(istart-1)) .ne. mptr) go to 90
            istart = istart - 1
            go to 10
          endif

      endif

 90   bsearch = istart

 99   return
      end

