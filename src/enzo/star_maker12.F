#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker12(nx, ny, nz,
     &                      d, dm, temp, u, v, w, cooltime,
     &                      dt, r, metal, dx, t, z, procnum, 
     &                      d1, x1, v1, t1,
     &                      nmax, xstart, ystart, zstart, ibuff, 
     &                      imetal, imethod, mintdyn,
     &                      odthresh, masseff, smthresh, level, 
     &                      np,
c     &                      rcut, gravcenter, np,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, tdp, tcp, metalf,
     &                      Bx, By, Bz)
c
c  CREATES STAR PARTICLES
c
c  Adapted by Elizabeth Tasker, July 2008
c  New feedback added, September 2008
c
c  modified1: enzo-dev (2.4) implementation BW, Feb 2016
c  modified1.1: enzo-dev (2.4) updated BW, Oct 2016 (fixed B units)
c 
c  Based on star_maker2.src (though with similar stochastic SF to star_maker1), 
c  this routine restricts star formation criteria to the density threshold; 
c  idea from Krumholz & Tan (2007) who argue a fixed SFR/t_ff in turbulant gas. 
c  Since gas must also be molecular, T < 3000 K. 
c
c  star_maker12:
c
c  If the gas is magnetically supercritical, forms stars as star_maker2
c
c  INPUTS:
c
c    d     - density field
c    dm    - dark matter field
c    temp  - temperature field
c    u,v,w - velocity fields
c    cooltime - cooling time in code units
c    r     - refinement field (non-zero if zone is further refined)
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
c    odthresh - overdensity threshold (some number * avg. density)
c               if negative, then use constant time rather than tdyn
c    masseff - gas-to-mass conversion efficiency ( 0<=masseff<=1 )
c    smthresh - star mass threshold (only creates stars with mass >
c        smthresh unless (random number) < starmass/smthresh )
c    mintdyn  - minimum dynamical time, in years
c    level - current level of refinement
c    procnum - processor number (for output)
c    gravcenter - Point source gravity position
c    rcut - causes radial cutoff for star particles (avoids central 3.5 kpc of disk). 
c    Bx, By, Bz - the magnetic field in code units.
c
c  OUTPUTS:
c
c    np   - number of particles created
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created
c    tcp      - creation time of particle
c    metalf   - metallicity fraction of particle
c    nmax     - particle array size specified by calling routine
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, nmax, np, level, imetal, imethod
      INTG_PREC procnum
c      INTG_PREC procnum, rcut
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), temp(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    r(nx,ny,nz), cooltime(nx,ny,nz), metal(nx,ny,nz)
c      R_PREC    gravcenter(3)
      R_PREC    dt, dx, z
      R_PREC    d1, x1, v1, t1
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(nmax), yp(nmax), zp(nmax)
      R_PREC    up(nmax), vp(nmax), wp(nmax)
      R_PREC    mp(nmax), tdp(nmax), tcp(nmax), metalf(nmax)
      R_PREC    odthresh, masseff, smthresh, mintdyn
      R_PREC    Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)

c
c      R_PREC   sformsum
c      save   sformsum
c      data   sformsum/0/
c
c  Locals:
c
      INTG_PREC  i, j, k, ii, clock
      R_PREC   div, tdyn, dtot
      R_PREC   pi, G, sndspdC
      R_PREC   isosndsp2, starmass, starfraction, bmass, jeanmass
      R_PREC   Bmag, Bregconst, magunits
      real*8 msolar
      parameter (pi=pi_val, G=GravConst, 
     &           sndspdC=1.3095e8_RKIND,
     &           msolar=SolarMass)
      character (len=90) :: filename

c
c  Stuff for the RNG:
c

      logical, save :: firstcall = .true.
      integer :: seed_size
      integer, dimension(:), allocatable :: seed

      R_PREC random, x, starradius
      external rand
c
      ii = 0
c
c  Initialize the RNG on the first call to the routine.
c

#define STOCHASTIC_STAR_FORMATION
c
#ifdef STOCHASTIC_STAR_FORMATION

      if (firstcall) then

         ! Get the size of the seed
         call random_seed(SIZE=seed_size)
         allocate(seed(seed_size))


         call system_clock(count=clock)
         seed = clock + 142 * [(i, i = 0,seed_size-1)]
         call random_seed(PUT=seed)

         deallocate(seed)
         firstcall = .false.
      endif
#endif

c
c  for each zone, : "star" particle is created if answers to all the
c  following questions are affirmative:
c
c    is this the finest level of refinement ?
c    is the density greater than a critical density ?
c    is the gas temperature lower than 3e3 (molecular)?
c

c      write(6,*) 'star_maker5 start'
c      write(*,*) 'star_maker5 start'
c      write(*,*) 'rcut = ', rcut
      
c      if (rcut > 0) write(6,*) 
c     & 'WARNING: radial cut-off in place. No stars at r < 2.5kpc'

c      write(*,*) 'star_maker11 start'
      write (filename, '( "SFcelldata_", I3.3, ".txt" )' ) procnum
      open (unit = 9000+procnum, file = filename, position="append")

      magunits = dsqrt(4._RKIND*pi*d1*v1**2)

      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff

c               write(*,*) d, dm, temp, u, v, w, cooltime //
c
c              1) finest level of refinement?
c
               if (r(i,j,k) .ne. 0._RKIND) goto 10
c
c              2) density greater than threshold
c
c               if (d(i,j,k) .lt. abs(odthresh)) goto 10
                
c              2b) density greater than critical mass-to-flux

c              Magnitude of B in physical units:
               Bmag = dsqrt(Bx(i,j,k)**2 + By(i,j,k)**2
     &               + Bz(i,j,k)**2)*magunits
               Bregconst = 0.126_RKIND/(dx*x1*sqrt(G))

               if (d(i,j,k)*d1 .lt. Bmag*Bregconst) goto 10
c               if (d(i,j,k)*d1 .lt. 0.5_RKIND*Bmag*Bregconst) goto 10
c               if (d(i,j,k)*d1 .lt. 2.0_RKIND*Bmag*Bregconst) goto 10


c              3) T < 3000 K

               if (temp(i,j,k) .gt. 3.e3_RKIND) goto 10 

c
c              4) Check to see if star is above threshold (given
c                 in units of solarmass)
c
               dtot = ( d(i,j,k) + dm(i,j,k) )*d1
               tdyn  = sqrt(3._RKIND*pi/32._RKIND/G/dtot)/t1

               if (odthresh .lt. 0._RKIND) 
     &               tdyn = sqrt(3._RKIND*pi/32._RKIND/G/abs(odthresh))
     &                      /t1
c               starfraction = min(masseff*dt/tdyn, 0.9_RKIND)
               starfraction = min(masseff*dt/tdyn, 0.5_RKIND)
               tdyn = max(tdyn, mintdyn*3.15e7_RKIND/t1)


#ifdef STOCHASTIC_STAR_FORMATION

c        if mass of star is less than mass threshold, only creates star
c        if (some random #) < starmass/massthresh 
c        and then create it with massthresh, unless not enough gas in cell.

      
               bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3 / msolar

               if (starfraction*bmass .lt. smthresh) then
c                  write(6,*) 'star_maker5: stochastic sf'
c                  set gas in cell thresh to 0.5
                  if (smthresh/bmass .gt. 0.5) goto 10
                  call random_number(random)
                  x = random
                  if (x .gt. (starfraction*bmass / smthresh) ) goto 10 
c                  starfraction = min(smthresh/bmass, 0.9)
c                  if (smthresh/bmass .gt. 0.9) goto 10
                  starfraction = smthresh/bmass

               endif

               write(*,*) 'I have beaten the probability ', 
     &                    starfraction*bmass

#else
c
c                 If star mass is less than threshold, then do make it
c
               if (starfraction*bmass .lt. smthresh) goto 10
#endif
c
c              Create a star particle
c
 
c               write(40,*) '# star_maker5 make star'
c               write(40,*) starfraction, bmass,smthresh,dt,tdyn   


               ii = ii + 1
               mp(ii)  = starfraction * d(i,j,k)
               tcp(ii) = t
               tdp(ii) = tdyn
               xp(ii) = xstart + (REAL(i,RKIND)-0.5_RKIND)*dx
               yp(ii) = ystart + (REAL(j,RKIND)-0.5_RKIND)*dx
               zp(ii) = zstart + (REAL(k,RKIND)-0.5_RKIND)*dx
               if (imethod .eq. 2) then
                  up(ii) = 0.5_RKIND*(u(i,j,k)+u(i+1,j,k))
                  vp(ii) = 0.5_RKIND*(v(i,j,k)+v(i,j+1,k))
                  wp(ii) = 0.5_RKIND*(w(i,j,k)+w(i,j,k+1))
               else
                  up(ii) = u(i,j,k)
                  vp(ii) = v(i,j,k)
                  wp(ii) = w(i,j,k)
               endif
c
c              Set the particle metal fraction
c
               if (imetal .eq. 1) then
                  metalf(ii) = metal(i,j,k)    ! in here metal is a fraction
c                  metalf(ii) = metal(i,j,k)/d(i,j,k)
               else
                  metalf(ii) = 0._RKIND
               endif

c         write out: star index, star mass, sim time, cell dyn time,
c                    cell density, temp, x, y, z, vx, vy, vz, Bx, By, Bz
c               write(200,2001) ii, starfraction*bmass, tcp(ii)*t1, 
               write(9000+procnum,2001) ii, starfraction*bmass, 
     &                    tcp(ii)*t1, tdp(ii)*t1, 
     &                    d(i,j,k)*d1, temp(i,j,k),
     &                    xp(ii)*x1, yp(ii)*x1, zp(ii)*x1, 
     &                    up(ii)*v1, vp(ii)*v1, wp(ii)*v1, 
     &                    Bx(i,j,k)*magunits, By(i,j,k)*magunits,
     &                    Bz(i,j,k)*magunits
 2001          format(i5,14(1pe10.3,1x))


c
c              Remove mass from grid
c
               d(i,j,k) = (1._RKIND - starfraction)*d(i,j,k)
c
c               write(7+procnum,1010) i,j,k, tdp(ii),tcp(ii),
c     &                       metalf(ii),d(i,j,k),dm(i,j,k),temp(i,j,k)
c               write(7+procnum,1000) level,bmass*starfraction,tcp(ii),
c     &                           tdp(ii)*t1,d(i,j,k)*d1,z,metalf(ii)

 1000          format(i5,1x,6(1pe10.3,1x))
 1010          format(3i5,10(1pe10.3,1x))
c
c              Do not generate more star particles than available
c
               if (ii .eq. nmax) goto 20
c
10          continue
c
            enddo
         enddo
      enddo
 20   continue
c	
      if (ii .ge. nmax) then
         write(6,*) 'star_maker11: reached max new particle count',nmax
         stop
      endif
      np = ii


      return
      end


cc=======================================================================
cc/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
cc
c      subroutine star_feedback11(nx, ny, nz,
c     &                      d, te, ge, u, v, w,
c     &                      idual, imethod, dt, r, dx, t,
c     &                      d1, x1, v1, t1, 
c     &                      npart, xstart, ystart, zstart,
c     &                      xp, yp, zp, up, vp, wp,
c     &                      mp, tcp,
c     &                      fbenergy, fbmassloss, fbstart, fbend, 
c     &                      nfb)
cc
cc  RELEASES "STAR" PARTICLE ENERGY
cc
cc  written by: Elizabeth Tasker
cc  date:       September, 2008
cc  modified1:
cc
cc  INPUTS:
cc
cc    d     - density field
cc    te,ge - total energy and gas energy fields
cc    u,v,w - velocity fields
cc    r     - refinement field (0 if zone is further refined)
cc    dt    - current timestep
cc    dx    - zone size (code units)
cc    t     - current time
cc    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
cc    nx,ny,nz - dimensions of field arrays
cc    idual    - dual energy flag
cc    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
cc
cc    x/y/z start - starting position of grid origin
cc    xp,yp,zp - positions of created particles
cc    up,vp,wp - velocities of created particles
cc    mp       - mass of new particles
cc    tcp      - creation time of particle (-1 if not a star particle)
cc    npart    - particle array size specified by calling routine
cc
cc    fbenergy - Energy/Mass/time from feedback.in (SNe + winds) as function of time
cc    fbmassloss - Mass loss / Mass / time "
cc    fbstart  - starting time for fbenergy
cc    fbend    - end time for fbenergy
cc    nfb      - number of bins
cc
cc  OUTPUTS:
cc    ge,e - modified field
cc
cc-----------------------------------------------------------------------
c       implicit none
cc-----------------------------------------------------------------------
cc
cc
cc  Arguments
cc
c      INTG_PREC nx, ny, nz, npart, idual, imethod, nfb
c      R_PREC    d(nx,ny,nz), te(nx,ny,nz)
c      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
c      R_PREC    r(nx,ny,nz), ge(nx,ny,nz)
c      R_PREC    dt, dx, z
c      R_PREC    d1, x1, v1, t1
c      P_PREC xstart, ystart, zstart, t
c      P_PREC xp(npart), yp(npart), zp(npart)
c      R_PREC    up(npart), vp(npart), wp(npart)
c      R_PREC    mp(npart), tcp(npart)
c      R_PREC    fbenergy, fbmassloss, fbstart, fbend
cc
cc  Locals
cc
c      INTG_PREC i, j, k, n, index
c      R_PREC    log_energy, energy, remainder, sptime
c      R_PREC    log_m_eject, m_eject, dratio
c    
cc
cc-----------------------------------------------------------------------
cc
cc     Loop over star particles
cc
cc      write(6,*) 'star_feedback5: start'
c      do n=1, npart
c         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0) then  
cc
cc
cc     interpolate from fbenergy the energy injects by star 
cc     particle at its current age
cc
cc
cc     index number of particle time in feedback array
cc       
c            sptime = max(tcp(n), fbstart)
c            sptime = min(tcp(n), fbend) 
cc
cc
c            index = min(nfb-1, 
c     &           max(1,int((sptime - fbstart)/(fbend- fbstart)*
c     &           real(nfb-1))+1))
cc     
c            remainder = ((sptime - fbstart)/(fbend- fbstart)*
c     &           real(nfb-1)+1)- real(index)
cc      
cc
c            log_energy = log(fbenergy(index)) + remainder*
c     &           (log(fbenergy(index+1)) - log(fbenergy(index)))
c
c            log_m_eject = log(fbmassloss(index)) + remainder*
c     &           (log(fbmassloss(index+1)) - log(fbmassloss(index)))
c
cc
cc                
cc
cc     Compute index of nearest cell
cc     
c            i = int((xp(n) - xstart)/dx) + 1
c            j = int((yp(n) - ystart)/dx) + 1
c            k = int((zp(n) - zstart)/dx) + 1
cc
c            if (i .lt. 1 .or. i .gt. nx .or. j .lt. 1 .or. j .gt. ny
c     &           .or. k .lt. 1 .or. k .gt. nz) then
c               write(6,*) 'warning: star particle out of grid',i,j,k
c               goto 100
c            endif
cc
cc 
c            energy = 10**log_energy * mp(n) * dt/ (d(i,j,k)+mp(n))
c
c            m_eject = 10**log_m_eject * mp(n) * dt
c
c            mp(n) = mp(n) - m_eject
cc
cc           
cc           Add energy to energy field
cc     
c
cc            write(6,*) 'te before: ', te(i,j,k)
c
c            dratio = d(i,j,k)/(d(i,j,k) + m_eject)
c            te(i,j,k) = te(i,j,k)*dratio + energy
c            if (idual .eq. 1) 
c     &           ge(i,j,k) = ge(i,j,k)*dratio + energy
cc
cc
cc            write(6,*) 'te after: ', te(i,j,k)
c
c
cc
cc              Mass and momentum feedback / conservation
cc
c            u(i,j,k) = u(i,j,k)*d(i,j,k) + m_eject * up(n)
c            v(i,j,k) = v(i,j,k)*d(i,j,k) + m_eject * vp(n)
c            w(i,j,k) = w(i,j,k)*d(i,j,k) + m_eject * wp(n)
c            d(i,j,k) = d(i,j,k) + m_eject
c            u(i,j,k) = u(i,j,k)/d(i,j,k)
c            v(i,j,k) = v(i,j,k)/d(i,j,k)
c            w(i,j,k) = w(i,j,k)/d(i,j,k)
cc
cc          If te is really total energy (and it is unless imethod=2),
cc             then just set this value
cc
c            if (imethod .ne. 2 .and. idual .eq. 1) then
c               te(i,j,k) = 0.5*(u(i,j,k)**2 + v(i,j,k)**2 + 
c     &                          w(i,j,k)**2) + ge(i,j,k)
c            endif
c
c         endif
cc
c 100        continue
cc     
c         enddo
cc
cc         write(6,*) 'star_feedback5: end'
c         return
c         end
