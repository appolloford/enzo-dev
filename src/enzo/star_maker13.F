#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c////////////////////////  SUBROUTINE STAR_MAKER \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_maker13(nx, ny, nz,
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
c  modified1: enzo-dev (2.4) implementation BW, Feb 2016 (star_maker12)
c  modified1: enzo-dev (2.5) implemented by Duncan Christie, Oct 2016 (this version)
c 
c  Based on star_maker2.src (though with similar stochastic SF to star_maker1), 
c  this routine restricts star formation criteria to the density threshold; 
c  idea from Krumholz & Tan (2007) who argue a fixed SFR/t_ff in turbulant gas. 
c  Since gas must also be molecular, T < 3000 K. 
c
c  Addtion for star_maker13:
c
c  If the gas is magneto-Jeans unstable, the routine behaves like star_maker12
c  If the gas is magneto-Jeans stable but Jeans unstable, it inserts a collapse retardation factor. 
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
c    rcut - causes radial cutoff for star particles (avoids central 2.5 kpc of disk). 
c    Bx, By, Bz - the magnetic field in code units.
c
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
c  Locals:
c
      INTG_PREC  i, j, k, ii, clock
      R_PREC   div, tdyn, dtot
      R_PREC   pi, G, sndspdC
      R_PREC   isosndsp2, starmass, starfraction, bmass, jeanmass
      R_PREC   Bmag, Bregconst, Bcrit, rho, cs, magunits,crf
      real*8 msolar
      parameter (pi=pi_val, G=GravConst, 
     &           sndspdC=1.3095e8_RKIND,
     &           msolar=SolarMass)

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

      magunits = dsqrt(4._RKIND*pi*d1*v1**2)

      do k=1+ibuff,nz-ibuff
         do j=1+ibuff,ny-ibuff
            do i=1+ibuff,nx-ibuff

c
c              1) finest level of refinement?
c
               if (r(i,j,k) .ne. 0._RKIND) goto 10

c
c              2) density greater than threshold
c
               if (d(i,j,k) .lt. abs(odthresh)) goto 10
                
c
c Calculate the critical magnetic field for dynamic collapse
c
                rho = d(i,j,k)*d1
                cs = dsqrt(1.38e-16_RKIND*temp(i,j,k)/(2.33_RKIND*1.6733e-24_RKIND))
                Bcrit = 4._RKIND*pi*rho*(dx*dx*x1*x1*rho*G/pi - cs*cs)
                
c
c There is potentially weird behaviour here if you want to create stars from gas that is Jeans-stable. It won't let you.
c

                if (Bcrit .lt. 0._RKIND) goto 10

                Bcrit = dsqrt(Bcrit)
                
c              Magnitude of B in physical units:
               Bmag = dsqrt(Bx(i,j,k)**2 + By(i,j,k)**2
     &               + Bz(i,j,k)**2)*magunits

c               Collapse retardation factor
               crf = 1._RKIND

              if (Bmag .gt. Bcrit) crf = 6.0 

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
               starfraction = min(masseff*dt/(crf*tdyn), 0.9_RKIND)
               tdyn = max(tdyn, mintdyn*3.15e7_RKIND/t1)


#ifdef STOCHASTIC_STAR_FORMATION

c        if mass of star is less than mass threshold, only creates star
c        if (some random #) < starmass/massthresh 
c        and then create it with massthresh, unless not enough gas in cell.

      
               bmass = d(i,j,k)*dble(d1)*dble(x1*dx)**3 / msolar

               if (starfraction*bmass .lt. smthresh) then
c                  write(6,*) 'star_maker5: stochastic sf'
                  call random_number(random)
                  x = random
                  if (x .gt. (starfraction*bmass / smthresh) ) goto 10 
                  starfraction = min(smthresh/bmass, 0.9)
               endif

               write(*,*) 'I have beaten the probability ', 
     &                    starfraction*bmass,crf

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


