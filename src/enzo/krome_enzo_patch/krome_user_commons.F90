
!############### MODULE ##############
module krome_user_commons
  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2021-02-02 01:45:04
  !  Changeset xxxxxxx
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  ! mock parameters, only for single grid model, change it!!
  real*8 :: gridsize, sputteringrate
  logical :: startr = .true.
  real*8 :: ebmaxh2=1.21d3,epsilon=0.01,ebmaxcrf=1.21d3,uvcreff=1.0d-3, &
      &  ebmaxcr=1.21d3,phi=1.0d5,ebmaxuvcr=1.0d4,uvy=0.1,h2form=0.0
  !Variables for self-shielding of CO and H2
  !dopw = doppler width (in s-1) of a typical transition
  !(assuming turbulent broadening with beta=3e5cms-1)
  !fosc = oscillator strength of a typical transition
  real*8 :: fosc = 1.d-2, dopw = 3.d10, xl = 1000.0, radw = 8.d7
  integer,parameter :: dimco=7, dimh2=6
  real*8 :: corates(7,6)=reshape((/0.000d+00, -1.408d-02, -1.099d-01, -4.400d-01,&
      &  -1.154d+00, -1.888d+00, -2.760d+00,&
      &  -8.539d-02, -1.015d-01, -2.104d-01, -5.608d-01,&
      &  -1.272d+00, -1.973d+00, -2.818d+00,&
      &  -1.451d-01, -1.612d-01, -2.708d-01, -6.273d-01,&
      &  -1.355d+00, -2.057d+00, -2.902d+00,&
      &  -4.559d-01, -4.666d-01, -5.432d-01, -8.665d-01,&
      &  -1.602d+00, -2.303d+00, -3.146d+00,&
      &  -1.303d+00, -1.312d+00, -1.367d+00, -1.676d+00,&
      &  -2.305d+00, -3.034d+00, -3.758d+00,&
      &  -3.883d+00, -3.888d+00, -3.936d+00, -4.197d+00,&
      &  -4.739d+00, -5.165d+00, -5.441d+00 /),shape(corates))
  real*8 :: y2r(7,6)
  real*8 :: ncogr(dimco) =(/12.0d+00, 13.0d+00, 14.0d+00, 15.0d+00,&
      &16.0d+00, 17.0d+00, 18.0d+00 /)
  real*8 :: nh2gr(dimh2)=(/18.0d+00, 19.0d+00, 20.0d+00, 21.0d+00,&
      &22.0d+00, 23.0d+00 /)
  real*8 :: nHAv(91,2) = reshape( (/ &
      -3.0d0, -2.9d0, -2.8d0, -2.7d0, -2.6d0, -2.5d0, -2.4d0, &
      -2.3d0, -2.2d0, -2.1d0, -2.0d0, -1.9d0, -1.8d0, -1.7d0, &
      -1.6d0, -1.5d0, -1.4d0, -1.3d0, -1.2d0, -1.1d0, -1.0d0, &
      -0.9d0, -0.8d0, -0.7d0, -0.6d0, -0.5d0, -0.4d0, -0.3d0, &
      -0.2d0, -0.1d0,  0.0d0,  0.1d0,  0.2d0,  0.3d0,  0.4d0, &
      0.5d0,  0.6d0,  0.7d0,  0.8d0,  0.9d0,  1.0d0,  1.1d0, &
      1.2d0,  1.3d0,  1.4d0,  1.5d0,  1.6d0,  1.7d0,  1.8d0, &
      1.9d0,  2.0d0,  2.1d0,  2.2d0,  2.3d0,  2.4d0,  2.5d0, &
      2.6d0,  2.7d0,  2.8d0,  2.9d0,  3.0d0,  3.1d0,  3.2d0, &
      3.3d0,  3.4d0,  3.5d0,  3.6d0,  3.7d0,  3.8d0,  3.9d0, &
      4.0d0,  4.1d0,  4.2d0,  4.3d0,  4.4d0,  4.5d0,  4.6d0, &
      4.7d0,  4.8d0,  4.9d0,  5.0d0,  5.1d0,  5.2d0,  5.3d0, &
      5.4d0,  5.5d0,  5.6d0,  5.7d0,  5.8d0,  5.9d0,  6.0d0, &
      3.85d-03, 3.90d-03, 3.97d-03, 4.06d-03, 4.18d-03, 4.32d-03, 4.48d-03, 4.69d-03, &
      4.96d-03, 5.29d-03, 5.71d-03, 6.22d-03, 6.85d-03, 7.65d-03, 8.65d-03, 9.91d-03, &
      1.13d-02, 1.30d-02, 1.51d-02, 1.78d-02, 2.12d-02, 2.33d-02, 2.60d-02, 2.93d-02, &
      3.36d-02, 3.89d-02, 4.41d-02, 5.07d-02, 5.90d-02, 6.95d-02, 8.27d-02, 9.47d-02, &
      1.10d-01, 1.29d-01, 1.53d-01, 1.83d-01, 2.02d-01, 2.25d-01, 2.55d-01, 2.92d-01, &
      3.39d-01, 3.76d-01, 4.22d-01, 4.81d-01, 5.55d-01, 6.48d-01, 7.06d-01, 7.78d-01, &
      8.68d-01, 9.83d-01, 1.13d+00, 1.21d+00, 1.31d+00, 1.44d+00, 1.60d+00, 1.80d+00, &
      1.93d+00, 2.08d+00, 2.28d+00, 2.54d+00, 2.85d+00, 3.04d+00, 3.28d+00, 3.58d+00, &
      3.95d+00, 4.42d+00, 4.70d+00, 5.05d+00, 5.48d+00, 6.03d+00, 6.72d+00, 7.12d+00, &
      7.61d+00, 8.24d+00, 9.03d+00, 1.00d+01, 1.06d+01, 1.13d+01, 1.22d+01, 1.33d+01, &
      1.47d+01, 1.54d+01, 1.64d+01, 1.76d+01, 1.92d+01, 2.11d+01, 2.22d+01, 2.35d+01, &
      2.52d+01, 2.73d+01, 3.00d+01 &
      /), shape(nHAv))

  ! reference: definition of h2col/ch2 and cocol
  ! h2col=0.5*abund(nh2,dstep)*density(dstep)*(cloudSize/real(points))
  ! cocol=0.5*abund(nco,dstep)*density(dstep)*(cloudSize/real(points))

  !user can add here the common variables needed
  !for rate calculation (e.g. optical depth, CR rate,
  !pressure, density, ...)
  real*8::tau,zrate,pah_size,gas_dust_ratio,krome_redshift
  real*8::krome_J21,krome_J21xr,user_tff,user_omega,user_nu
  real*8::krome_invdvdz !inverse of |velocity gradient| 1/abs(dv/dz)

  !$omp threadprivate(user_tff,krome_invdvdz)

contains

  subroutine set_gridsize(cellsize)
    implicit none
    real*8, intent(in) :: cellsize
    gridsize = cellsize
  end subroutine set_gridsize

  function get_gridsize()
    implicit none
    real*8 :: get_gridsize
    get_gridsize = gridsize
  end function get_gridsize

  subroutine set_sputtering(sputtering)
    implicit none
    real*8, intent(in) :: sputtering
    sputteringrate = sputtering
  end subroutine

  function get_sputtering()
    implicit none
    real*8 :: get_sputtering
    get_sputtering = sputteringrate
  end function

  !**********************
  !user can add here the functions he/she needs for
  !rate calculations (Kooij funcion provided as example)
  function kooij(kalpha,kbeta,kgamma,Tgas)
    real*8::kooij
    real*8,intent(in)::kalpha,kbeta,kgamma,Tgas
    kooij = kalpha*(Tgas/3d2)**kbeta*exp(-kgamma/Tgas)
  end function kooij

  !**********************
  !prototype for photo xsec kernel function.
  !this is called by xsec interpolator when resampling
  ! xsec produced by python into photobins at runtime
  function fxinterp(energy)
    implicit none
    real*8,intent(in)::energy
    real*8::fxinterp

    fxinterp = energy**2

  end function fxinterp

  ! The following functions must be included for uclchem network
  function get_Av(Hnuclei)
    implicit none
    real*8, intent(in) :: Hnuclei
    real*8 :: logH, get_Av
    integer :: i

    logH = log10(Hnuclei)
    if (logH >= 6.0d0) then
      get_Av = nHAv(91,2)
      return
    end if

    if (logH <= -3.0d0) then
      get_Av = nHAv(1,2)
      return
    end if

    do i = 1, 90
      if (nHAv(i,1) < logH .and. nHAv(i+1,1) >= logH) then
        get_Av = ( nHAv(i+1,2)*(logH-nHAv(i,1)) + nHAv(i,2)*(nHAv(i+1,1)-logH) ) &
            / ( nHAv(i+1,1) - nHAv(i,1) )
        return
      end if
    end do

  end function

  ! The following functions are copied from UCLCHEM v1.3

  !-----------------------------------------------------------------------------------------
  !h2d() and knrco() recalculate h2 and CO photodissociation using a self-shielding
  !treatment from van dishoeck and black (apj 334, p771 (1988))
  !All other functions are helpers for this purpose
  !Code by Serena Viti
  !-----------------------------------------------------------------------------------------

  function h2d(n)
    use krome_commons
    use krome_getphys
    implicit none
    real*8 :: n(:),ch2, taud, h2d, Hnuclei, av

    !h2col is h2 column density. Half total column density for uniform sphere.
    !Sum of shells for multidepth point (worked out in chem.f90.updateChemistry)

    !taud = opt. depth at line centre (assum. ortho:parah2=1)
    !pi**0.5 * e2 / (m(electr) * c) = 1.5e-2 cm2/s
    ! ch2 = 0.5*n(idx_H2)*get_Hnuclei(n(:))*gridsize
    ch2 = 0.5*n(idx_H2)*gridsize
    taud  = 0.5 * ch2 * 1.5e-2 * fosc / dopw
    Hnuclei = get_Hnuclei(n(:))
    av  = get_Av(Hnuclei)

    !c Here, the constant 5.1e-11 is by assuming rad is in Habing [double check]
    h2d = 5.1d-11 * user_rad * scat(xl, av) * fgk(taud)
  end function h2d

  function knrco(n)
    use krome_commons
    use krome_getphys
    implicit none
    real*8 :: n(:),ch2,ssf,lba,av,sca,chtot,knrco
    real*8 :: Hnuclei,h2col,cocol
    !double precision :: ssfco,lbar,scat,chtot

    ! h2col = 0.5*n(idx_H2)*get_Hnuclei(n(:))*gridsize
    ! cocol = 0.5*n(idx_H2)*get_Hnuclei(n(:))*gridsize
    h2col = 0.5*n(idx_H2)*gridsize
    cocol = 0.5*n(idx_CO)*gridsize
    !calculate photodissociation rates for co (species # nco; reaction
    !# nrco) according to van dishoeck and black (apj 334, p771 (1988))
    !cocol is the co column density (in cm-2); the scaling of the pdrate
    !by a factor of 1.8 is due to the change from draine's is uv field
    !to the habing field
    ssf = ssfco(n)
    lba = lbar(cocol,h2col)
    Hnuclei = get_Hnuclei(n(:))
    av  = get_Av(Hnuclei)
    sca = scat(lba, av)

    !The reason why rad is divided by 1.7 is that the alphas are for Draine and the rad is in
    !Habing units

    knrco = (2.d-10) * user_rad/1.7 * ssf * sca
  end function knrco

  function fgk(taudin)
    implicit none
    !calculates the line self shielding function
    !described in federman et al. apj vol.227 p.466.

    !--------------------------------------------------------------
    !input parameters
    !dopw : doppler   line width (in s-1) from sr setpmp
    !radw : radiative line width (in s-1) from sr setpmp
    !taud : optical depth at line center  from sr setequ
    !
    !program variables
    !fgk  : total self shielding function containing the radiative
    !       and the doppler contribution.
    !r    : parameter r  (eq. a2) in federman's paper
    !sj   : parameter jd (eq. a8) in federman's paper
    !       doppler contribution to self shielding function
    !sr   : parameter jr (eq. a9) in federman's paper
    !       radiative contribution to self shielding function
    !t    : parameter t1 (eq. a6) in federman's paper
    !u    : parameter u1 (eq. a6) in federman's paper
    !-----------------------------------------------------------

    !program variables type declaration
    real*8, intent(in) :: taudin
    real*8 :: r, sj, sr, t, u, taud, fgk
    !--------------------------------------------------------------

    !calculate wing contribution of self shielding function sr
    if (taudin .lt. 0.0d0) then
      taud = 0.0d0
    else
      taud = taudin
    endif

    if (radw .eq. 0.0d0) then
      sr = 0.0d0
    else
      r  = radw/(1.7724539d0*dopw)
      t  = 3.02d0 * ((r*1.0d+03)**(-0.064d0))
      u  = ((taud*r)**0.5d0)/t
      sr = ((0.78539816d0+u**2)**(-0.5d0))*r/t
    endif

    !calculate doppler contribution of self shielding function sj
    if (taud .eq. 0.0d0) then
      sj = 1.0d0
    else if (taud .lt. 2.0d0) then
      sj = exp(-0.6666667d0*taud)
    else if (taud .lt. 10.0d0) then
      sj = 0.638d0*taud**(-1.25d0)
    else if (taud .lt. 100.0d0) then
      sj = 0.505d0*taud**(-1.15d0)
    else
      sj = 0.344d0*taud**(-1.0667d0)
    endif

    !calculate total self shielding function fgk
    fgk = sj + sr
  end function fgk

  function scat(x1, av)
    use krome_commons
    implicit none
    !calculate the influence of dust extinction (g=0.8, omega=0.3)
    !wagenblast&hartquist, mnras237, 1019 (1989)

    !---------------------------------------------------------------------
    !         i/o variables type declaration
    !       scat   : factor describing the influence of grain scattering
    !                on the uv flux dependent on the total h number
    !                density and wavelength of the uv radiation
    !       x1      : wavelength (in angstrom)
    !       cdntot : total h number density (in cm-2)
    !
    !         program variables
    !       av     : visual extinction in magnitudes (cf. savage et al.,
    !                 1977 apj 216, p.291)
    !        expo   : exponent
    !        i      : loop index
    !        tl     : tau(lambda)
    !        tv     : tau(visual=5500a)
    !        xlamda : function which calculates tl/tv
    !        c(0)   : c(0) * exp(-k(0)*tau) : (rel.) intensity
    !                 decrease for 0<=tau<=1 caused by grain
    !                 scattering with g=0.8, omega=0.3
    !                 (approximation)
    !        c(i)   : sum (c(i) * exp(-k(i)*tau)) i=1,5  (rel.)
    !                 intensity decrease for 1<=tau<=oo caused by
    !                 grain scattering with g=0.8, omega=0.3.
    !                 (cf. flannery, roberge, and rybicki 1980,
    !                 apj 236, p.598).
    !        k(0)   : see c0
    !        k(i)   : see ci
    !---------------------------------------------------------------------

    integer :: i
    !   i/o variables type declaration
    real*8 :: cdntot,x1,av
    intent(in) :: x1, av
    !program variables type declaration
    real*8, dimension(6) :: c=(/1.0d0,2.006d0,-1.438d0,7.364d-01,-5.076d-01,-5.920d-02/)
    real*8, dimension(6) ::  k1=(/7.514d-01,8.490d-01,1.013d0,1.282d0,2.005d0,5.832d0/)
    real*8 :: expo, tl, tv, scat

    !calculate visual extinction
    !total column density of hydrogen nuclei / e(b-v) = 5.8e+21
    !atoms cm-2 mag-1 (bohlin, savage, and drake 1978; apj 224,132)
    !for lambda**-1 scattering : r = av / e(b-v) = 3.6 .

    !optical depth in the visual
    tv = av/ 1.086d0

    !make correction for the wavelength considered
    tl = tv * xlamda(x1)

    !calculate scat
    scat = 0.0d0
    if (tl.lt.1.0d0) then
      expo = k1(1)*tl
      if (expo.lt.35.0d0) then
        scat = c(1) * dexp(-expo)
      endif
    else
      do i=2,6
        expo = k1(i)*tl
        if (expo.lt.35.0d0) then
          scat = scat + c(i)*dexp(-expo)
        endif
      enddo
    endif

  end function scat

  function xlamda(x)
    implicit none
    !calculate  xlamda : =  tau(lambda) / tau(visual)
    !
    ! --------------------------------------------------------------
    !         i/o parameter
    !xlamda : =  tau(lambda) / tau(visual);
    !         tau(lambda) is the opt. depth for dust extinction at
    !         wavelength x (cf. b.d.savage and j.s.mathis, annual
    !         review of astronomy and astrophysics vol.17(1979),p.84)
    !x      : wavelength in angstrom
    ! --------------------------------------------------------------

    !i/o parameter type declaration
    real*8 :: x, xlamda
    intent(in) :: x

    if (x.lt.  910.0d0) then
      xlamda = 5.76d0
    else if (x.lt.  950.0d0) then
      xlamda = 5.76d0 - 1.45d-02*(x-  910.0d0)
    else if (x.lt. 1000.0d0) then
      xlamda = 5.18d0 - 1.06d-02*(x-  950.0d0)
    else if (x.lt. 1050.0d0) then
      xlamda = 4.65d0 - 9.68d-03*(x- 1000.0d0)
    else if (x.lt. 1110.0d0) then
      xlamda = 4.16d0 - 7.26d-03*(x- 1050.0d0)
    else if (x.lt. 1180.0d0) then
      xlamda = 3.73d0 - 4.61d-03*(x- 1110.0d0)
    else if (x.lt. 1250.0d0) then
      xlamda = 3.40d0 - 4.15d-03*(x- 1180.0d0)
    else if (x.lt. 1390.0d0) then
      xlamda = 3.11d0 - 2.67d-03*(x- 1250.0d0)
    else if (x.lt. 1490.0d0) then
      xlamda = 2.74d0 - 1.10d-03*(x- 1390.0d0)
    else if (x.lt. 1600.0d0) then
      xlamda = 2.63d0 - 8.80d-05*(x- 1490.0d0)
    else if (x.lt. 1700.0d0) then
      xlamda = 2.62d0 - 8.06d-04*(x- 1600.0d0)
    else if (x.lt. 1800.0d0) then
      xlamda = 2.54d0 - 3.87d-04*(x- 1700.0d0)
    else if (x.lt. 1900.0d0) then
      xlamda = 2.50d0 + 8.07d-04*(x- 1800.0d0)
    else if (x.lt. 2000.0d0) then
      xlamda = 2.58d0 + 2.00d-03*(x- 1900.0d0)
    else if (x.lt. 2100.0d0) then
      xlamda = 2.78d0 + 2.29d-03*(x- 2000.0d0)
    else if (x.lt. 2190.0d0) then
      xlamda = 3.01d0 + 1.22d-03*(x- 2100.0d0)
    else if (x.lt. 2300.0d0) then
      xlamda = 3.12d0 - 2.35d-03*(x- 2190.0d0)
    else if (x.lt. 2400.0d0) then
      xlamda = 2.86d0 - 2.81d-03*(x- 2300.0d0)
    else if (x.lt. 2500.0d0) then
      xlamda = 2.58d0 - 2.29d-03*(x- 2400.0d0)
    else if (x.lt. 2740.0d0) then
      xlamda = 2.35d0 - 1.46d-03*(x- 2500.0d0)
    else if (x.lt. 3440.0d0) then
      xlamda = 2.00d0 - 5.99d-04*(x- 2740.0d0)
    else if (x.lt. 4000.0d0) then
      xlamda = 1.58d0 - 2.88d-04*(x- 3440.0d0)
    else if (x.lt. 4400.0d0) then
      xlamda = 1.42d0 - 2.42d-04*(x- 4000.0d0)
    else if (x.lt. 5500.0d0) then
      xlamda = 1.32d0 - 2.93d-04*(x- 4400.0d0)
    else if (x.lt. 7000.0d0) then
      xlamda = 1.00d0 - 1.68d-04*(x- 5500.0d0)
    else if (x.lt. 9000.0d0) then
      xlamda = 0.75d0 - 1.32d-04*(x- 7000.0d0)
    else if (x.lt.12500.0d0) then
      xlamda = 0.48d0 - 5.81d-05*(x- 9000.0d0)
    else if (x.lt.22000.0d0) then
      xlamda = 0.28d0 - 1.66d-05*(x-12500.0d0)
    else if (x.lt.34000.0d0) then
      xlamda = 0.12d0 - 5.91d-06*(x-22000.0d0)
    else
      xlamda = 0.05d0 - 5.16d-11*(x-34000.0d0)
    endif

  end function xlamda

  function ssfco(n)
    use krome_commons
    use krome_getphys
    use krome_fit
    implicit none
    !calculates self-shielding factors for 12co transitions due to
    !12co self-shielding and h2 screening. values given in table 5
    !of van dishoeck and black, apj 334, p771 (1988) are used for
    !2dim spline interpolation. nco and nh2 are the 12co resp. h2
    !column densities for which the self-shielding factor is
    !interpolated.

    !common parameter
    !ssfcor : corates(7,6) : self-shielding factors for 12co
    !                     considering h2 line overlap according
    !                     to van dishoeck and black, apj 334, p771
    !                     (1988). logarithmic (base 10) values for
    !                     the self-shielding factors are stored.
    !                     1st index : variation of co column dens.
    !                     2nd index : variation of h2 column dens.
    !y2r(7,6) :   2nd derivative of rates from sr splie2
    !ncogr(7) :   12co column density grid (log10 values
    !             of column densities in cm-2)
    !nh2gr(6) :   h2 column density grid (log10 values
    !             of column densities in cm-2)
    !dimco :      dimension of ncogr
    !dimh2 :      dimension of nh2gr
    !startr :     is .true. when ssfcor is entered first
    !actual values are supplied by block data routine

    !program variables type declaration
    real*8 :: n(:), ch2, lognco, lognh2, ssfco, cocol
    intent(in) :: n
    if (startr) then
      call splie2_ucl(ncogr,nh2gr,corates,dimco,dimh2,y2r)
      startr = .false.
    endif

    ! ch2 = 0.5*n(idx_H2)*get_Hnuclei(n(:))*gridsize
    ! cocol = 0.5*n(idx_H2)*get_Hnuclei(n(:))*gridsize
    ch2 = 0.5*n(idx_H2)*gridsize
    cocol = 0.5*n(idx_CO)*gridsize

    lognco = dlog10(max(cocol,1.0d5))
    lognh2 = dlog10(max(ch2,1.0d10))

    if (lognco.lt.ncogr(1))      lognco = ncogr(1)
    if (lognh2.lt.nh2gr(1))      lognh2 = nh2gr(1)
    if (lognco.gt.ncogr(dimco))  lognco = ncogr(dimco)
    if (lognh2.gt.nh2gr(dimh2))  lognh2 = nh2gr(dimh2)

    call splin2_ucl(ncogr,nh2gr,corates,y2r,dimco,dimh2,lognco,&
        &               lognh2,ssfco)
    ssfco = 10.0d0**ssfco
  end function ssfco

  function lbar(u,w)
    implicit none
    !calculate lambda bar (in a) according to equ. 4 of van dishoeck
    !and black, apj 334, p771 (1988)
    ! --------------------------------------------------------------
    !       i/o parameter
    !       u : co column density in (cm-2)
    !       w : h2 column density in (cm-2)

    !        program variables
    !        lu : log10(co column density in cm-2)
    !        lw : log10(h2 column density in cm-2)

    !--------------------------------------------------------------
    !i/o parameter type declaration
    real*8 :: u, w, lu, lw, lbar

    lu = dlog10(dabs(u)+1.0d0)
    lw = dlog10(dabs(w)+1.0d0)

    lbar = (5675.0d0 - 200.6d0*lw) - (571.6d0 - 24.09d0*lw)*lu + &
        & (18.22d0 - 0.7664d0*lw)*lu**2

    !lbar represents the mean of the wavelengths of the 33
    !dissociating bands weighted by their fractional contribution
    !to the total rate of each depth. lbar cannot be larger than
    !the wavelength of band 33 (1076.1a) and not be smaller than
    !the wavelength of band 1 (913.6a).
    if (lbar.gt.1076.1d0)  lbar = 1076.1d0
    if (lbar.lt. 913.6d0)  lbar =  913.6d0
  end function lbar

  function desoh2(gama, mantle)
    use krome_commons
    use krome_getphys
    implicit none
    real*8 :: ebmaxh2 = 1.21d3
    real*8 :: gama, mantle, desoh2

    desoh2 = 0.0
    if (user_desorb == 0.d0 .or. user_h2desorb == 0.d0) return
    if (gama .le. ebmaxh2 .and. mantle .gt. 1d-30) desoh2 = 1.0
  end function desoh2

  function descr(gama, mantle)
    use krome_commons
    use krome_getphys
    implicit none
    real*8 :: ebmaxcr = 1.21d3
    real*8 :: gama, mantle, descr

    descr = 0.0
    if (user_desorb == 0.d0 .or. user_crdesorb == 0.d0) return
    if (gama .le. ebmaxcr .and. mantle .gt. 1d-30) descr = 1.0
  end function descr

  function deuvcr(gama, mantle)
    use krome_commons
    use krome_getphys
    implicit none
    real*8 :: ebmaxuvcr=1.0d4
    real*8 :: gama, mantle, deuvcr

    deuvcr = 0.0
    if (user_desorb == 0.d0 .or. user_uvcr == 0.d0) return
    if (gama .le. ebmaxuvcr .and. mantle .gt. 1d-30) deuvcr = 1.0
  end function deuvcr

  function freeze(Tgas)
    use krome_commons
    implicit none
    real*8 :: Tgas, freeze

    freeze = 0.0
    if (user_fr == 1.d0 .and. Tgas .le. 30.0) freeze = 1.0
  end function


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Mock subroutine for cshock that will sputter the ices based on Jimenez-Serra!
  ! 2008 paper.                                                                 !
  !
  ! Modified from UCLCEHM
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! function sputterrate(n)

  !   use krome_commons
  !   use krome_constants
  !   use krome_getphys
  !   implicit none
  !   real*8 :: sputterrate, n(:), mass(nspec), nh, s, coeff, mass_p, n_p, mantles, grainNumberDensity
  !   real*8, parameter :: GAS_DUST_NUMBER_RATIO=1.14d-12
  !   real*8, parameter :: grainRadius = 1.0d-5
  !   integer :: ispec
  !   integer :: projectiles(6)=(/idx_H2, idx_HE, idx_C, idx_O, idx_SI, idx_CO/)

  !   mass = get_mass()
  !   mantles = get_mantle(n(:))
  !   nh = get_Hnuclei(n(:))

  !   ! loop over projectile species and get rates of change of mantle for each, summing them
  !   sputterrate = 0.0
  !   do ispec = 1, size(projectiles) !!!! Make projectiles array in initialize
  !     ! projectile mass and number density
  !     mass_p = mass(projectiles(ispec))
  !     n_p = n(projectiles(ispec)) !* n_h
  !     ! leading coefficient
  !     coeff = dsqrt(8.0 * pi * boltzmann_erg * n(idx_Tgas) / mass_p) * grainRadius**2
  !     ! Variable relating mass and speed of projectile to energy, used in integration
  !     s = dsqrt(mass_p * gasvel**2 / (2.0*n(idx_Tgas)*boltzmann_erg))
  !     sputterrate = sputterrate + coeff * n_p * iceYieldRate(mass_p, n(idx_Tgas), s)
  !   end do

  !   grainNumberDensity = nh * GAS_DUST_NUMBER_RATIO
  !   ! Total rate/cm3 (ie released particles /cm3/s) is sputterRate (per grain) multiplied by grain number density
  !   sputterrate = sputterrate * grainNumberDensity
  !   sputterrate = min(sputterrate, mantles)


  ! end function

  function iceYield(specidx, vel)
    use krome_commons
    use krome_getphys
    implicit none

    real*8, parameter :: tarmass = 18.0 * 1.67353251819d-24 ! water ice
    real*8, parameter :: icebind = 0.53 * 1.6d-12
    real*8, parameter :: effic = 0.8
    real*8, parameter :: yieldconst = 8.3d-4

    integer, intent(in) :: specidx
    real*8, intent(in) :: vel
    real*8 :: yield, eta, eps, eps0, projmass, m(nspec)
    real*8 :: iceYield

    m(:) = get_mass()
    projmass = m(specidx)

    eta = 4.0 * effic * projmass * tarmass * (projmass+tarmass)**(-2.0)
    eps = eta * 0.5 * projmass * vel**2.0 / icebind
    eps0 = max(1.0, 4.0*eta)

    iceYield = 0.0
    if (eps > eps0) then
      iceYield = 2.0 * yieldconst * (eps-eps0)**2 / (1.0 + (eps/30.0)**(1.3333))
    endif

  end function iceYield


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Function calculates rate of change of ice mantle abundance of a species!
  !due to the impact of molecules of a given mass. actual rate is         !
  !proportional to projectile abundance                                   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function iceYieldRate(projectileMass, temp, s)
    use krome_constants
    implicit none
    real*8 :: iceYieldRate, projectileMass
    real*8 :: lowerLimit, upperLimit, temp, s, eps0, eta

    real*8, parameter :: iceBindingEnergy = 0.53*1.6d-12
    real*8, parameter :: targetMass = 18.0*1.67353251819d-24
    real*8, parameter :: iceYieldEfficiency = 0.8

    ! eta is effectively reduced mass of the collision
    eta = 4.0 * iceYieldEfficiency * projectileMass * targetMass * (projectileMass+targetMass)**(-2.0)
    eps0 = max(1.0, 4.0*eta)

    !Lower limit is xth in Jimenez-Serra et al. 2008
    lowerLimit = dsqrt(eps0 * iceBindingEnergy / (eta*boltzmann_erg*temp) )

    !Upper limit is just where the integrand goes to zero
    upperLimit = iceYieldIntegralLimit(lowerLimit, projectileMass, temp, s, eta)

    !calculate eq B.1 from Jimenez-Serra et al. 2008
    if ((upperlimit-lowerLimit) .gt. 1d-4) then
      !first get integral from Eq B.1 including 1/s factor
      iceYieldRate=trapezoidIntegrate(iceYieldIntegrand, lowerLimit, upperLimit, projectileMass, temp, s, eta)/s
    else
      iceYieldRate=0.0
    end if
  end function

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Function calculates integrand from Eq B.1 of Jimenez-Serra et al. 2008 !
  !                                                                       !
  !Inputs are mass of projectile and x. Returns value of integrand at x   !
  !allowing trapezium rule to integrate from xth to infinity              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function iceYieldIntegrand(x, mass_p, temp, s, eta)
    use krome_constants
    implicit none
    real*8 :: iceYieldIntegrand, x, mass_p, temp, s, eta, eps0, eps, yield

    real*8, parameter :: yieldConst = 8.3d-4
    real*8, parameter :: iceBindingEnergy = 0.53*1.6d-12

    ! epsilon is calculated from inmpact energy (Ep)
    eps0 = max(1.0, 4.0*eta)
    eps = (x**2) * boltzmann_erg * temp * eta / iceBindingEnergy

    ! this yield is for ice averaged over all angles. There's a different one for cores
    ! (Appendix B Jimenez-Serra 2008)
    ! it's 2 times the normal incidence yield, but there's a factor of 0.5 in integrand so we drop both
    yield = yieldConst * (eps-eps0)**2 / (1.0 + (eps/30.0)**(1.3333) )
    iceYieldIntegrand = yield * (x**2) * ( dexp(-(x-s)**2) - DEXP(-(x+s)**2) )
  end function iceYieldIntegrand

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Function to calculate the upper limit beyond which there's no point   !
  !evaluating the ice yield integrand. Ie trapezoids from upper limit to !
  !upperlimit+dx will have an area~0                                     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function iceYieldIntegralLimit(xth, mass_p, temp, s, eta)
    implicit none
    real*8 :: iceYieldIntegralLimit, xth, mass_p, temp, s, eta
    integer :: i
    i=1
    ! Take upperlimit to be half way between lowerLimit and 1000.
    iceYieldIntegralLimit = xth + (1d3-xth) * (0.5**i)
    ! decrease upper limit for as long as f(upperlimit) is <1.0e-20 and
    ! difference between lower and upper limit is not zero.
    do while ( &
          iceYieldIntegralLimit-xth .gt. 1.0d-3 .and. &
          iceYieldIntegrand(iceYieldIntegralLimit, mass_p, temp, s, eta) .lt. 1d-20 &
          )
      i = i+1
      iceYieldIntegralLimit = xth + (1d3-xth) * (0.5**i)
    end do
  end function iceYieldIntegralLimit

  ! Subroutine that calculates an integral using the trapezoidal method.
  function trapezoidIntegrate(func, xl, xu, mass_p, temp, s, eta)
    implicit none
    real*8 :: trapezoidIntegrate, func, xl, xu, mass_p, temp, s, eta, olds
    real*8, parameter :: tol = 1.0e-3
    integer :: j
    integer, parameter :: JMAX = 25
    external :: func

    olds = -1.0e30
    do j = 1, JMAX
      call trapzd(func, xl, xu, trapezoidIntegrate, j, mass_p, temp, s, eta)

      if (abs(trapezoidIntegrate-olds) .le. tol*abs(olds)) then
        return
      end if

      olds=trapezoidIntegrate
    end do
  end function

  subroutine trapzd(func, xl, xu, ret, iter, mass_p, temp, s, eta)
    real*8 :: func, xl, xu, ret, mass_p, temp, s, eta
    integer :: iter, it, j
    DOUBLE PRECISION del,sum,tnm,x
    external func
    if (iter .eq. 1) then
      ret = 0.5 * (xu-xl) * (func(xl, mass_p, temp, s, eta) + func(xu, mass_p, temp, s, eta))
    else
      it = 2**(iter-2)
      tnm = it
      del = (xu-xl)/tnm
      x = xl + 0.5*del
      sum = 0.0
      do j = 1, it
        sum = sum + func(x, mass_p, temp, s, eta)
        x = x+del
      end do
      ret = 0.5*(ret + (xu-xl) * sum / tnm)
    end if
    return
  end subroutine

end module krome_user_commons

