
!############### MODULE ##############
module krome_main

  integer::krome_call_to_fex
  !$omp threadprivate(krome_call_to_fex)

contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2019-12-02 01:47:36
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

  !********************************
  !KROME main (interface to the solver library)

  subroutine krome(x,Tgas,dt  )
    use krome_commons
    use krome_subs
    use krome_ode
    use krome_reduction
    use krome_dust
    use krome_getphys
    use krome_tabs
    implicit none
    real*8 :: Tgas,dt
    real*8 :: x(nmols)
    real*8 :: rhogas

    real*8::mass(nspec),n(nspec),tloc,xin
    real*8::rrmax,totmass,n_old(nspec),ni(nspec),invTdust(ndust)
    integer::icount,i,icount_max
    integer:: ierr

    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer::neq(1),itol,itask,istate,iopt,lrw,liw,mf
    integer::iwork(2621)
    real*8::atol(nspec),rtol(nspec)
    real*8::rwork(19792)
    logical::got_error,equil

    !****************************
    !init DLSODES (see DLSODES manual)
    call XSETF(0)!toggle solver verbosity
    got_error = .false.
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are scalar
    rtol(:) = 1.000000d-04 !relative tolerance
    atol(:) = 1.000000d-20 !absolute tolerance
    icount_max = 100 !maximum number of iterations

    itask = 1
    iopt = 0

    !MF=
    !  = 222 internal-generated JAC and sparsity
    !  = 121 user-provided JAC and internal generated sparsity
    !  =  22 internal-generated JAC but sparsity user-provided
    !  =  21 user-provided JAC and sparsity
    MF = 222
    !end init DLSODES
    !****************************

    ierr = 0 !error flag, zero==OK!
    n(:) = 0d0 !initialize densities

    n(1:nmols) = x(:)

    n(idx_Tgas) = Tgas !put temperature in the input array

    icount = 0 !count solver iterations
    istate = 1 !init solver state
    tloc = 0.d0 !set starting time

    !store initial values
    ni(:) = n(:)
    n_global(:) = n(:)

    n_old(:) = -1d99
    krome_call_to_fex = 0
    do
      icount = icount + 1
      !solve ODE
      CALL DLSODES(fex, NEQ(:), n(:), tloc, dt, &
          ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
          LIW, JES, MF)

      krome_call_to_fex = krome_call_to_fex + IWORK(12)
      !check DLSODES exit status
      if(istate==2) then
        exit !sucsessful integration
      elseif(istate==-1) then
        istate = 1 !exceeded internal max iterations
      elseif(istate==-5 .or. istate==-4) then
        istate = 3 !wrong sparsity recompute
      elseif(istate==-3) then
        n(:) = ni(:)
        istate = 1
      else
        got_error = .true.
      end if

      if(got_error.or.icount>icount_max) then
        if (krome_mpi_rank>0) then
          print *,krome_mpi_rank,"ERROR: wrong solver exit status!"
          print *,krome_mpi_rank,"istate:",istate
          print *,krome_mpi_rank,"iter count:",icount
          print *,krome_mpi_rank,"max iter count:",icount_max
          print *,krome_mpi_rank,"SEE KROME_ERROR_REPORT file"
        else
          print *,"ERROR: wrong solver exit status!"
          print *,"istate:",istate
          print *,"iter count:",icount
          print *,"max iter count:",icount_max
          print *,"SEE KROME_ERROR_REPORT file"
        end if
        call krome_dump(n(:), rwork(:), iwork(:), ni(:))
        stop
      end if

    end do

    !avoid negative species
    do i=1,nspec
      n(i) = max(n(i),0d0)
    end do

    !returns to user array
    x(:) = n(1:nmols)

    Tgas = n(idx_Tgas) !get new temperature

  end subroutine krome

  !*********************************
  !integrates to equilibrium using constant temperature
  subroutine krome_equilibrium(x,Tgas,verbosity)
    use krome_ode
    use krome_subs
    use krome_commons
    use krome_constants
    use krome_getphys
    use krome_tabs
    implicit none
    integer::mf,liw,lrw,itol,meth,iopt,itask,istate,neq(1)
    integer::i,imax
    integer,optional::verbosity
    integer::verbose
    real*8 :: Tgas
    real*8 :: x(nmols)
    real*8 :: rhogas
    real*8::tloc,n(nspec),mass(nspec),ni(nspec)
    real*8::dt,xin
    integer::iwork(2621)
    real*8::atol(nspec),rtol(nspec)
    real*8::rwork(19792)
    real*8::ertol,eatol,max_time,t_tot,ntot_tol,err_species
    logical::converged

    integer, save :: ncall=0
    integer, parameter :: ncall_print_frequency=20000
    integer :: ncallp
    integer::charges(nspec)
    real*8::masses(nspec)
    character*16::names(nspec)

    !set verbosity from argument
    verbose = 1 !default is verbose
    if(present(verbosity)) verbose = verbosity

    call XSETF(0)!toggle solver verbosity
    meth = 2
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are scalar
    rtol(:) = 1d-6 !relative tolerance
    atol(:) = 1d-20 !absolute tolerance

    ! Switches to decide when equilibrium has been reached
    ertol = 1d-5  ! relative min change in a species
    eatol = 1d-12 ! absolute min change in a species
    max_time=seconds_per_year*5d8 ! max time we will be integrating for

    !for DLSODES options see its manual
    iopt = 0
    itask = 1
    istate = 1

    mf = 222 !internally evaluated sparsity and jacobian
    tloc = 0d0 !initial time

    n(:) = 0d0 !initialize densities
    !copy into array
    n(nmols+1:) = 0d0
    n(1:nmols) = x(:)

    n(idx_Tgas) = Tgas

    !store previous values
    ni(:) = n(:)
    n_global(:) = ni(:)

    imax = 1000

    dt = seconds_per_year * 1d2
    t_tot = dt
    converged = .false.
    do while (.not. converged)
      do i=1,imax
        !solve ODE
        CALL DLSODES(fcn_tconst, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, jcn_dummy, MF)
        if(istate==2) then
          exit
        else
          istate=1
        end if
      end do
      !check errors
      if(istate.ne.2) then
        print *,"ERROR: no equilibrium found!"
        stop
      end if

      !avoid negative species
      do i=1,nspec
        n(i) = max(n(i),0d0)
      end do

      ! check if we have converged by comparing the error in any species with an relative abundance above eatol
      converged = maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) .lt. ertol &
          .or. t_tot .gt. max_time

      ! Increase integration time by a reasonable factor
      if(.not. converged) then
        dt = dt * 3.
        t_tot = t_tot + dt
        ni = n
        n_global = n
      endif
    enddo
    !returns to user array
    x(:) = n(1:nmols)

    if(t_tot > max_time .and. &
        maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) > 0.2 .and. verbose>0) then
    print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years.'
    print *, 'Tgas :', Tgas
    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a16)',"#","Name","m (g)","Chrg","  Current / Last"
    do i=1,nmols
      print '(I4,a10,E11.3,I5,2E14.6,E11.3)',i," "//names(i),masses(i),charges(i),n(i),ni(i),abs(n(i) - ni(i)) / max(n(i),eatol*sum(n(1:nmols)))
    end do
    print '(a30,2E14.6)'," sum",sum(n(1:nmols)),sum(ni(1:nmols))
    print *, 'Fractional error :', maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols))))
    print *, 'Absolute and relative floors:', eatol, ertol
  end if

  ! Print info ever so often
  !$omp critical
  ncall=ncall+1
  ncallp = ncall
  !$omp end critical

  if(modulo(ncallp,ncall_print_frequency)==0 .and. verbose>0) then
    print *, 'Found equilibrium for ', ncallp, ' cells.'
  end if

end subroutine krome_equilibrium

!********************
!dummy jacobian
subroutine jcn_dummy()
  implicit none
end subroutine jcn_dummy

!*******************
!dn/dt where dT/dt=0
subroutine fcn_tconst(n,tt,x,f)
  use krome_commons
  use krome_ode
  implicit none
  integer::n,ierr
  real*8::x(n),f(n),tt
  call fex(n,tt,x(:),f(:))
  f(idx_Tgas) = 0d0
end subroutine fcn_tconst

!*******************************
subroutine krome_dump(n,rwork,iwork,ni)
  use krome_commons
  use krome_subs
  use krome_tabs
  use krome_reduction
  use krome_ode
  use krome_getphys
  integer::fnum,i,iwork(:),idx(nrea),j
  real*8::n(:),rwork(:),rrmax,k(nrea),kmax,rperc,kperc,dn(nspec),tt,ni(:)
  character*16::names(nspec),FMTi,FMTr
  character*50::rnames(nrea),fname,prex
  integer,save::mx_dump=1000 ! max nr of reports before terminating
  fnum = 99
  if (krome_mpi_rank>0) then
    write(fname,'(a,i5.5)') "KROME_ERROR_REPORT_",krome_mpi_rank
  else
    fname = "KROME_ERROR_REPORT"
  endif
  open(fnum,FILE=trim(fname),status="replace")
  tt = 0d0
  names(:) = get_names()
  rnames(:) = get_rnames()
  call fex(nspec,tt,n(:),dn(:))

  write(fnum,*) "KROME ERROR REPORT"
  write(fnum,*)
  !SPECIES
  write(fnum,*) "Species abundances"
  write(fnum,*) "**********************"
  write(fnum,'(a5,a20,3a12)') "#","name","qty","dn/dt","ninit"
  write(fnum,*) "**********************"
  do i=1,nspec
    write(fnum,'(I5,a20,3E12.3e3)') i,names(i),n(i),dn(i),ni(i)
  end do
  write(fnum,*) "**********************"

  !F90 FRIENDLY RESTART
  write(fnum,*)
  write(fnum,*) "**********************"
  write(fnum,*) "F90-friendly species"
  write(fnum,*) "**********************"
  do i=1,nspec
    write(prex,'(a,i3,a)') "x(",i,") = "
    write(fnum,*) trim(prex),ni(i),"!"//names(i)
  end do

  write(fnum,*) "**********************"

  !RATE COEFFIECIENTS
  k(:) = coe_tab(n(:))
  idx(:) = idx_sort(k(:))
  kmax = maxval(k)
  write(fnum,*)
  write(fnum,*) "Rate coefficients (sorted) at Tgas",n(idx_Tgas)
  write(fnum,*) "**********************"
  write(fnum,'(a5,2a12,a10)') "#","k","k %","  name"
  write(fnum,*) "**********************"
  do j=1,nrea
    i = idx(j)
    kperc = 0.d0
    if(kmax>0.d0) kperc = k(i)*1d2/kmax
    write(fnum,'(I5,2E12.3e3,a2,a50)') i,k(i),kperc,"  ", rnames(i)
  end do
  write(fnum,*) "**********************"
  write(fnum,*)

  !FLUXES
  call load_arrays
  rrmax = fex_check(n(:), n(idx_Tgas))
  idx(:) = idx_sort(arr_flux(:))
  write(fnum,*)
  write(fnum,*) "Reaction magnitude (sorted) [k*n1*n2*n3*...]"
  write(fnum,*) "**********************"
  write(fnum,'(a5,2a12,a10)') "#","flux","flux %","  name"
  write(fnum,*) "**********************"
  do j=1,nrea
    i = idx(j)
    rperc = 0.d0
    if(rrmax>0.d0) rperc = arr_flux(i)*1d2/rrmax
    write(fnum,'(I5,2E12.3e3,a2,a50)') i,arr_flux(i),rperc,"  ",rnames(i)
  end do
  write(fnum,*) "**********************"
  write(fnum,*)

  !SOLVER
  FMTr = "(a30,E16.7e3)"
  FMTi = "(a30,I10)"
  write(fnum,*) "Solver-related information:"
  write(fnum,FMTr) "step size last",rwork(11)
  write(fnum,FMTr) "step size attempt",rwork(12)
  write(fnum,FMTr) "time current",rwork(13)
  write(fnum,FMTr) "tol scale factor",rwork(14)
  write(fnum,FMTi) "numeber of steps",iwork(11)
  write(fnum,FMTi) "call to fex",iwork(12)
  write(fnum,FMTi) "call to jex",iwork(13)
  write(fnum,FMTi) "last order used",iwork(14)
  write(fnum,FMTi) "order attempt",iwork(15)
  write(fnum,FMTi) "idx largest error",iwork(16)
  write(fnum,FMTi) "RWORK size required",iwork(17)
  write(fnum,FMTi) "IWORK size required",iwork(18)
  write(fnum,FMTi) "NNZ in Jac",iwork(19)
  write(fnum,FMTi) "extra fex to compute jac",iwork(20)
  write(fnum,FMTi) "number of LU decomp",iwork(21)
  write(fnum,FMTi) "base address in RWORK",iwork(22)
  write(fnum,FMTi) "base address of IAN",iwork(23)
  write(fnum,FMTi) "base address of JAN",iwork(24)
  write(fnum,FMTi) "NNZ in lower LU",iwork(25)
  write(fnum,FMTi) "NNZ in upper LU",iwork(21)
  write(fnum,*) "See DLSODES manual for further details on Optional Outputs"
  write(fnum,*)
  write(fnum,*) "END KROME ERROR REPORT"
  write(fnum,*)
  close(fnum)

  mx_dump = mx_dump - 1
  if (mx_dump==0) stop

end subroutine krome_dump

!********************************
subroutine krome_init()
  use krome_commons
  use krome_tabs
  use krome_subs
  use krome_reduction
  use krome_dust
  use krome_cooling
  use krome_photo
  use krome_fit

  !init phys common variables
  !$omp parallel
  phys_Tcmb = 2.73d0
  phys_zredshift = 0d0
  phys_orthoParaRatio = 3d0
  phys_metallicity = 0d0
  phys_Tfloor = 2.73d0
  !$omp end parallel

  !init metallicity default
  !assuming solar
  total_Z = 1d0

  !default D/D_sol = Z/Z_sol
  !assuming linear scaling
  dust2gas_ratio = total_Z

  !default broadening turubulence velocity
  broadeningVturb2 = 0d0

  !default clumping factor for
  ! H2 formation on dust by Jura/Gnedin
  clump_factor = 1d0

  !default for thermo toggle is ON
  !$omp parallel
  krome_thermo_toggle = 1
  !$omp end parallel

  !load arrays with ractants/products indexes
  call load_arrays()

  !initialize the table for exp(-a/T) function
  call init_exp_table()

  call load_parts()

  !init photo reactants indexes

  !get machine precision
  krome_epsilon = epsilon(0d0)

  !load verbatim reactions
  call loadReactionsVerbatim()

end subroutine krome_init

!****************************
function krome_get_coe(x,Tgas)
  !krome_get_coe: public interface to obtain rate coefficients
  use krome_commons
  use krome_subs
  use krome_tabs
  implicit none
  real*8 :: krome_get_coe(nrea), x(nmols), Tgas
  real*8::n(nspec)

  n(:) = 0d0
  n(1:nmols) = x(:)
  n(idx_Tgas) = Tgas
  krome_get_coe(:) = coe_tab(n(:))

end function krome_get_coe

!****************************
function krome_get_coeT(Tgas)
  !krome_get_coeT: public interface to obtain rate coefficients
  ! with argument Tgas only
  use krome_commons
  use krome_subs
  use krome_tabs
  implicit none
  real*8 :: krome_get_coeT(nrea),Tgas
  real*8::n(nspec)
  n(idx_Tgas) = Tgas
  krome_get_coeT(:) = coe_tab(n(:))
end function krome_get_coeT

end module krome_main
