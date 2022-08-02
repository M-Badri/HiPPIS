program main
!!
!! driver for BOMEX test case
!!

  implicit none

  !! 1D advection code used for debugging !!
  !!call advection1d()

  call bomex1()
  

end program main


subroutine bomex1()
!!
!! This subroutine is run the BOMEX example with the 
!! different confidgurations presented in the manuscript. 
!!  

  implicit none

  integer, parameter :: r8 = SELECTED_REAL_KIND(8)

  integer           :: nlevs			!! number of level used for simulation
  real(kind=r8)     :: cfl			!! clf condtion used
  integer           :: i, jj, k, ii		!! integer used for iteration
  character(len=16) :: snlevs			!! number of levels (string)
  character(len=16) :: scfl			!! cfl conditions (string)
  character(len=16) :: bomex_type		!! method used for advection in bomex (string)
  character(len=16) :: degree		        !! maximum polynomial degree used for interpolation (stringg)
  character(len=16) :: sst 			!! to choose the stencil selection procedure st=1 for ENO, st=2 for favoring symmetry, st=3 for favoring locality   
  character(len=16) :: seps0 			!! to specify eps for the PPI method
  character(len=16) :: seps1 			!! to specify eps for the PPI method
  character(len=16) :: mapping_type 		!! to specify the interpolation method Linear, Standard, PCHIP, DBI, and PPI
  


  nlevs = 600           !! 60, 180, 300, 600, 900, 1200 
  cfl = 0.1             !! 0.6, 0.3, 0.1, 0.01
  seps0 = "1.0e-5"      !! 1.0e-0, 1.0e-1, 2.0e-2, 2.0e-3, 2.0e-4, 2.0e-5, 2.0e-6, 0.0e-0
  seps1 = "1.0e-5"      !! 1.0e-0, 1.0e-1, 2.0e-2, 2.0e-3, 2.0e-4, 2.0e-5, 2.0e-6, 0.0e-0
  sst = "1"		!! 1, 2, 3
  bomex_type = "weno"   !! or "linear" or "cubic" or "weno"
  scfl = "_1"           !! CFL codition (cfl=0.1)
  snlevs = "600"        !! number of levels "nlevs"

    
    do jj=1, 3

      !! Set up different choices for stencil selction procedure
      !!   sst="1" ENO-like stencil selection procedure
      !!   sst="2" stencil selection that favor symetry
      !!   sst="3" stencil selection that favor locality
      if(jj==1)then
        sst ="1"
      elseif(jj==2)then
        sst = "2"
      elseif(jj==3)then
        sst = "3"
      endif

      write(*,*) 'nelevs=',nlevs, 'cfl =', cfl
      do ii=1,4
        !! Polynomial degree to be used for the simulation 
        if(ii.eq.1)then
          degree = "3"
        elseif(ii.eq.2)then
          degree = "5"
        elseif(ii.eq.3)then
          degree = "7"
        elseif(ii.eq.4)then
          degree = "9"
        elseif(ii.eq.5)then
          degree = "11"
        elseif(ii.eq.6)then
          degree = "13"
        endif

        mapping_type = "DBI" !! DBI (data-bounded interpolation)
        write(*,*) 'BOMEX simulation using DBI to map solution values between the physics and dynamics meshes. Max degree =', degree
        call bomex_mapping(nlevs, cfl, snlevs, scfl, bomex_type, mapping_type, degree, sst, seps0, seps1)

        mapping_type = "PPI" !! PPI (positivity preserving interpolation)
        write(*,*) 'BOMEX simulation using PPI to map solution values between the physics and dynamics meshes. &
                    Max degree =', degree
        call bomex_mapping(nlevs, cfl, snlevs, scfl, bomex_type, mapping_type, degree, sst, seps0, seps1)
      enddo

    enddo
    write(*,*) 'BOMEX simulation with the same mesh used for both physics and dynamics calculations'
    call bomex_no_mapping(nlevs, cfl, snlevs, scfl, bomex_type)

    write(*,*) 'BOMEX simulation using PCHIP to map solution values between the physics and dynamics meshes.'
    mapping_type = "PCHIP" 
    call bomex_mapping(nlevs, cfl, snlevs, scfl, bomex_type, mapping_type, degree, sst, seps0, seps1)

    !!! for interval I_{i} the stencil is V_4 = \{ x_{i-2}, x_{i-1}, x_{i}, x_{i+1}, x_{i+2}, x_{i+3} \}
    write(*,*) 'BOMEX simulation using a fifth order standar  polynomial interpolation to map solution values &
                between the physics and dynamics meshes.'
    mapping_type= "Standard" 
    call bomex_mapping(nlevs, cfl, snlevs, scfl, bomex_type, mapping_type, degree, sst, seps0, seps1)
   
    !!! for interval I_{i} the stencil is V_4 = \{ x_{i-2}, x_{i-1}, x_{i}, x_{i+1}, x_{i+2}, x_{i+3} \}
    write(*,*) 'BOMEX simulation using a fifth order standard polynomial interpolation with clipping &
                to map solution values between the physics and dynamics meshes.'
    mapping_type= "Clipping" 
    call bomex_mapping(nlevs, cfl, snlevs, scfl, bomex_type, mapping_type, degree, sst, seps0, seps1)

    !!write(*,*) 'BOMEX simulation using linear interpolation to map solution &
    !!            values between the physics and dynamics meshes.'
    !!mapping_type = "Linear" 
    !!call bomex_mapping(nlevs, cfl, snlevs, scfl, bomex_type, mapping_type, degree, sst, seps0, seps1)


  
end subroutine bomex1 

subroutine bomex_mapping(nz, cfl, snlevs, scfl, bomex_type, mapping_type, sdegree, sst, seps0, seps1)
!!
!! This subroutine run the BOMEX simulation with different meshes used for the 
!! dynamics and physics calculations. The dynamcis is calculated using 
!! uniformly-spaced points that indicate the boundary of each level. The physics 
!! mesh is calculated using the mid-point of each level. 
!!
!! When mapping from the physics to dynamics mesh the top and bottom boundry values 
!! using the results from the BOMEX simulation where the same mesh is used for both
!! the dynamics and physics mesh.  
!!
!! INPUT:
!! nz: number of levels
!! cfl: CFL conditions used for the simulation
!! snlevs: the number levels (string)
!! scfl: CFL conditions used for simulation (string)
!! bomex_type: method used for the advections in the dynamics. The options are ``Linear", ``Cubic", and ``WENO"
!! sdegree : maximu polynomial degree used for each interval (string)
!!

  use mod_bomex, only : bomex_init, set_bomex_forcing, bomex_ls_forcing
  use scm_physics, only : run_scm_physics
  use machine, only : r8 => real_kind
  use mod_adaptiveInterpolation, only : adaptiveInterpolation1D
  
  implicit none

  integer, intent(in)            :: nz				!! number of levels 
  real(kind=r8), intent(in)      :: cfl				!! CFL condition
  character(len=16), intent(in)  :: snlevs			!! number of levels (string)
  character(len=16), intent(in)  :: scfl			!! CFL condition (string)
  character(len=16), intent(in)  :: bomex_type			!! method used for adevction in BOMEX (Linear or WENO5)
  character(len=16), intent(in)  :: sdegree			!! maximum polynomial degree for each interval (string)
  character(len=16), intent(in)  :: mapping_type		!! type of interpolation method used Linear, PCHIP, DBI, PPI
  character(len=16), intent(in)  :: sst
  character(len=16), intent(in)  :: seps0
  character(len=16), intent(in)  :: seps1



  !! Variables specifict to interpolation methods !!
  integer                        :: degree			!! maximum polynomial degree for each interval
  integer                        :: limiter			!! type of interpolation
  integer			 :: st				!! to choose stencil construction procedure 
  integer, dimension(nz)         :: deg 			!! polynomial degree used for each interpolant
  real(kind=r8)                  :: eps0			!! eps
  real(kind=r8)                  :: eps1			!! eps

  !! Variable specific to BOMEX simulation !!
  integer                        :: n_steps			!! to track current time step
  integer                        :: start_map			
  integer                        :: i, k, kk, fid

  real(kind=r8), parameter       :: ztop = 3000.0_r8		!! top boundary in (m)
  real(kind=r8), dimension(nz)   :: u, v, th, prs, qv, qc, qr, rho, z
  real(kind=r8), dimension(nz+1) :: ui, vi, thi, prsi, prsi0, prsi2, qvi, qci, qri, rhoi,  zi
  real(kind=r8), dimension(nz+1) :: ug, vg, w, dthdt, dqvdt
  real(kind=r8), dimension(nz)   :: u2, v2, th2, prs2, qv2, qc2, qr2, rho2
  real(kind=r8), dimension(nz)   :: ug2, vg2, w2, dthdt2, dqvdt2
  real(kind=r8)                  :: ps2, heat2, evap2, stress2
  real(kind=r8)                  :: ps, dz, tf, t0, heat, evap, stress
  real(kind=r8)                  :: lat                         ! Only relevant for test=1
  real(kind=r8)                  :: precl
  real(kind=r8)                  :: dt				!! time step size
  real(kind=r8)                  :: dt_write			!! used to write simulation results to file
  logical                        :: lneg
  character(len=64)              :: fname			!! output file name 

  !! variablesspecific to PCHIP 
  integer                        :: nwk, nwki
  integer                        :: ierr
  real(kind=r8), dimension(nz+1)          :: fdi, d_tmpi
  real(kind=r8), dimension((nz+1)*2)      :: wki
  real(kind=r8), dimension(nz)            :: fd, d_tmp
  real(kind=r8), dimension(nz*2)          :: wk
  logical                                 :: spline 
  
  !! Variables specific to standard interpolation !!
  integer                                 :: ks, ke


  nwk = nz*2
  nwki = (nz+1)*2
  spline = .false.
  
  dt_write = 0.0_r8
  
  
  !! Use input information for build file name 
  !! For example  if bomex_type="weno", mapping_type = "PPI", snlevs="600", 
  !! scfl=_1, sdegree="5", st=1,eps0=_1, eps1=_1
  !!  bomexweno600cfl=_1PPI05st=1eps0=_01eps1=_01.dat
  if( mapping_type == "PPI" .or. mapping_type == "DBI") then
    fname =trim("bomex_data/")//trim("bomex")//trim(bomex_type)//trim(snlevs)//trim(scfl)&
         //trim(mapping_type)//trim(sdegree)//trim('st=')//trim(sst)&
         //trim("eps0=")//trim(seps0)//trim("eps1=")//trim(seps1)//trim(".dat")
  elseif( mapping_type == "DBI") then
    fname =trim("bomex_data/")//trim("bomex")//trim(bomex_type)//trim(snlevs)//trim(scfl)&
         //trim(mapping_type)//trim(sdegree)//trim("st=")//trim(sst)//trim(".dat")
  elseif(mapping_type == "Linear" .or. mapping_type == "Standard" .or. &
         mapping_type == "PCHIP" .or. mapping_type == "Clipping") then
    fname =trim("bomex_data/")//trim("bomex")//trim(bomex_type)//trim(snlevs)//trim(scfl)&
         //trim(mapping_type)//trim(".dat")
  else
    write(*,*) 'ERROR: Incorrect mappying_type the only option available for  &
               mapping types are Linear, Standard, Clipping, PCHIP, DBI and PPI'
  endif

  !! get maximum polynomial degree for each interval
  if(sdegree .eq. "1")then
    degree = 1
  elseif(sdegree .eq. "2")then
    degree = 2
  elseif(sdegree .eq. "3")then
    degree = 3
  elseif(sdegree .eq. "4")then
    degree = 4
  elseif(sdegree .eq. "5")then
    degree = 5
  elseif(sdegree .eq. "6")then
    degree = 6
  elseif(sdegree .eq. "7")then
    degree = 7
  elseif(sdegree .eq. "8")then
    degree = 8
  elseif(sdegree .eq. "9")then
    degree = 9
  elseif(sdegree .eq. "10")then
    degree = 10
  elseif(sdegree .eq. "11")then
    degree = 11
  elseif(sdegree .eq. "12")then
    degree = 12
  elseif(sdegree .eq. "13")then
    degree = 13
  elseif(sdegree .eq. "14")then
    degree = 14
  elseif(sdegree .eq. "15")then
    degree = 15
  elseif(sdegree .eq. "16")then
    degree = 16
  else
    write(*,*) "ERROR: Invalid degree ", sdegree
    write(*,*) "sdegree musbe between 1 and 16"
  endif
  
  !! get limiter
  if(mapping_type .eq. "DBI") then
    limiter = 1
  elseif(mapping_type .eq. "PPI")then
    limiter = 2
  !!else
  !!  write(*,*) "ERROR: Invalid limiter"
  !!  write(*,*) "limiter is a required parameter that must be DBI or PPI"
  endif

  !! set stencil construction procedure
  if(sst == "1")then
    st = 1
  elseif(sst == "2")then
    st = 2
  elseif(sst == "3")then
    st = 3
  endif

  !! Set eps0  
  if(seps0 == "1.0e-0")then
    eps0 = 1.0e-0
  elseif(seps0 == "1.0e-1")then
    eps0 = 1.0e-1
  elseif(seps0 == "1.0e-2")then
    eps0 = 1.0e-2
  elseif(seps0 == "1.0e-3")then
    eps0 = 1.0e-3
  elseif(seps0 == "1.0e-4")then
    eps0 = 1.0e-4
  elseif(seps0 == "1.0e-5")then
    eps0 = 1.0e-5
  elseif(seps0 == "1.0e-6")then
    eps0 = 1.0e-6
  elseif(seps0 == "0.0e-0")then
    eps0 = 0.0
  else
    Write(*,*)"WARNING: eps0 has not be set by the default values eps0=1.0e-2 will be used"
    eps0 = 1.0e-2
  endif

  !! Set eps1  
  if(seps1 == "1.0e-0")then
    eps1 = 1.0e-0
  elseif(seps1 == "1.0e-1")then
    eps1 = 1.0e-1
  elseif(seps1 == "1.0e-2")then
    eps1 = 1.0e-2
  elseif(seps1 == "1.0e-3")then
    eps1 = 1.0e-3
  elseif(seps1 == "1.0e-4")then
    eps1 = 1.0e-4
  elseif(seps1 == "1.0e-5")then
    eps1 = 1.0e-5
  elseif(seps1 == "1.0e-6")then
    eps1 = 1.0e-6
  elseif(seps1 == "0.0e-0")then
    eps1 = 0.0
  else
    Write(*,*)"RNING: eps1 has not be set by the default values eps1=1.0 will be used"
    eps1 = 1.0
  endif


  t0=0.0_r8
  tf=6*3600.0_r8
  dt = cfl *3000.0_r8/(real(nz, kind=r8))
  write(*,*) 'dt= ', dt
  lat = 0.0_r8

  start_map = 3600*6
  !!start_map = 15!!3600*3
  n_steps = 0
  
  !--- Evenly Spaced Grid ---!
  dz = ztop/nz
  zi(1) = 0.0
  do k=2,nz+1
    zi(k) = zi(k-1) + dz          
    z(k-1) = 0.5*(zi(k) + zi(k-1))
  end do


  ! Set bomex initial state
  call bomex_init(u,v,rho,th,qv,qc,qr,prs,z,ps,nz)

  call bomex_init(u2,v2,rho2,th2,qv2,qc2,qr2,prs2,z,ps2,nz)

  ! Only need prsi from this call
  call bomex_init(ui,vi,rhoi,thi,qvi,qci,qri,prsi,zi,ps,nz+1)
  prsi0 = prsi
  prsi2 = prsi

  call set_bomex_forcing(ug,vg,w,dthdt,dqvdt,heat,evap,stress,zi,nz+1)

  call set_bomex_forcing(ug2,vg2,w2,dthdt2,dqvdt2,heat2,evap2,stress2,z,nz)

  do while(t0<tf)
  
    n_steps= n_steps+1
    call bomex_ls_forcing(u2,v2,prs2,rho2,th2,qv2,qc2,qr2, &
                          ug2,vg2,w2,dthdt2,dqvdt2,z,dt,nz, bomex_type)
    call run_scm_physics(u2,v2,th2,qv2,qc2,qr2,prs2,prsi2,rho2,heat2,evap2,stress2,z,zi,dt,nz)

    !!if(n_steps > start_map )then
    !! "dynamics" caluclations 
    call bomex_ls_forcing(ui,vi,prsi,rhoi,thi,qvi,qci,qri, &
                          ug,vg,w,dthdt,dqvdt,zi,dt,nz+1, bomex_type)

    !!** Convert from "dynamics" grid to "physics" grid using the DBI or PPI method **!! 
    if(mapping_type == "PPI" .or. mapping_type == "DBI") then
      call adaptiveInterpolation1D(zi, ui, nz+1, z, u, nz, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(zi, vi, nz+1, z, v, nz, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(zi, prsi, nz+1, z, prs, nz, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(zi, rhoi, nz+1, z, rho, nz, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(zi, thi, nz+1, z, th, nz, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(zi, qvi, nz+1, z, qv, nz, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(zi, qci, nz+1, z, qc, nz, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(zi, qri, nz+1, z, qr, nz, degree, limiter, st, eps0, eps1, deg)

    !!** Convert from "dynamics" grid to "physics" grid using standard interpolation **!! 
    elseif(mapping_type == "Standard" .or. mapping_type == "Clipping") then
      do k=1, nz
        ks= max(1,k-2)
        ke= ks+5
        if(ke>nz+1)then
          ks=nz-4
          ke=nz+1
        endif
        call lagrangePolyVal(zi(ks:ke), ui(ks:ke), 6, z(k), u(k))
        call lagrangePolyVal(zi(ks:ke), vi(ks:ke), 6, z(k), v(k))
        call lagrangePolyVal(zi(ks:ke), rhoi(ks:ke), 6, z(k), rho(k))
        call lagrangePolyVal(zi(ks:ke), prsi(ks:ke), 6, z(k), prs(k))
        call lagrangePolyVal(zi(ks:ke), thi(ks:ke), 6, z(k), th(k))
        call lagrangePolyVal(zi(ks:ke), qvi(ks:ke), 6, z(k), qv(k))
        call lagrangePolyVal(zi(ks:ke), qci(ks:ke), 6, z(k), qc(k))
        call lagrangePolyVal(zi(ks:ke), qri(ks:ke), 6, z(k), qr(k))
      enddo
      if(mapping_type == "Clipping") then
        do k=1, nz
          if(qc(k) < 0.0 )then
            qc(k) = 0.0
          endif
        enddo
      endif

    !!** Convert from "dynamics" grid to "physics" grid using linear interpolation **!! 
    elseif(mapping_type == "Linear") then
      do k=1, nz
        u(k) = 0.5 * (ui(k+1) + ui(k))
        v(k) = 0.5 * (vi(k+1) + vi(k))
        prs(k) = 0.5 * (prsi(k+1) + prsi(k))
        rho(k) = 0.5 * (rhoi(k+1) + rhoi(k))
        th(k) = 0.5 * (thi(k+1) + thi(k))
        qv(k) = 0.5 * (qvi(k+1) + qvi(k))
        qc(k) = 0.5 * (qci(k+1) + qci(k))
        qr(k) = 0.5 * (qri(k+1) + qri(k))
      enddo

    !!** Convert from "dynamics" grid to "physics" grid using PCHIP **!! 
    elseif(mapping_type == "PCHIP") then
      call pchez(nz+1, zi, ui, d_tmpi, spline, wki, nwki, ierr)
      call pchev(nz+1, zi, ui, d_tmpi, nz, z, u, fd, ierr)

      call pchez(nz+1, zi, vi, d_tmpi, spline, wki, nwki, ierr)
      call pchev(nz+1, zi, vi, d_tmpi, nz, z, v, fd, ierr)

      call pchez(nz+1, zi, rhoi, d_tmpi, spline, wki, nwki, ierr)
      call pchev(nz+1, zi, rhoi, d_tmpi, nz, z, rho, fd, ierr)

      call pchez(nz+1, zi, thi, d_tmpi, spline, wki, nwki, ierr)
      call pchev(nz+1, zi, thi, d_tmpi, nz, z, th, fd, ierr)

      call pchez(nz+1, zi, prsi, d_tmpi, spline, wki, nwki, ierr)
      call pchev(nz+1, zi, prsi, d_tmpi, nz, z, prs, fd, ierr)

      call pchez(nz+1, zi, qvi, d_tmpi, spline, wki, nwki, ierr)
      call pchev(nz+1, zi, qvi, d_tmpi, nz, z, qv, fd, ierr)

      call pchez(nz+1, zi, qci, d_tmpi, spline, wki, nwki, ierr)
      call pchev(nz+1, zi, qci, d_tmpi, nz, z, qc, fd, ierr)

      call pchez(nz+1, zi, qri, d_tmpi, spline, wki, nwki, ierr)
      call pchev(nz+1, zi, qri, d_tmpi, nz, z, qr, fd, ierr)

    endif 

    
    !!if(n_steps >= start_map )then
    !!if (dt_write >=  30.0 .or. dt_write .eq. tf)then
    !! Wrire data to file at time t0 = 18000
    if (t0 <= 18000.0 .and. t0 + dt > 18000.0 ) then
      fid=10
      open(unit=fid,file=fname,status='unknown')
      do k=1,nz
        write(fid,'(9(3x,E12.5))') z(k),th(k),u(k),v(k),qv(k),qc(k),qr(k)
      end do
      close(fid)
    endif
    !!  dt_write = 0.0
    !!endif
    !!dt_write = dt_write + dt
    !!endif

    !!** "physics" calculations **!!
    call run_scm_physics(u,v,th,qv,qc,qr,prs,prsi0,rho,heat,evap,stress,z,zi,dt,nz)

    !!   the profile rho and prs are not interpolated back because they have not been modified by 
    !!   physics routines  

    !!** Convert from "physics" grid to "dynamics" grid using the DBI or the PPI method **!!
    if(mapping_type == "PPI" .or. mapping_type == "DBI") then
      call adaptiveInterpolation1D(z, u, nz, zi(2:nz), ui(2:nz), nz-1, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(z, v, nz, zi(2:nz), vi(2:nz), nz-1, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(z, th, nz, zi(2:nz), thi(2:nz), nz-1, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(z, qv, nz, zi(2:nz), qvi(2:nz), nz-1, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(z, qc, nz, zi(2:nz), qci(2:nz), nz-1, degree, limiter, st, eps0, eps1, deg)
      call adaptiveInterpolation1D(z, qr, nz, zi(2:nz), qri(2:nz), nz-1, degree, limiter, st, eps0, eps1, deg)

    !!** Convert from "physics" grid to "dynamics" grid using a standard interpolation **!!
    elseif(mapping_type == "Standard" .or. mapping_type == "Clipping") then
      do k=2, nz
        ks= max(1,k-2)
        ke= ks+5
        if(ke>nz)then
         ks=nz-5
         ke=nz
        endif
        call lagrangePolyVal(z(ks:ke), u(ks:ke), 6, zi(k), ui(k))
        call lagrangePolyVal(z(ks:ke), v(ks:ke), 6, zi(k), vi(k))
        call lagrangePolyVal(z(ks:ke), th(ks:ke), 6, zi(k), thi(k))
        call lagrangePolyVal(z(ks:ke), qv(ks:ke), 6, zi(k), qvi(k))
        call lagrangePolyVal(z(ks:ke), qc(ks:ke), 6, zi(k), qci(k))
        call lagrangePolyVal(z(ks:ke), qr(ks:ke), 6, zi(k), qri(k))
      enddo
      if(mapping_type == "Clipping") then
        do k=2, nz
          if(qc(k) < 0.0)then
            qc(k) = 0.0
          endif
        enddo
      endif

    !!** Convert from "physics" grid to "dynamics" grid using linear interpolation **!!
    elseif(mapping_type == "Linear") then
      do k=2, nz
        ui(k) = 0.5 * (u(k) + u(k-1))
        vi(k) = 0.5 * (v(k) + v(k-1))
        thi(k) = 0.5 * (th(k) + th(k-1))
        qvi(k) = 0.5 * (qv(k) + qv(k-1))
        qci(k) = 0.5 * (qc(k) + qc(k-1))
        qri(k) = 0.5 * (qr(k) + qr(k-1))
      enddo

    !!** Convert from "physics" grid to "dynamics" grid using PCHIP **!!
    elseif(mapping_type == "PCHIP") then
      call pchez(nz, z, u, d_tmp, spline, wk, nwk, ierr)
      call pchev(nz, z, u, d_tmp, nz-1, zi(2:nz), ui(2:nz), fdi, ierr)

      call pchez(nz, z, v, d_tmp, spline, wk, nwk, ierr)
      call pchev(nz, z, v, d_tmp, nz-1, zi(2:nz), vi(2:nz), fdi, ierr)

      call pchez(nz, z, th, d_tmp, spline, wk, nwk, ierr)
      call pchev(nz, z, th, d_tmp, nz-1, zi(2:nz), thi(2:nz), fdi, ierr)

      call pchez(nz, z, qv, d_tmp, spline, wk, nwk, ierr)
      call pchev(nz, z, qv, d_tmp, nz-1, zi(2:nz), qvi(2:nz), fdi, ierr)

      call pchez(nz, z, qc, d_tmp, spline, wk, nwk, ierr)
      call pchev(nz, z, qc, d_tmp, nz-1, zi(2:nz), qci(2:nz), fdi, ierr)

      call pchez(nz, z, qr, d_tmp, spline, wk, nwk, ierr)
      call pchev(nz, z, qr, d_tmp, nz-1, zi(2:nz), qri(2:nz), fdi, ierr)
    endif
 

    !! Extrapolate bottom values using solution from  target profile
    !! where the same mesh is used for both the dynamics and physics 
    !! calculations.
    ui(1)   = u2(1)     
    vi(1)   = v2(1)     
    thi(1)  = th2(1)    
    qvi(1)  = qv2(1)    
    qci(1)  = qc2(1)    
    qri(1)  = qr2(1)    

    !! Extrapolate top values using solution from  target profile
    !! where the same mesh is used for both the dynamics and physics 
    !! calculations.
    ui(nz+1)    = u2(nz)     
    vi(nz+1)    = v2(nz)     
    thi(nz+1)   = th2(nz)  
    qvi(nz+1)   = qv2(nz)  
    qci(nz+1)   = qc2(nz)  
    qri(nz+1)   = qr2(nz)  

    !! move to next time step
    t0 = t0+dt
    
  end do

end subroutine bomex_mapping

subroutine bomex_no_mapping(nz, cfl, snlevs, scfl, bomex_type)
!!
!! This subroutine run the BOMEX simulation with same mesh used for the 
!! dynamics and physics calculations. 
!!
!! INPUT:
!! nz: number of levels
!! cfl: CFL conditions used for the simulation
!! snlevs: the number levels (string)
!! scfl: CFL conditions used for simulation (string)
!! bomex_type: method used for the advections in the dynamics. The options are ``Linear", ``Cubic", and ``WENO"
!!


  use mod_bomex, only : bomex_init, set_bomex_forcing, bomex_ls_forcing
  use scm_physics, only : run_scm_physics
  use machine, only : r8 => real_kind
  
  
  implicit none
  
 
  integer, intent(in)           :: nz   !! nz= nlevs number of levels
  real(kind=r8), intent(in)     :: cfl  !! CFL condition used to calculate dt 
  character(len=16), intent(in) :: snlevs
  character(len=16), intent(in) :: scfl
  character(len=16), intent(in) :: bomex_type
  character(len=32)             :: fname
  real(kind=r8), parameter      :: ztop = 3000.0_r8
  real(kind=r8), dimension(nz)  :: u, v, th, prs, qv, qc, qr, rho, z
  real(kind=r8), dimension(nz+1):: ui, vi, thi, prsi, qvi, qci, qri, rhoi,  zi
  real(kind=r8), dimension(nz)  :: ug, vg, w, dthdt, dqvdt
  
  real(kind=r8) :: lat   ! Only relevant for test=1
  real(kind=r8) :: precl
  real(kind=r8) :: dt
  
  real(kind=r8) :: ps, dz, tf, t0, heat, evap, stress
  integer       :: i, k, kk, fid
  
  real(kind=r8)                   :: dt_write
  
  integer                         :: n_steps, start_map
  

  
  dt_write = 0.0_r8
  
  !! get filename 
  fname = trim("bomex_data/")//trim("bomex")//trim(bomex_type)//trim(snlevs)//trim(scfl)//trim(".dat")
  
  t0=0.0_r8         !! start time
  tf=6*3600.0_r8    !! end time 
  
  dt = cfl *3000.0_r8/(real(nz, kind=r8))
  lat = 0.0_r8
  
  n_steps = 0
  start_map = 3600*6

  !--- Evenly Spaced Grid ---!
  dz = ztop/real(nz, kind=r8)
  zi(1) = 0.0
  do k=2,nz+1
    zi(k) = zi(k-1) + dz          
    z(k-1) = 0.5*(zi(k) + zi(k-1))
  end do
  
  
  ! Set bomex initial state
  call bomex_init(u,v,rho,th,qv,qc,qr,prs,z,ps,nz)
  
  ! Only need prsi from this call
  call bomex_init(ui,vi,rhoi,thi,qvi,qci,qri,prsi,zi,ps,nz+1)
    
  call set_bomex_forcing(ug,vg,w,dthdt,dqvdt,heat,evap,stress,z,nz)
  

  do while(t0<tf)

    n_steps = n_steps+1

    !!** dynamics calculations **!!
    call bomex_ls_forcing(u,v,prs,rho,th,qv,qc,qr, &
                          ug,vg,w,dthdt,dqvdt,z,dt,nz, bomex_type)


    !!!! write to file every 15 min
    !!if(n_steps >= start_map)then
    !!if (dt_write >=  30.0 .or. dt_write .eq. tf .or. n_steps .eq. start_map)then
    if (t0 <= 18000.0 .and. t0 + dt > 18000.0 ) then
      fid=10
      open(unit=fid,file=fname,status='unknown')
      do k=1,nz
        write(fid,'(9(3x,E12.5))') z(k),th(k),u(k),v(k),qv(k),qc(k),qr(k)
      end do
      close(fid)
    endif
    !!  dt_write = 0.0
    !!endif
    !!dt_write = dt_write + dt
    !!endif


    !!** "physics" calulation.
    call run_scm_physics(u,v,th,qv,qc,qr,prs,prsi,rho,heat,evap,stress,z,zi,dt,nz)

    t0 = t0+dt
  end do

end subroutine bomex_no_mapping

subroutine f ( n, t, y, ydot )

!*****************************************************************************80
!
!! F evaluates the right hand sides of the ODE's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) t
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ydot(n)

  ydot(1) = y(2)
  ydot(2) = -y(1)

  return
end

subroutine advection1d()
!!
!! This example to solve the hyperbolic equation q_{t} + w q_{x} = 0 
!! use an upwind WENO method
!!
!! q^{t+1}_{i} = q^{t}_{i} - w/dx ( q^{t}_{i + 1/2} - q^{t}_{i-1/2} )
!!
!! Where the q^{t}_{i + 1/2} and  q^{t}_{i-1/2} are obtained using weno recontructions
!!
!!

  use mod_bomex

  implicit none

  integer, parameter    :: n = 201
  integer               :: k, fid
  real(kind=8)          :: a, b, dx, dt, t, t_old
  real(kind=8)          :: x(n), q_new(n)
  real(kind=8)          :: q_old(n), q1(n), q2(n)
  real(kind=8)          :: fluxm(n), fluxp(n)

  !! create th mesh 
  a = -1.0;		!! left bound
  b = 1.0;		!! right bound
  dx = (b-a) / real(n-1, kind=8)
  do k=1, n
    x(k) = a + dx * real(k-1, kind=8)
    q_new(k) = 1.0 / (1.0 + exp(-200.0*x(k)))
  enddo
  dt = 0.001;
  t=0;
  t_old=t

  fid=10
  open(unit=fid,file='advection1d.dat',status='unknown')
  do k=1,n
  write(fid,'(2(3x,E12.5))') x(k),q_new(k)
  write(6,'(2(3x,E12.5))') x(k),q_new(k)
  end do
  close(fid)
  write(6,'(A)') '---------------------------------------------'


  do while( t<0.5) 
    q_old = q_new; 
    q1 = q_old;
    q2 = q_old;
    call computeflux(q_old, fluxm, fluxp, n);   !! compute upwind flux
    do k=4,n-3
      q1(k) = q_old(k) + dt/dx * (fluxp(k+1)-fluxp(k));
      !q1(k) = q_old(k) - dt/dx * (fluxm(k)-fluxm(k-1));
      !write(*,*) 'k=',k, 'q1(k)=', q1(k)
    enddo
    call computeflux(q1, fluxm, fluxp, n);   !! compute upwind flux
    !print*, 'fluxm =', fluxm
    do k=4,n-3
      q2(k) = 1.0/4.0*(3.0*q_old(k) + q1(k) + dt/dx * (fluxp(k+1)-fluxp(k)));
      !q2(k) = 1.0/4.0*(3.0*q_old(k) + q1(k) - dt/dx * (fluxm(k)-fluxm(k-1)));
    enddo
    call computeflux(q2, fluxm, fluxp, n);   !!compute upwind flux
    do k=4,n-3
      q_new(k) = 1.0/3.0*(q_old(k) + 2.0*q2(k) + 2.0*dt/dx * (fluxp(k+1)-fluxp(k)));
      !q_new(k) = 1.0/3.0*(q_old(k) + 2.0*q2(k) - 2.0*dt/dx * (fluxm(k)-fluxm(k-1)));
    enddo
    if(abs(t-t_old-0.01) <= 1e-5) then
      write(*,*) 'Write to file t = ', t
      fid=10
      open(unit=fid,file='advection1d.dat',status='old', access='append')
      do k=1,n
      write(fid,'(2(3x,E12.5))') x(k),q_new(k)
      end do
      close(fid)
      t_old = t
    endif
    t = t + dt
    !!write(*,*) 't =', t
  enddo

 
end subroutine advection1d


subroutine lagrangePolyVal(x, y, n, xout, yout)
!!
!! contruct and evaluate lagrange polynomial at xout
!!
!! INPUT
!! x: input points to be used for the lagrange interpolation
!! y: data values associated to the input points in x
!! n: number of points
!! xout: point to interpolate to
!!
!! OUTPUT
!! yout: interpolated value
!!

  integer, intent(in)                   :: n
  real(kind=8), intent(in)              :: x(n), y(n), xout
  real(kind=8), intent(out)             :: yout
  integer                               :: i, j
  real(kind=8)                          :: tmp1, tmp2

  tmp2 = 0.0
  do i=1, n
    tmp1 = 1.0
    do j=1, n
      if(i /= j) then
       tmp1 = tmp1 * (xout - x(j)) / (x(i)-x(j))
      endif 
    enddo
    tmp2 = tmp2 + y(i) * tmp1
  enddo
  yout = tmp2

end subroutine



