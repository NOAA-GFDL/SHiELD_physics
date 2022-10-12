
!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Atmos Drivers project.
!*
!* This is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* It is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module atmos_model_mod
!-----------------------------------------------------------------------
!<OVERVIEW>
!  Driver for the atmospheric model, contains routines to advance the
!  atmospheric model state by one time step.
!</OVERVIEW>

!<DESCRIPTION>
!     This version of atmos_model_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the atmospheric model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Most atmospheric processes (dynamics,radiation,etc.) are performed
!     in the down routine. The up routine finishes the vertical diffusion
!     and computes moisture related terms (convection,large-scale condensation,
!     and precipitation).

!     The boundary variables needed by other component models for coupling
!     are contained in a derived data type. A variable of this derived type
!     is returned when initializing the atmospheric model. It is used by other
!     routines in this module and by coupling routines. The contents of
!     this derived type should only be modified by the atmospheric model.

!</DESCRIPTION>

use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
use mpp_mod,            only: mpp_clock_end, CLOCK_COMPONENT, MPP_CLOCK_SYNC
use mpp_mod,            only: mpp_min, mpp_max, mpp_error, mpp_chksum
use mpp_domains_mod,    only: domain2d
use mpp_mod,            only: mpp_get_current_pelist_name
#ifdef INTERNAL_FILE_NML
use mpp_mod,            only: input_nml_file
#else
use fms_mod,            only: open_namelist_file
#endif
use fms_mod,            only: file_exist, error_mesg
use fms_mod,            only: close_file, write_version_number, stdlog, stdout
use fms_mod,            only: clock_flag_default
use fms_mod,            only: check_nml_error
use diag_manager_mod,   only: diag_send_complete_instant
use time_manager_mod,   only: time_type, get_time, get_date, &
                              operator(+), operator(-)
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_number_tracers, get_tracer_names
use xgrid_mod,          only: grid_box_type
use atmosphere_mod,     only: atmosphere_init
use atmosphere_mod,     only: atmosphere_restart
use atmosphere_mod,     only: atmosphere_end
use atmosphere_mod,     only: atmosphere_state_update
use atmosphere_mod,     only: atmos_phys_driver_statein
use atmosphere_mod,     only: atmosphere_control_data
use atmosphere_mod,     only: atmosphere_resolution, atmosphere_domain
use atmosphere_mod,     only: atmosphere_grid_bdry, atmosphere_grid_ctr
use atmosphere_mod,     only: atmosphere_dynamics, atmosphere_diag_axes
use atmosphere_mod,     only: atmosphere_etalvls, atmosphere_hgt
!rab use atmosphere_mod,     only: atmosphere_tracer_postinit
use atmosphere_mod,     only: atmosphere_diss_est, atmosphere_nggps_diag
use atmosphere_mod,     only: atmosphere_scalar_field_halo
use atmosphere_mod,     only: set_atmosphere_pelist
use atmosphere_mod,     only: atmosphere_coarse_graining_parameters
use atmosphere_mod,     only: atmosphere_coarse_diag_axes
use atmosphere_mod,     only: atmosphere_coarsening_strategy
use atmosphere_mod,     only: Atm, mygrid
use block_control_mod,  only: block_control_type, define_blocks_packed
use IPD_typedefs,       only: IPD_init_type, IPD_control_type, &
                              IPD_data_type, IPD_diag_type,    &
                              IPD_restart_type, kind_phys
use IPD_driver,         only: IPD_initialize, IPD_setup_step, &
                              IPD_radiation_step,             &
                              IPD_physics_step1,              &
                              IPD_physics_step2, IPD_physics_end
#ifdef STOCHY 
use stochastic_physics, only: init_stochastic_physics,         &
                              run_stochastic_physics
use stochastic_physics_sfc, only: run_stochastic_physics_sfc
#endif
use FV3GFS_io_mod,      only: FV3GFS_restart_read, FV3GFS_restart_write, &
                              FV3GFS_IPD_checksum,                       &
                              gfdl_diag_register, gfdl_diag_output, &
                              FV3GFS_restart_write_coarse, FV3GFS_diag_register_coarse, &
                              sfc_data_override
use FV3GFS_io_mod,      only: register_diag_manager_controlled_diagnostics, register_coarse_diag_manager_controlled_diagnostics
use FV3GFS_io_mod,      only: send_diag_manager_controlled_diagnostic_data
use fv_iau_mod,         only: iau_external_data_type,getiauforcing,iau_initialize
use module_ocean,       only: ocean_init
!-----------------------------------------------------------------------

implicit none
private

public update_atmos_radiation_physics
public update_atmos_model_state
public update_atmos_model_dynamics
public atmos_model_init, atmos_model_end, atmos_data_type
public atmos_model_restart
!-----------------------------------------------------------------------

!<PUBLICTYPE >
 type atmos_data_type
     type (domain2d)               :: domain             ! domain decomposition
     type (domain2d)               :: domain_for_read    ! domain decomposition for reads
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     real,                 pointer, dimension(:,:) :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
     real,                 pointer, dimension(:,:) :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
     real(kind=kind_phys), pointer, dimension(:,:) :: lon      => null() ! local longitude axis grid box centers in radians.
     real(kind=kind_phys), pointer, dimension(:,:) :: lat      => null() ! local latitude axis grid box centers in radians.
     type (time_type)              :: Time               ! current time
     type (time_type)              :: Time_step          ! atmospheric time step.
     type (time_type)              :: Time_init          ! reference time.
     integer                       :: iau_offset         ! iau running window length
     integer, pointer              :: pelist(:) =>null() ! pelist where atmosphere is running.
     logical                       :: pe                 ! current pe.
     type(grid_box_type)           :: grid               ! hold grid information needed for 2nd order conservative flux exchange
                                                         ! to calculate gradient on cubic sphere grid.
     integer                       :: layout(2)          ! computer task laytout
     logical                       :: regional           ! true if domain is regional
     logical                       :: bounded_domain     ! true if domain is bounded
     real(kind=8), pointer, dimension(:) :: ak
     real(kind=8), pointer, dimension(:) :: bk
     real(kind=8), pointer, dimension(:,:,:) :: layer_hgt
     real(kind=8), pointer, dimension(:,:,:) :: level_hgt
     real(kind=kind_phys), pointer, dimension(:,:) :: dx
     real(kind=kind_phys), pointer, dimension(:,:) :: dy
     real(kind=8), pointer, dimension(:,:) :: area
     type(domain2d)                :: coarse_domain      ! domain decomposition of the coarse grid
     logical                       :: write_coarse_restart_files  ! whether to write coarse restart files
     logical                       :: write_only_coarse_intermediate_restarts  ! whether to write only coarse intermediate restart files
     character(len=64)             :: coarsening_strategy  ! Strategy for coarse-graining diagnostics and restart files
end type atmos_data_type
!</PUBLICTYPE >

integer :: fv3Clock, getClock, overrideClock, updClock, setupClock, radClock, physClock

!-----------------------------------------------------------------------
integer :: blocksize       = 1
logical :: chksum_debug    = .false.
logical :: dycore_only     = .false.
logical :: debug           = .false.
logical :: sync            = .false.
logical :: first_time_step = .false.
logical :: fprint          = .true.
logical :: enforce_rst_cksum = .true. ! enforce or override data integrity restart checksums
real, dimension(4096) :: fdiag = 0. ! xic: TODO: this is hard coded, space can run out in some cases. Should make it allocatable.
logical :: fdiag_override = .false. ! lmh: if true overrides fdiag and fhzer: all quantities are zeroed out
                                    ! after every calcluation, output interval and accumulation/avg/max/min
                                    ! are controlled by diag_manager, fdiag controls output interval only
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, first_time_step, fdiag, fprint, &
                           fdiag_override, enforce_rst_cksum
type (time_type) :: diag_time, diag_time_fhzero
logical :: fdiag_fix = .false.

!--- concurrent and decoupled radiation and physics variables
!----------------
!  IPD containers
!----------------
type(IPD_control_type)              :: IPD_Control
type(IPD_data_type),    allocatable :: IPD_Data(:)  ! number of blocks
type(IPD_diag_type)                 :: IPD_Diag(250)
type(IPD_restart_type)              :: IPD_Restart

!--------------
! IAU container
!--------------
type(iau_external_data_type)        :: IAU_Data

!-----------------
!  Block container
!-----------------
type (block_control_type), target   :: Atm_block

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

real(kind=kind_phys), parameter :: zero = 0.0_kind_phys

contains

!#######################################################################
! <SUBROUTINE NAME="update_radiation_physics">
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down".
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_radiation_physics (Atmos)
!   </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

subroutine update_atmos_radiation_physics (Atmos)
!-----------------------------------------------------------------------
  type (atmos_data_type), intent(in) :: Atmos
!--- local variables---
    integer :: nb, jdat(8)
    integer :: nthrds

#ifdef OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif

    if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "statein driver"
!--- get atmospheric state from the dynamic core
    call set_atmosphere_pelist()
    call mpp_clock_begin(getClock)
    if (IPD_control%do_skeb) call atmosphere_diss_est (IPD_control%skeb_npass) !  do smoothing for SKEB
    call atmos_phys_driver_statein (IPD_data, Atm_block)
    call mpp_clock_end(getClock)

!--- get varied surface data
    call mpp_clock_begin(overrideClock)
    call sfc_data_override (Atmos%Time, IPD_data, Atm_block, IPD_Control)
    call mpp_clock_end(overrideClock)

!--- if dycore only run, set up the dummy physics output state as the input state
    if (dycore_only) then
      do nb = 1,Atm_block%nblks
        IPD_Data(nb)%Stateout%gu0 = IPD_Data(nb)%Statein%ugrs
        IPD_Data(nb)%Stateout%gv0 = IPD_Data(nb)%Statein%vgrs
        IPD_Data(nb)%Stateout%gt0 = IPD_Data(nb)%Statein%tgrs
        IPD_Data(nb)%Stateout%gq0 = IPD_Data(nb)%Statein%qgrs
      enddo
    else
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "setup step"

!--- update IPD_Control%jdat(8)
      jdat(:) = 0
      call get_date (Atmos%Time, jdat(1), jdat(2), jdat(3),  &
                                 jdat(5), jdat(6), jdat(7))
      IPD_Control%jdat(:) = jdat(:)
!--- execute the IPD atmospheric setup step
      call mpp_clock_begin(setupClock)
      call IPD_setup_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)

#ifdef STOCHY
!--- call stochastic physics pattern generation / cellular automata
      if (IPD_Control%do_sppt .OR. IPD_Control%do_shum .OR. IPD_Control%do_skeb .OR. IPD_Control%do_sfcperts) then
         call run_stochastic_physics(IPD_Control, IPD_Data(:)%Grid, IPD_Data(:)%Coupling, nthrds)
      end if
#endif

      call mpp_clock_end(setupClock)

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "radiation driver"
!--- execute the IPD atmospheric radiation subcomponent (RRTM)
      call mpp_clock_begin(radClock)
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_radiation_step (IPD_Control, IPD_Data(nb), IPD_Diag, IPD_Restart)
      enddo
      call mpp_clock_end(radClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'RADIATION STEP  ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "physics driver"
!--- execute the IPD atmospheric physics step1 subcomponent (main physics driver)
      call mpp_clock_begin(physClock)
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_physics_step1 (IPD_Control, IPD_Data(nb), IPD_Diag, IPD_Restart)
      enddo
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP1   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "stochastic physics driver"
!--- execute the IPD atmospheric physics step2 subcomponent (stochastic physics driver)
      call mpp_clock_begin(physClock)
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_physics_step2 (IPD_Control, IPD_Data(nb), IPD_Diag, IPD_Restart)
      enddo
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP2   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif
      call getiauforcing(IPD_Control,IAU_data)
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "end of radiation and physics step"
    endif

!-----------------------------------------------------------------------
 end subroutine update_atmos_radiation_physics
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step, iau_offset)

#ifdef OPENMP
  use omp_lib
#endif
  use mpp_mod, only: mpp_npes

  type (atmos_data_type), intent(inout) :: Atmos
  type (time_type), intent(in) :: Time_init, Time, Time_step
  integer, intent(in) :: iau_offset
!--- local variables ---
  integer :: unit, ntdiag, ntfamily, i, j, k
  integer :: mlon, mlat, nlon, nlat, nlev, sec, dt, sec_prev
  integer :: ierr, io, logunit
  integer :: idx, tile_num
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: blk, ibs, ibe, jbs, jbe
  real(kind=kind_phys) :: dt_phys
  real, allocatable :: q(:,:,:,:), p_half(:,:,:)
  character(len=80) :: control
  character(len=64) :: filename, filename2, pelist_name
  character(len=132) :: text
  logical :: p_hydro, hydro, fexist
  logical, save :: block_message = .true.
  type(IPD_init_type) :: Init_parm
  integer :: bdat(8), cdat(8)
  integer :: ntracers
  integer :: kdt_prev
  character(len=32), allocatable, target :: tracer_names(:)
  integer :: coarse_diagnostic_axes(4)
  integer :: nthrds
  !-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step
   Atmos % iau_offset = iau_offset
   call get_time (Atmos % Time_step, sec)
   call get_time (Atmos%Time - Atmos%Time_init, sec_prev)
   dt_phys = real(sec)      ! integer seconds
   kdt_prev = int(sec_prev / dt_phys)

   logunit = stdlog()

!-----------------------------------------------------------------------
! initialize atmospheric model -----

!---------- initialize atmospheric dynamics -------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area, IAU_Data)

   IF ( file_exist('input.nml')) THEN
#ifdef INTERNAL_FILE_NML
      read(input_nml_file, nml=atmos_model_nml, iostat=io)
      ierr = check_nml_error(io, 'atmos_model_nml')
#else
      unit = open_namelist_file ( )
      ierr=1
      do while (ierr /= 0)
         read  (unit, nml=atmos_model_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'atmos_model_nml')
      enddo
 10     call close_file (unit)
#endif
   endif
!-----------------------------------------------------------------------
   call atmosphere_resolution (nlon, nlat, global=.false.)
   call atmosphere_resolution (mlon, mlat, global=.true.)
   call alloc_atmos_data_type (nlon, nlat, Atmos)
   call atmosphere_domain (Atmos%domain, Atmos%domain_for_read, Atmos%layout, Atmos%regional, &
                           Atmos%bounded_domain)
   call atmosphere_diag_axes (Atmos%axes)
   call atmosphere_etalvls (Atmos%ak, Atmos%bk, flip=.true.)
   call atmosphere_grid_bdry (Atmos%lon_bnd, Atmos%lat_bnd, global=.false.)
   call atmosphere_grid_ctr (Atmos%lon, Atmos%lat)
   call atmosphere_hgt (Atmos%layer_hgt, 'layer', relative=.false., flip=.true.)
   call atmosphere_hgt (Atmos%level_hgt, 'level', relative=.false., flip=.true.)
   call atmosphere_coarse_graining_parameters(Atmos%coarse_domain, Atmos%write_coarse_restart_files, &
                                              Atmos%write_only_coarse_intermediate_restarts)
   call atmosphere_coarsening_strategy(Atmos%coarsening_strategy)

!-----------------------------------------------------------------------
!--- before going any further check definitions for 'blocks'
!-----------------------------------------------------------------------
   call atmosphere_control_data (isc, iec, jsc, jec, nlev, p_hydro, hydro, tile_num)
   call define_blocks_packed ('atmos_model', Atm_block, isc, iec, jsc, jec, nlev, &
                              blocksize, block_message)

   allocate(IPD_Data(Atm_block%nblks))

#ifdef OPENMP
   nthrds = omp_get_max_threads()
#else
   nthrds = 1
#endif

!--- update IPD_Control%jdat(8)
   bdat(:) = 0
   call get_date (Time_init, bdat(1), bdat(2), bdat(3),  &
                             bdat(5), bdat(6), bdat(7))
   cdat(:) = 0
   call get_date (Time,      cdat(1), cdat(2), cdat(3),  &
                             cdat(5), cdat(6), cdat(7))
   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
   allocate (tracer_names(ntracers))
   do i = 1, ntracers
     call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
   enddo
!--- setup IPD Init_parm
   Init_parm%me              =  mpp_pe()
   Init_parm%master          =  mpp_root_pe()
   Init_parm%tile_num        =  tile_num
   Init_parm%isc             =  isc
   Init_parm%jsc             =  jsc
   Init_parm%nx              =  nlon
   Init_parm%ny              =  nlat
   Init_parm%levs            =  nlev
   Init_parm%cnx             =  mlon
   Init_parm%cny             =  mlat
   Init_parm%gnx             =  Init_parm%cnx*4
   Init_parm%gny             =  Init_parm%cny*2
   Init_parm%nlunit          =  9999
   Init_parm%logunit         =  logunit
   Init_parm%bdat(:)         =  bdat(:)
   Init_parm%cdat(:)         =  cdat(:)
   Init_parm%dt_dycore       =  dt_phys
   Init_parm%dt_phys         =  dt_phys
   Init_parm%iau_offset      =  Atmos%iau_offset
   Init_parm%blksz           => Atm_block%blksz
   Init_parm%ak              => Atmos%ak
   Init_parm%bk              => Atmos%bk
   Init_parm%xlon            => Atmos%lon
   Init_parm%xlat            => Atmos%lat
   Init_parm%area            => Atmos%area
   Init_parm%tracer_names    => tracer_names

#ifdef INTERNAL_FILE_NML
   allocate(Init_parm%input_nml_file, mold=input_nml_file)
   Init_parm%input_nml_file  => input_nml_file
   Init_parm%fn_nml='using internal file'
#else
   pelist_name=mpp_get_current_pelist_name()
   Init_parm%fn_nml='input_'//trim(pelist_name)//'.nml'
   inquire(FILE=Init_parm%fn_nml, EXIST=fexist)
   if (.not. fexist ) then
      Init_parm%fn_nml='input.nml'
   endif
#endif

   call IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Init_parm)

#ifdef STOCHY
   if (IPD_Control%do_sppt .OR. IPD_Control%do_shum .OR. IPD_Control%do_skeb .OR. IPD_Control%do_sfcperts) then
      ! Initialize stochastic physics
      call init_stochastic_physics(IPD_Control, Init_parm, mpp_npes(), nthrds)
      if (mpp_pe() == mpp_root_pe()) print *,'do_skeb=',IPD_Control%do_skeb
   end if

   if (IPD_Control%do_sfcperts) then
      ! Get land surface perturbations here (move to GFS_time_vary
      ! step if wanting to update each time-step)
      call run_stochastic_physics_sfc(IPD_Control, IPD_Data(:)%Grid, IPD_Data(:)%Coupling)
   end if
#endif

   Atm(mygrid)%flagstruct%do_diss_est = IPD_Control%do_skeb

!  initialize the IAU module
   call iau_initialize (IPD_Control,IAU_data,Init_parm)

   IPD_Control%kdt_prev = kdt_prev

!--- initialize slab ocean model or mixed layer ocean model
#ifdef INTERNAL_FILE_NML
   if (IPD_Control%do_ocean) call ocean_init (IPD_Control, Init_parm%logunit, input_nml_file)
#else
   if (IPD_Control%do_ocean) call ocean_init (IPD_Control, Init_parm%logunit)
#endif

   Init_parm%blksz           => null()
   Init_parm%ak              => null()
   Init_parm%bk              => null()
   Init_parm%xlon            => null()
   Init_parm%xlat            => null()
   Init_parm%area            => null()
   Init_parm%tracer_names    => null()
   deallocate (tracer_names)

   !--- update tracers in FV3 with any initialized during the physics/radiation init phase
!rab   call atmosphere_tracer_postinit (IPD_Data, Atm_block)

   call atmosphere_nggps_diag (Time, init=.true.)
   call gfdl_diag_register (Time, IPD_Data(:)%Sfcprop, IPD_Data(:)%IntDiag, IPD_Data%Cldprop, &
        Atm_block, Atmos%axes, IPD_Control%nfxr, IPD_Control%ldiag3d, &
        IPD_Control%nkld, IPD_Control%levs)
   call register_diag_manager_controlled_diagnostics(Time, IPD_Data(:)%IntDiag, Atm_block%nblks, Atmos%axes)
   if (Atm(mygrid)%coarse_graining%write_coarse_diagnostics) then
       call atmosphere_coarse_diag_axes(coarse_diagnostic_axes)
       call FV3GFS_diag_register_coarse(Time, coarse_diagnostic_axes)
       call register_coarse_diag_manager_controlled_diagnostics(Time, coarse_diagnostic_axes)
    endif
   if (.not. dycore_only) &
      call FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, IPD_Control, Atmos%domain_for_read, enforce_rst_cksum)
      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'RESTART READ  ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

   !--- set the initial diagnostic timestamp
   diag_time = Time
   if (Atmos%iau_offset > zero) then
     call get_time (Atmos%Time - Atmos%Time_init, sec)
     if (sec < Atmos%iau_offset*3600) then
       diag_time = Atmos%Time_init
       diag_time_fhzero = Atmos%Time
     endif
   endif

   !---- print version number to logfile ----

   call write_version_number ( version, tagname )
   !--- write the namelist to a log file
   if (mpp_pe() == mpp_root_pe()) then
      unit = stdlog( )
      write (unit, nml=atmos_model_nml)
      call close_file (unit)
   endif

   !--- get fdiag
#ifdef GFS_PHYS
!--- check fdiag to see if it is an interval or a list
   if (fdiag_override) then
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---OVERRIDING fdiag: USING SETTINGS IN diag_table for GFS PHYSICS DIAGS"
      IPD_Control%fhzero = dt_phys / 3600.
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---fhzero IS SET TO dt_atmos: ALL DIAGNOSTICS ARE SINGLE-STEP"
   else
      if (nint(fdiag(2)) == 0) then
         fdiag_fix = .true.
         do i = 2, size(fdiag,1)
            fdiag(i) = fdiag(1) * i
         enddo
      endif
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---fdiag",fdiag(1:40)
   endif
#endif

   setupClock = mpp_clock_id( 'GFS Step Setup        ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   overrideClock = mpp_clock_id( 'GFS Override          ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   radClock   = mpp_clock_id( 'GFS Radiation         ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   physClock  = mpp_clock_id( 'GFS Physics           ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   getClock   = mpp_clock_id( 'Dynamics get state    ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   updClock   = mpp_clock_id( 'Dynamics update state ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   if (sync) then
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default+MPP_CLOCK_SYNC, grain=CLOCK_COMPONENT )
   else
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   endif

!-----------------------------------------------------------------------
end subroutine atmos_model_init
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_dynamics"
!
! <OVERVIEW>
subroutine update_atmos_model_dynamics (Atmos)
! run the atmospheric dynamics to advect the properties
  type (atmos_data_type), intent(in) :: Atmos

    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call atmosphere_dynamics (Atmos%Time)
    call mpp_clock_end(fv3Clock)

end subroutine update_atmos_model_dynamics
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_state"
!
! <OVERVIEW>
subroutine update_atmos_model_state (Atmos)
! to update the model state after all concurrency is completed
  type (atmos_data_type), intent(inout) :: Atmos
!--- local variables
  integer :: isec,seconds,isec_fhzero
  real(kind=kind_phys) :: time_int, time_intfull
  integer :: is, ie, js, je, kt

    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call mpp_clock_begin(updClock)
    call atmosphere_state_update (Atmos%Time, IPD_Data, IAU_Data, Atm_block)
    call mpp_clock_end(updClock)
    call mpp_clock_end(fv3Clock)

    if (chksum_debug) then
      if (mpp_pe() == mpp_root_pe()) print *,'UPDATE STATE    ', IPD_Control%kdt, IPD_Control%fhour
      call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
    endif

!------ advance time ------
    Atmos % Time = Atmos % Time + Atmos % Time_step

    call atmosphere_control_data(is, ie, js, je, kt)
    call send_diag_manager_controlled_diagnostic_data(Atmos%Time, &
       Atm_block, IPD_Data, IPD_Control%nx, IPD_Control%ny, IPD_Control%levs, &
       Atm(mygrid)%coarse_graining%write_coarse_diagnostics, &
       real(Atm(mygrid)%delp(is:ie,js:je,:), kind=kind_phys), &
       Atmos%coarsening_strategy, real(Atm(mygrid)%ptop, kind=kind_phys))

    call get_time (Atmos%Time - diag_time, isec)
    call get_time (Atmos%Time - Atmos%Time_init, seconds)

    time_int = real(isec)
    if (ANY(nint(fdiag(:)*3600.0) == seconds) .or. (fdiag_fix .and. mod(seconds, nint(fdiag(1)*3600.0)) .eq. 0) .or. (IPD_Control%kdt == 1 .and. first_time_step) ) then
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---isec,seconds",isec,seconds
      if (mpp_pe() == mpp_root_pe()) write(6,*) ' gfs diags time since last bucket empty: ',time_int/3600.,'hrs'
      call atmosphere_nggps_diag(Atmos%Time)
    endif
    if (ANY(nint(fdiag(:)*3600.0) == seconds) .or. (fdiag_fix .and. mod(seconds, nint(fdiag(1)*3600.0)) .eq. 0) .or. (IPD_Control%kdt == 1 .and. first_time_step)) then
      if(Atmos%iau_offset > zero) then
        if( time_int - Atmos%iau_offset*3600. > zero ) then
          time_int = time_int - Atmos%iau_offset*3600.
        else if(seconds == Atmos%iau_offset*3600) then
          call get_time (Atmos%Time - diag_time_fhzero, isec_fhzero)
          time_int = real(isec_fhzero)
          if (mpp_pe() == mpp_root_pe()) write(6,*) "---iseczero",isec_fhzero
        endif
      endif
      time_intfull = real(seconds)
      if(Atmos%iau_offset > zero) then
        if( time_intfull - Atmos%iau_offset*3600. > zero) then
          time_intfull = time_intfull - Atmos%iau_offset*3600.
        endif
      endif
      call gfdl_diag_output(Atmos%Time, Atm_block, IPD_Data, IPD_Control%nx, IPD_Control%ny, fprint, &
                            IPD_Control%levs, 1, 1, 1.d0, time_int, time_intfull, &
                            IPD_Control%fhswr, IPD_Control%fhlwr, &
                            mod(seconds, nint(fdiag(1)*3600.0)) .eq. 0, &
                            Atm(mygrid)%coarse_graining%write_coarse_diagnostics,&
                            real(Atm(mygrid)%delp(is:ie,js:je,:), kind=kind_phys), &
                            Atmos%coarsening_strategy, real(Atm(mygrid)%ptop, kind=kind_phys))
      call diag_send_complete_instant (Atmos%Time)
      if (mod(isec,nint(3600*IPD_Control%fhzero)) == 0) diag_time = Atmos%Time
    endif

 end subroutine update_atmos_model_state
! </SUBROUTINE>



!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_end (Atmos)
  type (atmos_data_type), intent(inout) :: Atmos
!---local variables
  integer :: idx

    call IPD_physics_end (IPD_Control)

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----

    call atmosphere_end (Atmos % Time, Atmos%grid)
    if (.not. dycore_only) then
       call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
            IPD_Control, Atmos%domain)
       if (Atmos%write_coarse_restart_files) then
          call FV3GFS_restart_write_coarse(IPD_Data, IPD_Restart, Atm_block, &
            IPD_Control, Atmos%coarse_domain)
       endif
    endif

end subroutine atmos_model_end

! </SUBROUTINE>
!#######################################################################
! <SUBROUTINE NAME="atmos_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine atmos_model_restart(Atmos, timestamp)
  type (atmos_data_type),   intent(inout) :: Atmos
  character(len=*),  intent(in)           :: timestamp

    call atmosphere_restart(timestamp)
    if (.not. dycore_only) then
       if (.not. Atmos%write_only_coarse_intermediate_restarts) then
          call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
               IPD_Control, Atmos%domain, timestamp)
       endif
       if (Atmos%write_coarse_restart_files) then
          call FV3GFS_restart_write_coarse(IPD_Data, IPD_Restart, Atm_block, &
               IPD_Control, Atmos%coarse_domain, timestamp)
       endif
    endif
end subroutine atmos_model_restart
! </SUBROUTINE>

!#######################################################################
!#######################################################################
! <SUBROUTINE NAME="atmos_data_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the atmos_data_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the atmos_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, atm)
! </TEMPLATE>

! <IN NAME="Atm" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields in the atmos_data_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine atmos_data_type_chksum(id, timestep, atm)
type(atmos_data_type), intent(in) :: atm
    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    integer :: n, outunit

100 format("CHECKSUM::",A32," = ",Z20)
101 format("CHECKSUM::",A16,a,'%',a," = ",Z20)

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(Atmos_data_type):: ', id, timestep
  write(outunit,100) ' atm%lon_bnd                ', mpp_chksum(atm%lon_bnd               )
  write(outunit,100) ' atm%lat_bnd                ', mpp_chksum(atm%lat_bnd               )
  write(outunit,100) ' atm%lon                    ', mpp_chksum(atm%lon                   )
  write(outunit,100) ' atm%lat                    ', mpp_chksum(atm%lat                   )

end subroutine atmos_data_type_chksum

! </SUBROUTINE>

  subroutine alloc_atmos_data_type (nlon, nlat, Atmos)
   integer, intent(in) :: nlon, nlat
   type(atmos_data_type), intent(inout) :: Atmos
    allocate ( Atmos % lon_bnd  (nlon+1,nlat+1), &
               Atmos % lat_bnd  (nlon+1,nlat+1), &
               Atmos % lon      (nlon,nlat),     &
               Atmos % lat      (nlon,nlat)      )

  end subroutine alloc_atmos_data_type

  subroutine dealloc_atmos_data_type (Atmos)
   type(atmos_data_type), intent(inout) :: Atmos
    deallocate (Atmos%lon_bnd, &
                Atmos%lat_bnd, &
                Atmos%lon,     &
                Atmos%lat      )
  end subroutine dealloc_atmos_data_type

end module atmos_model_mod
