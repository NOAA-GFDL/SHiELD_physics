module FV3GFS_io_mod

!-----------------------------------------------------------------------
!    gfs_physics_driver_mod defines the GFS physics routines used by
!    the GFDL FMS system to obtain tendencies and boundary fluxes due
!    to the physical parameterizations and processes that drive
!    atmospheric time tendencies for use by other components, namely
!    the atmospheric dynamical core.
!
!    NOTE: This module currently supports only the operational GFS
!          parameterizations as of September 2015.  Further development
!          is needed to support the full suite of physical
!          parameterizations present in the GFS physics package.
!-----------------------------------------------------------------------
!
!--- FMS/GFDL modules
  use block_control_mod,  only: block_control_type
  use mpp_mod,            only: mpp_error, mpp_pe, mpp_root_pe, &
                                mpp_chksum, NOTE, FATAL, mpp_get_current_pelist_name
  use mpp_domains_mod,    only: domain2d, mpp_get_compute_domain
  use fms_mod,            only: stdout
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, unlimited,      &
                                open_file, close_file, register_field, &
                                register_axis, register_restart_field, &
                                register_variable_attribute,           &
                                read_restart, write_restart,           &
                                get_global_io_domain_indices,          &
                                dimension_exists, write_data
  use time_manager_mod,   only: time_type
  use data_override_mod,  only: data_override
  use diag_manager_mod,   only: register_diag_field, send_data
  use fv_mp_mod,          only: is_master, mp_reduce_sum, mp_reduce_min, mp_reduce_max
!
!--- GFS physics modules
  use machine,            only: kind_phys
!--- variables needed for calculating 'sncovr'
  use namelist_soilveg,   only: salp_data, snupx

!
! --- variables needed for Noah MP init
!
  use noahmp_tables,      only: laim_table,saim_table,sla_table,      &
                                bexp_table,smcmax_table,smcwlt_table, &
                                dwsat_table,dksat_table,psisat_table, &
                                isurban_table,isbarren_table,         &
                                isice_table,iswater_table

!
!--- GFS_typedefs
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_diag_type, GFS_grid_type
  use GFS_typedefs,       only: GFS_cldprop_type
  use ozne_def,           only: oz_coeff
!
!--- IPD typdefs
  use IPD_typedefs,       only: IPD_control_type, IPD_data_type, &
                                IPD_restart_type
!--- GFS physics constants
  use physcons,           only: pi => con_pi, RADIUS => con_rerth, rd => con_rd
!--- needed for dq3dt output
  use ozne_def,           only: oz_coeff
!--- needed for cold-start capability to initialize q2m
  use gfdl_cld_mp_mod,    only: wqs, qs_init
  use coarse_graining_mod, only: block_mode, block_upsample, block_min, block_max, block_sum, weighted_block_average
  use coarse_graining_mod, only: MODEL_LEVEL, PRESSURE_LEVEL
  use coarse_graining_mod, only: vertical_remapping_requirements, get_coarse_array_bounds
  use coarse_graining_mod, only: vertically_remap_field, mask_area_weights
!
!-----------------------------------------------------------------------
  implicit none
  private

  !--- public interfaces ---
  public  FV3GFS_restart_read, FV3GFS_restart_write, FV3GFS_restart_write_coarse
  public  FV3GFS_IPD_checksum
  public  gfdl_diag_register, gfdl_diag_output
  public  FV3GFS_diag_register_coarse, register_diag_manager_controlled_diagnostics
  public  register_coarse_diag_manager_controlled_diagnostics
  public  send_diag_manager_controlled_diagnostic_data
  public  sfc_data_override
  public  gfdl_diag_type, Diag, Diag_diag_manager_controlled

  !--- GFDL filenames
  character(len=32)  :: fn_oro = 'oro_data.nc'
  character(len=32)  :: fn_srf = 'sfc_data.nc'
  character(len=32)  :: fn_phy = 'phy_data.nc'
  character(len=32)  :: fn_ifsSST = 'ifs_sst_data.nc'

  !--- GFDL FMS netcdf restart data types
  type(FmsNetcdfDomainFile_t) :: Oro_restart
  type(FmsNetcdfDomainFile_t) :: Sfc_restart
  type(FmsNetcdfDomainFile_t) :: Phy_restart
  type(FmsNetcdfDomainFile_t) :: ifsSST_restart

  !--- GFDL FMS restart containers
  character(len=32),    allocatable,         dimension(:)       :: oro_name2, sfc_name2, sfc_name3
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: oro_var2, sfc_var2, phy_var2
  real(kind=kind_phys), allocatable, target, dimension(:,:)     :: ifsSST
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3, phy_var3
  !--- Noah MP restart containers
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3sn,sfc_var3eq,sfc_var3zn

  ! Coarse graining
  real(kind=kind_phys), allocatable, target, dimension(:,:,:) :: sfc_var2_coarse
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3_coarse
  type(FmsNetcdfDomainFile_t) :: Sfc_restart_coarse

  integer :: isco, ieco, jsco, jeco, levo

!-RAB
  type data_subtype
    real(kind=kind_phys), dimension(:),   pointer :: var2 => NULL()
    real(kind=kind_phys), dimension(:,:), pointer :: var3 => NULL()
    real(kind=kind_phys), dimension(:),   pointer :: var21 => NULL()
  end type data_subtype
  !--- data type definition for use with GFDL FMS diagnostic manager until write component is working
  type gfdl_diag_type
    integer :: id = -1
    integer :: axes = -1
    logical :: time_avg = .false.
    character(len=64)    :: time_avg_kind = ''
    character(len=64)    :: mod_name = ''
    character(len=128)   :: name = ''
    character(len=128)   :: desc = ''
    character(len=64)    :: unit = ''
    character(len=64)    :: mask = ''
    character(len=64)    :: intpl_method = ''
    real(kind=kind_phys) :: cnvfac
    type(data_subtype), dimension(:), allocatable :: data

    ! Add an attribute that specifies the coarse-graining method for the
    ! variable.  By default we will set this as unspecified and raise an error
    ! if a user asks to coarse-grain a variable that does not have a supported
    ! method for coarse-graining.  Currently supported methods are:
    !
    ! 'area_weighted'
    ! 'mass_weighted'
    !
    ! In the future we may want to support more methods, e.g. for the land
    ! surface variables, which may require masking.
    character(len=64) :: coarse_graining_method = 'unspecified'
!rab    real(kind=kind_phys), dimension(:),   pointer :: var2 => NULL()
!rab    real(kind=kind_phys), dimension(:),   pointer :: var21 => NULL()
   end type gfdl_diag_type
   real(kind=kind_phys) :: zhour
!
   integer :: tot_diag_idx = 0
   integer, parameter :: DIAG_SIZE = 500
   real(kind=kind_phys), parameter :: missing_value = 9.99e20
   type(gfdl_diag_type), dimension(DIAG_SIZE) :: Diag, Diag_coarse, Diag_diag_manager_controlled, Diag_diag_manager_controlled_coarse
!-RAB


!--- miscellaneous other variables
  logical :: module_is_initialized = .FALSE.

  character(len=64) :: AREA_WEIGHTED = 'area_weighted'
  character(len=64) :: MASKED_AREA_WEIGHTED = 'masked_area_weighted'
  character(len=64) :: MASS_WEIGHTED = 'mass_weighted'
  character(len=64) :: MODE = 'mode'

  CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!--------------------
! FV3GFS_restart_read
!--------------------
  subroutine FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, Model, fv_domain, enforce_rst_cksum)
    type(IPD_data_type),      intent(inout) :: IPD_Data(:)
    type(IPD_restart_type),   intent(inout) :: IPD_Restart
    type(block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),   intent(inout) :: Model
    type(domain2d),           intent(in)    :: fv_domain
    logical,                  intent(in)    :: enforce_rst_cksum

    !--- read in surface data from chgres
    call sfc_prop_restart_read (IPD_Data%Sfcprop, Atm_block, Model, fv_domain, enforce_rst_cksum)

    !--- read in
    if (Model%sfc_override) call sfc_prop_override  (IPD_Data%Sfcprop, IPD_Data%Grid, Atm_block, Model, fv_domain)

    !--- read in physics restart data
    call phys_restart_read (IPD_Restart, Atm_block, Model, fv_domain, enforce_rst_cksum)

  end subroutine FV3GFS_restart_read

!---------------------
! FV3GFS_restart_write
!---------------------
  subroutine FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, Model, fv_domain, timestamp)
    type(IPD_data_type),         intent(inout) :: IPD_Data(:)
    type(IPD_restart_type),      intent(inout) :: IPD_Restart
    type(block_control_type),    intent(in)    :: Atm_block
    type(IPD_control_type),      intent(in)    :: Model
    type(domain2d),              intent(in)    :: fv_domain
    character(len=32), optional, intent(in)    :: timestamp

    !--- read in surface data from chgres
    call sfc_prop_restart_write (IPD_Data%Sfcprop, Atm_block, Model, fv_domain, timestamp)

    !--- read in physics restart data
    call phys_restart_write (IPD_Restart, Atm_block, Model, fv_domain, timestamp)

  end subroutine FV3GFS_restart_write

  subroutine FV3GFS_restart_write_coarse (IPD_Data, IPD_Restart, Atm_block, Model, coarse_domain, timestamp)
    type(IPD_data_type),         intent(inout) :: IPD_Data(:)
    type(IPD_restart_type),      intent(inout) :: IPD_Restart
    type(block_control_type),    intent(in)    :: Atm_block
    type(IPD_control_type),      intent(in)    :: Model
    type(domain2d),              intent(in)    :: coarse_domain
    character(len=32), optional, intent(in)    :: timestamp

    if (present(timestamp)) then
      call sfc_prop_restart_write_coarse (IPD_Data%Sfcprop, Atm_block, Model, &
            coarse_domain, IPD_Data%Grid, timestamp)
    else
      call sfc_prop_restart_write_coarse (IPD_Data%Sfcprop, Atm_block, Model, &
            coarse_domain, IPD_Data%Grid)
    endif
  end subroutine FV3GFS_restart_write_coarse

!--------------------
! FV3GFS_IPD_checksum
!--------------------
 subroutine FV3GFS_IPD_checksum (Model, IPD_Data, Atm_block)
   !--- interface variables
   type(IPD_control_type),    intent(in) :: Model
   type(IPD_data_type),       intent(in) :: IPD_Data(:)
   type (block_control_type), intent(in) :: Atm_block
   !--- local variables
   integer :: outunit, j, i, ix, nb, isc, iec, jsc, jec, lev, ct, l, ntr
   integer :: nsfcprop2d, idx_opt
   real(kind=kind_phys), allocatable :: temp2d(:,:,:)
   real(kind=kind_phys), allocatable :: temp3d(:,:,:,:)
   character(len=32) :: name

   isc = Model%isc
   iec = Model%isc+Model%nx-1
   jsc = Model%jsc
   jec = Model%jsc+Model%ny-1
   lev = Model%levs

   ntr = size(IPD_Data(1)%Statein%qgrs,3)

   if(Model%lsm == Model%lsm_noahmp) then
     nsfcprop2d = 154
   else
     nsfcprop2d = 100
   endif

   allocate (temp2d(isc:iec,jsc:jec,nsfcprop2d+Model%ntot3d+Model%nctp))
   allocate (temp3d(isc:iec,jsc:jec,1:lev,17+Model%ntot3d+2*ntr))

   temp2d = 0.
   temp3d = 0.

   do j=jsc,jec
     do i=isc,iec
       nb = Atm_block%blkno(i,j)
       ix = Atm_block%ixp(i,j)
       !--- statein pressure
       temp2d(i,j, 1) = IPD_Data(nb)%Statein%pgr(ix)
       temp2d(i,j, 2) = IPD_Data(nb)%Sfcprop%slmsk(ix)
       temp2d(i,j, 3) = IPD_Data(nb)%Sfcprop%tsfc(ix)
       temp2d(i,j, 4) = IPD_Data(nb)%Sfcprop%tisfc(ix)
       temp2d(i,j, 5) = IPD_Data(nb)%Sfcprop%snowd(ix)
       temp2d(i,j, 6) = IPD_Data(nb)%Sfcprop%zorl(ix)
       temp2d(i,j, 7) = IPD_Data(nb)%Sfcprop%fice(ix)
       temp2d(i,j, 8) = IPD_Data(nb)%Sfcprop%hprime(ix,1)
       temp2d(i,j, 9) = IPD_Data(nb)%Sfcprop%sncovr(ix)
       temp2d(i,j,10) = IPD_Data(nb)%Sfcprop%snoalb(ix)
       temp2d(i,j,11) = IPD_Data(nb)%Sfcprop%alvsf(ix)
       temp2d(i,j,12) = IPD_Data(nb)%Sfcprop%alnsf(ix)
       temp2d(i,j,13) = IPD_Data(nb)%Sfcprop%alvwf(ix)
       temp2d(i,j,14) = IPD_Data(nb)%Sfcprop%alnwf(ix)
       temp2d(i,j,15) = IPD_Data(nb)%Sfcprop%facsf(ix)
       temp2d(i,j,16) = IPD_Data(nb)%Sfcprop%facwf(ix)
       temp2d(i,j,17) = IPD_Data(nb)%Sfcprop%slope(ix)
       temp2d(i,j,18) = IPD_Data(nb)%Sfcprop%shdmin(ix)
       temp2d(i,j,19) = IPD_Data(nb)%Sfcprop%shdmax(ix)
       temp2d(i,j,20) = IPD_Data(nb)%Sfcprop%tg3(ix)
       temp2d(i,j,21) = IPD_Data(nb)%Sfcprop%vfrac(ix)
       temp2d(i,j,22) = IPD_Data(nb)%Sfcprop%vtype(ix)
       temp2d(i,j,23) = IPD_Data(nb)%Sfcprop%stype(ix)
       temp2d(i,j,24) = IPD_Data(nb)%Sfcprop%uustar(ix)
       temp2d(i,j,25) = IPD_Data(nb)%Sfcprop%oro(ix)
       temp2d(i,j,26) = IPD_Data(nb)%Sfcprop%oro_uf(ix)
       temp2d(i,j,27) = IPD_Data(nb)%Sfcprop%hice(ix)
       temp2d(i,j,28) = IPD_Data(nb)%Sfcprop%weasd(ix)
       temp2d(i,j,29) = IPD_Data(nb)%Sfcprop%canopy(ix)
       temp2d(i,j,30) = IPD_Data(nb)%Sfcprop%ffmm(ix)
       temp2d(i,j,31) = IPD_Data(nb)%Sfcprop%ffhh(ix)
       temp2d(i,j,32) = IPD_Data(nb)%Sfcprop%f10m(ix)
       temp2d(i,j,33) = IPD_Data(nb)%Sfcprop%tprcp(ix)
       temp2d(i,j,34) = IPD_Data(nb)%Sfcprop%srflag(ix)
       temp2d(i,j,35) = IPD_Data(nb)%Sfcprop%slc(ix,1)
       temp2d(i,j,36) = IPD_Data(nb)%Sfcprop%slc(ix,2)
       temp2d(i,j,37) = IPD_Data(nb)%Sfcprop%slc(ix,3)
       temp2d(i,j,38) = IPD_Data(nb)%Sfcprop%slc(ix,4)
       temp2d(i,j,39) = IPD_Data(nb)%Sfcprop%smc(ix,1)
       temp2d(i,j,40) = IPD_Data(nb)%Sfcprop%smc(ix,2)
       temp2d(i,j,41) = IPD_Data(nb)%Sfcprop%smc(ix,3)
       temp2d(i,j,42) = IPD_Data(nb)%Sfcprop%smc(ix,4)
       temp2d(i,j,43) = IPD_Data(nb)%Sfcprop%stc(ix,1)
       temp2d(i,j,44) = IPD_Data(nb)%Sfcprop%stc(ix,2)
       temp2d(i,j,45) = IPD_Data(nb)%Sfcprop%stc(ix,3)
       temp2d(i,j,46) = IPD_Data(nb)%Sfcprop%stc(ix,4)
       temp2d(i,j,47) = IPD_Data(nb)%Sfcprop%t2m(ix)
       temp2d(i,j,48) = IPD_Data(nb)%Sfcprop%q2m(ix)
       temp2d(i,j,49) = IPD_Data(nb)%Coupling%nirbmdi(ix)
       temp2d(i,j,50) = IPD_Data(nb)%Coupling%nirdfdi(ix)
       temp2d(i,j,51) = IPD_Data(nb)%Coupling%visbmdi(ix)
       temp2d(i,j,52) = IPD_Data(nb)%Coupling%visdfdi(ix)
       temp2d(i,j,53) = IPD_Data(nb)%Coupling%nirbmui(ix)
       temp2d(i,j,54) = IPD_Data(nb)%Coupling%nirdfui(ix)
       temp2d(i,j,55) = IPD_Data(nb)%Coupling%visbmui(ix)
       temp2d(i,j,56) = IPD_Data(nb)%Coupling%visdfui(ix)
       temp2d(i,j,57) = IPD_Data(nb)%Coupling%sfcdsw(ix)
       temp2d(i,j,58) = IPD_Data(nb)%Coupling%sfcnsw(ix)
       temp2d(i,j,59) = IPD_Data(nb)%Coupling%sfcdlw(ix)
       temp2d(i,j,60) = IPD_Data(nb)%Grid%xlon(ix)
       temp2d(i,j,61) = IPD_Data(nb)%Grid%xlat(ix)
       temp2d(i,j,62) = IPD_Data(nb)%Grid%xlat_d(ix)
       temp2d(i,j,63) = IPD_Data(nb)%Grid%sinlat(ix)
       temp2d(i,j,64) = IPD_Data(nb)%Grid%coslat(ix)
       temp2d(i,j,65) = IPD_Data(nb)%Grid%area(ix)
       temp2d(i,j,66) = IPD_Data(nb)%Grid%dx(ix)
       if (Model%ntoz > 0) then
         temp2d(i,j,67) = IPD_Data(nb)%Grid%ddy_o3(ix)
       endif
       if (Model%h2o_phys) then
         temp2d(i,j,68) = IPD_Data(nb)%Grid%ddy_h(ix)
       endif
       temp2d(i,j,69) = IPD_Data(nb)%Cldprop%cv(ix)
       temp2d(i,j,70) = IPD_Data(nb)%Cldprop%cvt(ix)
       temp2d(i,j,71) = IPD_Data(nb)%Cldprop%cvb(ix)
       temp2d(i,j,72) = IPD_Data(nb)%Radtend%sfalb(ix)
       temp2d(i,j,73) = IPD_Data(nb)%Radtend%coszen(ix)
       temp2d(i,j,74) = IPD_Data(nb)%Radtend%tsflw(ix)
       temp2d(i,j,75) = IPD_Data(nb)%Radtend%semis(ix)
       temp2d(i,j,76) = IPD_Data(nb)%Radtend%coszdg(ix)
       temp2d(i,j,77) = IPD_Data(nb)%Radtend%sfcfsw(ix)%upfxc
       temp2d(i,j,78) = IPD_Data(nb)%Radtend%sfcfsw(ix)%upfx0
       temp2d(i,j,79) = IPD_Data(nb)%Radtend%sfcfsw(ix)%dnfxc
       temp2d(i,j,80) = IPD_Data(nb)%Radtend%sfcfsw(ix)%dnfx0
       temp2d(i,j,81) = IPD_Data(nb)%Radtend%sfcflw(ix)%upfxc
       temp2d(i,j,82) = IPD_Data(nb)%Radtend%sfcflw(ix)%upfx0
       temp2d(i,j,83) = IPD_Data(nb)%Radtend%sfcflw(ix)%dnfxc
       temp2d(i,j,84) = IPD_Data(nb)%Radtend%sfcflw(ix)%dnfx0

        idx_opt = 85
       if (Model%lsm == Model%lsm_noahmp) then
        temp2d(i,j,idx_opt) = IPD_Data(nb)%Sfcprop%snowxy(ix)
        temp2d(i,j,idx_opt+1) = IPD_Data(nb)%Sfcprop%tvxy(ix)
        temp2d(i,j,idx_opt+2) = IPD_Data(nb)%Sfcprop%tgxy(ix)
        temp2d(i,j,idx_opt+3) = IPD_Data(nb)%Sfcprop%canicexy(ix)
        temp2d(i,j,idx_opt+4) = IPD_Data(nb)%Sfcprop%canliqxy(ix)
        temp2d(i,j,idx_opt+5) = IPD_Data(nb)%Sfcprop%eahxy(ix)
        temp2d(i,j,idx_opt+6) = IPD_Data(nb)%Sfcprop%tahxy(ix)
        temp2d(i,j,idx_opt+7) = IPD_Data(nb)%Sfcprop%cmxy(ix)
        temp2d(i,j,idx_opt+8) = IPD_Data(nb)%Sfcprop%chxy(ix)
        temp2d(i,j,idx_opt+9) = IPD_Data(nb)%Sfcprop%fwetxy(ix)
        temp2d(i,j,idx_opt+10) = IPD_Data(nb)%Sfcprop%sneqvoxy(ix)
        temp2d(i,j,idx_opt+11) = IPD_Data(nb)%Sfcprop%alboldxy(ix)
        temp2d(i,j,idx_opt+12) = IPD_Data(nb)%Sfcprop%qsnowxy(ix)
        temp2d(i,j,idx_opt+13) = IPD_Data(nb)%Sfcprop%wslakexy(ix)
        temp2d(i,j,idx_opt+14) = IPD_Data(nb)%Sfcprop%zwtxy(ix)
        temp2d(i,j,idx_opt+15) = IPD_Data(nb)%Sfcprop%waxy(ix)
        temp2d(i,j,idx_opt+16) = IPD_Data(nb)%Sfcprop%wtxy(ix)
        temp2d(i,j,idx_opt+17) = IPD_Data(nb)%Sfcprop%lfmassxy(ix)
        temp2d(i,j,idx_opt+18) = IPD_Data(nb)%Sfcprop%rtmassxy(ix)
        temp2d(i,j,idx_opt+19) = IPD_Data(nb)%Sfcprop%stmassxy(ix)
        temp2d(i,j,idx_opt+20) = IPD_Data(nb)%Sfcprop%woodxy(ix)
        temp2d(i,j,idx_opt+21) = IPD_Data(nb)%Sfcprop%stblcpxy(ix)
        temp2d(i,j,idx_opt+22) = IPD_Data(nb)%Sfcprop%fastcpxy(ix)
        temp2d(i,j,idx_opt+23) = IPD_Data(nb)%Sfcprop%xsaixy(ix)
        temp2d(i,j,idx_opt+24) = IPD_Data(nb)%Sfcprop%xlaixy(ix)
        temp2d(i,j,idx_opt+25) = IPD_Data(nb)%Sfcprop%taussxy(ix)
        temp2d(i,j,idx_opt+26) = IPD_Data(nb)%Sfcprop%smcwtdxy(ix)
        temp2d(i,j,idx_opt+27) = IPD_Data(nb)%Sfcprop%deeprechxy(ix)
        temp2d(i,j,idx_opt+28) = IPD_Data(nb)%Sfcprop%rechxy(ix)

        temp2d(i,j,idx_opt+29) = IPD_Data(nb)%Sfcprop%snicexy(ix,-2)
        temp2d(i,j,idx_opt+30) = IPD_Data(nb)%Sfcprop%snicexy(ix,-1)
        temp2d(i,j,idx_opt+31) = IPD_Data(nb)%Sfcprop%snicexy(ix,0)
        temp2d(i,j,idx_opt+32) = IPD_Data(nb)%Sfcprop%snliqxy(ix,-2)
        temp2d(i,j,idx_opt+33) = IPD_Data(nb)%Sfcprop%snliqxy(ix,-1)
        temp2d(i,j,idx_opt+34) = IPD_Data(nb)%Sfcprop%snliqxy(ix,0)
        temp2d(i,j,idx_opt+35) = IPD_Data(nb)%Sfcprop%tsnoxy(ix,-2)
        temp2d(i,j,idx_opt+36) = IPD_Data(nb)%Sfcprop%tsnoxy(ix,-1)
        temp2d(i,j,idx_opt+37) = IPD_Data(nb)%Sfcprop%tsnoxy(ix,0)
        temp2d(i,j,idx_opt+38) = IPD_Data(nb)%Sfcprop%smoiseq(ix,1)
        temp2d(i,j,idx_opt+39) = IPD_Data(nb)%Sfcprop%smoiseq(ix,2)
        temp2d(i,j,idx_opt+40) = IPD_Data(nb)%Sfcprop%smoiseq(ix,3)
        temp2d(i,j,idx_opt+41) = IPD_Data(nb)%Sfcprop%smoiseq(ix,4)
        temp2d(i,j,idx_opt+42) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,-2)
        temp2d(i,j,idx_opt+43) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,-1)
        temp2d(i,j,idx_opt+44) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,0)
        temp2d(i,j,idx_opt+45) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,1)
        temp2d(i,j,idx_opt+46) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,2)
        temp2d(i,j,idx_opt+47) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,3)
        temp2d(i,j,idx_opt+48) = IPD_Data(nb)%Sfcprop%zsnsoxy(ix,4)
        temp2d(i,j,idx_opt+49) = IPD_Data(nb)%Sfcprop%albdvis(ix)
        temp2d(i,j,idx_opt+50) = IPD_Data(nb)%Sfcprop%albdnir(ix)
        temp2d(i,j,idx_opt+51) = IPD_Data(nb)%Sfcprop%albivis(ix)
        temp2d(i,j,idx_opt+52) = IPD_Data(nb)%Sfcprop%albinir(ix)
        temp2d(i,j,idx_opt+53) = IPD_Data(nb)%Sfcprop%emiss(ix)
        idx_opt = 139
       endif

       if (Model%nstf_name(1) > 0) then
         temp2d(i,j,idx_opt) = IPD_Data(nb)%Sfcprop%tref(ix)
         temp2d(i,j,idx_opt+1) = IPD_Data(nb)%Sfcprop%z_c(ix)
         temp2d(i,j,idx_opt+2) = IPD_Data(nb)%Sfcprop%c_0(ix)
         temp2d(i,j,idx_opt+3) = IPD_Data(nb)%Sfcprop%c_d(ix)
         temp2d(i,j,idx_opt+4) = IPD_Data(nb)%Sfcprop%w_0(ix)
         temp2d(i,j,idx_opt+5) = IPD_Data(nb)%Sfcprop%w_d(ix)
         temp2d(i,j,idx_opt+6) = IPD_Data(nb)%Sfcprop%xt(ix)
         temp2d(i,j,idx_opt+7) = IPD_Data(nb)%Sfcprop%xs(ix)
         temp2d(i,j,idx_opt+8) = IPD_Data(nb)%Sfcprop%xu(ix)
         temp2d(i,j,idx_opt+9) = IPD_Data(nb)%Sfcprop%xz(ix)
         temp2d(i,j,idx_opt+10) = IPD_Data(nb)%Sfcprop%zm(ix)
         temp2d(i,j,idx_opt+11) = IPD_Data(nb)%Sfcprop%xtts(ix)
         temp2d(i,j,idx_opt+12) = IPD_Data(nb)%Sfcprop%xzts(ix)
         temp2d(i,j,idx_opt+13) = IPD_Data(nb)%Sfcprop%ifd(ix)
         temp2d(i,j,idx_opt+14) = IPD_Data(nb)%Sfcprop%dt_cool(ix)
         temp2d(i,j,idx_opt+15) = IPD_Data(nb)%Sfcprop%qrain(ix)
       endif

       do l = 1,Model%ntot2d
         temp2d(i,j,nsfcprop2d+l) = IPD_Data(nb)%Tbd%phy_f2d(ix,l)
       enddo

       do l = 1,Model%nctp
         temp2d(i,j,nsfcprop2d+Model%ntot2d+l) = IPD_Data(nb)%Tbd%phy_fctd(ix,l)
       enddo

       temp3d(i,j,:, 1) = IPD_Data(nb)%Statein%phii(ix,1:lev)
       temp3d(i,j,:, 2) = IPD_Data(nb)%Statein%prsi(ix,1:lev)
       temp3d(i,j,:, 3) = IPD_Data(nb)%Statein%prsik(ix,1:lev)
       temp3d(i,j,:, 4) = IPD_Data(nb)%Statein%phil(ix,:)
       temp3d(i,j,:, 5) = IPD_Data(nb)%Statein%prsl(ix,:)
       temp3d(i,j,:, 6) = IPD_Data(nb)%Statein%prslk(ix,:)
       temp3d(i,j,:, 7) = IPD_Data(nb)%Statein%ugrs(ix,:)
       temp3d(i,j,:, 8) = IPD_Data(nb)%Statein%vgrs(ix,:)
       temp3d(i,j,:, 9) = IPD_Data(nb)%Statein%vvl(ix,:)
       temp3d(i,j,:,10) = IPD_Data(nb)%Statein%tgrs(ix,:)
       temp3d(i,j,:,11) = IPD_Data(nb)%Stateout%gu0(ix,:)
       temp3d(i,j,:,12) = IPD_Data(nb)%Stateout%gv0(ix,:)
       temp3d(i,j,:,13) = IPD_Data(nb)%Stateout%gt0(ix,:)
       temp3d(i,j,:,14) = IPD_Data(nb)%Radtend%htrsw(ix,:)
       temp3d(i,j,:,15) = IPD_Data(nb)%Radtend%htrlw(ix,:)
       temp3d(i,j,:,16) = IPD_Data(nb)%Radtend%swhc(ix,:)
       temp3d(i,j,:,17) = IPD_Data(nb)%Radtend%lwhc(ix,:)
       do l = 1,Model%ntot3d
         temp3d(i,j,:,17+l) = IPD_Data(nb)%Tbd%phy_f3d(ix,:,l)
       enddo
       do l = 1,ntr
         temp3d(i,j,:,17+Model%ntot3d+l)     = IPD_Data(nb)%Statein%qgrs(ix,:,l)
         temp3d(i,j,:,17+Model%ntot3d+ntr+l) = IPD_Data(nb)%Stateout%gq0(ix,:,l)
       enddo
     enddo
   enddo

   outunit = stdout()
   do i = 1,nsfcprop2d+Model%ntot2d+Model%nctp
     write (name, '(i3.3,3x,4a)') i, ' 2d '
     write(outunit,100) name, mpp_chksum(temp2d(:,:,i:i))
   enddo
   do i = 1,17+Model%ntot3d+2*ntr
     write (name, '(i2.2,3x,4a)') i, ' 3d '
     write(outunit,100) name, mpp_chksum(temp3d(:,:,:,i:i))
   enddo
100 format("CHECKSUM::",A32," = ",Z20)

   end subroutine FV3GFS_IPD_checksum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
! register_sfc_prop_restart_vars
!----------------------------------------------------------------------
!    creates and populates a data type for surface variables which is
!    then used to "register" restart variables with the GFDL FMS
!    restart subsystem.
!
!    calls:  register_restart_field
!
!----------------------------------------------------------------------
   subroutine register_sfc_prop_restart_vars(Model, nx, ny, nvar_s2m, action)
    type(IPD_control_type), intent(in)  :: Model
    integer,                intent(in)  :: nx
    integer,                intent(in)  :: ny
    integer,                intent(out) :: nvar_s2m
    character(len=*),       intent(in)  :: action  !< alloc, read, write

    !--- local variables
    integer :: is, ie
    integer :: lsoil, num
    integer :: nvar_s2mp, nvar_s2o
    integer :: nvar_s3,  nvar_s3mp
    logical :: opt
    character(len=8) :: dim_names_2d(3), dim_names_3d(4)
    integer, allocatable, dimension(:) :: buffer
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()

    if (Model%cplflx) then ! needs more variables
      nvar_s2m = 34
    else
      nvar_s2m = 32
    endif
    nvar_s2o = 18
    nvar_s3  = 3

    if (Model%lsm == Model%lsm_noahmp) then
      nvar_s2mp = 40       !mp 2D
      nvar_s3mp = 5        !mp 3D
    else
      nvar_s2mp = 0        !mp 2D
      nvar_s3mp = 0        !mp 3D
    endif

    if (Model%use_ifs_ini_sst) then
      allocate(ifsSST(nx,ny))
    endif

    if (.not. allocated(sfc_name2)) then
      !--- allocate the various containers needed for restarts
      allocate(sfc_name2(nvar_s2m+nvar_s2o+nvar_s2mp))
      allocate(sfc_name3(nvar_s3+nvar_s3mp))

      allocate(sfc_var2(nx,ny,nvar_s2m+nvar_s2o+nvar_s2mp))
      allocate(sfc_var3(nx,ny,Model%lsoil,nvar_s3))
      sfc_var2 = -9999._kind_phys
      sfc_var3 = -9999._kind_phys
!
      if (Model%lsm == Model%lsm_noahmp) then
        allocate(sfc_var3sn(nx,ny,-2:0,4:6))
        allocate(sfc_var3eq(nx,ny,1:4,7:7))
        allocate(sfc_var3zn(nx,ny,-2:4,8:8))
        sfc_var3sn = -9999._kind_phys
        sfc_var3eq = -9999._kind_phys
        sfc_var3zn = -9999._kind_phys
      end if

      !--- names of the 2D variables to save
      sfc_name2(1)  = 'slmsk'
      sfc_name2(2)  = 'tsea'    !tsfc
      sfc_name2(3)  = 'sheleg'  !weasd
      sfc_name2(4)  = 'tg3'
      sfc_name2(5)  = 'zorl'
      sfc_name2(6)  = 'alvsf'
      sfc_name2(7)  = 'alvwf'
      sfc_name2(8)  = 'alnsf'
      sfc_name2(9)  = 'alnwf'
      sfc_name2(10) = 'facsf'
      sfc_name2(11) = 'facwf'
      sfc_name2(12) = 'vfrac'
      sfc_name2(13) = 'canopy'
      sfc_name2(14) = 'f10m'
      sfc_name2(15) = 't2m'
      sfc_name2(16) = 'q2m'
      sfc_name2(17) = 'vtype'
      sfc_name2(18) = 'stype'
      sfc_name2(19) = 'uustar'
      sfc_name2(20) = 'ffmm'
      sfc_name2(21) = 'ffhh'
      sfc_name2(22) = 'hice'
      sfc_name2(23) = 'fice'
      sfc_name2(24) = 'tisfc'
      sfc_name2(25) = 'tprcp'
      sfc_name2(26) = 'srflag'
      sfc_name2(27) = 'snwdph'  !snowd
      sfc_name2(28) = 'shdmin'
      sfc_name2(29) = 'shdmax'
      sfc_name2(30) = 'slope'
      sfc_name2(31) = 'snoalb'
      !--- variables below here are optional
      sfc_name2(32) = 'sncovr'
      if(Model%cplflx) then
        sfc_name2(33) = 'tsfcl'   !temp on land portion of a cell
        sfc_name2(34) = 'zorll'   !zorl on land portion of a cell
      end if

      !--- NSSTM inputs only needed when (nstf_name(1) > 0) .and. (nstf_name(2)) == 0)
      sfc_name2(nvar_s2m+1)  = 'tref'
      sfc_name2(nvar_s2m+2)  = 'z_c'
      sfc_name2(nvar_s2m+3)  = 'c_0'
      sfc_name2(nvar_s2m+4)  = 'c_d'
      sfc_name2(nvar_s2m+5)  = 'w_0'
      sfc_name2(nvar_s2m+6)  = 'w_d'
      sfc_name2(nvar_s2m+7)  = 'xt'
      sfc_name2(nvar_s2m+8)  = 'xs'
      sfc_name2(nvar_s2m+9)  = 'xu'
      sfc_name2(nvar_s2m+10) = 'xv'
      sfc_name2(nvar_s2m+11) = 'xz'
      sfc_name2(nvar_s2m+12) = 'zm'
      sfc_name2(nvar_s2m+13) = 'xtts'
      sfc_name2(nvar_s2m+14) = 'xzts'
      sfc_name2(nvar_s2m+15) = 'd_conv'
      sfc_name2(nvar_s2m+16) = 'ifd'
      sfc_name2(nvar_s2m+17) = 'dt_cool'
      sfc_name2(nvar_s2m+18) = 'qrain'
      !
      ! Only needed when Noah MP LSM is used - 39 2D
      !
      if (Model%lsm == Model%lsm_noahmp) then
        sfc_name2(nvar_s2m+19) = 'snowxy'
        sfc_name2(nvar_s2m+20) = 'tvxy'
        sfc_name2(nvar_s2m+21) = 'tgxy'
        sfc_name2(nvar_s2m+22) = 'canicexy'
        sfc_name2(nvar_s2m+23) = 'canliqxy'
        sfc_name2(nvar_s2m+24) = 'eahxy'
        sfc_name2(nvar_s2m+25) = 'tahxy'
        sfc_name2(nvar_s2m+26) = 'cmxy'
        sfc_name2(nvar_s2m+27) = 'chxy'
        sfc_name2(nvar_s2m+28) = 'fwetxy'
        sfc_name2(nvar_s2m+29) = 'sneqvoxy'
        sfc_name2(nvar_s2m+30) = 'alboldxy'
        sfc_name2(nvar_s2m+31) = 'qsnowxy'
        sfc_name2(nvar_s2m+32) = 'wslakexy'
        sfc_name2(nvar_s2m+33) = 'zwtxy'
        sfc_name2(nvar_s2m+34) = 'waxy'
        sfc_name2(nvar_s2m+35) = 'wtxy'
        sfc_name2(nvar_s2m+36) = 'lfmassxy'
        sfc_name2(nvar_s2m+37) = 'rtmassxy'
        sfc_name2(nvar_s2m+38) = 'stmassxy'
        sfc_name2(nvar_s2m+39) = 'woodxy'
        sfc_name2(nvar_s2m+40) = 'stblcpxy'
        sfc_name2(nvar_s2m+41) = 'fastcpxy'
        sfc_name2(nvar_s2m+42) = 'xsaixy'
        sfc_name2(nvar_s2m+43) = 'xlaixy'
        sfc_name2(nvar_s2m+44) = 'taussxy'
        sfc_name2(nvar_s2m+45) = 'smcwtdxy'
        sfc_name2(nvar_s2m+46) = 'deeprechxy'
        sfc_name2(nvar_s2m+47) = 'rechxy'
        sfc_name2(nvar_s2m+48) = 'drainncprv'
        sfc_name2(nvar_s2m+49) = 'draincprv'
        sfc_name2(nvar_s2m+50) = 'dsnowprv'
        sfc_name2(nvar_s2m+51) = 'dgraupelprv'
        sfc_name2(nvar_s2m+52) = 'diceprv'
        sfc_name2(nvar_s2m+53) = 'albdvis'
        sfc_name2(nvar_s2m+54) = 'albdnir'
        sfc_name2(nvar_s2m+55) = 'albivis'
        sfc_name2(nvar_s2m+56) = 'albinir'
        sfc_name2(nvar_s2m+57) = 'emiss'
        sfc_name2(nvar_s2m+58) = 'scolor'
      endif

      !--- names of the 3D variables to save
      sfc_name3(1) = 'stc'
      sfc_name3(2) = 'smc'
      sfc_name3(3) = 'slc'
      !--- Noah MP
      if (Model%lsm == Model%lsm_noahmp) then
        sfc_name3(4) = 'snicexy'
        sfc_name3(5) = 'snliqxy'
        sfc_name3(6) = 'tsnoxy'
        sfc_name3(7) = 'smoiseq'
        sfc_name3(8) = 'zsnsoxy'
      endif
    endif  ! if not allocated

    if (trim(action) == "alloc") then
      return
    elseif (trim(action) == "read") then
      !--- register the axes for restarts
      if (dimension_exists(Sfc_restart, "xaxis_1")) then
        call register_axis(Sfc_restart, "xaxis_1", "X")
        call register_axis(Sfc_restart, "yaxis_1", "Y")
        call register_axis(Sfc_restart, "zaxis_1", dimension_length=4)
        call register_axis(Sfc_restart, "zaxis_2", dimension_length=3)
        call register_axis(Sfc_restart, "zaxis_3", dimension_length=7)
        call register_axis(Sfc_restart, "Time", unlimited)
        dim_names_2d(1) = "xaxis_1"
        dim_names_2d(2) = "yaxis_1"
        dim_names_2d(3) = "Time"
        dim_names_3d(1) = "xaxis_1"
        dim_names_3d(2) = "yaxis_1"
        dim_names_3d(3) = "zaxis_1"     ! to be reset as needed when variables are registered
        dim_names_3d(4) = "Time"
      else
        call register_axis(Sfc_restart, 'lon', 'X')
        call register_axis(Sfc_restart, 'lat', 'Y')
        call register_axis(Sfc_restart, 'lsoil', dimension_length=Model%lsoil)
        dim_names_2d(1) = "lat"
        dim_names_2d(2) = "lon"
        dim_names_3d(1) = "lat"
        dim_names_3d(2) = "lon"
        dim_names_3d(3) = "lsoil"
      endif

    elseif (trim(action) == "write") then

      !--- register the axes for restarts
      call register_axis(Sfc_restart, 'xaxis_1', 'X')
      call register_field(Sfc_restart, 'xaxis_1', 'double', (/'xaxis_1'/))
      call register_variable_attribute(Sfc_restart, 'xaxis_1', 'cartesian_axis', 'X', str_len=1)
      call get_global_io_domain_indices(Sfc_restart, 'xaxis_1', is, ie, indices=buffer)
      call write_data(Sfc_restart, "xaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Sfc_restart, 'yaxis_1', 'Y')
      call register_field(Sfc_restart, 'yaxis_1', 'double', (/'yaxis_1'/))
      call register_variable_attribute(Sfc_restart, 'yaxis_1', 'cartesian_axis', 'Y', str_len=1)
      call get_global_io_domain_indices(Sfc_restart, 'yaxis_1', is, ie, indices=buffer)
      call write_data(Sfc_restart, "yaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Sfc_restart, 'zaxis_1', dimension_length=Model%lsoil)
      call register_field(Sfc_restart, 'zaxis_1', 'double', (/'zaxis_1'/))
      call register_variable_attribute(Sfc_restart, 'zaxis_1', 'cartesian_axis', 'Z', str_len=1)
      allocate( buffer(Model%lsoil) )
      do lsoil=1, Model%lsoil
         buffer(lsoil) = lsoil
      end do
      call write_data(Sfc_restart, 'zaxis_1', buffer)
      deallocate(buffer)

      call register_axis(Sfc_restart, 'Time', unlimited)
      call register_field(Sfc_restart, 'Time', 'double', (/'Time'/))
      call register_variable_attribute(Sfc_restart, 'Time', 'cartesian_axis', 'T', str_len=1)
      call write_data(Sfc_restart, 'Time', 1)

      if (Model%lsm == Model%lsm_noahmp) then
        call register_axis(Sfc_restart, 'zaxis_2', dimension_length=3)
        call register_field(Sfc_restart, 'zaxis_2', 'double', (/'zaxis_2'/))
        call register_variable_attribute(Sfc_restart, 'zaxis_2', 'cartesian_axis', 'Z', str_len=1)
        allocate( buffer(3) )
        do lsoil=-2,0
           buffer(lsoil+3) = lsoil
        end do
        call write_data(Sfc_restart, 'zaxis_2', buffer)
        deallocate(buffer)

        call register_axis(Sfc_restart, 'zaxis_3', dimension_length=7)
        call register_field(Sfc_restart, 'zaxis_3', 'double', (/'zaxis_3'/))
        call register_variable_attribute(Sfc_restart, 'zaxis_3', 'cartesian_axis', 'Z', str_len=1)
        allocate( buffer(7) )
        do lsoil=-2,4
           buffer(lsoil+3) = lsoil
        end do
        call write_data(Sfc_restart, 'zaxis_3', buffer)
        deallocate(buffer)
      endif   ! if (lsm_noahmp)

      dim_names_2d(1) = "xaxis_1"
      dim_names_2d(2) = "yaxis_1"
      dim_names_2d(3) = "Time"
      dim_names_3d(1) = "xaxis_1"
      dim_names_3d(2) = "yaxis_1"
      dim_names_3d(3) = "zaxis_1"     ! to be reset as needed when variables are registered
      dim_names_3d(4) = "Time"

    else  ! error case

      call mpp_error(FATAL,"FV3GFS_io::register_sfc_prop_restart_vars action not found")

    endif  ! end of if (read)

    !--- register IFS SST
    if (Model%use_ifs_ini_sst) then
      var2_p => ifsSST
      opt = .false.
      call register_restart_field(ifsSST_restart, 'sst', var2_p, dim_names_2d, is_optional=opt)
      nullify(var2_p)
    endif

    !--- register the 2D fields
    do num = 1,nvar_s2m
      var2_p => sfc_var2(:,:,num)
      opt = .false.
      if (trim(sfc_name2(num)) == 'sncovr'.or.trim(sfc_name2(num)) == 'tsfcl'.or.trim(sfc_name2(num)) == 'zorll') opt = .true.
      call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dim_names_2d, is_optional=opt)
      nullify(var2_p)
    enddo

    !--- register NSST variables
    if (Model%nstf_name(1) > 0) then
      opt = .true.
      if (Model%nstf_name(2) == 0) opt = .false.
      do num = nvar_s2m+1,nvar_s2m+nvar_s2o
        var2_p => sfc_var2(:,:,num)
        call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dim_names_2d, is_optional=opt)
        nullify(var2_p)
      enddo
    endif

    !--- Noah MP register only necessary only lsm = 2, not necessary has values
    if (nvar_s2mp > 0) then
      opt = .true.
      do num = nvar_s2m+nvar_s2o+1,nvar_s2m+nvar_s2o+nvar_s2mp
        var2_p => sfc_var2(:,:,num)
        call register_restart_field(Sfc_restart, sfc_name2(num), var2_p, dim_names_2d, is_optional=opt)
        nullify(var2_p)
      enddo
    endif ! noahmp

    !--- register the 3D fields
    do num = 1,nvar_s3
      var3_p => sfc_var3(:,:,:,num)
      dim_names_3d(3) = "zaxis_1"
      call register_restart_field(Sfc_restart, sfc_name3(num), var3_p, dim_names_3d)
      nullify(var3_p)
    enddo

    !--- register the NOAH-MP 3D fields
    if (Model%lsm == Model%lsm_noahmp) then
      opt = .true.
      do num = nvar_s3+1,nvar_s3+3
        var3_p => sfc_var3sn(:,:,:,num)
        dim_names_3d(3) = "zaxis_2"
        call register_restart_field(Sfc_restart, sfc_name3(num), var3_p, dim_names_3d, is_optional=opt)
        nullify(var3_p)
      enddo

      var3_p => sfc_var3eq(:,:,:,7)
      dim_names_3d(3) = "zaxis_1"
      call register_restart_field(Sfc_restart, sfc_name3(7), var3_p, dim_names_3d, is_optional=opt)
      nullify(var3_p)

      var3_p => sfc_var3zn(:,:,:,8)
      dim_names_3d(3) = "zaxis_3"
      call register_restart_fIeld(Sfc_restart, sfc_name3(8), var3_p, dim_names_3d, is_optional=opt)
      nullify(var3_p)
    endif   !mp

   end subroutine register_sfc_prop_restart_vars

!----------------------------------------------------------------------
! sfc_prop_restart_read
!----------------------------------------------------------------------
!    calls a routine to "register" restart variables with the GFDL FMS
!    restart subsystem.
!
!    calls a GFDL FMS routine to restore the data from a restart file.
!    calculates sncovr if it is not present in the restart file.
!
!    calls:  open_file, register_sfc_prop_restart_vars, read_restart,
!            close_file
!
!    opens:  oro_data.tile?.nc, sfc_data.tile?.nc
!
!----------------------------------------------------------------------
  subroutine sfc_prop_restart_read (Sfcprop, Atm_block, Model, fv_domain, enforce_rst_cksum)
    !--- interface variable definitions
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop(:)
    type (block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),    intent(inout) :: Model
    type (domain2d),           intent(in)    :: fv_domain
    logical,                   intent(in)    :: enforce_rst_cksum
    !--- local variables
    integer :: i, j, k, ix, lsoil, num, nb
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: nvar_o2, nvar_s2m
    integer :: isnow
    logical :: opt
    character(len=8) :: dim_names_2d(3)
    real(kind=kind_phys), pointer, dimension(:,:) :: var2_p => NULL()

    character(len=64) :: fname
    !--- local variables for sncovr calculation
    integer :: vegtyp

    real(kind=kind_phys) :: rsnow, tem
    !--- Noah MP
    integer              :: soiltyp,ns,imon,iter,imn
    real(kind=kind_phys) :: masslai, masssai,snd
    real(kind=kind_phys) :: ddz,expon,aa,bb,smc,func,dfunc,dx
    real(kind=kind_phys) :: bexp, smcmax, smcwlt,dwsat,dksat,psisat
    real(kind=kind_phys) :: emg, emv

    real(kind=kind_phys), dimension(-2:0) :: dzsno
    real(kind=kind_phys), dimension(-2:4) :: dzsnso

    real(kind=kind_phys), dimension(4), save :: zsoil,dzs
    data dzs   /0.1,0.3,0.6,1.0/
    data zsoil /-0.1,-0.4,-1.0,-2.0/


    nvar_o2  = 19

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)

    !--- Open the restart file and associate it with the Oro_restart fileobject
    fname='INPUT/'//trim(fn_oro)
    if (open_file(Oro_restart, fname, "read", fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)) then
      call register_axis(Oro_restart, "lat", "y")
      call register_axis(Oro_restart, "lon", "x")
      dim_names_2d(1) = "lat"
      dim_names_2d(2) = "lon"

      !--- OROGRAPHY FILE
      if (.not. allocated(oro_name2)) then
      !--- allocate the various containers needed for orography data
        allocate(oro_name2(nvar_o2))
        allocate(oro_var2(nx,ny,nvar_o2))
        oro_var2 = -9999._kind_phys

        oro_name2(1)  = 'stddev'     ! hprime(ix,1)
        oro_name2(2)  = 'convexity'  ! hprime(ix,2)
        oro_name2(3)  = 'oa1'        ! hprime(ix,3)
        oro_name2(4)  = 'oa2'        ! hprime(ix,4)
        oro_name2(5)  = 'oa3'        ! hprime(ix,5)
        oro_name2(6)  = 'oa4'        ! hprime(ix,6)
        oro_name2(7)  = 'ol1'        ! hprime(ix,7)
        oro_name2(8)  = 'ol2'        ! hprime(ix,8)
        oro_name2(9)  = 'ol3'        ! hprime(ix,9)
        oro_name2(10) = 'ol4'        ! hprime(ix,10)
        oro_name2(11) = 'theta'      ! hprime(ix,11)
        oro_name2(12) = 'gamma'      ! hprime(ix,12)
        oro_name2(13) = 'sigma'      ! hprime(ix,13)
        oro_name2(14) = 'elvmax'     ! hprime(ix,14)
        oro_name2(15) = 'orog_filt'  ! oro
        oro_name2(16) = 'orog_raw'   ! oro_uf
        oro_name2(17) = 'land_frac'  ! land fraction [0:1]
        !--- variables below here are optional
        oro_name2(18) = 'lake_frac'  ! lake fraction [0:1]
        oro_name2(19) = 'lake_depth' ! lake depth(m)
      !--- register the 2D fields
        do num = 1,nvar_o2
          var2_p => oro_var2(:,:,num)
          opt = .false.
          if (trim(oro_name2(num)) == 'lake_frac' .or. trim(oro_name2(num)) == 'lake_depth') opt = .true.
          call register_restart_field(Oro_restart, oro_name2(num), var2_p, dim_names_2d, is_optional=opt)
          nullify(var2_p)
        enddo
      endif

      !--- read the orography restart/data
      call mpp_error(NOTE,'reading topographic/orographic information from INPUT/oro_data.tile*.nc')
      call read_restart(Oro_restart, ignore_checksum=enforce_rst_cksum)
      call close_file(Oro_restart)

      !--- copy data into GFS containers
      do nb = 1, Atm_block%nblks
        !--- 2D variables
        do ix = 1, Atm_block%blksz(nb)
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          !--- stddev
          Sfcprop(nb)%hprim(ix)     = oro_var2(i,j,1)
          !--- hprime(1:14)
          Sfcprop(nb)%hprime(ix,1)  = oro_var2(i,j,1)
          Sfcprop(nb)%hprime(ix,2)  = oro_var2(i,j,2)
          Sfcprop(nb)%hprime(ix,3)  = oro_var2(i,j,3)
          Sfcprop(nb)%hprime(ix,4)  = oro_var2(i,j,4)
          Sfcprop(nb)%hprime(ix,5)  = oro_var2(i,j,5)
          Sfcprop(nb)%hprime(ix,6)  = oro_var2(i,j,6)
          Sfcprop(nb)%hprime(ix,7)  = oro_var2(i,j,7)
          Sfcprop(nb)%hprime(ix,8)  = oro_var2(i,j,8)
          Sfcprop(nb)%hprime(ix,9)  = oro_var2(i,j,9)
          Sfcprop(nb)%hprime(ix,10) = oro_var2(i,j,10)
          Sfcprop(nb)%hprime(ix,11) = oro_var2(i,j,11)
          Sfcprop(nb)%hprime(ix,12) = oro_var2(i,j,12)
          Sfcprop(nb)%hprime(ix,13) = oro_var2(i,j,13)
          Sfcprop(nb)%hprime(ix,14) = oro_var2(i,j,14)
          !--- oro
          Sfcprop(nb)%oro(ix)       = oro_var2(i,j,15)
          !--- oro_uf
          Sfcprop(nb)%oro_uf(ix)    = oro_var2(i,j,16)
          Sfcprop(nb)%landfrac(ix)  = oro_var2(i,j,17) !land frac [0:1]
          Sfcprop(nb)%lakefrac(ix)  = oro_var2(i,j,18) !lake frac [0:1]
        enddo
      enddo

      if (nint(oro_var2(1,1,18)) == -9999._kind_phys) then ! lakefrac doesn't exist in the restart, need to create it
        if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - will computing lakefrac')
        Model%frac_grid = .false.
      else
        Model%frac_grid = .true.
      endif

      if (Model%me == Model%master ) write(0,*)' resetting Model%frac_grid=',Model%frac_grid

      !--- deallocate containers
      deallocate(oro_name2, oro_var2)

    else ! no file - cold_start (no way yet to create orography on-the-fly)

      call mpp_error(NOTE,'No INPUT/oro_data.tile*.nc orographic data found; setting to 0')
      !--- copy data into GFS containers
      do nb = 1, Atm_block%nblks
        !--- 2D variables
        do ix = 1, Atm_block%blksz(nb)
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          !--- stddev
          Sfcprop(nb)%hprim(ix)      = 0.0
          !--- hprime(1:14)
          Sfcprop(nb)%hprime(ix,1:14)  = 0.0
          !--- oro
          Sfcprop(nb)%oro(ix)        = 0.0
          !--- oro_uf
          Sfcprop(nb)%oro_uf(ix)     = 0.0
        enddo
      enddo

    endif

    if (Model%use_ifs_ini_sst) then
      !--- Open the restart file and associate it with the ifsSST_restart fileobject
      fname='INPUT/'//trim(fn_ifsSST)
      if (open_file(ifsSST_restart, fname, "read", fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)) then
        !--- read the IFS SST restart/data
        call mpp_error(NOTE,'reading ifs SST data from INPUT/ifsSST_data.tile*.nc')
        call read_restart(ifsSST_restart, ignore_checksum=enforce_rst_cksum)
        call close_file(ifsSST_restart)
      else
        call mpp_error(FATAL,'No ifs SST data.')
      endif
    endif

    !--- Open the restart file and associate it with the Sfc_restart fileobject
    fname='INPUT/'//trim(fn_srf)
    if (open_file(Sfc_restart, fname, "read", fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)) then

      !--- register any variables needed for the restart read
      call register_sfc_prop_restart_vars(Model, nx, ny, nvar_s2m, action="read")

      !--- read the surface restart/data
      call mpp_error(NOTE,'reading surface properties data from INPUT/sfc_data.tile*.nc')
      call read_restart(Sfc_restart, ignore_checksum=enforce_rst_cksum)
      call close_file(Sfc_restart)

      !--- place the data into the block GFS containers
      do nb = 1, Atm_block%nblks
         do ix = 1, Atm_block%blksz(nb)
           i = Atm_block%index(nb)%ii(ix) - isc + 1
           j = Atm_block%index(nb)%jj(ix) - jsc + 1

!--- 2D variables
!    ------------
           Sfcprop(nb)%slmsk(ix)  = sfc_var2(i,j,1)    !--- slmsk
           if (Model%use_ifs_ini_sst) then
              Sfcprop(nb)%tsfco(ix)  = ifsSST(i,j)    !--- tsfc (sst in ifsSST file)
           else
              Sfcprop(nb)%tsfco(ix)  = sfc_var2(i,j,2)    !--- tsfc (tsea in sfc file)
           endif
           Sfcprop(nb)%weasd(ix)  = sfc_var2(i,j,3)    !--- weasd (sheleg in sfc file)
           Sfcprop(nb)%tg3(ix)    = sfc_var2(i,j,4)    !--- tg3
           Sfcprop(nb)%zorlo(ix)  = sfc_var2(i,j,5)    !--- zorl on ocean
           Sfcprop(nb)%alvsf(ix)  = sfc_var2(i,j,6)    !--- alvsf
           Sfcprop(nb)%alvwf(ix)  = sfc_var2(i,j,7)    !--- alvwf
           Sfcprop(nb)%alnsf(ix)  = sfc_var2(i,j,8)    !--- alnsf
           Sfcprop(nb)%alnwf(ix)  = sfc_var2(i,j,9)    !--- alnwf
           Sfcprop(nb)%facsf(ix)  = sfc_var2(i,j,10)   !--- facsf
           Sfcprop(nb)%facwf(ix)  = sfc_var2(i,j,11)   !--- facwf
           Sfcprop(nb)%vfrac(ix)  = sfc_var2(i,j,12)   !--- vfrac
           Sfcprop(nb)%canopy(ix) = sfc_var2(i,j,13)   !--- canopy
           Sfcprop(nb)%f10m(ix)   = sfc_var2(i,j,14)   !--- f10m
           Sfcprop(nb)%t2m(ix)    = sfc_var2(i,j,15)   !--- t2m
           Sfcprop(nb)%q2m(ix)    = sfc_var2(i,j,16)   !--- q2m
           Sfcprop(nb)%vtype(ix)  = sfc_var2(i,j,17)   !--- vtype
           Sfcprop(nb)%stype(ix)  = sfc_var2(i,j,18)   !--- stype
           Sfcprop(nb)%uustar(ix) = sfc_var2(i,j,19)   !--- uustar
           Sfcprop(nb)%ffmm(ix)   = sfc_var2(i,j,20)   !--- ffmm
           Sfcprop(nb)%ffhh(ix)   = sfc_var2(i,j,21)   !--- ffhh
           Sfcprop(nb)%hice(ix)   = sfc_var2(i,j,22)   !--- hice
           Sfcprop(nb)%fice(ix)   = sfc_var2(i,j,23)   !--- fice
           Sfcprop(nb)%tisfc(ix)  = sfc_var2(i,j,24)   !--- tisfc
           Sfcprop(nb)%tprcp(ix)  = sfc_var2(i,j,25)   !--- tprcp
           Sfcprop(nb)%srflag(ix) = sfc_var2(i,j,26)   !--- srflag
           Sfcprop(nb)%snowd(ix)  = sfc_var2(i,j,27)   !--- snowd (snwdph in the file)
           Sfcprop(nb)%shdmin(ix) = sfc_var2(i,j,28)   !--- shdmin
           Sfcprop(nb)%shdmax(ix) = sfc_var2(i,j,29)   !--- shdmax
           Sfcprop(nb)%slope(ix)  = sfc_var2(i,j,30)   !--- slope
           Sfcprop(nb)%snoalb(ix) = sfc_var2(i,j,31)   !--- snoalb
           Sfcprop(nb)%sncovr(ix) = sfc_var2(i,j,32)   !--- sncovr
           if(Model%cplflx) then
             Sfcprop(nb)%tsfcl(ix)  = sfc_var2(i,j,33) !--- sfcl  (temp on land portion of a cell)
             Sfcprop(nb)%zorll(ix)  = sfc_var2(i,j,34) !--- zorll (zorl on land portion of a cell)
           end if
           !
           !--- NSSTM variables
           if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 1)) then
             !--- nsstm tref
             Sfcprop(nb)%tref(ix)    = Sfcprop(nb)%tsfco(ix)
             Sfcprop(nb)%xz(ix)      = 30.0d0
           endif
           if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 0)) then
             Sfcprop(nb)%tref(ix)    = sfc_var2(i,j,nvar_s2m+1)  !--- nsstm tref
             Sfcprop(nb)%z_c(ix)     = sfc_var2(i,j,nvar_s2m+2)  !--- nsstm z_c
             Sfcprop(nb)%c_0(ix)     = sfc_var2(i,j,nvar_s2m+3)  !--- nsstm c_0
             Sfcprop(nb)%c_d(ix)     = sfc_var2(i,j,nvar_s2m+4)  !--- nsstm c_d
             Sfcprop(nb)%w_0(ix)     = sfc_var2(i,j,nvar_s2m+5)  !--- nsstm w_0
             Sfcprop(nb)%w_d(ix)     = sfc_var2(i,j,nvar_s2m+6)  !--- nsstm w_d
             Sfcprop(nb)%xt(ix)      = sfc_var2(i,j,nvar_s2m+7)  !--- nsstm xt
             Sfcprop(nb)%xs(ix)      = sfc_var2(i,j,nvar_s2m+8)  !--- nsstm xs
             Sfcprop(nb)%xu(ix)      = sfc_var2(i,j,nvar_s2m+9)  !--- nsstm xu
             Sfcprop(nb)%xv(ix)      = sfc_var2(i,j,nvar_s2m+10) !--- nsstm xv
             Sfcprop(nb)%xz(ix)      = sfc_var2(i,j,nvar_s2m+11) !--- nsstm xz
             Sfcprop(nb)%zm(ix)      = sfc_var2(i,j,nvar_s2m+12) !--- nsstm zm
             Sfcprop(nb)%xtts(ix)    = sfc_var2(i,j,nvar_s2m+13) !--- nsstm xtts
             Sfcprop(nb)%xzts(ix)    = sfc_var2(i,j,nvar_s2m+14) !--- nsstm xzts
             Sfcprop(nb)%d_conv(ix)  = sfc_var2(i,j,nvar_s2m+15) !--- nsstm d_conv
             Sfcprop(nb)%ifd(ix)     = sfc_var2(i,j,nvar_s2m+16) !--- nsstm ifd
             Sfcprop(nb)%dt_cool(ix) = sfc_var2(i,j,nvar_s2m+17) !--- nsstm dt_cool
             Sfcprop(nb)%qrain(ix)   = sfc_var2(i,j,nvar_s2m+18) !--- nsstm qrain
           endif

! Noah MP
! -------
           if (Model%lsm == Model%lsm_noahmp) then
             Sfcprop(nb)%snowxy(ix)     = sfc_var2(i,j,nvar_s2m+19)
             Sfcprop(nb)%tvxy(ix)       = sfc_var2(i,j,nvar_s2m+20)
             Sfcprop(nb)%tgxy(ix)       = sfc_var2(i,j,nvar_s2m+21)
             Sfcprop(nb)%canicexy(ix)   = sfc_var2(i,j,nvar_s2m+22)
             Sfcprop(nb)%canliqxy(ix)   = sfc_var2(i,j,nvar_s2m+23)
             Sfcprop(nb)%eahxy(ix)      = sfc_var2(i,j,nvar_s2m+24)
             Sfcprop(nb)%tahxy(ix)      = sfc_var2(i,j,nvar_s2m+25)
             Sfcprop(nb)%cmxy(ix)       = sfc_var2(i,j,nvar_s2m+26)
             Sfcprop(nb)%chxy(ix)       = sfc_var2(i,j,nvar_s2m+27)
             Sfcprop(nb)%fwetxy(ix)     = sfc_var2(i,j,nvar_s2m+28)
             Sfcprop(nb)%sneqvoxy(ix)   = sfc_var2(i,j,nvar_s2m+29)
             Sfcprop(nb)%alboldxy(ix)   = sfc_var2(i,j,nvar_s2m+30)
             Sfcprop(nb)%qsnowxy(ix)    = sfc_var2(i,j,nvar_s2m+31)
             Sfcprop(nb)%wslakexy(ix)   = sfc_var2(i,j,nvar_s2m+32)
             Sfcprop(nb)%zwtxy(ix)      = sfc_var2(i,j,nvar_s2m+33)
             Sfcprop(nb)%waxy(ix)       = sfc_var2(i,j,nvar_s2m+34)
             Sfcprop(nb)%wtxy(ix)       = sfc_var2(i,j,nvar_s2m+35)
             Sfcprop(nb)%lfmassxy(ix)   = sfc_var2(i,j,nvar_s2m+36)
             Sfcprop(nb)%rtmassxy(ix)   = sfc_var2(i,j,nvar_s2m+37)
             Sfcprop(nb)%stmassxy(ix)   = sfc_var2(i,j,nvar_s2m+38)
             Sfcprop(nb)%woodxy(ix)     = sfc_var2(i,j,nvar_s2m+39)
             Sfcprop(nb)%stblcpxy(ix)   = sfc_var2(i,j,nvar_s2m+40)
             Sfcprop(nb)%fastcpxy(ix)   = sfc_var2(i,j,nvar_s2m+41)
             Sfcprop(nb)%xsaixy(ix)     = sfc_var2(i,j,nvar_s2m+42)
             Sfcprop(nb)%xlaixy(ix)     = sfc_var2(i,j,nvar_s2m+43)
             Sfcprop(nb)%taussxy(ix)    = sfc_var2(i,j,nvar_s2m+44)
             Sfcprop(nb)%smcwtdxy(ix)   = sfc_var2(i,j,nvar_s2m+45)
             Sfcprop(nb)%deeprechxy(ix) = sfc_var2(i,j,nvar_s2m+46)
             Sfcprop(nb)%rechxy(ix)     = sfc_var2(i,j,nvar_s2m+47)
             Sfcprop(nb)%drainncprv(ix) = sfc_var2(i,j,nvar_s2m+48)
             Sfcprop(nb)%draincprv(ix)  = sfc_var2(i,j,nvar_s2m+49)
             Sfcprop(nb)%dsnowprv(ix)   = sfc_var2(i,j,nvar_s2m+50)
             Sfcprop(nb)%dgraupelprv(ix)= sfc_var2(i,j,nvar_s2m+51)
             Sfcprop(nb)%diceprv(ix)    = sfc_var2(i,j,nvar_s2m+52)
             Sfcprop(nb)%albdvis(ix)    = sfc_var2(i,j,nvar_s2m+53)
             Sfcprop(nb)%albdnir(ix)    = sfc_var2(i,j,nvar_s2m+54)
             Sfcprop(nb)%albivis(ix)    = sfc_var2(i,j,nvar_s2m+55)
             Sfcprop(nb)%albinir(ix)    = sfc_var2(i,j,nvar_s2m+56)
             Sfcprop(nb)%emiss(ix)      = sfc_var2(i,j,nvar_s2m+57)
             Sfcprop(nb)%scolor(ix)     = sfc_var2(i,j,nvar_s2m+58)
           endif


           !--- 3D variables
           do lsoil = 1,Model%lsoil
             Sfcprop(nb)%stc(ix,lsoil) = sfc_var3(i,j,lsoil,1)   !--- stc
             Sfcprop(nb)%smc(ix,lsoil) = sfc_var3(i,j,lsoil,2)   !--- smc
             Sfcprop(nb)%slc(ix,lsoil) = sfc_var3(i,j,lsoil,3)   !--- slc
           enddo

           if (Model%lsm == Model%lsm_noahmp) then
             do lsoil = -2, 0
               Sfcprop(nb)%snicexy(ix,lsoil) = sfc_var3sn(i,j,lsoil,4)
               Sfcprop(nb)%snliqxy(ix,lsoil) = sfc_var3sn(i,j,lsoil,5)
               Sfcprop(nb)%tsnoxy(ix,lsoil)  = sfc_var3sn(i,j,lsoil,6)
             enddo

             do lsoil = 1, 4
               Sfcprop(nb)%smoiseq(ix,lsoil)  = sfc_var3eq(i,j,lsoil,7)
             enddo

             do lsoil = -2, 4
               Sfcprop(nb)%zsnsoxy(ix,lsoil)  = sfc_var3zn(i,j,lsoil,8)
             enddo
           endif

        enddo   !ix
      enddo    !nb

      if (Model%use_ifs_ini_sst) deallocate (ifsSST)

      call mpp_error(NOTE, 'gfs_driver:: - after put to container ')
! so far: At cold start everything is 9999.0, warm start snowxy has values
!         but the 3D of snow fields are not available because not allocated yet.
!         ix,nb loops may be consolidate with the Noah MP isnowxy init
!         restore traditional vars first,we need some of them to init snow fields
!         snow depth to actual snow layers; so we can allocate and register
!         note zsnsoxy is from -2:4 - isnowxy is from 0:-2, but we need
!         exact snow layers to pass 3D fields correctly, snow layers are
!         different fro grid to grid, we have to init point by point/grid.
!         It has to be done after the weasd is available
!         sfc_var2(1,1,32) is the first; we need this to allocate snow related fields

      !--- if sncovr does not exist in the restart, need to create it
      if (nint(sfc_var2(1,1,32)) == -9999) then
        if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver::surface_props_input - computing sncovr')
        !--- compute sncovr from existing variables
        !--- code taken directly from read_fix.f
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            Sfcprop(nb)%sncovr(ix) = 0.0
            if (Sfcprop(nb)%slmsk(ix) > 0.001) then
              vegtyp = Sfcprop(nb)%vtype(ix)
              if (vegtyp == 0) vegtyp = 7
              rsnow  = 0.001*Sfcprop(nb)%weasd(ix)/snupx(vegtyp)
              if (0.001*Sfcprop(nb)%weasd(ix) < snupx(vegtyp)) then
                Sfcprop(nb)%sncovr(ix) = 1.0 - (exp(-salp_data*rsnow) - rsnow*exp(-salp_data))
              else
                Sfcprop(nb)%sncovr(ix) = 1.0
              endif
            endif
          enddo
        enddo
      endif

      if(Model%frac_grid) then ! 3-way composite
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            tem = 1.0 - Sfcprop(nb)%landfrac(ix) - Sfcprop(nb)%fice(ix)
            Sfcprop(nb)%zorl(ix) = Sfcprop(nb)%zorll(ix) * Sfcprop(nb)%landfrac(ix) &
                                 + Sfcprop(nb)%zorll(ix) * Sfcprop(nb)%fice(ix)     & !zorl ice = zorl land
                                 + Sfcprop(nb)%zorlo(ix) * tem
            Sfcprop(nb)%tsfc(ix) = Sfcprop(nb)%tsfcl(ix) * Sfcprop(nb)%landfrac(ix) &
                                 + Sfcprop(nb)%tisfc(ix) * Sfcprop(nb)%fice(ix)     &
                                 + Sfcprop(nb)%tsfco(ix) * tem
          enddo
        enddo
      else     ! in this case ice fracion is fraction of water fraction
        do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            !--- specify tsfcl/zorll from existing variable tsfco/zorlo
            Sfcprop(nb)%tsfcl(ix) = Sfcprop(nb)%tsfco(ix)
            Sfcprop(nb)%zorll(ix) = Sfcprop(nb)%zorlo(ix)
            Sfcprop(nb)%zorl(ix)  = Sfcprop(nb)%zorlo(ix)
            Sfcprop(nb)%tsfc(ix)  = Sfcprop(nb)%tsfco(ix)
            if (Sfcprop(nb)%slmsk(ix) > 1.9) then
              Sfcprop(nb)%landfrac(ix) = 0.0
            else
              Sfcprop(nb)%landfrac(ix) = Sfcprop(nb)%slmsk(ix)
            endif
          enddo
        enddo
      endif ! if (Model%frac_grid)

      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          if (Sfcprop(nb)%lakefrac(ix) > 0.0) then
            Sfcprop(nb)%oceanfrac(ix) = 0.0 ! lake & ocean don't coexist in a cell
          else
            Sfcprop(nb)%oceanfrac(ix) = 1.0 - Sfcprop(nb)%landfrac(ix)  !LHS:ocean frac [0:1]
          endif

        enddo
      enddo

      if (Model%lsm == Model%lsm_noahmp) then
        if (sfc_var2(1,1,nvar_s2m+19) < -9990. ) then
          !--- initialize NOAH MP properties
          if (Model%me == Model%master ) call mpp_error(NOTE, 'gfs_driver:: - Cold start Noah MP ')

          do nb = 1, Atm_block%nblks
            do ix = 1, Atm_block%blksz(nb)

              Sfcprop(nb)%tvxy(ix)     = missing_value
              Sfcprop(nb)%tgxy(ix)     = missing_value
              Sfcprop(nb)%tahxy(ix)    = missing_value
              Sfcprop(nb)%canicexy(ix) = missing_value
              Sfcprop(nb)%canliqxy(ix) = missing_value
              Sfcprop(nb)%eahxy(ix)    = missing_value
              Sfcprop(nb)%cmxy(ix)     = missing_value
              Sfcprop(nb)%chxy(ix)     = missing_value
              Sfcprop(nb)%fwetxy(ix)   = missing_value
              Sfcprop(nb)%sneqvoxy(ix) = missing_value
              Sfcprop(nb)%alboldxy(ix) = missing_value
              Sfcprop(nb)%qsnowxy(ix)  = missing_value
              Sfcprop(nb)%wslakexy(ix) = missing_value
              Sfcprop(nb)%taussxy(ix)  = missing_value
              Sfcprop(nb)%waxy(ix)     = missing_value
              Sfcprop(nb)%wtxy(ix)     = missing_value
              Sfcprop(nb)%zwtxy(ix)    = missing_value
              Sfcprop(nb)%xlaixy(ix)   = missing_value
              Sfcprop(nb)%xsaixy(ix)   = missing_value

              Sfcprop(nb)%lfmassxy(ix) = missing_value
              Sfcprop(nb)%stmassxy(ix) = missing_value
              Sfcprop(nb)%rtmassxy(ix) = missing_value
              Sfcprop(nb)%woodxy(ix)   = missing_value
              Sfcprop(nb)%stblcpxy(ix) = missing_value
              Sfcprop(nb)%fastcpxy(ix) = missing_value
              Sfcprop(nb)%smcwtdxy(ix) = missing_value
              Sfcprop(nb)%deeprechxy(ix) = missing_value
              Sfcprop(nb)%rechxy(ix)     = missing_value
              Sfcprop(nb)%drainncprv(ix) = 0.
              Sfcprop(nb)%draincprv(ix)  = 0.
              Sfcprop(nb)%dsnowprv(ix)   = 0.
              Sfcprop(nb)%dgraupelprv(ix)= 0.
              Sfcprop(nb)%diceprv(ix)    = 0.
              Sfcprop(nb)%albdvis(ix)    = missing_value
              Sfcprop(nb)%albdnir(ix)    = missing_value
              Sfcprop(nb)%albivis(ix)    = missing_value
              Sfcprop(nb)%albinir(ix)    = missing_value
              Sfcprop(nb)%emiss(ix)      = missing_value
              Sfcprop(nb)%scolor(ix)     = 0.0

              Sfcprop(nb)%snowxy (ix)   = missing_value
              Sfcprop(nb)%snicexy(ix, -2:0) = missing_value
              Sfcprop(nb)%snliqxy(ix, -2:0) = missing_value
              Sfcprop(nb)%tsnoxy (ix, -2:0) = missing_value
              Sfcprop(nb)%smoiseq(ix,  1:4) = missing_value
              Sfcprop(nb)%zsnsoxy(ix, -2:4) = missing_value

              if (Sfcprop(nb)%slmsk(ix) > 0.01) then

                Sfcprop(nb)%tvxy(ix)     = Sfcprop(nb)%tsfcl(ix)
                Sfcprop(nb)%tgxy(ix)     = Sfcprop(nb)%tsfcl(ix)
                Sfcprop(nb)%tahxy(ix)    = Sfcprop(nb)%tsfcl(ix)

                if (Sfcprop(nb)%snowd(ix) > 0.01 .and. Sfcprop(nb)%tsfcl(ix) > 273.15 ) Sfcprop(nb)%tvxy(ix)  = 273.15
                if (Sfcprop(nb)%snowd(ix) > 0.01 .and. Sfcprop(nb)%tsfcl(ix) > 273.15 ) Sfcprop(nb)%tgxy(ix)  = 273.15
                if (Sfcprop(nb)%snowd(ix) > 0.01 .and. Sfcprop(nb)%tsfcl(ix) > 273.15 ) Sfcprop(nb)%tahxy(ix) = 273.15

                Sfcprop(nb)%canicexy(ix) = 0.0
                Sfcprop(nb)%canliqxy(ix) = Sfcprop(nb)%canopy(ix)

                Sfcprop(nb)%eahxy(ix)    = 2000.0

      !---  eahxy = psfc*qv/(0.622+qv); qv is mixing ratio, converted from sepcific
      !     humidity specific humidity /(1.0 - specific humidity)

                Sfcprop(nb)%cmxy(ix)     = 0.0
                Sfcprop(nb)%chxy(ix)     = 0.0
                Sfcprop(nb)%fwetxy(ix)   = 0.0
                Sfcprop(nb)%sneqvoxy(ix) = Sfcprop(nb)%weasd(ix)     ! mm
                Sfcprop(nb)%alboldxy(ix) = 0.65
                Sfcprop(nb)%qsnowxy(ix)  = 0.0

      !---  if (Sfcprop(nb)%srflag(ix) > 0.001) Sfcprop(nb)%qsnowxy(ix) = Sfcprop(nb)%tprcp(ix)/Model%dtp
      !     already set to 0.0
                Sfcprop(nb)%wslakexy(ix) = 0.0
                Sfcprop(nb)%taussxy(ix)  = 0.0

                Sfcprop(nb)%albdvis(ix)  = 0.2
                Sfcprop(nb)%albdnir(ix)  = 0.2
                Sfcprop(nb)%albivis(ix)  = 0.2
                Sfcprop(nb)%albinir(ix)  = 0.2

                Sfcprop(nb)%waxy(ix)     = 4900.0
                Sfcprop(nb)%wtxy(ix)     = Sfcprop(nb)%waxy(ix)
                Sfcprop(nb)%zwtxy(ix)    = (25.0 + 2.0) - Sfcprop(nb)%waxy(ix) / 1000.0 /0.2
                !
                vegtyp                   = Sfcprop(nb)%vtype(ix)
                if (vegtyp == 0) vegtyp = 7
                imn                      = Model%idate(2)

                if ((vegtyp == isbarren_table) .or. (vegtyp == isice_table) .or. (vegtyp == isurban_table) .or. &
                  & (vegtyp == iswater_table)) then

                  Sfcprop(nb)%xlaixy(ix)   = 0.0
                  Sfcprop(nb)%xsaixy(ix)   = 0.0

                  Sfcprop(nb)%lfmassxy(ix) = 0.0
                  Sfcprop(nb)%stmassxy(ix) = 0.0
                  Sfcprop(nb)%rtmassxy(ix) = 0.0

                  Sfcprop(nb)%woodxy   (ix) = 0.0
                  Sfcprop(nb)%stblcpxy (ix) = 0.0
                  Sfcprop(nb)%fastcpxy (ix) = 0.0

                else


                  Sfcprop(nb)%xlaixy(ix)   = max(laim_table(vegtyp, imn),0.05)
                  Sfcprop(nb)%xsaixy(ix)   = max(Sfcprop(nb)%xlaixy(ix)*0.1,0.05)

                  masslai                  = 1000.0 / max(sla_table(vegtyp),1.0)
                  Sfcprop(nb)%lfmassxy(ix) = Sfcprop(nb)%xlaixy(ix)*masslai
                  masssai                  = 1000.0 / 3.0
                  Sfcprop(nb)%stmassxy(ix) = Sfcprop(nb)%xsaixy(ix)* masssai

                  Sfcprop(nb)%rtmassxy(ix) = 500.0

                  Sfcprop(nb)%woodxy  (ix) = 500.0
                  Sfcprop(nb)%stblcpxy(ix) = 1000.0
                  Sfcprop(nb)%fastcpxy(ix) = 1000.0

                endif  ! non urban ...

                emv = 1. - exp(-(Sfcprop(nb)%xlaixy(ix)+Sfcprop(nb)%xsaixy(ix))/1.0) ! Ignore snow-buried sai and lai during the initialization
                emg = 0.97*(1.-Sfcprop(nb)%sncovr(ix)) + 1.0*Sfcprop(nb)%sncovr(ix)
                Sfcprop(nb)%emiss(ix) = Sfcprop(nb)%vfrac(ix) * ( emg*(1-emv) + emv + emv*(1-emv)*(1-emg) ) + (1-Sfcprop(nb)%vfrac(ix)) * emg

                if ( vegtyp == isice_table )  then
                  do lsoil = 1,Model%lsoil
                    Sfcprop(nb)%stc(ix,lsoil) = min(Sfcprop(nb)%stc(ix,lsoil),min(Sfcprop(nb)%tg3(ix),263.15))
                    Sfcprop(nb)%smc(ix,lsoil) = 1
                    Sfcprop(nb)%slc(ix,lsoil) = 0
                  enddo
                endif

                snd   = Sfcprop(nb)%snowd(ix)/1000.0  ! go to m from snwdph

                if (Sfcprop(nb)%weasd(ix) /= 0.0 .and. snd == 0.0 ) then
                  snd = Sfcprop(nb)%weasd(ix)*0.005
                  Sfcprop(nb)%snowd(ix) = snd*1000.0
                endif

                ! cap SNOW at 2000, maintain density
                if (Sfcprop(nb)%weasd(ix) > 2000.0 ) then
                  snd = snd * 2000.0 / Sfcprop(nb)%weasd(ix)
                  Sfcprop(nb)%weasd(ix) = 2000.0
                  Sfcprop(nb)%snowd(ix) = snd*1000.0
                endif

                if (vegtyp == 15) then                      ! land ice in MODIS/IGBP
                  if ( Sfcprop(nb)%weasd(ix) < 0.1) then
                    Sfcprop(nb)%weasd(ix) = 0.1
                    snd                   = 0.01
                  endif
                endif

                if (snd < 0.025 ) then
                  Sfcprop(nb)%snowxy(ix)   = 0.0
                  dzsno(-2:0)              = 0.0
                elseif (snd >= 0.025 .and. snd <= 0.05 ) then
                  Sfcprop(nb)%snowxy(ix)   = -1.0
                  dzsno(0)                 = snd
                elseif (snd > 0.05 .and. snd <= 0.10 ) then
                  Sfcprop(nb)%snowxy(ix)   = -2.0
                  dzsno(-1)                = 0.5*snd
                  dzsno(0)                 = 0.5*snd
                elseif (snd > 0.10 .and. snd <= 0.25 ) then
                  Sfcprop(nb)%snowxy(ix)   = -2.0
                  dzsno(-1)                = 0.05
                  dzsno(0)                 = snd - 0.05
                elseif (snd > 0.25 .and. snd <= 0.45 ) then
                  Sfcprop(nb)%snowxy(ix)   = -3.0
                  dzsno(-2)                = 0.05
                  dzsno(-1)                = 0.5*(snd-0.05)
                  dzsno(0)                 = 0.5*(snd-0.05)
                elseif (snd > 0.45) then
                  Sfcprop(nb)%snowxy(ix)   = -3.0
                  dzsno(-2)                = 0.05
                  dzsno(-1)                = 0.20
                  dzsno(0)                 = snd - 0.05 - 0.20
                else
                  call mpp_error(FATAL, 'problem with the logic assigning snow layers.')
                endif

      !---  Now we have the snowxy field
      !     snice + snliq + tsno allocation and compute them from what we have
                Sfcprop(nb)%tsnoxy(ix,-2:0)  = 0.0
                Sfcprop(nb)%snicexy(ix,-2:0) = 0.0
                Sfcprop(nb)%snliqxy(ix,-2:0) = 0.0
                Sfcprop(nb)%zsnsoxy(ix,-2:4) = 0.0

                isnow = nint(Sfcprop(nb)%snowxy(ix))+1    ! snowxy <=0.0, dzsno >= 0.0

                do ns = isnow , 0
                  Sfcprop(nb)%tsnoxy(ix,ns)  = Sfcprop(nb)%tgxy(ix)
                  Sfcprop(nb)%snliqxy(ix,ns) = 0.0
                  Sfcprop(nb)%snicexy(ix,ns) = 1.00 * dzsno(ns) * Sfcprop(nb)%weasd(ix)/snd
                enddo
      !
      !--- zsnsoxy, all negative ?
      !
                do ns = isnow, 0
                  dzsnso(ns) = -dzsno(ns)
                enddo

                do ns = 1 , 4
                  dzsnso(ns) = -dzs(ns)
                enddo
      !
      !--- Assign to zsnsoxy
      !
                Sfcprop(nb)%zsnsoxy(ix,isnow) = dzsnso(isnow)
                do ns = isnow+1,4
                  Sfcprop(nb)%zsnsoxy(ix,ns) = Sfcprop(nb)%zsnsoxy(ix,ns-1) + dzsnso(ns)
                enddo

      !
      !--- smoiseq
      !    Init water table related quantities here
      !
                soiltyp  = Sfcprop(nb)%stype(ix)


                if (soiltyp == 0) then
                  Sfcprop(nb)%stype(ix) = 16
                  soiltyp  = Sfcprop(nb)%stype(ix)
                endif

                bexp   = bexp_table(soiltyp)
                smcmax = smcmax_table(soiltyp)
                smcwlt = smcwlt_table(soiltyp)
                dwsat  = dwsat_table(soiltyp)
                dksat  = dksat_table(soiltyp)
                psisat = -psisat_table(soiltyp)

                if (vegtyp == isurban_table) then
                  smcmax = 0.45
                  smcwlt = 0.40
                endif

                if ((bexp > 0.0) .and. (smcmax > 0.0) .and. (-psisat > 0.0 )) then
                  do ns = 1, Model%lsoil
                    if ( ns == 1 )then
                      ddz = -zsoil(ns+1) * 0.5
                    elseif ( ns < Model%lsoil ) then
                      ddz = ( zsoil(ns-1) - zsoil(ns+1) ) * 0.5
                    else
                      ddz = zsoil(ns-1) - zsoil(ns)
                    endif
      !
      !--- Use newton-raphson method to find eq soil moisture
      !
                    expon = bexp + 1.
                    aa    = dwsat / ddz
                    bb    = dksat / smcmax ** expon

                    smc = 0.5 * smcmax

                    do iter = 1, 100
                      func  = (smc - smcmax) * aa +  bb * smc ** expon
                      dfunc = aa + bb * expon * smc ** bexp
                      dx    = func / dfunc
                      smc   = smc - dx
                      if ( abs (dx) < 1.e-6) exit
                    enddo                               ! iteration
                    Sfcprop(nb)%smoiseq(ix,ns) = min(max(smc,1.e-4),smcmax*0.99)
                  enddo                                 ! ddz soil layer
                else                                    ! bexp <= 0.0
                  Sfcprop(nb)%smoiseq(ix,1:4) = smcmax
                endif                                   ! end the bexp condition

                Sfcprop(nb)%smcwtdxy(ix)   = smcmax
                Sfcprop(nb)%deeprechxy(ix) = 0.0
                Sfcprop(nb)%rechxy(ix)     = 0.0

                ! Use a default value of 4 for the soil color category over
                ! land when cold starting. Note this will get overridden during
                ! cycling if soil color data is provided. If soil color data is
                ! not provided then the soil color will remain 4 over land.
                Sfcprop(nb)%scolor(ix)     = 4.0
              endif !end if slmsk>0.01 (land only)

            enddo ! ix
          enddo  ! nb
        endif
      endif !if Noah MP cold start ends

    else   !--- ELSE of IF (open_file(fn_srf) ...

      !--- Noah MP define arbitrary value (number layers of snow) to indicate
      !    coldstart(sfcfile doesn't include noah mp fields) or not

      call mpp_error(NOTE,'No INPUT/sfc_data.tile*.nc surface data found; cold-starting land surface')
      !Need a namelist for options:
      ! 1. choice of sst (uniform, profiles) --- ML0 should relax to this
      ! 2. Choice of veg, soil type with certain soil T,q,ql
      ! How to fix day of year (for astronomy)?
      !--- place the data into the block GFS containers
      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          !--- 2D variables
          !--- slmsk
          Sfcprop(nb)%slmsk(ix)  = 0.
          !--- tsfc (tsea in sfc file)
          Sfcprop(nb)%tsfc(ix)   = Model%Ts0 ! should specify some latitudinal profile
          !--- weasd (sheleg in sfc file)
          Sfcprop(nb)%weasd(ix)  = 0.0
          !--- tg3
          Sfcprop(nb)%tg3(ix)    = 290. !generic value, probably not good; real value latitude-dependent
          !--- zorl
          Sfcprop(nb)%zorl(ix)   = 0.1 ! typical ocean value; different values for different land surfaces (use a lookup table?)
          !--- alvsf
          Sfcprop(nb)%alvsf(ix)  = 0.06
          !--- alvwf
          Sfcprop(nb)%alvwf(ix)  = 0.06
          !--- alnsf
          Sfcprop(nb)%alnsf(ix)  = 0.06
          !--- alnwf
          Sfcprop(nb)%alnwf(ix)  = 0.06
          !--- facsf
          Sfcprop(nb)%facsf(ix)  = 0.0
          !--- facwf
          Sfcprop(nb)%facwf(ix)  = 0.0
          !--- vfrac
          Sfcprop(nb)%vfrac(ix)  = 0.0
          !--- canopy
          Sfcprop(nb)%canopy(ix) = 0.0
          !--- f10m
          Sfcprop(nb)%f10m(ix)   = 0.9
          !--- t2m
          Sfcprop(nb)%t2m(ix)    = Sfcprop(nb)%tsfc(ix)
          !--- q2m
          Sfcprop(nb)%q2m(ix)    = 0.0 ! initially dry atmosphere?
          !--- vtype
          Sfcprop(nb)%vtype(ix)  = 0.0
          !--- stype
          Sfcprop(nb)%stype(ix)  = 0.0
          !--- uustar
          Sfcprop(nb)%uustar(ix) = 0.5
          !--- ffmm
          Sfcprop(nb)%ffmm(ix)   = 10.
          !--- ffhh
          Sfcprop(nb)%ffhh(ix)   = 10.
          !--- hice
          Sfcprop(nb)%hice(ix)   = 0.0
          !--- fice
          Sfcprop(nb)%fice(ix)   = 0.0
          !--- tisfc
          Sfcprop(nb)%tisfc(ix)  = Sfcprop(nb)%tsfc(ix)
          !--- tprcp
          Sfcprop(nb)%tprcp(ix)  = 0.0
          !--- srflag
          Sfcprop(nb)%srflag(ix) = 0.0
          !--- snowd (snwdph in the file)
          Sfcprop(nb)%snowd(ix)  = 0.0
          !--- shdmin
          Sfcprop(nb)%shdmin(ix) = 0.0 !this and the next depend on the surface type
          !--- shdmax
          Sfcprop(nb)%shdmax(ix) = 0.0
          !--- slope
          Sfcprop(nb)%slope(ix)  = 0.0 ! also land-surface dependent
          !--- snoalb
          Sfcprop(nb)%snoalb(ix) = 0.0
          !--- sncovr
          Sfcprop(nb)%sncovr(ix) = 0.0
          !
          !--- NSSTM variables
          if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 1)) then
             !--- nsstm tref
             Sfcprop(nb)%tref(ix)    = Sfcprop(nb)%tsfc(ix)
             Sfcprop(nb)%xz(ix)      = 30.0d0
          endif
          if ((Model%nstf_name(1) > 0) .and. (Model%nstf_name(2) == 0)) then
             !return an error
             call mpp_error(FATAL, 'cold-starting does not support NSST.')
          endif

          !--- 3D variables
          ! these are all set to ocean values.
          !--- stc
          Sfcprop(nb)%stc(ix,:) = Sfcprop(nb)%tsfc(ix)
          !--- smc
          Sfcprop(nb)%smc(ix,:) = 1.0
          !--- slc
          Sfcprop(nb)%slc(ix,:) = 1.0
        enddo
      enddo
      !--- end of file not existing cold-start case

    endif  !--- END of IF (open_file(fn_srf) ...


  end subroutine sfc_prop_restart_read

  subroutine compute_surface_type_fraction(Atm_block, isc, jsc, nx, ny, diagnostic, result)
    type(block_control_type), intent(in) :: Atm_block
    integer, intent(in) :: isc, jsc, nx, ny
    type(gfdl_diag_type), intent(in) :: diagnostic
    real(kind=kind_phys), intent(out) :: result(nx, ny)

    integer :: i, j, ii, jj, nb, ix, surface_type_code

    select case(trim(diagnostic%name))
      case ('ocean_fraction')
        surface_type_code = 0
      case ('land_fraction')
        surface_type_code = 1
      case ('sea_ice_fraction')
        surface_type_code = 2
    end select

    do j = 1, ny
      jj = j + jsc - 1
      do i = 1, nx
          ii = i + isc - 1
          nb = Atm_block%blkno(ii,jj)
          ix = Atm_block%ixp(ii,jj)
          if (nint(diagnostic%data(nb)%var2(ix)) .eq. surface_type_code) then
            result(i,j) = 1.0
          else
            result(i,j) = 0.0
          endif
      enddo
    enddo
  end subroutine compute_surface_type_fraction

  subroutine sfc_data_override(Time, IPD_data, Atm_block, Model)

    implicit none
    !--- interface variable definitions
    type(time_type),           intent(in)    :: Time
    type(IPD_data_type),       intent(inout) :: IPD_Data(:)
    type (block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),    intent(in)    :: Model
    !--- local variables
    integer :: i, j, ix, nb
    integer :: isc, iec, jsc, jec

    logical :: used
    real, allocatable :: sst(:,:), ci(:,:)

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec

    if (Model%use_ext_sst) then

        ! Here is a sample data_table that will enable reading in
        ! external SSTs and sea-ice from an external file.
        !
        !"ATM", "sst", "sst", "INPUT/ec_sst.nc", "bilinear", 1.0
        !"ATM", "ci", "ci", "INPUT/ec_sst.nc", "bilinear", 1.0

        allocate(sst(isc:iec,jsc:jec))
        allocate(ci(isc:iec,jsc:jec))
        call data_override('ATM', 'sst', sst, Time, override=used)
        if (.not. used) then
           call mpp_error(FATAL, " SST dataset not specified in data_table.")
        endif
        call data_override('ATM', 'ci', ci, Time, override=used)
        if (.not. used) then
           call mpp_error(NOTE, " Sea ice fraction dataset not specified in data_table. No override will occur.")
           ci(:,:) = -999.
        endif
        do nb = 1, Atm_block%nblks
           do ix = 1, Atm_block%blksz(nb)
              i = Atm_block%index(nb)%ii(ix)
              j = Atm_block%index(nb)%jj(ix)
              IPD_Data(nb)%Statein%sst(ix) = sst(i,j)
              IPD_Data(nb)%Statein%ci(ix) = ci(i,j)
           enddo
        enddo
        deallocate(sst)
        deallocate(ci)

    endif

  end subroutine sfc_data_override

  subroutine sfc_prop_override(Sfcprop, Grid, Atm_block, Model, fv_domain)

    implicit none
    !--- interface variable definitions
    type(GFS_sfcprop_type),    intent(inout) :: Sfcprop(:)
    type(GFS_grid_type),       intent(inout) :: Grid(:)
    type (block_control_type), intent(in)    :: Atm_block
    type(IPD_control_type),    intent(in)    :: Model
    type (domain2d),           intent(in)    :: fv_domain
    !--- local variables
    integer :: i, j, k, ix, lsoil, num, nb
    integer :: isc, iec, jsc, jec, npz, nx, ny, ios

    logical :: ideal_sst = .false.
    real(kind=kind_phys) :: sst_max = 300.
    real(kind=kind_phys) :: sst_min = 271.14 ! -2c --> sea ice
    integer :: sst_profile = 0

    logical :: ideal_land = .false.
    !Assuming modern veg/soil types
    ! sample Amazon settings; values for OKC in comments
    integer :: vegtype = 2 ! 12
    integer :: soiltype = 9 ! 8
    real(kind=kind_phys) :: vegfrac = 0.8 ! 0.25 -- 0.75
    real(kind=kind_phys) :: zorl = 265 ! 15
    !uniform soil temperature and moisture for now
    real(kind=kind_phys) :: stc = 300. ! 310.
    real(kind=kind_phys) :: smc = 0.4 ! wet season vs. 0.08 dry ! 0.2 okc highly variable and patchy

    namelist /sfc_prop_override_nml/ &
         ideal_sst, sst_max, sst_profile, & !Aquaplanet SST options
         ideal_land, vegtype, soiltype, & ! idealized soil/veg options
         vegfrac, zorl, stc, smc

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)

    !--- read the sfc_prop_override namelist
    read(Model%input_nml_file, nml=sfc_prop_override_nml, iostat=ios)

    call qs_init

    call mpp_error(NOTE, "Calling sfc_prop_override")

    if (ideal_sst) then
       do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
             i = Atm_block%index(nb)%ii(ix) - isc + 1
             j = Atm_block%index(nb)%jj(ix) - jsc + 1
             !--- slmsk
             Sfcprop(nb)%slmsk(ix)  = 0.0
             !--- tsfc (tsea in sfc file)
             select case (sst_profile)
                case (0)
                   Sfcprop(nb)%tsfc(ix)   = sst_max
                case (1) ! symmetric
                   Sfcprop(nb)%tsfc(ix)   = sst_min + (sst_max - sst_min)*Grid(nb)%coslat(ix)
                case default
                   call mpp_error(FATAL, "value of sst_profile not defined.")
             end select
             !--- zorl
             Sfcprop(nb)%zorl(ix)   = zorl
             !--- vfrac
             Sfcprop(nb)%vfrac(ix)  = 0.0
             if (Sfcprop(nb)%tsfc(ix) <= sst_min) then
                !--- hice
                Sfcprop(nb)%hice(ix)   = 1.0
                !--- fice
                Sfcprop(nb)%fice(ix)   = 1.0
                Sfcprop(nb)%tsfc(ix)   = sst_min
             else
                !--- hice
                Sfcprop(nb)%hice(ix)   = 0.0
                !--- fice
                Sfcprop(nb)%fice(ix)   = 0.0
             endif
             !--- tisfc
             Sfcprop(nb)%tisfc(ix)  = Sfcprop(nb)%tsfc(ix)
             !--- t2m ! slt. unstable
             Sfcprop(nb)%t2m(ix)    = Sfcprop(nb)%t2m(ix) * 0.98
             !--- q2m ! use RH = 98% and assume ps = 1000 mb
             Sfcprop(nb)%q2m(ix)    = wqs (Sfcprop(nb)%t2m(ix), 1.e5/rd/Sfcprop(nb)%t2m(ix), Sfcprop(nb)%q2m(ix))
             !--- vtype
             Sfcprop(nb)%vtype(ix)  = 0
             !--- stype
             Sfcprop(nb)%stype(ix)  = 0
             !Override MLO properties also
             if (Model%do_ocean) then
                Sfcprop(nb)%ts_clim_iano(ix) = Sfcprop(nb)%tsfc(ix)
                Sfcprop(nb)%tsclim(ix) = Sfcprop(nb)%tsfc(ix)
                Sfcprop(nb)%ts_som(ix) = Sfcprop(nb)%tsfc(ix)
             endif
          enddo
       enddo


    elseif (ideal_land) then
       do nb = 1, Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
             i = Atm_block%index(nb)%ii(ix) - isc + 1
             j = Atm_block%index(nb)%jj(ix) - jsc + 1
             !--- slmsk
             Sfcprop(nb)%slmsk(ix)  = 1.0
             !--- tsfc (tsea in sfc file)
             Sfcprop(nb)%tsfc(ix)   = stc
             !--- weasd (sheleg in sfc file)
             Sfcprop(nb)%weasd(ix)  = 0.0 ! snow
             !--- tg3
             Sfcprop(nb)%tg3(ix)    = stc ! simple approach
             !--- zorl
             Sfcprop(nb)%zorl(ix)   = zorl
             !--- vfrac
             Sfcprop(nb)%vfrac(ix)  = vegfrac
             !--- canopy
             Sfcprop(nb)%canopy(ix) = 0.0 !this quantity is quite variable
             !--- t2m
             Sfcprop(nb)%t2m(ix)    = stc * 0.98 !slt unstable
             !--- q2m ! use RH = 98%
             Sfcprop(nb)%q2m(ix)    = wqs (Sfcprop(nb)%t2m(ix), 1.e5/rd/Sfcprop(nb)%t2m(ix), Sfcprop(nb)%q2m(ix))
             !--- vtype
             Sfcprop(nb)%vtype(ix)  = vegtype
             !--- stype
             Sfcprop(nb)%stype(ix)  = soiltype
             !--- hice
             Sfcprop(nb)%hice(ix)   = 0.0
             !--- fice
             Sfcprop(nb)%fice(ix)   = 0.0
             !--- tisfc
             Sfcprop(nb)%tisfc(ix)  = stc
             !--- snowd (snwdph in the file)
             Sfcprop(nb)%snowd(ix)  = 0.0
             !--- snoalb
             Sfcprop(nb)%snoalb(ix) = 0.5
             !--- sncovr
             Sfcprop(nb)%sncovr(ix) = 0.0
             !--- 3D variables
             do lsoil = 1,Model%lsoil
                !--- stc
                Sfcprop(nb)%stc(ix,lsoil) = stc
                !--- smc
                Sfcprop(nb)%smc(ix,lsoil) = smc
                !--- slc = smc
                Sfcprop(nb)%slc(ix,lsoil) = smc
             enddo

          enddo
       enddo

    endif


  end subroutine sfc_prop_override


!----------------------------------------------------------------------
! sfc_prop_restart_write
!----------------------------------------------------------------------
!    routine to write out GFS surface restarts via the GFDL FMS restart
!    subsystem.
!    takes an optional argument to append timestamps for intermediate
!    restarts.
!
!    calls:  register_restart_field, save_restart
!----------------------------------------------------------------------
  subroutine sfc_prop_restart_write (Sfcprop, Atm_block, Model, fv_domain, timestamp)
    !--- interface variable definitions
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: Atm_block
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- local variables
    integer :: i, j, k, nb, ix, lsoil, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: nvar_s2m
    character(len=64) :: fname
    character(len=32) :: fn_srf = 'sfc_data.nc'

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)

    !--- Open the restart file and associate it with the Sfc_restart fileobject
    if (present(timestamp)) then
      fname='RESTART/'//trim(timestamp)//'.'//trim(fn_srf)
    else
      fname='RESTART/'//trim(fn_srf)
    endif

    if (open_file(Sfc_restart, fname, "overwrite", fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)) then

      !--- register the variables needed for the restart write
      call register_sfc_prop_restart_vars(Model, nx, ny, nvar_s2m, action="write")

      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          !--- 2D variables
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          sfc_var2(i,j,1)  = Sfcprop(nb)%slmsk(ix) !--- slmsk
          sfc_var2(i,j,2)  = Sfcprop(nb)%tsfc(ix)  !--- tsfc (tsea in sfc file)
          sfc_var2(i,j,3)  = Sfcprop(nb)%weasd(ix) !--- weasd (sheleg in sfc file)
          sfc_var2(i,j,4)  = Sfcprop(nb)%tg3(ix)   !--- tg3
          sfc_var2(i,j,5)  = Sfcprop(nb)%zorl(ix)  !--- zorl
          sfc_var2(i,j,6)  = Sfcprop(nb)%alvsf(ix) !--- alvsf
          sfc_var2(i,j,7)  = Sfcprop(nb)%alvwf(ix) !--- alvwf
          sfc_var2(i,j,8)  = Sfcprop(nb)%alnsf(ix) !--- alnsf
          sfc_var2(i,j,9)  = Sfcprop(nb)%alnwf(ix) !--- alnwf
          sfc_var2(i,j,10) = Sfcprop(nb)%facsf(ix) !--- facsf
          sfc_var2(i,j,11) = Sfcprop(nb)%facwf(ix) !--- facwf
          sfc_var2(i,j,12) = Sfcprop(nb)%vfrac(ix) !--- vfrac
          sfc_var2(i,j,13) = Sfcprop(nb)%canopy(ix)!--- canopy
          sfc_var2(i,j,14) = Sfcprop(nb)%f10m(ix)  !--- f10m
          sfc_var2(i,j,15) = Sfcprop(nb)%t2m(ix)   !--- t2m
          sfc_var2(i,j,16) = Sfcprop(nb)%q2m(ix)   !--- q2m
          sfc_var2(i,j,17) = Sfcprop(nb)%vtype(ix) !--- vtype
          sfc_var2(i,j,18) = Sfcprop(nb)%stype(ix) !--- stype
          sfc_var2(i,j,19) = Sfcprop(nb)%uustar(ix)!--- uustar
          sfc_var2(i,j,20) = Sfcprop(nb)%ffmm(ix)  !--- ffmm
          sfc_var2(i,j,21) = Sfcprop(nb)%ffhh(ix)  !--- ffhh
          sfc_var2(i,j,22) = Sfcprop(nb)%hice(ix)  !--- hice
          sfc_var2(i,j,23) = Sfcprop(nb)%fice(ix)  !--- fice
          sfc_var2(i,j,24) = Sfcprop(nb)%tisfc(ix) !--- tisfc
          sfc_var2(i,j,25) = Sfcprop(nb)%tprcp(ix) !--- tprcp
          sfc_var2(i,j,26) = Sfcprop(nb)%srflag(ix)!--- srflag
          sfc_var2(i,j,27) = Sfcprop(nb)%snowd(ix) !--- snowd (snwdph in the file)
          sfc_var2(i,j,28) = Sfcprop(nb)%shdmin(ix)!--- shdmin
          sfc_var2(i,j,29) = Sfcprop(nb)%shdmax(ix)!--- shdmax
          sfc_var2(i,j,30) = Sfcprop(nb)%slope(ix) !--- slope
          sfc_var2(i,j,31) = Sfcprop(nb)%snoalb(ix)!--- snoalb
          sfc_var2(i,j,32) = Sfcprop(nb)%sncovr(ix)!--- sncovr
          if (Model%cplflx) then
            sfc_var2(i,j,33) = Sfcprop(nb)%tsfcl(ix) !--- tsfcl (temp on land)
            sfc_var2(i,j,34) = Sfcprop(nb)%zorll(ix) !--- zorll (zorl on land)
          end if
          !--- NSSTM variables
          if (Model%nstf_name(1) > 0) then
            sfc_var2(i,j,nvar_s2m+1) = Sfcprop(nb)%tref(ix)    !--- nsstm tref
            sfc_var2(i,j,nvar_s2m+2) = Sfcprop(nb)%z_c(ix)     !--- nsstm z_c
            sfc_var2(i,j,nvar_s2m+3) = Sfcprop(nb)%c_0(ix)     !--- nsstm c_0
            sfc_var2(i,j,nvar_s2m+4) = Sfcprop(nb)%c_d(ix)     !--- nsstm c_d
            sfc_var2(i,j,nvar_s2m+5) = Sfcprop(nb)%w_0(ix)     !--- nsstm w_0
            sfc_var2(i,j,nvar_s2m+6) = Sfcprop(nb)%w_d(ix)     !--- nsstm w_d
            sfc_var2(i,j,nvar_s2m+7) = Sfcprop(nb)%xt(ix)      !--- nsstm xt
            sfc_var2(i,j,nvar_s2m+8) = Sfcprop(nb)%xs(ix)      !--- nsstm xs
            sfc_var2(i,j,nvar_s2m+9) = Sfcprop(nb)%xu(ix)      !--- nsstm xu
            sfc_var2(i,j,nvar_s2m+10) = Sfcprop(nb)%xv(ix)     !--- nsstm xv
            sfc_var2(i,j,nvar_s2m+11) = Sfcprop(nb)%xz(ix)     !--- nsstm xz
            sfc_var2(i,j,nvar_s2m+12) = Sfcprop(nb)%zm(ix)     !--- nsstm zm
            sfc_var2(i,j,nvar_s2m+13) = Sfcprop(nb)%xtts(ix)   !--- nsstm xtts
            sfc_var2(i,j,nvar_s2m+14) = Sfcprop(nb)%xzts(ix)   !--- nsstm xzts
            sfc_var2(i,j,nvar_s2m+15) = Sfcprop(nb)%d_conv(ix) !--- nsstm d_conv
            sfc_var2(i,j,nvar_s2m+16) = Sfcprop(nb)%ifd(ix)    !--- nsstm ifd
            sfc_var2(i,j,nvar_s2m+17) = Sfcprop(nb)%dt_cool(ix)!--- nsstm dt_cool
            sfc_var2(i,j,nvar_s2m+18) = Sfcprop(nb)%qrain(ix)  !--- nsstm qrain
          endif

  ! Noah MP
          if (Model%lsm == Model%lsm_noahmp) then
            sfc_var2(i,j,nvar_s2m+19) = Sfcprop(nb)%snowxy(ix)
            sfc_var2(i,j,nvar_s2m+20) = Sfcprop(nb)%tvxy(ix)
            sfc_var2(i,j,nvar_s2m+21) = Sfcprop(nb)%tgxy(ix)
            sfc_var2(i,j,nvar_s2m+22) = Sfcprop(nb)%canicexy(ix)
            sfc_var2(i,j,nvar_s2m+23) = Sfcprop(nb)%canliqxy(ix)
            sfc_var2(i,j,nvar_s2m+24) = Sfcprop(nb)%eahxy(ix)
            sfc_var2(i,j,nvar_s2m+25) = Sfcprop(nb)%tahxy(ix)
            sfc_var2(i,j,nvar_s2m+26) = Sfcprop(nb)%cmxy(ix)
            sfc_var2(i,j,nvar_s2m+27) = Sfcprop(nb)%chxy(ix)
            sfc_var2(i,j,nvar_s2m+28) = Sfcprop(nb)%fwetxy(ix)
            sfc_var2(i,j,nvar_s2m+29) = Sfcprop(nb)%sneqvoxy(ix)
            sfc_var2(i,j,nvar_s2m+30) = Sfcprop(nb)%alboldxy(ix)
            sfc_var2(i,j,nvar_s2m+31) = Sfcprop(nb)%qsnowxy(ix)
            sfc_var2(i,j,nvar_s2m+32) = Sfcprop(nb)%wslakexy(ix)
            sfc_var2(i,j,nvar_s2m+33) = Sfcprop(nb)%zwtxy(ix)
            sfc_var2(i,j,nvar_s2m+34) = Sfcprop(nb)%waxy(ix)
            sfc_var2(i,j,nvar_s2m+35) = Sfcprop(nb)%wtxy(ix)
            sfc_var2(i,j,nvar_s2m+36) = Sfcprop(nb)%lfmassxy(ix)
            sfc_var2(i,j,nvar_s2m+37) = Sfcprop(nb)%rtmassxy(ix)
            sfc_var2(i,j,nvar_s2m+38) = Sfcprop(nb)%stmassxy(ix)
            sfc_var2(i,j,nvar_s2m+39) = Sfcprop(nb)%woodxy(ix)
            sfc_var2(i,j,nvar_s2m+40) = Sfcprop(nb)%stblcpxy(ix)
            sfc_var2(i,j,nvar_s2m+41) = Sfcprop(nb)%fastcpxy(ix)
            sfc_var2(i,j,nvar_s2m+42) = Sfcprop(nb)%xsaixy(ix)
            sfc_var2(i,j,nvar_s2m+43) = Sfcprop(nb)%xlaixy(ix)
            sfc_var2(i,j,nvar_s2m+44) = Sfcprop(nb)%taussxy(ix)
            sfc_var2(i,j,nvar_s2m+45) = Sfcprop(nb)%smcwtdxy(ix)
            sfc_var2(i,j,nvar_s2m+46) = Sfcprop(nb)%deeprechxy(ix)
            sfc_var2(i,j,nvar_s2m+47) = Sfcprop(nb)%rechxy(ix)
            sfc_var2(i,j,nvar_s2m+48) = Sfcprop(nb)%drainncprv(ix)
            sfc_var2(i,j,nvar_s2m+49) = Sfcprop(nb)%draincprv(ix)
            sfc_var2(i,j,nvar_s2m+50) = Sfcprop(nb)%dsnowprv(ix)
            sfc_var2(i,j,nvar_s2m+51) = Sfcprop(nb)%dgraupelprv(ix)
            sfc_var2(i,j,nvar_s2m+52) = Sfcprop(nb)%diceprv(ix)
            sfc_var2(i,j,nvar_s2m+53) = Sfcprop(nb)%albdvis(ix)
            sfc_var2(i,j,nvar_s2m+54) = Sfcprop(nb)%albdnir(ix)
            sfc_var2(i,j,nvar_s2m+55) = Sfcprop(nb)%albivis(ix)
            sfc_var2(i,j,nvar_s2m+56) = Sfcprop(nb)%albinir(ix)
            sfc_var2(i,j,nvar_s2m+57) = Sfcprop(nb)%emiss(ix)
            sfc_var2(i,j,nvar_s2m+58) = Sfcprop(nb)%scolor(ix)
          endif

          !--- 3D variables
          do lsoil = 1,Model%lsoil
            sfc_var3(i,j,lsoil,1) = Sfcprop(nb)%stc(ix,lsoil) !--- stc
            sfc_var3(i,j,lsoil,2) = Sfcprop(nb)%smc(ix,lsoil) !--- smc
            sfc_var3(i,j,lsoil,3) = Sfcprop(nb)%slc(ix,lsoil) !--- slc
          enddo
  ! 5 Noah MP 3D
          if (Model%lsm == Model%lsm_noahmp) then

            do lsoil = -2,0
              sfc_var3sn(i,j,lsoil,4) = Sfcprop(nb)%snicexy(ix,lsoil)
              sfc_var3sn(i,j,lsoil,5) = Sfcprop(nb)%snliqxy(ix,lsoil)
              sfc_var3sn(i,j,lsoil,6) = Sfcprop(nb)%tsnoxy(ix,lsoil)
            enddo

            do lsoil = 1,Model%lsoil
              sfc_var3eq(i,j,lsoil,7)  = Sfcprop(nb)%smoiseq(ix,lsoil)
            enddo

            do lsoil = -2,4
              sfc_var3zn(i,j,lsoil,8)  = Sfcprop(nb)%zsnsoxy(ix,lsoil)
            enddo

          endif  ! Noah MP
        enddo
      enddo

      call write_restart(Sfc_restart)
      call close_file(Sfc_restart)
    endif
  end subroutine sfc_prop_restart_write

  subroutine sfc_prop_restart_write_coarse(Sfcprop, Atm_block, Model, coarse_domain, Grid, timestamp)
    type(GFS_sfcprop_type),      intent(in) :: Sfcprop(:)
    type(block_control_type),    intent(in) :: atm_block
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: coarse_domain
    type(GFS_grid_type),         intent(in) :: Grid(:)
    character(len=32), optional, intent(in) :: timestamp

    integer :: i, j, k, nb, ix, lsoil, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: nvar2m, nvar2o, nvar3, nvar_s2m

    integer :: is_coarse, ie_coarse, js_coarse, je_coarse, nx_coarse, ny_coarse
    character(len=64) :: fname
    character(len=32) :: fn_srf_coarse = 'sfc_data_coarse.nc'
    real(kind=kind_phys), allocatable, dimension(:,:) :: area, &
      dominant_sfc_type, dominant_vtype, dominant_stype, &
      tisfc_area_average, only_area_weighted_zorl, &
      only_area_weighted_canopy, coarsened_area_times_fice, &
      coarsened_area_times_sncovr, coarsened_area_times_vfrac
    logical, allocatable, dimension(:,:) :: sfc_type_mask, sfc_and_vtype_mask, sfc_and_stype_mask
    real(kind=kind_phys) :: FREEZING, VTYPE_LAND_ICE, STYPE_LAND_ICE, SHDMIN_CANOPY_THRESHOLD

    if (Model%cplflx .or. (Model%lsm .eq. Model%lsm_noahmp) .or. (Model%nstf_name(1) > 0)) then
      call mpp_error(FATAL, 'Coarse graining strategy not defined for land surface model configuration')
    endif

    nvar2m = 32
    nvar3 = 3

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)

    call mpp_get_compute_domain(coarse_domain, is_coarse, ie_coarse, js_coarse, je_coarse)
    nx_coarse = ie_coarse - is_coarse + 1
    ny_coarse = je_coarse - js_coarse + 1

    allocate(area(nx,ny))
    do nb = 1, Atm_block%nblks
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - isc + 1
        j = Atm_block%index(nb)%jj(ix) - jsc + 1
        area(i,j) = Grid(nb)%area(ix)
      enddo
    enddo

    allocate(dominant_sfc_type(nx,ny))
    allocate(dominant_vtype(nx,ny))
    allocate(dominant_stype(nx,ny))
    allocate(sfc_type_mask(nx,ny))
    allocate(sfc_and_vtype_mask(nx,ny))
    allocate(sfc_and_stype_mask(nx,ny))
    allocate(tisfc_area_average(nx_coarse,nx_coarse))
    allocate(only_area_weighted_zorl(nx_coarse,nx_coarse))
    allocate(only_area_weighted_canopy(nx_coarse,nx_coarse))
    allocate(coarsened_area_times_fice(nx_coarse,nx_coarse))
    allocate(coarsened_area_times_sncovr(nx_coarse,nx_coarse))
    allocate(coarsened_area_times_vfrac(nx_coarse,nx_coarse))

    if (.not. allocated(sfc_var2_coarse)) then
      allocate(sfc_var2_coarse(nx_coarse,ny_coarse,nvar2m))
      allocate(sfc_var3_coarse(nx_coarse,ny_coarse,Model%lsoil,nvar3))
    endif


    !--- Open the restart file and associate it with the Sfc_restart fileobject
    if (present(timestamp)) then
      fname='RESTART/'//trim(timestamp)//'.'//trim(fn_srf_coarse)
    else
      fname='RESTART/'//trim(fn_srf_coarse)
    endif

    if (open_file(Sfc_restart_coarse, fname, "overwrite", coarse_domain, is_restart=.true., dont_add_res_to_filename=.true.)) then

      !--- register surface properties for coarse restart
      call register_sfc_prop_restart_vars(Model, nx, ny, nvar_s2m, action="alloc")

      call register_coarse_sfc_prop_restart_fields(Model, sfc_var2_coarse, sfc_var3_coarse, nvar2m, nvar3)

      do nb = 1, Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          sfc_var2(i,j,1)  = Sfcprop(nb)%slmsk(ix)
          sfc_var2(i,j,2)  = Sfcprop(nb)%tsfc(ix)
          sfc_var2(i,j,3)  = Sfcprop(nb)%weasd(ix)
          sfc_var2(i,j,4)  = Sfcprop(nb)%tg3(ix)
          sfc_var2(i,j,5)  = Sfcprop(nb)%zorl(ix)
          sfc_var2(i,j,6)  = Sfcprop(nb)%alvsf(ix)
          sfc_var2(i,j,7)  = Sfcprop(nb)%alvwf(ix)
          sfc_var2(i,j,8)  = Sfcprop(nb)%alnsf(ix)
          sfc_var2(i,j,9)  = Sfcprop(nb)%alnwf(ix)
          sfc_var2(i,j,10) = Sfcprop(nb)%facsf(ix)
          sfc_var2(i,j,11) = Sfcprop(nb)%facwf(ix)
          sfc_var2(i,j,12) = Sfcprop(nb)%vfrac(ix)
          sfc_var2(i,j,13) = Sfcprop(nb)%canopy(ix)
          sfc_var2(i,j,14) = Sfcprop(nb)%f10m(ix)
          sfc_var2(i,j,15) = Sfcprop(nb)%t2m(ix)
          sfc_var2(i,j,16) = Sfcprop(nb)%q2m(ix)
          sfc_var2(i,j,17) = Sfcprop(nb)%vtype(ix)
          sfc_var2(i,j,18) = Sfcprop(nb)%stype(ix)
          sfc_var2(i,j,19) = Sfcprop(nb)%uustar(ix)
          sfc_var2(i,j,20) = Sfcprop(nb)%ffmm(ix)
          sfc_var2(i,j,21) = Sfcprop(nb)%ffhh(ix)
          sfc_var2(i,j,22) = Sfcprop(nb)%hice(ix)
          sfc_var2(i,j,23) = Sfcprop(nb)%fice(ix)
          sfc_var2(i,j,24) = Sfcprop(nb)%tisfc(ix)
          sfc_var2(i,j,25) = Sfcprop(nb)%tprcp(ix)
          sfc_var2(i,j,26) = Sfcprop(nb)%srflag(ix)
          sfc_var2(i,j,27) = Sfcprop(nb)%snowd(ix)
          sfc_var2(i,j,28) = Sfcprop(nb)%shdmin(ix)
          sfc_var2(i,j,29) = Sfcprop(nb)%shdmax(ix)
          sfc_var2(i,j,30) = Sfcprop(nb)%slope(ix)
          sfc_var2(i,j,31) = Sfcprop(nb)%snoalb(ix)
          sfc_var2(i,j,32) = Sfcprop(nb)%sncovr(ix)
          do lsoil = 1,Model%lsoil
            sfc_var3(i,j,lsoil,1) = Sfcprop(nb)%stc(ix,lsoil)
            sfc_var3(i,j,lsoil,2) = Sfcprop(nb)%smc(ix,lsoil)
            sfc_var3(i,j,lsoil,3) = Sfcprop(nb)%slc(ix,lsoil)
          enddo
        enddo
      enddo

      ! Coarse grain all the variables

      ! First coarse-grain the land surface type and upsample it back to the native resolution
      call block_mode(sfc_var2(:,:,1), sfc_var2_coarse(:,:,1))
      call block_upsample(sfc_var2_coarse(:,:,1), dominant_sfc_type)
      sfc_type_mask = (dominant_sfc_type .eq. sfc_var2(:,:,1))

      ! Then coarse-grain the vegetation and soil types and upsample them too
      call block_mode(sfc_var2(:,:,17), sfc_type_mask, sfc_var2_coarse(:,:,17))
      call block_upsample(sfc_var2_coarse(:,:,17), dominant_vtype)
      call block_mode(sfc_var2(isc:iec,jsc:jec,18), sfc_type_mask, sfc_var2_coarse(:,:,18))
      call block_upsample(sfc_var2_coarse(:,:,18), dominant_stype)

      sfc_and_vtype_mask = (sfc_type_mask .and. (dominant_vtype .eq. sfc_var2(:,:,17)))
      sfc_and_stype_mask = (sfc_type_mask .and. (dominant_stype .eq. sfc_var2(:,:,18)))

      ! Take the area weighted mean over full blocks for the surface temperature
      call weighted_block_average(area, sfc_var2(:,:,2), sfc_var2_coarse(:,:,2))

      ! Take the area weighted average over the dominant surface type for tg3
      call weighted_block_average(area, sfc_var2(:,:,4), sfc_type_mask, sfc_var2_coarse(:,:,4))

      ! Take the area weighted average over the dominant surface type for vfrac
      call weighted_block_average(area, sfc_var2(:,:,12), sfc_type_mask, sfc_var2_coarse(:,:,12))

      ! Take the area and vfrac weighted average over the dominant surface and vegetation type for zorl and canopy
      call weighted_block_average(area * sfc_var2(:,:,12), sfc_var2(:,:,5), sfc_and_vtype_mask, &
           sfc_var2_coarse(:,:,5))
      call weighted_block_average(area * sfc_var2(:,:,12), sfc_var2(:,:,13), sfc_and_vtype_mask, &
           sfc_var2_coarse(:,:,13))

      ! Also compute a simple area weighted average over the dominant surface and
      ! vegetation type for zorl and canopy; this will be used in the event that
      ! the sum of vfrac is equal to zero.
      call weighted_block_average(area, sfc_var2(:,:,5), sfc_and_vtype_mask, &
           only_area_weighted_zorl)
      call weighted_block_average(area, sfc_var2(:,:,13), sfc_and_vtype_mask, &
           only_area_weighted_canopy)

      call block_sum(area * sfc_var2(:,:,12), sfc_and_vtype_mask, coarsened_area_times_vfrac)

      ! If the dominant surface type is ocean or sea-ice then just use the
      ! area weighted average over the dominant surface and vegetation type for zorl or canopy.
      where (coarsened_area_times_vfrac .eq. 0.0)
         sfc_var2_coarse(:,:,5) = only_area_weighted_zorl
         sfc_var2_coarse(:,:,13) = only_area_weighted_canopy
      endwhere

      ! Take the area weighted average of the albedo variables
      call weighted_block_average(area, sfc_var2(:,:,6), sfc_var2_coarse(:,:,6))
      call weighted_block_average(area, sfc_var2(:,:,7), sfc_var2_coarse(:,:,7))
      call weighted_block_average(area, sfc_var2(:,:,8), sfc_var2_coarse(:,:,8))
      call weighted_block_average(area, sfc_var2(:,:,9), sfc_var2_coarse(:,:,9))
      call weighted_block_average(area, sfc_var2(:,:,10), sfc_var2_coarse(:,:,10))
      call weighted_block_average(area, sfc_var2(:,:,11), sfc_var2_coarse(:,:,11))

      ! Take the area weighted average of f10, t2m, q2m, uustar, ffmm, and ffhh
      call weighted_block_average(area, sfc_var2(:,:,14), sfc_var2_coarse(:,:,14))
      call weighted_block_average(area, sfc_var2(:,:,15), sfc_var2_coarse(:,:,15))
      call weighted_block_average(area, sfc_var2(:,:,16), sfc_var2_coarse(:,:,16))
      call weighted_block_average(area, sfc_var2(:,:,19), sfc_var2_coarse(:,:,19))
      call weighted_block_average(area, sfc_var2(:,:,20), sfc_var2_coarse(:,:,20))
      call weighted_block_average(area, sfc_var2(:,:,21), sfc_var2_coarse(:,:,21))

      ! Take the area weighted average over the dominant surface type for fice
      call weighted_block_average(area, sfc_var2(:,:,23), sfc_type_mask, sfc_var2_coarse(:,:,23))

      ! Compute the area weighted average of tpcrp
      call weighted_block_average(area, sfc_var2(:,:,25), sfc_var2_coarse(:,:,25))

      ! Take the mode for srflag
      call block_mode(sfc_var2(:,:,26), sfc_var2_coarse(:,:,26))

      ! Take the area weighted average for snow depth
      call weighted_block_average(area, sfc_var2(:,:,27), sfc_var2_coarse(:,:,27))

      ! Take the min and max over the dominant sfc type for shdmin and shdmax
      call block_min(sfc_var2(:,:,28), sfc_type_mask, sfc_var2_coarse(:,:,28))
      call block_max(sfc_var2(:,:,29), sfc_type_mask, sfc_var2_coarse(:,:,29))

      ! Take the masked block mode over the dominant surface type for slope
      call block_mode(sfc_var2(:,:,30), sfc_type_mask, sfc_var2_coarse(:,:,30))

      ! Take the block maximum for the snoalb
      call block_max(sfc_var2(:,:,31), sfc_type_mask, sfc_var2_coarse(:,:,31))

      ! Take the area weighted average over the dominant surface type for sncovr
      call weighted_block_average(area, sfc_var2(:,:,32), sfc_type_mask, sfc_var2_coarse(:,:,32))

      ! For sheleg take the area and sncovr weighted average; zero out any regions where the snow cover fraction is zero over the block.
      call weighted_block_average(area * sfc_var2(:,:,32), sfc_var2(:,:,3), sfc_var2_coarse(:,:,3))
      call block_sum(area * sfc_var2(:,:,32), coarsened_area_times_sncovr)
      where (coarsened_area_times_sncovr .eq. 0.0)
         sfc_var2_coarse(:,:,3) = 0.0
      endwhere

      ! Do something similar for hice
      call weighted_block_average(area * sfc_var2(:,:,23), sfc_var2(:,:,22), sfc_var2_coarse(:,:,22))
      call block_sum(area * sfc_var2(:,:,23), coarsened_area_times_fice)
      where (coarsened_area_times_fice .eq. 0.0)
         sfc_var2_coarse(:,:,22) = 0.0
      endwhere

      ! Over sea ice compute the area and ice fraction weighted average of tisfc; over all
      ! other surfaces use just the area weighted average of tisfc.
      call weighted_block_average(area * sfc_var2(:,:,23), sfc_var2(:,:,24), sfc_type_mask, sfc_var2_coarse(:,:,24))
      call weighted_block_average(area, sfc_var2(:,:,24), sfc_type_mask, tisfc_area_average)
      where (sfc_var2_coarse(:,:,1) .lt. 2.0)
         sfc_var2_coarse(:,:,24) = tisfc_area_average
      endwhere

      ! Apply corrections to 2D variables based on surface_chgres.F90
      FREEZING = 273.16
      VTYPE_LAND_ICE = 15.0
      STYPE_LAND_ICE = 16.0
      SHDMIN_CANOPY_THRESHOLD = 0.011

      ! Correction (1)
      ! Clip tsea and tg3 at 273.16 K if a cell contains land ice.
      where ((sfc_var2_coarse(:,:,2) .gt. FREEZING) .and. (sfc_var2_coarse(:,:,17) .eq. VTYPE_LAND_ICE))
         sfc_var2_coarse(:,:,2) = FREEZING
      endwhere
      where ((sfc_var2_coarse(:,:,4) .gt. FREEZING) .and. (sfc_var2_coarse(:,:,17) .eq. VTYPE_LAND_ICE))
         sfc_var2_coarse(:,:,4) = FREEZING
      endwhere

      ! Correction (2)
      ! If a cell contains land ice, make sure the soil type is ice.
      where (sfc_var2_coarse(:,:,17) .eq. VTYPE_LAND_ICE)
         sfc_var2_coarse(:,:,18) = STYPE_LAND_ICE
      endwhere

      ! Correction (3)
      ! If a cell does not contain vegetation, i.e. if shdmin < 0.011,
      ! then set the canopy moisture content to zero.
      where (sfc_var2_coarse(:,:,28) .lt. SHDMIN_CANOPY_THRESHOLD)
         sfc_var2_coarse(:,:,13) = 0.0
      endwhere

      ! Correction (4)
      ! If a cell contains land ice, then shdmin is set to zero.
      where (sfc_var2_coarse(:,:,17) .eq. VTYPE_LAND_ICE)
         sfc_var2_coarse(:,:,28) = 0.0
      endwhere

      ! For the 3D variables (all soil properties) take the area weighted average
      ! over the dominant surface and soil type.
      do num = 1,nvar3
        call weighted_block_average(area, sfc_var3(:,:,:,num), sfc_and_stype_mask, Model%lsoil, sfc_var3_coarse(:,:,:,num))
      enddo

      call write_restart(Sfc_restart_coarse)
      call close_file(Sfc_restart_coarse)
    endif


  end subroutine sfc_prop_restart_write_coarse


  subroutine register_coarse_sfc_prop_restart_fields(Model, var2, var3, nvar2, nvar3)
    type(IPD_control_type), intent(in) :: Model
    real(kind=kind_phys), target, intent(inout) :: var2(:,:,:)
    real(kind=kind_phys), target, intent(inout) :: var3(:,:,:,:)
    integer, intent(in) :: nvar2, nvar3

    integer :: is, ie, num, lsoil
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()
    integer, allocatable, dimension(:) :: buffer
    character(len=8) :: dim_names_2d(3), dim_names_3d(4)

    !--- register the axes for restarts
    call register_axis(Sfc_restart_coarse, 'xaxis_1', 'X')
    call register_field(Sfc_restart_coarse, 'xaxis_1', 'double', (/'xaxis_1'/))
    call register_variable_attribute(Sfc_restart_coarse, 'xaxis_1', 'cartesian_axis', 'X', str_len=1)
    call get_global_io_domain_indices(Sfc_restart_coarse, 'xaxis_1', is, ie, indices=buffer)
    call write_data(Sfc_restart_coarse, "xaxis_1", buffer)
    deallocate(buffer)

    call register_axis(Sfc_restart_coarse, 'yaxis_1', 'Y')
    call register_field(Sfc_restart_coarse, 'yaxis_1', 'double', (/'yaxis_1'/))
    call register_variable_attribute(Sfc_restart_coarse, 'yaxis_1', 'cartesian_axis', 'Y', str_len=1)
    call get_global_io_domain_indices(Sfc_restart_coarse, 'yaxis_1', is, ie, indices=buffer)
    call write_data(Sfc_restart_coarse, "yaxis_1", buffer)
    deallocate(buffer)

    call register_axis(Sfc_restart_coarse, 'zaxis_1', dimension_length=Model%lsoil)
    call register_field(Sfc_restart_coarse, 'zaxis_1', 'double', (/'zaxis_1'/))
    call register_variable_attribute(Sfc_restart_coarse, 'zaxis_1', 'cartesian_axis', 'Z', str_len=1)
    allocate( buffer(Model%lsoil) )
    do lsoil=1, Model%lsoil
       buffer(lsoil) = lsoil
    end do
    call write_data(Sfc_restart_coarse, 'zaxis_1', buffer)
    deallocate(buffer)

    call register_axis(Sfc_restart_coarse, 'Time', unlimited)
    call register_field(Sfc_restart_coarse, 'Time', 'double', (/'Time'/))
    call register_variable_attribute(Sfc_restart_coarse, 'Time', 'cartesian_axis', 'T', str_len=1)
    call write_data(Sfc_restart_coarse, 'Time', 1)

    !--- Assign dimensions to array for use in register_restart_field
    dim_names_2d(1) = "xaxis_1"
    dim_names_2d(2) = "yaxis_1"
    dim_names_2d(3) = "Time"
    dim_names_3d(1) = "xaxis_1"
    dim_names_3d(2) = "yaxis_1"
    dim_names_3d(3) = "zaxis_1"
    dim_names_3d(4) = "Time"

    do num = 1,nvar2
      var2_p => var2(:,:,num)
      call register_restart_fIeld(Sfc_restart_coarse, sfc_name2(num), var2_p, dim_names_2d)
      nullify(var2_p)
    enddo

    do num = 1,nvar3
      var3_p => var3(:,:,:,num)
      call register_restart_field(Sfc_restart_coarse, sfc_name3(num), var3_p, dim_names_3d)
      nullify(var3_p)
    enddo

 end subroutine register_coarse_sfc_prop_restart_fields

!----------------------------------------------------------------------
! phys_restart_read
!----------------------------------------------------------------------
!    creates and populates a data type which is then used to "register"
!    restart variables with the GFDL FMS restart subsystem.
!    calls a GFDL FMS routine to restore the data from a restart file.
!    calculates sncovr if it is not present in the restart file.
!
!    calls:  register_restart_field, restart_state, free_restart
!
!    opens:  phys_data.tile?.nc
!
!----------------------------------------------------------------------
  subroutine phys_restart_read (IPD_Restart, Atm_block, Model, fv_domain, enforce_rst_cksum)
    !--- interface variable definitions
    type(IPD_restart_type),      intent(in) :: IPD_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    logical,                     intent(in) :: enforce_rst_cksum
    !--- local variables
    integer :: i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: nvar2d, nvar3d
    character(len=64) :: fname
    character(len=8) :: dim_names_2d(3), dim_names_3d(4)
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()


    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)
    nvar2d = IPD_Restart%num2d
    nvar3d = IPD_Restart%num3d

    !--- Assign dimensions to array for use in register_restart_field
    dim_names_2d(1) = "xaxis_1"
    dim_names_2d(2) = "yaxis_1"
    dim_names_2d(3) = "Time"
    dim_names_3d(1) = "xaxis_1"
    dim_names_3d(2) = "yaxis_1"
    dim_names_3d(3) = "zaxis_1"
    dim_names_3d(4) = "Time"

    !--- Open the restart file and associate it with the Phy_restart fileobject
    fname='INPUT/'//trim(fn_phy)
    if (open_file(Phy_restart, fname, "read", fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)) then

      !--- register the axes for restarts
      call register_axis(Phy_restart, "xaxis_1", "X")
      call register_axis(Phy_restart, "yaxis_1", "Y")
      call register_axis(Phy_restart, "zaxis_1", npz)
      call register_axis(Phy_restart, "Time", unlimited)

      !--- register the restart fields
      if (.not. allocated(phy_var2)) then
        allocate (phy_var2(nx,ny,nvar2d))
        allocate (phy_var3(nx,ny,npz,nvar3d))
        phy_var2 = 0.0_kind_phys
        phy_var3 = 0.0_kind_phys

        do num = 1,nvar2d
          var2_p => phy_var2(:,:,num)
          call register_restart_field (Phy_restart, trim(IPD_Restart%name2d(num)), &
                                       var2_p, dim_names_2d, is_optional=.true.)
          nullify(var2_p)
        enddo
        do num = 1,nvar3d
          var3_p => phy_var3(:,:,:,num)
          call register_restart_field (Phy_restart, trim(IPD_restart%name3d(num)), &
                                       var3_p, dim_names_3d, is_optional=.true.)
          nullify(var3_p)
        enddo
      endif

      !--- read the surface restart/data
      call mpp_error(NOTE,'reading physics restart data from INPUT/phy_data.tile*.nc')
      call read_restart(Phy_restart, ignore_checksum=enforce_rst_cksum)
      call close_file(Phy_restart)
    else
      call mpp_error(NOTE,'No physics restarts - cold starting physical parameterizations')
      return
    endif

    !--- place the data into the block GFS containers
    !--- phy_var* variables
    do num = 1,nvar2d
      do nb = 1,Atm_block%nblks
        do ix = 1, Atm_block%blksz(nb)
          i = Atm_block%index(nb)%ii(ix) - isc + 1
          j = Atm_block%index(nb)%jj(ix) - jsc + 1
          IPD_Restart%data(nb,num)%var2p(ix) = phy_var2(i,j,num)
        enddo
      enddo
    enddo
    do num = 1,nvar3d
      do nb = 1,Atm_block%nblks
        do k=1,npz
          do ix = 1, Atm_block%blksz(nb)
            i = Atm_block%index(nb)%ii(ix) - isc + 1
            j = Atm_block%index(nb)%jj(ix) - jsc + 1
            IPD_Restart%data(nb,num)%var3p(ix,k) = phy_var3(i,j,k,num)
          enddo
        enddo
      enddo
    enddo

    if (allocated(phy_var2)) deallocate(phy_var2)
    if (allocated(phy_var3)) deallocate(phy_var3)

  end subroutine phys_restart_read


!----------------------------------------------------------------------
! phys_restart_write
!----------------------------------------------------------------------
!    routine to write out GFS surface restarts via the GFDL FMS restart
!    subsystem.
!    takes an optional argument to append timestamps for intermediate
!    restarts.
!
!    calls:  register_restart_field, save_restart
!----------------------------------------------------------------------
  subroutine phys_restart_write (IPD_Restart, Atm_block, Model, fv_domain, timestamp)
    !--- interface variable definitions
    type(IPD_restart_type),      intent(in) :: IPD_Restart
    type(block_control_type),    intent(in) :: Atm_block
    type(IPD_control_type),      intent(in) :: Model
    type(domain2d),              intent(in) :: fv_domain
    character(len=32), optional, intent(in) :: timestamp
    !--- local variables
    integer :: is, ie, i, j, k, nb, ix, num
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: nvar2d, nvar3d
    integer, allocatable, dimension(:) :: buffer
    character(len=64) :: fname
    character(len=8) :: dim_names_2d(3), dim_names_3d(4)
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()

    isc = Atm_block%isc
    iec = Atm_block%iec
    jsc = Atm_block%jsc
    jec = Atm_block%jec
    npz = Atm_block%npz
    nx = (iec - isc + 1)
    ny = (jec - jsc + 1)
    nvar2d = IPD_Restart%num2d
    nvar3d = IPD_Restart%num3d

    !--- Assign dimensions to array for use in register_restart_field
    dim_names_2d(1) = "xaxis_1"
    dim_names_2d(2) = "yaxis_1"
    dim_names_2d(3) = "Time"
    dim_names_3d(1) = "xaxis_1"
    dim_names_3d(2) = "yaxis_1"
    dim_names_3d(3) = "zaxis_1"
    dim_names_3d(4) = "Time"

    !--- Open the restart file and associate it with the Phy_restart fileobject
    if (present(timestamp)) then
      fname='RESTART/'//trim(timestamp)//'.'//trim(fn_phy)
    else
      fname='RESTART/'//trim(fn_phy)
    endif

    if (open_file(Phy_restart, fname, "overwrite", fv_domain, is_restart=.true., dont_add_res_to_filename=.true.)) then

      !--- register the axes for restarts
      call register_axis(Phy_restart, "xaxis_1", "X")
      call register_field(Phy_restart, 'xaxis_1', 'double', (/'xaxis_1'/))
      call register_variable_attribute(Phy_restart, 'xaxis_1', 'cartesian_axis', 'X', str_len=1)
      call get_global_io_domain_indices(Phy_restart, 'xaxis_1', is, ie, indices=buffer)
      call write_data(Phy_restart, "xaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Phy_restart, "yaxis_1", "Y")
      call register_field(Phy_restart, 'yaxis_1', 'double', (/'yaxis_1'/))
      call register_variable_attribute(Phy_restart, 'yaxis_1', 'cartesian_axis', 'Y', str_len=1)
      call get_global_io_domain_indices(Phy_restart, 'yaxis_1', is, ie, indices=buffer)
      call write_data(Phy_restart, "yaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Phy_restart, "zaxis_1", npz)
      call register_field(Phy_restart, 'zaxis_1', 'double', (/'zaxis_1'/))
      call register_variable_attribute(Phy_restart, 'zaxis_1', 'cartesian_axis', 'Z', str_len=1)
      allocate( buffer(npz) )
      do i=1, npz
         buffer(i)=i
      end do
      call write_data(Phy_restart, "zaxis_1", buffer)
      deallocate(buffer)

      call register_axis(Phy_restart, "Time", unlimited)
      call register_field(Phy_restart, 'Time', 'double', (/'Time'/))
      call register_variable_attribute(Phy_restart, 'Time', 'cartesian_axis', 'T', str_len=1)
      call write_data(Phy_restart, "Time", 1)

      !--- register the restart fields
      if (.not. allocated(phy_var2)) then
        allocate (phy_var2(nx,ny,nvar2d))
        allocate (phy_var3(nx,ny,npz,nvar3d))
        phy_var2 = 0.0_kind_phys
        phy_var3 = 0.0_kind_phys

        do num = 1,nvar2d
          var2_p => phy_var2(:,:,num)
          call register_restart_field (Phy_restart, trim(IPD_Restart%name2d(num)), &
                                       var2_p, dim_names_2d)
          nullify(var2_p)
        enddo
        do num = 1,nvar3d
          var3_p => phy_var3(:,:,:,num)
          call register_restart_field (Phy_restart, trim(IPD_restart%name3d(num)), &
                                       var3_p, dim_names_3d)
          nullify(var3_p)
        enddo
      endif

      !--- 2D variables
      do num = 1,nvar2d
        do nb = 1,Atm_block%nblks
          do ix = 1, Atm_block%blksz(nb)
            i = Atm_block%index(nb)%ii(ix) - isc + 1
            j = Atm_block%index(nb)%jj(ix) - jsc + 1
            phy_var2(i,j,num) = IPD_Restart%data(nb,num)%var2p(ix)
          enddo
        enddo
      enddo
      !--- 3D variables
      do num = 1,nvar3d
        do nb = 1,Atm_block%nblks
          do k=1,npz
            do ix = 1, Atm_block%blksz(nb)
              i = Atm_block%index(nb)%ii(ix) - isc + 1
              j = Atm_block%index(nb)%jj(ix) - jsc + 1
              phy_var3(i,j,k,num) = IPD_Restart%data(nb,num)%var3p(ix,k)
            enddo
          enddo
        enddo
      enddo

      call write_restart(Phy_restart)
      call close_file(Phy_restart)

      if (allocated(phy_var2)) deallocate (phy_var2)
      if (allocated(phy_var3)) deallocate (phy_var3)
    endif

  end subroutine phys_restart_write

  subroutine register_diag_manager_controlled_diagnostics(Time, Sfcprop, IntDiag, Model, nblks, axes)
    type(time_type), intent(in) :: Time
    type(Gfs_sfcprop_type), intent(in) :: Sfcprop(:)
    type(GFS_diag_type), intent(in) :: IntDiag(:)
    type(IPD_control_type), intent(in) :: Model
    integer, intent(in) :: nblks
    integer, intent(in) :: axes(4)

    integer :: nb
    integer :: index = 1

    if (Model%ldiag3d) then
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_air_temperature_due_to_longwave_heating'
      Diag_diag_manager_controlled(index)%desc = 'temperature tendency due to longwave radiation'
      Diag_diag_manager_controlled(index)%unit = 'K/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%t_dt(:,:,1)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_air_temperature_due_to_shortwave_heating'
      Diag_diag_manager_controlled(index)%desc = 'temperature tendency due to shortwave radiation'
      Diag_diag_manager_controlled(index)%unit = 'K/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%t_dt(:,:,2)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_air_temperature_due_to_turbulence'
      Diag_diag_manager_controlled(index)%desc = 'temperature tendency due to turbulence scheme'
      Diag_diag_manager_controlled(index)%unit = 'K/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%t_dt(:,:,3)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_air_temperature_due_to_deep_convection'
      Diag_diag_manager_controlled(index)%desc = 'temperature tendency due to deep convection'
      Diag_diag_manager_controlled(index)%unit = 'K/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%t_dt(:,:,4)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_air_temperature_due_to_shallow_convection'
      Diag_diag_manager_controlled(index)%desc = 'temperature tendency due to shallow convection'
      Diag_diag_manager_controlled(index)%unit = 'K/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%t_dt(:,:,5)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_air_temperature_due_to_microphysics'
      Diag_diag_manager_controlled(index)%desc = 'temperature tendency due to micro-physics'
      Diag_diag_manager_controlled(index)%unit = 'K/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%t_dt(:,:,6)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_air_temperature_due_to_dissipation_of_gravity_waves'
      Diag_diag_manager_controlled(index)%desc = 'temperature tendency due to gravity wave drag'
      Diag_diag_manager_controlled(index)%unit = 'K/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%t_dt(:,:,7)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky'
      Diag_diag_manager_controlled(index)%desc = 'temperature tendency due to clear sky longwave radiation'
      Diag_diag_manager_controlled(index)%unit = 'K/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%t_dt(:,:,8)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky'
      Diag_diag_manager_controlled(index)%desc = 'temperature tendency due to clear sky shortwave radiation'
      Diag_diag_manager_controlled(index)%unit = 'K/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%t_dt(:,:,9)
      enddo

   ! Vertically integrated instantaneous temperature tendency diagnostics
      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_air_temperature_due_to_longwave_heating'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated temperature tendency due to longwave radiation'
      Diag_diag_manager_controlled(index)%unit = 'W/m**2'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%t_dt_int(:,1)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_air_temperature_due_to_shortwave_heating'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated temperature tendency due to shortwave radiation'
      Diag_diag_manager_controlled(index)%unit = 'W/m**2'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%t_dt_int(:,2)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_air_temperature_due_to_turbulence'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated temperature tendency due to turbulence scheme'
      Diag_diag_manager_controlled(index)%unit = 'W/m**2'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%t_dt_int(:,3)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_air_temperature_due_to_deep_convection'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated temperature tendency due to deep convection'
      Diag_diag_manager_controlled(index)%unit = 'W/m**2'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%t_dt_int(:,4)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_air_temperature_due_to_shallow_convection'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated temperature tendency due to shallow convection'
      Diag_diag_manager_controlled(index)%unit = 'W/m**2'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%t_dt_int(:,5)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_air_temperature_due_to_microphysics'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated temperature tendency due to micro-physics'
      Diag_diag_manager_controlled(index)%unit = 'W/m**2'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%t_dt_int(:,6)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_air_temperature_due_to_dissipation_of_gravity_waves'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated temperature tendency due to gravity wave drag'
      Diag_diag_manager_controlled(index)%unit = 'W/m**2'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%t_dt_int(:,7)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated temperature tendency due to clear sky longwave radiation'
      Diag_diag_manager_controlled(index)%unit = 'W/m**2'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%t_dt_int(:,8)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated temperature tendency due to clear sky shortwave radiation'
      Diag_diag_manager_controlled(index)%unit = 'W/m**2'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%t_dt_int(:,9)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_specific_humidity_due_to_turbulence'
      Diag_diag_manager_controlled(index)%desc = 'water vapor tendency due to turbulence scheme'
      Diag_diag_manager_controlled(index)%unit = 'kg/kg/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%q_dt(:,:,1)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_specific_humidity_due_to_deep_convection'
      Diag_diag_manager_controlled(index)%desc = 'water vapor tendency due to deep convection'
      Diag_diag_manager_controlled(index)%unit = 'kg/kg/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%q_dt(:,:,2)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_specific_humidity_due_to_shallow_convection'
      Diag_diag_manager_controlled(index)%desc = 'water vapor tendency due to shallow convection'
      Diag_diag_manager_controlled(index)%unit = 'kg/kg/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%q_dt(:,:,3)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_specific_humidity_due_to_microphysics'
      Diag_diag_manager_controlled(index)%desc = 'water vapor tendency due to microphysics'
      Diag_diag_manager_controlled(index)%unit = 'kg/kg/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%q_dt(:,:,4)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'tendency_of_specific_humidity_due_to_change_in_atmosphere_mass'
      Diag_diag_manager_controlled(index)%desc = 'residual water vapor tendency'
      Diag_diag_manager_controlled(index)%unit = 'kg/kg/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'mass_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%q_dt(:,:,5)
      enddo

      ! Vertically integrated instantaneous specific humidity tendency diagnostics
      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_specific_humidity_due_to_turbulence'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated water vapor tendency due to turbulence scheme'
      Diag_diag_manager_controlled(index)%unit = 'kg/m**2/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%q_dt_int(:,1)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_specific_humidity_due_to_deep_convection'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated water vapor tendency due to deep convection'
      Diag_diag_manager_controlled(index)%unit = 'kg/m**2/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%q_dt_int(:,2)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_specific_humidity_due_to_shallow_convection'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated water vapor tendency due to shallow convection'
      Diag_diag_manager_controlled(index)%unit = 'kg/m**2/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%q_dt_int(:,3)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_specific_humidity_due_to_microphysics'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated water vapor tendency due to microphysics'
      Diag_diag_manager_controlled(index)%unit = 'kg/m**2/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%q_dt_int(:,4)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 2
      Diag_diag_manager_controlled(index)%name = 'vertically_integrated_tendency_of_specific_humidity_due_to_change_in_atmosphere_mass'
      Diag_diag_manager_controlled(index)%desc = 'vertically integrated residual water vapor tendency'
      Diag_diag_manager_controlled(index)%unit = 'kg/m**2/s'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%q_dt_int(:,5)
      enddo

      index = index + 1
      Diag_diag_manager_controlled(index)%axes = 3
      Diag_diag_manager_controlled(index)%name = 'co2'
      Diag_diag_manager_controlled(index)%desc = 'carbon dioxide concentration'
      Diag_diag_manager_controlled(index)%unit = 'volume mixing ratio'
      Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
      Diag_diag_manager_controlled(index)%coarse_graining_method = AREA_WEIGHTED
      allocate (Diag_diag_manager_controlled(index)%data(nblks))
      do nb = 1,nblks
         Diag_diag_manager_controlled(index)%data(nb)%var3 => IntDiag(nb)%co2(:,:)
      enddo

   endif

   index = index + 1
   Diag_diag_manager_controlled(index)%axes = 0
   Diag_diag_manager_controlled(index)%name = 'global_mean_co2'
   Diag_diag_manager_controlled(index)%desc = 'global mean carbon dioxide concentration'
   Diag_diag_manager_controlled(index)%unit = 'volume mixing ratio'
   Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
   Diag_diag_manager_controlled(index)%coarse_graining_method = AREA_WEIGHTED
   allocate (Diag_diag_manager_controlled(index)%data(nblks))
   do nb = 1,nblks
       Diag_diag_manager_controlled(index)%data(nb)%var2 => IntDiag(nb)%column_moles_co2_per_square_meter
       Diag_diag_manager_controlled(index)%data(nb)%var21 => IntDiag(nb)%column_moles_dry_air_per_square_meter
   enddo

   index = index + 1
   Diag_diag_manager_controlled(index)%axes = 2
   Diag_diag_manager_controlled(index)%name = 'ocean_fraction'
   Diag_diag_manager_controlled(index)%desc = 'fraction of grid cell classified as ocean type'
   Diag_diag_manager_controlled(index)%unit = 'fraction'
   Diag_diag_manager_controlled(index)%mod_name = 'gfs_sfc'
   Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
   allocate (Diag_diag_manager_controlled(index)%data(nblks))
   do nb = 1,nblks
     Diag_diag_manager_controlled(index)%data(nb)%var2 => Sfcprop(nb)%slmsk(:)
   enddo

   index = index + 1
   Diag_diag_manager_controlled(index)%axes = 2
   Diag_diag_manager_controlled(index)%name = 'land_fraction'
   Diag_diag_manager_controlled(index)%desc = 'fraction of grid cell classified as land type'
   Diag_diag_manager_controlled(index)%unit = 'fraction'
   Diag_diag_manager_controlled(index)%mod_name = 'gfs_sfc'
   Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
   allocate (Diag_diag_manager_controlled(index)%data(nblks))
   do nb = 1,nblks
     Diag_diag_manager_controlled(index)%data(nb)%var2 => Sfcprop(nb)%slmsk(:)
   enddo

   index = index + 1
   Diag_diag_manager_controlled(index)%axes = 2
   Diag_diag_manager_controlled(index)%name = 'sea_ice_fraction'
   Diag_diag_manager_controlled(index)%desc = 'fraction of grid cell classified as sea ice type'
   Diag_diag_manager_controlled(index)%unit = 'fraction'
   Diag_diag_manager_controlled(index)%mod_name = 'gfs_sfc'
   Diag_diag_manager_controlled(index)%coarse_graining_method = 'area_weighted'
   allocate (Diag_diag_manager_controlled(index)%data(nblks))
   do nb = 1,nblks
     Diag_diag_manager_controlled(index)%data(nb)%var2 => Sfcprop(nb)%slmsk(:)
   enddo

   index = index + 1
   Diag_diag_manager_controlled(index)%axes = 2
   Diag_diag_manager_controlled(index)%name = 'mixed_layer_depth'
   Diag_diag_manager_controlled(index)%desc = 'ocean mixed layer depth'
   Diag_diag_manager_controlled(index)%unit = 'm'
   Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
   Diag_diag_manager_controlled(index)%coarse_graining_method = AREA_WEIGHTED
   allocate (Diag_diag_manager_controlled(index)%data(nblks))
   do nb = 1,nblks
     Diag_diag_manager_controlled(index)%data(nb)%var2 => Sfcprop(nb)%mld(:)
   enddo

   index = index + 1
   Diag_diag_manager_controlled(index)%axes = 2
   Diag_diag_manager_controlled(index)%name = 'prescribed_mixed_layer_depth'
   Diag_diag_manager_controlled(index)%desc = 'prescribed ocean mixed layer depth'
   Diag_diag_manager_controlled(index)%unit = 'm'
   Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
   Diag_diag_manager_controlled(index)%coarse_graining_method = AREA_WEIGHTED
   allocate (Diag_diag_manager_controlled(index)%data(nblks))
   do nb = 1,nblks
     Diag_diag_manager_controlled(index)%data(nb)%var2 => Sfcprop(nb)%mldclim(:)
   enddo

   index = index + 1
   Diag_diag_manager_controlled(index)%axes = 2
   Diag_diag_manager_controlled(index)%name = 'prescribed_qflux'
   Diag_diag_manager_controlled(index)%desc = 'prescribed ocean Q-flux'
   Diag_diag_manager_controlled(index)%unit = 'W/m**2'
   Diag_diag_manager_controlled(index)%mod_name = 'gfs_phys'
   Diag_diag_manager_controlled(index)%coarse_graining_method = AREA_WEIGHTED
   allocate (Diag_diag_manager_controlled(index)%data(nblks))
   do nb = 1,nblks
     Diag_diag_manager_controlled(index)%data(nb)%var2 => Sfcprop(nb)%qfluxadj(:)
   enddo

   do index = 1, DIAG_SIZE
      if (trim(Diag_diag_manager_controlled(index)%name) .eq. '') exit  ! No need to populate non-existent diagnostics
        if (Diag_diag_manager_controlled(index)%axes .gt. 0) then
          Diag_diag_manager_controlled(index)%id = register_diag_field(trim(Diag_diag_manager_controlled(index)%mod_name), &
              & trim(Diag_diag_manager_controlled(index)%name),  &
              & axes(1:Diag_diag_manager_controlled(index)%axes), Time, trim(Diag_diag_manager_controlled(index)%desc), &
              & trim(Diag_diag_manager_controlled(index)%unit), missing_value=real(missing_value))
        else
          ! Scalar diagnostics are registered without any axes, so must be handled differently.
          Diag_diag_manager_controlled(index)%id = register_diag_field(trim(Diag_diag_manager_controlled(index)%mod_name), &
              & trim(Diag_diag_manager_controlled(index)%name), Time, trim(Diag_diag_manager_controlled(index)%desc), &
              & trim(Diag_diag_manager_controlled(index)%unit), missing_value=real(missing_value))
        endif
   enddo
  end subroutine register_diag_manager_controlled_diagnostics

!-------------------------------------------------------------------------
!--- gfdl_diag_register ---
!-------------------------------------------------------------------------
!    creates and populates a data type which is then used to "register"
!    GFS physics diagnostic variables with the GFDL FMS diagnostic manager.
!    includes short & long names, units, conversion factors, etc.
!    there is no copying of data, but instead a clever use of pointers.
!    calls a GFDL FMS routine to register diagnositcs and compare against
!    the diag_table to determine what variables are to be output.
!
!    calls:  register_diag_field
!-------------------------------------------------------------------------
!    Current sizes
!    13+NFXR - radiation
!    76+pl_coeff - physics
!-------------------------------------------------------------------------
  subroutine gfdl_diag_register(Time, Sfcprop, Gfs_diag, Model, Cldprop, Atm_block, axes)
    use physcons,  only: con_g
!--- subroutine interface variable definitions
    type(time_type),           intent(in) :: Time
    type(Gfs_sfcprop_type),    intent(in) :: Sfcprop(:)
    type(GFS_diag_type),       intent(in) :: Gfs_diag(:)
    type(IPD_control_type),    intent(in) :: Model
    type(GFS_cldprop_type),    intent(in) :: Cldprop(:)
    type (block_control_type), intent(in) :: Atm_block
    integer, dimension(4),     intent(in) :: axes
!--- local variables
    integer :: idx, num, nb, nblks, nx, ny, k
    integer, allocatable :: blksz(:)
    character(len=2) :: xtra
    real(kind=kind_phys), parameter :: cn_one = 1._kind_phys
    real(kind=kind_phys), parameter :: cn_100 = 100._kind_phys
    real(kind=kind_phys), parameter :: cn_th  = 1000._kind_phys
    real(kind=kind_phys), parameter :: cn_hr  = 3600._kind_phys

    nblks = Atm_block%nblks
    allocate (blksz(nblks))
    blksz(:) = Atm_block%blksz(:)

    isco   = Atm_block%isc
    ieco   = Atm_block%iec
    jsco   = Atm_block%jsc
    jeco   = Atm_block%jec
    levo   = Model%levs

    Diag(:)%id = -99
    Diag(:)%axes = -99
    Diag(:)%cnvfac = cn_one
    Diag(:)%time_avg = .FALSE.
    Diag(:)%time_avg_kind = ''
    Diag(:)%mask = ''
    Diag(:)%coarse_graining_method = 'area_weighted'
    Diag(:)%intpl_method = 'nearest_stod'

    idx = 0

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ALBDOsfc'
    Diag(idx)%desc = 'surface albedo'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%mask = 'positive_flux'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2  => Gfs_diag(nb)%fluxr(:,3)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,4)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'USWRFsfc'
    Diag(idx)%desc = 'Interval-averaged unadjusted upward shortwave flux at the surface'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,3)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DSWRFsfc'
    Diag(idx)%desc = 'Interval-averaged unadjusted downward shortwave flux at the surface'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,4)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DLWRFsfc'
    Diag(idx)%desc = 'Interval-averaged unadjusted downward longwave flux at the surface'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_lw'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,19)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ULWRFsfc'
    Diag(idx)%desc = 'Interval-averaged unadjusted upward longwave flux at the surface'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_lw'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,20)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'duvb_ave'
    Diag(idx)%desc = 'UV-B Downward Solar Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,21)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cduvb_ave'
    Diag(idx)%desc = 'Clear sky UV-B Downward Solar Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,22)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'vbdsf_ave'
    Diag(idx)%desc = 'Visible Beam Downward Solar Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,24)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'vddsf_ave'
    Diag(idx)%desc = 'Visible Diffuse Downward Solar Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,25)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'nbdsf_ave'
    Diag(idx)%desc = 'Near IR Beam Downward Solar Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,26)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'nddsf_ave'
    Diag(idx)%desc = 'Near IR Diffuse Downward Solar Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,27)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'csulf_avetoa'
    Diag(idx)%desc = 'Clear Sky Upward Long Wave Flux at toa'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_lw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,28)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'csusf_avetoa'
    Diag(idx)%desc = 'Clear Sky Upward Short Wave Flux at toa'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,29)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'csdlf_ave'
    Diag(idx)%desc = 'Clear Sky Downward Long Wave Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_lw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,30)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'csusf_ave'
    Diag(idx)%desc = 'Clear Sky Upward Short Wave Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,31)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'csdsf_ave'
    Diag(idx)%desc = 'Clear Sky Downward Short Wave Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,32)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'csulf_ave'
    Diag(idx)%desc = 'Clear Sky Upward Long Wave Flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_lw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,33)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DSWRFtoa'
    Diag(idx)%desc = 'top of atmos downward shortwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,23)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'USWRFtoa'
    Diag(idx)%desc = 'top of atmos upward shortwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_sw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,2)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ULWRFtoa'
    Diag(idx)%desc = 'top of atmos upward longwave flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_lw'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,1)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDCclm'
    Diag(idx)%desc = 'atmos column total cloud cover [%]'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,17)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDCbndcl'
    Diag(idx)%desc = 'boundary layer cloud layer total cloud cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,18)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDCcnvcl'
    Diag(idx)%desc = 'convective cloud layer total cloud cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Cldprop(nb)%cv(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'PREScnvclt'
    Diag(idx)%desc = 'pressure at convective cloud top level'
    Diag(idx)%unit = 'pa'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%mask = 'cldmask'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Cldprop(nb)%cvt(:)
      Diag(idx)%data(nb)%var21 => Cldprop(nb)%cv(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'PREScnvclb'
    Diag(idx)%desc = 'pressure at convective cloud bottom level'
    Diag(idx)%unit = 'pa'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%mask = 'cldmask'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Cldprop(nb)%cvb(:)
      Diag(idx)%data(nb)%var21 => Cldprop(nb)%cv(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDChcl'
    Diag(idx)%desc = 'high cloud level total cloud cover [%]'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'PRES_avehct'
    Diag(idx)%desc = 'pressure high cloud top level'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    Diag(idx)%mask = "cldmask_ratio"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,8)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'PRES_avehcb'
    Diag(idx)%desc = 'pressure high cloud bottom level'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    Diag(idx)%mask = "cldmask_ratio"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,11)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TEMP_avehct'
    Diag(idx)%desc = 'temperature high cloud top level'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    Diag(idx)%mask = "cldmask_ratio"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,14)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDCmcl'
    Diag(idx)%desc = 'mid cloud level total cloud cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,6)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'PRES_avemct'
    Diag(idx)%desc = 'pressure middle cloud top level'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    Diag(idx)%mask = "cldmask_ratio"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,9)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,6)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'PRES_avemcb'
    Diag(idx)%desc = 'pressure middle cloud bottom level'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    Diag(idx)%mask = "cldmask_ratio"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,12)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,6)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TEMP_avemct'
    Diag(idx)%desc = 'temperature middle cloud top level'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    Diag(idx)%mask = "cldmask_ratio"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,15)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,6)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TCDClcl'
    Diag(idx)%desc = 'low cloud level total cloud cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_100
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,7)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'PRES_avelct'
    Diag(idx)%desc = 'pressure low cloud top level'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    Diag(idx)%mask = "cldmask_ratio"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,10)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,7)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'PRES_avelcb'
    Diag(idx)%desc = 'pressure low cloud bottom level'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    Diag(idx)%mask = "cldmask_ratio"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2  => Gfs_diag(nb)%fluxr(:,13)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,7)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'TEMP_avelct'
    Diag(idx)%desc = 'temperature low cloud top level'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'rad_swlw_min'
    Diag(idx)%mask = "cldmask_ratio"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2  => Gfs_diag(nb)%fluxr(:,16)
      Diag(idx)%data(nb)%var21 => Gfs_diag(nb)%fluxr(:,7)
    enddo

!--- accumulated diagnostics ---
    do num = 1,Model%nfxr
      write (xtra,'(I2.2)') num
      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'fluxr_'//trim(xtra)
      Diag(idx)%desc = 'fluxr diagnostic '//trim(xtra)//' - GFS radiation'
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%fluxr(:,num)
      enddo
    enddo

!--- averaged diagnostics ---
    do num = 1,Model%nkld
      write (xtra,'(I2.2)') num
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'cloud_'//trim(xtra)
      Diag(idx)%desc = 'cloud diagnostic '//trim(xtra)//' - GFS radiation'
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cloud(:,:,num)
      enddo
    enddo

!--- the next two appear to be appear to be coupling fields in gloopr
!--- each has four elements
!rab    do num = 1,4
!rab      write (xtra,'(I1)') num
!rab      idx = idx + 1
!rab      Diag(idx)%axes = 2
!rab      Diag(idx)%name = 'dswcmp_'//trim(xtra)
!rab      Diag(idx)%desc = 'dswcmp dagnostic '//trim(xtra)//' - GFS radiation'
!rab      Diag(idx)%unit = 'XXX'
!rab      Diag(idx)%mod_name = 'gfs_phys'
!rab      do nb = 1,nblks
!rab        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dswcmp(:,num)
!rab      enddo
!rab    enddo
!rab
!rab    do num = 1,4
!rab      write (xtra,'(I1)') num
!rab      idx = idx + 1
!rab      Diag(idx)%axes = 2
!rab      Diag(idx)%name = 'uswcmp_'//trim(xtra)
!rab      Diag(idx)%desc = 'uswcmp dagnostic '//trim(xtra)//' - GFS radiation'
!rab      Diag(idx)%unit = 'XXX'
!rab      Diag(idx)%mod_name = 'gfs_phys'
!rab      do nb = 1,nblks
!rab        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%uswcmp(:,num)
!rab      enddo
!rab    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_upfxc'
    Diag(idx)%desc = 'total sky upward sw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topfsw(:)%upfxc
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_dnfxc'
    Diag(idx)%desc = 'total sky downward sw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topfsw(:)%dnfxc
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sw_upfx0'
    Diag(idx)%desc = 'clear sky upward sw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topfsw(:)%upfx0
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'lw_upfxc'
    Diag(idx)%desc = 'total sky upward lw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topflw(:)%upfxc
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'lw_upfx0'
    Diag(idx)%desc = 'clear sky upward lw flux at toa - GFS radiation'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%topflw(:)%upfx0
    enddo

!--- physics accumulated diagnostics ---
    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ssrun_acc'
    Diag(idx)%desc = 'surface storm water runoff - GFS lsm'
    Diag(idx)%unit = 'kg/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%srunoff(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'evbs_ave'
    Diag(idx)%desc = 'Direct Evaporation from Bare Soil - GFS lsm'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%evbsa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'evcw_ave'
    Diag(idx)%desc = 'Canopy water evaporation - GFS lsm'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%evcwa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snohf_ave'
    Diag(idx)%desc = 'Snow Phase Change Heat Flux - GFS lsm'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%snohfa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'trans_ave'
    Diag(idx)%desc = 'transpiration - GFS lsm'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%transa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sbsno_ave'
    Diag(idx)%desc = 'Sublimation (evaporation from snow) - GFS lsm'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%sbsnoa(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snowc_ave'
    Diag(idx)%desc = 'snow cover - GFS lsm'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%cnvfac = cn_100
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%snowca(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'soilm'
    Diag(idx)%desc = 'total column soil moisture content [kg/m**2]'
    Diag(idx)%unit = 'kg/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%mask = "land_only"
    Diag(idx)%coarse_graining_method = MASKED_AREA_WEIGHTED
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2  => Gfs_diag(nb)%soilm(:)
      Diag(idx)%data(nb)%var21 => Sfcprop(nb)%slmsk(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tmpmin2m'
    Diag(idx)%desc = 'min temperature at 2m height'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tmpmin(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tmpmax2m'
    Diag(idx)%desc = 'max temperature at 2m height'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tmpmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dusfc'
    Diag(idx)%desc = 'surface zonal momentum flux [N/m**2]'
    Diag(idx)%unit = 'N/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dusfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dvsfc'
    Diag(idx)%desc = 'surface meridional momentum flux'
    Diag(idx)%unit = 'N/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dvsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'shtfl_ave'
    Diag(idx)%desc = 'surface sensible heat flux'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dtsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'lhtfl_ave'
    Diag(idx)%desc = 'surface latent heat flux [W/m**2]'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dqsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totprcp_ave'
    Diag(idx)%desc = 'surface precipitation rate [kg/m**2/s]'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'full'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totprcp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totprcpb_ave'
    Diag(idx)%desc = 'bucket surface precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totprcpb(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'gflux_ave'
    Diag(idx)%desc = 'surface ground heat flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
!    Diag(idx)%mask = "land_ice_only"
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2  => Gfs_diag(nb)%gflux(:)
!      Diag(idx)%data(nb)%var21 => Sfcprop(nb)%slmsk(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DSWRF'
    Diag(idx)%desc = 'Interval-averaged zenith-angle-adjusted downward shortwave flux at the surface'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    ! This diagnostic is always meant to refer to what the model felt, so it
    ! refers to the override flux when overriding and the RRTMG flux when not.
    if (Model%override_surface_radiative_fluxes) then
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dswsfc_override(:)
      enddo
    else
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dswsfc(:)
      enddo
    endif

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'USWRF'
    Diag(idx)%desc = 'Interval-averaged zenith-angle-adjusted upward shortwave flux at the surface'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    ! This diagnostic is always meant to refer to what the model felt, so it
    ! refers to the override flux when overriding and the RRTMG flux when not.
    if (Model%override_surface_radiative_fluxes) then
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%uswsfc_override(:)
      enddo
    else
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%uswsfc(:)
      enddo
    endif

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DLWRF'
    Diag(idx)%desc = 'Interval-averaged surface-temperature-adjusted downward longwave flux at the surface'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    ! This diagnostic is always meant to refer to what the model felt, so it
    ! refers to the override flux when overriding and the RRTMG flux when not.
    if (Model%override_surface_radiative_fluxes) then
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dlwsfc_override(:)
      enddo
    else
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dlwsfc(:)
      enddo
    endif

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ULWRF'
    Diag(idx)%desc = 'Interval-averaged surface-temperature-adjusted upward longwave flux at the surface'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ulwsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sunsd_acc'
    Diag(idx)%desc = 'sunshine duration time'
    Diag(idx)%unit = 's'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%suntim(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'watr_acc'
    Diag(idx)%desc = 'total water runoff'
    Diag(idx)%unit = 'kg/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%runoff(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'pevpr_ave'
    Diag(idx)%desc = 'averaged potential evaporation rate'
    Diag(idx)%unit = 'W/M**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ep(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cwork_ave'
    Diag(idx)%desc = 'cloud work function (valid only with sas)'
    Diag(idx)%unit = 'J/kg'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cldwrk(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'u-gwd_ave'
    Diag(idx)%desc = 'surface zonal gravity wave stress'
    Diag(idx)%unit = 'N/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dugwd(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'v-gwd_ave'
    Diag(idx)%desc = 'surface meridional gravity wave stress'
    Diag(idx)%unit = 'N/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dvgwd(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'psmean'
    Diag(idx)%desc = 'surface pressure'
    Diag(idx)%unit = 'kPa'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%psmean(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cnvprcp_ave'
    Diag(idx)%desc = 'averaged surface convective precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'full'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cnvprcp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cnvprcpb_ave'
    Diag(idx)%desc = 'averaged bucket surface convective precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cnvprcpb(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cnvprcp'
    Diag(idx)%desc = 'surface convective precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cnvprcp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'spfhmin2m'
    Diag(idx)%desc = 'minimum specific humidity at 2m height'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%spfhmin(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'spfhmax2m'
    Diag(idx)%desc = 'maximum specific humidity at 2m height'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%spfhmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'u10mmax'
    Diag(idx)%desc = 'maximum (magnitude) u-wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'vector_bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%u10mmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'v10mmax'
    Diag(idx)%desc = 'maximum (magnitude) v-wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'vector_bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%v10mmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'wind10mmax'
    Diag(idx)%desc = 'maximum wind speed'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%wind10mmax(:)
    enddo

!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'u10max'
!    Diag(idx)%desc = 'hourly maximum (magnitude) u-wind'
!    Diag(idx)%unit = 'm/s'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'vector_bilinear'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%u10max(:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'v10max'
!    Diag(idx)%desc = 'hourly maximum (magnitude) v-wind'
!    Diag(idx)%unit = 'm/s'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'vector_bilinear'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%v10max(:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'spd10max'
!    Diag(idx)%desc = 'hourly maximum wind speed'
!    Diag(idx)%unit = 'm/s'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'bilinear'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%spd10max(:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 't02max'
!    Diag(idx)%desc = 'max hourly 2m Temperature'
!    Diag(idx)%unit = 'K'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'bilinear'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%t02max(:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 't02min'
!    Diag(idx)%desc = 'min hourly 2m Temperature'
!    Diag(idx)%unit = 'K'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'bilinear'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%t02min(:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'rh02max'
!    Diag(idx)%desc = 'max hourly 2m RH'
!    Diag(idx)%unit = '%'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'bilinear'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%rh02max(:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'rh02min'
!    Diag(idx)%desc = 'min hourly 2m RH'
!    Diag(idx)%unit = '%'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'bilinear'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%rh02min(:)
!    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'rain'
    Diag(idx)%desc = 'total rain at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%rain(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'rainc'
    Diag(idx)%desc = 'convective rain at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%rainc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ice'
    Diag(idx)%desc = 'ice fall at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snow'
    Diag(idx)%desc = 'snow fall at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%snow(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'graupel'
    Diag(idx)%desc = 'graupel fall at this time step'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%graupel(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totice_ave'
    Diag(idx)%desc = 'surface ice precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'full'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'toticeb_ave'
    Diag(idx)%desc = 'bucket surface ice precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%toticeb(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totsnw_ave'
    Diag(idx)%desc = 'surface snow precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'full'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totsnw(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totsnwb_ave'
    Diag(idx)%desc = 'bucket surface snow precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totsnwb(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totgrp_ave'
    Diag(idx)%desc = 'surface graupel precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%time_avg_kind = 'full'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totgrp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'totgrpb_ave'
    Diag(idx)%desc = 'bucket surface graupel precipitation rate'
    Diag(idx)%unit = 'kg/m**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_th
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%totgrpb(:)
    enddo

!--- physics instantaneous diagnostics ---
    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'u10m'
    Diag(idx)%desc = '10 meter u wind [m/s]'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'vector_bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%u10m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'v10m'
    Diag(idx)%desc = '10 meter v wind [m/s]'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'vector_bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%v10m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dpt2m'
    Diag(idx)%desc = '2 meter dew point temperature [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dpt2m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hgt_hyblev1'
    Diag(idx)%desc = 'layer 1 height'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%zlvl(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'psurf'
    Diag(idx)%desc = 'surface pressure [Pa]'
    Diag(idx)%unit = 'Pa'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%psurf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hpbl'
    Diag(idx)%desc = 'surface planetary boundary layer height [m]'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%hpbl(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hgamt'
    Diag(idx)%desc = 'ysu counter-gradient heat flux factor'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%hgamt(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hfxpbl'
    Diag(idx)%desc = 'ysu entrainment heat flux factor'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%hfxpbl(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'xmb_shal'
    Diag(idx)%desc = 'cloud base mass flux from mass-flux shal cnv'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%xmb_shal(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tfac_shal'
    Diag(idx)%desc = 'Tadv/Tcnv factor from  mass-flux shal cnv'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tfac_shal(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sigma_shal'
    Diag(idx)%desc = 'updraft fractional area from mass-flux shal cnv'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%sigma_shal(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'pwat'
    Diag(idx)%desc = 'atmos column precipitable water [kg/m**2]'
    Diag(idx)%unit = 'kg/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%pwat(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tmp_hyblev1'
    Diag(idx)%desc = 'layer 1 temperature'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%t1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'spfh_hyblev1'
    Diag(idx)%desc = 'layer 1 specific humidity'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%q1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ugrd_hyblev1'
    Diag(idx)%desc = 'layer 1 zonal wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'vector_bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%u1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'vgrd_hyblev1'
    Diag(idx)%desc = 'layer 1 meridional wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'vector_bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%v1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sfexc'
    Diag(idx)%desc = 'Exchange Coefficient'
    Diag(idx)%unit = 'kg/m2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%chh(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'acond'
    Diag(idx)%desc = 'Aerodynamic conductance'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cmm(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DLWRFI'
    Diag(idx)%desc = 'Instantaneous surface-temperature-adjusted downward longwave flux at the surface'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    ! This diagnostic is always meant to refer to what the model felt, so it
    ! refers to the override flux when overriding and the RRTMG flux when not.
    if (Model%override_surface_radiative_fluxes) then
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dlwsfci_override(:)
      enddo
    else
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dlwsfci(:)
      enddo
    endif

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ULWRFI'
    Diag(idx)%desc = 'Instantaneous surface-temperature-adjusted upward longwave flux at the surface'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ulwsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'DSWRFI'
    Diag(idx)%desc = 'Instantaneous zenith-angle-adjusted downward shortwave flux at the surface'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    ! This diagnostic is always meant to refer to what the model felt, so it
    ! refers to the override flux when overriding and the RRTMG flux when not.
    if (Model%override_surface_radiative_fluxes) then
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dswsfci_override(:)
      enddo
    else
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dswsfci(:)
      enddo
    endif

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'USWRFI'
    Diag(idx)%desc = 'Instantaneous zenith-angle-adjusted upward shortwave flux at the surface'
    Diag(idx)%unit = 'w/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    ! This diagnostic is always meant to refer to what the model felt, so it
    ! refers to the override flux when overriding and the RRTMG flux when not.
    if (Model%override_surface_radiative_fluxes) then
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%uswsfci_override(:)
      enddo
    else
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%uswsfci(:)
      enddo
    endif

    if (Model%override_surface_radiative_fluxes) then
      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'DLWRF_from_rrtmg'
      Diag(idx)%desc = 'Interval-averaged native RRTMG surface-temperature-adjusted downward longwave flux at the surface'
      Diag(idx)%unit = 'W/m**2'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%cnvfac = cn_one
      Diag(idx)%time_avg = .TRUE.
      Diag(idx)%intpl_method = 'bilinear'
      Diag(idx)%coarse_graining_method = 'area_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dlwsfc(:)
      enddo

      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'DLWRFI_from_rrtmg'
      Diag(idx)%desc = 'Instantaneous native RRTMG surface-temperature-adjusted downward longwave flux at the surface'
      Diag(idx)%unit = 'W/m**2'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%intpl_method = 'bilinear'
      Diag(idx)%coarse_graining_method = 'area_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dlwsfci(:)
      enddo

      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'ULWRF_from_rrtmg'
      Diag(idx)%desc = 'Interval-averaged native RRTMG surface-temperature-adjusted upward longwave flux at the surface'
      Diag(idx)%unit = 'W/m**2'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%cnvfac = cn_one
      Diag(idx)%time_avg = .TRUE.
      Diag(idx)%intpl_method = 'bilinear'
      Diag(idx)%coarse_graining_method = 'area_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ulwsfc(:)
      enddo

      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'ULWRFI_from_rrtmg'
      Diag(idx)%desc = 'Instantaneous native RRTMG surface-temperature-adjusted upward longwave flux at the surface'
      Diag(idx)%unit = 'W/m**2'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%intpl_method = 'bilinear'
      Diag(idx)%coarse_graining_method = 'area_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ulwsfci(:)
      enddo

      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'DSWRF_from_rrtmg'
      Diag(idx)%desc = 'Interval-averaged native RRTMG zenith-angle-adjusted downward shortwave flux at the surface'
      Diag(idx)%unit = 'W/m**2'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%cnvfac = cn_one
      Diag(idx)%time_avg = .TRUE.
      Diag(idx)%intpl_method = 'bilinear'
      Diag(idx)%coarse_graining_method = 'area_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dswsfc(:)
      enddo

      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'DSWRFI_from_rrtmg'
      Diag(idx)%desc = 'Instantaneous native RRTMG zenith-angle-adjusted downward shortwave flux at the surface'
      Diag(idx)%unit = 'W/m**2'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%intpl_method = 'bilinear'
      Diag(idx)%coarse_graining_method = 'area_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dswsfci(:)
      enddo

      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'USWRF_from_rrtmg'
      Diag(idx)%desc = 'Interval-averaged native RRTMG zenith-angle-adjusted upward shortwave flux at the surface'
      Diag(idx)%unit = 'W/m**2'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%cnvfac = cn_one
      Diag(idx)%time_avg = .TRUE.
      Diag(idx)%intpl_method = 'bilinear'
      Diag(idx)%coarse_graining_method = 'area_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%uswsfc(:)
      enddo

      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'USWRFI_from_rrtmg'
      Diag(idx)%desc = 'Instantaneous native RRTMG zenith-angle-adjusted upward shortwave flux at the surface'
      Diag(idx)%unit = 'W/m**2'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%intpl_method = 'bilinear'
      Diag(idx)%coarse_graining_method = 'area_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%uswsfci(:)
      enddo
    endif

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dusfci'
    Diag(idx)%desc = 'instantaneous u component of surface stress'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dusfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'dvsfci'
    Diag(idx)%desc = 'instantaneous v component of surface stress'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dvsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'shtfl'
    Diag(idx)%desc = 'instantaneous surface sensible heat flux'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dtsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'lhtfl'
    Diag(idx)%desc = 'instantaneous surface latent heat flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%dqsfci(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'gfluxi'
    Diag(idx)%desc = 'instantaneous surface ground heat flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%gfluxi(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'pevpr'
    Diag(idx)%desc = 'instantaneous surface potential evaporation'
    Diag(idx)%unit = 'W/M**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%epi(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'wilt'
    Diag(idx)%desc = 'wiltimg point (volumetric)'
    Diag(idx)%unit = 'Proportion'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%smcwlt2(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'fldcp'
    Diag(idx)%desc = 'Field Capacity (volumetric)'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%smcref2(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'wet1'
    Diag(idx)%desc = 'normalized soil wetness'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%wet1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cpofp'
    Diag(idx)%desc = 'Percent frozen precipitation'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%intpl_method = 'bilinear'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%sr(:)
    enddo

#if defined (USE_COSP)
!--- 2D diagnostic variables from the CFMIP Observation Simulator Package (COSP), Linjiong Zhou

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cltisccp'
    Diag(idx)%desc = 'ISCCP Total Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cltisccp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'meantbisccp'
    Diag(idx)%desc = 'ISCCP all-sky 10.5 micron brightness temperature / toa_brightness_temperature'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%meantbisccp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'meantbclrisccp'
    Diag(idx)%desc = 'ISCCP clear-sky 10.5 micron brightness temperature / toa_brightness_temperature_assuming_clear_sky'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%meantbclrisccp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'pctisccp'
    Diag(idx)%desc = 'ISCCP Mean Cloud Top Pressure / air_pressure_at_cloud_top'
    Diag(idx)%unit = 'hPa'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%pctisccp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tauisccp'
    Diag(idx)%desc = 'ISCCP Mean Optical Depth / atmosphere_optical_thickness_due_to_cloud'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%tauisccp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'albisccp'
    Diag(idx)%desc = 'ISCCP Mean Cloud Albedo / cloud_albedo'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%albisccp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'misr_meanztop'
    Diag(idx)%desc = 'MISR Mean Cloud Top Height / cloud_top_altitude'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%misr_meanztop(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'misr_cldarea'
    Diag(idx)%desc = 'MISR cloud cover / cloud_area_fraction'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%misr_cldarea(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cltmodis'
    Diag(idx)%desc = 'MODIS Total Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cltmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clwmodis'
    Diag(idx)%desc = 'MODIS Liquid Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clwmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'climodis'
    Diag(idx)%desc = 'MODIS Ice Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%climodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clhmodis'
    Diag(idx)%desc = 'MODIS High Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clhmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clmmodis'
    Diag(idx)%desc = 'MODIS Mid Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clmmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cllmodis'
    Diag(idx)%desc = 'MODIS Low Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cllmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tautmodis'
    Diag(idx)%desc = 'MODIS Total Cloud Optical Thickness / atmosphere_optical_thickness_due_to_cloud'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%tautmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tauwmodis'
    Diag(idx)%desc = 'MODIS Liquid Cloud Optical Thickness / atmosphere_optical_thickness_due_to_cloud'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%tauwmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tauimodis'
    Diag(idx)%desc = 'MODIS Ice Cloud Optical Thickness / atmosphere_optical_thickness_due_to_cloud'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%tauimodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tautlogmodis'
    Diag(idx)%desc = 'MODIS Total Cloud Optical Thickness (Log10 Mean) / atmosphere_optical_thickness_due_to_cloud'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%tautlogmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tauwlogmodis'
    Diag(idx)%desc = 'MODIS Liquid Cloud Optical Thickness (Log10 Mean) / atmosphere_optical_thickness_due_to_cloud'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%tauwlogmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tauilogmodis'
    Diag(idx)%desc = 'MODIS Ice Cloud Optical Thickness (Log10 Mean) / atmosphere_optical_thickness_due_to_cloud'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%tauilogmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'reffclwmodis'
    Diag(idx)%desc = 'MODIS Liquid Cloud Particle Size / effective_radius_of_cloud_liquid_water_particle'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%reffclwmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'reffclimodis'
    Diag(idx)%desc = 'MODIS Ice Cloud Particle Size / effective_radius_of_cloud_liquid_water_particle'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%reffclimodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'pctmodis'
    Diag(idx)%desc = 'MODIS Cloud Top Pressure / air_pressure_at_cloud_top'
    Diag(idx)%unit = 'hPa'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%pctmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'lwpmodis'
    Diag(idx)%desc = 'MODIS Cloud Liquid Water Path / atmosphere_cloud_liquid_water_content'
    Diag(idx)%unit = 'kg m-2'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%lwpmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'iwpmodis'
    Diag(idx)%desc = 'MODIS Cloud Ice Water Path / atmosphere_mass_content_of_cloud_ice'
    Diag(idx)%unit = 'kg m-2'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%iwpmodis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cltlidarradar'
    Diag(idx)%desc = 'CALIPSO and CloudSat Total Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cltlidarradar(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cllcalipsoice'
    Diag(idx)%desc = 'CALIPSO Ice Low Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cllcalipsoice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clmcalipsoice'
    Diag(idx)%desc = 'CALIPSO Ice Mid Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clmcalipsoice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clhcalipsoice'
    Diag(idx)%desc = 'CALIPSO Ice High Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clhcalipsoice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cltcalipsoice'
    Diag(idx)%desc = 'CALIPSO Ice Total Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cltcalipsoice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cllcalipsoliq'
    Diag(idx)%desc = 'CALIPSO Liquid Low Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cllcalipsoliq(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clmcalipsoliq'
    Diag(idx)%desc = 'CALIPSO Liquid Mid Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clmcalipsoliq(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clhcalipsoliq'
    Diag(idx)%desc = 'CALIPSO Liquid High Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clhcalipsoliq(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cltcalipsoliq'
    Diag(idx)%desc = 'CALIPSO Liquid Total Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cltcalipsoliq(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cllcalipsoun'
    Diag(idx)%desc = 'CALIPSO Undefined-Phase Low Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cllcalipsoun(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clmcalipsoun'
    Diag(idx)%desc = 'CALIPSO Undefined-Phase Mid Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clmcalipsoun(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clhcalipsoun'
    Diag(idx)%desc = 'CALIPSO Undefined-Phase High Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clhcalipsoun(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cltcalipsoun'
    Diag(idx)%desc = 'CALIPSO Undefined-Phase Total Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cltcalipsoun(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cllcalipso'
    Diag(idx)%desc = 'CALIPSO Low Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cllcalipso(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clmcalipso'
    Diag(idx)%desc = 'CALIPSO Mid Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clmcalipso(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clhcalipso'
    Diag(idx)%desc = 'CALIPSO High Level Cloud Fraction / cloud_area_fraction_in_atmosphere_layer'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clhcalipso(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cltcalipso'
    Diag(idx)%desc = 'CALIPSO Total Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cltcalipso(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clopaquecalipso'
    Diag(idx)%desc = 'CALIPSO Opaque Cloud Cover / opaque_cloud_cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clopaquecalipso(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clthincalipso'
    Diag(idx)%desc = 'CALIPSO Thin Cloud Cover / thin_cloud_cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clthincalipso(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clzopaquecalipso'
    Diag(idx)%desc = 'CALIPSO z_opaque Altitude / z_opaque'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clzopaquecalipso(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clopaquetemp'
    Diag(idx)%desc = 'CALIPSO Opaque Cloud Temperature / opaque_cloud_temperature'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clopaquetemp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clthintemp'
    Diag(idx)%desc = 'CALIPSO Thin Cloud Temperature / thin_cloud_temperature'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clthintemp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clzopaquetemp'
    Diag(idx)%desc = 'CALIPSO z_opaque Temperature / z_opaque_temperature'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clzopaquetemp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clopaquemeanz'
    Diag(idx)%desc = 'CALIPSO Opaque Cloud Altitude / opaque_cloud_altitude'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clopaquemeanz(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clthinmeanz'
    Diag(idx)%desc = 'CALIPSO Thin Cloud Altitude / thin_cloud_altitude'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clthinmeanz(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clthinemis'
    Diag(idx)%desc = 'CALIPSO Thin Cloud Emissivity / thin_cloud_emissivity'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clthinemis(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clopaquemeanzse'
    Diag(idx)%desc = 'CALIPSO Opaque Cloud Altitude with respect to SE / opaque_cloud_altitude_se'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clopaquemeanzse(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clthinmeanzse'
    Diag(idx)%desc = 'CALIPSO Thin Cloud Altitude with respect to SE / thin_cloud_altitude_se'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clthinmeanzse(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clzopaquecalipsose'
    Diag(idx)%desc = 'CALIPSO z_opaque Altitude with respect to SE / z_opaque_se'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clzopaquecalipsose(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cllgrLidar532'
    Diag(idx)%desc = 'GROUND LIDAR Low Level Cloud Cover / grLidar532_low_cloud_cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cllgrLidar532(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clmgrLidar532'
    Diag(idx)%desc = 'GROUND LIDAR Mid Level Cloud Cover / grLidar532_mid_cloud_cover'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clmgrLidar532(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clhgrLidar532'
    Diag(idx)%desc = 'GROUND LIDAR High Level Cloud Cover / grLidar532_high_cloud_cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clhgrLidar532(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cltgrLidar532'
    Diag(idx)%desc = 'GROUND LIDAR Total Cloud Cover / grLidar532_total_cloud_cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cltgrLidar532(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cllatlid'
    Diag(idx)%desc = 'ATLID Low Level Cloud Cover / atlid_low_cloud_cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cllatlid(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clmatlid'
    Diag(idx)%desc = 'ATLID Mid Level Cloud Cover / atlid_mid_cloud_cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clmatlid(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'clhatlid'
    Diag(idx)%desc = 'ATLID High Level Cloud Cover /  atlid_high_cloud_cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%clhatlid(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cltatlid'
    Diag(idx)%desc = 'ATLID Total Cloud Cover / atlid_total_cloud_cover'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cltatlid(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag0'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag0'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag0(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag1'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag1'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag1(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag2'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag2'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag2(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag3'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag3'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag3(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag4'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag4'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag4(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag5'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag5'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag5(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag6'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag6'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag6(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag7'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag7'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag7(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag8'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag8'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag8(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ptcloudsatflag9'
    Diag(idx)%desc = 'Cloudsat precipitation cover for flag9'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%ptcloudsatflag9(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cloudsatpia'
    Diag(idx)%desc = 'Cloudsat path integrated attenuation'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cloudsatpia(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cloudsat_tcc'
    Diag(idx)%desc = 'CloudSat Total Cloud Fraction / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cloudsat_tcc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'cloudsat_tcc2'
    Diag(idx)%desc = 'CloudSat Total Cloud Fraction (no 1km) / cloud_area_fraction'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%cloudsat_tcc2(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'npdfcld'
    Diag(idx)%desc = '# of Non-Precipitating Clouds / number_of_slwc_nonprecip'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%npdfcld(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'npdfdrz'
    Diag(idx)%desc = '# of Drizzling Clouds / number_of_slwc_drizzle'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%npdfdrz(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'npdfrain'
    Diag(idx)%desc = '# of Precipitating Clouds / number_of_slwc_precip'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%npdfrain(:)
    enddo
#endif

#if defined (COSP_OFFLINE)
!--- 2D/3D variables for the offline CFMIP Observation Simulator Package (COSP), Linjiong Zhou

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'skt'
    Diag(idx)%desc = 'Skin temperature'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%skt(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'surfelev'
    Diag(idx)%desc = 'Surface Elevation'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%surfelev(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'landmask'
    Diag(idx)%desc = 'Land/sea mask'
    Diag(idx)%unit = '0/1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%landmask(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'sunlit'
    Diag(idx)%desc = 'Sunlit flag'
    Diag(idx)%unit = 'none'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%cosp%sunlit(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'p'
    Diag(idx)%desc = 'Model pressure levels'
    Diag(idx)%unit = 'pa'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%p(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'ph'
    Diag(idx)%desc = 'Moddel pressure at half levels'
    Diag(idx)%unit = 'pa'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%ph(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'zlev'
    Diag(idx)%desc = 'Model level height'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%zlev(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'zlev_half'
    Diag(idx)%desc = 'Model level height at half-levels'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%zlev_half(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'T'
    Diag(idx)%desc = 'Temperature'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%T(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'sh'
    Diag(idx)%desc = 'Specific humidity'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%sh(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'tca'
    Diag(idx)%desc = 'Total cloud fraction'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%tca(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'cca'
    Diag(idx)%desc = 'Convective cloud fraction'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%cca(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'u_wind'
    Diag(idx)%desc = 'U-component of wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%u_wind(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'v_wind'
    Diag(idx)%desc = 'V-component of wind'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%v_wind(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'mr_lsliq'
    Diag(idx)%desc = 'Mass mixing ratio for stratiform cloud liquid'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%mr_lsliq(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'mr_lsice'
    Diag(idx)%desc = 'Mass mixing ratio for stratiform cloud ice'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%mr_lsice(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'mr_ccliq'
    Diag(idx)%desc = 'Mass mixing ratio for convective cloud liquid'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%mr_ccliq(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'mr_ccice'
    Diag(idx)%desc = 'Mass mixing ratio for convective cloud ice'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%mr_ccice(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'mr_ozone'
    Diag(idx)%desc = 'Mass mixing ratio for ozone'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%mr_ozone(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'fl_lsrain'
    Diag(idx)%desc = 'Precipitation flux (rain) for stratiform cloud'
    Diag(idx)%unit = 'kg/m^2/s'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%fl_lsrain(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'fl_lssnow'
    Diag(idx)%desc = 'Precipitation flux (snow) for stratiform cloud'
    Diag(idx)%unit = 'kg/m^2/s'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%fl_lssnow(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'fl_lsgrpl'
    Diag(idx)%desc = 'Precipitation flux (groupel) for stratiform cloud'
    Diag(idx)%unit = 'kg/m^2/s'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%fl_lsgrpl(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'fl_ccrain'
    Diag(idx)%desc = 'Precipitation flux (rain) for convective cloud'
    Diag(idx)%unit = 'kg/m^2/s'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%fl_ccrain(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'fl_ccsnow'
    Diag(idx)%desc = 'Precipitation flux (snow) for convective cloud'
    Diag(idx)%unit = 'kg/m^2/s'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%fl_ccsnow(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'dtau_s'
    Diag(idx)%desc = '0.67micron optical depth (stratiform cloud)'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%dtau_s(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'dtau_c'
    Diag(idx)%desc = '0.67micron optical depth (convective cloud)'
    Diag(idx)%unit = '1'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%dtau_c(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'dem_s'
    Diag(idx)%desc = '11micron emissivity (stratiform cloud)'
    Diag(idx)%unit = 'none'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%dem_s(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'dem_c'
    Diag(idx)%desc = '11microm emissivity (convective cloud)'
    Diag(idx)%unit = 'none'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%dem_c(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'Reff_LSCLIQ'
    Diag(idx)%desc = 'Subcolumn effective radius for stratiform cloud liquid'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%Reff_LSCLIQ(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'Reff_LSCICE'
    Diag(idx)%desc = 'Subcolumn effective radius for stratiform cloud ice'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%Reff_LSCICE(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'Reff_LSRAIN'
    Diag(idx)%desc = 'Subcolumn effective radius for stratiform rain'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%Reff_LSRAIN(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'Reff_LSSNOW'
    Diag(idx)%desc = 'Subcolumn effective radius for stratiform snow'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%Reff_LSSNOW(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'Reff_LSGRPL'
    Diag(idx)%desc = 'Subcolumn effective radius for stratiform graupel'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'cosp'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%cosp%Reff_LSGRPL(:,:)
    enddo

#endif

!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'crain_ave'
!    Diag(idx)%desc = 'averaged categorical rain'
!    Diag(idx)%unit = 'number'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'bilinear'
!    Diag(idx)%cnvfac = cn_one
!    Diag(idx)%time_avg = .TRUE.
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tdomr(:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'csnow_ave'
!    Diag(idx)%desc = 'averaged categorical snow'
!    Diag(idx)%unit = 'number'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'bilinear'
!    Diag(idx)%cnvfac = cn_one
!    Diag(idx)%time_avg = .TRUE.
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tdoms(:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'cfrzr_ave'
!    Diag(idx)%desc = 'averaged categorical freezing rain'
!    Diag(idx)%unit = 'number'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'bilinear'
!    Diag(idx)%cnvfac = cn_one
!    Diag(idx)%time_avg = .TRUE.
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tdomzr(:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'cicep_ave'
!    Diag(idx)%desc = 'averaged categorical sleet'
!    Diag(idx)%unit = 'number'
!    Diag(idx)%mod_name = 'gfs_phys'
!    Diag(idx)%intpl_method = 'bilinear'
!    Diag(idx)%cnvfac = cn_one
!    Diag(idx)%time_avg = .TRUE.
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tdomip(:)
!    enddo
!
!--- three-dimensional variables that need to be handled special when writing
    if (Model%ldiag3d) then

    do num = 1,6
      write (xtra,'(I1)') num
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'dt3dt_'//trim(xtra)
      Diag(idx)%desc = 'temperature change due to physics '//trim(xtra)//''
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      Diag(idx)%coarse_graining_method = 'mass_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%dt3dt(:,:,num)
      enddo
    enddo

    do num = 1,5+oz_coeff
      write (xtra,'(I1)') num
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'dq3dt_'//trim(xtra)
      Diag(idx)%desc = 'moisture change due to physics '//trim(xtra)//''
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      Diag(idx)%coarse_graining_method = 'mass_weighted'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%dq3dt(:,:,num)
      enddo
    enddo

    do num = 1,4
      write (xtra,'(I1)') num
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'du3dt_'//trim(xtra)
      Diag(idx)%desc = 'u momentum change due to physics '//trim(xtra)//''
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%du3dt(:,:,num)
      enddo
    enddo

    do num = 1,4
      write (xtra,'(I1)') num
      idx = idx + 1
      Diag(idx)%axes = 3
      Diag(idx)%name = 'dv3dt_'//trim(xtra)
      Diag(idx)%desc = 'v momentum change due to physics '//trim(xtra)//''
      Diag(idx)%unit = 'XXX'
      Diag(idx)%mod_name = 'gfs_phys'
      Diag(idx)%time_avg = .TRUE.
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
         Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%dv3dt(:,:,num)
      enddo
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'dkt_pbl'
    Diag(idx)%desc = 'instantaneous heat diffusion coefficient'
    Diag(idx)%unit = 'm**2/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%dkt(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'flux_cg'
    Diag(idx)%desc = 'instantaneous counter-gradient heat flux in ysu'
    Diag(idx)%unit = 'K*m/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%flux_cg(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'flux_en'
    Diag(idx)%desc = 'instantaneous entrainment heat flux in ysu'
    Diag(idx)%unit = 'K*m/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%flux_en(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'wu2_shal'
    Diag(idx)%desc = 'updraft velocity square from shallow convection'
    Diag(idx)%unit = 'm**2/s**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%wu2_shal(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'eta_shal'
    Diag(idx)%desc = 'normalized mass flux from shallow convection'
    Diag(idx)%unit = 'non-dim'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%eta_shal(:,:)
    enddo

!    idx = idx + 1
!    Diag(idx)%axes = 3
!    Diag(idx)%name = 'refl_10cm'
!    Diag(idx)%desc = 'Radar reflectivity'
!    Diag(idx)%unit = 'dBz'
!    Diag(idx)%mod_name = 'gfs_phys'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%refl_10cm(:,:)
!    enddo
!
!    idx = idx + 1
!    Diag(idx)%axes = 3
!    Diag(idx)%name = 'cnvw'
!    Diag(idx)%desc = 'subgrid scale convective cloud water'
!    Diag(idx)%unit = 'kg/kg'
!    Diag(idx)%mod_name = 'gfs_phys'
!    allocate (Diag(idx)%data(nblks))
!    if( Model%ncnvw > 0 ) then
!      do nb = 1,nblks
!        Diag(idx)%data(nb)%var3 => Tbd(nb)%phy_f3d(:,:,Model%ncnvw)
!      enddo
!    endif

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'diss_est'
    Diag(idx)%desc = 'dissipation rate for skeb'
    Diag(idx)%unit = 'none'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%diss_est(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'skebu_wts'
    Diag(idx)%desc = 'perturbation velocity'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%skebu_wts(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'skebv_wts'
    Diag(idx)%desc = 'perturbation velocity'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%skebv_wts(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'zmtnblck'
    Diag(idx)%desc = 'level of dividing streamline'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%zmtnblck(:)
    enddo

!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'refdmax'
!    Diag(idx)%desc = 'max hourly 1-km agl reflectivity'
!    Diag(idx)%unit = 'dBZ'
!    Diag(idx)%mod_name = 'gfs_phys'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%refdmax(:)
!    enddo
!    idx = idx + 1
!    Diag(idx)%axes = 2
!    Diag(idx)%name = 'refdmax263k'
!    Diag(idx)%desc = 'max hourly -10C reflectivity'
!    Diag(idx)%unit = 'dBZ'
!    Diag(idx)%mod_name = 'gfs_phys'
!    allocate (Diag(idx)%data(nblks))
!    do nb = 1,nblks
!      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%refdmax263k(:)
!    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'sppt_wts'
    Diag(idx)%desc = 'perturbation velocity'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%sppt_wts(:,:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 3
    Diag(idx)%name = 'shum_wts'
    Diag(idx)%desc = 'perturbation velocity'
    Diag(idx)%unit = 'm/s'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%shum_wts(:,:)
    enddo

!!$    idx = idx + 1
!!$    Diag(idx)%axes = 3
!!$    !Requires lgocart = .T.
!!$    Diag(idx)%name = 'dqdt_v'
!!$    Diag(idx)%desc = 'instantaneous total moisture tendency'
!!$    Diag(idx)%unit = 'XXX'
!!$    Diag(idx)%mod_name = 'gfs_phys'
!!$    allocate (Diag(idx)%data(nblks))
!!$    do nb = 1,nblks
!!$       Diag(idx)%data(nb)%var3 => Gfs_diag(nb)%Diag(nb)%dqdt_v(:,:,num)
!!$    enddo

!
!--- prognostic variable tendencies (T, u, v, sph, clwmr, o3)
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'dtemp_dt'
!rab    Diag(idx)%desc = 'GFS radiation/physics temperature tendency'
!rab    Diag(idx)%unit = 'K/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'du_dt'
!rab    Diag(idx)%desc = 'GFS radiation/physics horizontal wind component tendency'
!rab    Diag(idx)%unit = 'm/s/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'dv_dt'
!rab    Diag(idx)%desc = 'GFS radiation/physics meridional wind component tendency'
!rab    Diag(idx)%unit = 'm/s/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'dsphum_dt'
!rab    Diag(idx)%desc = 'GFS radiation/physics specific humidity tendency'
!rab    Diag(idx)%unit = 'kg/kg/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'dclwmr_dt'
!rab    Diag(idx)%desc = 'GFS radiation/radiation cloud water mixing ratio tendency'
!rab    Diag(idx)%unit = 'kg/kg/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'
!rab
!rab    idx = idx + 1
!rab    Diag(idx)%axes = 3
!rab    Diag(idx)%name = 'do3mr_dt'
!rab    Diag(idx)%desc = 'GFS radiation/radiation ozone mixing ratio tendency'
!rab    Diag(idx)%unit = 'kg/kg/s'
!rab    Diag(idx)%mod_name = 'gfs_phys'

    endif

!--- Surface diagnostics in gfs_sfc
    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'alnsf'
    Diag(idx)%desc = 'mean nir albedo with strong cosz dependency'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%alnsf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'alnwf'
    Diag(idx)%desc = 'mean nir albedo with weak cosz dependency'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%alnwf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'alvsf'
    Diag(idx)%desc = 'mean vis albedo with strong cosz dependency'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%alvsf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'alvwf'
    Diag(idx)%desc = 'mean vis albedo with weak cosz dependency'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%alvwf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'canopy'
    Diag(idx)%desc = 'canopy water (cnwat in gfs data)'
    Diag(idx)%unit = '%'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%canopy(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'f10m'
    Diag(idx)%desc = '10-meter wind speed divided by lowest model wind speed'
    Diag(idx)%unit = 'N/A'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%f10m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'facsf'
    Diag(idx)%desc = 'fractional coverage with strong cosz dependency'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%facsf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'facwf'
    Diag(idx)%desc = 'fractional coverage with weak cosz dependency'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%facwf(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ffhh'
    Diag(idx)%desc = 'fh parameter from PBL scheme'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%ffhh(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ffmm'
    Diag(idx)%desc = 'fm parameter from PBL scheme'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%ffmm(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'uustar'
    Diag(idx)%desc = 'uustar surface frictional wind'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%uustar(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'slope'
    Diag(idx)%desc = 'surface slope type'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%slope(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'fice'
    Diag(idx)%desc = 'surface ice concentration (ice=1; no ice=0) [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%fice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'hice'
    Diag(idx)%desc = 'sea ice thickness (icetk in gfs_data)'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%hice(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snoalb'
    Diag(idx)%desc = 'maximum snow albedo in fraction (salbd?? in gfs data)'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%snoalb(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'shdmax'
    Diag(idx)%desc = 'maximum fractional coverage of green vegetation'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%shdmax(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'shdmin'
    Diag(idx)%desc = 'minimum fractional coverage of green vegetation'
    Diag(idx)%unit = 'XXX'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%shdmin(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snowd'
    Diag(idx)%desc = 'surface snow depth [m]'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%cnvfac = cn_one/cn_th
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%snowd(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'snow_cover'
    Diag(idx)%desc = 'snow cover area fraction'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%sncovr(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'soil_color'
    Diag(idx)%desc = 'soil color category'
    Diag(idx)%unit = 'none'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%scolor(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'crain'
    Diag(idx)%desc = 'instantaneous categorical rain'
    Diag(idx)%unit = 'number'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%cnvfac = cn_one
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%srflag(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'stype'
    Diag(idx)%desc = 'soil type in integer 1-9'
    Diag(idx)%unit = 'N/A'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stype(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'q2m'
    Diag(idx)%desc = '2m specific humidity [kg/kg]'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%q2m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 't2m'
    Diag(idx)%desc = '2m temperature [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%t2m(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tsfc'
    Diag(idx)%desc = 'surface temperature [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%tsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'qsfc'
    Diag(idx)%desc = 'surface specific humidity [kg/kg]'
    Diag(idx)%unit = 'kg/kg'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%qsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tg3'
    Diag(idx)%desc = 'deep soil temperature'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%tg3(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tisfc'
    Diag(idx)%desc = 'surface temperature over ice fraction'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%tisfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tprcp'
    Diag(idx)%desc = 'total precipitation'
    Diag(idx)%unit = 'kg/m**2'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%tprcp(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'vtype'
    Diag(idx)%desc = 'vegetation type in integer 1-13'
    Diag(idx)%unit = 'number'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%vtype(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'weasd'
    Diag(idx)%desc = 'surface snow water equivalent [kg/m**2]'
    Diag(idx)%unit = 'kg/m**2'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%weasd(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'HGTsfc'
    Diag(idx)%desc = 'surface geopotential height [gpm]'
    Diag(idx)%unit = 'gpm'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%cnvfac = con_g
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%oro(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SLMSKsfc'
    Diag(idx)%desc = 'sea-land-ice mask (0-sea, 1-land, 2-ice)'
    Diag(idx)%unit = 'N/A'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%coarse_graining_method = MODE
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%slmsk(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ZORLsfc'
    Diag(idx)%desc = 'surface roughness [m]'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%cnvfac = cn_one/cn_100
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%zorl(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'VFRACsfc'
    Diag(idx)%desc = 'vegetation fraction'
    Diag(idx)%unit = 'N/A'
    Diag(idx)%mod_name = 'gfs_sfc'
    Diag(idx)%cnvfac = cn_100
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%vfrac(:)
    enddo

    do num = 1,4
      write (xtra,'(I1)') num
      idx = idx + 1
      Diag(idx)%axes = 2
      Diag(idx)%name = 'slc_'//trim(xtra)
      Diag(idx)%desc = 'liquid soil mositure at layer-'//trim(xtra)
      Diag(idx)%unit = 'xxx'
      Diag(idx)%mod_name = 'gfs_sfc'
      allocate (Diag(idx)%data(nblks))
      do nb = 1,nblks
        Diag(idx)%data(nb)%var2 => Sfcprop(nb)%slc(:,num)
      enddo
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILW1'
    Diag(idx)%desc = 'volumetric soil moisture 0-10cm [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%smc(:,1)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILW2'
    Diag(idx)%desc = 'volumetric soil moisture 10-40cm [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%smc(:,2)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILW3'
    Diag(idx)%desc = 'volumetric soil moisture 40-100cm [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%smc(:,3)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILW4'
    Diag(idx)%desc = 'volumetric soil moisture 100-200cm [fraction]'
    Diag(idx)%unit = 'fraction'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%smc(:,4)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILT1'
    Diag(idx)%desc = 'soil temperature 0-10cm [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stc(:,1)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILT2'
    Diag(idx)%desc = 'soil temperature 10-40cm [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stc(:,2)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILT3'
    Diag(idx)%desc = 'soil temperature 40-100cm [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stc(:,3)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'SOILT4'
    Diag(idx)%desc = 'soil temperature 100-200cm [K]'
    Diag(idx)%unit = 'K'
    Diag(idx)%mod_name = 'gfs_sfc'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Sfcprop(nb)%stc(:,4)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'netflxsfc'
    Diag(idx)%desc = 'net surface heat flux [W/m**2]'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%netflxsfc(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'qflux_restore'
    Diag(idx)%desc = 'restoring flux'
    Diag(idx)%unit = 'W/m**2'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%coarse_graining_method = AREA_WEIGHTED
    Diag(idx)%time_avg = .TRUE.
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%qflux_restore(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'tclim_iano'
    Diag(idx)%desc = 'climatological SST plus initial anomaly'
    Diag(idx)%unit = 'degree C'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%tclim_iano(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'MLD'
    Diag(idx)%desc = 'Interval-average ocean mixed layer depth'
    Diag(idx)%unit = 'm'
    Diag(idx)%mod_name = 'gfs_phys'
    Diag(idx)%cnvfac = cn_one
    Diag(idx)%time_avg = .TRUE.
    Diag(idx)%coarse_graining_method = AREA_WEIGHTED
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%mld(:)
    enddo

    idx = idx + 1
    Diag(idx)%axes = 2
    Diag(idx)%name = 'ps_dt'
    Diag(idx)%desc = 'surface pressure tendency'
    Diag(idx)%unit = 'Pa/3hr'
    Diag(idx)%mod_name = 'gfs_phys'
    allocate (Diag(idx)%data(nblks))
    do nb = 1,nblks
      Diag(idx)%data(nb)%var2 => Gfs_diag(nb)%ps_dt(:)
    enddo

    tot_diag_idx = idx

    if (idx > DIAG_SIZE) then
      call mpp_error(FATAL, 'gfs_driver::gfs_diag_register - need to increase DIAG_SIZE')
    endif

    do idx = 1,tot_diag_idx
      if (diag(idx)%axes == -99) then
        call mpp_error(FATAL, 'gfs_driver::gfs_diag_register - attempt to register an undefined variable')
      endif
      Diag(idx)%id = register_diag_field (trim(Diag(idx)%mod_name), trim(Diag(idx)%name),  &
                                          axes(1:Diag(idx)%axes), Time, trim(Diag(idx)%desc), &
                                          trim(Diag(idx)%unit), missing_value=real(missing_value))
    enddo

  end subroutine gfdl_diag_register

  subroutine populate_coarse_diag_type(diagnostic, coarse_diagnostic)
    type(gfdl_diag_type), intent(in) :: diagnostic
    type(gfdl_diag_type), intent(inout) :: coarse_diagnostic

    ! We leave the data attribute empty for these, because we will coarsen it
    ! directly from the data attribute in the full resolution version of each
    ! diagnostic.
    coarse_diagnostic%axes = diagnostic%axes
    coarse_diagnostic%time_avg = diagnostic%time_avg
    coarse_diagnostic%mod_name = diagnostic%mod_name
    coarse_diagnostic%name = trim(diagnostic%name) // '_coarse'
    coarse_diagnostic%desc = diagnostic%desc
    coarse_diagnostic%unit = diagnostic%unit
    coarse_diagnostic%cnvfac = diagnostic%cnvfac
    coarse_diagnostic%coarse_graining_method = diagnostic%coarse_graining_method
  end subroutine populate_coarse_diag_type

  subroutine fv3gfs_diag_register_coarse(Time, coarse_axes)
    type(time_type), intent(in) :: Time
    integer, intent(in) :: coarse_axes(4)

    integer :: index

    do index = 1, DIAG_SIZE
       if (Diag(index)%name .eq. '') exit  ! No need to populate non-existent coarse diagnostics
       call populate_coarse_diag_type(Diag(index), Diag_coarse(index))
       Diag_coarse(index)%id = register_diag_field( &
            trim(Diag_coarse(index)%mod_name), trim(Diag_coarse(index)%name),  &
            coarse_axes(1:Diag_coarse(index)%axes), Time, trim(Diag_coarse(index)%desc), &
            trim(Diag_coarse(index)%unit), missing_value=real(missing_value))
    enddo
  end subroutine fv3gfs_diag_register_coarse

  subroutine register_coarse_diag_manager_controlled_diagnostics(Time, coarse_axes)
    type(time_type), intent(in) :: Time
    integer, intent(in) :: coarse_axes(4)

    integer :: index

    do index = 1, DIAG_SIZE
       if (Diag(index)%name .eq. '') exit  ! No need to populate non-existent coarse diagnostics
       call populate_coarse_diag_type(Diag_diag_manager_controlled(index), Diag_diag_manager_controlled_coarse(index))
       Diag_diag_manager_controlled_coarse(index)%id = register_diag_field( &
            trim(Diag_diag_manager_controlled_coarse(index)%mod_name), trim(Diag_diag_manager_controlled_coarse(index)%name),  &
            coarse_axes(1:Diag_diag_manager_controlled_coarse(index)%axes), Time, trim(Diag_diag_manager_controlled_coarse(index)%desc), &
            trim(Diag_diag_manager_controlled_coarse(index)%unit), missing_value=real(missing_value))
    enddo
  end subroutine register_coarse_diag_manager_controlled_diagnostics

  subroutine send_diag_manager_controlled_diagnostic_data(Time, Atm_block, IPD_Data, nx, ny, levs, &
    write_coarse_diagnostics, delp, coarsening_strategy, ptop)
    type(time_type),           intent(in) :: Time
    type(block_control_type),  intent(in) :: Atm_block
    type(IPD_data_type),       intent(in) :: IPD_Data(:)
    integer,                   intent(in) :: nx, ny, levs
    logical,                   intent(in) :: write_coarse_diagnostics
    real(kind=kind_phys),      intent(in) :: delp(isco:ieco,jsco:jeco,1:levo)
    character(len=64),         intent(in) :: coarsening_strategy
    real(kind=kind_phys),      intent(in) :: ptop

    logical :: require_area, require_masked_area, require_mass, require_masked_mass, require_vertical_remapping
    real(kind=kind_phys), allocatable :: area(:,:)
    real(kind=kind_phys), allocatable :: mass(:,:,:), phalf(:,:,:), phalf_coarse_on_fine(:,:,:)
    real(kind=kind_phys), allocatable :: masked_area(:,:,:)

    real(kind=kind_phys) :: scalar
    real(kind=kind_phys) :: var2d(nx, ny)
    real(kind=kind_phys) :: var3d(nx, ny, levs)
    integer :: i, j, ii, jj, k, isc, jsc, ix, nb, index, used

    isc   = atm_block%isc
    jsc   = atm_block%jsc

    if (write_coarse_diagnostics) then
      call determine_required_coarse_graining_weights(Diag_diag_manager_controlled_coarse, coarsening_strategy, &
                & require_area, require_masked_area, require_mass, require_vertical_remapping)
      if (.not. require_vertical_remapping) then
        if (require_area) then
          allocate(area(nx, ny))
          call get_area(Atm_block, IPD_Data, nx, ny, area)
        endif
        if (require_mass) then
          allocate(mass(nx, ny, levs))
          call get_mass(Atm_block, IPD_Data, delp, nx, ny, levs, mass)
        endif
      else
        allocate(area(nx, ny))
        allocate(phalf(nx, ny, levs + 1))
        allocate(phalf_coarse_on_fine(nx, ny, levs + 1))
        allocate(masked_area(nx, ny, levs))
        call get_area(Atm_block, IPD_Data, nx, ny, area)
        call vertical_remapping_requirements(delp, area, ptop, phalf, phalf_coarse_on_fine)
        call mask_area_weights(area, phalf, phalf_coarse_on_fine, masked_area)
      endif
    endif

    do index = 1, DIAG_SIZE
      if (trim(Diag_diag_manager_controlled(index)%name) .eq. '') exit
      if (Diag_diag_manager_controlled(index)%id .gt. 0 .or. Diag_diag_manager_controlled_coarse(index)%id .gt. 0) then
        if (Diag_diag_manager_controlled(index)%axes .eq. 2) then
          if (trim(Diag_diag_manager_controlled(index)%name) .eq. 'ocean_fraction' .or. &
              trim(Diag_diag_manager_controlled(index)%name) .eq. 'land_fraction' .or. &
              trim(Diag_diag_manager_controlled(index)%name) .eq. 'sea_ice_fraction') then
            call compute_surface_type_fraction(Atm_block, isc, jsc, nx, ny, Diag_diag_manager_controlled(index), var2d)
          else
            do j = 1, ny
              jj = j + jsc - 1
              do i = 1, nx
                  ii = i + isc - 1
                  nb = Atm_block%blkno(ii,jj)
                  ix = Atm_block%ixp(ii,jj)
                  var2d(i,j) = Diag_diag_manager_controlled(index)%data(nb)%var2(ix)
              enddo
            enddo
          endif
          if (Diag_diag_manager_controlled(index)%id > 0) then
            used = send_data(Diag_diag_manager_controlled(index)%id, var2d, Time)
          endif
          if (Diag_diag_manager_controlled_coarse(index)%id > 0) then
            call store_data2D_coarse(Diag_diag_manager_controlled_coarse(index)%id, Diag_diag_manager_controlled_coarse(index)%name, &
                & Diag_diag_manager_controlled_coarse(index)%coarse_graining_method, nx, ny, var2d, area, Time)
          endif
        elseif (Diag_diag_manager_controlled(index)%axes .eq. 3) then
          do k=1, levs
            do j = 1, ny
              jj = j + jsc - 1
              do i = 1, nx
                  ii = i + isc - 1
                  nb = Atm_block%blkno(ii,jj)
                  ix = Atm_block%ixp(ii,jj)
                  var3d(i,j,k) = Diag_diag_manager_controlled(index)%data(nb)%var3(ix,levs - k + 1)
              enddo
            enddo
          enddo
          if (Diag_diag_manager_controlled(index)%id .gt. 0) then
            used = send_data(Diag_diag_manager_controlled(index)%id, var3d, Time)
          endif
          if (Diag_diag_manager_controlled_coarse(index)%id > 0) then
            if (trim(coarsening_strategy) .eq. MODEL_LEVEL) then
              call store_data3D_coarse_model_level(Diag_diag_manager_controlled_coarse(index)%id, &
                  & Diag_diag_manager_controlled_coarse(index)%name, &
                  & Diag_diag_manager_controlled_coarse(index)%coarse_graining_method, &
                  & nx, ny, levs, var3d, area, mass, Time)
            elseif (trim(coarsening_strategy) .eq. PRESSURE_LEVEL) then
              call store_data3D_coarse_pressure_level(Diag_diag_manager_controlled_coarse(index)%id, &
                  & Diag_diag_manager_controlled_coarse(index)%name, &
                  & Diag_diag_manager_controlled_coarse(index)%coarse_graining_method, &
                  & nx, ny, levs, var3d, phalf, phalf_coarse_on_fine, masked_area, Time, ptop)
            else
              call mpp_error(FATAL, 'Invalid coarse-graining strategy provided.')
            endif
          endif
        elseif (trim(Diag_diag_manager_controlled(index)%name) .eq. 'global_mean_co2') then
          if (Diag_diag_manager_controlled(index)%id > 0) then
             call compute_global_mean_co2(Atm_block, IPD_Data, nx, ny, Diag_diag_manager_controlled(index), scalar)
             used = send_data(Diag_diag_manager_controlled(index)%id, scalar, Time)
          endif
          if (Diag_diag_manager_controlled_coarse(index)%id > 0) then
             call mpp_error(FATAL, 'global_mean_co2_coarse is not a valid diagnostic; use global_mean_co2 instead.')
          endif
        endif
      endif
    enddo
  end subroutine send_diag_manager_controlled_diagnostic_data

!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!--- gfs_diag_output ---
!-------------------------------------------------------------------------
!    routine to transfer the diagnostic data to the GFDL FMS diagnostic
!    manager for eventual output to the history files.
!
!    calls:  send_data
!-------------------------------------------------------------------------
  subroutine gfdl_diag_output(Time, Atm_block, IPD_Data, nx, ny, fprint, &
                             levs, ntcw, ntoz, dt, time_int, time_intfull, &
                             fhswr, fhlwr, &
                             prt_stats, write_coarse_diagnostics, delp, &
                             coarsening_strategy, ptop)
!--- subroutine interface variable definitions
    logical :: fprint
    type(time_type),           intent(in) :: Time
    type (block_control_type), intent(in) :: Atm_block
    type(IPD_data_type),       intent(in) :: IPD_Data(:)
    integer,                   intent(in) :: nx, ny, levs, ntcw, ntoz
    real(kind=kind_phys),      intent(in) :: dt
    real(kind=kind_phys),      intent(in) :: time_int
    real(kind=kind_phys),      intent(in) :: time_intfull
    real(kind=kind_phys),      intent(in) :: fhswr, fhlwr
    logical,                   intent(in) :: prt_stats
    logical,                   intent(in) :: write_coarse_diagnostics
    real(kind=kind_phys),      intent(in) :: delp(isco:ieco,jsco:jeco,1:levo)
    character(len=64),         intent(in) :: coarsening_strategy
    real(kind=kind_phys),      intent(in) :: ptop
!--- local variables
    integer :: i, j, k, idx, nblks, nb, ix, ii, jj, kflip
    integer :: is_in, js_in, isc, jsc
    character(len=2) :: xtra
    real(kind=kind_phys), dimension(nx*ny) :: var2p
    real(kind=kind_phys), dimension(nx*ny,levs) :: var3p
    real(kind=kind_phys), dimension(nx,ny) :: var2, area, lat, lon, one, landmask, seamask, icemask
    real(kind=kind_phys), dimension(nx,ny,levs) :: var3
    real(kind=kind_phys) :: rdt, rtime_int, rtime_intfull, lcnvfac
    real(kind=kind_phys) :: rtime_radsw, rtime_radlw
    logical :: used

    ! Local variables required for coarse-grianing
    logical :: require_area, require_masked_area, require_mass, require_masked_mass, require_vertical_remapping
    real(kind=kind_phys), allocatable :: mass(:,:,:), phalf(:,:,:), phalf_coarse_on_fine(:,:,:)
    real(kind=kind_phys), allocatable :: masked_area(:,:,:)

     nblks = Atm_block%nblks
     rdt = 1.0d0/dt
     rtime_int = 1.0d0/time_int
     rtime_intfull = 1.0d0/time_intfull
     rtime_radsw   = 1.0d0/fhswr
     rtime_radlw   = 1.0d0/fhlwr

     isc = Atm_block%isc
     jsc = Atm_block%jsc
     is_in = Atm_block%isc
     js_in = Atm_block%jsc

     !Metrics
     do j = 1, ny
        jj = j + jsc -1
        do i = 1, nx
           ii = i + isc -1
           nb = Atm_block%blkno(ii,jj)
           ix = Atm_block%ixp(ii,jj)
           area(i,j) = IPD_Data(nb)%Grid%area(ix)
           lat(i,j)  = IPD_Data(nb)%Grid%xlat(ix)
           lon(i,j)  = IPD_Data(nb)%Grid%xlon(ix)
           one(i,j)  = 1.
           landmask(i,j) = IPD_Data(nb)%Sfcprop%slmsk(ix)
           seamask(i,j)  = 1. - landmask(i,j)
           icemask(i,j)  = landmask(i,j) - 1.
        enddo
     enddo

     if (write_coarse_diagnostics) then
        call determine_required_coarse_graining_weights(diag_coarse, coarsening_strategy, require_area, &
                & require_masked_area, require_mass, require_vertical_remapping)
        if (.not. require_vertical_remapping) then
          if (require_mass) then
             allocate(mass(nx, ny, levs))
             call get_mass(Atm_block, IPD_Data, delp, nx, ny, levs, mass)
          endif
       else
          allocate(phalf(nx, ny, levs + 1))
          allocate(phalf_coarse_on_fine(nx, ny, levs + 1))
          allocate(masked_area(nx, ny, levs))
          call vertical_remapping_requirements(delp, area, ptop, phalf, phalf_coarse_on_fine)
          call mask_area_weights(area, phalf, phalf_coarse_on_fine, masked_area)
       endif
    endif

     do idx = 1,tot_diag_idx
       if ((Diag(idx)%id > 0) .or. (diag_coarse(idx)%id > 0)) then
         lcnvfac = Diag(idx)%cnvfac
         if (Diag(idx)%time_avg) then
           if ( trim(Diag(idx)%time_avg_kind) == 'full' ) then
             lcnvfac = lcnvfac*rtime_intfull
           else if ( trim(Diag(idx)%time_avg_kind) == 'rad_lw' ) then
             lcnvfac = lcnvfac*min(rtime_radlw,rtime_int)
           else if ( trim(Diag(idx)%time_avg_kind) == 'rad_sw' ) then
             lcnvfac = lcnvfac*min(rtime_radsw,rtime_int)
           else if ( trim(Diag(idx)%time_avg_kind) == 'rad_swlw_min' ) then
             lcnvfac = lcnvfac*min(max(rtime_radsw,rtime_radlw),rtime_int)
           else
             lcnvfac = lcnvfac*rtime_int
           endif
         endif
         if (Diag(idx)%axes == 2) then
           if (trim(diag(idx)%mask) == 'positive_flux') then
             !--- albedos are actually a ratio of two radiation surface properties
             var2(1:nx,1:ny) = 0._kind_phys
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                 if (Diag(idx)%data(nb)%var21(ix) > 0._kind_phys) &
                   var2(i,j) = max(0._kind_phys,min(1._kind_phys,Diag(idx)%data(nb)%var2(ix)/Diag(idx)%data(nb)%var21(ix)))*lcnvfac
               enddo
             enddo
           elseif (trim(Diag(idx)%mask) == 'land_ice_only') then
             !--- need to "mask" gflux to output valid data over land/ice only
             var2(1:nx,1:ny) = missing_value
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                  if (Diag(idx)%data(nb)%var21(ix) /= 0) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
               enddo
             enddo
           elseif (trim(Diag(idx)%mask) == 'land_only') then
             !--- need to "mask" soilm to have value only over land
             var2(1:nx,1:ny) = missing_value
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                 if (Diag(idx)%data(nb)%var21(ix) == 1) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
               enddo
             enddo
           elseif (trim(Diag(idx)%mask) == 'cldmask') then
             !--- need to "mask" soilm to have value only over land
             var2(1:nx,1:ny) = missing_value
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                 if (Diag(idx)%data(nb)%var21(ix)*100. > 0.5) var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
               enddo
             enddo
           elseif (trim(Diag(idx)%mask) == 'cldmask_ratio') then
             !--- need to "mask" soilm to have value only over land
             var2(1:nx,1:ny) = missing_value
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                 if (Diag(idx)%data(nb)%var21(ix)*100.*lcnvfac > 0.5) var2(i,j) = Diag(idx)%data(nb)%var2(ix)/ &
                     Diag(idx)%data(nb)%var21(ix)
               enddo
             enddo
           elseif (trim(Diag(idx)%mask) == '') then
             do j = 1, ny
               jj = j + jsc -1
               do i = 1, nx
                 ii = i + isc -1
                 nb = Atm_block%blkno(ii,jj)
                 ix = Atm_block%ixp(ii,jj)
                 var2(i,j) = Diag(idx)%data(nb)%var2(ix)*lcnvfac
               enddo
             enddo
           endif
           if (Diag(idx)%id > 0) then
              used=send_data(Diag(idx)%id, var2, Time)
            endif
            if (Diag_coarse(idx)%id > 0) then
               call store_data2D_coarse(Diag_coarse(idx)%id, Diag_coarse(idx)%name, &
                    Diag_coarse(idx)%coarse_graining_method, nx, ny, var2, area, Time)
            endif

           !!!! Accumulated diagnostics --- lmh 19 sep 17
           if (fprint .and. prt_stats) then
           select case (trim(Diag(idx)%name))
           case('totprcpb_ave')
              call prt_gb_nh_sh_us('Total Precip (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, one, 86400.)
              call prt_gb_nh_sh_us('Land Precip  (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, landmask, 86400.)
              call prt_gb_nh_sh_us('Ocean Precip (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, seamask, 86400.)
              call prt_gb_nh_sh_us('SeaIce Precip (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, icemask, 86400.)
           case('cnvprcpb_ave')
              call prt_gb_nh_sh_us('Total Convective Precip (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, one, 86400.)
              call prt_gb_nh_sh_us('Land Convective Precip  (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, landmask, 86400.)
              call prt_gb_nh_sh_us('Ocean Convective Precip (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, seamask, 86400.)
              call prt_gb_nh_sh_us('SeaIce Convective Precip (mm/d)', 1, nx, 1, ny, var2, area, lon, lat, icemask, 86400.)
           case('totsnwb_ave')
              call prt_gb_nh_sh_us('Total Snowfall (9:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, one, 777600.)
              call prt_gb_nh_sh_us('Land Snowfall  (9:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, landmask, 777600.)
              call prt_gb_nh_sh_us('Ocean Snowfall (9:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, seamask, 777600.)
              call prt_gb_nh_sh_us('SeaIce Snowfall (9:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, icemask, 777600.)
!           case('totgrp') ! Tiny??
!              call prt_gb_nh_sh_us('Total Icefall (2:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, one, 172800.)
!              call prt_gb_nh_sh_us('Land Icefall  (2:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, landmask, 172800.)
!              call prt_gb_nh_sh_us('Ocean Icefall (2:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, seamask, 172800.)
!              call prt_gb_nh_sh_us('SeaIce Icefall (2:1 mm/d)', 1, nx, 1, ny, var2, area, lon, lat, icemask, 172800.)
           case('lhtfl_ave')
              call prt_gb_nh_sh_us('Total sfc LH flux ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land sfc LH flux  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean sfc LH flux ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce sfc LH flux ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('shtfl_ave')
              call prt_gb_nh_sh_us('Total sfc SH flux ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land sfc SH flux  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean sfc SH flux ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce sfc SH flux ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('hpbl')
              call prt_gb_nh_sh_us('Total pbl height ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land pbl height  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean pbl height ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce pbl height ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('dusfc')
              call prt_gb_nh_sh_us('Total u-wind stress ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land u-wind stress  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean u-wind stress ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce u-wind stress ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('dvsfc')
              call prt_gb_nh_sh_us('Total v-wind stress ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land v-wind stress  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean v-wind stress ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce v-wind stress ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('DSWRFtoa')
              call prt_gb_nh_sh_us('TOA SW down ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
           case('USWRFtoa')
              call prt_gb_nh_sh_us('TOA SW up ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
           case('ULWRFtoa')
              call prt_gb_nh_sh_us('TOA LW up ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
           case('u10m')
              call prt_gb_nh_sh_us('Total 10-m u avg ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land 10-m u avg  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean 10-m u avg ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce 10-m u avg ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('v10m')
              call prt_gb_nh_sh_us('Total 10-m v avg ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land 10-m v avg  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean 10-m v avg ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce 10-m v avg ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('acond')
              call prt_gb_nh_sh_us('Total momentum exchange coefficient ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land momentum exchange coefficient  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean momentum exchange coefficient ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce momentum exchange coefficient ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('sfexc')
              call prt_gb_nh_sh_us('Total thermal exchange coefficient ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land thermal exchange coefficient  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean thermal exchange coefficient ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce thermal exchange coefficient ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('ffmm')
              call prt_gb_nh_sh_us('Total ffmm for PBL ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land ffmm for PBL  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean ffmm for PBL ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce ffmm for PBL ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('ffhh')
              call prt_gb_nh_sh_us('Total ffhh for PBL ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land ffhh for PBL  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean ffhh for PBL ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce ffhh for PBL ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('ZORLsfc')
              call prt_gb_nh_sh_us('Total surface roughness ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land surface roughness  ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean surface roughness ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce surface roughness ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('q2m')
              call prt_gb_nh_sh_us('Total 2-m Q avg ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land 2-m Q avg ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean 2-m Q avg ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce 2-m Q avg ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
           case('t2m')
              call prt_gb_nh_sh_us('Total 2-m T avg ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land 2-m T avg ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean 2-m T avg ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce 2-m T avg ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
              call prt_gb_nh_sh_us('2-m T max ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MAX')
              call prt_gb_nh_sh_us('2-m T min ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MIN')
           case('tsfc')
              call prt_gb_nh_sh_us('Total sfc T avg ', 1, nx, 1, ny, var2, area, lon, lat, one, 1.)
              call prt_gb_nh_sh_us('Land sfc T avg ', 1, nx, 1, ny, var2, area, lon, lat, landmask, 1.)
              call prt_gb_nh_sh_us('Ocean sfc T avg ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1.)
              call prt_gb_nh_sh_us('SeaIce sfc T avg ', 1, nx, 1, ny, var2, area, lon, lat, icemask, 1.)
              call prt_gb_nh_sh_us('sfc T max ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MAX')
              call prt_gb_nh_sh_us('sfc T min ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MIN')
              call prt_gb_nh_sh_us('SST max ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1., 'MAX')
              call prt_gb_nh_sh_us('SST min ', 1, nx, 1, ny, var2, area, lon, lat, seamask, 1., 'MIN')
           case('ps_dt')
              call prt_gb_nh_sh_us('ps_dt max ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MAX')
              call prt_gb_nh_sh_us('ps_dt min ', 1, nx, 1, ny, var2, area, lon, lat, one, 1., 'MIN')
           end select
           endif
         elseif (Diag(idx)%axes == 3) then
           !--- dt3dt variables ---- restored 16 feb 18 lmh
            do k=1,levs
            kflip=levs+1-k
            do j=1,ny
            jj = j + jsc -1
            do i=1,nx
               ii = i + isc - 1
               nb = Atm_block%blkno(ii,jj)
               ix = Atm_block%ixp(ii,jj)
               var3(i,j,k) = Diag(idx)%data(nb)%var3(ix,kflip)*lcnvfac
            enddo
            enddo
            enddo
            if (Diag(idx)%id > 0) then
               used=send_data(Diag(idx)%id, var3, Time)
            endif
            if (Diag_coarse(idx)%id > 0) then
               if (trim(coarsening_strategy) .eq. MODEL_LEVEL) then
                  call store_data3D_coarse_model_level(Diag_coarse(idx)%id, Diag_coarse(idx)%name, &
                       Diag_coarse(idx)%coarse_graining_method, &
                       nx, ny, levo, var3, area, mass, Time)
               elseif (trim(coarsening_strategy) .eq. PRESSURE_LEVEL) then
                  call store_data3D_coarse_pressure_level(Diag_coarse(idx)%id, Diag_coarse(idx)%name, &
                       Diag_coarse(idx)%coarse_graining_method, &
                       nx, ny, levo, var3, phalf, phalf_coarse_on_fine, masked_area, Time, ptop)
               else
                  call mpp_error(FATAL, 'Invalid coarse-graining strategy provided.')
               endif
            endif

#ifdef JUNK
           !--- dq3dt variables
           do num = 1,5+Mdl_parms%pl_coeff
             write(xtra,'(i1)') num
             if (trim(Diag(idx)%name) == 'dq3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dq3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
             endif
           enddo
           !--- du3dt and dv3dt variables
           do num = 1,4
             write(xtra,'(i1)') num
             if (trim(Diag(idx)%name) == 'du3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%du3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
             endif
             if (trim(Diag(idx)%name) == 'dv3dt_'//trim(xtra)) then
               var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dv3dt(1:ngptc,levs:1:-1,num:num), (/nx,ny,levs/))
               used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
             endif
           enddo
           if (trim(Diag(idx)%name) == 'dqdt_v') then
             var3(1:nx,1:ny,1:levs) = RESHAPE(Gfs_diag%dqdt_v(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- temperature tendency
           if (trim(Diag(idx)%name) == 'dtemp_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%tgrs(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gt0(1:ngptc,levs:1:-1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- horizontal wind component tendency
           if (trim(Diag(idx)%name) == 'du_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%ugrs(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gu0(1:ngptc,levs:1:-1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- meridional wind component tendency
           if (trim(Diag(idx)%name) == 'dv_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%vgrs(1:ngptc,levs:1:-1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gv0(1:ngptc,levs:1:-1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- specific humidity tendency
           if (trim(Diag(idx)%name) == 'dsphum_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%qgrs(1:ngptc,levs:1:-1,1:1), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gq0(1:ngptc,levs:1:-1,1:1), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- cloud water mixing ration tendency
           if (trim(Diag(idx)%name) == 'dclwmr_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%qgrs(1:ngptc,levs:1:-1,ntcw:ntcw), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gq0(1:ngptc,levs:1:-1,ntcw:ntcw), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
           endif
           !--- ozone mixing ration tendency
           if (trim(Diag(idx)%name) == 'do3mr_dt') then
             var3(1:nx,1:ny,1:levs) =  RESHAPE(Statein%qgrs(1:ngptc,levs:1:-1,ntoz:ntoz), (/nx,ny,levs/))
             var3(1:nx,1:ny,1:levs) = (RESHAPE(Stateout%gq0(1:ngptc,levs:1:-1,ntoz:ntoz), (/nx,ny,levs/))  &
                                        - var3(1:nx,1:ny,1:levs))*rdt
             used=send_data(Diag(idx)%id, var3, Time, is_in=is_in, js_in=js_in, ks_in=1)
          endif
#endif
       endif
    endif
 enddo


  end subroutine gfdl_diag_output

 subroutine compute_global_mean_co2(Atm_block, IPD_Data, nx, ny, Diag, global_mean_co2)
    type (block_control_type), intent(in) :: Atm_block
    type(IPD_data_type),       intent(in) :: IPD_Data(:)
    integer, intent(in) :: nx, ny
    type(gfdl_diag_type), intent(in) :: Diag
    real(kind=kind_phys), intent(out) :: global_mean_co2

    real(kind=kind_phys) :: moles_dry_air, moles_co2, area
    integer :: j, jj, i, ii, nb, ix, isc, jsc

    moles_dry_air = 0.0
    moles_co2 = 0.0

    isc = Atm_block%isc
    jsc = Atm_block%jsc

    do j = 1, ny
       jj = j + jsc - 1
       do i = 1, nx
          ii = i + isc - 1
          nb = Atm_block%blkno(ii,jj)
          ix = Atm_block%ixp(ii,jj)
          area = IPD_Data(nb)%Grid%area(ix)
          moles_dry_air = moles_dry_air + area * Diag%data(nb)%var21(ix)
          moles_co2 = moles_co2 + area * Diag%data(nb)%var2(ix)
        enddo
     enddo

    call mp_reduce_sum(moles_dry_air)
    call mp_reduce_sum(moles_co2)

    global_mean_co2 = moles_co2 / moles_dry_air
 end subroutine compute_global_mean_co2

!-------------------------------------------------------------------------
 subroutine prt_gb_nh_sh_us(qname, is,ie, js,je, a2, area, lon, lat, mask, fac, operation_in) !Prints averages/sums, or maxes/mins
  use physcons,    pi=>con_pi
  character(len=*), intent(in)::  qname
  integer, intent(in):: is, ie, js, je
  real(kind=kind_phys), intent(in), dimension(is:ie, js:je):: a2
  real(kind=kind_phys), intent(in), dimension(is:ie, js:je):: area, lon, lat, mask
  real, intent(in) :: fac
  character(len=*), intent(in), OPTIONAL :: operation_in
! Local:
  real(kind=kind_phys), parameter:: rad2deg = 180./pi
  real(kind=kind_phys) :: slat, slon
  real(kind=kind_phys):: t_eq, t_nh, t_sh, t_gb, t_us
  real(kind=kind_phys):: area_eq, area_nh, area_sh, area_gb, area_us
  integer:: i,j
  character(len=100) :: diagstr
  character(len=20)  :: diagstr1
  character(len=3)  :: operation

  if (present(operation_in)) then
     operation = operation_in(1:3)
  else
     operation = 'SUM'
  endif

     if (operation == "MAX") then
        t_eq =-1.e14   ;    t_nh =-1.e14;    t_sh =-1.e14;    t_gb =-1.e14;    t_us =-1.e14
        area_eq = 0.   ; area_nh = 0.   ; area_sh = 0.   ; area_gb = 0.   ; area_us = 0.
        do j=js,je
        do i=is,ie
           if (mask(i,j) <= 1.e-6) cycle

           slat = lat(i,j)*rad2deg
           slon = lon(i,j)*rad2deg
           area_gb = 1.
           t_gb = max(t_gb,a2(i,j))
           if( (slat>-20. .and. slat<20.) ) then
                area_eq = 1.
                t_eq = max(t_eq,a2(i,j))
           elseif( slat>=20. .and. slat<80. ) then
                area_nh = 1.
                t_nh = max(t_nh,a2(i,j))
           elseif( slat<=-20. .and. slat>-80. ) then
                area_sh = 1.
                t_sh = max(t_sh,a2(i,j))
           endif
           if ( slat>25.  .and. slat<50. .and. &
                slon>235. .and. slon<300. ) then
              area_us = 1.
              t_us = max(t_us,a2(i,j))
           endif
        enddo
        enddo

        call mp_reduce_max(   t_gb)
        call mp_reduce_max(   t_nh)
        call mp_reduce_max(   t_sh)
        call mp_reduce_max(   t_eq)
        call mp_reduce_max(   t_us)
     elseif (operation == "MIN") then
        t_eq = 1.e14   ;    t_nh = 1.e14;    t_sh = 1.e14;    t_gb = 1.e14;    t_us = 1.e14
        area_eq = 0.   ; area_nh = 0.   ; area_sh = 0.   ; area_gb = 0.   ; area_us = 0.
        do j=js,je
        do i=is,ie
           if (mask(i,j) <= 1.e-6) cycle

           slat = lat(i,j)*rad2deg
           slon = lon(i,j)*rad2deg
           area_gb = 1.
           t_gb = min(t_gb,a2(i,j))
           if( (slat>-20. .and. slat<20.) ) then
                area_eq = 1.
                t_eq = min(t_eq,a2(i,j))
           elseif( slat>=20. .and. slat<80. ) then
                area_nh = 1.
                t_nh = min(t_nh,a2(i,j))
           elseif( slat<=-20. .and. slat>-80. ) then
                area_sh = 1.
                t_sh = min(t_sh,a2(i,j))
           endif
           if ( slat>25.  .and. slat<50. .and. &
                slon>235. .and. slon<300. ) then
              area_us = 1.
              t_us = min(t_us,a2(i,j))
           endif
        enddo
        enddo

        call mp_reduce_min(   t_gb)
        call mp_reduce_min(   t_nh)
        call mp_reduce_min(   t_sh)
        call mp_reduce_min(   t_eq)
        call mp_reduce_min(   t_us)
     else
        t_eq = 0.   ;    t_nh = 0.;    t_sh = 0.;    t_gb = 0.;    t_us = 0.
        area_eq = 0.; area_nh = 0.; area_sh = 0.; area_gb = 0.; area_us = 0.
        operation = 'SUM'
        do j=js,je
        do i=is,ie
           slat = lat(i,j)*rad2deg
           slon = lon(i,j)*rad2deg
           area_gb = area_gb + area(i,j)*mask(i,j)
           t_gb = t_gb + a2(i,j)*area(i,j)*mask(i,j)
           if( (slat>-20. .and. slat<20.) ) then
                area_eq = area_eq + area(i,j)*mask(i,j)
                t_eq = t_eq + a2(i,j)*area(i,j)*mask(i,j)
           elseif( slat>=20. .and. slat<80. ) then
                area_nh = area_nh + area(i,j)*mask(i,j)
                t_nh = t_nh + a2(i,j)*area(i,j)*mask(i,j)
           elseif( slat<=-20. .and. slat>-80. ) then
                area_sh = area_sh + area(i,j)*mask(i,j)
                t_sh = t_sh + a2(i,j)*area(i,j)*mask(i,j)
           endif
           if ( slat>25.  .and. slat<50. .and. &
                slon>235. .and. slon<300. ) then
              area_us = area_us + area(i,j)*mask(i,j)
              t_us = t_us + a2(i,j)*area(i,j)*mask(i,j)
           endif
        enddo
        enddo

        call mp_reduce_sum(area_gb)
        call mp_reduce_sum(   t_gb)
        call mp_reduce_sum(area_nh)
        call mp_reduce_sum(   t_nh)
        call mp_reduce_sum(area_sh)
        call mp_reduce_sum(   t_sh)
        call mp_reduce_sum(area_eq)
        call mp_reduce_sum(   t_eq)
        call mp_reduce_sum(area_us)
        call mp_reduce_sum(   t_us)
     endif

     diagstr = trim(qname) // ' ' // trim(mpp_get_current_pelist_name()) // ' '
     !if (area_gb < 1.) then
     !   diagstr1 = ''
     !elseif( area_gb <= 4.*pi*RADIUS*RADIUS*.98) then
     !   write(diagstr1,101) 'Grid', t_gb/area_gb*fac
     !else
     if (area_gb <= 1.e-6) return
     write(diagstr1,101) 'GB', t_gb/area_gb*fac
     !endif
     diagstr = trim(diagstr) // trim(diagstr1)
     if (area_nh <= 1.e-6 ) then
        diagstr1 = ''
     else
        write(diagstr1,101) 'NH', t_nh/area_nh*fac
     endif
     diagstr = trim(diagstr) // trim(diagstr1)
     if (area_sh <= 1.e-6 ) then
        diagstr1 = ''
     else
        write(diagstr1,101) 'SH', t_sh/area_sh*fac
     endif
     diagstr = trim(diagstr) // trim(diagstr1)
     if (area_eq <= 1.e-6) then
        diagstr1 = ''
     else
        write(diagstr1,101) 'EQ', t_eq/area_eq*fac
     endif
     diagstr = trim(diagstr) // trim(diagstr1)
     if (area_us <= 1.e-6) then
        diagstr1 = ''
     else
        write(diagstr1,101) 'US', t_us/area_us*fac
     endif
     diagstr = trim(diagstr) // trim(diagstr1)

     if (is_master()) write(*,'(A)') trim(diagstr)

101  format(3x, A, ': ', F7.2)

   end subroutine prt_gb_nh_sh_us

 subroutine determine_required_coarse_graining_weights(coarse_diag, coarsening_strategy, require_area, &
          & require_masked_area, require_mass, require_vertical_remapping)
   type(gfdl_diag_type), intent(in) :: coarse_diag(:)
   character(len=64), intent(in) :: coarsening_strategy
   logical, intent(out) :: require_area, require_masked_area, require_mass, require_vertical_remapping

   require_area = any(coarse_diag%id .gt. 0 .and. coarse_diag%coarse_graining_method .eq. AREA_WEIGHTED)
   require_mass = any(coarse_diag%id .gt. 0 .and. coarse_diag%coarse_graining_method .eq. MASS_WEIGHTED)

   if (trim(coarsening_strategy) .eq. PRESSURE_LEVEL) then
     require_masked_area = any(coarse_diag%id .gt. 0 .and. coarse_diag%axes .eq. 3 .and. &
                             & coarse_diag%coarse_graining_method .eq. AREA_WEIGHTED)
     require_vertical_remapping = any(coarse_diag%id .gt. 0 .and. coarse_diag%axes .eq. 3)
   else
     require_masked_area = .false.
     require_vertical_remapping = .false.
   endif
 end subroutine determine_required_coarse_graining_weights

 subroutine get_area(Atm_block, IPD_Data, nx, ny, area)
   type(block_control_type), intent(in) :: Atm_block
   type(IPD_data_type), intent(in) :: IPD_Data(:)
   integer, intent(in) :: nx, ny
   real(kind=kind_phys), intent(out) :: area(1:nx,1:ny)

   integer :: i, ii, j, jj, block_number, column
   do j = 1, ny
      jj = j + jsco - 1
      do i = 1, nx
         ii = i + isco - 1
         block_number = Atm_block%blkno(ii,jj)
         column = Atm_block%ixp(ii,jj)
         area(i,j) = IPD_Data(block_number)%Grid%area(column)
      enddo
   enddo
 end subroutine get_area

 subroutine get_mass(Atm_block, IPD_Data, delp, nx, ny, nz, mass)
   type(block_control_type), intent(in) :: Atm_block
   type(IPD_data_type), intent(in) :: IPD_Data(:)
   integer, intent(in) :: nx, ny, nz
   real(kind=kind_phys), intent(in) :: delp(1:nx,1:ny,1:nz)
   real(kind=kind_phys), intent(out) :: mass(1:nx,1:ny,1:nz)

   integer :: i, ii, j, jj, k, block_number, column, isc, jsc
   real(kind=kind_phys) :: area_value

   do k = 1, nz
      do j = 1, ny
         jj = j + jsco - 1
         do i = 1, nx
            ii = i + isco - 1
            block_number = Atm_block%blkno(ii,jj)
            column = Atm_block%ixp(ii,jj)
            area_value = IPD_Data(block_number)%Grid%area(column)
            mass(i,j,k) = area_value * delp(i,j,k)
         enddo
      enddo
   enddo
 end subroutine get_mass

 subroutine store_data2D_coarse(id, name, method, nx, ny, full_resolution_field, area, Time)
   integer, intent(in) :: id
   character(len=64), intent(in) :: name
   character(len=64), intent(in) :: method
   integer, intent(in) :: nx, ny
   real(kind=kind_phys), intent(in) :: full_resolution_field(1:nx,1:ny)
   real(kind=kind_phys), intent(in) :: area(1:nx,1:ny)
   type(time_type), intent(in) :: Time

   real(kind=kind_phys), allocatable :: coarse(:,:)
   character(len=128) :: message
   integer :: is_coarse, ie_coarse, js_coarse, je_coarse, nx_coarse, ny_coarse
   logical :: used

   call get_coarse_array_bounds(is_coarse, ie_coarse, js_coarse, je_coarse)
   nx_coarse = ie_coarse - is_coarse + 1
   ny_coarse = je_coarse - js_coarse + 1

   allocate(coarse(nx_coarse, ny_coarse))

   if (method .eq. AREA_WEIGHTED) then
      call weighted_block_average(area, full_resolution_field, coarse)
   elseif (method .eq. MASKED_AREA_WEIGHTED) then
      call weighted_block_average(area, full_resolution_field, full_resolution_field .ne. missing_value, coarse)
   elseif (method .eq. MODE) then
      call block_mode(full_resolution_field, coarse)
   elseif (method .eq. MASS_WEIGHTED) then
      message = 'mass_weighted is not a valid coarse_graining_method for 2D variable ' // trim(name)
      call mpp_error(FATAL, message)
   else
      message = 'A valid coarse_graining_method must be specified for ' // trim(name)
      call mpp_error(FATAL, message)
   endif
   used = send_data(id, coarse, Time)
 end subroutine store_data2D_coarse

 subroutine store_data3D_coarse_model_level(id, name, method, nx, ny, nz, full_resolution_field, &
      area, mass, Time)
   integer, intent(in) :: id
   character(len=64), intent(in) :: name
   character(len=64), intent(in) :: method
   integer, intent(in) :: nx, ny, nz
   real(kind=kind_phys), intent(in) :: full_resolution_field(1:nx,1:ny,1:nz)
   real(kind=kind_phys), intent(in) :: area(1:nx,1:ny)
   real(kind=kind_phys), intent(in) :: mass(1:nx,1:ny,1:nz)
   type(time_type), intent(in) :: Time

   real(kind=kind_phys), allocatable :: coarse(:,:,:)
   character(len=128) :: message
   integer :: is_coarse, ie_coarse, js_coarse, je_coarse, nx_coarse, ny_coarse
   logical :: used

   call get_coarse_array_bounds(is_coarse, ie_coarse, js_coarse, je_coarse)
   nx_coarse = ie_coarse - is_coarse + 1
   ny_coarse = je_coarse - js_coarse + 1

   allocate(coarse(nx_coarse, ny_coarse, nz))

   if (method .eq. AREA_WEIGHTED) then
      call weighted_block_average(area, full_resolution_field, coarse)
   elseif (method .eq. MASKED_AREA_WEIGHTED) then
      message = 'Masked area-weighted coarse-graining is not currently implemented for 3D variables'
      call mpp_error(FATAL, message)
   elseif (method .eq. MASS_WEIGHTED) then
      call weighted_block_average(mass, full_resolution_field, coarse)
   elseif (method .eq. MODE) then
      message = 'Block mode coarse-graining is not currently implemented for 3D variables'
      call mpp_error(FATAL, message)
   else
      message = 'A valid coarse_graining_method must be specified for ' // trim(name)
      call mpp_error(FATAL, message)
   endif
   used = send_data(id, coarse, Time)
 end subroutine store_data3D_coarse_model_level

 subroutine store_data3D_coarse_pressure_level(id, name, method, nx, ny, nz, full_resolution_field, &
        & phalf, phalf_coarse_on_fine, masked_area, Time, ptop)
    integer, intent(in) :: id
    character(len=64), intent(in) :: name
    character(len=64), intent(in) :: method
    integer, intent(in) :: nx, ny, nz
    real(kind=kind_phys), intent(in) :: full_resolution_field(1:nx,1:ny,1:nz)
    real(kind=kind_phys), intent(in) :: phalf(1:nx,1:ny,1:nz + 1)
    real(kind=kind_phys), intent(in) :: phalf_coarse_on_fine(1:nx,1:ny,1:nz + 1)
    real(kind=kind_phys), intent(in) :: masked_area(1:nx,1:ny,1:nz)
    type(time_type), intent(in) :: Time
    real(kind=kind_phys), intent(in) :: ptop

    real(kind=kind_phys), allocatable :: remapped(:,:,:), coarse(:,:,:)
    character(len=128) :: message
    integer :: is_coarse, ie_coarse, js_coarse, je_coarse, nx_coarse, ny_coarse
    logical :: used

    call get_coarse_array_bounds(is_coarse, ie_coarse, js_coarse, je_coarse)
    nx_coarse = ie_coarse - is_coarse + 1
    ny_coarse = je_coarse - js_coarse + 1

    allocate(remapped(nx, ny, nz))
    allocate(coarse(nx_coarse, ny_coarse, nz))

    call vertically_remap_field(phalf, full_resolution_field, phalf_coarse_on_fine, ptop, remapped)

    ! AREA_WEIGHTED and MASS_WEIGHTED are equivalent in pressure level coarse-graining
    if (method .eq. AREA_WEIGHTED .or. method .eq. MASS_WEIGHTED) then
      call weighted_block_average(masked_area, remapped, coarse)
    elseif (method .eq. MASKED_AREA_WEIGHTED) then
      message = 'Masked area-weighted coarse-graining is not currently implemented for 3D variables'
      call mpp_error(FATAL, message)
    elseif (method .eq. MODE) then
      message = 'Block mode coarse-graining is not currently implemented for 3D variables'
      call mpp_error(FATAL, message)
    else
      message = 'A valid coarse_graining_method must be specified for ' // trim(name)
      call mpp_error(FATAL, message)
    endif
    used = send_data(id, coarse, Time)
end subroutine store_data3D_coarse_pressure_level

end module FV3GFS_io_mod



