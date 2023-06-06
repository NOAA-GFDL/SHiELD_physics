module cosp2_test

  use physcons,            ONLY: grav => con_g
  use GFS_typedefs,        ONLY: GFS_control_type, GFS_diag_type,                         &
                                 GFS_statein_type, GFS_stateout_type, GFS_sfcprop_type,   &
                                 GFS_radtend_type, GFS_init_type
  
  implicit none

  public :: cosp2_offline

contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  ! SUBROUTINE cosp2_offline
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine cosp2_offline (Model, Statein, Stateout, Sfcprop, Radtend, Diag, Init_parm)

  implicit none

  type(GFS_init_type), intent(in) :: Init_parm
  type (GFS_control_type), intent (in) :: Model
  type (GFS_statein_type), intent (in) :: Statein(:)
  type (GFS_stateout_type), intent (in) :: Stateout(:)
  type (GFS_sfcprop_type), intent (in) :: Sfcprop(:)
  type (GFS_radtend_type), intent (in) :: Radtend(:)

  type (GFS_diag_type), intent (inout) :: Diag(:)

  integer :: nb, Nlevels, nblks

  nblks = size(Init_parm%blksz)

  Nlevels = Model%levs

  do nb = 1, nblks

    Diag(nb)%cosp%p           = Statein(nb)%prsl
    Diag(nb)%cosp%ph          = Statein(nb)%prsi(:,1:Nlevels)
    Diag(nb)%cosp%zlev        = Statein(nb)%phil / grav
    Diag(nb)%cosp%zlev_half   = Statein(nb)%phii(:,1:Nlevels) / grav
    Diag(nb)%cosp%T           = Stateout(nb)%gt0
    Diag(nb)%cosp%sh          = Stateout(nb)%gq0(:,:,1)
    Diag(nb)%cosp%tca         = Stateout(nb)%gq0(:,:,Model%ntclamt)
    Diag(nb)%cosp%cca         = 0
    Diag(nb)%cosp%mr_lsliq    = Stateout(nb)%gq0(:,:,Model%ntcw)
    Diag(nb)%cosp%mr_lsice    = Stateout(nb)%gq0(:,:,Model%ntiw)
    Diag(nb)%cosp%mr_ccliq    = 0.0
    Diag(nb)%cosp%mr_ccice    = 0.0
    Diag(nb)%cosp%fl_lsrain   = Diag(nb)%pfr / 86400.
    Diag(nb)%cosp%fl_lssnow   = Diag(nb)%pfs / 86400.
    Diag(nb)%cosp%fl_lsgrpl   = Diag(nb)%pfg / 86400.
    Diag(nb)%cosp%fl_ccrain   = 0.0
    Diag(nb)%cosp%fl_ccsnow   = 0.0
    Diag(nb)%cosp%Reff_LSCLIQ = Diag(nb)%reff(:,:,1) * 1.e-6
    Diag(nb)%cosp%Reff_LSCICE = Diag(nb)%reff(:,:,2) * 1.e-6
    Diag(nb)%cosp%Reff_LSRAIN = Diag(nb)%reff(:,:,3) * 1.e-6
    Diag(nb)%cosp%Reff_LSSNOW = Diag(nb)%reff(:,:,4) * 1.e-6
    Diag(nb)%cosp%Reff_LSGRPL = Diag(nb)%reff(:,:,5) * 1.e-6
    Diag(nb)%cosp%dtau_s      = Diag(nb)%ctau(:,:,1)
    Diag(nb)%cosp%dtau_c      = 0.0
    Diag(nb)%cosp%dem_s       = Diag(nb)%ctau(:,:,2)
    Diag(nb)%cosp%dem_c       = 0.0
    Diag(nb)%cosp%skt         = Sfcprop(nb)%tsfc
    Diag(nb)%cosp%landmask    = 1-abs(Sfcprop(nb)%slmsk-1)
    Diag(nb)%cosp%mr_ozone    = Stateout(nb)%gq0(:,:,Model%ntoz)
    Diag(nb)%cosp%u_wind      = Stateout(nb)%gu0
    Diag(nb)%cosp%v_wind      = Stateout(nb)%gv0
    Diag(nb)%cosp%sunlit      = ceiling(Radtend(nb)%coszen)
    Diag(nb)%cosp%surfelev    = Sfcprop(nb)%oro

  enddo
 
  end subroutine cosp2_offline

 end module cosp2_test

