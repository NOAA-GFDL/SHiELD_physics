module physics_abstraction_layer

  use GFS_typedefs,    only: init_type        =>  GFS_init_type,     &
                             control_type     =>  GFS_control_type,  &
                             statein_type     =>  GFS_statein_type,  &
                             stateout_type    =>  GFS_stateout_type, &
                             sfcprop_type     =>  GFS_sfcprop_type,  &
                             coupling_type    =>  GFS_coupling_type, &
                             grid_type        =>  GFS_grid_type,     &
                             statemid_type    =>  GFS_statemid_type, &
                             cldprop_type     =>  GFS_cldprop_type,  &
                             radtend_type     =>  GFS_radtend_type,  &
                             intdiag_type     =>  GFS_diag_type,     &
                             overrides_type   =>  GFS_overrides_type

  use GFS_driver,      only: initialize       =>  GFS_initialize,       &
                             time_vary_step   =>  GFS_time_vary_step,   &
                             radiation_step1  =>  GFS_radiation_driver, &
                             physics_step1_down    =>  GFS_physics_driver_down,   &
                             physics_step1_up      =>  GFS_physics_driver_up,   &
                             physics_step2    =>  GFS_stochastic_driver,&
                             physics_end      =>  GFS_physics_end

!----------------------
!  public physics types
!----------------------
  public  init_type
  public  control_type
  public  statein_type
  public  stateout_type
  public  sfcprop_type
  public  coupling_type
  public  grid_type
  public  statemid_type
  public  cldprop_type
  public  radtend_type
  public  intdiag_type

!--------------------------
!  public physics functions
!--------------------------
  public  initialize
  public  time_vary_step
  public  radiation_step1
  public  physics_step1
  public  physics_step2
  public  physics_end

CONTAINS

end module physics_abstraction_layer
