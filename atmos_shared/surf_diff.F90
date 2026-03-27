module surf_diff_mod

! Mostly taken from vert_diff in atmos_phys
! could be replaced by the original routines
! when merging SHiELD and AMx
! Scientific description can be found in the 
! document vert_diff.tech.ps dated March 2001
! joseph.mouallem@noaa.gov
! March 2026


use constants_mod, only:  GRAV, RDGAS, RVGAS, CP_AIR
use IPD_typedefs,           only: IPD_data_type, kind_phys
use fv_iau_mod,             only: IAU_external_data_type
use block_control_mod,      only: block_control_type

implicit none

private

public :: surf_diff_type, compute_nu, compute_e, vert_diff_down_2

  real(kind=kind_phys), parameter :: d608 = (RVGAS-RDGAS)/RDGAS
  logical :: use_virtual_temp_vert_diff = .true.
  logical :: do_mcm_plev                = .false.



type surf_diff_type

  real, pointer, dimension(:,:) :: dtmass  => NULL(),   &
                                   dflux_t => NULL(),   &
                                   delta_t => NULL(),   &
                                   delta_u => NULL(),   &
                                   delta_v => NULL()
  real, pointer, dimension(:,:,:) :: tdt_dyn => NULL(), &
                                     qdt_dyn => NULL(), &
                                     dgz_dyn => NULL(), &
                                     ddp_dyn => NULL(), &
                                     tdt_rad => NULL()   !miz

  real, pointer, dimension(:,:,:) :: dflux_tr => NULL(),& ! tracer flux tendency
                                     delta_tr => NULL()   ! tracer tendency
end type surf_diff_type
!
!
contains
!
!

! inspired by atmos_phys/atmos_param/vert_diff/vert_diff.F90
subroutine compute_nu (diff, p_half, p_full, z_full, t, q, nu)

!-----------------------------------------------------------------------
real,    intent(in)    , dimension(:,:,:) :: diff, p_half, p_full, &
                                             z_full, t, q
real,    intent(out)   , dimension(:,:,:) :: nu

real, dimension(size(t,1),size(t,2),size(t,3)) :: rho_half, tt
integer :: nlev
!-----------------------------------------------------------------------

nlev = size(nu,3)

if ( use_virtual_temp_vert_diff ) then
  tt = t * (1.0 + d608*q)           ! virtual temperature
else
  tt = t ! Take out virtual temperature effect here to mimic supersource
endif

rho_half(:,:,2:nlev) =  &         ! density at half levels
      2.0*p_half(:,:,2:nlev)/(RDGAS*(tt(:,:,2:nlev)+tt(:,:,1:nlev-1)))

if(do_mcm_plev) then
  nu(:,:,2:nlev) = GRAV*rho_half(:,:,2:nlev)*rho_half(:,:,2:nlev)*diff(:,:,2:nlev)/ &
                    (p_full(:,:,2:nlev)-p_full(:,:,1:nlev-1))
else
  nu(:,:,2:nlev) = rho_half(:,:,2:nlev)*diff(:,:,2:nlev) /  &
                    (z_full(:,:,1:nlev-1)-z_full(:,:,2:nlev))
endif
!-----------------------------------------------------------------------

end subroutine compute_nu

subroutine compute_e (delt, mu, nu, e, a, b, c, g)

!-----------------------------------------------------------------------

real,    intent(in)                       :: delt
real,    intent(in)    , dimension(:,:,:) :: mu, nu
real,    intent(out)   , dimension(:,:,:) :: e, a, b, c, g

integer :: k, nlev

!-----------------------------------------------------------------------

 nlev = size(mu,3)

 a(:,:,1:nlev-1) = - mu(:,:,1:nlev-1)*nu(:,:,2:nlev)*delt
 a(:,:,nlev    ) =   0.0
 c(:,:,2:nlev  ) = - mu(:,:,2:nlev  )*nu(:,:,2:nlev)*delt
 c(:,:,1       ) =   0.0

 b = 1.0 - a - c

 e(:,:,1)   =   - a(:,:,1)/b(:,:,1)
 do  k= 2, nlev - 1
    g(:,:,k) = 1.0/(b(:,:,k) + c(:,:,k)*e(:,:,k-1))
    e(:,:,k) = - a(:,:,k)*g(:,:,k)
 enddo

!-----------------------------------------------------------------------

end subroutine compute_e

subroutine compute_f (dt_xi, b, c, g, f)

!-----------------------------------------------------------------------
real,    intent(in)    , dimension(:,:,:) :: dt_xi, b, c, g
real,    intent(out)   , dimension(:,:,:) :: f

integer :: k
!-----------------------------------------------------------------------

 f(:,:,1) =   dt_xi(:,:,1)/b(:,:,1)

 do  k = 2, size(b,3)-1
    f(:,:,k) = (dt_xi(:,:,k) - c(:,:,k)*f(:,:,k-1))*g(:,:,k)
 enddo

!-----------------------------------------------------------------------

end subroutine compute_f

subroutine explicit_tend (mu, nu, xi, dt_xi)

!-----------------------------------------------------------------------

real,    intent(in)    , dimension(:,:,:) :: mu, nu, xi
real,    intent(inout) , dimension(:,:,:) :: dt_xi

real, dimension(size(xi,1),size(xi,2),size(xi,3)) :: fluxx

integer :: nlev

!-----------------------------------------------------------------------

 nlev = size(mu,3)

 fluxx(:,:,1)      = 0.0
 fluxx(:,:,2:nlev) = nu(:,:,2:nlev)*(xi(:,:,2:nlev) - xi(:,:,1:nlev-1))

 dt_xi(:,:,1:nlev-1) = dt_xi(:,:,1:nlev-1) +  &
    mu(:,:,1:nlev-1)*(fluxx(:,:,2:nlev) - fluxx(:,:,1:nlev-1))
 dt_xi(:,:,nlev)     = dt_xi(:,:,nlev) - mu(:,:,nlev)*fluxx(:,:,nlev)

!-----------------------------------------------------------------------

end subroutine explicit_tend

subroutine vert_diff_down_2 &
      (delt, mu, nu, xi_1, xi_2, dt_xi_1, dt_xi_2, e, f_1, f_2, &
       mu_delt_n, nu_n, e_n1, f_1_delt_n1, f_2_delt_n1,         &
       delta_1_n, delta_2_n, kbot)

!-----------------------------------------------------------------------

real,    intent(in)                       :: delt
real,    intent(in)    , dimension(:,:,:) :: mu, nu, xi_1, xi_2
real,    intent(in)    , dimension(:,:,:) :: dt_xi_1, dt_xi_2
real,    intent(out)   , dimension(:,:,:) :: e, f_1, f_2
real,    intent(out)   , dimension(:,:)   :: mu_delt_n, nu_n, e_n1,    &
                                             f_1_delt_n1, f_2_delt_n1, &
                                             delta_1_n, delta_2_n

integer, intent(in),    dimension(:,:), optional :: kbot

real, dimension(size(xi_1,1),size(xi_1,2),size(xi_1,3)) :: a, b, c, g, &
                                                      dt_xi_11, dt_xi_22

integer :: i, j, kb, nlev

!-----------------------------------------------------------------------

! local copy of input 
  dt_xi_11 = dt_xi_1
  dt_xi_22 = dt_xi_2

 call explicit_tend (mu, nu, xi_1, dt_xi_11)
 call explicit_tend (mu, nu, xi_2, dt_xi_22)

 call compute_e (delt, mu, nu, e, a, b, c, g)

 call compute_f (dt_xi_11, b, c, g, f_1)
 call compute_f (dt_xi_22, b, c, g, f_2)

 if (present(kbot)) then
    do j=1,size(xi_1,2)
    do i=1,size(xi_1,1)
        kb = kbot(i,j)
        mu_delt_n(i,j)  =      mu(i,j,kb  )*delt
             nu_n(i,j)  =      nu(i,j,kb  )
            e_n1(i,j)  =       e(i,j,kb-1)
     f_1_delt_n1(i,j)  =     f_1(i,j,kb-1)*delt
     f_2_delt_n1(i,j)  =     f_2(i,j,kb-1)*delt
        delta_1_n(i,j)  = dt_xi_11(i,j,kb  )*delt
        delta_2_n(i,j)  = dt_xi_22(i,j,kb  )*delt
    enddo
    enddo
 else
        nlev = size(mu,3)
        mu_delt_n(:,:)  =      mu(:,:,nlev  )*delt
             nu_n(:,:)  =      nu(:,:,nlev  )
            e_n1(:,:)  =       e(:,:,nlev-1)
     f_1_delt_n1(:,:)  =     f_1(:,:,nlev-1)*delt
     f_2_delt_n1(:,:)  =     f_2(:,:,nlev-1)*delt
        delta_1_n(:,:)  = dt_xi_11(:,:,nlev  )*delt
        delta_2_n(:,:)  = dt_xi_22(:,:,nlev  )*delt
 endif

!-----------------------------------------------------------------------

end subroutine vert_diff_down_2

end module
