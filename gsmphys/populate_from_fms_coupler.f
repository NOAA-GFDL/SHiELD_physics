! This routine is a combination of sfc_ocean_coupled and sfc_diff_coupled
! It basically populate the needed variables from the full coupler output
! such as eavp, hflx, qsurf, stress to drive the physics bypassing all the
! old caluclation at the surface from the physics ocean model, land model, 
! ice and others (these are located within the do iter=1,2 loop in the 
! GFS_physics_driver)

! Joseph M.
! 09/05/2025

      subroutine populate_from_fms_coupler
!  ---  inputs from sfc_ocean_coupled:
     &     ( im, ps, u1, v1, t1, q1, tskin, cm, ch,                     &
     &       prsl1, prslki, ddvel, maxevap,                             &
     &       shflx, lhflx, qsfc, ! kgao: shflx and lhflx from coupler
! ---- additional input from sfc_diff_coupled
     &       z1, snwdph, z0rl, ztrl, rb, stress, fm, fh, ustar,         &
     &       wind, fm10, fh2, tsurf,                                    &
!  ---  outputs:
     &       qsurf, cmm, chh, gflux, evap, hflx, ep                     &
     &     )

      use machine , only : kind_phys
      use funcphys, only : fpvs
      use physcons, only : cp => con_cp, rd => con_rd, eps => con_eps,  &
     &                     epsm1 => con_epsm1,                          &
     &                     rvrdm1 => con_fvirt
!
      implicit none
!
!  ---  constant parameters:
      real (kind=kind_phys), parameter :: cpinv  = 1.0/cp               &

!  ---  inputs:
      integer, intent(in) :: im

      real (kind=kind_phys), dimension(im), intent(in) :: ps, u1, v1,   &
     &      t1, q1, tskin, cm, ch, prsl1, prslki, ddvel, maxevap

      real (kind=kind_phys), dimension(im), intent(in) :: shflx,        &
     &        lhflx, qsfc ! kgao

      real(kind=kind_phys), dimension(im):: z1, z0rl, ztrl, rb
     &,                                    stress, fm, fh, ustar, wind
     &,                                    fm10, fh2, tsurf, snwdph

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(out) :: qsurf,       &
     &       cmm, chh, gflux, evap, hflx, ep 

!  ---  locals:

      real(kind=kind_phys) :: q0, qss_local, rch, rho,  tem, qs1, thv1,
     &               tvs,    z0, zt, z0max, ztmax, tem1,   tem2

      integer :: i

!!!!!!!!!!!!!!!!!!!!!!!!!
! FROM sfc_coupled_ocean
!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, im

          wind(i)     = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))                &
     &                 + max( 0.0, min( ddvel(i), 30.0 ) ), 1.0)

          q0       = max( q1(i), 1.0e-8 )
          rho      = prsl1(i) / (rd*t1(i)*(1.0 + rvrdm1*q0))

          qss_local      = fpvs( tskin(i) )
          qss_local      = eps*qss_local / (ps(i) + epsm1*qss_local)

          evap(i)  = 0.0
          hflx(i)  = 0.0
          ep(i)    = 0.0
          gflux(i) = 0.0

          rch      = rho * cp * ch(i) * wind(i)
          cmm(i)   = cm(i) * wind(i)
          chh(i)   = rho * ch(i) * wind(i)

          tem      = 1.0 / rho
          qsurf(i) = qss_local

          hflx(i)  = shflx(i) * tem * cpinv
          evap(i)  = lhflx(i) * tem
          qsurf(i) = qsfc(i)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     FROM surf_diff_coupled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
          tem1    = 1.0 + rvrdm1 * max(q1(i),1.e-8)
          thv1    = t1(i) * prslki(i) * tem1
          tvs     = 0.5 * (tsurf(i)+tskin(i)) * tem1
          qs1     = fpvs(t1(i))
          qs1     = max(1.0e-8, eps * qs1 / (prsl1(i) + epsm1 * qs1))

          ! --- z0/zt from coupler
          z0      = 0.01 * z0rl(i)
          zt      = 0.01 * ztrl(i)

          z0max   = max(1.0e-6, min(z0,z1(i)))
          ztmax   = max(zt,1.0e-6)

          ! --- call similarity
          ! kgao: use z0/zt to do sfc diagnosis
          call monin_obukhov_similarity
     &     (z1(i), snwdph(i), thv1, wind(i), z0max, ztmax, tvs,
     &      rb(i), fm(i), fh(i), fm10(i), fh2(i),
     &      cm(i), ch(i), tem1, tem2) !stress(i), ustar(i))

          ! kgao: use ustar from coupler to get stress
          stress(i) =  ustar(i) * ustar(i)

      enddo           ! end of do i=1,im loop

      return

      end subroutine populate_from_fms_coupler


