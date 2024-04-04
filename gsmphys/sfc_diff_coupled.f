      subroutine sfc_diff_coupled(im,ps,u1,v1,t1,q1,z1,
     &                    snwdph,tskin,z0rl,ztrl,cm,ch,rb,
     &                    prsl1,prslki,islimsk,
     &                    stress,fm,fh,
     &                    ustar,wind,ddvel,fm10,fh2,
     &                    sigmaf,vegtype,shdmax,ivegsrc,
     &                    tsurf,flag_iter) !,redrag,
!     &                    z0s_max)
!     &                    do_z0_moon, do_z0_hwrf15, do_z0_hwrf17,
!     &                    do_z0_hwrf17_hwonly, wind_th_hwrf)

      use machine , only : kind_phys
      use funcphys, only : fpvs    
      use physcons, grav => con_g,       cp => con_cp
     &,             rvrdm1 => con_fvirt, rd => con_rd
     &,             eps => con_eps, epsm1 => con_epsm1

      implicit none

! ---  input/output

      integer              im, ivegsrc

      real(kind=kind_phys), dimension(im)::ps,  u1, v1, t1, q1, z1
     &,                                    tskin, z0rl, ztrl, cm, ch, rb
     &,                                    prsl1, prslki, stress
     &,                                    fm, fh, ustar, wind, ddvel
     &,                                    fm10, fh2, sigmaf, shdmax
     &,                                    tsurf, snwdph
      integer, dimension(im)             ::vegtype, islimsk

      logical   flag_iter(im)
!      logical   redrag        
!      logical   do_z0_moon, do_z0_hwrf15, do_z0_hwrf17 ! kgao 
!     &,         do_z0_hwrf17_hwonly                    ! kgao 

! ---  local

      integer   i
!
      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,qs1,
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2, rat,
     &                     thv1,   tvs,    z1i,    z0, zt, z0max, ztmax,
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,
     &                     hl110,  hlt,    hltinf, olinf,
     &                     restar, czilc,  tem1,   tem2,
     &                     u10m, v10m, ws10m, ws10m_moon,         !kgao
     &                     z0_1, zt_1, fm1, fh1, ustar_1, ztmax_1 !kgao
!

!      real(kind=kind_phys),intent(in   ) :: z0s_max, wind_th_hwrf ! kgao 

      real(kind=kind_phys), parameter ::
     &              charnock=.014, ca=.4 
     &,             vis=1.4e-5, rnu=1.51e-5, visi=1.0/vis
     &,             log01=log(0.01), log05=log(0.05), log07=log(0.07)
     &,             ztmin1=-999.0

!================================================ 
! Main program starts here
!================================================ 

      do i=1,im

        if(flag_iter(i)) then 
 
! --- get variables at model lowest layer and surface (water/ice/land) 

          wind(i) = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))
     &                + max(0.0, min(ddvel(i), 30.0)), 1.0)
          tem1    = 1.0 + rvrdm1 * max(q1(i),1.e-8)
          thv1    = t1(i) * prslki(i) * tem1
          tvs     = 0.5 * (tsurf(i)+tskin(i)) * tem1
          qs1     = fpvs(t1(i))
          qs1     = max(1.0e-8, eps * qs1 / (prsl1(i) + epsm1 * qs1))

          !(sea/land/ice mask =0/1/2)
          if(islimsk(i) == 1 .or. islimsk(i) == 2) then ! over land or sea ice 

!================================================ 
! if over land or sea ice:
!    step 1 - get z0/zt 
!    step 2 - call similarity
!================================================
          
! --- get surface roughness for momentum (z0)

            z0      = 0.01 * z0rl(i) 
            z0max   = max(1.0e-6, min(z0,z1(i)))

            !xubin's new z0  over land and sea ice
            tem1 = 1.0 - shdmax(i) ! shdmax is max vegetation area fraction
            tem2 = tem1 * tem1
            tem1 = 1.0  - tem2
          
            if( ivegsrc == 1 ) then

              if (vegtype(i) == 10) then
                z0max = exp( tem2*log01 + tem1*log07 )
              elseif (vegtype(i) == 6) then
                z0max = exp( tem2*log01 + tem1*log05 )
              elseif (vegtype(i) == 7) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01
              elseif (vegtype(i) == 16) then
!               z0max = exp( tem2*log01 + tem1*log01 )
                z0max = 0.01
              else
                z0max = exp( tem2*log01 + tem1*log(z0max) )
              endif

            elseif (ivegsrc == 2 ) then

                if (vegtype(i) == 7) then
                  z0max = exp( tem2*log01 + tem1*log07 )
                elseif (vegtype(i) == 8) then
                  z0max = exp( tem2*log01 + tem1*log05 )
                elseif (vegtype(i) == 9) then
!                 z0max = exp( tem2*log01 + tem1*log01 )
                  z0max = 0.01
                elseif (vegtype(i) == 11) then
!                 z0max = exp( tem2*log01 + tem1*log01 )
                  z0max = 0.01
                else
                  z0max = exp( tem2*log01 + tem1*log(z0max) )
                endif

            z0max  = max(z0max,1.0e-6)

            endif

! --- get surface roughness for heat (zt)

!           czilc = 10.0 ** (- (0.40/0.07) * z0) ! let czilc depend on canopy height 
            czilc = 0.8

            tem1 = 1.0 - sigmaf(i)
            ztmax = z0max*exp( - tem1*tem1
     &                         * czilc*ca*sqrt(ustar(i)*(0.01/1.5e-05)))

            ztmax  = max(ztmax,1.0e-6)

! --- call similarity

            call monin_obukhov_similarity
     &       (z1(i), snwdph(i), thv1, wind(i), z0max, ztmax, tvs,
     &        rb(i), fm(i), fh(i), fm10(i), fh2(i),
     &        cm(i), ch(i), stress(i), ustar(i))

          elseif (islimsk(i) == 0) then

!================================================ 
! if over water 
! - redesigned by Kun Gao for coupling with MOM6
! - not updating z0, zt and ustar over ocean 
!================================================

            ! --- z0/zt from coupler
            z0      = 0.01 * z0rl(i)
            zt      = 0.01 * ztrl(i)

            z0max   = max(1.0e-6, min(z0,z1(i)))
            ztmax   = max(zt,1.0e-6)

            ! --- call similarity
            ! kgao: use z0/zt to do sfc diagnosis
            call monin_obukhov_similarity
     &       (z1(i), snwdph(i), thv1, wind(i), z0max, ztmax, tvs,
     &        rb(i), fm(i), fh(i), fm10(i), fh2(i),
     &        cm(i), ch(i), tem1, tem2) !stress(i), ustar(i))

            ! kgao: use ustar from coupler to get stress
            stress(i) =  ustar(i) * ustar(i)

          endif       ! end of if(islimsk) loop
        endif         ! end of if(flagiter) loop
      enddo           ! end of do i=1,im loop

      return
      end subroutine 
