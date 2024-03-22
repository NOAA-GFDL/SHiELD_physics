      subroutine sfc_diff_gfdl_coupled(im,ps,u1,v1,t1,q1,z1,
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
! redesigned by Kun Gao for coupling with MOM6
! now do not update z0rl and ustar over ocean points 
!================================================

            ! --- get z0/zt
            z0      = 0.01 * z0rl(i) ! kgao: from coupler
            zt      = 0.01 * ztrl(i)

            z0max   = max(1.0e-6, min(z0,z1(i)))
            ztmax   = max(zt,1.0e-6)

            ! --- call similarity
            call monin_obukhov_similarity
     &       (z1(i), snwdph(i), thv1, wind(i), z0max, ztmax, tvs,
     &        rb(i), fm(i), fh(i), fm10(i), fh2(i),
     &        cm(i), ch(i), tem1, tem2) !stress(i), ustar(i))

            ! kgao: use ustar from coupler to update stress
            stress(i) =  ustar(i) * ustar(i)

            ! kgao: get zt as in old sfc_diff.f
            ! this should be removed once using ztrl from coupler  
            ustar_1 = sqrt(grav * z0max / charnock)
            restar = max(ustar_1*z0max*visi, 0.000001)
            rat    = min(7.0, 2.67 * sqrt(sqrt(restar)) - 2.57)
            zt     = z0max * exp(-rat) ! zeng, zhao and dickinson 1997 (eq 25)

            ! kgao: do not update z0 
            !z0rl(i) = 100.0 * z0max
            ztrl(i) = 100.0 * ztmax

          endif       ! end of if(islimsk) loop
        endif         ! end of if(flagiter) loop
      enddo           ! end of do i=1,im loop

      return
      end subroutine 

! =======================================================================

      subroutine monin_obukhov_similarity
     &     ( z1, snwdph, thv1, wind, z0max, ztmax, tvs,
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar)

! --- input 
! z1     - lowest model level height 
! snwdph - surface snow thickness 
! wind   - wind speed at lowest model layer
! thv1   - virtual potential temp at lowest model layer
! tvs    - surface temp
! z0max  - surface roughness length for momentum
! ztmax  - surface roughness length for heat
!
! --- output 
! rb        - a bulk richardson number
! fm, fh    - similarity function defined at lowest model layer 
! fm10, fh2 - similarity function defined at 10m (for momentum) and 2m (for heat)
! cm, ch    - surface exchange coefficients for momentum and heat
! stress    - surface wind stress
! ustar     - surface frictional velocity 

      use machine , only : kind_phys
      use physcons, grav => con_g

!  ---  inputs:
      real(kind=kind_phys), intent(in) :: 
     &       z1, snwdph, thv1, wind, z0max, ztmax, tvs

!  ---  outputs:
      real(kind=kind_phys), intent(out) ::
     &       rb, fm, fh, fm10, fh2, cm, ch, stress, ustar

!  ---  locals:

      real(kind=kind_phys), parameter :: alpha=5., a0=-3.975
     &,             a1=12.32, alpha4=4.0*alpha
     &,             b1=-7.755,  b2=6.041,  alpha2=alpha+alpha, beta=1.0
     &,             a0p=-7.941, a1p=24.75, b1p=-8.705, b2p=7.899
     &,             ztmin1=-999.0, ca=.4

      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2,
     &                     z1i,
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,
     &                     hl110,  hlt,    hltinf, olinf,
     &                     tem1,   tem2,   ztmax1

          z1i = 1.0 / z1

          tem1   = z0max/z1
          if (abs(1.0-tem1) > 1.0e-6) then
            ztmax1 = - beta*log(tem1)/(alpha2*(1.-tem1))
          else
            ztmax1 = 99.0
          endif
          if( z0max < 0.05 .and. snwdph < 10.0 ) ztmax1 = 99.0

!
!  compute stability indices (rb and hlinf)
!
          dtv     = thv1 - tvs
          adtv    = max(abs(dtv),0.001)
          dtv     = sign(1.,dtv) * adtv
          rb      = max(-5000.0, (grav+grav) * dtv * z1
     &            / ((thv1 + tvs) * wind * wind))
          tem1    = 1.0 / z0max
          tem2    = 1.0 / ztmax
          fm      = log((z0max+z1)  * tem1)
          fh      = log((ztmax+z1)  * tem2)
          fm10    = log((z0max+10.) * tem1)
          fh2     = log((ztmax+2.)  * tem2)
          hlinf   = rb * fm * fm / fh
          hlinf   = min(max(hlinf,ztmin1),ztmax1)
!
!  stable case
!
          if (dtv >= 0.0) then
            hl1 = hlinf
            if(hlinf > .25) then
              tem1   = hlinf * z1i
              hl0inf = z0max * tem1
              hltinf = ztmax * tem1
              aa     = sqrt(1. + alpha4 * hlinf)
              aa0    = sqrt(1. + alpha4 * hl0inf)
              bb     = aa
              bb0    = sqrt(1. + alpha4 * hltinf)
              pm     = aa0 - aa + log( (aa + 1.)/(aa0 + 1.) )
              ph     = bb0 - bb + log( (bb + 1.)/(bb0 + 1.) )
              fms    = fm - pm
              fhs    = fh - ph
              hl1    = fms * fms * rb / fhs
              hl1    = min(max(hl1, ztmin1), ztmax1)
            endif
!
!  second iteration
!
            tem1  = hl1 * z1i
            hl0   = z0max * tem1
            hlt   = ztmax * tem1
            aa    = sqrt(1. + alpha4 * hl1)
            aa0   = sqrt(1. + alpha4 * hl0)
            bb    = aa
            bb0   = sqrt(1. + alpha4 * hlt)
            pm    = aa0 - aa + log( (1.0+aa)/(1.0+aa0) )
            ph    = bb0 - bb + log( (1.0+bb)/(1.0+bb0) )
            hl110 = hl1 * 10. * z1i
            hl110 = min(max(hl110, ztmin1), ztmax1)
            aa    = sqrt(1. + alpha4 * hl110)
            pm10  = aa0 - aa + log( (1.0+aa)/(1.0+aa0) )
            hl12  = (hl1+hl1) * z1i
            hl12  = min(max(hl12,ztmin1),ztmax1)
!           aa    = sqrt(1. + alpha4 * hl12)
            bb    = sqrt(1. + alpha4 * hl12)
            ph2   = bb0 - bb + log( (1.0+bb)/(1.0+bb0) )
!
!  unstable case - check for unphysical obukhov length
!
          else                          ! dtv < 0 case
            olinf = z1 / hlinf
            tem1  = 50.0 * z0max
            if(abs(olinf) <= tem1) then
              hlinf = -z1 / tem1
              hlinf = min(max(hlinf,ztmin1),ztmax1)
            endif
!
!  get pm and ph
!
            if (hlinf >= -0.5) then
              hl1   = hlinf
              pm    = (a0  + a1*hl1)  * hl1   / (1.+ (b1+b2*hl1)  *hl1)
              ph    = (a0p + a1p*hl1) * hl1   / (1.+ (b1p+b2p*hl1)*hl1)
              hl110 = hl1 * 10. * z1i
              hl110 = min(max(hl110, ztmin1), ztmax1)
              pm10  = (a0 + a1*hl110) * hl110 / (1.+(b1+b2*hl110)*hl110)
              hl12  = (hl1+hl1) * z1i
              hl12  = min(max(hl12, ztmin1), ztmax1)
              ph2   = (a0p + a1p*hl12) * hl12 / (1.+(b1p+b2p*hl12)*hl12)
            else                       ! hlinf < 0.05
              hl1   = -hlinf
              tem1  = 1.0 / sqrt(hl1)
              pm    = log(hl1) + 2. * sqrt(tem1) - .8776
              ph    = log(hl1) + .5 * tem1 + 1.386
!             pm    = log(hl1) + 2.0 * hl1 ** (-.25) - .8776
!             ph    = log(hl1) + 0.5 * hl1 ** (-.5) + 1.386
              hl110 = hl1 * 10. * z1i
              hl110 = min(max(hl110, ztmin1), ztmax1)
              pm10  = log(hl110) + 2.0 / sqrt(sqrt(hl110)) - .8776
!             pm10  = log(hl110) + 2. * hl110 ** (-.25) - .8776
              hl12  = (hl1+hl1) * z1i
              hl12  = min(max(hl12, ztmin1), ztmax1)
              ph2   = log(hl12) + 0.5 / sqrt(hl12) + 1.386
!             ph2   = log(hl12) + .5 * hl12 ** (-.5) + 1.386
            endif

          endif          ! end of if (dtv >= 0 ) then loop
!
!  finish the exchange coefficient computation to provide fm and fh
!
          fm        = fm - pm
          fh        = fh - ph
          fm10      = fm10 - pm10
          fh2       = fh2 - ph2
          cm        = ca * ca / (fm * fm)
          ch        = ca * ca / (fm * fh)
          tem1      = 0.00001/z1
          cm        = max(cm, tem1)
          ch        = max(ch, tem1)
          stress    = cm * wind * wind
          ustar     = sqrt(stress)

      return
      end subroutine monin_obukhov_similarity
