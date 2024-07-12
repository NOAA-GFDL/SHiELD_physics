      !NOTE:
      ! This routine contains parameters added by Sofar Ocean to support
      ! additional coupling between the atmosphere, waves, and
      ! ocean models at high temporal and spatial resolutions.
      !
      ! Edits were made in 2023 by:
      ! Stephen G. Penny, Sofar Ocean (steve.penny@sofarocean.com)
      ! and
      ! Christie Hegermiller, Sofar Ocean

      subroutine sfc_diff_gfdl(im,ps,u1,v1,t1,q1,z1,
     &                    snwdph,tskin,z0rl,ztrl,cm,ch,rb,
     &                    prsl1,prslki,islimsk,
     &                    stress,fm,fh,
     &                    charnock,                             ! Sofar added Spring 2023
     &                    rhoa,                                 ! Sofar added 9/22/23
     &                    u10m_array,v10m_array,                ! Sofar added 11/17/23
     &                    u10n,v10n,                            ! Sofar added 9/22/23
     &                    fm10_neutral,                         ! Sofar added 10/19/23
     &                    ustar,wind,ddvel,fm10,fh2,
     &                    sigmaf,vegtype,shdmax,ivegsrc,
     &                    tsurf,flag_iter,redrag,
     &                    z0s_max,
     &                    do_z0_moon, do_z0_hwrf15, do_z0_hwrf17,
     &                    do_z0_hwrf17_hwonly, wind_th_hwrf)

! oct 2019 - a clean and updated version by Kun Gao at GFDL (Kun.Gao@noaa.gov)

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
     &,                                    fm, fh
     &,                                    charnock                     ! Sofar added Spring 2023
     &,                                    rhoa, u10n, v10n             ! Sofar added 9/22/23
     &,                                    u10m_array, v10m_array       ! Sofar added 11/17/23
     &,                                    ustar, wind 
     &,                                    ddvel
     &,                                    fm10, fh2, sigmaf, shdmax
     &,                                    tsurf, snwdph
     &,                                    fm_neutral, fm10_neutral     ! Sofar added Spring 2023
      real(kind=kind_phys) :: ws1, ws10n                                ! Sofar added 10/19/23
      integer, dimension(im)             ::vegtype, islimsk

      logical   flag_iter(im)
      logical   redrag        
      logical   do_z0_moon, do_z0_hwrf15, do_z0_hwrf17 ! kgao 
     &,         do_z0_hwrf17_hwonly                    ! kgao 

! ---  local

      integer   i
!
      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,qs1,
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2, rat,
     &                     thv1,   tvs,    z1i,    z0, zt, z0max, ztmax,
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,
     &                     tv1,                                   ! Sofar added 9/22/23
     &                     hl110,  hlt,    hltinf, olinf,
     &                     restar, czilc,  tem1,   tem2,
     &                     u10m, v10m, ws10m, ws10m_moon,         !kgao
     &                     z0_1, zt_1, fm1, fh1, ustar_1, ztmax_1 !kgao
!

      real(kind=kind_phys),intent(in   ) :: z0s_max, wind_th_hwrf ! kgao 

      real(kind=kind_phys), parameter ::
! CH     &              charnock=.014, ca=.4 
     &              ca=.4 
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
          tv1     = t1(i) * tem1                        ! Sofar added 9/22/23
          rhoa(i)   = prsl1(i) / (rd*tv1)               ! Sofar added 9/22/23
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
     &        fm_neutral(i), fm10_neutral(i),                          !(ADDED by Sofar)
     &        cm(i), ch(i), stress(i), ustar(i))

          elseif (islimsk(i) == 0) then ! over water

!================================================ 
! if over water (redesigned by Kun Gao)
!    iteration 1 
!         step 1 get z0/zt from previous step
!         step 2 call similarity
!    iteration 2 
!         step 1 update z0/zt 
!         step 2 call similarity 
!
! CH comment: Iteration 2 can be dropped because we have updated
! z0rl with wave-dependent roughness. STEVE: does this still apply (10/19/2023)?
!================================================

! === iteration 1

            ! --- get z0/zt
            z0      = 0.01 * z0rl(i)  
            zt      = 0.01 * ztrl(i)

            z0max   = max(1.0e-6, min(z0,z1(i)))                       !STEVE: note that z1 is the layer 1 height (m)
            ztmax   = max(zt,1.0e-6)

            ! --- call similarity
            call monin_obukhov_similarity
     &       (z1(i), snwdph(i), thv1, wind(i), z0max, ztmax, tvs,
     &        rb(i), fm(i), fh(i), fm10(i), fh2(i),
     &        fm_neutral(i), fm10_neutral(i),                          !(ADDED by Sofar)
     &        cm(i), ch(i), stress(i), ustar(i))

! === iteration 2

            ! --- get z0/zt following the old sfc_diff.f 
            z0 = (charnock(i) / grav) * ustar(i) * ustar(i)
            if (redrag) then
               z0 = max(min(z0, z0s_max), 1.e-7)
            else
               z0 = max(min(z0,.1), 1.e-7)
            endif

            ! zt calculations copied from old sfc_diff.f
            !ustar(i) = sqrt(grav * z0 / charnock(i))
            !restar = max(ustar(i)*z0max*visi, 0.000001)
            !rat    = min(7.0, 2.67 * sqrt(sqrt(restar)) - 2.57)
            !ztmax  = z0max * exp(-rat)

            ustar_1 = sqrt(grav * z0 / charnock(i))
            restar = max(ustar_1*z0max*visi, 0.000001)
            rat    = min(7.0, 2.67 * sqrt(sqrt(restar)) - 2.57)
            zt     = z0max * exp(-rat) ! zeng, zhao and dickinson 1997 (eq 25)

            ! --- update z0/zt with new options
            ! only z0 options in the following
            ! will add zt options in the future

            ! ---------------------------------------------- - Sofar (start)
            ! Compute neutral wind and u/v components at 10 m height
            ! Following IFS calculation, p.52:
            ! https://www.ecmwf.int/sites/default/files/elibrary/2021/20198-ifs-documentation-cy47r3-part-vi-physical-processes.pdf
            ! In the open sea, with fetch at least 5 km, rougness length is about 0.0002 (m)
            ! The equation below is for what ECMWF assumes for: z0M < 0.03 
            ! This may need to be updated for larger waves (z0M > 0.03)

            ! These didn't work:
!           u10n = u1(i) * fm10_neutral(i) / fm_neutral(i) !(ADDED by Sofar) CH: I don't think this is correct yet, unfortunately.
!           v10n = v1(i) * fm10_neutral(i) / fm_neutral(i) !(ADDED by Sofar) 
!!          u10n = u1(i) * fm10_neutral(i) / fm !_neutral(i) !(ADDED by Sofar)
!!          v10n = v1(i) * fm10_neutral(i) / fm !_neutral(i) !(ADDED by Sofar) SGP: see https://github.com/wavespotter/EarthSystemModel/issues/27

            ! Trying to match ECMWF closely:
            ! Compute lowest level wind speed
!           ws1 = sqrt(u1(i)*u1(i) + v1(i)*v1(i))
            ! Compute 10-meter equivalent neutral wind speed
            ! example: u10n = u1(i) * (ws10n/ws1) * fm10_neutral(i) / fm(i)
!           ws10n = ws1 * fm10_neutral(i) / fm_neutral(i) !(ADDED by Sofar)
!           ws10n = ws1 * fm10(i) / fm(i) !(ADDED by Sofar) TEST
!           ! Convert to components  by making use of the wind direction from the lowest model level
!           u10n = (u1(i) / ws1) * ws10n
!           v10n = (v1(i) / ws1) * ws10n

            ! Compute neutral winds for output to Sfcprop
            u10n(i) = u1(i) * fm10_neutral(i) / fm(i) 
            v10n(i) = v1(i) * fm10_neutral(i) / fm(i)

            ! Compute standard 10m winds for output to Sfcprop
            u10m_array(i) = u1(i) * fm10(i) / fm(i)
            v10m_array(i) = v1(i) * fm10(i) / fm(i)

            ! ---------------------------------------------- - Sofar (end)
              
            u10m = u1(i) * fm10(i) / fm(i)
            v10m = v1(i) * fm10(i) / fm(i)
            ws10m = sqrt(u10m*u10m + v10m*v10m)

            if (do_z0_hwrf15) then 
            ! option 1: HWRF15, originally developed by URI/GFDL 
                call cal_z0_hwrf15(ws10m, z0)
                call cal_zt_hwrf15(ws10m, zt)

            elseif (do_z0_hwrf17) then 
            ! option 2: HWRF17
                call cal_z0_hwrf17(ws10m, z0)
                call cal_zt_hwrf17(ws10m, zt)

            elseif (do_z0_hwrf17_hwonly) then 
            ! option 3: HWRF17 under high wind only
                if (ws10m > wind_th_hwrf) then
                   call cal_z0_hwrf17(ws10m, z0)
                   z0 = max(min(z0, z0s_max), 1.e-7) ! must apply limiter here
                endif

            elseif (do_z0_moon) then 
            ! option 4: Moon et al 2007 under high winds (same as in HiRAM)
               ws10m_moon = 2.458 + ustar(i)*(20.255-0.56*ustar(i))  ! Eq(7) Moon et al. 2007 
               if ( ws10m_moon > 20. ) then
                  call cal_z0_moon(ws10m_moon, z0, charnock(i))
                  z0 = max(min(z0, z0s_max), 1.e-7) ! must apply limiter here
               endif
            endif

            z0max = max(z0,1.0e-6)
            ztmax  = max(zt,1.0e-6)

            ! --- call similarity
            call monin_obukhov_similarity
     &       (z1(i), snwdph(i), thv1, wind(i), z0max, ztmax, tvs,
     &        rb(i), fm(i), fh(i), fm10(i), fh2(i),
     &        fm_neutral(i), fm10_neutral(i),                          !(ADDED by Sofar)
     &        cm(i), ch(i), stress(i), ustar(i))

            z0rl(i) = 100.0 * z0max
            ztrl(i) = 100.0 * ztmax

          endif       ! end of if(islimsk) loop
        endif         ! end of if(flagiter) loop
      enddo           ! end of do i=1,im loop

      return
      end subroutine sfc_diff_gfdl

! =======================================================================

      subroutine cal_z0_hwrf15(ws10m, z0)
      ! coded by Kun Gao (Kun.Gao@noaa.gov)
      ! originally developed by URI/GFDL
      use machine , only : kind_phys
      real(kind=kind_phys) :: ws10m, z0
      real(kind=kind_phys), parameter ::
     &  a0=-8.367276172397277e-12
     &, a1=1.7398510865876079e-09
     &, a2=-1.331896578363359e-07
     &, a3=4.507055294438727e-06
     &, a4=-6.508676881906914e-05
     &, a5=0.00044745137674732834
     &, a6=-0.0010745704660847233
     &, b0=2.1151080765239772e-13
     &, b1=-3.2260663894433345e-11
     &, b2=-3.329705958751961e-10
     &, b3=1.7648562021709124e-07
     &, b4=7.107636825694182e-06
     &, b5=-0.0013914681964973246
     &, b6=0.0406766967657759

      if (ws10m <= 5.0) then
         z0 = 0.0185/9.8*(7.59e-4*ws10m**2+2.46e-2*ws10m)**2
      elseif (ws10m > 5.0 .and. ws10m <= 10.) then
         z0 = 0.00000235*(ws10m**2-25.)+3.805129199617346e-05
      elseif (ws10m > 10.0 .and. ws10m <= 60.) then
         z0 = a6 + a5*ws10m + a4*ws10m**2 + a3*ws10m**3
     &      + a2*ws10m**4 + a1*ws10m**5 + a0*ws10m**6
      else
         z0 = b6 + b5*ws10m + b4*ws10m**2 + b3*ws10m**3
     &      + b2*ws10m**4 + b1*ws10m**5 + b0*ws10m**6
      endif

      end subroutine cal_z0_hwrf15

      subroutine cal_zt_hwrf15(ws10m, zt)
      ! coded by Kun Gao (Kun.Gao@noaa.gov)
      ! originally developed by URI/GFDL
      use machine , only : kind_phys
      real(kind=kind_phys) :: ws10m, zt
      real(kind=kind_phys), parameter ::
     &  a0 = 2.51715926619e-09
     &, a1 = -1.66917514012e-07
     &, a2 = 4.57345863551e-06
     &, a3 = -6.64883696932e-05
     &, a4 = 0.00054390175125
     &, a5 = -0.00239645231325
     &, a6 = 0.00453024927761
     &, b0 = -1.72935914649e-14
     &, b1 = 2.50587455802e-12
     &, b2 = -7.90109676541e-11
     &, b3 = -4.40976353607e-09
     &, b4 = 3.68968179733e-07
     &, b5 = -9.43728336756e-06
     &, b6 = 8.90731312383e-05
     &, c0 = 4.68042680888e-14
     &, c1 = -1.98125754931e-11
     &, c2 = 3.41357133496e-09
     &, c3 = -3.05130605309e-07
     &, c4 = 1.48243563819e-05
     &, c5 = -0.000367207751936
     &, c6 = 0.00357204479347

      if (ws10m <= 7.0) then
         zt = 0.0185/9.8*(7.59e-4*ws10m**2+2.46e-2*ws10m)**2
      elseif (ws10m > 7.0 .and. ws10m <= 15.) then
         zt = a6 + a5*ws10m + a4*ws10m**2 + a3*ws10m**3
     &      + a2*ws10m**4 + a1*ws10m**5 + a0*ws10m**6
      elseif (ws10m > 15.0 .and. ws10m <= 60.) then
         zt = b6 + b5*ws10m + b4*ws10m**2 + b3*ws10m**3
     &      + b2*ws10m**4 + b1*ws10m**5 + b0*ws10m**6
      else
         zt = c6 + c5*ws10m + c4*ws10m**2 + c3*ws10m**3
     &      + c2*ws10m**4 + c1*ws10m**5 + c0*ws10m**6
      endif
      end subroutine cal_zt_hwrf15

! =======================================================================

      subroutine cal_z0_hwrf17(ws10m, z0)
      ! coded by Kun Gao (Kun.Gao@noaa.gov)
      use machine , only : kind_phys
      real(kind=kind_phys) :: ws10m, z0
      real(kind=kind_phys), parameter ::
     &  p13=-1.296521881682694e-02
     &, p12= 2.855780863283819e-01
     &, p11=-1.597898515251717e+00
     &, p10=-8.396975715683501e+00
     &, p25= 3.790846746036765e-10
     &, p24= 3.281964357650687e-09
     &, p23= 1.962282433562894e-07
     &, p22=-1.240239171056262e-06
     &, p21=1.739759082358234e-07
     &, p20=2.147264020369413e-05
     &, p35=1.840430200185075e-07
     &, p34=-2.793849676757154e-05
     &, p33=1.735308193700643e-03
     &, p32=-6.139315534216305e-02
     &, p31=1.255457892775006e+00
     &, p30=-1.663993561652530e+01
     &, p40=4.579369142033410e-04

      if (ws10m <= 6.5) then
         z0 = exp( p10 + p11*ws10m + p12*ws10m**2 + p13*ws10m**3)
      elseif (ws10m > 6.5 .and. ws10m <= 15.7) then
         z0 = p25*ws10m**5 + p24*ws10m**4 + p23*ws10m**3 
     &      + p22*ws10m**2 + p21*ws10m + p20
      elseif (ws10m > 15.7 .and. ws10m <= 53.) then
         z0 = exp( p35*ws10m**5 + p34*ws10m**4 + p33*ws10m**3 
     &      + p32*ws10m**2 + p31*ws10m + p30 )
      else
         z0 = p40
      endif
      end subroutine cal_z0_hwrf17

      subroutine cal_zt_hwrf17(ws10m, zt)
      ! coded by Kun Gao (Kun.Gao@noaa.gov)
      use machine , only : kind_phys
      real(kind=kind_phys) :: ws10m, zt
      real(kind=kind_phys), parameter  :: p00 =  1.100000000000000e-04,
     &      p15 = -9.144581627678278e-10, p14 =  7.020346616456421e-08,
     &      p13 = -2.155602086883837e-06, p12 =  3.333848806567684e-05,
     &      p11 = -2.628501274963990e-04, p10 =  8.634221567969181e-04,
     &      p25 = -8.654513012535990e-12, p24 =  1.232380050058077e-09,
     &      p23 = -6.837922749505057e-08, p22 =  1.871407733439947e-06,
     &      p21 = -2.552246987137160e-05, p20 =  1.428968311457630e-04,
     &      p35 =  3.207515102100162e-12, p34 = -2.945761895342535e-10,
     &      p33 =  8.788972147364181e-09, p32 = -3.814457439412957e-08,
     &      p31 = -2.448983648874671e-06, p30 =  3.436721779020359e-05,
     &      p45 = -3.530687797132211e-11, p44 =  3.939867958963747e-09,
     &      p43 = -1.227668406985956e-08, p42 = -1.367469811838390e-05,
     &      p41 =  5.988240863928883e-04, p40 = -7.746288511324971e-03,
     &      p56 = -1.187982453329086e-13, p55 =  4.801984186231693e-11,
     &      p54 = -8.049200462388188e-09, p53 =  7.169872601310186e-07,
     &      p52 = -3.581694433758150e-05, p51 =  9.503919224192534e-04,
     &      p50 = -1.036679430885215e-02,
     &      p60 =  4.751256171799112e-05

      if (ws10m >= 0.0 .and. ws10m < 5.9 ) then
         zt = p00
      elseif (ws10m >= 5.9 .and. ws10m <= 15.4) then
         zt = p10 + ws10m * (p11 + ws10m * (p12 + ws10m * (p13
     &               + ws10m * (p14 + ws10m * p15))))
      elseif (ws10m > 15.4 .and. ws10m <= 21.6) then
         zt = p20 + ws10m * (p21 + ws10m * (p22 + ws10m * (p23
     &               + ws10m * (p24 + ws10m * p25))))
      elseif (ws10m > 21.6 .and. ws10m <= 42.2) then
         zt = p30 + ws10m * (p31 + ws10m * (p32 + ws10m * (p33
     &               + ws10m * (p34 + ws10m * p35))))
      elseif ( ws10m > 42.2 .and. ws10m <= 53.3) then
         zt = p40 + ws10m * (p41 + ws10m * (p42 + ws10m * (p43
     &               + ws10m * (p44 + ws10m * p45))))
      elseif ( ws10m > 53.3 .and. ws10m <= 80.0) then
         zt = p50 + ws10m * (p51 + ws10m * (p52 + ws10m * (p53
     &               + ws10m * (p54 + ws10m * (p55 + ws10m * p56)))))
      elseif ( ws10m > 80.0) then
         zt = p60
      endif
      end subroutine cal_zt_hwrf17

! =======================================================================

      subroutine cal_z0_moon(ws10m, z0, charnock)
      ! coded by Kun Gao (Kun.Gao@noaa.gov)
      use machine , only : kind_phys
      use physcons, grav => con_g

      real(kind=kind_phys) :: ws10m, z0
      real(kind=kind_phys) :: ustar_th, z0_adj 

      real(kind=kind_phys), parameter ::
!     &          charnock=.014
     &          wind_th_moon = 20. 
     &,         a = 0.56
     &,         b = -20.255
     &,         c = wind_th_moon - 2.458

      ustar_th = (-b-sqrt(b*b-4*a*c))/(2*a)

      z0_adj = 0.001*(0.085*wind_th_moon - 0.58) - 
     &        (charnock/grav)*ustar_th*ustar_th 

      z0 = 0.001*(0.085*ws10m - 0.58) - z0_adj  ! Eq(8b) Moon et al. 2007 modified by kgao

      end subroutine cal_z0_moon

! =======================================================================

      subroutine monin_obukhov_similarity
     &     ( z1, snwdph, thv1, wind, z0max, ztmax, tvs,
     &       rb, fm, fh, fm10, fh2, 
     &       fm_neutral, fm10_neutral,                                 !(ADDED by Sofar)
     &       cm, ch, stress, ustar)

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
! fm_neutral, fm10_neutral - similarity functions for neutral profile  !(ADDED by Sofar)
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
     &       rb, fm, fh, fm10, fh2, 
     &       fm_neutral, fm10_neutral,                                 !(ADDED by Sofar)
     &       cm, ch, stress, ustar

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

          fm_neutral = fm                                              !(ADDED by Sofar)
          fm10_neutral = fm10                                          !(ADDED by Sofar)
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
