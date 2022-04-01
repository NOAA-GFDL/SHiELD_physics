      subroutine sfc_diff(im,ps,u1,v1,t1,q1,z1,
     &                    snwdph,tskin,z0rl,cm,ch,rb,
     &                    prsl1,prslki,islimsk,
     &                    stress,fm,fh,
     &                    ustar,wind,ddvel,fm10,fh2,
     &                    sigmaf,vegtype,shdmax,ivegsrc,
     &                    tsurf,flag_iter,redrag,czilc,
     &                    z0s_max,
     &                    do_z0_moon, do_z0_hwrf15, do_z0_hwrf17,
     &                    do_z0_hwrf17_hwonly, wind_th_hwrf)
!
      use machine , only : kind_phys
      use funcphys, only : fpvs    
      use physcons, grav => con_g,       cp => con_cp
     &,             rvrdm1 => con_fvirt, rd => con_rd
     &,             eps => con_eps, epsm1 => con_epsm1

      implicit none
!
      integer              im, ivegsrc
      real(kind=kind_phys), dimension(im) :: ps,  u1, v1, t1, q1, z1
     &,                                      tskin, z0rl, cm,  ch, rb
     &,                                      prsl1, prslki, stress
     &,                                      fm, fh, ustar, wind, ddvel
     &,                                      fm10, fh2, sigmaf, shdmax
     &,                                      tsurf, snwdph
      integer, dimension(im)              ::  vegtype, islimsk
      real(kind=kind_phys) czilc

      logical   flag_iter(im) ! added by s.lu
      logical   redrag        ! reduced drag coeff. flag for high wind over sea (j.han)
      logical   do_z0_moon, do_z0_hwrf15, do_z0_hwrf17 
     &,         do_z0_hwrf17_hwonly ! added by kgao 
!
!     locals
!
      integer   i
!
      real(kind=kind_phys) aa,     aa0,    bb,     bb0, dtv,   adtv,qs1,
     &                     hl1,    hl12,   pm,     ph,  pm10,  ph2, rat,
     &                     thv1,   tvs,    z1i,    z0,  z0max, ztmax,
     &                     fms,    fhs,    hl0,    hl0inf, hlinf,
     &                     hl110,  hlt,    hltinf, olinf,
     &                     restar, tem1,   tem2, ztmax1,
     &                     z0_adj, wind_th_moon, ustar_th, a,b,c, !kgao 
     &                     u10m, v10m, ws10m !kgao
!

      real(kind=kind_phys),intent(in   ) :: z0s_max, wind_th_hwrf ! kgao 

      real(kind=kind_phys), parameter ::
     &              charnock=.014, ca=.4  ! ca - von karman constant
     &,             alpha=5.,   a0=-3.975, a1=12.32, alpha4=4.0*alpha
     &,             b1=-7.755,  b2=6.041,  alpha2=alpha+alpha, beta=1.0
     &,             a0p=-7.941, a1p=24.75, b1p=-8.705, b2p=7.899
     &,             vis=1.4e-5, rnu=1.51e-5, visi=1.0/vis

     &,             log01=log(0.01), log05=log(0.05), log07=log(0.07)
     &,             ztmin1=-999.0
! following is added by kgao 
     &,         bs0=-8.367276172397277e-12 
     &,         bs1=1.7398510865876079e-09
     &,         bs2=-1.331896578363359e-07
     &,         bs3=4.507055294438727e-06
     &,         bs4=-6.508676881906914e-05
     &,         bs5=0.00044745137674732834
     &,         bs6=-0.0010745704660847233
     &,         cf0=2.1151080765239772e-13
     &,         cf1=-3.2260663894433345e-11
     &,         cf2=-3.329705958751961e-10
     &,         cf3=1.7648562021709124e-07
     &,         cf4=7.107636825694182e-06
     &,         cf5=-0.0013914681964973246
     &,         cf6=0.0406766967657759
     &,         p13=-1.296521881682694e-02
     &,         p12= 2.855780863283819e-01
     &,         p11=-1.597898515251717e+00
     &,         p10=-8.396975715683501e+00
     &,         p25= 3.790846746036765e-10
     &,         p24= 3.281964357650687e-09
     &,         p23= 1.962282433562894e-07
     &,         p22=-1.240239171056262e-06
     &,         p21=1.739759082358234e-07
     &,         p20=2.147264020369413e-05
     &,         p35=1.840430200185075e-07
     &,         p34=-2.793849676757154e-05
     &,         p33=1.735308193700643e-03
     &,         p32=-6.139315534216305e-02
     &,         p31=1.255457892775006e+00
     &,         p30=-1.663993561652530e+01
     &,         p40=4.579369142033410e-04

!     parameter (charnock=.014,ca=.4)!c ca is the von karman constant
!     parameter (alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041)
!     parameter (a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5)

!     real(kind=kind_phys) aa1,bb1,bb2,cc,cc1,cc2,arnu
!     parameter (aa1=-1.076,bb1=.7045,cc1=-.05808)
!     parameter (bb2=-.1954,cc2=.009999)
!     parameter (arnu=.135*rnu)
!
!    z0s_max=.196e-2 for u10_crit=25 m/s
!    z0s_max=.317e-2 for u10_crit=30 m/s
!    z0s_max=.479e-2 for u10_crit=35 m/s
!
! mbek -- toga-coare flux algorithm
!     parameter (rnu=1.51e-5,arnu=0.11*rnu)
!
!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals, wind is wind speed, 
!  surface roughness length is converted to m from cm
!
      do i=1,im
        if(flag_iter(i)) then 
          wind(i) = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))
     &                + max(0.0, min(ddvel(i), 30.0)), 1.0)
          tem1    = 1.0 + rvrdm1 * max(q1(i),1.e-8)
          thv1    = t1(i) * prslki(i) * tem1
          tvs     = 0.5 * (tsurf(i)+tskin(i)) * tem1
          qs1     = fpvs(t1(i))
          qs1     = max(1.0e-8, eps * qs1 / (prsl1(i) + epsm1 * qs1))

          z0      = 0.01 * z0rl(i)
          z0max   = max(1.0e-6, min(z0,z1(i)))
          z1i     = 1.0 / z1(i)

!  compute stability dependent exchange coefficients
!  this portion of the code is presently suppressed
!

          if(islimsk(i) == 0) then            ! over ocean
            ustar(i) = sqrt(grav * z0 / charnock)

!**  test xubin's new z0

!           ztmax  = z0max

            restar = max(ustar(i)*z0max*visi, 0.000001)

!           restar = log(restar)
!           restar = min(restar,5.)
!           restar = max(restar,-5.)
!           rat    = aa1 + (bb1 + cc1*restar) * restar
!           rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
!  rat taken from zeng, zhao and dickinson 1997

            rat    = min(7.0, 2.67 * sqrt(sqrt(restar)) - 2.57)
            ztmax  = z0max * exp(-rat)

          else                                ! over land and sea ice
!** xubin's new z0  over land and sea ice
            tem1 = 1.0 - shdmax(i)
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

            endif
            z0max = max(z0max,1.0e-6)
!
!           czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height dependance of czil
!            czilc = 0.8

            tem1 = 1.0 - sigmaf(i)
            ztmax = z0max*exp( - tem1*tem1
     &                         * czilc*ca*sqrt(ustar(i)*(0.01/1.5e-05)))

          endif       ! end of if(islimsk(i) == 0) then

          ztmax  = max(ztmax,1.0e-6)
          tem1   = z0max/z1(i)
          if (abs(1.0-tem1) > 1.0e-6) then
            ztmax1 = - beta*log(tem1)/(alpha2*(1.-tem1))
          else
            ztmax1 = 99.0
          endif
          if( z0max < 0.05 .and. snwdph(i) < 10.0 ) ztmax1 = 99.0


!  compute stability indices (rb and hlinf)

          dtv     = thv1 - tvs
          adtv    = max(abs(dtv),0.001)
          dtv     = sign(1.,dtv) * adtv
          rb(i)   = max(-5000.0, (grav+grav) * dtv * z1(i)
     &            / ((thv1 + tvs) * wind(i) * wind(i)))
          tem1    = 1.0 / z0max
          tem2    = 1.0 / ztmax
          fm(i)   = log((z0max+z1(i)) * tem1)
          fh(i)   = log((ztmax+z1(i)) * tem2)
          fm10(i) = log((z0max+10.)   * tem1)
          fh2(i)  = log((ztmax+2.)    * tem2)
          hlinf   = rb(i) * fm(i) * fm(i) / fh(i)
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
              fms    = fm(i) - pm
              fhs    = fh(i) - ph
              hl1    = fms * fms * rb(i) / fhs
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
            olinf = z1(i) / hlinf
            tem1  = 50.0 * z0max
            if(abs(olinf) <= tem1) then
              hlinf = -z1(i) / tem1
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
          fm(i)     = fm(i) - pm
          fh(i)     = fh(i) - ph
          fm10(i)   = fm10(i) - pm10
          fh2(i)    = fh2(i) - ph2
          cm(i)     = ca * ca / (fm(i) * fm(i))
          ch(i)     = ca * ca / (fm(i) * fh(i))
          tem1      = 0.00001/z1(i)
          cm(i) = max(cm(i), tem1)
          ch(i) = max(ch(i), tem1)
          stress(i) = cm(i) * wind(i) * wind(i)
          ustar(i)  = sqrt(stress(i))
!
!  update z0 over ocean
!
          if(islimsk(i) == 0) then

            z0 = (charnock / grav) * ustar(i) * ustar(i)

! mbek -- toga-coare flux algorithm
!           z0 = (charnock / grav) * ustar(i)*ustar(i) +  arnu/ustar(i)
!  new implementation of z0
!           cc = ustar(i) * z0 / rnu
!           pp = cc / (1. + cc)
!           ff = grav * arnu / (charnock * ustar(i) ** 3)
!           z0 = arnu / (ustar(i) * ff ** pp)

! -------------------------- modify z0 by kgao

! diagnose 10m wind (same as sfc_diag.f)
 
            u10m = u1(i) * fm10(i) / fm(i)
            v10m = v1(i) * fm10(i) / fm(i)
            ws10m = sqrt(u10m*u10m + v10m*v10m) 

! option - URI/GFDL (HWRF 2015)
! note there is discontinuity at 10m/s in original formulation
! needs to be fixed

            if (do_z0_hwrf15) then
              if (ws10m <= 5.0) then
                 z0 = 0.0185/9.8*(7.59e-4*ws10m**2+2.46e-2*ws10m)**2
              elseif (ws10m > 5.0 .and. ws10m <= 10.) then
                 z0 = 0.00000235*(ws10m**2-25.)+3.805129199617346e-05
              elseif (ws10m > 10.0 .and. ws10m <= 60.) then
                 z0 = bs6 + bs5*ws10m + bs4*ws10m**2 + bs3*ws10m**3 
     &              + bs2*ws10m**4 + bs1*ws10m**5 + bs0*ws10m**6
              else
                 z0 = cf6 + cf5*ws10m + cf4*ws10m**2 + cf3*ws10m**3
     &              + cf2*ws10m**4 + cf1*ws10m**5 + cf0*ws10m**6
              endif
            endif

! option - HWRF 2017

            if (do_z0_hwrf17) then
              if (ws10m <= 6.5) then
                z0 = exp( p10 + p11*ws10m + p12*ws10m**2 + p13*ws10m**3)
              elseif (ws10m > 6.5 .and. ws10m <= 15.7) then
                z0 = p25*ws10m**5 + p24*ws10m**4 + p23*ws10m**3 
     &             + p22*ws10m**2 + p21*ws10m + p20
              elseif (ws10m > 15.7 .and. ws10m <= 53.) then
                z0 = exp( p35*ws10m**5 + p34*ws10m**4 + p33*ws10m**3 
     &                  + p32*ws10m**2 + p31*ws10m + p30 )
              else
                z0 = p40
              endif
            endif

! option - GFS (low wind) + HWRF 2017 (high wind)

            if (do_z0_hwrf17_hwonly) then

              if (ws10m > wind_th_hwrf .and. ws10m <= 53.) then
                z0 = exp( p35*ws10m**5 + p34*ws10m**4 + p33*ws10m**3
     &                  + p32*ws10m**2 + p31*ws10m + p30 )
              elseif (ws10m > 53.) then
                z0 = p40
              endif

            endif

! option - GFS (low wind) + Moon et al (high wind)

            if (do_z0_moon) then
              wind_th_moon = 20. 
              a = 0.56
              b = -20.255
              c = wind_th_moon - 2.458
              ustar_th = (-b-sqrt(b*b-4*a*c))/(2*a)

              z0_adj = 0.001*(0.085*wind_th_moon - 0.58) - 
     &                 (charnock/grav)*ustar_th*ustar_th 

              ws10m = 2.458 + ustar(i)*(20.255-0.56*ustar(i))  ! Eq(7) Moon et al. 2007 
              if ( ws10m > wind_th_moon ) then                 ! No modification in low wind conditions 
                 z0 = 0.001*(0.085*ws10m - 0.58) - z0_adj      ! Eq(8b) Moon et al. 2007 modified by kgao
              endif
            endif

! ----------------------------  modify z0 end

            if (redrag) then
              z0rl(i) = 100.0 * max(min(z0, z0s_max), 1.e-7)
            else
              z0rl(i) = 100.0 * max(min(z0,.1), 1.e-7)
            endif
          endif
        endif                ! end of if(flagiter) loop
      enddo

      return
      end
