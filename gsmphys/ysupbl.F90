!-------------------------------------------------------------------------------
!
   subroutine ysupbl(ix,im,km,ndiff,ntcw,ntiw,                                 &
                  vtnp,utnp,ttnp,qtnp,                                         &
                  ux,vx,tx,qx,p2di,p2d,pi2d,                                   &
                  phi2di,psfcpa,                                               &
                  htrsw,htrlw,xmu,                                             &
                  z0rl,ust,hpbl,psim,psih,                                     & 
                  islmsk,heat,evap,wspd,br,                                    & 
                  dusfc,dvsfc,dtsfc,dqsfc,                                     &
                  dt,kpbl1d,u10,v10,                                           &
                  kinver,xkzm_m_in,xkzm_h_in,xkzm_s,xkzminv,                   &
                  dspheat,ent_fac,dkt,pfac_q,brcr_ub,rlam,afac,bfac)
!-------------------------------------------------------------------------------
   use machine, only : kind_phys
!   use mpp_mod, only: mpp_pe
   use physcons, cp=>con_cp, g=>con_g, rovcp=>con_rocp, rd=>con_rd, rv=>con_rv,&
                 rovg=>con_rog, ep1=>con_fvirt, ep2=>con_eps, xlv=>con_hvap
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!     this code is a revised vertical diffusion package ("ysupbl")
!     with a nonlocal turbulent mixing in the pbl after "mrfpbl".
!     the ysupbl (hong et al. 2006) is based on the study of noh
!     et al.(2003) and accumulated realism of the behavior of the
!     troen and mahrt (1986) concept implemented by hong and pan(1996).
!     the major ingredient of the ysupbl is the inclusion of an explicit
!     treatment of the entrainment processes at the entrainment layer.
!     this routine uses an implicit approach for vertical flux
!     divergence and does not require "miter" timesteps.
!     it includes vertical diffusion in the stable atmosphere
!     and moist vertical diffusion in clouds.
!
!     mrfpbl:
!     coded by song-you hong (ncep), implemented by jimy dudhia (ncar)
!              fall 1996
!
!     ysupbl:
!     coded by song-you hong (yonsei university) and implemented by
!              song-you hong (yonsei university) and jimy dudhia (ncar)
!              summer 2002
!
!     further modifications :
!              an enhanced stable layer mixing, april 2008
!              ==> increase pbl height when sfc is stable (hong 2010)
!              pressure-level diffusion, april 2009
!              ==> negligible differences
!              implicit forcing for momentum with clean up, july 2009
!              ==> prevents model blowup when sfc layer is too low
!              incresea of lamda, maximum (30, 0.1 x del z) feb 2010
!              ==> prevents model blowup when delz is extremely large
!              revised prandtl number at surface, peggy lemone, feb 2010
!              ==> increase kh, decrease mixing due to counter-gradient term
!              revised thermal, shin et al. mon. wea. rev. , songyou hong, aug 2011
!              ==> reduce the thermal strength when z1 < 0.1 h
!              revised prandtl number for free convection, dudhia, mar 2012
!              ==> pr0 = 1 + bke (=0.272) when newtral, kh is reduced
!              minimum kzo = 0.01, lo = min (30m,delz), hong, mar 2012
!              ==> weaker mixing when stable, and les resolution in vertical
!              gz1oz0 is removed, and phim phih are ln(z1/z0)-phim,h, hong, mar 2012
!              ==> consider thermal z0 when differs from mechanical z0
!              a bug fix in wscale computation in stable bl, sukanta basu, jun 2012
!              ==> wscale becomes small with height, and less mixing in stable bl
!              revision in background diffusion (kzo), jan 2016
!              ==> kzo = 0.1 for momentum and = 0.01 for mass to account for
!                  internal wave mixing of large et al. (1994), songyou hong, feb 2016
!              ==> alleviate superious excessive mixing when delz is large
!              Added minimum diffusion (a bound, not additional background) and
!              dissipative heating from the GFS EDMF scheme. Hailey Shin, feb 2018.
!
!     references:
!
!        hong (2010) quart. j. roy. met. soc
!        hong, noh, and dudhia (2006), mon. wea. rev.
!        hong and pan (1996), mon. wea. rev.
!        noh, chun, hong, and raasch (2003), boundary layer met.
!        troen and mahrt (1986), boundary layer met.
!
!-------------------------------------------------------------------------------
!
   integer,parameter :: nqci   = 2
   integer,parameter :: imvdif = 1
   integer,parameter :: ysu_topdown_pblmix  = 1
   real(kind=kind_phys),parameter :: karman = 0.4
   real(kind=kind_phys),parameter :: xkzminm = 0.1,xkzminh = 0.01
   real(kind=kind_phys),parameter :: xkzmin = 0.01,xkzmax = 1000.
   real(kind=kind_phys),parameter :: rimin = -100.
!   real(kind=kind_phys),parameter :: rlam = 30.,prmin = 0.25,prmax = 4.
   real(kind=kind_phys),parameter :: prmin = 0.25,prmax = 4.
!   real(kind=kind_phys),parameter :: brcr_ub = 0.0,brcr_sb = 0.25
   real(kind=kind_phys),parameter :: brcr_sb = 0.25
   real(kind=kind_phys),parameter :: cori = 1.e-4
!   real(kind=kind_phys),parameter :: afac = 6.8,bfac = 6.8
   real(kind=kind_phys),parameter :: pfac = 2.0 !,pfac_q = 2.0
   real(kind=kind_phys),parameter :: phifac = 8.,sfcfrac = 0.1
   real(kind=kind_phys),parameter :: d1 = 0.02, d2 = 0.05, d3 = 0.001
   real(kind=kind_phys),parameter :: h1 = 0.33333335, h2 = 0.6666667
   real(kind=kind_phys),parameter :: zfmin = 1.e-8
   real(kind=kind_phys),parameter :: aphi5 = 5.,aphi16 = 16.
   real(kind=kind_phys),parameter :: tmin=1.e-2
   real(kind=kind_phys),parameter :: gamcrt = 3.,gamcrq = 2.e-3
   real(kind=kind_phys),parameter :: xka = 2.4e-5
   real(kind=kind_phys),parameter :: rcl=1.0
   real(kind=kind_phys),parameter :: dw2min=0.0001
!
   integer,intent(in   ) :: ix,im,km,ndiff,ntcw,ntiw
!
   real(kind=kind_phys),intent(in   ) :: dt
   real(kind=kind_phys),intent(in   ) :: xkzm_m_in,xkzm_h_in,xkzm_s,xkzminv,ent_fac,pfac_q
   real(kind=kind_phys),intent(in   ) :: brcr_ub,rlam,afac,bfac
!
   real(kind=kind_phys),dimension(   1:ix ,   1:km  )                 , & !! Statein%ugrs (ix,km)
                        intent(in   ) ::                            ux, & !! 
                                                                    vx, & !!
                                                                    tx, & !! Radtend%htrsw (ix,km)
                                                                 htrsw, & !!
                                                                 htrlw
!
   real(kind=kind_phys),dimension(   1:ix ,   1:km *ndiff )           , & !! 
                        intent(in   ) ::                            qx
!
   real(kind=kind_phys),dimension(   1:im ,   1:km  )                 , & !! dvdt (im,km)
                        intent(inout) ::                          vtnp, & !!
                                                                  utnp, & !!
                                                                  ttnp
   real(kind=kind_phys),dimension(   1:im ,   1:km *ndiff )           , &
                        intent(inout) ::                          qtnp
!
   real(kind=kind_phys),dimension(   1:ix ,   1:km +1 )               , & !! Statein%prsi
                        intent(in   ) ::                          p2di, &
                                                                phi2di
!
   real(kind=kind_phys),dimension(   1:ix ,   1:km  )                 , & !! Statein%prsl
                        intent(in   ) ::                           p2d, & !! Statein%prslk
                                                                  pi2d
!
   real(kind=kind_phys),dimension(   1:im  )                          , & !! Sfcprop%uustar
                        intent(in   ) ::                           ust, & !! Sfcprop%zorl
                                                                  z0rl, &
                                                                   xmu
   real(kind=kind_phys),dimension(   1:im  )                          , &
                        intent(in   ) ::                          heat, &
                                                                  evap
   real(kind=kind_phys),dimension(   1:im  )                          , & !! rb
                        intent(in   ) ::                            br, & !! Sfcprop%ffmm
                                                                  psim, & !! Sfcprop%ffhh
                                                                  psih, & !! wind
                                                                  wspd, &
                                                                psfcpa
!
   real(kind=kind_phys),dimension(   1:im  )                          , & !! Diag%hpbl
                        intent(inout) ::                          hpbl, & !! dusfc1, dvsfc1
                                                           dusfc,dvsfc, & !! dtsfc1, dqsfc1
                                                           dtsfc,dqsfc
!
   integer,dimension(   1:im  ),intent(in   ) ::                islmsk, kinver ! kinver = levs for most purposes
   integer,dimension(   1:im  ),intent(out  ) ::                kpbl1d
   logical,intent(in   ) ::                                    dspheat
!
! local vars
!
   real(kind=kind_phys),dimension(   1:im  ) ::                  xland, &
                                                                   hfx, &
                                                                   qfx
!
   real(kind=kind_phys),dimension(   1:im  ) ::                    hol, &
                                                                   znt
   real(kind=kind_phys),dimension(   1:im ,   1:km+1  ) ::          zq, &
                                                               p2diorg
!
   real(kind=kind_phys),dimension(   1:im ,   1:km  ) ::                &
                                                        thx,thvx,thlix, &
                                                                   del, &
                                                                   dza, &
                                                                   dzq, &
                                                                 xkzom, &
                                                                 xkzoh, &
                                                                    za
   real(kind=kind_phys),dimension(   1:im ,   1:km  ) ::      rthraten 
!
   real(kind=kind_phys),dimension(   1:im  ) ::                         &
                                                                  rhox, &
                                                                govrth, &
                                                           zl1,thermal, &
                                                                wscale, &
                                                           hgamt,hgamq, &
                                                             brdn,brup, &
                                                             phim,phih, &
                                                                 prpbl, &
                                                       wspd1,thermalli
!
   real(kind=kind_phys),dimension(   1:im ,   1:km  ) ::     xkzm,xkzh, &
                                                                 f1,f2, &
                                                                 r1,r2, &
                                                                 ad,au, &
                                                                    cu, &
                                                                    al, &
                                                                  xkzq, &
                                                                  zfac, &
                                                                 rhox2, &
                                                                hgamt2
!
   real(kind=kind_phys),dimension(   1:im  )                          , & !! Diag%u10m
                        intent(in   ) ::                           u10, & !! Diag%v10m
                                                                   v10
   real(kind=kind_phys),dimension(   1:im  ) ::                         &
                                                                  brcr, &
                                                                 sflux, &
                                                                  zol1, &
                                                             brcr_sbro
!
   real(kind=kind_phys),dimension(   1:im ,   1:km , ndiff) ::   r3,f3
   real(kind=kind_phys),dimension(   1:ix ,   1:km , nqci ) ::    qxci
   integer, dimension(   1:im  ) ::                       kpbl,kpblold
   integer, dimension(   1:im  ) ::                                kx1
!
   logical, dimension(   1:im  ) ::                             pblflg, &
                                                                sfcflg, & ! br <= 0 
                                                                stable, &
                                                              cloudflg

   real(kind=kind_phys),dimension(   1:im, 1:km-1), intent(OUT), OPTIONAL :: dkt
! Local:
   real(kind=kind_phys),dimension(   1:im ,   1:km  ) :: diss
! SJL
   real(kind=kind_phys),dimension(   1:ix ,   1:km  ) :: p2m

   logical :: definebrup
!
   integer :: n,i,k,l,ic,is,kk
   integer :: klpbl, ktrace1, ktrace2, ktrace3
   integer :: its,ite,kts,kte
!
!
   real(kind=kind_phys) :: dt2,rdt,spdk2,fm,fh,hol1
   real(kind=kind_phys) :: gamfac,vpert,prnum,prnum0
   real(kind=kind_phys) :: ss,ri,qmean,tmean,alph,chi,zk,rl2,dk,sri
   real(kind=kind_phys) :: brint,dtodsd,dtodsu,rdz
   real(kind=kind_phys) :: dsdzt,dsdzq,dsdz2,rlamdz
   real(kind=kind_phys) :: utend,vtend,ttend,qtend
   real(kind=kind_phys) :: dtstep,govrthv
   real(kind=kind_phys) :: cont, conq, conw, conwrc
!
   real(kind=kind_phys),dimension(   1:im ,   1:km  ) ::       wscalek, &
                                                              wscalek2
   real(kind=kind_phys),dimension(   1:im  ) ::                  wstar, &
                                                                 delta
   real(kind=kind_phys),dimension(   1:im ,   1:km  ) ::   xkzml,xkzhl, &
                                                        zfacent,entfac
   real(kind=kind_phys),dimension(   1:im  ) ::                   ust3, &
                                                                wstar3, &
                                                              wstar3_2, &
                                                           hgamu,hgamv, &
                                                               wm2, we, &
                                                                bfxpbl, &
                                                         hfxpbl,qfxpbl, &
                                                         ufxpbl,vfxpbl, &
                                                                 dthvx
!
   real(kind=kind_phys),dimension(   1:im  ) ::                uox,vox
   real(kind=kind_phys) :: ptem
   real(kind=kind_phys) :: xkzm_m,xkzm_h,xkzm_fac
   real(kind=kind_phys) :: dw2,shr2,ti,bf,tem,tem2
!
   real(kind=kind_phys) :: prnumfac,bfx0,hfx0,qfx0,delb,dux,dvx,        &
            dsdzu,dsdzv,wm3,dthx,dqx,wspd10,ross,tem1,dsig,tvcon,conpr, &
            prfac,prfac2,phim8z,radsum,tmp1,templ,rvls,temps,ent_eff,   &
            rcldb,bruptmp,radflux
!
!-------------------------------------------------------------------------------
!
   its = 1
   ite = im
   kts = 1
   kte = km
   klpbl = kte
!
   cont=cp/g
   conq=xlv/g
   conw=1./g
   conwrc = conw*sqrt(rcl)
   conpr = bfac*karman*sfcfrac
!
   xkzm_fac = 0.5
   xkzm_m   = xkzm_m_in*xkzm_fac
   xkzm_h   = xkzm_h_in*xkzm_fac
!
!  k-start index for tracer diffusion
!
   ktrace1 = 0
   ktrace2 = 0 + kte*(ntcw-1)
   ktrace3 = 0 + kte*(ntiw-1)
!
   do k = kts,kte
     do i = its,ite
       rthraten(i,k) = (htrsw(i,k)*xmu(i)+htrlw(i,k))/pi2d(i,k) !converts temp/s to theta/s
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       qxci(i,k,1:nqci) = 0.0
       if(ntcw > 0) qxci(i,k,1) = qx(i,ktrace2+k)
       if(ntiw > 0) qxci(i,k,2) = qx(i,ktrace3+k)
     enddo
   enddo
!
   do i = its,ite
     znt(i) = 0.01*z0rl(i)
   enddo
!
   do k = kts,kte
     do i = its,ite
       thx(i,k) = tx(i,k)/pi2d(i,k)
!      thlix(i,k) = (tx(i,k)-xlv*qxci(i,k,1)/cp-2.834E6*qxci(i,k,2)/cp)/pi2d(i,k)
       thlix(i,k) = (tx(i,k)-(xlv*qxci(i,k,1)+2.834E6*qxci(i,k,2))/cp) / pi2d(i,k)
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       tvcon = (1.+ep1*qx(i,k))
       thvx(i,k) = thx(i,k)*tvcon
     enddo
   enddo
!
   do i = its,ite
     tvcon = (1.+ep1*qx(i,1))
     rhox(i) = psfcpa(i)/(rd*tx(i,1)*tvcon)
     govrth(i) = g/thx(i,1)
   enddo
!
!-----compute the height of full- and half-sigma levels above ground
!     level, and the layer thicknesses.
!
   do k = kts,kte
     do i = its,ite
       zq(i,k) = phi2di(i,k)/g
       if(k.eq.kte) zq(i,k+1) =  phi2di(i,k+1)/g
       tvcon = (1.+ep1*qx(i,k))
       rhox2(i,k) = p2d(i,k)/(rd*tx(i,k)*tvcon)
! SJL
       p2m(i,k) = 0.5*(p2di(i,k) + p2di(i,k+1))
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       za(i,k) = 0.5*(zq(i,k)+zq(i,k+1))
       dzq(i,k) = zq(i,k+1)-zq(i,k)
       del(i,k) = p2di(i,k)-p2di(i,k+1)
       p2diorg(i,k) = p2di(i,k)
       if(k.eq.kte) p2diorg(i,k+1) = p2di(i,k+1)
     enddo
   enddo
!
   do i = its,ite
     dza(i,1) = za(i,1)
   enddo
!
   do k = kts+1,kte
     do i = its,ite
       dza(i,k) = za(i,k)-za(i,k-1)
     enddo
   enddo
!
   do i = its,ite
      uox(i) = 0.0
      vox(i) = 0.0
      xland(i) = 1
      if(islmsk(i).eq.0) xland(i) = 2
      hfx(i) = heat(i)*rhox(i)*cp
      qfx(i) = evap(i)*rhox(i)
   enddo

!
!
!-----initialize vertical tendencies
!
   do i = its,ite
!    wspd1(i) = sqrt( (ux(i,1)-uox(i))*(ux(i,1)-uox(i)) + (vx(i,1)-vox(i))*(vx(i,1)-vox(i)) )+1.e-9
     wspd1(i) = max(1.e-9, sqrt((ux(i,1)-uox(i))**2 + (vx(i,1)-vox(i))**2))
   enddo
!
!---- compute vertical diffusion
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     compute preliminary variables
!
   dtstep = dt
!  dt2 = 2.*dtstep
   dt2 = dtstep
   rdt = 1./dt2
!
   do i = its,ite
     bfxpbl(i) = 0.0
     hfxpbl(i) = 0.0
     qfxpbl(i) = 0.0
     ufxpbl(i) = 0.0
     vfxpbl(i) = 0.0
     hgamu(i)  = 0.0
     hgamv(i)  = 0.0
     delta(i)  = 0.0
     wstar3_2(i) =  0.0
   enddo
!
   do k = kts,klpbl
     do i = its,ite
       wscalek(i,k) = 0.0
       wscalek2(i,k) = 0.0
     enddo
   enddo
!
   do k = kts,klpbl
     do i = its,ite
       zfac(i,k) = 0.0
     enddo
   enddo
!
!<--- background vertical diffusivity
!
!   do k = kts,klpbl-1
!     do i = its,ite
!       xkzom(i,k) = xkzminm
!       xkzoh(i,k) = xkzminh
!     enddo
!   enddo
!
   do i = its,ite
     kx1(i) = 1
   enddo
!
   do k = kts,kte
     do i = its,ite
       xkzoh(i,k) = 0.0
       xkzom(i,k) = 0.0
     enddo
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       if (k.lt.kinver(i)) then
         ptem       = p2di(i,k+1)/p2di(i,kts)
         tem1       = 1.0 - ptem
         tem1       = tem1 * tem1 * 10.0
         xkzoh(i,k) = xkzm_h * min(1.0, exp(-tem1))
         if (ptem.ge.xkzm_s) then
            !constant xkzm_m below xkzm_s (default 1.0)
           xkzom(i,k) = xkzm_m
           kx1(i)     = k + 1
         else
           if (k.eq.kx1(i) .and. k.gt.kts) then
             tem1 = 1.0 - p2di(i,k+1)/p2di(i,k)
           else
             tem1 = 1.0 - p2di(i,k+1)/p2di(i,kts)
           endif
           tem1 = tem1 * tem1 * 5.0
           !xkzm_m * exp(-5.0*(1.-dp_rat)**2) < xkzm_m
           xkzom(i,k) = xkzm_m * min(1.0, exp(-tem1))
         endif
       endif
     enddo
   enddo
!
   do k = kts,kte/2
     do i = its,ite
       if(zq(i,k+1) .gt. 250.) then
         tem1 = (tx(i,k+1)-tx(i,k))/(za(i,k+1)-za(i,k))
         if(tem1 .gt. 1.e-5) then
           xkzoh(i,k) = min(xkzoh(i,k),xkzminv)
         endif
       endif
     enddo
   enddo
!
   do i = its,ite
     dusfc(i) = 0.
     dvsfc(i) = 0.
     dtsfc(i) = 0.
     dqsfc(i) = 0.
   enddo
!
   do i = its,ite
     hgamt(i)  = 0.
     hgamq(i)  = 0.
     wscale(i) = 0.
     kpbl(i)   = 1
     hpbl(i)   = zq(i,1)
     zl1(i)    = za(i,1)
     thermal(i)= thvx(i,1)
     thermalli(i) = thlix(i,1)
     pblflg(i) = .true.
     sfcflg(i) = .true.
!    sflux(i) = hfx(i)/rhox(i)/cp + qfx(i)/rhox(i)*ep1*thx(i,1)
     sflux(i) = hfx(i)/(rhox(i)*cp) + qfx(i)/rhox(i)*ep1*thx(i,1)
     if(br(i).gt.0.0) sfcflg(i) = .false.
   enddo
!
!     compute the first guess of pbl height
!
   do i = its,ite
     stable(i) = .false.
     brup(i) = br(i)
     brcr(i) = brcr_ub
   enddo
!
   do k = 2,klpbl
     do i = its,ite
       if(.not.stable(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
!        brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         brup(i) = (thvx(i,k)-thermal(i))*g*za(i,k)/(thvx(i,1)*spdk2)
         kpbl(i) = k
         stable(i) = brup(i).gt.brcr(i)
       endif
     enddo
   enddo
!
   do i = its,ite
     k = kpbl(i)
     if(brdn(i).ge.brcr(i))then
       brint = 0.
     elseif(brup(i).le.brcr(i))then
       brint = 1.
     else
       brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
     endif
     hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
     if(hpbl(i).lt.zq(i,2)) kpbl(i) = 1
     if(kpbl(i).le.1) pblflg(i) = .false.
   enddo
!
   do i = its,ite
     fm = psim(i)
     fh = psih(i)
     zol1(i) = max(br(i)*fm*fm/fh,rimin)
     if(sfcflg(i))then
       zol1(i) = min(zol1(i),-zfmin)
     else
       zol1(i) = max(zol1(i),zfmin)
     endif
     hol1 = zol1(i)*hpbl(i)/zl1(i)*sfcfrac
     if(sfcflg(i))then
       phim(i) = (1.-aphi16*hol1)**(-1./4.)
       phih(i) = (1.-aphi16*hol1)**(-1./2.)
       bfx0 = max(sflux(i),0.)
!      hfx0 = max(hfx(i)/rhox(i)/cp,0.)
       hfx0 = max(hfx(i)/(rhox(i)*cp), 0.) ! not used
       qfx0 = max(ep1*thx(i,1)*qfx(i)/rhox(i),0.) ! not used
       wstar3(i) = (govrth(i)*bfx0*hpbl(i))
       wstar(i) = (wstar3(i))**h1
     else
       phim(i) = (1.+aphi5*hol1)
       phih(i) = phim(i)
       wstar(i)  = 0.
       wstar3(i) = 0.
     endif
     ust3(i)   = ust(i)**3.
     wscale(i) = (ust3(i)+phifac*karman*wstar3(i)*0.5)**h1
     wscale(i) = min(wscale(i),ust(i)*aphi16)
     wscale(i) = max(wscale(i),ust(i)/aphi5)
   enddo
!
!     compute the surface variables for pbl height estimation
!     under unstable conditions
!
   do i = its,ite
     if(sfcflg(i).and.sflux(i).gt.0.0)then
!      gamfac   = bfac/rhox(i)/wscale(i)
       gamfac   = bfac / (rhox(i)*wscale(i))
       hgamt(i) = min(gamfac*hfx(i)/cp,gamcrt)
       hgamq(i) = min(gamfac*qfx(i),gamcrq)
       vpert = (hgamt(i)+ep1*thx(i,1)*hgamq(i))/bfac*afac
       thermal(i) = thermal(i)+max(vpert,0.)*min(za(i,1)/(sfcfrac*hpbl(i)),1.0)
       thermalli(i)= thermalli(i)+max(vpert,0.)*min(za(i,1)/(sfcfrac*hpbl(i)),1.0)
       hgamt(i) = max(hgamt(i),0.0)
       hgamq(i) = max(hgamq(i),0.0)
!      brint    = -15.9*ust(i)*ust(i)/wspd(i)*wstar3(i)/(wscale(i)**4.)
       brint    = -15.9*ust(i)*ust(i)*wstar3(i)/(wspd(i)*wscale(i)**4.)
       hgamu(i) = brint*ux(i,1)
       hgamv(i) = brint*vx(i,1)
     else
       pblflg(i) = .false.
     endif
   enddo
!
!     enhance the pbl height by considering the thermal
!
   do i = its,ite
     if(pblflg(i))then
       kpbl(i) = 1
       hpbl(i) = zq(i,1)
     endif
   enddo
!
   do i = its,ite
     if(pblflg(i))then
       stable(i) = .false.
       brup(i) = br(i)
       brcr(i) = brcr_ub
     endif
   enddo
!
   do k = 2,klpbl
     do i = its,ite
       if(.not.stable(i).and.pblflg(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
!        brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         brup(i) = (thvx(i,k)-thermal(i))*g*za(i,k)/(thvx(i,1)*spdk2)
         kpbl(i) = k
         stable(i) = brup(i).gt.brcr(i)
       endif
     enddo
   enddo
!
!     enhance pbl by theta-li
!
   if (ysu_topdown_pblmix.eq.1)then
     do i = its,ite
        kpblold(i) = kpbl(i)
        definebrup=.false.
        do k = kpblold(i), kte-1
           spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
!          bruptmp = (thlix(i,k)-thermalli(i))*(g*za(i,k)/thlix(i,1))/spdk2
           bruptmp = (thlix(i,k)-thermalli(i))*g*za(i,k)/(thlix(i,1)*spdk2)
           stable(i) = bruptmp.ge.brcr(i)
           if (definebrup) then
           kpbl(i) = k
           brup(i) = bruptmp
           definebrup=.false.
           endif
           if (.not.stable(i)) then !overwrite brup brdn values
           brdn(i)=bruptmp
           definebrup=.true.
           pblflg(i)=.true.
           endif
        enddo
     enddo
   endif

   do i = its,ite
     if(pblflg(i)) then
       k = kpbl(i)
       if(brdn(i).ge.brcr(i))then
         brint = 0.
       elseif(brup(i).le.brcr(i))then
         brint = 1.
       else
         brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
       endif
       hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
       if(hpbl(i).lt.zq(i,2)) kpbl(i) = 1
       if(kpbl(i).le.1) pblflg(i) = .false.
     endif
   enddo
!
!     stable boundary layer
!
   do i = its,ite
     if((.not.sfcflg(i)).and.hpbl(i).lt.zq(i,2)) then
       brup(i) = br(i)
       stable(i) = .false.
     else
       stable(i) = .true.
     endif
   enddo
!
   do i = its,ite
     if((.not.stable(i)).and.((xland(i)-1.5).ge.0))then
       wspd10 = u10(i)*u10(i) + v10(i)*v10(i)
       wspd10 = sqrt(wspd10)
       ross = wspd10 / (cori*znt(i))
       brcr_sbro(i) = min(0.16*(1.e-7*ross)**(-0.18),.3)
     endif
   enddo
!
   do i = its,ite
     if(.not.stable(i))then
       if((xland(i)-1.5).ge.0)then
         brcr(i) = brcr_sbro(i)
       else
         brcr(i) = brcr_sb
       endif
     endif
   enddo
!
   do k = 2,klpbl
     do i = its,ite
       if(.not.stable(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
!        brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         brup(i) = (thvx(i,k)-thermal(i))*g*za(i,k)/(thvx(i,1)*spdk2)
         kpbl(i) = k
         stable(i) = brup(i).gt.brcr(i)
       endif
     enddo
   enddo
!
   do i = its,ite
     if((.not.sfcflg(i)).and.hpbl(i).lt.zq(i,2)) then
       k = kpbl(i)
       if(brdn(i).ge.brcr(i))then
         brint = 0.
       elseif(brup(i).le.brcr(i))then
         brint = 1.
       else
         brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
       endif
       hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
       if(hpbl(i).lt.zq(i,2)) kpbl(i) = 1
       if(kpbl(i).le.1) pblflg(i) = .false.
     endif
   enddo
!
!     estimate the entrainment parameters
!
   do i = its,ite
     cloudflg(i)=.false. 
     if(pblflg(i)) then
       k = kpbl(i) - 1
       wm3       = wstar3(i) + 5. * ust3(i)
       wm2(i)    = wm3**h2
!       bfxpbl(i) = -ent_fac*thvx(i,1)/g*wm3/hpbl(i)
!       bfxpbl(i) = -0.15*thvx(i,1)/g*wm3/hpbl(i)
       bfxpbl(i) = -ent_fac*thvx(i,1)*wm3/(g*hpbl(i))
       dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),tmin)
       we(i) = max(bfxpbl(i)/dthvx(i),-sqrt(wm2(i)))
       if((qxci(i,k,1)+qxci(i,k,2)).gt.0.01e-3.and.ysu_topdown_pblmix.eq.1)then
           if ( kpbl(i) .ge. 2) then
                cloudflg(i)=.true. 
                templ=thlix(i,k)*(p2di(i,k+1)/100000)**rovcp
                !rvls is ws at full level
                rvls=100.*6.112*EXP(17.67*(templ-273.16)/(templ-29.65))*(ep2/p2di(i,k+1))
                temps=templ + ((qx(i,k)+qxci(i,k,1))-rvls)/(cp/xlv  + &
                ep2*xlv*rvls/(rd*templ**2))
                rvls=100.*6.112*EXP(17.67*(temps-273.15)/(temps-29.65))*(ep2/p2di(i,k+1))
                rcldb=max((qx(i,k)+qxci(i,k,1))-rvls,0.)
                !entrainment efficiency
                dthvx(i)  = (thlix(i,k+2)+thx(i,k+2)*ep1*(qx(i,k+2)+qxci(i,k+2,1))) &
                          - (thlix(i,k) + thx(i,k)  *ep1*(qx(i,k)  +qxci(i,k  ,1)))
                dthvx(i)  = max(dthvx(i),0.1)
                tmp1      = xlv/cp * rcldb/(pi2d(i,k)*dthvx(i))
                ent_eff   = 0.2 * 8. * tmp1 +0.2

                radsum=0.
                do kk = 1,kpbl(i)-1
                   radflux=rthraten(i,kk)*pi2d(i,kk) !converts theta/s to temp/s
                   radflux=radflux*cp/g*(p2diORG(i,kk)-p2diORG(i,kk+1)) ! converts temp/s to W/m^2
                   if (radflux < 0.0 ) radsum=abs(radflux)+radsum
                enddo
                radsum=max(radsum,0.0)

                !recompute entrainment from sfc thermals
                bfx0 = max(max(sflux(i),0.0)-radsum/rhox2(i,k)/cp,0.)
                bfx0 = max(sflux(i),0.0)
                wm3 = (govrth(i)*bfx0*hpbl(i))+5. * ust3(i)
                wm2(i)    = wm3**h2
!                bfxpbl(i) = -ent_fac*thvx(i,1)/g*wm3/hpbl(i)
!                bfxpbl(i) = -0.15*thvx(i,1)/g*wm3/hpbl(i)
                bfxpbl(i) = -ent_fac*thvx(i,1)*wm3/(g*hpbl(i))
                dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),tmin)
                we(i) = max(bfxpbl(i)/dthvx(i),-sqrt(wm2(i)))

                !entrainment from PBL top thermals
                bfx0 = max(radsum/rhox2(i,k)/cp-max(sflux(i),0.0),0.)
                bfx0 = max(radsum/rhox2(i,k)/cp,0.)
                wm3       = (g/thvx(i,k)*bfx0*hpbl(i)) ! this is wstar3(i)
                wm2(i)    = wm2(i)+wm3**h2
                bfxpbl(i) = - ent_eff * bfx0
                dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),0.1)
                we(i) = we(i) + max(bfxpbl(i)/dthvx(i),-sqrt(wm3**h2))

                !wstar3_2
                bfx0 = max(radsum/rhox2(i,k)/cp,0.)
                wstar3_2(i) =  (g/thvx(i,k)*bfx0*hpbl(i))
                !recompute hgamt 
                wscale(i) = (ust3(i)+phifac*karman*(wstar3(i)+wstar3_2(i))*0.5)**h1
                wscale(i) = min(wscale(i),ust(i)*aphi16)
                wscale(i) = max(wscale(i),ust(i)/aphi5)
                gamfac   = bfac/rhox(i)/wscale(i)
                hgamt(i) = min(gamfac*hfx(i)/cp,gamcrt)
                hgamq(i) = min(gamfac*qfx(i),gamcrq)
                gamfac   = bfac/rhox2(i,k)/wscale(i)
                hgamt2(i,k) = min(gamfac*radsum/cp,gamcrt)
                hgamt(i) = max(hgamt(i),0.0) + max(hgamt2(i,k),0.0)
                brint    = -15.9*ust(i)*ust(i)/wspd(i)*(wstar3(i)+wstar3_2(i))/(wscale(i)**4.)
                hgamu(i) = brint*ux(i,1)
                hgamv(i) = brint*vx(i,1)
           endif
       endif
       prpbl(i) = 1.0
       dthx  = max(thx(i,k+1)-thx(i,k),tmin)
       dqx   = min(qx(i,k+1)-qx(i,k),0.0)
       hfxpbl(i) = we(i)*dthx
       qfxpbl(i) = we(i)*dqx
!
       dux = ux(i,k+1)-ux(i,k)
       dvx = vx(i,k+1)-vx(i,k)
       if(dux.gt.tmin) then
         ufxpbl(i) = max(prpbl(i)*we(i)*dux,-ust(i)*ust(i))
       elseif(dux.lt.-tmin) then
         ufxpbl(i) = min(prpbl(i)*we(i)*dux,ust(i)*ust(i))
       else
         ufxpbl(i) = 0.0
       endif
       if(dvx.gt.tmin) then
         vfxpbl(i) = max(prpbl(i)*we(i)*dvx,-ust(i)*ust(i))
       elseif(dvx.lt.-tmin) then
         vfxpbl(i) = min(prpbl(i)*we(i)*dvx,ust(i)*ust(i))
       else
         vfxpbl(i) = 0.0
       endif
       delb  = govrth(i)*d3*hpbl(i)
       delta(i) = min(d1*hpbl(i) + d2*wm2(i)/delb,100.)
     endif
   enddo
!
   do k = kts,klpbl
     do i = its,ite
       if(pblflg(i).and.k.ge.kpbl(i))then
         entfac(i,k) = ((zq(i,k+1)-hpbl(i))/delta(i))**2.
       else
         entfac(i,k) = 1.e30
       endif
     enddo
   enddo
!
!     compute diffusion coefficients below pbl
!
   do k = kts,klpbl
     do i = its,ite
       if(k.lt.kpbl(i)) then
         zfac(i,k) = min(max((1.-(zq(i,k+1)-zl1(i))/(hpbl(i)-zl1(i))),zfmin),1.)
         zfacent(i,k) = (1.-zfac(i,k))**3.
         wscalek(i,k) = (ust3(i)+phifac*karman*wstar3(i)*(1.-zfac(i,k)))**h1
         wscalek2(i,k) = (phifac*karman*wstar3_2(i)*(zfac(i,k)))**h1
         if(sfcflg(i)) then
           prfac = conpr
           prfac2 = 15.9*(wstar3(i)+wstar3_2(i))/ust3(i)/(1.+4.*karman*(wstar3(i)+wstar3_2(i))/ust3(i))
           prnumfac = -3.*(max(zq(i,k+1)-sfcfrac*hpbl(i),0.))**2./hpbl(i)**2.
         else
           prfac = 0.
           prfac2 = 0.
           prnumfac = 0.
           phim8z = 1.+aphi5*zol1(i)*zq(i,k+1)/zl1(i)
           wscalek(i,k) = ust(i)/phim8z
           wscalek(i,k) = max(wscalek(i,k),0.001)
         endif
         prnum0 = (phih(i)/phim(i)+prfac)
         prnum0 = max(min(prnum0,prmax),prmin)
           xkzm(i,k) = wscalek(i,k) *karman*    zq(i,k+1)      *    zfac(i,k)**pfac+ &
                       wscalek2(i,k)*karman*(hpbl(i)-zq(i,k+1))*(1-zfac(i,k))**pfac
         !Do not include xkzm at kpbl-1 since it changes entrainment
         if (k.eq.kpbl(i)-1.and.cloudflg(i).and.we(i).lt.0.0) then
           xkzm(i,k) = 0.0
         endif
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzq(i,k) = xkzm(i,k)/prnum*zfac(i,k)**(pfac_q-pfac)
         prnum0 = prnum0/(1.+prfac2*karman*sfcfrac)
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzh(i,k) = xkzm(i,k)/prnum
         xkzm(i,k) = max(xkzm(i,k),xkzom(i,k))
         xkzh(i,k) = max(xkzh(i,k),xkzoh(i,k))
         xkzq(i,k) = max(xkzq(i,k),xkzoh(i,k))
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         xkzq(i,k) = min(xkzq(i,k),xkzmax)
       endif
     enddo
   enddo
!
!     compute diffusion coefficients over pbl (free atmosphere)
!
!   if (lprnt) then
!      i = im/2
!      write(mpp_pe()+1000,*) kpbl(i), hpbl(i), pblflg(i), sfcflg(i)
!   endif
   do k = kts,kte-1
     do i = its,ite
        !prnum = 0. ! DEBUG
       if(k.ge.kpbl(i)) then
!        ss = ((ux(i,k+1)-ux(i,k))*(ux(i,k+1)-ux(i,k))                         &
!             +(vx(i,k+1)-vx(i,k))*(vx(i,k+1)-vx(i,k)))                        &
!             /(dza(i,k+1)*dza(i,k+1))+1.e-9
         ss = ((ux(i,k+1)-ux(i,k))**2 + (vx(i,k+1)-vx(i,k))**2) / dza(i,k+1)**2
         ss = max(ss, 1.e-8)
         govrthv = g/(0.5*(thvx(i,k+1)+thvx(i,k)))
         ri = govrthv*(thvx(i,k+1)-thvx(i,k))/(ss*dza(i,k+1))
         if(imvdif.eq.1.and.ndiff.ge.3)then
           if((qxci(i,k,1)+qxci(i,k,2)).gt.0.01e-3.and.(qxci(i           &
             ,k+1,1)+qxci(i,k+1,2)).gt.0.01e-3)then
!      in cloud
             qmean = 0.5*(qx(i,k)+qx(i,k+1))
             tmean = 0.5*(tx(i,k)+tx(i,k+1))
             alph  = xlv*qmean/rd/tmean
             chi   = xlv*xlv*qmean/cp/rv/tmean/tmean
!            ri    = (1.+alph)*(ri-g*g/ss/tmean/cp*((chi-alph)/(1.+chi)))
             ri = (1.+alph) * (ri - g*g/(ss*tmean*cp)*(chi-alph)/(1.+chi))
           endif
         endif
         zk = karman*zq(i,k+1)
         rlamdz = min(max(0.1*dza(i,k+1),rlam),300.) ! was constant 150 in EDMF
         rlamdz = min(dza(i,k+1),rlamdz)
         rl2 = (zk*rlamdz/(rlamdz+zk))**2
         dk = rl2*sqrt(ss)
         if(ri.lt.0.)then
! unstable regime (same as in EDMF?)
           ri = max(ri, rimin)
           sri = sqrt(-ri)
           xkzm(i,k) = dk*(1+8.*(-ri)/(1+1.746*sri))
           xkzh(i,k) = dk*(1+8.*(-ri)/(1+1.286*sri))
         else
! stable regime
           xkzh(i,k) = dk/(1+5.*ri)**2
           prnum = 1.0+2.1*ri ! set to 1 above PBL in EDMF
           prnum = min(prnum,prmax)
           xkzm(i,k) = xkzh(i,k)*prnum
         endif
!
         xkzm(i,k) = max(xkzm(i,k),xkzom(i,k))
         xkzh(i,k) = max(xkzh(i,k),xkzoh(i,k))
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         xkzml(i,k) = xkzm(i,k)
         xkzhl(i,k) = xkzh(i,k)
       endif
     enddo
   enddo
!
!     compute tridiagonal matrix elements for heat
!
   do k = kts,kte
     do i = its,ite
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
     enddo
   enddo
!
   do i = its,ite
     ad(i,1) = 1.
     f1(i,1) = thx(i,1)-300.+hfx(i)/cont/del(i,1)*dt2
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
! SJL       dsig   = p2d(i,k)-p2d(i,k+1)
       dsig   = p2m(i,k) - p2m(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*rdz
       if(pblflg(i).and.k.lt.kpbl(i)) then
         dsdzt = tem1*(-hgamt(i)*xkzh(i,k)/hpbl(i)-hfxpbl(i)*zfacent(i,k))
         f1(i,k)   = f1(i,k)+dtodsd*dsdzt
         f1(i,k+1) = thx(i,k+1)-300.-dtodsu*dsdzt
       elseif(pblflg(i).and.k.ge.kpbl(i).and.entfac(i,k).lt.4.6) then
         xkzh(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
         xkzh(i,k) = sqrt(xkzh(i,k)*xkzhl(i,k))
         xkzh(i,k) = max(xkzh(i,k),xkzoh(i,k))
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         f1(i,k+1) = thx(i,k+1)-300.
       else
         f1(i,k+1) = thx(i,k+1)-300.
       endif
       tem1   = dsig*xkzh(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
     enddo
   enddo
!
! copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
     enddo
   enddo
!
   call tridin_ysu(al,ad,cu,r1,au,f1,its,ite,kts,kte,1)
!
!     recover tendencies of heat
!
   do k = kte,kts,-1
     do i = its,ite
       ttend = (f1(i,k)-thx(i,k)+300.)*rdt*pi2d(i,k)
       ttnp(i,k) = ttnp(i,k)+ttend 
       dtsfc(i) = dtsfc(i)+ttend*cont*del(i,k)/pi2d(i,k)
     enddo
   enddo

   if (present(dkt)) then
      do k=kts,kte-1
      do i=its,ite
         dkt(i,k) = xkzh(i,k)
      enddo
      enddo
   endif
!
!     compute tridiagonal matrix elements for moisture, clouds, and gases
!
   do k = kts,kte
     do i = its,ite
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
     enddo
   enddo
!
   do ic = 1,ndiff
     do i = its,ite
       do k = kts,kte
         f3(i,k,ic) = 0.
       enddo
     enddo
   enddo
!
   do i = its,ite
     ad(i,1) = 1.
     f3(i,1,1) = qx(i,1)+qfx(i)*g/del(i,1)*dt2
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * kte
       do i = its,ite
         f3(i,1,ic) = qx(i,1+is)
       enddo
     enddo
   endif
!
   do k = kts,kte-1
     do i = its,ite
       if(k.ge.kpbl(i)) then
         xkzq(i,k) = xkzh(i,k)
       endif
     enddo
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
! SJL       dsig   = p2d(i,k)-p2d(i,k+1)
       dsig   = p2m(i,k) - p2m(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*rdz
       if(pblflg(i).and.k.lt.kpbl(i)) then
         dsdzq = tem1*(-qfxpbl(i)*zfacent(i,k))
         f3(i,k,1) = f3(i,k,1)+dtodsd*dsdzq
         f3(i,k+1,1) = qx(i,k+1)-dtodsu*dsdzq
       elseif(pblflg(i).and.k.ge.kpbl(i).and.entfac(i,k).lt.4.6) then
         xkzq(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
         xkzq(i,k) = sqrt(xkzq(i,k)*xkzhl(i,k))
         xkzq(i,k) = max(xkzq(i,k),xkzoh(i,k))
         xkzq(i,k) = min(xkzq(i,k),xkzmax)
         f3(i,k+1,1) = qx(i,k+1)
       else
         f3(i,k+1,1) = qx(i,k+1)
       endif
       tem1   = dsig*xkzq(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
!      exch_hx(i,k+1) = xkzh(i,k)
     enddo
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * kte
       do k = kts,kte-1
         do i = its,ite
           f3(i,k+1,ic) = qx(i,k+1+is)
         enddo
       enddo
     enddo
   endif
!
! copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
     enddo
   enddo
!
   do ic = 1,ndiff
     do k = kts,kte
       do i = its,ite
         r3(i,k,ic) = f3(i,k,ic)
       enddo
     enddo
   enddo
!
!     solve tridiagonal problem for moisture, clouds, and gases
!
   call tridin_ysu(al,ad,cu,r3,au,f3,its,ite,kts,kte,ndiff)
!
!     recover tendencies of heat and moisture
!
   do k = kte,kts,-1
     do i = its,ite
       qtend = (f3(i,k,1)-qx(i,k))*rdt
       qtnp(i,k) = qtnp(i,k)+qtend
       dqsfc(i) = dqsfc(i)+qtend*conq*del(i,k)
     enddo
   enddo


!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * kte
       do k = kte,kts,-1
         do i = its,ite
           qtend = (f3(i,k,ic)-qx(i,k+is))*rdt
           qtnp(i,k+is) = qtnp(i,k+is)+qtend
         enddo
       enddo
     enddo
   endif
!
!     compute tridiagonal matrix elements for momentum
!
   do i = its,ite
     do k = kts,kte
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
       f2(i,k) = 0.
     enddo
   enddo
!
   do i = its,ite
! paj: ctopo=1 if topo_wind=0 (default)
! mchen add this line to make sure NMM can still work with YSU PBL
!<---hns
!     if(present(ctopo)) then
!       ad(i,1) = 1.+ctopo(i)*ust(i)**2/wspd1(i)*rhox(i)*g/del(i,1)*dt2         &
!        *(wspd1(i)/wspd(i))**2
!     else
       ad(i,1) = 1.+ust(i)**2/wspd1(i)*rhox(i)*g/del(i,1)*dt2                  &
        *(wspd1(i)/wspd(i))**2
!     endif
!--->hns
     f1(i,1) = ux(i,1)+uox(i)*ust(i)**2*g/del(i,1)*dt2/wspd1(i)
     f2(i,1) = vx(i,1)+vox(i)*ust(i)**2*g/del(i,1)*dt2/wspd1(i)
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
! SJL       dsig   = p2d(i,k)-p2d(i,k+1)
       dsig   = p2m(i,k) - p2m(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*rdz
       if(pblflg(i).and.k.lt.kpbl(i))then
         dsdzu     = tem1*(-hgamu(i)*xkzm(i,k)/hpbl(i)-ufxpbl(i)*zfacent(i,k))
         dsdzv     = tem1*(-hgamv(i)*xkzm(i,k)/hpbl(i)-vfxpbl(i)*zfacent(i,k))
         f1(i,k)   = f1(i,k)+dtodsd*dsdzu
         f1(i,k+1) = ux(i,k+1)-dtodsu*dsdzu
         f2(i,k)   = f2(i,k)+dtodsd*dsdzv
         f2(i,k+1) = vx(i,k+1)-dtodsu*dsdzv
       elseif(pblflg(i).and.k.ge.kpbl(i).and.entfac(i,k).lt.4.6) then
         xkzm(i,k) = prpbl(i)*xkzh(i,k)
         xkzm(i,k) = sqrt(xkzm(i,k)*xkzml(i,k))
         xkzm(i,k) = max(xkzm(i,k),xkzom(i,k))
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         f1(i,k+1) = ux(i,k+1)
         f2(i,k+1) = vx(i,k+1)
       else
         f1(i,k+1) = ux(i,k+1)
         f2(i,k+1) = vx(i,k+1)
       endif
       tem1   = dsig*xkzm(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
     enddo
   enddo
!
! copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
       r2(i,k) = f2(i,k)
     enddo
   enddo
!
!     solve tridiagonal problem for momentum
!
   call tridi1n(al,ad,cu,r1,r2,au,f1,f2,its,ite,kts,kte,1)
!
!     recover tendencies of momentum
!
   do k = kte,kts,-1
     do i = its,ite
       utend = (f1(i,k)-ux(i,k))*rdt
       vtend = (f2(i,k)-vx(i,k))*rdt
       utnp(i,k) = utnp(i,k)+utend
       vtnp(i,k) = vtnp(i,k)+vtend
       dusfc(i) = dusfc(i) + utend*conwrc*del(i,k)
       dvsfc(i) = dvsfc(i) + vtend*conwrc*del(i,k)
     enddo
   enddo
!
! paj: ctopo2=1 if topo_wind=0 (default)
!
!<---hns
!   do i = its,ite
!     if(present(ctopo).and.present(ctopo2)) then ! mchen for NMM
!       u10(i) = ctopo2(i)*u10(i)+(1-ctopo2(i))*ux(i,1)
!       v10(i) = ctopo2(i)*v10(i)+(1-ctopo2(i))*vx(i,1)
!     endif !mchen
!   enddo
!--->hns
!
!---- end of vertical diffusion
!
!
!   compute tke dissipation rate
!
   if(dspheat) then
!
      do k = kts,kte-1
        do i = its,ite
          rdz  = 1./dza(i,k+1)
          dw2  = (ux(i,k)-ux(i,k+1))**2.+(vx(i,k)-vx(i,k+1))**2.
          shr2 = max(dw2,dw2min)*rdz*rdz
          ti   = 2./(tx(i,k)+tx(i,k+1))
          bf   = (thvx(i,k+1)-thvx(i,k))*rdz
          diss(i,k) = xkzm(i,k)*shr2-g*ti*xkzh(i,k)*bf
        enddo
      enddo
!
!     add dissipative heating at the first model layer
!
      do i = its,ite
         tem   = govrth(i)*sflux(i)
         tem1  = tem + ust(i)**2.*wspd1(i)/za(i,1)
         tem2  = 0.5 * (tem1+diss(i,1))
         tem2  = max(tem2, 0.)
         ttend = tem2 / cp
         ttnp(i,1) = ttnp(i,1)+0.5*ttend
      enddo
!
!     add dissipative heating above the first model layer
!
      do k = kts+1,kte-1
        do i = its,ite
          tem = 0.5 * (diss(i,k-1)+diss(i,k))
          tem  = max(tem, 0.)
          ttend = tem / cp
          ttnp(i,k) = ttnp(i,k) + 0.5*ttend
        enddo
      enddo
!
   endif
!
   do i = its,ite
     kpbl1d(i) = kpbl(i)
   enddo
!
!
   end subroutine ysupbl
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine tridi1n(cl,cm,cu,r1,r2,au,f1,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
   use machine, only : kind_phys
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in ) :: its,ite, kts,kte, nt
!
   real(kind=kind_phys), dimension( its:ite, kts+1:kte+1 )            , &
         intent(in   )  ::                                          cl
!
   real(kind=kind_phys), dimension( its:ite, kts:kte )                , &
         intent(in   )  ::                                          cm, &
                                                                    r1
   real(kind=kind_phys), dimension( its:ite, kts:kte,nt )             , &
         intent(in   )  ::                                          r2
!
   real(kind=kind_phys), dimension( its:ite, kts:kte )                , &
         intent(inout)  ::                                          au, &
                                                                    cu, &
                                                                    f1
   real(kind=kind_phys), dimension( its:ite, kts:kte,nt )             , &
         intent(inout)  ::                                          f2
!
   real(kind=kind_phys)    :: fk
   integer :: i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
   l = ite
   n = kte
!
   do i = its,l
     fk = 1./cm(i,1)
     au(i,1) = fk*cu(i,1)
     f1(i,1) = fk*r1(i,1)
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./cm(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
!
   do k = kts+1,n-1
     do i = its,l
       fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
       au(i,k) = fk*cu(i,k)
       f1(i,k) = fk*(r1(i,k)-cl(i,k)*f1(i,k-1))
     enddo
   enddo
!
   do it = 1,nt
     do k = kts+1,n-1
       do i = its,l
!        fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
!        f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
         f2(i,k,it) = (r2(i,k,it)-cl(i,k)*f2(i,k-1,it)) / (cm(i,k)-cl(i,k)*au(i,k-1))
       enddo
     enddo
   enddo
!
   do i = its,l
!    fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
!    f1(i,n) = fk*(r1(i,n)-cl(i,n)*f1(i,n-1))
     f1(i,n) = (r1(i,n)-cl(i,n)*f1(i,n-1)) / (cm(i,n)-cl(i,n)*au(i,n-1))
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
       f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
     enddo
   enddo
!
   do k = n-1,kts,-1
     do i = its,l
       f1(i,k) = f1(i,k)-au(i,k)*f1(i,k+1)
     enddo
   enddo
!
   do it = 1,nt
     do k = n-1,kts,-1
       do i = its,l
         f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
       enddo
     enddo
   enddo
!
   end subroutine tridi1n
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine tridin_ysu(cl,cm,cu,r2,au,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
   use machine, only : kind_phys
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in ) ::     its,ite, kts,kte, nt
!
   real(kind=kind_phys), dimension( its:ite, kts+1:kte+1 )            , &
         intent(in   )  ::                                          cl
!
   real(kind=kind_phys), dimension( its:ite, kts:kte )                , &
         intent(in   )  ::                                          cm
   real(kind=kind_phys), dimension( its:ite, kts:kte,nt )             , &
         intent(in   )  ::                                          r2
!
   real(kind=kind_phys), dimension( its:ite, kts:kte )                , &
         intent(inout)  ::                                          au, &
                                                                    cu
   real(kind=kind_phys), dimension( its:ite, kts:kte,nt )             , &
         intent(inout)  ::                                          f2
!
   real(kind=kind_phys) :: fk
   integer :: i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
   l = ite
   n = kte
!
   do it = 1,nt
     do i = its,l
       fk = 1./cm(i,1)
       au(i,1) = fk*cu(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
!
   do it = 1,nt
     do k = kts+1,n-1
       do i = its,l
         fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
         au(i,k) = fk*cu(i,k)
         f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
       enddo
     enddo
   enddo
!
   do it = 1,nt
     do i = its,l
!      fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
!      f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
       f2(i,n,it) = (r2(i,n,it)-cl(i,n)*f2(i,n-1,it))/(cm(i,n)-cl(i,n)*au(i,n-1))
     enddo
   enddo
!
   do it = 1,nt
     do k = n-1,kts,-1
       do i = its,l
         f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
       enddo
     enddo
   enddo
!
   end subroutine tridin_ysu
!-------------------------------------------------------------------------------
