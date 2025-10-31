!!!!! ==================================================================  !!!!!
!  subroutine 'satmedmfvdifq.f' computes subgrid vertical turbulence mixing
!  using scale-aware TKE-based moist eddy-diffusion mass-flux (EDMF) parameterization
!
!  --- Overview
!
!  Originally developed by Jongil Han at NOAA/NCEP/EMC 
!
!  1) For the convective boundary layer, the scheme adopts
!  EDMF parameterization (Siebesma et al., 2007) to take
!  into account nonlocal transport by large eddies (mfpbltq.f).
!
!  2) A new mass-flux parameterization for stratocumulus-top-induced turbulence
!  mixing has been introduced (previously, it was eddy diffusion form)
!  [mfscu.f].
!
!  3) For local turbulence mixing, a TKE closure model is used.
!
!  --- Updates
!
!  1) May 2019 by Jongil Han (EMC)
!    goals: to have better low-level inversion, 
!           to reduce the cold bias in lower troposphere,
!           to reduce the negative wind speed bias in upper troposphere
!  changes: reduce the minimum and maximum characteristic mixing lengths,
!           reduce core downdraft and updraft fractions,
!           change of updraft top height calculation,
!           reduce the background diffusivity with increasing surface layer stability (for inversion)

!  2) Jul 2019 by Kun Gao (GFDL; kun.gao@noaa.gov)
!      goal: to allow for tke advection
!    change: rearange tracers (q1) and their tendencies (rtg)
!            tke no longer needs to be the last tracer
!  3) Nov 2019 by Kun Gao
!     turn off non-local mixing for hydrometers to avoid unphysical negative values 
!  4) Jan 2020 by Kun Gao 
!     add rlmn2 parameter (set to 10.) to be consistent with EMC's version 
!  5) Jun 2020 by Kun Gao
!     a) disable the upper-limter on background diff. in inversion layer
!        over land points to be consistent with EMC's version
!     b) use different xkzm_m,xkzm_h for land, ocean and sea ice points
!     c) add option for turning off HB19 formula for surface backgroud diff. (do_dk_hb19)
!  
!  6) Jul 2020 from Jongil Han: significant revisions to improve SCu
!     a) revised xkzo and rlmnz in inversion layer
!     b) limited updraft overshooting
!
!  7) Aug 2025: Joseph Mouallem (GFDL)
!     broken into down and up sweeps for the temperature, moisture and tracers
!
!----------------------------------------------------------------------

      subroutine satmedmfvdifq_down(ix,im,km,ntrac,ntcw,ntiw,ntke,
     &     dv,du,tdt,rtg_in,u1,v1,t1,q1_in,
     &     swh,hlw,xmu,garea,zvfun,islimsk,
     &     psk,rbsoil,zorl,u10m,v10m,fm,fh,
     &     tsea,heat,evap,stress,spd1,kpbl,
     &     prsi,del,prsl,prslk,phii,phil,delt,
     &     dspheat,dusfc,dvsfc,dtsfc,dqsfc,hpbl,
     &     kinver,xkzm_mo,xkzm_ho,xkzm_ml,xkzm_hl,xkzm_mi,xkzm_hi,
     &     xkzm_s,xkzinv,rlmx,zolcru,cs0,
     &     do_dk_hb19,xkgdx,dspfac,bl_upfr,bl_dnfr,
     &     l2_diag_opt,use_lup_only,l1l2_blend_opt,use_l1_sfc,
     &     use_tke_ent_det,use_shear,
     &     dkt_out, flux_up, flux_dn, elm,
!    to be passed for the upward sweep
     &     au_out, f1_out, f2_out, diss_out) 

      use machine  , only : kind_phys
      use funcphys , only : fpvs
      use physcons, grav => con_g, rd => con_rd, cp => con_cp
     &,             rv => con_rv, hvap => con_hvap
     &,             hfus => con_hfus, fv => con_fvirt
     &,             eps => con_eps, epsm1 => con_epsm1
!
      implicit none
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      real(kind=kind_phys) au_out(im,km-1)
      real(kind=kind_phys) f1_out(im,km), f2_out(im,km*(ntrac-1)),
     &                     diss_out(im,km-1) 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      integer ix, im, km, ntrac, ntcw, ntiw, ntke, ntcw_new
      integer l2_diag_opt, l1l2_blend_opt
      integer kpbl(im), kinver(im), islimsk(im)
!
      real(kind=kind_phys) delt, xkzm_mo, xkzm_ho, xkzm_s, dspfac,
     &                     bl_upfr, bl_dnfr, xkzm_ml, xkzm_hl,
     &                     xkzm_mi, xkzm_hi
      real(kind=kind_phys) dv(im,km),     du(im,km),
     &                     tdt(im,km),    rtg(im,km,ntrac),
     &                     u1(ix,km),     v1(ix,km),
     &                     t1(ix,km),     q1(ix,km,ntrac),
     &                     swh(ix,km),    hlw(ix,km),
     &                     xmu(im),       garea(im),
     &                     zvfun(im),
     &                     psk(ix),       rbsoil(im),
     &                     zorl(im),      tsea(im),
     &                     u10m(im),      v10m(im),
     &                     fm(im),        fh(im),
     &                     evap(im),      heat(im),
     &                     stress(im),    spd1(im),
     &                     prsi(ix,km+1), del(ix,km),
     &                     prsl(ix,km),   prslk(ix,km),
     &                     phii(ix,km+1), phil(ix,km),
     &                     dusfc(im),     dvsfc(im),
     &                     dtsfc(im),     dqsfc(im),
     &                     hpbl(im),
     &                     q1_in(ix,km,ntrac),  
     &                     rtg_in(im,km,ntrac)
! kgao note - q1 and rtg are local var now 
!
      logical dspheat, do_dk_hb19, use_lup_only, use_l1_sfc,
     &        use_tke_ent_det, use_shear

      real(kind=kind_phys)::dkt_out(im,km),flux_up(im,km),flux_dn(im,km)
!
!----------------------------------------------------------------------
!***
!***  local variables
!***
      integer i,is,k,kk,n,ndt,km1,kmpbl,kmscu,ntrac1,kps
      integer lcld(im),kcld(im),krad(im),mrad(im)
      integer kx1(im), kpblx(im)
!
      real(kind=kind_phys) tke(im,km),  tkeh(im,km-1), e2(im,0:km)
!
      real(kind=kind_phys) theta(im,km),thvx(im,km),  thlvx(im,km),
     &                     qlx(im,km),  thetae(im,km),thlx(im,km),
     &                     slx(im,km),  svx(im,km),   qtx(im,km),
     &                     tvx(im,km),  pix(im,km),   radx(im,km-1),
     &                     dku(im,km-1),dkt(im,km-1), dkq(im,km-1),
     &                     cku(im,km-1),ckt(im,km-1)
!
      real(kind=kind_phys) plyr(im,km), rhly(im,km),  cfly(im,km),
     &                     qstl(im,km)
!
      real(kind=kind_phys) dtdz1(im), gdx(im),
     &                     phih(im),  phim(im),    
     &                     phims(im), prn(im,km-1),
     &                     rbdn(im),  rbup(im),    thermal(im),
     &                     ustar(im), wstar(im),   hpblx(im),
     &                     ust3(im),  wst3(im),
     &                     z0(im),    crb(im),
     &                     hgamt(im), hgamq(im),
     &                     wscale(im),vpert(im),
     &                     zol(im),   sflux(im),
     &                     tx1(im),   tx2(im),
     &                     vez0fun(im), tkemean(im), sumx(im)
!
      real(kind=kind_phys) vegflo, vegfup, z0lo, z0up, vc0, zc0, csmf
!
      real(kind=kind_phys) radmin(im)
!
      real(kind=kind_phys) zi(im,km+1),  zl(im,km),   zm(im,km),
     &                     xkzo(im,km),xkzmo(im,km),
     &                     xkzm_hx(im),  xkzm_mx(im),
     &                     ri(im,km-1),  tkmnz(im,km-1),
     &                     rdzt(im,km-1),rlmnz(im,km),
     &                     al(im,km-1),  ad(im,km),   au(im,km-1),
     &                     f1(im,km),    f2(im,km*(ntrac-1))
!
      real(kind=kind_phys) elm(im,km),   ele(im,km),
     &                     ckz(im,km),   chz(im,km),  frik(im),
     &                     diss(im,km-1),prod(im,km-1), 
     &                     bf(im,km-1),  shr2(im,km-1), wush(im,km),
     &                     xlamue(im,km-1), xlamde(im,km-1),
     &                     gotvx(im,km), rlam(im,km-1)
!
!   variables for updrafts (thermals)
!
      real(kind=kind_phys) tcko(im,km),  qcko(im,km,ntrac),
     &                     ucko(im,km),  vcko(im,km),
     &                     buou(im,km),  xmf(im,km)
!
!   variables for stratocumulus-top induced downdrafts
!
      real(kind=kind_phys) tcdo(im,km),  qcdo(im,km,ntrac),
     &                     ucdo(im,km),  vcdo(im,km),
     &                     buod(im,km),  xmfd(im,km)
!
      logical  pblflg(im), sfcflg(im), flg(im)
      logical  scuflg(im), pcnvflg(im)
      logical  mlenflg
!
!  pcnvflg: true for unstable pbl
!
      real(kind=kind_phys) aphi16,  aphi5,
     &                     wfac,    cfac,
     &                     gamcrt,  gamcrq, sfcfrac,
     &                     conq,    cont,   conw,
     &                     dsdz2,   dsdzt,  dkmax,
     &                     dsig,    dt2,    dtodsd,
     &                     dtodsu,  g,      factor, dz,
     &                     gocp,    gravi,  zol1,   zolcru,
     &                     buop,    shrp,   dtn,
     &                     prnum,   prmax,  prmin,  prtke,
     &                     prscu,   pr0,
     &                     dw2,     dw2min, zk,     
     &                     elmfac,  elefac, dspmax,
     &                     alp,     clwt,   cql,
     &                     f0,      robn,   crbmin, crbmax,
     &                     es,      qs,     value,  onemrh,
     &                     cfh,     gamma,  elocp,  el2orc,
     &                     epsi,    beta,   chx,    cqx,
     &                     rdt,     rdz,    qmin,   qlmin,
     &                     rimin,   rbcr,   rbint,  tdzmin,
     &                     rlmn,    rlmn1,  rlmn2,  
     &                     rlmx,    elmx,
     &                     ttend,   utend,  vtend,  qtend,
     &                     zfac,    zfmin,  vk,     spdk2,
     &                     tkmin,   tkminx, xkgdx,  xkzinv,
     &                     zlup,    zldn,   bsum,   cs0,
     &                     tem,     tem1,   tem2,   tem3,
     &                     ptem,    ptem0,  ptem1,  ptem2
!
      real(kind=kind_phys) ck0, ck1, ch0, ch1, ce0, rchck
!
      real(kind=kind_phys) qlcr, zstblmax, hcrinv
!
      real(kind=kind_phys) h1 
!!
      parameter(gravi=1.0/grav)
      parameter(g=grav)
      parameter(gocp=g/cp)
      parameter(cont=cp/g,conq=hvap/g,conw=1.0/g)  ! for del in pa
!     parameter(cont=1000.*cp/g,conq=1000.*hvap/g,conw=1000./g) !kpa
      parameter(elocp=hvap/cp,el2orc=hvap*hvap/(rv*cp))
      parameter(wfac=7.0,cfac=4.5)
      parameter(gamcrt=3.,gamcrq=0.,sfcfrac=0.1)
      parameter(vk=0.4,rimin=-100.)
!      parameter(rbcr=0.25,zolcru=-0.02,tdzmin=1.e-3)
      parameter(rbcr=0.25,tdzmin=1.e-3)
      parameter(rlmn=30.,rlmn1=5.,rlmn2=10.)
!      parameter(rlmx=300.,elmx=300.)
      parameter(prmin=0.25,prmax=4.0)
      parameter(pr0=1.0,prtke=1.0,prscu=0.67)
      parameter(f0=1.e-4,crbmin=0.15,crbmax=0.35)
      parameter(tkmin=1.e-9,tkminx=0.2,dspmax=10.0)
      parameter(qmin=1.e-8,qlmin=1.e-12,zfmin=1.e-8)
      parameter(aphi5=5.,aphi16=16.)
      parameter(elmfac=1.0,elefac=1.0,cql=100.)
      parameter(dw2min=1.e-4,dkmax=1000.)!,xkgdx=5000.)
      parameter(qlcr=3.5e-5,zstblmax=2500.) !,xkzinv=0.1)
      parameter(h1=0.33333333,hcrinv=250.)
      parameter(ck0=0.4,ck1=0.15,ch0=0.4,ch1=0.15)
!     parameter(ce0=0.4,cs0=0.5)
      parameter(ce0=0.4)
      parameter(rchck=1.5,ndt=20)
      parameter(vegflo=0.1,vegfup=1.0,z0lo=0.1,z0up=1.0)
      parameter(vc0=1.0,zc0=1.0)
      parameter(csmf=0.5)
!
!************************************************************************
      elmx = rlmx
      dt2  = delt
      rdt = 1. / dt2
      dkt_out = 0.
      flux_up = 0.
      flux_dn = 0. 
!
! kgao note (jul 2019) 
! the code was originally written assuming ntke=ntrac
! in this version ntke does not need to be equal to ntrac
! in the following we rearrange q1 (and rtg) so that tke is the last tracer
!
      !if(ntrac >= 3 ) then
        if(ntke == ntrac) then ! tke is the last tracer
          q1(:,:,:)  = q1_in(:,:,:)
          rtg(:,:,:) = rtg_in(:,:,:)
        else                   ! tke is not
          do kk = 1, ntke-1
             q1(:,:,kk)  = q1_in(:,:,kk)
             rtg(:,:,kk) = rtg_in(:,:,kk)
          enddo
          do kk = ntke+1, ntrac
             q1(:,:,kk-1)  = q1_in(:,:,kk)
             rtg(:,:,kk-1) = rtg_in(:,:,kk)
          enddo
          q1(:,:,ntrac)  = q1_in(:,:,ntke)
          rtg(:,:,ntrac) = rtg_in(:,:,ntke)
        endif
      !endif
!
      ntrac1 = ntrac - 1
      km1 = km - 1
      kmpbl = km / 2
      kmscu = km / 2
!
      do k=1,km
        do i=1,im
          zi(i,k) = phii(i,k) * gravi
          zl(i,k) = phil(i,k) * gravi
          xmf(i,k) = 0.
          xmfd(i,k) = 0.
          buou(i,k) = 0.
          buod(i,k) = 0.
          ckz(i,k) = ck1
          chz(i,k) = ch1
          rlmnz(i,k) = rlmn
        enddo
      enddo
      do i=1,im
        frik(i) = 1.0
      enddo
      do i=1,im
        zi(i,km+1) = phii(i,km+1) * gravi
      enddo
      do k=1,km
        do i=1,im
          zm(i,k) = zi(i,k+1)
        enddo
      enddo
! horizontal grid size
      do i=1,im
        gdx(i) = sqrt(garea(i))
      enddo
!
      do k=1,km
        do i=1,im
          tke(i,k) = max(q1_in(i,k,ntke), tkmin) ! tke at layer centers
        enddo
      enddo
      do k=1,km1
        do i=1,im
          tkeh(i,k) = 0.5 * (tke(i,k) + tke(i,k+1)) ! tke at interfaces
        enddo
      enddo
!
      do k = 1,km1
        do i=1,im
          rdzt(i,k) = 1.0 / (zl(i,k+1) - zl(i,k))
          prn(i,k)  = pr0
        enddo
      enddo
!
! set background diffusivities as a function of
!  horizontal grid size with xkzm_h & xkzm_m for gdx >= 25km
!  and 0.01 for gdx=5m, i.e.,
!  xkzm_hx = 0.01 + (xkzm_h - 0.01)/(xkgdx-5.) * (gdx-5.)
!  xkzm_mx = 0.01 + (xkzm_h - 0.01)/(xkgdx-5.) * (gdx-5.)
!
      do i=1,im
        kx1(i) = 1
        tx1(i) = 1.0 / prsi(i,1)
        tx2(i) = tx1(i)

        ! kgao change - set surface value of background diff (dk) below

        !if(gdx(i) >= xkgdx) then
        !  xkzm_hx(i) = xkzm_h
        !  xkzm_mx(i) = xkzm_m
        !else
        !  tem  = 1. / (xkgdx - 5.)
        !  tem1 = (xkzm_h - 0.01) * tem
        !  tem2 = (xkzm_m - 0.01) * tem
        !  ptem = gdx(i) - 5.
        !  xkzm_hx(i) = 0.01 + tem1 * ptem
        !  xkzm_mx(i) = 0.01 + tem2 * ptem
        !endif

        if (do_dk_hb19) then               ! use eq43 in HB2019

          if(gdx(i) >= xkgdx) then         ! resolution coarser than xkgdx
            if( islimsk(i) == 1 ) then     ! land points
              xkzm_hx(i) = xkzm_hl
              xkzm_mx(i) = xkzm_ml
            elseif ( islimsk(i) == 2 ) then! sea ice points
              xkzm_hx(i) = xkzm_hi
              xkzm_mx(i) = xkzm_mi
            else                           ! ocean points
              xkzm_hx(i) = xkzm_ho
              xkzm_mx(i) = xkzm_mo
            endif
          else                             ! resolution finer than xkgdx
            tem  = 1. / (xkgdx - 5.)
            if ( islimsk(i) == 1 ) then    ! land points
              tem1 = (xkzm_hl - 0.01) * tem
              tem2 = (xkzm_ml - 0.01) * tem
            elseif ( islimsk(i) == 2 ) then! sea ice points
              tem1 = (xkzm_hi - 0.01) * tem
              tem2 = (xkzm_mi - 0.01) * tem
            else                           ! ocean points
              tem1 = (xkzm_ho - 0.01) * tem
              tem2 = (xkzm_mo - 0.01) * tem
            endif
            ptem = gdx(i) - 5.
            xkzm_hx(i) = 0.01 + tem1 * ptem
            xkzm_mx(i) = 0.01 + tem2 * ptem
          endif

        else ! use values in the namelist; no res dependency

          if ( islimsk(i) == 1 ) then     ! land points
              xkzm_hx(i) = xkzm_hl
              xkzm_mx(i) = xkzm_ml
          elseif ( islimsk(i) == 2 ) then ! sea ice points
              xkzm_hx(i) = xkzm_hi
              xkzm_mx(i) = xkzm_mi
          else                            ! ocean points
              xkzm_hx(i) = xkzm_ho
              xkzm_mx(i) = xkzm_mo
          endif
        endif
      enddo

      do k = 1,km
        do i=1,im
          xkzo(i,k)  = 0.0
          xkzmo(i,k) = 0.0
          if (k < kinver(i)) then
!                                  vertical background diffusivity
            ptem      = prsi(i,k+1) * tx1(i)
            tem1      = 1.0 - ptem
            tem2      = tem1 * tem1 * 10.0
            tem2      = min(1.0, exp(-tem2))
            xkzo(i,k) = xkzm_hx(i) * tem2
!
            ptem      = prsl(i,k) * tx1(i)
            tem1      = 1.0 - ptem
            tem2      = tem1 * tem1 * 2.5
            tem2      = min(1.0, exp(-tem2))
            rlmnz(i,k)= rlmn * tem2
            rlmnz(i,k)= max(rlmnz(i,k), rlmn1)
!                                 vertical background diffusivity for momentum
            if (ptem >= xkzm_s) then
              xkzmo(i,k) = xkzm_mx(i)
              kx1(i)     = k + 1
            else
              if (k == kx1(i) .and. k > 1) tx2(i) = 1.0 / prsi(i,k)
              tem1 = 1.0 - prsi(i,k+1) * tx2(i)
              tem1 = tem1 * tem1 * 5.0
              xkzmo(i,k) = xkzm_mx(i) * min(1.0, exp(-tem1))
            endif
          endif
        enddo
      enddo
!
      do i = 1,im
         z0(i)    = 0.01 * zorl(i)
         dusfc(i) = 0.
         dvsfc(i) = 0.
         dtsfc(i) = 0.
         dqsfc(i) = 0.
         kpbl(i) = 1
         hpbl(i) = 0.
         kpblx(i) = 1
         hpblx(i) = 0.
         pblflg(i)= .true.
         sfcflg(i)= .true.
         if(rbsoil(i) > 0.) sfcflg(i) = .false.
         pcnvflg(i)= .false.
         scuflg(i)= .true.
         if(scuflg(i)) then
           radmin(i)= 0.
           mrad(i)  = km1
           krad(i)  = 1
           lcld(i)  = km1
           kcld(i)  = km1
         endif
      enddo
!           
! compute a function for green vegetation fraction and surface roughness
!       
      do i = 1,im
        !tem = (sigmaf(i) - vegflo) / (vegfup - vegflo)
        !tem = min(max(tem, 0.), 1.)
        !tem1 = sqrt(tem)
        !ptem = (z0(i) - z0lo) / (z0up - z0lo)
        !ptem = min(max(ptem, 0.), 1.)
        vez0fun(i) = 1. !(1. + vc0 * tem1) * (1. + zc0 * ptem)
      enddo
!
      do k=1,km
        do i=1,im
          pix(i,k)   = psk(i) / prslk(i,k)
          theta(i,k) = t1(i,k) * pix(i,k)
          if(ntiw > 0) then
            tem = max(q1_in(i,k,ntcw),qlmin)
            tem1 = max(q1_in(i,k,ntiw),qlmin)
            qlx(i,k) = tem + tem1
            ptem = hvap*tem + (hvap+hfus)*tem1
            slx(i,k)   = cp * t1(i,k) + phil(i,k) - ptem
          else
            qlx(i,k) = max(q1_in(i,k,ntcw),qlmin)
            slx(i,k)   = cp * t1(i,k) + phil(i,k) - hvap*qlx(i,k)
          endif
          tem2       = 1.+fv*max(q1(i,k,1),qmin)-qlx(i,k)
          thvx(i,k)  = theta(i,k) * tem2
          tvx(i,k)   = t1(i,k) * tem2
          qtx(i,k) = max(q1(i,k,1),qmin)+qlx(i,k)
          thlx(i,k)  = theta(i,k) - pix(i,k)*elocp*qlx(i,k)
          thlvx(i,k) = thlx(i,k) * (1. + fv * qtx(i,k))
          svx(i,k)   = cp * tvx(i,k)
          ptem1      = elocp * pix(i,k) * max(q1(i,k,1),qmin)
          thetae(i,k)= theta(i,k) +  ptem1
          gotvx(i,k) = g / tvx(i,k)
        enddo
      enddo
!
! compute an empirical cloud fraction based on 
!    Xu & Randall's (1996,JAS) study
!
      do k = 1, km
        do i = 1, im
          plyr(i,k)   = 0.01 * prsl(i,k)   ! pa to mb (hpa)
!  --- ...  compute relative humidity
          es  = 0.01 * fpvs(t1(i,k))       ! fpvs in pa
          qs  = max(qmin, eps * es / (plyr(i,k) + epsm1*es))
          rhly(i,k) = max(0.0, min(1.0, max(qmin, q1(i,k,1))/qs))
          qstl(i,k) = qs
        enddo
      enddo
!
      do k = 1, km
        do i = 1, im
          cfly(i,k) = 0.
          clwt = 1.0e-6 * (plyr(i,k)*0.001)
          if (qlx(i,k) > clwt) then
            onemrh= max(1.e-10, 1.0-rhly(i,k))
            tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)
            tem1  = cql / tem1
            value = max(min( tem1*qlx(i,k), 50.0), 0.0)
            tem2  = sqrt(sqrt(rhly(i,k)))
            cfly(i,k) = min(max(tem2*(1.0-exp(-value)), 0.0), 1.0)
          endif
        enddo
      enddo
!
!  compute buoyancy modified by clouds
!
      do k = 1, km1
        do i = 1, im
          tem  = 0.5 * (svx(i,k) + svx(i,k+1))
          tem1 = 0.5 * (t1(i,k) + t1(i,k+1))
          tem2 = 0.5 * (qstl(i,k) + qstl(i,k+1))
          cfh  = min(cfly(i,k+1),0.5*(cfly(i,k)+cfly(i,k+1)))
          alp  = g / tem
          gamma = el2orc * tem2 / (tem1**2)
          epsi  = tem1 / elocp
          beta  = (1. + gamma*epsi*(1.+fv)) / (1. + gamma)
          chx   = cfh * alp * beta + (1. - cfh) * alp
          cqx   = cfh * alp * hvap * (beta - epsi)
          cqx   = cqx + (1. - cfh) * fv * g
          ptem1 = (slx(i,k+1)-slx(i,k))*rdzt(i,k)
          ptem2 = (qtx(i,k+1)-qtx(i,k))*rdzt(i,k)
          bf(i,k) = chx * ptem1 + cqx * ptem2
        enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k=1,km1
        do i=1,im
          dku(i,k)  = 0.
          dkt(i,k)  = 0.
          dkq(i,k)  = 0.
          cku(i,k)  = 0.
          ckt(i,k)  = 0.
          tem       = zi(i,k+1)-zi(i,k)
          radx(i,k) = tem*(swh(i,k)*xmu(i)+hlw(i,k))
        enddo
      enddo
!
      do i = 1,im
         sflux(i)  = heat(i) + evap(i)*fv*theta(i,1)
         if(.not.sfcflg(i) .or. sflux(i) <= 0.) pblflg(i)=.false.
      enddo
!
!  compute critical bulk richardson number 
!
      do i = 1,im
        if(pblflg(i)) then
!         thermal(i) = thvx(i,1)
          thermal(i) = thlvx(i,1)
          crb(i) = rbcr
        else
          thermal(i) = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
          tem = sqrt(u10m(i)**2+v10m(i)**2)
          tem = max(tem, 1.)
          robn = tem / (f0 * z0(i))
          tem1 = 1.e-7 * robn
          crb(i) = 0.16 * (tem1 ** (-0.18))
          crb(i) = max(min(crb(i), crbmax), crbmin)
        endif
      enddo
!
      do i=1,im
         dtdz1(i)  = dt2 / (zi(i,2)-zi(i,1))
      enddo
!
      do i=1,im
         ustar(i) = sqrt(stress(i))
      enddo
!
!  compute buoyancy (bf) and winshear square
!
      do k = 1, km1
      do i = 1, im
         rdz  = rdzt(i,k)
!        bf(i,k) = gotvx(i,k)*(thvx(i,k+1)-thvx(i,k))*rdz
         dw2  = (u1(i,k)-u1(i,k+1))**2
     &        + (v1(i,k)-v1(i,k+1))**2
         shr2(i,k) = max(dw2,dw2min)*rdz*rdz
         ri(i,k) = max(bf(i,k)/shr2(i,k),rimin)
      enddo
      enddo
!
! find pbl height based on bulk richardson number (mrf pbl scheme)
!   and also for diagnostic purpose
!
      do i=1,im
         flg(i) = .false.
         rbup(i) = rbsoil(i)
      enddo
!
      do k = 1, kmpbl
      do i = 1, im
        if(.not.flg(i)) then
          rbdn(i) = rbup(i)
          spdk2   = max((u1(i,k)**2+v1(i,k)**2),1.)
!         rbup(i) = (thvx(i,k)-thermal(i))*
!    &              (g*zl(i,k)/thvx(i,1))/spdk2
          rbup(i) = (thlvx(i,k)-thermal(i))*
     &              (g*zl(i,k)/thlvx(i,1))/spdk2
          kpblx(i) = k
          flg(i)  = rbup(i) > crb(i)
        endif
      enddo
      enddo
      do i = 1,im
        if(kpblx(i) > 1) then
          k = kpblx(i)
          if(rbdn(i) >= crb(i)) then
            rbint = 0.
          elseif(rbup(i) <= crb(i)) then
            rbint = 1.
          else
            rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
          endif
          hpblx(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
          if(hpblx(i) < zi(i,kpblx(i))) kpblx(i)=kpblx(i)-1
        else
          hpblx(i) = zl(i,1)
          kpblx(i) = 1
        endif
        hpbl(i) = hpblx(i)
        kpbl(i) = kpblx(i)
        if(kpbl(i) <= 1) pblflg(i)=.false.
      enddo
!
!  compute mean tke within pbl
!
      do i = 1, im
        sumx(i) = 0.
        tkemean(i) = 0.
      enddo
      do k = 1, kmpbl
      do i = 1, im
        if(k < kpbl(i)) then
          dz = zi(i,k+1) - zi(i,k)
          tkemean(i) = tkemean(i) + tke(i,k) * dz
          sumx(i) = sumx(i) + dz
        endif
      enddo
      enddo
      do i = 1, im
        if(tkemean(i) > 0. .and. sumx(i) > 0.) then
          tkemean(i) = tkemean(i) / sumx(i)
        endif
      enddo
!
!  compute wind shear term as a sink term for updraft and downdraft
!  velocity
!
      kps = max(kmpbl, kmscu)
      do k = 2, kps
      do i = 1, im
        dz = zi(i,k+1) - zi(i,k)
        tem = (0.5*(u1(i,k-1)-u1(i,k+1))/dz)**2
        tem1 = tem+(0.5*(v1(i,k-1)-v1(i,k+1))/dz)**2
        wush(i,k) = csmf * sqrt(tem1)
      enddo
      enddo
!
!     compute similarity parameters
!
      do i=1,im
         zol(i) = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
         if(sfcflg(i)) then
           zol(i) = min(zol(i),-zfmin)
         else
           zol(i) = max(zol(i),zfmin)
         endif
!
         zol1 = zol(i)*sfcfrac*hpbl(i)/zl(i,1)
         if(sfcflg(i)) then
           tem     = 1.0 / (1. - aphi16*zol1)
           phih(i) = sqrt(tem)
           phim(i) = sqrt(phih(i))
           tem1    = 1.0 / (1. - aphi16*zol(i))
           phims(i) = sqrt(sqrt(tem1))
         else
           phim(i) = 1. + aphi5*zol1
           phih(i) = phim(i)
           phims(i) = 1. + aphi5*zol(i)
         endif
      enddo
!
      do i=1,im
        if(pblflg(i)) then
          if(zol(i) < zolcru) then
            pcnvflg(i) = .true.
          endif
          wst3(i) = gotvx(i,1)*sflux(i)*hpbl(i)
          wstar(i)= wst3(i)**h1
          ust3(i) = ustar(i)**3.
          wscale(i)=(ust3(i)+wfac*vk*wst3(i)*sfcfrac)**h1
          ptem = ustar(i)/aphi5
          wscale(i) = max(wscale(i),ptem)
        endif
      enddo
!
! compute a thermal excess
!
      do i = 1,im
         if(pcnvflg(i)) then
           hgamt(i) = heat(i)/wscale(i)
           hgamq(i) = evap(i)/wscale(i)
           vpert(i) = hgamt(i) + hgamq(i)*fv*theta(i,1)
           vpert(i) = max(vpert(i),0.)
           tem = min(cfac*vpert(i),gamcrt)
           thermal(i)= thermal(i) + tem !jih jul2020
         endif
      enddo
!
!  enhance the pbl height by considering the thermal excess
!     (overshoot pbl top) -- jih jul2020
!
      do i=1,im
         flg(i)  = .true.
         if(pcnvflg(i)) then
           flg(i)  = .false.
           rbup(i) = rbsoil(i)
         endif
      enddo
      do k = 2, kmpbl
      do i = 1, im
        if(.not.flg(i)) then
          rbdn(i) = rbup(i)
          spdk2   = max((u1(i,k)**2+v1(i,k)**2),1.)
          rbup(i) = (thlvx(i,k)-thermal(i))*
     &              (g*zl(i,k)/thlvx(i,1))/spdk2
          kpbl(i) = k
          flg(i)  = rbup(i) > crb(i)
        endif
      enddo
      enddo
      do i = 1,im
        if(pcnvflg(i)) then
           k = kpbl(i)
           if(rbdn(i) >= crb(i)) then
             rbint = 0.
           elseif(rbup(i) <= crb(i)) then
             rbint = 1.
           else
             rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
           endif
           hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
           if(hpbl(i) < zi(i,kpbl(i))) then
             kpbl(i) = kpbl(i) - 1
           endif
           if(kpbl(i) <= 1) then
              pcnvflg(i) = .false.
              pblflg(i) = .false.
           endif
        endif
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  look for stratocumulus
!
      do i=1,im
         flg(i)  = scuflg(i)
      enddo
      do k = 1, km1
        do i=1,im
          if(flg(i).and.zl(i,k) >= zstblmax) then
             lcld(i)=k
             flg(i)=.false.
          endif
      enddo
      enddo
      do i = 1, im
        flg(i)=scuflg(i)
      enddo
      do k = kmscu,1,-1
      do i = 1, im
        if(flg(i) .and. k <= lcld(i)) then
          if(qlx(i,k) >= qlcr) then
             kcld(i)=k
             flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, im
        if(scuflg(i) .and. kcld(i)==km1) scuflg(i)=.false.
      enddo
!
      do i = 1, im
        flg(i)=scuflg(i)
      enddo
      do k = kmscu,1,-1
      do i = 1, im
        if(flg(i) .and. k <= kcld(i)) then
          if(qlx(i,k) >= qlcr) then
            if(radx(i,k) < radmin(i)) then
              radmin(i)=radx(i,k)
              krad(i)=k
            endif
          else
            flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, im
        if(scuflg(i) .and. krad(i) <= 1) scuflg(i)=.false.
        if(scuflg(i) .and. radmin(i)>=0.) scuflg(i)=.false.
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute components for mass flux mixing by large thermals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k = 1, km
        do i = 1, im
          if(pcnvflg(i)) then
            tcko(i,k) = t1(i,k)
            ucko(i,k) = u1(i,k)
            vcko(i,k) = v1(i,k)
          endif
          if(scuflg(i)) then
            tcdo(i,k) = t1(i,k)
            ucdo(i,k) = u1(i,k)
            vcdo(i,k) = v1(i,k)
          endif
        enddo
      enddo
      do kk = 1, ntrac1
      do k = 1, km
        do i = 1, im
          if(pcnvflg(i)) then
            qcko(i,k,kk) = q1(i,k,kk)
          endif
          if(scuflg(i)) then
            qcdo(i,k,kk) = q1(i,k,kk)
          endif
        enddo
      enddo
      enddo
! kgao note - change ntcw if q1 is rearranged
      if (ntke > ntcw) then
         ntcw_new = ntcw
      else
         ntcw_new = ntcw-1
      endif
! EDMF parameterization Siebesma et al.(2007) 
      call mfpbltq(im,ix,km,kmpbl,ntcw_new,ntrac1,dt2,
     &    pcnvflg,zl,zm,q1,t1,u1,v1,plyr,pix,thlx,thvx,
     &    gdx,hpbl,kpbl,vpert,buou,
     &    use_shear,wush,
     &    use_tke_ent_det,tkemean,vez0fun,xmf,
     &    tcko,qcko,ucko,vcko,xlamue,bl_upfr)
! mass-flux parameterization for stratocumulus-top-induced turbulence mixing
      call mfscuq(im,ix,km,kmscu,ntcw_new,ntrac1,dt2,
     &    scuflg,zl,zm,q1,t1,u1,v1,plyr,pix,
     &    thlx,thvx,thlvx,gdx,thetae,
     &    krad,mrad,radmin,buod,
     &    use_shear,wush,
     &    use_tke_ent_det,tkemean,vez0fun,xmfd,
     &    tcdo,qcdo,ucdo,vcdo,xlamde,bl_dnfr)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   compute prandtl number and exchange coefficient varying with height
!
      do k = 1, kmpbl
        do i = 1, im
          if(k < kpbl(i)) then
            tem = phih(i)/phim(i)
            ptem = sfcfrac*hpbl(i)
            tem1 = max(zi(i,k+1)-ptem, 0.)
            tem2 = tem1 / (hpbl(i) - ptem)
            if(pcnvflg(i)) then
              tem = min(tem, pr0)
              prn(i,k) = tem + (pr0 - tem) * tem2
            else
              tem = max(tem, pr0)
              prn(i,k) = tem
            endif
            prn(i,k) = min(prn(i,k),prmax)
            prn(i,k) = max(prn(i,k),prmin)
!
            ckz(i,k) = ck0 + (ck1 - ck0) * tem2
            ckz(i,k) = max(min(ckz(i,k), ck0), ck1)
            chz(i,k) = ch0 + (ch1 - ch0) * tem2
            chz(i,k) = max(min(chz(i,k), ch0), ch1)
!
          endif
        enddo
      enddo
!
! Above a threshold height (hcrinv), the background vertical diffusivities & mixing length 
!    in the inversion layers are set to much smaller values (xkzinv & rlmn2)
!
! Below the threshold height (hcrinv), the background vertical diffusivities & mixing length 
!    in the inversion layers are increased with increasing roughness length & vegetation fraction
!
      do k = 1,km1
        do i=1,im
          if(zi(i,k+1) > hcrinv) then
            tem1 = tvx(i,k+1)-tvx(i,k)
            if(tem1 >= 0. .and. islimsk(i) == 0) then ! kgao note: only apply limiter over ocean points
              xkzo(i,k)  = min(xkzo(i,k), xkzinv)
              xkzmo(i,k) = min(xkzmo(i,k), xkzinv)
              rlmnz(i,k) = min(rlmnz(i,k), rlmn2)
            endif
          else
            tem1 = tvx(i,k+1)-tvx(i,k)
            if(tem1 > 0.) then
              ptem = xkzo(i,k) * zvfun(i)
              xkzo(i,k) = min(max(ptem, xkzinv), xkzo(i,k))
              ptem = xkzmo(i,k) * zvfun(i)
              xkzmo(i,k) = min(max(ptem, xkzinv), xkzmo(i,k))
              ptem = rlmnz(i,k) * zvfun(i)
              rlmnz(i,k) = min(max(ptem, rlmn2), rlmnz(i,k))
            endif
          endif
        enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute an asymtotic mixing length

      do k = 1, km1
        do i = 1, im

          if (l2_diag_opt == 0) then
          ! kgao 12/08/2023: original method as in Han and Bretherton 2019
          ! but additionally considers shear effect
          zlup = 0.0
          bsum = 0.0
          mlenflg = .true.
          do n = k, km1
            if(mlenflg) then
              dz = zl(i,n+1) - zl(i,n)
              tem3=((u1(i,n+1)-u1(i,n))/dz)**2
              tem3=tem3+((v1(i,n+1)-v1(i,n))/dz)**2
              tem3=cs0*sqrt(tem3)*sqrt(tke(i,k))
              ptem = (gotvx(i,n)*(thvx(i,n+1)-thvx(i,k))+tem3)*dz
              bsum = bsum + ptem
              zlup = zlup + dz
              if(bsum >= tke(i,k)) then
                if(ptem >= 0.) then
                  tem2 = max(ptem, zfmin)
                else
                  tem2 = min(ptem, -zfmin)
                endif
                ptem1 = (bsum - tke(i,k)) / tem2
                zlup = zlup - ptem1 * dz 
                zlup = max(zlup, 0.)
                mlenflg = .false.
              endif
            endif
          enddo
          zldn = 0.0
          bsum = 0.0
          mlenflg = .true.
          do n = k, 1, -1
            if(mlenflg) then
              if(n == 1) then
                dz = zl(i,1)
                tem1 = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
                !jih jul2020
                tem3 = (u1(i,1)/dz)**2
                tem3 = tem3+(v1(i,1)/dz)**2
                tem3 = cs0*sqrt(tem3)*sqrt(tke(i,1))
              else
                dz = zl(i,n) - zl(i,n-1)
                tem1 = thvx(i,n-1)
!               tem1 = thlvx(i,n-1)
                tem3 = ((u1(i,n)-u1(i,n-1))/dz)**2
                tem3 = tem3+((v1(i,n)-v1(i,n-1))/dz)**2
                tem3 = cs0*sqrt(tem3)*sqrt(tke(i,k))
              endif
              ptem = (gotvx(i,n)*(thvx(i,k)-tem1)+tem3)*dz
              bsum = bsum + ptem
              zldn = zldn + dz
              if(bsum >= tke(i,k)) then
                if(ptem >= 0.) then
                  tem2 = max(ptem, zfmin)
                else
                  tem2 = min(ptem, -zfmin)
                endif
                ptem1 = (bsum - tke(i,k)) / tem2
                zldn = zldn - ptem1 * dz 
                zldn = max(zldn, 0.)
                mlenflg = .false.
              endif
            endif
          enddo

          else if (l2_diag_opt == 1) then
          ! kgao 12/08/2023: a new method for diagnosing L2
          zlup = 0.0
          mlenflg = .true.
          e2(i,k) = max(2.*tke(i,k), 0.001)
          do n = k, km1
            if(mlenflg) then
              dz = zl(i,n+1) - zl(i,n)
              tem1 = 2.*gotvx(i,n+1)*(thvx(i,k)-thvx(i,n+1))
              tem2 = cs0*sqrt(e2(i,n))*sqrt(shr2(i,n))
              e2(i,n+1) = e2(i,n) + (tem1 - tem2) * dz
              zlup = zlup + dz
              if(e2(i,n+1) < 0.) then
                ptem = e2(i,n+1) / (e2(i,n+1) - e2(i,n))
                zlup = zlup - ptem * dz
                zlup = max(zlup, 0.)
                mlenflg = .false.
              endif
            endif
          enddo
          zldn = 0.0
          mlenflg = .true.
          do n = k, 1, -1
            if(mlenflg) then
              if(n == 1) then
                dz = zl(i,1)
                tem = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
                tem1 = 2.*gotvx(i,n)*(tem-thvx(i,k))
                tem2 = ustar(i)*phims(i)/(vk*dz)
                tem2 = cs0*sqrt(e2(i,n))*tem2
                e2(i,n-1) = e2(i,n) + (tem1 - tem2) * dz
              else
                dz = zl(i,n) - zl(i,n-1)
                tem1 = 2.*gotvx(i,n-1)*(thvx(i,n-1)-thvx(i,k))
                tem2 = cs0*sqrt(e2(i,n))*sqrt(shr2(i,n-1))
                e2(i,n-1) = e2(i,n) + (tem1 - tem2) * dz
              endif
              zldn = zldn + dz
              if(e2(i,n-1) < 0.) then
                ptem = e2(i,n-1) / (e2(i,n-1) - e2(i,n))
                zldn = zldn - ptem * dz
                zldn = max(zldn, 0.)
                mlenflg = .false.
              endif
            endif
          enddo
          endif ! end-if of l2_diag_opt
!
          tem = 0.5 * (zi(i,k+1)-zi(i,k))
          tem1 = min(tem, rlmnz(i,k))
!

          ! kgao 08/29/23: add option to use L_up as L2 
          ! zldn is strongly limited by the layer height for the near-surface levels
          ! it is not physical to use limiter below because zk already considers this factor
          if (use_lup_only) then
             ptem2 = zlup
          else
             ptem2 = min(zlup,zldn)
          endif

          rlam(i,k) = elmfac * ptem2
          rlam(i,k) = max(rlam(i,k), tem1)
          rlam(i,k) = min(rlam(i,k), rlmx)
!
          ptem2 = sqrt(zlup*zldn)

          ele(i,k) = elefac * ptem2
          ele(i,k) = max(ele(i,k), tem1)
          ele(i,k) = min(ele(i,k), elmx)
!
        enddo
      enddo
!
      do k = 1, km1
        do i = 1, im
          tem = vk * zl(i,k)
          if (zol(i) < 0.) then
            ptem = 1. - 100. * zol(i)
            ptem1 = ptem**0.2
            zk = tem * ptem1
          elseif (zol(i) >= 1.) then
            zk = tem / 3.7
          else
            ptem = 1. + 2.7 * zol(i)
            zk = tem / ptem
          endif 

          ! kgao 12/08/2023: introduce multiple l1 and l2 blending options
          if (l1l2_blend_opt == 0) then
            ! original as in HK19
            elm(i,k) = zk*rlam(i,k)/(rlam(i,k)+zk)

          else if ( l1l2_blend_opt == 1 ) then
            ! HAFA method as in wang et al 2023 WaF; use zk as elm within surface layer
            tem = 1.
            if ( sfcflg(i) .and. hpbl(i) > 200.) then
               tem1 = min(100., hpbl(i)*0.05)                        ! sfc layer height
               if (zl(i,k) < tem1) then                              ! for layers below sfc layer
                  tem = 0.
               elseif ( zl(i,k) >= tem1 .and. zl(i,k) < 2*tem1) then ! transition layers 
                  tem = (zl(i,k) - tem1) / tem1
               endif
            endif
            elm(i,k) = 1. / ( 1./zk + tem * 1./rlam(i,k) )

          else if (l1l2_blend_opt == 2) then
            ! HAFB blending method
            elm(i,k) = sqrt( 1.0/( 1.0/(zk**2)+1.0/(rlam(i,k)**2) ) )
          endif

          ! kgao 12/08/2023: use l1 as l at the lowest layer
          if (use_l1_sfc) then 
             if(k == 1) elm(i,k)=zk
          endif

          dz = zi(i,k+1) - zi(i,k)
          tem = max(gdx(i),dz)
          elm(i,k) = min(elm(i,k), tem)
          ele(i,k) = min(ele(i,k), tem)
!
        enddo
      enddo
      do i = 1, im
        elm(i,km) = elm(i,km1)
        ele(i,km) = ele(i,km1)
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute eddy diffusivities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k = 1, km1
        do i = 1, im
           xkzo(i,k) = 0.5 * (xkzo(i,k) + xkzo(i,k+1))
           xkzmo(i,k) = 0.5 * (xkzmo(i,k) + xkzmo(i,k+1))
        enddo
      enddo
      do k = 1, km1
        do i = 1, im
           tem = 0.5 * (elm(i,k) + elm(i,k+1))
           tem = tem * sqrt(tkeh(i,k))
           if(k < kpbl(i)) then
             if(pcnvflg(i)) then
               dku(i,k) = ckz(i,k) * tem
               dkt(i,k) = dku(i,k) / prn(i,k)
             else
               if(ri(i,k) < 0.) then ! unstable regime
                 dku(i,k) = ckz(i,k) * tem
                 dkt(i,k) = dku(i,k) / prn(i,k)
               else             ! stable regime
                 dkt(i,k) = chz(i,k) * tem
                 dku(i,k) = dkt(i,k) * prn(i,k)
               endif
             endif
           else
              if(ri(i,k) < 0.) then ! unstable regime
                dku(i,k) = ck1 * tem
                dkt(i,k) = rchck * dku(i,k)
              else             ! stable regime
                dkt(i,k) = ch1 * tem
                prnum = 1.0 + 2.1 * ri(i,k)
                prnum = min(prnum,prmax)
                dku(i,k) = dkt(i,k) * prnum
              endif
           endif
!
           if(scuflg(i)) then
             if(k >= mrad(i) .and. k < krad(i)) then
                tem1 = ckz(i,k) * tem
                ptem1 = tem1 / prscu
                dku(i,k) = max(dku(i,k), tem1)
                dkt(i,k) = max(dkt(i,k), ptem1)
             endif
           endif
!
           dkq(i,k) = prtke * dkt(i,k)
!
           dkt(i,k) = min(dkt(i,k),dkmax)
           dkt(i,k) = max(dkt(i,k),xkzo(i,k))
           dkq(i,k) = min(dkq(i,k),dkmax)
           dkq(i,k) = max(dkq(i,k),xkzo(i,k))
           dku(i,k) = min(dku(i,k),dkmax)
           dku(i,k) = max(dku(i,k),xkzmo(i,k))
!
        enddo
      enddo
!
!  compute a minimum TKE deduced from background diffusivity for momentum.
!
      do k = 1, km1
        do i = 1, im
          if(k == 1) then
            tem = ckz(i,1)
            tem1 = 0.5 * xkzmo(i,1)
          else
            tem = 0.5 * (ckz(i,k-1) + ckz(i,k))
            tem1 = 0.5 * (xkzmo(i,k-1) + xkzmo(i,k))
          endif
          ptem = tem1 / (tem * elm(i,k))
          tkmnz(i,k) = ptem * ptem
          tkmnz(i,k) = min(tkmnz(i,k), tkminx)
          tkmnz(i,k) = max(tkmnz(i,k), tkmin)
        enddo
      enddo
! kgao
      do k=1,km1
        do i=1,im
           dkt_out(i,k) = dkt(i,k)
       enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute buoyancy and shear productions of tke 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k = 1, km1
        do i = 1, im
          if (k == 1) then
            tem = -dkt(i,1) * bf(i,1)
!           if(pcnvflg(i)) then
!             ptem1 = xmf(i,1) * buou(i,1)
!           else
              ptem1 = 0.
!           endif
            if(scuflg(i) .and. mrad(i) == 1) then
              ptem2 = xmfd(i,1) * buod(i,1)
            else
              ptem2 = 0.
            endif
            tem = tem + ptem1 + ptem2
            buop = 0.5 * (gotvx(i,1) * sflux(i) + tem)
!
            tem1 = dku(i,1) * shr2(i,1)
!
            tem = (u1(i,2)-u1(i,1))*rdzt(i,1)
!           if(pcnvflg(i)) then
!             ptem = xmf(i,1) * tem
!             ptem1 = 0.5 * ptem * (u1(i,2)-ucko(i,2))
!           else
              ptem1 = 0.
!           endif
            if(scuflg(i) .and. mrad(i) == 1) then
              ptem = ucdo(i,1)+ucdo(i,2)-u1(i,1)-u1(i,2)
              ptem = 0.5 * tem * xmfd(i,1) * ptem
            else
              ptem = 0.
            endif
            ptem1 = ptem1 + ptem
!
            tem = (v1(i,2)-v1(i,1))*rdzt(i,1)
!           if(pcnvflg(i)) then
!             ptem = xmf(i,1) * tem
!             ptem2 = 0.5 * ptem * (v1(i,2)-vcko(i,2))
!           else
              ptem2 = 0.
!           endif
            if(scuflg(i) .and. mrad(i) == 1) then
              ptem = vcdo(i,1)+vcdo(i,2)-v1(i,1)-v1(i,2)
              ptem = 0.5 * tem * xmfd(i,1) * ptem
            else
              ptem = 0.
            endif
            ptem2 = ptem2 + ptem
!
!           tem2 = stress(i)*spd1(i)/zl(i,1)
            tem2 = stress(i)*ustar(i)*phim(i)/(vk*zl(i,1))
            shrp = 0.5 * (tem1 + ptem1 + ptem2 + tem2)
          else
            tem1 = -dkt(i,k-1) * bf(i,k-1)
            tem2 = -dkt(i,k) * bf(i,k)
            tem  = 0.5 * (tem1 + tem2)
            if(pcnvflg(i) .and. k <= kpbl(i)) then
              ptem = 0.5 * (xmf(i,k-1) + xmf(i,k))
              ptem1 = ptem * buou(i,k)
            else
              ptem1 = 0.
            endif
            if(scuflg(i)) then
              if(k >= mrad(i) .and. k < krad(i)) then
                ptem0 = 0.5 * (xmfd(i,k-1) + xmfd(i,k))
                ptem2 = ptem0 * buod(i,k)
              else
                ptem2 = 0.
              endif
            else
              ptem2 = 0.
            endif
            buop = tem + ptem1 + ptem2
!
            tem1 = dku(i,k-1) * shr2(i,k-1)
            tem2 = dku(i,k) * shr2(i,k)
            tem  = 0.5 * (tem1 + tem2)
            tem1 = (u1(i,k+1)-u1(i,k))*rdzt(i,k)
            tem2 = (u1(i,k)-u1(i,k-1))*rdzt(i,k-1)
            if(pcnvflg(i) .and. k <= kpbl(i)) then
              ptem = xmf(i,k) * tem1 + xmf(i,k-1) * tem2
              ptem1 = 0.5 * ptem * (u1(i,k)-ucko(i,k))
            else
              ptem1 = 0.
            endif
            if(scuflg(i)) then
              if(k >= mrad(i) .and. k < krad(i)) then
                ptem0 = xmfd(i,k) * tem1 + xmfd(i,k-1) * tem2
                ptem2 = 0.5 * ptem0 * (ucdo(i,k)-u1(i,k))
              else
                ptem2 = 0.
              endif
            else
              ptem2 = 0.
            endif
            shrp = tem + ptem1 + ptem2
            tem1 = (v1(i,k+1)-v1(i,k))*rdzt(i,k)
            tem2 = (v1(i,k)-v1(i,k-1))*rdzt(i,k-1)
            if(pcnvflg(i) .and. k <= kpbl(i)) then
              ptem = xmf(i,k) * tem1 + xmf(i,k-1) * tem2
              ptem1 = 0.5 * ptem * (v1(i,k)-vcko(i,k))
            else
              ptem1 = 0.
            endif
            if(scuflg(i)) then
              if(k >= mrad(i) .and. k < krad(i)) then
                ptem0 = xmfd(i,k) * tem1 + xmfd(i,k-1) * tem2
                ptem2 = 0.5 * ptem0 * (vcdo(i,k)-v1(i,k))
              else
                ptem2 = 0.
              endif
            else
              ptem2 = 0.
            endif
            shrp = shrp + ptem1 + ptem2
          endif
          prod(i,k) = buop + shrp
        enddo
      enddo
!
!----------------------------------------------------------------------
!     first predict tke due to tke production & dissipation(diss) 
!
      dtn = dt2 / float(ndt)
      do n = 1, ndt
      do k = 1,km1
        do i=1,im
           tem = sqrt(tke(i,k))
           ptem = ce0 / ele(i,k)
           diss(i,k) = ptem * tke(i,k) * tem
           tem1 = prod(i,k) + tke(i,k) / dtn
           diss(i,k)=max(min(diss(i,k), tem1), 0.)
           tke(i,k) = tke(i,k) + dtn * (prod(i,k)-diss(i,k))! no diffusion yet
!          tke(i,k) = max(tke(i,k), tkmin)
           tke(i,k) = max(tke(i,k), tkmnz(i,k))
        enddo
      enddo
      enddo
!
!     compute updraft & downdraft properties for tke
!
      do k = 1, km
        do i = 1, im
          if(pcnvflg(i)) then
! kgao change
!            qcko(i,k,ntke) = tke(i,k)
            qcko(i,k,ntrac) = tke(i,k)
          endif
          if(scuflg(i)) then
! kgao change
!            qcdo(i,k,ntke) = tke(i,k)
            qcdo(i,k,ntrac) = tke(i,k)
          endif
        enddo
      enddo
      do k = 2, kmpbl
        do i = 1, im
          if (pcnvflg(i) .and. k <= kpbl(i)) then
             dz   = zl(i,k) - zl(i,k-1)
             tem  = 0.5 * xlamue(i,k-1) * dz
             factor = 1. + tem
! kgao change
!             qcko(i,k,ntke)=((1.-tem)*qcko(i,k-1,ntke)+tem*
!     &                (tke(i,k)+tke(i,k-1)))/factor
             qcko(i,k,ntrac)=((1.-tem)*qcko(i,k-1,ntrac)+tem*
     &                (tke(i,k)+tke(i,k-1)))/factor
          endif
        enddo
      enddo
      do k = kmscu, 1, -1
        do i = 1, im
          if (scuflg(i) .and. k < krad(i)) then
            if(k >= mrad(i)) then
              dz = zl(i,k+1) - zl(i,k)
              tem  = 0.5 * xlamde(i,k) * dz
              factor = 1. + tem
! kgao change
!              qcdo(i,k,ntke)=((1.-tem)*qcdo(i,k+1,ntke)+tem*
!     &                 (tke(i,k)+tke(i,k+1)))/factor
              qcdo(i,k,ntrac)=((1.-tem)*qcdo(i,k+1,ntrac)+tem*
     &                 (tke(i,k)+tke(i,k+1)))/factor
            endif
          endif
        enddo
      enddo
!
!----------------------------------------------------------------------
!     compute tridiagonal matrix elements for turbulent kinetic energy
!
      do i=1,im
         ad(i,1) = 1.0
         f1(i,1) = tke(i,1)
      enddo
!
      do k = 1,km1
        do i=1,im
          dtodsd  = dt2/del(i,k)
          dtodsu  = dt2/del(i,k+1)
          dsig    = prsl(i,k)-prsl(i,k+1)
          rdz     = rdzt(i,k)
          tem1    = dsig * dkq(i,k) * rdz
          dsdz2   = tem1 * rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
          ad(i,k) = ad(i,k)-au(i,k)
          ad(i,k+1)= 1.-al(i,k)
          tem2    = dsig * rdz
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             ptem      = 0.5 * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             tem       = tke(i,k) + tke(i,k+1)
! kgao change
!             ptem      = qcko(i,k,ntke) + qcko(i,k+1,ntke)
             ptem      = qcko(i,k,ntrac) + qcko(i,k+1,ntrac) 
             f1(i,k)   = f1(i,k)-(ptem-tem)*ptem1
             f1(i,k+1) = tke(i,k+1)+(ptem-tem)*ptem2
          else
             f1(i,k+1) = tke(i,k+1)
          endif
!
          if(scuflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              ptem      = 0.5 * tem2 * xmfd(i,k)
              ptem1     = dtodsd * ptem
              ptem2     = dtodsu * ptem
              tem       = tke(i,k) + tke(i,k+1)
! kgao change
!              ptem      = qcdo(i,k,ntke) + qcdo(i,k+1,ntke)
              ptem      = qcdo(i,k,ntrac) + qcdo(i,k+1,ntrac)
              f1(i,k)   = f1(i,k) + (ptem - tem) * ptem1
              f1(i,k+1) = f1(i,k+1) - (ptem - tem) * ptem2
            endif
          endif
!
        enddo
      enddo
c
c     solve tridiagonal problem for tke
c
      call tridit(im,km,1,al,ad,au,f1,au,f1)
c
c     recover tendency of tke
c
      do k = 1,km
         do i = 1,im
! fix negative tke 
           f1(i,k) = max(f1(i,k), tkmin)
! kgao change
!            qtend = (f1(i,k)-q1(i,k,ntke))*rdt
!            rtg(i,k,ntke) = rtg(i,k,ntke)+qtend
            qtend = (f1(i,k)-q1(i,k,ntrac))*rdt
            rtg(i,k,ntrac) = rtg(i,k,ntrac)+qtend
         enddo
      enddo
c
c     compute tridiagonal matrix elements for heat and moisture (and other tracers, except tke)
c
      do i=1,im
         ad(i,1) = 1.
         f1(i,1) = t1(i,1)   + dtdz1(i) * heat(i)
         f2(i,1) = q1(i,1,1) + dtdz1(i) * evap(i)
      enddo
      if(ntrac1 >= 2) then
        do kk = 2, ntrac1
          is = (kk-1) * km
          do i = 1, im
            f2(i,1+is) = q1(i,1,kk)
          enddo
        enddo
      endif
c
      do k = 1,km1
        do i = 1,im
          dtodsd  = dt2/del(i,k)
          dtodsu  = dt2/del(i,k+1)
          dsig    = prsl(i,k)-prsl(i,k+1)
          rdz     = rdzt(i,k)
          tem1    = dsig * dkt(i,k) * rdz
          dsdzt   = tem1 * gocp
          dsdz2   = tem1 * rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
          ad(i,k) = ad(i,k)-au(i,k)
          ad(i,k+1)= 1.-al(i,k)
          tem2    = dsig * rdz
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             ptem      = 0.5 * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             tem       = t1(i,k) + t1(i,k+1)
             ptem      = tcko(i,k) + tcko(i,k+1)
             f1(i,k)   = f1(i,k)+dtodsd*dsdzt -(ptem-tem)*ptem1
             f1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt +(ptem-tem)*ptem2
             ! kgao - updraft mass flux
             flux_up(i,k) = xmf(i,k) !0.5*(ptem-tem)*xmf(i,k)

             tem       = q1(i,k,1) + q1(i,k+1,1)
             ptem      = qcko(i,k,1) + qcko(i,k+1,1)
             f2(i,k)   = f2(i,k) - (ptem - tem) * ptem1
             f2(i,k+1) = q1(i,k+1,1) + (ptem - tem) * ptem2
          else
             f1(i,k)   = f1(i,k)+dtodsd*dsdzt
             f1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt
             f2(i,k+1) = q1(i,k+1,1)
          endif
!
          if(scuflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              ptem      = 0.5 * tem2 * xmfd(i,k)
              ptem1     = dtodsd * ptem
              ptem2     = dtodsu * ptem
              ptem      = tcdo(i,k) + tcdo(i,k+1)
              tem       = t1(i,k) + t1(i,k+1)
              f1(i,k)   = f1(i,k) + (ptem - tem) * ptem1
              f1(i,k+1) = f1(i,k+1) - (ptem - tem) * ptem2
              ! kgao - downdraft mass flux
              flux_dn(i,k) = xmfd(i,k) !-0.5*(ptem-tem)*xmfd(i,k)

              tem       = q1(i,k,1) + q1(i,k+1,1)
              ptem      = qcdo(i,k,1) + qcdo(i,k+1,1)
              f2(i,k)   = f2(i,k) + (ptem - tem) * ptem1
              f2(i,k+1) = f2(i,k+1) - (ptem - tem) * ptem2
            endif
          endif
        enddo
      enddo
!
      if(ntrac1 >= 2) then
        do kk = 2, ntrac1
          is = (kk-1) * km
          do k = 1, km1
            do i = 1, im
              if(pcnvflg(i) .and. k < kpbl(i)) then
                dtodsd = dt2/del(i,k)
                dtodsu = dt2/del(i,k+1)
                dsig  = prsl(i,k)-prsl(i,k+1)
                tem   = dsig * rdzt(i,k)
                ptem  = 0.5 * tem * xmf(i,k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                tem1  = qcko(i,k,kk) + qcko(i,k+1,kk)
                tem2  = q1(i,k,kk) + q1(i,k+1,kk)
                ! kgao note - turn off non-local mixing 
                f2(i,k+is) = f2(i,k+is) !- (tem1 - tem2) * ptem1
                f2(i,k+1+is)= q1(i,k+1,kk) !+ (tem1 - tem2) * ptem2
              else
                f2(i,k+1+is) = q1(i,k+1,kk)
              endif
!
              if(scuflg(i)) then
                if(k >= mrad(i) .and. k < krad(i)) then
                  dtodsd = dt2/del(i,k)
                  dtodsu = dt2/del(i,k+1)
                  dsig  = prsl(i,k)-prsl(i,k+1)
                  tem   = dsig * rdzt(i,k)
                  ptem  = 0.5 * tem * xmfd(i,k)
                  ptem1 = dtodsd * ptem
                  ptem2 = dtodsu * ptem
                  tem1  = qcdo(i,k,kk) + qcdo(i,k+1,kk)
                  tem2  = q1(i,k,kk) + q1(i,k+1,kk)
                  ! kgao note - turn off non-local mixing 
                  f2(i,k+is)  = f2(i,k+is) !+ (tem1 - tem2) * ptem1
                  f2(i,k+1+is)= f2(i,k+1+is) !- (tem1 - tem2) * ptem2
                endif
              endif
!
            enddo
          enddo
        enddo
      endif


      
c
c     solve tridiagonal problem for heat and moisture

c     Original:
c     call tridin(im,km,ntrac1,al,ad,au,f1,f2,au,f1,f2)

c     Modified: downward sweep only
      call tridin_down(im,km,ntrac1,al,ad,au,f1,f2,au,f1,f2)
      
c     Save upper diag and RHS for upward sweep 
      au_out(:,:)=au(:,:) 
      f1_out(:,:)=f1(:,:)
      f2_out(:,:)=f2(:,:)
      diss_out(:,:)=diss(:,:)

      ! to be rearranged in the upward sweep
      rtg_in(:,:,:)=rtg(:,:,:)


! Commented out, to be done in the upward sweep
! after updated all levels
!
! kgao note - rearrange tracer tendencies 
!
      !if(ntrac >= 3 ) then 
!        if(ntke == ntrac) then ! tke is the last tracer
!          rtg_in(:,:,:) = rtg(:,:,:) 
!        else                   ! tke is not
!          do kk = 1, ntke-1
!             rtg_in(:,:,kk) = rtg(:,:,kk)
!          enddo
!          rtg_in(:,:,ntke) = rtg(:,:,ntrac)
!          do kk = ntke+1, ntrac
!             rtg_in(:,:,kk) = rtg(:,:,kk-1)
!          enddo
!        endif
      !endif
!
!     add tke dissipative heating to temperature tendency
!
!      if(dspheat) then
!      do k = 1,km1
!        do i = 1,im
!!         tem = min(diss(i,k), dspmax)
!!         ttend = tem / cp
!          ttend = diss(i,k) / cp
!          tdt(i,k) = tdt(i,k) + dspfac * ttend
!        enddo
!      enddo
!      endif
c
c     compute tridiagonal matrix elements for momentum
c
      do i=1,im
         ad(i,1) = 1.0 + dtdz1(i) * stress(i) / spd1(i)
         f1(i,1) = u1(i,1)
         f2(i,1) = v1(i,1)
      enddo
c
      do k = 1,km1
        do i=1,im
          dtodsd  = dt2/del(i,k)
          dtodsu  = dt2/del(i,k+1)
          dsig    = prsl(i,k)-prsl(i,k+1)
          rdz     = rdzt(i,k)
          tem1    = dsig * dku(i,k) * rdz
          dsdz2   = tem1*rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
          ad(i,k) = ad(i,k)-au(i,k)
          ad(i,k+1)= 1.-al(i,k)
          tem2    = dsig * rdz
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             ptem      = 0.5 * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             tem       = u1(i,k) + u1(i,k+1)
             ptem      = ucko(i,k) + ucko(i,k+1)
             f1(i,k)   = f1(i,k) - (ptem - tem) * ptem1
             f1(i,k+1) = u1(i,k+1) + (ptem - tem) * ptem2
             tem       = v1(i,k) + v1(i,k+1)
             ptem      = vcko(i,k) + vcko(i,k+1)
             f2(i,k)   = f2(i,k) - (ptem - tem) * ptem1
             f2(i,k+1) = v1(i,k+1) + (ptem - tem) * ptem2
          else
             f1(i,k+1) = u1(i,k+1)
             f2(i,k+1) = v1(i,k+1)
          endif
!
          if(scuflg(i)) then
            if(k >= mrad(i) .and. k < krad(i)) then
              ptem      = 0.5 * tem2 * xmfd(i,k)
              ptem1     = dtodsd * ptem
              ptem2     = dtodsu * ptem
              tem       = u1(i,k) + u1(i,k+1)
              ptem      = ucdo(i,k) + ucdo(i,k+1)
              f1(i,k)   = f1(i,k) + (ptem - tem) *ptem1
              f1(i,k+1) = f1(i,k+1) - (ptem - tem) *ptem2
              tem       = v1(i,k) + v1(i,k+1)
              ptem      = vcdo(i,k) + vcdo(i,k+1)
              f2(i,k)   = f2(i,k) + (ptem - tem) * ptem1
              f2(i,k+1) = f2(i,k+1) - (ptem - tem) * ptem2
            endif
          endif
!
        enddo
      enddo
c
c     solve tridiagonal problem for momentum
c
      call tridi2(im,km,al,ad,au,f1,f2,au,f1,f2)
c
c     recover tendencies of momentum
c
      do k = 1,km
         do i = 1,im
            utend = (f1(i,k)-u1(i,k))*rdt
            vtend = (f2(i,k)-v1(i,k))*rdt
            du(i,k)  = du(i,k)+utend
            dv(i,k)  = dv(i,k)+vtend
            dusfc(i) = dusfc(i)+conw*del(i,k)*utend
            dvsfc(i) = dvsfc(i)+conw*del(i,k)*vtend
         enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  pbl height for diagnostic purpose
!
      do i = 1, im
         hpbl(i) = hpblx(i)
         kpbl(i) = kpblx(i)
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
      end



! the upward sweep for heat, moisture and tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine satmedmfvdifq_up(ix,im,km,ntrac,ntke,
     &     tdt,rtg_in,q1,t1,dspheat,dspfac,
     &     del,delt,
     &     dtsfc,dqsfc,
     &     au_in,f1_in,f2_in,diss_in)
!
      use machine  , only : kind_phys
      use physcons, grav => con_g, rd => con_rd, cp => con_cp
     &,             rv => con_rv, hvap => con_hvap
!
      implicit none
!
      logical dspheat
      integer ix, im, km, km1, ntrac,ntke
      real(kind=kind_phys) delt, dspfac
      real(kind=kind_phys) tdt(im,km), rtg_in(im,km,ntrac)
      real(kind=kind_phys) q1(ix,km,ntrac), t1(ix,km)
      real(kind=kind_phys) del(ix,km)
      real(kind=kind_phys) dtsfc(im), dqsfc(im)
      
      real(kind=kind_phys) rtg(im,km,ntrac)
      ! Input arrays from part 1
      real(kind=kind_phys) au_in(im,km-1), diss_in(im,km-1)
      real(kind=kind_phys) f1_in(im,km), f2_in(im,km*(ntrac-1))
!
      ! Local variables
      integer i, k, kk, is, ntrac1
      real(kind=kind_phys) rdt, ttend, qtend
      real(kind=kind_phys) cont, conq
!
      parameter(cont=cp/grav,conq=hvap/grav)
!
      rdt = 1. / delt
      km1=km-1
      ntrac1= ntrac - 1
      rtg=rtg_in

      ! Perform upward sweep for heat, moisture and tracers
      call tridin_up(im,km,ntrac1,au_in,f1_in,f2_in)

      ! Apply tendencies for heat and moisture
      do k = 1,km
         do i = 1,im
            ttend      = (f1_in(i,k)-t1(i,k))*rdt
            qtend      = (f2_in(i,k)-q1(i,k,1))*rdt
            tdt(i,k)   = tdt(i,k)+ttend
            rtg(i,k,1) = rtg(i,k,1)+qtend
            dtsfc(i)   = dtsfc(i)+cont*del(i,k)*ttend
            dqsfc(i)   = dqsfc(i)+conq*del(i,k)*qtend
         enddo
      enddo

      ! Apply tendencies for other tracers
      if(ntrac1 >= 2) then
        do kk = 2, ntrac1
          is = (kk-1) * km
          do k = 1, km
            do i = 1, im
              qtend = (f2_in(i,k+is)-q1(i,k,kk))*rdt
              rtg(i,k,kk) = rtg(i,k,kk)+qtend
            enddo
          enddo
        enddo
      endif

!!!!! moved from _down - Joseph
!!!!!! kgao note - rearrange tracer tendencies 
      !if(ntrac >= 3 ) then 
        if(ntke == ntrac) then ! tke is the last tracer
          rtg_in(:,:,:) = rtg(:,:,:)
        else                   ! tke is not
          do kk = 1, ntke-1
             rtg_in(:,:,kk) = rtg(:,:,kk)
          enddo
          rtg_in(:,:,ntke) = rtg(:,:,ntrac)
          do kk = ntke+1, ntrac
             rtg_in(:,:,kk) = rtg(:,:,kk-1)
          enddo
        endif
      !endif
!
!     add tke dissipative heating to temperature tendency
!
      if(dspheat) then
      do k = 1,km1
        do i = 1,im
!         tem = min(diss(i,k), dspmax)
!         ttend = tem / cp
          ttend = diss_in(i,k) / cp
          tdt(i,k) = tdt(i,k) + dspfac * ttend
        enddo
      enddo
      endif

      return
      end
