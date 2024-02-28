!>  \file som_mlm.f90
!!  This file contains routines for a Slab Ocean Model (SOM) and
!!  also a Mixed Layer Ocean Model (MLM)
!
!!  Contacted Baoqiang Xiang at baoqiang.xiang@noaa.gov

!  ==========================================================  !!!!!
!                          'module_ocean' description          !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!    this module sets up SST using a slab ocean model (SOM) or         !
!    mixed layer ocean model (MLM)                                                                  !
!    in the module, the externally callabe subroutines are :           !
!                                                                      !
!      'ocean_init'   -- initialization SOM  by setting some namelists   !
!                                                                      !
!      'update_ocean' -- update SST with the combined effect of net      !
!                      surface heat flux and the nudging term          !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!


!========================================!
      module module_ocean                !
!........................................!
!
      use physcons,          only : con_tice
      use physparam,         only : kind_phys
      use GFS_typedefs,      only : GFS_control_type, GFS_grid_type
!      use constants_mod,     only : omega, grav
!
      implicit   none
      private

      public  ocean_init, update_ocean
!
      real (kind=kind_phys)   :: width_buffer, minmld,  &
                                 cpwater, rhowater, omega, grav
      parameter(minmld       = 10.)  ! minimum mixed layer depth
      parameter(width_buffer = 15.)  ! the width of a buffer band where SST is determined by both SOM/MLM
                                     ! and climatology (or climatology plus initial anomaly)
      parameter(cpwater      = 4000.)
      parameter(rhowater     = 1000.)
      parameter(omega        = 7.292e-5)
      parameter(grav         = 9.80)

!    namelist variables
      character(len=24)     :: ocean_option       = 'SOM'     ! option to set ocean mixed layer depth (MLD)
                                                              ! using either 'SOM' or 'MLM'
      character(len=24)     :: mld_option         = 'obs'     ! option to set ocean mixed layer depth (MLD)
                                                              ! using either 'obs' or 'const'
      real(kind=kind_phys)  :: mld_obs_ratio      = 1.        ! tuning parameter for observed MLD
      real(kind=kind_phys)  :: stress_ratio       = 1.        ! how much of wind stress is applied in the mixed layer
      integer               :: restore_method     = 1         ! option 1: nudging toward observational climatology
                                                              ! option 2: nudging toward observational climatology plus
                                                              !           initial anomaly with a decay time scale of FTSFS (90 days)
                                                              ! option 3: nudging toward observed SST
      logical               :: use_old_mlm        = .false.   ! if true: very similar to WRF model
      logical               :: use_rain_flux      = .false.   ! considering the rainfall induced surface flux
      logical               :: use_qflux          = .false.   ! considering the qflux correction
      logical               :: do_mld_restore     = .false.   ! restoring MLD toward observed climatology
      real(kind=kind_phys)  :: const_mld          = 40.       ! constant ocean MLD (meter)
      real(kind=kind_phys)  :: Gam                = 0.14      ! ocean temp lapese rate at the bottom of MLD (degree per m)
      real(kind=kind_phys)  :: eps_day            = 10.       ! damping time scale of ocean current (days)
      real(kind=kind_phys)  :: sst_restore_tscale = 3.        ! restoring time scale for sst (day)
      real(kind=kind_phys)  :: mld_restore_tscale = 1.        ! restoring time scale for mld (day)
      real(kind=kind_phys)  :: start_lat          = -30.      ! latitude starting from? Note that this value should not be smaller than -maxlat.
      real(kind=kind_phys)  :: end_lat            = 30.       ! latitude ending with? Note that this value should not be bigger than maxlat.
      real(kind=kind_phys)  :: tday1              = 3.        !
      real(kind=kind_phys)  :: tday2              = 10.       !
      real(kind=kind_phys)  :: sst_restore_tscale1= 3.        ! restoring time scale for sst during the period from 1 to tday1
      real(kind=kind_phys)  :: sst_restore_tscale2= 10.       ! restoring time scale for sst for the period beyond tday2
      real(kind=kind_phys)  :: mld_restore_tscale1= 3.        ! restoring time scale for mld during the period from 1 to tday1
      real(kind=kind_phys)  :: mld_restore_tscale2= 10.       ! restoring time scale for mld for the period beyond tday2
                                                              ! beyond the latitude bands (start_lat:end_lat), using climatological SST or
                                                              ! climatological SST plus initial anomaly
      logical               :: use_tvar_restore_sst  = .false.! using time varying restoring time scale for sst
      logical               :: use_tvar_restore_mld  = .false.! using time varying restoring time scale for mld
      real(kind=kind_phys)  :: maxlat = 60.                   ! maximum latitudinal extent of the SOM/MLM. Generally set to <= 60, though
                                                              ! it can be useful to set to > 60 if the desired start_lat and end_lat are
                                                              ! poleward of 60 degrees. If set to > 60, it is recommended to also set
                                                              ! gfs_physics_nml.disable_radiation_quasi_sea_ice to .true. to prevent
                                                              ! an unphysical quasi-ice-albedo feedback from occurring. Set to 90 along
                                                              ! with ocean_nml.start_lat = -90 and ocean_nml.end_lat = 90 to enable
                                                              ! running with a global interactive ocean.

      namelist /ocean_nml/   &
       ocean_option, mld_option, mld_obs_ratio, stress_ratio, restore_method,  &
       use_old_mlm, use_rain_flux, use_qflux, do_mld_restore, const_mld, Gam,  &
       eps_day, sst_restore_tscale, mld_restore_tscale, start_lat, end_lat,    &
       tday1, tday2, sst_restore_tscale1, sst_restore_tscale2, mld_restore_tscale1, &
       mld_restore_tscale2, use_tvar_restore_sst, use_tvar_restore_mld, maxlat

! =================
      contains
! =================


!-----------------------------------
      subroutine ocean_init                                               &
     &     ( Model, logunit, input_nml_file )!  ---  inputs:

!  ===================================================================  !
!                                                                       !
!  this program is the initialization program for SOM/MLM model             !
!                                                                       !
! usage:         call ocean_init                                          !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!
      implicit none

!  ---  inputs:
      type (GFS_control_type),    intent(in)  :: Model
      integer,    intent(in)  :: logunit
      character (len = *), optional, intent (in) :: input_nml_file(:)

!  ---  outputs: ( none )

!  ---  locals:
      integer    :: ios
      logical    :: exists
!
!===> ...  begin here
!
!--- read in the namelist
!      inquire (file=trim(Model%fn_nml), exist=exists)
!      if (.not. exists) then
!      write(6,*) 'ocean_namelist_read:: namelist file: ',trim(Model%fn_nml),' does not exist'
!      stop
!      else
!      open (unit=Model%nlunit, file=Model%fn_nml, READONLY, status='OLD', iostat=ios)
!      endif
!      rewind(Model%nlunit)
!      read (Model%nlunit, nml=ocean_nml)
!      close (Model%nlunit)

#ifdef INTERNAL_FILE_NML
        read(input_nml_file, nml=ocean_nml)
#else
!       print *,' in sfcsub nlunit=',nlunit,' me=',me,' ialb=',ialb
       inquire (file=trim(Model%fn_nml), exist=exists)
       if (.not. exists) then
        write(6,*) 'ocean_namelist_read:: namelist file: ',trim(Model%fn_nml),' does not exist'
        stop
       else
        open (unit=Model%nlunit, file=Model%fn_nml, READONLY, status='OLD', iostat=ios)
       endif
       rewind(Model%nlunit)
       read (Model%nlunit,ocean_nml)
       close (Model%nlunit)
#endif

      if (start_lat < -maxlat) then
       write(*,*) 'start_lat should not be smaller than', -maxlat
       call abort
      endif
      if (end_lat > maxlat) then
       write(*,*) 'end_lat should not be larger than', maxlat
       call abort
      endif

      if (restore_method == 3 .and. .not. Model%use_ext_sst) then
         write(6,*) ' som_mlm::ocean_init(): Cannot use restore_method == 3'
         write(6,*) '                        unless external SST provided '
         write(6,*) '                        (use_ext_sst = .true.). Stop.'
         call abort
      endif

!--- write namelist to log file ---
      if (Model%me == Model%master) then
       write(logunit, *) "============================================="
       write(logunit, *) "Slab (or Mixed Layer) Ocean Model"
       write(logunit, nml=ocean_nml)
      endif
!
      return
!...................................
      end subroutine ocean_init
!-----------------------------------
!
      subroutine update_ocean                                                      &
       (im, dtp, Grid, islmsk, kdt, kdt_prev, netflxsfc, taux, tauy, rain, tair,   &
        qflux_restore, qflux_adj, mldclim, tsclim, ts_clim_iano, ts_obs, ts_som,   &
        tsfc, tml, tml0, mld, mld0, huml, hvml, tmoml, tmoml0, iau_offset)  
    
!  ===================================================================  !
!                                                                       !
!  this program computes the updated SST based on a simple SOM/MLM model    !
!  Within start_lat - end_lat, SST is determined by SOM/MLM
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer, intent(in)                              :: im, kdt, kdt_prev, iau_offset
      real,    intent(in)                              :: dtp  ! model time step
      type (GFS_grid_type), intent(in)                 :: Grid
      integer, dimension(:), intent(in)                :: islmsk
      real (kind=kind_phys), dimension(:), intent(in)  ::    &
           netflxsfc,     & ! net surface heat flux
           taux,          &
           tauy,          &
           rain,          &
           tair,          & ! lowest model level temp
           qflux_adj,     & ! lowest model level temp
           mldclim,       & ! ocean MLD
           tsclim,        & ! observed climatological SST
           ts_clim_iano,  & ! observed climatological SST plus initial anomaly
           ts_obs           ! observed SST (for simulation)

!  ---  inoutputs
      real (kind=kind_phys), dimension(:), intent(inout) ::   &
           ts_som,        &
           tsfc,          &
           tml,           &
           tml0,          &
           mld,           &
           mld0,          &
           huml,          &
           hvml,          &
           tmoml,         &
           tmoml0

      real (kind=kind_phys), dimension(:), intent(out)   ::   &
           qflux_restore  ! restoring flux for diagnosis purpose

!  ---  locals:
      real (kind=kind_phys)                              ::   &
           lat, mlcp, mldc, taut, taum, &
           alphat,alpham, bufzs,        &
           bufzn, fcor, c1, c2, r1, r2
      real (kind=kind_phys), dimension (size(tsfc,1))    :: tsfc1, tsfc2
      real (kind=kind_phys), dimension (size(tsfc,1))    :: qsfc
      integer :: i
      real (kind=kind_phys)                              ::   &
           tmlp, mldp, humlp, hvmlp, mldn, tmln, tmomln, fday, tem
!
!===> ...  begin here
!
      if (iau_offset > 0 .and. kdt_prev > 0) then
        fday = (kdt - kdt_prev) * dtp / 86400.
      else
        fday = kdt * dtp / 86400.
      endif
!
      qsfc = 0.
      if (use_tvar_restore_sst) then
       if (fday < tday1) then
        taut = sst_restore_tscale1*86400.
       elseif (fday >= tday1 .and. fday < tday2 ) then
        tem = (sst_restore_tscale2 - sst_restore_tscale1)/(tday2-tday1)
        taut = (tem*(fday-tday1) + sst_restore_tscale1) *86400.
       else
        taut = sst_restore_tscale2*86400.
       endif
      else
       taut = sst_restore_tscale*86400.
      endif
      alphat = 1. + dtp/taut
!
      if (use_tvar_restore_mld) then
       if (fday < tday1) then
        taum = mld_restore_tscale1*86400.
       elseif (fday >= tday1 .and. fday < tday2 ) then
        tem = (mld_restore_tscale2 - mld_restore_tscale1)/(tday2-tday1)
        taum = (tem*(fday-tday1) + mld_restore_tscale1) *86400.
       else
        taum = mld_restore_tscale2*86400.
       endif
      else
       taum = mld_restore_tscale*86400.
      endif
      alpham = 1. + dtp/taum

!    two buffer zones:
!    (bufzs - start_lat)  and (end_lat - bufzn)
      bufzs = max(-maxlat - 0.0001, start_lat - width_buffer)  
      bufzn = min( maxlat + 0.0001, end_lat + width_buffer)     
!      
      if (kdt == 1 .or. (iau_offset > 0 .and. kdt-kdt_prev == 1)) then
       do i=1, im
        ts_som(i)=tsfc(i)
       enddo
       if (ocean_option == "MLM") then
        do i=1, im
        tml(i)=tsfc(i)
        tml0(i)=tsfc(i)
        enddo
        do i=1, im
        mld(i)=mldclim(i)
        mld0(i)=mld(i)
        huml(i)=0.
        hvml(i)=0.
        tmoml(i)=tsfc(i)-5.
        tmoml0(i)=tmoml(i)
        enddo
       endif
      endif

      ! Reset the slab ocean temperature or mixed layer ocean properties in
      ! any sea ice grid cells.  This is required in the (rare) instances
      ! when one is running with the SOM/MLM with a latitudinal extent such
      ! that grid cells transition from ocean to sea ice and back.
      where (islmsk .eq. 2)
       ts_som = con_tice
      endwhere

      if (ocean_option == "MLM") then
       where(islmsk .eq. 2)
        tml = con_tice
        tml0 = con_tice
        mld = mldclim
        mld0 = mld
        huml = 0.
        hvml = 0.
        tmoml = con_tice - 5.
        tmoml0 = tmoml
       endwhere
      endif
!
      if (use_rain_flux) then
       do i=1,im
        if ( islmsk(i) ==0 ) then
         qsfc(i) = netflxsfc (i) + cpwater* rain(i)/dtp*( tair(i)-tsfc(i) )
        endif
       enddo
      else
       qsfc = netflxsfc
      endif

      if (use_qflux) then
       do i=1,im
         qsfc(i) = qsfc (i) + qflux_adj (i)
       enddo
      endif

      do i = 1, im
!
       if (mld_option == 'const') then
        mlcp = const_mld * rhowater * cpwater    ! rho*Cp*mld
        mldc = const_mld
       elseif (mld_option == 'obs') then
        mlcp =  max(minmld, mld_obs_ratio* mldclim(i)) * rhowater *cpwater   ! rho*Cp*mld
        mldc =  mld_obs_ratio* mldclim(i)
       else
        write(*,*) ' mld_option can only be const or obs now'
        call abort
       endif

       fcor = 2 * omega * sin (Grid%xlat(i))

       if ( islmsk(i) ==0 ) then
        if (ocean_option == "SOM") then
         mld(i)   =  mldc
        elseif (ocean_option == "MLM") then
         tmlp    =  tml(i)
         mldp    =  mld(i)
         humlp   =  huml(i)
         hvmlp   =  hvml(i)
         tmln    =  tml0(i)
!
!         tmomln  =  tmoml0(i)
         tmomln  =  tmoml(i)
!
         mldn    =  mld0(i)
         call MLM1D(dtp, fcor, taum, alpham, qsfc(i), taux(i), tauy(i),     &
             tmlp, tmln, tmomln, mldp, mldn, mldc, humlp, hvmlp)
        endif !end ocean_option

        select case (restore_method)
        case(1)
           tsfc2(i) = tsclim(i)
        case(2)
           tsfc2(i) = ts_clim_iano(i)
        case (3)
           tsfc2(i) = ts_obs(i)
        case default
           !call mpp_error(FATAL, 'restore_method not implemented')
           print*, 'restore_method = ', restore_method, ' not implemented'
           stop 121
        end select

        select case (ocean_option)
        case("SOM")
          if (use_qflux) then
           tsfc1(i) = ts_som(i) + qsfc(i)/mlcp*dtp
          else
           tsfc1(i) = (ts_som(i) + qsfc(i)/mlcp*dtp + tsfc2(i)/taut*dtp ) / alphat
          endif
        case("MLM")
          tsfc1(i) = (tmlp + tsfc2(i)/taut*dtp)/alphat
          tml(i)   =  tsfc1(i)
          mld(i)   =  mldp
          huml(i)  =  humlp
          hvml(i)  =  hvmlp
          tmoml(i) =  tml(i) - 5. ! not used
        case default
           !call mpp_error(FATAL, "ocean_option must be SOM or MLM")
           print*, 'ocean_option must be SOM or MLM; ocean_option set to ', ocean_option
           stop 122
        end select
        qflux_restore(i) = (ts_clim_iano(i) - tsfc1(i)) * mlcp / taut  ! for diagnosis purpose only

        ts_som (i) = tsfc1 (i)
       endif  ! end islmsk
      enddo

      do i = 1, im
       if (islmsk(i) == 0 ) then
        lat = Grid%xlat(i) * 57.29578
        c1 = min(1.0, abs((lat -bufzs) / (start_lat-bufzs)) )
!        r1 = (exp(c1**interp_order)-1.)/(exp(1.0)-1.)
        c2 = min(1.0, abs((bufzn - lat) / (bufzn - end_lat)) )
!        r2 = (exp(c2**interp_order)-1.)/(exp(1.0)-1.)
        if (lat >= start_lat .and. lat<= end_lat ) then
         tsfc(i) = tsfc1(i)
        elseif (lat >= bufzs .and. lat < start_lat) then  ! the first buffer zone
!         tsfc(i) = c1 * tsfc1(i) + (1.-c1) * tsfc2(i)
         tsfc(i) = c1 * ts_som(i) + (1.-c1) * tsfc2(i)
        elseif (lat > end_lat .and. lat <= bufzn) then    ! the second buffer zone
!         tsfc(i) = c2 * tsfc1(i) + (1.-c2) * tsfc2(i)
         tsfc(i) = c2 * ts_som(i) + (1.-c2) * tsfc2(i)
        else
         tsfc(i) = tsfc2(i)
        endif
       endif !end islmsk
      enddo
!
      return
!...................................
      end subroutine update_ocean
!-----------------------------------

      subroutine MLM1D(dt, F, taum, alpham, qsfc, taux, tauy,           &
                       tml, tml0, tmoml, H, H0, HC, huml, hvml)
!----------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------
!
!  SUBROUTINE OCEANML CALCULATES THE SEA SURFACE TEMPERATURE
!  FROM A SIMPLE OCEAN MIXED LAYER MODEL BASED ON
!  (Pollard, Rhines and Thompson (1973).
!
!-- DT          time step (second)
!-- F           Coriolis parameter
!-- taum        MLD restoring time scale
!-- alpham      MLD restoring parameter
!-- qsfc        net surface heat flux
!-- taux        wind stress at zonal direction
!-- tauy        wind stress at meridional direction
!-- tml         ocean mixed layer temperature (K)
!-- tml0        ocean mixed layer temperature (K) at initial time or previous time step
!-- tmoml       top 200 m ocean mean temperature (K) at initial time or previous time step
!-- H           ocean mixed layer depth (m)
!-- H0          ocean mixed layer depth (m) at initial time or nudged MLD toward climatology
!-- HC          climatological or constant ocean mixed layer depth (m)
!-- huml        ocean mixed layer u component of wind
!-- hvml        ocean mixed layer v component of wind
!
!  Note: Part of the code for this subroutine is from WRF model
!----------------------------------------------------------------

      REAL,    INTENT(INOUT)    :: tml, H, huml, hvml

      REAL,    INTENT(IN   )    :: dt, F, taum, alpham, qsfc, taux, tauy,      &
                                   tml0, tmoml, H0, HC
! Local
      REAL ::  alp, BV2, A1, A2, A3, B2, u, v,  &
           hu1, hv1, hu2, hv2, q, hold,         &
           hsqrd, thp, taux2, tauy2, fdt, damp

      hu1=huml
      hv1=hvml
      fdt = f * dt
!
      alp=max((tml-273.15)*1.e-5, 1.e-6)
      BV2=alp*grav*Gam
      thp=tml0-Gam*(h-h0)
      if (use_old_mlm) then
       A1=(tml-thp)*h - 0.5*Gam*h*h
      else
       A1=(tml-tml0)*h + 0.5*Gam*(h-h0)*abs(h-h0)
      endif
      if(h.ne.0.)then
        u=hu1/h
        v=hv1/h
      else
        u=0.
        v=0.
      endif

! determine how much of wind stress is applied to mixed layer
       taux2 = taux * stress_ratio
       tauy2 = tauy * stress_ratio
       q=qsfc/(rhowater*cpwater)
! note: forward-backward coriolis force for effective time-centering
       if (use_old_mlm) then
!        hu2=hu1+dt*( f*hv1 + taux2/rhowater - damp*hu1)
!        hv2=hv1+dt*(-f*hu2 + tauy2/rhowater - damp*hv1)
         damp = 1. / 86400./eps_day
        hu2=( hu1+dt*( f*hv1 + taux2/rhowater ) )/(1.0+damp*dt)
        hv2=( hv1+dt*(-f*hu2 + tauy2/rhowater ) )/(1.0+damp*dt)
       else
        hu2=( (1-fdt**2/4.)*hu1+fdt*hv1+taux2/rhowater*dt+f*dt**2/2./rhowater*tauy2 ) / &
             (1.+fdt**2/4.)
        hv2=hv1+tauy2/rhowater*dt-fdt/2.*(hu2+hu1)
       endif
! consider the flux effect
       A2 = A1+q*dt
       A3 = A1+q*dt - 0.5*Gam*h0**2

       huml=hu2
       hvml=hv2

       hold=h
       B2=hu2*hu2+hv2*hv2
       if (use_old_mlm) then
        hsqrd=-A2/Gam + sqrt(A2*A2/(Gam*Gam) + 2.*B2/BV2)
       else
        hsqrd=-A3/Gam + sqrt(A3*A3/(Gam*Gam) + 2.*B2/BV2)
       endif
       h=sqrt(max(hsqrd,0.0))
       h=min(h, 500.0)

!       write(0,*) 'test0',h,hc,taum,alpham,dt
       if(do_mld_restore) then
         h  = (h + HC/taum*dt)/alpham
       endif
!       write(0,*) 'test1',h,hc,taum,alpham,dt
! limit to posit ive h change
!       if (use_old_mlm) then
!        if(h.lt.hold) h=hold
!       else
!        if(h.lt.hold) h=h0
!       endif
! no change unless tml is warmer than layer mean temp tmol or tsk-5 (see omlinit)
       if(tml.ge.tmoml .and. h.ne.0.)then

! no change unless tml is warmer than layer mean temp tmoml or tsk-5 (see omlinit)
         if(tml.ge.tmoml)then
          if (use_old_mlm) then
! if MLD does not deepen, we only consider the surface heat flux effect
           if (h <= hold) then
            tml=max(tml + q*dt/h, tmoml)
           else
            tml=max(tml0 - Gam*(h-h0) + 0.5*Gam*h + A2/h, tmoml)
           endif
          else
           tml=max(tml0 -0.5* Gam*(h-h0)*abs(h-h0)/h + A2/h, tmoml)
          endif
         else
           tml=tmoml
         endif
         u=hu2/h
         v=hv2/h
       else
         tml=tml0
         u=0.
         v=0.
       endif
!
       tml = max (273.15, tml)

      end subroutine MLM1D

      end module module_ocean

!=========================================
