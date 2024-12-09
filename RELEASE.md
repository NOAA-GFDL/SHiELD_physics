# RELEASE NOTES for FV3 202411: Summary
FV3-202411-public --- November 2024
Lucas Harris, GFDL lucas.harris@noaa.gov

This version has been tested with:  
FV3 Dynamical Core release FV3-202411-public from https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere  
FMS release 2024.03 from https://github.com/NOAA-GFDL/FMS  
FMS Coupler release 2024.03.01 from https://github.com/NOAA-GFDL/FMScoupler  
Atmos Drivers release FV3-202411-public from https://github.com/NOAA-GFDL/atmos_drivers  

This release includes the following:
- Support for running SHiELD in the full FMSCoupler (Joseph)
- GFDL Microphysics (Linjiong)
  - Replaced the hardcoded qa with aerosol data
  - Replaced the hardcoded hydrostatic with the input hydrostatic variable
  - Added two namelist options (fast_fr_mlt and fast_dep_sub) to control freezing/melting and deposition/sublimation in the fast microphysics.
  - Included a missing term in the energy conservation formula (credit: Tristan Abbott).  May affect  prediction of processes depending strongly on microphysics. Compile the model with -DENG_CNV_OLD to revert this change.
  - Added diagnostic cloud content and cloud effective radii of all cloud hydrometeors (qc*, re*).
  - Added diagnostic microphysical process rate (mpp*).
  - Removed unused Keihl et al. (1994) cloud water effective radius diagnosis.
  - A new COSP output variable named clisccp was added.
- Surface processes and PBL (Kun, Linjiong, Kai, Spencer)
  - Brought EMC updates from dev/emc branch related to TKE-EDMF, deep and shallow convection schemes
  - Added a constraint to the PBL bottom layer to prevent negative water vapor. 
  - Revise the surface process to prevent excessive downward latent heat transfer 
  - Exposed additional tunable parameters in TKE-EDMF PBL scheme
  - Fix miscalculated albedo over the ocean when ialb = 2
  - Slab ocean bug fixes and enhancements
  - Upgraded land-ice related processes in Noah-MP based on version 4.5 of the official Noah-MP repository and some updates from NCAR/ccpp-physics, mostly from EMC.
    - Computing the snow height within the snowwater_glacier and snowwater subroutines
    - Limiting the compaction of snow in various places
    - Calculating the equilibrium state of snow only if snow exists. The maximum water depth becomes a tunable parameter and its default value remains 5 m.
    - Prescribed soil color capability
    - Prescribed snow albedo capability
- Convection schemes (Kun, Linjiong)
  - Brought EMC updates from dev/emc branch
  - Merged samfshalcnv_gfdl.f and samfshalcnv.f.  Use gfs_physics_nml::limit_shal_conv = .true. to recover samfshalcnv_gfdl. 
  - Created the shallow convection option 5 based on option 3, to turn off the shallow convection scheme if the diagnosed cloud depth or cloud top exceeds certain critical values (set by cthk_shal and top_shal, respectively)
  - Added a parameter dxcrtas to the deep convection; for grid-spacing less than this (default 8000 m) the AS's quasi-equilibrium assumption is considered invalid.
- Radiation (Kai)
  - Added a namelist option, gfs_physics_nml::fco2_scaling to uniformly scale CO2
- Diagnostics (Spencer)
  - Added pressure_level_extrapolate, blended-area-weighted, and simplest model_level_area_weighted coarse-graining strategies.
  - Implemented area-weighted mean over the dominant surface type for the snoalb, shdmin, and shdmax surface fields.
  - Fixed a bug that affected coarse-graining the soil type field and some other soil properties downstream.
  - Additional CO2 and surface-type fraction diagnostics.
- Software updates (Spencer, Lucas, Joseph)
  - Updates for running SHiELD in the FV3net Python wrapper
  - Bugfix in ozphys.F for ldiag3d = .true..
  - Cleaned up obscure diagnostics in radiation_astronomy.f.
  - Namelist option gfs_physics_nml::landseaprt to allow turning off land/sea partitioned global intervals in stdout.
  - Added variable-by-variable defined missing_value to gfdl_diag_type.
  - Add additional atmos driver files to run SHiELD and SHiEMOM with the full FMScoupler, mimicking the structure of AM4 atmos_phys.
  - Several diagnostic bug fixes.


# RELEASE NOTES for FV3 202305: Summary
FV3-202305-public --- May 2023
Lucas Harris, GFDL lucas.harris@noaa.gov

This version has been tested with the FV3 Dynamical Core release 202305
and with FMS release 2023.01 from https://github.com/NOAA-GFDL/FMS

- COSP diagnostics (Linjiong)
  - Cleaned up the interface of the COSP and call COSP only at the diagnostic time step
  - Added functionality for COSP to run offline
- Modified shallow convection (Kun):
  - Introduced a new file mfshalcnv_gfdl.f, which will be used if imfshalcnv is set to 4. In the new scheme, if the diagnosed cloud depth or cloud top exceeds certain critical values (set by cthk_shal and top_shal, respectively), shallow convection will not be called.
  - Introduced more diagnostic fields related to the modified shallow convection scheme
- GFDL MP (Linjiong)
  - Removed unused 3d microphysics diagnostics to save time and memory
  - Update gfdl_mp_nml reading code to avoid model crash for absent nml
  - Added options to sub-cycling condensation evaporation (nconds), control timescale or evaporation (do_evap_timescale), and delay condensation and evaporation (delay_cond_evap)
  - Removed grid size in GFDL MP energy and mass calculation
- Fixed potential memory access outside of allocated arrays when NOAHMP is turned on (Kai-Yuan)
- Added check for TKE tracer if TKE_EDMF scheme is used and added default values for undefined values (Lucas)
- Add a function to use IFS initial SST for short-term forecast (Jan-Huey)
  - Namelist parameter: use_ifs_ini_sst
  - The IFS sst data on 6 tiles (ifsSST_data_tile*.nc) need to be put in the INPUT dir


# RELEASE NOTES for FV3 202210: Summary
FV3-202210-public --- October 2022
Lucas Harris, GFDL lucas.harris@noaa.gov

This version has been tested with the FV3 Dynamical Core release 202210
and with FMS release 2022.03 from https://github.com/NOAA-GFDL/FMS

This release includes the following:
- Release of the GFDL Microphysics Version 3
- Fix for first_time_step bug
- Fix for bug in which global physics diagnostic messages were printed out every timestep.
- Fix specification of time-averaged flux variable
- Fix segmentation fault when writing coarse-grained surface restart files
- COSP Implementation
- NOAH MP update
- Introduce namelist parameter, Ts0, that is used to specify the surface temperature
- Added option to fix solar declination for doubly periodic experiments
- Diagnostic totprcp_ave has been renamed to totprcpb_ave
