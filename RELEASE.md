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
