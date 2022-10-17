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
