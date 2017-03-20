!> @file radiation_model_mod.f90
!------------------------------------------------------------------------------!
! This file is part of PALM.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2016 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: radiation_model_mod.f90 2012 2016-09-19 17:31:38Z kanani $
!
! 2011 2016-09-19 17:29:57Z kanani
! Removed CALL of auxiliary SUBROUTINE get_usm_info,
! flag urban_surface is now defined in module control_parameters.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added calculation of solar directional vector for new urban surface 
! model,
! accounted for urban_surface model in radiation_check_parameters,
! correction of comments for zenith angle.
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1976 2016-07-27 13:28:04Z maronga
! Output of 2D/3D/masked data is now directly done within this module. The
! radiation schemes have been simplified for better usability so that
! rad_lw_in, rad_lw_out, rad_sw_in, and rad_sw_out are available independent of
! the radiation code used.
! 
! 1856 2016-04-13 12:56:17Z maronga
! Bugfix: allocation of rad_lw_out for radiation_scheme = 'clear-sky'
! 
! 1853 2016-04-11 09:00:35Z maronga
! Added routine for radiation_scheme = constant.
!  
! 1849 2016-04-08 11:33:18Z hoffmann 
! Adapted for modularization of microphysics
!
! 1826 2016-04-07 12:01:39Z maronga
! Further modularization.
! 
! 1788 2016-03-10 11:01:04Z maronga
! Added new albedo class for pavements / roads.
!
! 1783 2016-03-06 18:36:17Z raasch
! palm-netcdf-module removed in order to avoid a circular module dependency,
! netcdf-variables moved to netcdf-module, new routine netcdf_handle_error_rad
! added
!
! 1757 2016-02-22 15:49:32Z maronga
! Added parameter unscheduled_radiation_calls. Bugfix: interpolation of sounding
! profiles for pressure and temperature above the LES domain.
! 
! 1709 2015-11-04 14:47:01Z maronga
! Bugfix: set initial value for rrtm_lwuflx_dt to zero, small formatting
! corrections
! 
! 1701 2015-11-02 07:43:04Z maronga
! Bugfixes: wrong index for output of timeseries, setting of nz_snd_end
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added option for spin-up runs without radiation (skip_time_do_radiation). Bugfix
! in calculation of pressure profiles. Bugfix in calculation of trace gas profiles.
! Added output of radiative heating rates.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1606 2015-06-29 10:43:37Z maronga
! Added preprocessor directive __netcdf to allow for compiling without netCDF.
! Note, however, that RRTMG cannot be used without netCDF.
! 
! 1590 2015-05-08 13:56:27Z maronga
! Bugfix: definition of character strings requires same length for all elements
! 
! 1587 2015-05-04 14:19:01Z maronga
! Added albedo class for snow
! 
! 1585 2015-04-30 07:05:52Z maronga
! Added support for RRTMG
! 
! 1571 2015-03-12 16:12:49Z maronga
! Added missing KIND attribute. Removed upper-case variable names
! 
! 1551 2015-03-03 14:18:16Z maronga
! Added support for data output. Various variables have been renamed. Added
! interface for different radiation schemes (currently: clear-sky, constant, and
! RRTM (not yet implemented).
! 
! 1496 2014-12-02 17:25:50Z maronga
! Initial revision
! 
!
! Description:
! ------------
!> Radiation models and interfaces
!> @todo move variable definitions used in radiation_init only to the subroutine
!>       as they are no longer required after initialization.
!> @todo Output of full column vertical profiles used in RRTMG
!> @todo Output of other rrtm arrays (such as volume mixing ratios)
!> @todo Adapt for use with topography
!>
!> @note Many variables have a leading dummy dimension (0:0) in order to
!>       match the assume-size shape expected by the RRTMG model.
!------------------------------------------------------------------------------!
 MODULE radiation_model_mod
 
    USE arrays_3d,                                                             &
        ONLY:  dzw, hyp, pt, q, ql, zu, zw

    USE cloud_parameters,                                                      &
        ONLY:  cp, l_d_cp, rho_l

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, cloud_physics, g, initializing_actions,         &
               large_scale_forcing, lsf_surf, phi, pt_surface, rho_surface,    &
               surface_pressure, time_since_reference_point

    USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb_s_inner, nzb, nzt

    USE kinds

    USE microphysics_mod,                                                      &
        ONLY:  nc_const, sigma_gc

#if defined ( __netcdf )
    USE NETCDF
#endif

#if defined ( __rrtmg )
    USE parrrsw,                                                               &
        ONLY:  naerec, nbndsw

    USE parrrtm,                                                               &
        ONLY:  nbndlw

    USE rrtmg_lw_init,                                                         &
        ONLY:  rrtmg_lw_ini

    USE rrtmg_sw_init,                                                         &
        ONLY:  rrtmg_sw_ini

    USE rrtmg_lw_rad,                                                          &
        ONLY:  rrtmg_lw

    USE rrtmg_sw_rad,                                                          &
        ONLY:  rrtmg_sw
#endif



    IMPLICIT NONE

    CHARACTER(10) :: radiation_scheme = 'clear-sky' ! 'constant', 'clear-sky', or 'rrtmg'

!
!-- Predefined Land surface classes (albedo_type) after Briegleb (1992)
    CHARACTER(37), DIMENSION(0:17), PARAMETER :: albedo_type_name = (/      &
                                   'user defined                         ', & !  0 
                                   'ocean                                ', & !  1
                                   'mixed farming, tall grassland        ', & !  2
                                   'tall/medium grassland                ', & !  3 
                                   'evergreen shrubland                  ', & !  4
                                   'short grassland/meadow/shrubland     ', & !  5
                                   'evergreen needleleaf forest          ', & !  6
                                   'mixed deciduous evergreen forest     ', & !  7
                                   'deciduous forest                     ', & !  8
                                   'tropical evergreen broadleaved forest', & !  9
                                   'medium/tall grassland/woodland       ', & ! 10
                                   'desert, sandy                        ', & ! 11 
                                   'desert, rocky                        ', & ! 12 
                                   'tundra                               ', & ! 13
                                   'land ice                             ', & ! 14 
                                   'sea ice                              ', & ! 15 
                                   'snow                                 ', & ! 16
                                   'pavement/roads                       '  & ! 17
                                                         /)

    INTEGER(iwp) :: albedo_type  = 5,    & !< Albedo surface type (default: short grassland)
                    day,                 & !< current day of the year
                    day_init     = 172,  & !< day of the year at model start (21/06)
                    dots_rad     = 0       !< starting index for timeseries output 

    LOGICAL ::  unscheduled_radiation_calls = .TRUE., & !< flag parameter indicating whether additional calls of the radiation code are allowed
                constant_albedo = .FALSE.,            & !< flag parameter indicating whether the albedo may change depending on zenith
                force_radiation_call = .FALSE.,       & !< flag parameter for unscheduled radiation calls
                lw_radiation = .TRUE.,                & !< flag parameter indicating whether longwave radiation shall be calculated
                radiation = .FALSE.,                  & !< flag parameter indicating whether the radiation model is used
                sun_up    = .TRUE.,                   & !< flag parameter indicating whether the sun is up or down
                sw_radiation = .TRUE.,                 & !< flag parameter indicing whether shortwave radiation shall be calculated
                sun_direction = .FALSE.                 !< flag parameter indicing whether solar direction shall be calculated


    REAL(wp), PARAMETER :: d_seconds_hour  = 0.000277777777778_wp,  & !< inverse of seconds per hour (1/3600)
                           d_hours_day    = 0.0416666666667_wp,     & !< inverse of hours per day (1/24)
                           sigma_sb       = 5.67037321E-8_wp,       & !< Stefan-Boltzmann constant
                           solar_constant = 1368.0_wp                 !< solar constant at top of atmosphere

    REAL(wp) :: albedo = 9999999.9_wp,           & !< NAMELIST alpha
                albedo_lw_dif = 9999999.9_wp,    & !< NAMELIST aldif
                albedo_lw_dir = 9999999.9_wp,    & !< NAMELIST aldir
                albedo_sw_dif = 9999999.9_wp,    & !< NAMELIST asdif
                albedo_sw_dir = 9999999.9_wp,    & !< NAMELIST asdir
                decl_1,                          & !< declination coef. 1
                decl_2,                          & !< declination coef. 2
                decl_3,                          & !< declination coef. 3
                dt_radiation = 0.0_wp,           & !< radiation model timestep
                emissivity = 0.98_wp,            & !< NAMELIST surface emissivity
                lambda = 0.0_wp,                 & !< longitude in degrees
                lon = 0.0_wp,                    & !< longitude in radians
                lat = 0.0_wp,                    & !< latitude in radians
                net_radiation = 0.0_wp,          & !< net radiation at surface
                skip_time_do_radiation = 0.0_wp, & !< Radiation model is not called before this time
                sky_trans,                       & !< sky transmissivity
                time_radiation = 0.0_wp,         & !< time since last call of radiation code
                time_utc,                        & !< current time in UTC
                time_utc_init = 43200.0_wp         !< UTC time at model start (noon)

    REAL(wp), DIMENSION(0:0) ::  zenith,         & !< cosine of solar zenith angle
                                 sun_dir_lat,    & !< solar directional vector in latitudes
                                 sun_dir_lon       !< solar directional vector in longitudes

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
                alpha,                       & !< surface broadband albedo (used for clear-sky scheme)
                rad_lw_out_change_0,         & !< change in LW out due to change in surface temperature
                rad_net,                     & !< net radiation at the surface
                rad_net_av                     !< average of rad_net

!
!-- Land surface albedos for solar zenith angle of 60° after Briegleb (1992)     
!-- (shortwave, longwave, broadband):   sw,      lw,      bb,
    REAL(wp), DIMENSION(0:2,1:17), PARAMETER :: albedo_pars = RESHAPE( (/& 
                                   0.06_wp, 0.06_wp, 0.06_wp,            & !  1
                                   0.09_wp, 0.28_wp, 0.19_wp,            & !  2
                                   0.11_wp, 0.33_wp, 0.23_wp,            & !  3
                                   0.11_wp, 0.33_wp, 0.23_wp,            & !  4
                                   0.14_wp, 0.34_wp, 0.25_wp,            & !  5
                                   0.06_wp, 0.22_wp, 0.14_wp,            & !  6
                                   0.06_wp, 0.27_wp, 0.17_wp,            & !  7
                                   0.06_wp, 0.31_wp, 0.19_wp,            & !  8
                                   0.06_wp, 0.22_wp, 0.14_wp,            & !  9
                                   0.06_wp, 0.28_wp, 0.18_wp,            & ! 10
                                   0.35_wp, 0.51_wp, 0.43_wp,            & ! 11
                                   0.24_wp, 0.40_wp, 0.32_wp,            & ! 12
                                   0.10_wp, 0.27_wp, 0.19_wp,            & ! 13
                                   0.90_wp, 0.65_wp, 0.77_wp,            & ! 14
                                   0.90_wp, 0.65_wp, 0.77_wp,            & ! 15
                                   0.95_wp, 0.70_wp, 0.82_wp,            & ! 16
                                   0.08_wp, 0.08_wp, 0.08_wp             & ! 17
                                 /), (/ 3, 17 /) )

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: &
                        rad_lw_cs_hr,                  & !< longwave clear sky radiation heating rate (K/s)
                        rad_lw_cs_hr_av,               & !< average of rad_lw_cs_hr
                        rad_lw_hr,                     & !< longwave radiation heating rate (K/s)
                        rad_lw_hr_av,                  & !< average of rad_sw_hr
                        rad_lw_in,                     & !< incoming longwave radiation (W/m2)
                        rad_lw_in_av,                  & !< average of rad_lw_in
                        rad_lw_out,                    & !< outgoing longwave radiation (W/m2)
                        rad_lw_out_av,                 & !< average of rad_lw_out
                        rad_sw_cs_hr,                  & !< shortwave clear sky radiation heating rate (K/s)
                        rad_sw_cs_hr_av,               & !< average of rad_sw_cs_hr
                        rad_sw_hr,                     & !< shortwave radiation heating rate (K/s)
                        rad_sw_hr_av,                  & !< average of rad_sw_hr
                        rad_sw_in,                     & !< incoming shortwave radiation (W/m2)
                        rad_sw_in_av,                  & !< average of rad_sw_in
                        rad_sw_out,                    & !< outgoing shortwave radiation (W/m2)
                        rad_sw_out_av                    !< average of rad_sw_out


!
!-- Variables and parameters used in RRTMG only
#if defined ( __rrtmg )
    CHARACTER(LEN=12) :: rrtm_input_file = "RAD_SND_DATA" !< name of the NetCDF input file (sounding data)


!
!-- Flag parameters for RRTMGS (should not be changed)
    INTEGER(iwp), PARAMETER :: rrtm_idrv     = 1, & !< flag for longwave upward flux calculation option (0,1)
                               rrtm_inflglw  = 2, & !< flag for lw cloud optical properties (0,1,2)
                               rrtm_iceflglw = 0, & !< flag for lw ice particle specifications (0,1,2,3)
                               rrtm_liqflglw = 1, & !< flag for lw liquid droplet specifications
                               rrtm_inflgsw  = 2, & !< flag for sw cloud optical properties (0,1,2)
                               rrtm_iceflgsw = 0, & !< flag for sw ice particle specifications (0,1,2,3)
                               rrtm_liqflgsw = 1    !< flag for sw liquid droplet specifications

!
!-- The following variables should be only changed with care, as this will 
!-- require further setting of some variables, which is currently not
!-- implemented (aerosols, ice phase).
    INTEGER(iwp) :: nzt_rad,           & !< upper vertical limit for radiation calculations
                    rrtm_icld = 0,     & !< cloud flag (0: clear sky column, 1: cloudy column)
                    rrtm_iaer = 0        !< aerosol option flag (0: no aerosol layers, for lw only: 6 (requires setting of rrtm_sw_ecaer), 10: one or more aerosol layers (not implemented)

    INTEGER(iwp) :: nc_stat !< local variable for storin the result of netCDF calls for error message handling

    LOGICAL :: snd_exists = .FALSE.      !< flag parameter to check whether a user-defined input files exists

    REAL(wp), PARAMETER :: mol_mass_air_d_wv = 1.607793_wp !< molecular weight dry air / water vapor

    REAL(wp), DIMENSION(:), ALLOCATABLE :: hyp_snd,     & !< hypostatic pressure from sounding data (hPa)
                                           q_snd,       & !< specific humidity from sounding data (kg/kg) - dummy at the moment
                                           rrtm_tsfc,   & !< dummy array for storing surface temperature
                                           t_snd          !< actual temperature from sounding data (hPa)

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: aldif,          & !< longwave diffuse albedo solar angle of 60°
                                             aldir,          & !< longwave direct albedo solar angle of 60°
                                             asdif,          & !< shortwave diffuse albedo solar angle of 60°
                                             asdir,          & !< shortwave direct albedo solar angle of 60°
                                             rrtm_ccl4vmr,   & !< CCL4 volume mixing ratio (g/mol)
                                             rrtm_cfc11vmr,  & !< CFC11 volume mixing ratio (g/mol)
                                             rrtm_cfc12vmr,  & !< CFC12 volume mixing ratio (g/mol)
                                             rrtm_cfc22vmr,  & !< CFC22 volume mixing ratio (g/mol)
                                             rrtm_ch4vmr,    & !< CH4 volume mixing ratio
                                             rrtm_cicewp,    & !< in-cloud ice water path (g/m²)
                                             rrtm_cldfr,     & !< cloud fraction (0,1)
                                             rrtm_cliqwp,    & !< in-cloud liquid water path (g/m²)
                                             rrtm_co2vmr,    & !< CO2 volume mixing ratio (g/mol)
                                             rrtm_emis,      & !< surface emissivity (0-1)    
                                             rrtm_h2ovmr,    & !< H2O volume mixing ratio
                                             rrtm_n2ovmr,    & !< N2O volume mixing ratio
                                             rrtm_o2vmr,     & !< O2 volume mixing ratio
                                             rrtm_o3vmr,     & !< O3 volume mixing ratio
                                             rrtm_play,      & !< pressure layers (hPa, zu-grid)
                                             rrtm_plev,      & !< pressure layers (hPa, zw-grid)
                                             rrtm_reice,     & !< cloud ice effective radius (microns)
                                             rrtm_reliq,     & !< cloud water drop effective radius (microns)
                                             rrtm_tlay,      & !< actual temperature (K, zu-grid)
                                             rrtm_tlev,      & !< actual temperature (K, zw-grid)
                                             rrtm_lwdflx,    & !< RRTM output of incoming longwave radiation flux (W/m2)
                                             rrtm_lwdflxc,   & !< RRTM output of outgoing clear sky longwave radiation flux (W/m2) 
                                             rrtm_lwuflx,    & !< RRTM output of outgoing longwave radiation flux (W/m2)
                                             rrtm_lwuflxc,   & !< RRTM output of incoming clear sky longwave radiation flux (W/m2)
                                             rrtm_lwuflx_dt, & !< RRTM output of incoming clear sky longwave radiation flux (W/m2)
                                             rrtm_lwuflxc_dt,& !< RRTM output of outgoing clear sky longwave radiation flux (W/m2)
                                             rrtm_lwhr,      & !< RRTM output of longwave radiation heating rate (K/d)
                                             rrtm_lwhrc,     & !< RRTM output of incoming longwave clear sky radiation heating rate (K/d)
                                             rrtm_swdflx,    & !< RRTM output of incoming shortwave radiation flux (W/m2)
                                             rrtm_swdflxc,   & !< RRTM output of outgoing clear sky shortwave radiation flux (W/m2) 
                                             rrtm_swuflx,    & !< RRTM output of outgoing shortwave radiation flux (W/m2)
                                             rrtm_swuflxc,   & !< RRTM output of incoming clear sky shortwave radiation flux (W/m2)
                                             rrtm_swhr,      & !< RRTM output of shortwave radiation heating rate (K/d)
                                             rrtm_swhrc        !< RRTM output of incoming shortwave clear sky radiation heating rate (K/d) 

!
!-- Definition of arrays that are currently not used for calling RRTMG (due to setting of flag parameters)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rad_lw_cs_in,   & !< incoming clear sky longwave radiation (W/m2) (not used)
                                                rad_lw_cs_out,  & !< outgoing clear sky longwave radiation (W/m2) (not used)
                                                rad_sw_cs_in,   & !< incoming clear sky shortwave radiation (W/m2) (not used)
                                                rad_sw_cs_out,  & !< outgoing clear sky shortwave radiation (W/m2) (not used)
                                                rrtm_aldif,     & !< surface albedo for longwave diffuse radiation
                                                rrtm_aldir,     & !< surface albedo for longwave direct radiation
                                                rrtm_asdif,     & !< surface albedo for shortwave diffuse radiation
                                                rrtm_asdir,     & !< surface albedo for shortwave direct radiation
                                                rrtm_lw_tauaer, & !< lw aerosol optical depth
                                                rrtm_lw_taucld, & !< lw in-cloud optical depth
                                                rrtm_sw_taucld, & !< sw in-cloud optical depth
                                                rrtm_sw_ssacld, & !< sw in-cloud single scattering albedo
                                                rrtm_sw_asmcld, & !< sw in-cloud asymmetry parameter
                                                rrtm_sw_fsfcld, & !< sw in-cloud forward scattering fraction
                                                rrtm_sw_tauaer, & !< sw aerosol optical depth
                                                rrtm_sw_ssaaer, & !< sw aerosol single scattering albedo
                                                rrtm_sw_asmaer, & !< sw aerosol asymmetry parameter
                                                rrtm_sw_ecaer     !< sw aerosol optical detph at 0.55 microns (rrtm_iaer = 6 only)

#endif

    INTERFACE radiation_check_data_output
       MODULE PROCEDURE radiation_check_data_output
    END INTERFACE radiation_check_data_output

    INTERFACE radiation_check_data_output_pr
       MODULE PROCEDURE radiation_check_data_output_pr
    END INTERFACE radiation_check_data_output_pr
  
    INTERFACE radiation_check_parameters
       MODULE PROCEDURE radiation_check_parameters
    END INTERFACE radiation_check_parameters
  
    INTERFACE radiation_clearsky
       MODULE PROCEDURE radiation_clearsky
    END INTERFACE radiation_clearsky
 
    INTERFACE radiation_constant
       MODULE PROCEDURE radiation_constant
    END INTERFACE radiation_constant
 
    INTERFACE radiation_control
       MODULE PROCEDURE radiation_control
    END INTERFACE radiation_control

    INTERFACE radiation_3d_data_averaging
       MODULE PROCEDURE radiation_3d_data_averaging
    END INTERFACE radiation_3d_data_averaging

    INTERFACE radiation_data_output_2d
       MODULE PROCEDURE radiation_data_output_2d
    END INTERFACE radiation_data_output_2d

    INTERFACE radiation_data_output_3d
       MODULE PROCEDURE radiation_data_output_3d
    END INTERFACE radiation_data_output_3d

    INTERFACE radiation_data_output_mask
       MODULE PROCEDURE radiation_data_output_mask
    END INTERFACE radiation_data_output_mask

    INTERFACE radiation_define_netcdf_grid
       MODULE PROCEDURE radiation_define_netcdf_grid
    END INTERFACE radiation_define_netcdf_grid

    INTERFACE radiation_header
       MODULE PROCEDURE radiation_header
    END INTERFACE radiation_header 
  
    INTERFACE radiation_init
       MODULE PROCEDURE radiation_init
    END INTERFACE radiation_init

    INTERFACE radiation_parin
       MODULE PROCEDURE radiation_parin
    END INTERFACE radiation_parin
    
    INTERFACE radiation_rrtmg
       MODULE PROCEDURE radiation_rrtmg
    END INTERFACE radiation_rrtmg

    INTERFACE radiation_tendency
       MODULE PROCEDURE radiation_tendency
       MODULE PROCEDURE radiation_tendency_ij
    END INTERFACE radiation_tendency

    INTERFACE radiation_read_restart_data
       MODULE PROCEDURE radiation_read_restart_data
    END INTERFACE radiation_read_restart_data

    INTERFACE radiation_last_actions
       MODULE PROCEDURE radiation_last_actions
    END INTERFACE radiation_last_actions

    SAVE

    PRIVATE

!
!-- Public functions / NEEDS SORTING
    PUBLIC radiation_check_data_output, radiation_check_data_output_pr,        &
           radiation_check_parameters, radiation_control,                      &
           radiation_header, radiation_init, radiation_parin,                  &
           radiation_3d_data_averaging, radiation_tendency,                    &
           radiation_data_output_2d, radiation_data_output_3d,                 &
           radiation_define_netcdf_grid, radiation_last_actions,               &
           radiation_read_restart_data, radiation_data_output_mask
    
!
!-- Public variables and constants / NEEDS SORTING
    PUBLIC dots_rad, dt_radiation, force_radiation_call,                       &
           rad_net, rad_net_av, radiation, radiation_scheme, rad_lw_in,        &
           rad_lw_in_av, rad_lw_out, rad_lw_out_av, rad_lw_out_change_0,       &
           rad_lw_cs_hr, rad_lw_cs_hr_av, rad_lw_hr, rad_lw_hr_av, rad_sw_in,  &
           rad_sw_in_av, rad_sw_out, rad_sw_out_av, rad_sw_cs_hr,              &
           rad_sw_cs_hr_av, rad_sw_hr, rad_sw_hr_av, sigma_sb,                 &
           skip_time_do_radiation, time_radiation, unscheduled_radiation_calls,&
           zenith, calc_zenith, sun_direction, sun_dir_lat, sun_dir_lon,       &
           day_init, time_utc_init


#if defined ( __rrtmg )
    PUBLIC rrtm_aldif, rrtm_aldir, rrtm_asdif, rrtm_asdir
#endif

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine controls the calls of the radiation schemes
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_control
 
 
       IMPLICIT NONE


       SELECT CASE ( TRIM( radiation_scheme ) )

          CASE ( 'constant' )
             CALL radiation_constant
          
          CASE ( 'clear-sky' )
             CALL radiation_clearsky
       
          CASE ( 'rrtmg' )
             CALL radiation_rrtmg

          CASE DEFAULT

       END SELECT


    END SUBROUTINE radiation_control

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_check_data_output( var, unit, i, ilen, k )
 
 
       USE control_parameters,                                                 &
           ONLY: data_output, message_string

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  unit     !< 
       CHARACTER (LEN=*) ::  var      !<

       INTEGER(iwp) :: i
       INTEGER(iwp) :: ilen
       INTEGER(iwp) :: k

       SELECT CASE ( TRIM( var ) )

         CASE ( 'rad_lw_cs_hr', 'rad_lw_hr', 'rad_sw_cs_hr', 'rad_sw_hr' )
             IF (  .NOT.  radiation  .OR.  radiation_scheme /= 'rrtmg' )  THEN
                message_string = '"output of "' // TRIM( var ) // '" requi' // &
                                 'res radiation = .TRUE. and ' //              &
                                 'radiation_scheme = "rrtmg"'
                CALL message( 'check_parameters', 'PA0406', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'W/m2'     

          CASE ( 'rad_net*', 'rrtm_aldif*', 'rrtm_aldir*', 'rrtm_asdif*',      &
                 'rrtm_asdir*' )
             IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
                message_string = 'illegal value for data_output: "' //         &
                                 TRIM( var ) // '" & only 2d-horizontal ' //   &
                                 'cross sections are allowed for this value'
                CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
             ENDIF
             IF (  .NOT.  radiation  .OR.  radiation_scheme /= "rrtmg" )  THEN
                IF ( TRIM( var ) == 'rrtm_aldif*'  .OR.                        &
                     TRIM( var ) == 'rrtm_aldir*'  .OR.                        &
                     TRIM( var ) == 'rrtm_asdif*'  .OR.                        &
                     TRIM( var ) == 'rrtm_asdir*'      )                       &
                THEN
                   message_string = 'output of "' // TRIM( var ) // '" require'&
                                    // 's radiation = .TRUE. and radiation_sch'&
                                    // 'eme = "rrtmg"'
                   CALL message( 'check_parameters', 'PA0409', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

             IF ( TRIM( var ) == 'rad_net*'      ) unit = 'W/m2'   
             IF ( TRIM( var ) == 'rrtm_aldif*'   ) unit = ''   
             IF ( TRIM( var ) == 'rrtm_aldir*'   ) unit = '' 
             IF ( TRIM( var ) == 'rrtm_asdif*'   ) unit = '' 
             IF ( TRIM( var ) == 'rrtm_asdir*'   ) unit = '' 

          CASE DEFAULT
             unit = 'illegal'

       END SELECT


    END SUBROUTINE radiation_check_data_output

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for radiation model
!------------------------------------------------------------------------------!  
    SUBROUTINE radiation_check_data_output_pr( variable, var_count, unit, dopr_unit )
 
       USE arrays_3d,                                                          &
           ONLY: zu

       USE control_parameters,                                                 &
           ONLY: data_output_pr, message_string

       USE indices

       USE profil_parameter

       USE statistics

       IMPLICIT NONE
   
       CHARACTER (LEN=*) ::  unit      !< 
       CHARACTER (LEN=*) ::  variable  !< 
       CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit
 
       INTEGER(iwp) ::  user_pr_index !< 
       INTEGER(iwp) ::  var_count     !< 

       SELECT CASE ( TRIM( variable ) )
       
         CASE ( 'rad_net' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme == 'constant' )&
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 101
                dopr_unit  = 'W/m2'
                hom(:,2,101,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_in' )
             IF ( (  .NOT.  radiation)  .OR.  radiation_scheme == 'constant' ) &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 102
                dopr_unit  = 'W/m2'
                hom(:,2,102,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit 
             ENDIF

          CASE ( 'rad_lw_out' )
             IF ( (  .NOT. radiation )  .OR.  radiation_scheme == 'constant' ) &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 103
                dopr_unit  = 'W/m2'
                hom(:,2,103,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit   
             ENDIF

          CASE ( 'rad_sw_in' )
             IF ( (  .NOT. radiation )  .OR.  radiation_scheme == 'constant' ) &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 104
                dopr_unit  = 'W/m2'
                hom(:,2,104,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_out')
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme == 'constant' )&
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme = "constant"'
                CALL message( 'check_parameters', 'PA0408', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 105
                dopr_unit  = 'W/m2'
                hom(:,2,105,:)  = SPREAD( zw, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_cs_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 106
                dopr_unit  = 'K/h'
                hom(:,2,106,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_lw_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 107
                dopr_unit  = 'K/h'
                hom(:,2,107,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_cs_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 108
                dopr_unit  = 'K/h'
                hom(:,2,108,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF

          CASE ( 'rad_sw_hr' )
             IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )   &
             THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(var_count) ) // ' is' // &
                                 'not available for radiation = .FALSE. or ' //&
                                 'radiation_scheme /= "rrtmg"'
                CALL message( 'check_parameters', 'PA0413', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(var_count) = 109
                dopr_unit  = 'K/h'
                hom(:,2,109,:)  = SPREAD( zu, 2, statistic_regions+1 )
                unit = dopr_unit
             ENDIF


          CASE DEFAULT
             unit = 'illegal'

       END SELECT


    END SUBROUTINE radiation_check_data_output_pr
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_check_parameters

       USE control_parameters,                                                 &
           ONLY: message_string, topography, urban_surface
                 
    
       IMPLICIT NONE
       

       IF ( radiation_scheme /= 'constant'   .AND.                             &
            radiation_scheme /= 'clear-sky'  .AND.                             &
            radiation_scheme /= 'rrtmg' )  THEN
          message_string = 'unknown radiation_scheme = '//                     &
                           TRIM( radiation_scheme )
          CALL message( 'check_parameters', 'PA0405', 1, 2, 0, 6, 0 )
       ELSEIF ( radiation_scheme == 'rrtmg' )  THEN
#if ! defined ( __rrtmg )
          message_string = 'radiation_scheme = "rrtmg" requires ' //           & 
                           'compilation of PALM with pre-processor ' //        &
                           'directive -D__rrtmg'
          CALL message( 'check_parameters', 'PA0407', 1, 2, 0, 6, 0 )
#endif
#if defined ( __rrtmg ) && ! defined( __netcdf )
          message_string = 'radiation_scheme = "rrtmg" requires ' //           & 
                           'the use of NetCDF (preprocessor directive ' //     &
                           '-D__netcdf'
          CALL message( 'check_parameters', 'PA0412', 1, 2, 0, 6, 0 )
#endif

       ENDIF

       IF ( albedo_type == 0  .AND.  albedo == 9999999.9_wp  .AND.             &
            radiation_scheme == 'clear-sky')  THEN
          message_string = 'radiation_scheme = "clear-sky" in combination' //  & 
                           'with albedo_type = 0 requires setting of albedo'// &
                           ' /= 9999999.9'
          CALL message( 'check_parameters', 'PA0410', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( albedo_type == 0  .AND.  radiation_scheme == 'rrtmg'  .AND.        &
          (    albedo_lw_dif == 9999999.9_wp .OR. albedo_lw_dir == 9999999.9_wp&
          .OR. albedo_sw_dif == 9999999.9_wp .OR. albedo_sw_dir == 9999999.9_wp& 
          ) ) THEN
          message_string = 'radiation_scheme = "rrtmg" in combination' //      & 
                           'with albedo_type = 0 requires setting of ' //      &
                           'albedo_lw_dif /= 9999999.9' //                     &
                           'albedo_lw_dir /= 9999999.9' //                     &
                           'albedo_sw_dif /= 9999999.9 and' //                 &
                           'albedo_sw_dir /= 9999999.9'
          CALL message( 'check_parameters', 'PA0411', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    The following paramter check is temporarily extended by the urban_surface
!--    flag, until a better solution comes up to omit this check in case of
!--    urban surface model is used.
       IF ( topography /= 'flat'  .AND.  .NOT.  urban_surface )  THEN
          message_string = 'radiation scheme cannot be used ' //               & 
                           'in combination with  topography /= "flat"'
          CALL message( 'check_parameters', 'PA0414', 1, 2, 0, 6, 0 )
       ENDIF
 
    END SUBROUTINE radiation_check_parameters 
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_init
    
       IMPLICIT NONE

!
!--    Allocate array for storing the surface net radiation
       IF ( .NOT. ALLOCATED ( rad_net ) )  THEN
          ALLOCATE ( rad_net(nysg:nyng,nxlg:nxrg) )
          rad_net = 0.0_wp
       ENDIF

!
!--    Allocate array for storing the surface net radiation
       IF ( .NOT. ALLOCATED ( rad_lw_out_change_0 ) )  THEN
          ALLOCATE ( rad_lw_out_change_0(nysg:nyng,nxlg:nxrg) )
          rad_lw_out_change_0 = 0.0_wp
       ENDIF

!
!--    Fix net radiation in case of radiation_scheme = 'constant'
       IF ( radiation_scheme == 'constant' )  THEN
          rad_net = net_radiation
!          radiation = .FALSE.
!
!--    Calculate orbital constants
       ELSE
          decl_1 = SIN(23.45_wp * pi / 180.0_wp)
          decl_2 = 2.0_wp * pi / 365.0_wp
          decl_3 = decl_2 * 81.0_wp
          lat    = phi * pi / 180.0_wp
          lon    = lambda * pi / 180.0_wp
       ENDIF

       IF ( radiation_scheme == 'clear-sky'  .OR.                              &
            radiation_scheme == 'constant')  THEN

          ALLOCATE ( alpha(nysg:nyng,nxlg:nxrg) )

          IF ( .NOT. ALLOCATED ( rad_sw_in ) )  THEN
             ALLOCATE ( rad_sw_in(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_sw_out ) )  THEN
             ALLOCATE ( rad_sw_out(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_in_av ) )  THEN
             ALLOCATE ( rad_sw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_sw_out_av ) )  THEN
             ALLOCATE ( rad_sw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in ) )  THEN
             ALLOCATE ( rad_lw_in(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_lw_out ) )  THEN
             ALLOCATE ( rad_lw_out(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in_av ) )  THEN
             ALLOCATE ( rad_lw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( .NOT. ALLOCATED ( rad_lw_out_av ) )  THEN
             ALLOCATE ( rad_lw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
          ENDIF

          rad_sw_in  = 0.0_wp
          rad_sw_out = 0.0_wp
          rad_lw_in  = 0.0_wp
          rad_lw_out = 0.0_wp

!
!--       Overwrite albedo if manually set in parameter file
          IF ( albedo_type /= 0 .AND. albedo == 9999999.9_wp )  THEN
             albedo = albedo_pars(2,albedo_type)
          ENDIF
    
          alpha = albedo
  
!
!--    Initialization actions for RRTMG
       ELSEIF ( radiation_scheme == 'rrtmg' )  THEN
#if defined ( __rrtmg )
!
!--       Allocate albedos
          ALLOCATE ( rrtm_aldif(0:0,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( rrtm_aldir(0:0,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( rrtm_asdif(0:0,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( rrtm_asdir(0:0,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( aldif(nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( aldir(nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( asdif(nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( asdir(nysg:nyng,nxlg:nxrg) )

          IF ( albedo_type /= 0 )  THEN
             IF ( albedo_lw_dif == 9999999.9_wp )  THEN
                albedo_lw_dif = albedo_pars(0,albedo_type)
                albedo_lw_dir = albedo_lw_dif
             ENDIF
             IF ( albedo_sw_dif == 9999999.9_wp )  THEN
                albedo_sw_dif = albedo_pars(1,albedo_type)
                albedo_sw_dir = albedo_sw_dif
             ENDIF
          ENDIF

          aldif(:,:) = albedo_lw_dif
          aldir(:,:) = albedo_lw_dir
          asdif(:,:) = albedo_sw_dif
          asdir(:,:) = albedo_sw_dir
!
!--       Calculate initial values of current (cosine of) the zenith angle and 
!--       whether the sun is up
          CALL calc_zenith     
!
!--       Calculate initial surface albedo
          IF ( .NOT. constant_albedo )  THEN
             CALL calc_albedo
          ELSE
             rrtm_aldif(0,:,:) = aldif(:,:)
             rrtm_aldir(0,:,:) = aldir(:,:)
             rrtm_asdif(0,:,:) = asdif(:,:) 
             rrtm_asdir(0,:,:) = asdir(:,:)   
          ENDIF

!
!--       Allocate surface emissivity
          ALLOCATE ( rrtm_emis(0:0,1:nbndlw+1) )
          rrtm_emis = emissivity

!
!--       Allocate 3d arrays of radiative fluxes and heating rates
          IF ( .NOT. ALLOCATED ( rad_sw_in ) )  THEN
             ALLOCATE ( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_in = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_in_av ) )  THEN
             ALLOCATE ( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_out ) )  THEN
             ALLOCATE ( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_out = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_out_av ) )  THEN
             ALLOCATE ( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_hr ) )  THEN
             ALLOCATE ( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_hr_av ) )  THEN
             ALLOCATE ( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_hr_av = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_cs_hr ) )  THEN
             ALLOCATE ( rad_sw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_cs_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_sw_cs_hr_av ) )  THEN
             ALLOCATE ( rad_sw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_sw_cs_hr_av = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in ) )  THEN
             ALLOCATE ( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_in     = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_in_av ) )  THEN
             ALLOCATE ( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_out ) )  THEN
             ALLOCATE ( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
            rad_lw_out    = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_out_av ) )  THEN
             ALLOCATE ( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_hr ) )  THEN
             ALLOCATE ( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_hr_av ) )  THEN
             ALLOCATE ( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_hr_av = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_cs_hr ) )  THEN
             ALLOCATE ( rad_lw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_cs_hr = 0.0_wp
          ENDIF

          IF ( .NOT. ALLOCATED ( rad_lw_cs_hr_av ) )  THEN
             ALLOCATE ( rad_lw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             rad_lw_cs_hr_av = 0.0_wp
          ENDIF

          ALLOCATE ( rad_sw_cs_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( rad_sw_cs_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_sw_cs_in  = 0.0_wp
          rad_sw_cs_out = 0.0_wp

          ALLOCATE ( rad_lw_cs_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ALLOCATE ( rad_lw_cs_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_lw_cs_in  = 0.0_wp
          rad_lw_cs_out = 0.0_wp

!
!--       Allocate dummy array for storing surface temperature
          ALLOCATE ( rrtm_tsfc(1) )

!
!--       Initialize RRTMG
          IF ( lw_radiation )  CALL rrtmg_lw_ini ( cp )
          IF ( sw_radiation )  CALL rrtmg_sw_ini ( cp )

!
!--       Set input files for RRTMG
          INQUIRE(FILE="RAD_SND_DATA", EXIST=snd_exists) 
          IF ( .NOT. snd_exists )  THEN
             rrtm_input_file = "rrtmg_lw.nc"
          ENDIF

!
!--       Read vertical layers for RRTMG from sounding data
!--       The routine provides nzt_rad, hyp_snd(1:nzt_rad),
!--       t_snd(nzt+2:nzt_rad), rrtm_play(1:nzt_rad), rrtm_plev(1_nzt_rad+1), 
!--       rrtm_tlay(nzt+2:nzt_rad), rrtm_tlev(nzt+2:nzt_rad+1)
          CALL read_sounding_data

!
!--       Read trace gas profiles from file. This routine provides
!--       the rrtm_ arrays (1:nzt_rad+1)
          CALL read_trace_gas_data
#endif
       ENDIF

!
!--    Perform user actions if required
       CALL user_init_radiation

!
!--    Calculate radiative fluxes at model start
       IF (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

          SELECT CASE ( radiation_scheme )
             CASE ( 'rrtmg' )
                CALL radiation_rrtmg
             CASE ( 'clear-sky' )
                CALL radiation_clearsky
             CASE ( 'constant' )
                CALL radiation_constant
             CASE DEFAULT
          END SELECT

       ENDIF

       RETURN

    END SUBROUTINE radiation_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> A simple clear sky radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_clearsky


       IMPLICIT NONE

       INTEGER(iwp) :: i, j, k   !< loop indices
       REAL(wp)     :: exn,   &  !< Exner functions at surface
                       exn1,  &  !< Exner functions at first grid level
                       pt1       !< potential temperature at first grid level

!
!--    Calculate current zenith angle
       CALL calc_zenith

!
!--    Calculate sky transmissivity
       sky_trans = 0.6_wp + 0.2_wp * zenith(0)

!
!--    Calculate value of the Exner function
       exn = (surface_pressure / 1000.0_wp )**0.286_wp
!
!--    Calculate radiation fluxes and net radiation (rad_net) for each grid 
!--    point
       DO i = nxlg, nxrg
          DO j = nysg, nyng
             k = nzb_s_inner(j,i)

             exn1 = (hyp(k+1) / 100000.0_wp )**0.286_wp

             rad_sw_in(0,j,i)  = solar_constant * sky_trans * zenith(0)
             rad_sw_out(0,j,i) = alpha(j,i) * rad_sw_in(0,j,i)
             rad_lw_out(0,j,i) = emissivity * sigma_sb * (pt(k,j,i) * exn)**4

             IF ( cloud_physics )  THEN
                pt1 = pt(k+1,j,i) + l_d_cp / exn1 * ql(k+1,j,i)
                rad_lw_in(0,j,i)  = 0.8_wp * sigma_sb * (pt1 * exn1)**4
             ELSE
                rad_lw_in(0,j,i)  = 0.8_wp * sigma_sb * (pt(k+1,j,i) * exn1)**4
             ENDIF

             rad_net(j,i) = rad_sw_in(0,j,i) - rad_sw_out(0,j,i)               &
                            + rad_lw_in(0,j,i) - rad_lw_out(0,j,i)


             rad_lw_out_change_0(j,i) = 3.0_wp * sigma_sb * emissivity         &
                                        * (pt(k,j,i) * exn) ** 3

          ENDDO
       ENDDO

    END SUBROUTINE radiation_clearsky


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This scheme keeps the prescribed net radiation constant during the run
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_constant


       IMPLICIT NONE

       INTEGER(iwp) :: i, j, k   !< loop indices
       REAL(wp)     :: exn,   &  !< Exner functions at surface
                       exn1,  &  !< Exner functions at first grid level
                       pt1       !< potential temperature at first grid level

!
!--    Calculate value of the Exner function
       exn = (surface_pressure / 1000.0_wp )**0.286_wp
!
!--    Prescribe net radiation and estimate the remaining radiative fluxes
       DO i = nxlg, nxrg
          DO j = nysg, nyng
             k = nzb_s_inner(j,i)

             rad_net(j,i)      = net_radiation

             exn1 = (hyp(k+1) / 100000.0_wp )**0.286_wp

             IF ( cloud_physics )  THEN
                pt1 = pt(k+1,j,i) + l_d_cp / exn1 * ql(k+1,j,i)
                rad_lw_in(0,j,i)  = 0.8_wp * sigma_sb * (pt1 * exn1)**4
             ELSE
                rad_lw_in(0,j,i)  = 0.8_wp * sigma_sb * (pt(k+1,j,i) * exn1)**4
             ENDIF

             rad_lw_out(0,j,i) = emissivity * sigma_sb * (pt(k,j,i) * exn)**4

             rad_sw_in(0,j,i) = ( rad_net(j,i) - rad_lw_in(0,j,i)              &
                                  + rad_lw_out(0,j,i) )                        &
                                  / ( 1.0_wp - alpha(j,i) )

          ENDDO
       ENDDO

    END SUBROUTINE radiation_constant

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_header ( io )


       IMPLICIT NONE
 
       INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file
    

       
!
!--    Write radiation model header
       WRITE( io, 3 )

       IF ( radiation_scheme == "constant" )  THEN
          WRITE( io, 4 ) net_radiation
       ELSEIF ( radiation_scheme == "clear-sky" )  THEN
          WRITE( io, 5 )
       ELSEIF ( radiation_scheme == "rrtmg" )  THEN
          WRITE( io, 6 )
          IF ( .NOT. lw_radiation )  WRITE( io, 10 )
          IF ( .NOT. sw_radiation )  WRITE( io, 11 )
       ENDIF 

       IF ( albedo_type == 0 )  THEN
          WRITE( io, 7 ) albedo
       ELSE
          WRITE( io, 8 ) TRIM( albedo_type_name(albedo_type) )
       ENDIF
       IF ( constant_albedo )  THEN
          WRITE( io, 9 )
       ENDIF
       
       IF ( radiation .AND. radiation_scheme /= 'constant' )  THEN
          WRITE ( io, 1 )  lambda
          WRITE ( io, 2 )  day_init, time_utc_init
       ENDIF

       WRITE( io, 12 ) dt_radiation
 

 1 FORMAT ('    Geograph. longitude            :   lambda = ',F4.1,' degr')
 2 FORMAT ('    Day of the year at model start :   day_init = ',I3   &
            /'    UTC time at model start        :   time_utc_init = ',F7.1' s')
 3 FORMAT (//' Radiation model information:'/                                  &
              ' ----------------------------'/)
 4 FORMAT ('    --> Using constant net radiation: net_radiation = ', F6.2,        &
           // 'W/m**2')
 5 FORMAT ('    --> Simple radiation scheme for clear sky is used (no clouds,',   &
                   ' default)')
 6 FORMAT ('    --> RRTMG scheme is used')
 7 FORMAT (/'    User-specific surface albedo: albedo =', F6.3)
 8 FORMAT (/'    Albedo is set for land surface type: ', A)
 9 FORMAT (/'    --> Albedo is fixed during the run')
10 FORMAT (/'    --> Longwave radiation is disabled')
11 FORMAT (/'    --> Shortwave radiation is disabled.')
12 FORMAT  ('    Timestep: dt_radiation = ', F6.2, '  s')


    END SUBROUTINE radiation_header
   

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &radiation_par for radiation model
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_parin


       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file 
       
       NAMELIST /radiation_par/   albedo, albedo_type, albedo_lw_dir,          &
                                  albedo_lw_dif, albedo_sw_dir, albedo_sw_dif, &
                                  constant_albedo, day_init, dt_radiation,     &
                                  lambda, lw_radiation, net_radiation,         &
                                  radiation_scheme, skip_time_do_radiation,    &
                                  sw_radiation, time_utc_init,                 &
                                  unscheduled_radiation_calls
       
       line = ' '
       
!
!--    Try to find radiation model package
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&radiation_par' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, radiation_par )

!
!--    Set flag that indicates that the radiation model is switched on
       radiation = .TRUE.

 10    CONTINUE
       

    END SUBROUTINE radiation_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Implementation of the RRTMG radiation_scheme
!------------------------------------------------------------------------------!
    SUBROUTINE radiation_rrtmg

       USE indices,                                                            &
           ONLY:  nbgp

       USE particle_attributes,                                                &
           ONLY:  grid_particles, number_of_particles, particles,              &
                  particle_advection_start, prt_count

       IMPLICIT NONE

#if defined ( __rrtmg )

       INTEGER(iwp) :: i, j, k, n !< loop indices

       REAL(wp)     ::  s_r2, &   !< weighted sum over all droplets with r^2
                        s_r3      !< weighted sum over all droplets with r^3

!
!--    Calculate current (cosine of) zenith angle and whether the sun is up
       CALL calc_zenith     
!
!--    Calculate surface albedo
       IF ( .NOT. constant_albedo )  THEN
          CALL calc_albedo
       ENDIF

!
!--    Prepare input data for RRTMG

!
!--    In case of large scale forcing with surface data, calculate new pressure
!--    profile. nzt_rad might be modified by these calls and all required arrays
!--    will then be re-allocated
       IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
          CALL read_sounding_data
          CALL read_trace_gas_data
       ENDIF
!
!--    Loop over all grid points
       DO i = nxl, nxr
          DO j = nys, nyn

!
!--          Prepare profiles of temperature and H2O volume mixing ratio
             rrtm_tlev(0,nzb+1) = pt(nzb,j,i) * ( surface_pressure             &
                                                  / 1000.0_wp )**0.286_wp

             DO k = nzb+1, nzt+1
                rrtm_tlay(0,k) = pt(k,j,i) * ( (hyp(k) ) / 100000.0_wp         &
                                 )**0.286_wp + l_d_cp * ql(k,j,i)
                rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * (q(k,j,i) - ql(k,j,i))

             ENDDO

!
!--          Avoid temperature/humidity jumps at the top of the LES domain by 
!--          linear interpolation from nzt+2 to nzt+7
             DO k = nzt+2, nzt+7
                rrtm_tlay(0,k) = rrtm_tlay(0,nzt+1)                            &
                              + ( rrtm_tlay(0,nzt+8) - rrtm_tlay(0,nzt+1) )    &
                              / ( rrtm_play(0,nzt+8) - rrtm_play(0,nzt+1) )    &
                              * ( rrtm_play(0,k) - rrtm_play(0,nzt+1) )

                rrtm_h2ovmr(0,k) = rrtm_h2ovmr(0,nzt+1)                        &
                              + ( rrtm_h2ovmr(0,nzt+8) - rrtm_h2ovmr(0,nzt+1) )&
                              / ( rrtm_play(0,nzt+8)   - rrtm_play(0,nzt+1)   )&
                              * ( rrtm_play(0,k) - rrtm_play(0,nzt+1) )

             ENDDO

!--          Linear interpolate to zw grid
             DO k = nzb+2, nzt+8
                rrtm_tlev(0,k)   = rrtm_tlay(0,k-1) + (rrtm_tlay(0,k) -        &
                                   rrtm_tlay(0,k-1))                           &
                                   / ( rrtm_play(0,k) - rrtm_play(0,k-1) )     &
                                   * ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
             ENDDO


!
!--          Calculate liquid water path and cloud fraction for each column.
!--          Note that LWP is required in g/m² instead of kg/kg m.
             rrtm_cldfr  = 0.0_wp
             rrtm_reliq  = 0.0_wp
             rrtm_cliqwp = 0.0_wp
             rrtm_icld   = 0

             DO k = nzb+1, nzt+1
                rrtm_cliqwp(0,k) =  ql(k,j,i) * 1000.0_wp *                    &
                                    (rrtm_plev(0,k) - rrtm_plev(0,k+1))        &
                                    * 100.0_wp / g 

                IF ( rrtm_cliqwp(0,k) > 0.0_wp )  THEN
                   rrtm_cldfr(0,k) = 1.0_wp
                   IF ( rrtm_icld == 0 )  rrtm_icld = 1

!
!--                Calculate cloud droplet effective radius
                   IF ( cloud_physics )  THEN
                      rrtm_reliq(0,k) = 1.0E6_wp * ( 3.0_wp * ql(k,j,i)        &
                                        * rho_surface                          &
                                        / ( 4.0_wp * pi * nc_const * rho_l )   &
                                        )**0.33333333333333_wp                 &
                                        * EXP( LOG( sigma_gc )**2 )

                   ELSEIF ( cloud_droplets )  THEN
                      number_of_particles = prt_count(k,j,i)

                      IF (number_of_particles <= 0)  CYCLE
                      particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                      s_r2 = 0.0_wp
                      s_r3 = 0.0_wp

                      DO  n = 1, number_of_particles
                         IF ( particles(n)%particle_mask )  THEN
                            s_r2 = s_r2 + particles(n)%radius**2 * &
                                   particles(n)%weight_factor
                            s_r3 = s_r3 + particles(n)%radius**3 * &
                                   particles(n)%weight_factor
                         ENDIF
                      ENDDO

                      IF ( s_r2 > 0.0_wp )  rrtm_reliq(0,k) = s_r3 / s_r2

                   ENDIF

!
!--                Limit effective radius
                   IF ( rrtm_reliq(0,k) > 0.0_wp )  THEN
                      rrtm_reliq(0,k) = MAX(rrtm_reliq(0,k),2.5_wp)
                      rrtm_reliq(0,k) = MIN(rrtm_reliq(0,k),60.0_wp)
                  ENDIF
                ENDIF
             ENDDO

!
!--          Set surface temperature
             rrtm_tsfc = pt(nzb,j,i) * (surface_pressure / 1000.0_wp )**0.286_wp

             IF ( lw_radiation )  THEN
               CALL rrtmg_lw( 1, nzt_rad      , rrtm_icld    , rrtm_idrv      ,&
               rrtm_play       , rrtm_plev    , rrtm_tlay    , rrtm_tlev      ,&
               rrtm_tsfc       , rrtm_h2ovmr  , rrtm_o3vmr   , rrtm_co2vmr    ,&
               rrtm_ch4vmr     , rrtm_n2ovmr  , rrtm_o2vmr   , rrtm_cfc11vmr  ,&
               rrtm_cfc12vmr   , rrtm_cfc22vmr, rrtm_ccl4vmr , rrtm_emis      ,&
               rrtm_inflglw    , rrtm_iceflglw, rrtm_liqflglw, rrtm_cldfr     ,&
               rrtm_lw_taucld  , rrtm_cicewp  , rrtm_cliqwp  , rrtm_reice     ,& 
               rrtm_reliq      , rrtm_lw_tauaer,                               &
               rrtm_lwuflx     , rrtm_lwdflx  , rrtm_lwhr  ,                   &
               rrtm_lwuflxc    , rrtm_lwdflxc , rrtm_lwhrc ,                   &
               rrtm_lwuflx_dt  ,  rrtm_lwuflxc_dt )

!
!--             Save fluxes
                DO k = nzb, nzt+1
                   rad_lw_in(k,j,i)  = rrtm_lwdflx(0,k)
                   rad_lw_out(k,j,i) = rrtm_lwuflx(0,k)
                ENDDO

!
!--             Save heating rates (convert from K/d to K/h)
                DO k = nzb+1, nzt+1
                   rad_lw_hr(k,j,i)     = rrtm_lwhr(0,k)  * d_hours_day
                   rad_lw_cs_hr(k,j,i)  = rrtm_lwhrc(0,k) * d_hours_day
                ENDDO

!
!--             Save change in LW heating rate
                rad_lw_out_change_0(j,i) = rrtm_lwuflx_dt(0,nzb)

             ENDIF

             IF ( sw_radiation .AND. sun_up )  THEN
                CALL rrtmg_sw( 1, nzt_rad      , rrtm_icld  , rrtm_iaer       ,&
               rrtm_play       , rrtm_plev    , rrtm_tlay  , rrtm_tlev        ,&
               rrtm_tsfc       , rrtm_h2ovmr  , rrtm_o3vmr , rrtm_co2vmr      ,&
               rrtm_ch4vmr     , rrtm_n2ovmr  , rrtm_o2vmr , rrtm_asdir(:,j,i),&
               rrtm_asdif(:,j,i), rrtm_aldir(:,j,i), rrtm_aldif(:,j,i), zenith,&
               0.0_wp          , day          , solar_constant,   rrtm_inflgsw,&
               rrtm_iceflgsw   , rrtm_liqflgsw, rrtm_cldfr , rrtm_sw_taucld   ,&
               rrtm_sw_ssacld  , rrtm_sw_asmcld, rrtm_sw_fsfcld, rrtm_cicewp  ,&
               rrtm_cliqwp     , rrtm_reice   , rrtm_reliq , rrtm_sw_tauaer   ,&
               rrtm_sw_ssaaer     , rrtm_sw_asmaer  , rrtm_sw_ecaer ,          &
               rrtm_swuflx     , rrtm_swdflx  , rrtm_swhr  ,                   &
               rrtm_swuflxc    , rrtm_swdflxc , rrtm_swhrc )
  
!
!--             Save fluxes
                DO k = nzb, nzt+1
                   rad_sw_in(k,j,i)  = rrtm_swdflx(0,k)
                   rad_sw_out(k,j,i) = rrtm_swuflx(0,k)
                ENDDO

!
!--             Save heating rates (convert from K/d to K/s)
                DO k = nzb+1, nzt+1
                   rad_sw_hr(k,j,i)     = rrtm_swhr(0,k)  * d_hours_day
                   rad_sw_cs_hr(k,j,i)  = rrtm_swhrc(0,k) * d_hours_day
                ENDDO

             ENDIF

!
!--          Calculate surface net radiation
             rad_net(j,i) = rad_sw_in(nzb,j,i) - rad_sw_out(nzb,j,i)           &
                            + rad_lw_in(nzb,j,i) - rad_lw_out(nzb,j,i)

          ENDDO
       ENDDO

       CALL exchange_horiz( rad_lw_in,  nbgp )
       CALL exchange_horiz( rad_lw_out, nbgp )
       CALL exchange_horiz( rad_lw_hr,    nbgp )
       CALL exchange_horiz( rad_lw_cs_hr, nbgp )

       CALL exchange_horiz( rad_sw_in,  nbgp )
       CALL exchange_horiz( rad_sw_out, nbgp ) 
       CALL exchange_horiz( rad_sw_hr,    nbgp )
       CALL exchange_horiz( rad_sw_cs_hr, nbgp )

       CALL exchange_horiz_2d( rad_net, nbgp )
       CALL exchange_horiz_2d( rad_lw_out_change_0, nbgp )
#endif

    END SUBROUTINE radiation_rrtmg


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the cosine of the zenith angle (variable is called zenith)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_zenith

       IMPLICIT NONE

       REAL(wp) ::  declination,  & !< solar declination angle
                    hour_angle      !< solar hour angle
!
!--    Calculate current day and time based on the initial values and simulation
!--    time
       day = day_init + INT(FLOOR( (time_utc_init + time_since_reference_point)    &
                               / 86400.0_wp ), KIND=iwp)
       time_utc = MOD((time_utc_init + time_since_reference_point), 86400.0_wp)


!
!--    Calculate solar declination and hour angle   
       declination = ASIN( decl_1 * SIN(decl_2 * REAL(day, KIND=wp) - decl_3) )
       hour_angle  = 2.0_wp * pi * (time_utc / 86400.0_wp) + lon - pi

!
!--    Calculate cosine of solar zenith angle
       zenith(0) = SIN(lat) * SIN(declination) + COS(lat) * COS(declination)      &
                                            * COS(hour_angle)
       zenith(0) = MAX(0.0_wp,zenith(0))

!
!--    Calculate solar directional vector
       IF ( sun_direction )  THEN
!--       Direction in longitudes equals to sin(solar_azimuth) * sin(zenith)
          sun_dir_lon(0) = -SIN(hour_angle) * COS(declination)
!--       Direction in latitues equals to cos(solar_azimuth) * sin(zenith)
          sun_dir_lat(0) = SIN(declination) * COS(lat) - COS(hour_angle) &
                              * COS(declination) * SIN(lat)
       ENDIF

!
!--    Check if the sun is up (otheriwse shortwave calculations can be skipped)
       IF ( zenith(0) > 0.0_wp )  THEN
          sun_up = .TRUE.
       ELSE
          sun_up = .FALSE.
       END IF

    END SUBROUTINE calc_zenith

#if defined ( __rrtmg ) && defined ( __netcdf )
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates surface albedo components based on Briegleb (1992) and 
!> Briegleb et al. (1986)
!------------------------------------------------------------------------------!
    SUBROUTINE calc_albedo

        IMPLICIT NONE

        IF ( sun_up )  THEN
!
!--        Ocean
           IF ( albedo_type == 1 )  THEN
              rrtm_aldir(0,:,:) = 0.026_wp / ( zenith(0)**1.7_wp + 0.065_wp )  &
                                  + 0.15_wp * ( zenith(0) - 0.1_wp )           &
                                            * ( zenith(0) - 0.5_wp )           &
                                            * ( zenith(0) - 1.0_wp )
              rrtm_asdir(0,:,:) = rrtm_aldir(0,:,:)
!
!--        Snow
           ELSEIF ( albedo_type == 16 )  THEN
              IF ( zenith(0) < 0.5_wp )  THEN
                 rrtm_aldir(0,:,:) = 0.5_wp * (1.0_wp - aldif)                 &
                                     * ( 3.0_wp / (1.0_wp + 4.0_wp             &
                                     * zenith(0))) - 1.0_wp
                 rrtm_asdir(0,:,:) = 0.5_wp * (1.0_wp - asdif)                 &
                                     * ( 3.0_wp / (1.0_wp + 4.0_wp             &
                                     * zenith(0))) - 1.0_wp

                 rrtm_aldir(0,:,:) = MIN(0.98_wp, rrtm_aldir(0,:,:))
                 rrtm_asdir(0,:,:) = MIN(0.98_wp, rrtm_asdir(0,:,:))
              ELSE
                 rrtm_aldir(0,:,:) = aldif
                 rrtm_asdir(0,:,:) = asdif
              ENDIF
!
!--        Sea ice
           ELSEIF ( albedo_type == 15 )  THEN
                 rrtm_aldir(0,:,:) = aldif
                 rrtm_asdir(0,:,:) = asdif

!
!--        Asphalt
           ELSEIF ( albedo_type == 17 )  THEN
                 rrtm_aldir(0,:,:) = aldif
                 rrtm_asdir(0,:,:) = asdif
!
!--        Land surfaces
           ELSE
              SELECT CASE ( albedo_type )

!
!--              Surface types with strong zenith dependence
                 CASE ( 1, 2, 3, 4, 11, 12, 13 )
                    rrtm_aldir(0,:,:) = aldif * 1.4_wp /                       &
                                        (1.0_wp + 0.8_wp * zenith(0))
                    rrtm_asdir(0,:,:) = asdif * 1.4_wp /                       &
                                        (1.0_wp + 0.8_wp * zenith(0))
!
!--              Surface types with weak zenith dependence
                 CASE ( 5, 6, 7, 8, 9, 10, 14 )
                    rrtm_aldir(0,:,:) = aldif * 1.1_wp /                       &
                                        (1.0_wp + 0.2_wp * zenith(0))
                    rrtm_asdir(0,:,:) = asdif * 1.1_wp /                       &
                                        (1.0_wp + 0.2_wp * zenith(0))

                 CASE DEFAULT

              END SELECT
           ENDIF
!
!--        Diffusive albedo is taken from Table 2
           rrtm_aldif(0,:,:) = aldif
           rrtm_asdif(0,:,:) = asdif

        ELSE

           rrtm_aldir(0,:,:) = 0.0_wp
           rrtm_asdir(0,:,:) = 0.0_wp
           rrtm_aldif(0,:,:) = 0.0_wp
           rrtm_asdif(0,:,:) = 0.0_wp
        ENDIF
    END SUBROUTINE calc_albedo

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read sounding data (pressure and temperature) from RADIATION_DATA.
!------------------------------------------------------------------------------!
    SUBROUTINE read_sounding_data

       IMPLICIT NONE

       INTEGER(iwp) :: id,           & !< NetCDF id of input file
                       id_dim_zrad,  & !< pressure level id in the NetCDF file
                       id_var,       & !< NetCDF variable id
                       k,            & !< loop index
                       nz_snd,       & !< number of vertical levels in the sounding data
                       nz_snd_start, & !< start vertical index for sounding data to be used
                       nz_snd_end      !< end vertical index for souding data to be used

       REAL(wp) :: t_surface           !< actual surface temperature

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  hyp_snd_tmp, & !< temporary hydrostatic pressure profile (sounding)
                                               t_snd_tmp      !< temporary temperature profile (sounding)

!
!--    In case of updates, deallocate arrays first (sufficient to check one
!--    array as the others are automatically allocated). This is required
!--    because nzt_rad might change during the update
       IF ( ALLOCATED ( hyp_snd ) )  THEN
          DEALLOCATE( hyp_snd )
          DEALLOCATE( t_snd )
          DEALLOCATE( q_snd  )
          DEALLOCATE ( rrtm_play )
          DEALLOCATE ( rrtm_plev )
          DEALLOCATE ( rrtm_tlay )
          DEALLOCATE ( rrtm_tlev )

          DEALLOCATE ( rrtm_h2ovmr )
          DEALLOCATE ( rrtm_cicewp )
          DEALLOCATE ( rrtm_cldfr )
          DEALLOCATE ( rrtm_cliqwp )
          DEALLOCATE ( rrtm_reice )
          DEALLOCATE ( rrtm_reliq )
          DEALLOCATE ( rrtm_lw_taucld )
          DEALLOCATE ( rrtm_lw_tauaer )

          DEALLOCATE ( rrtm_lwdflx  )
          DEALLOCATE ( rrtm_lwdflxc )
          DEALLOCATE ( rrtm_lwuflx  )
          DEALLOCATE ( rrtm_lwuflxc )
          DEALLOCATE ( rrtm_lwuflx_dt )
          DEALLOCATE ( rrtm_lwuflxc_dt )
          DEALLOCATE ( rrtm_lwhr  )
          DEALLOCATE ( rrtm_lwhrc )

          DEALLOCATE ( rrtm_sw_taucld )
          DEALLOCATE ( rrtm_sw_ssacld )
          DEALLOCATE ( rrtm_sw_asmcld )
          DEALLOCATE ( rrtm_sw_fsfcld )
          DEALLOCATE ( rrtm_sw_tauaer )
          DEALLOCATE ( rrtm_sw_ssaaer )
          DEALLOCATE ( rrtm_sw_asmaer ) 
          DEALLOCATE ( rrtm_sw_ecaer )   
 
          DEALLOCATE ( rrtm_swdflx  )
          DEALLOCATE ( rrtm_swdflxc )
          DEALLOCATE ( rrtm_swuflx  )
          DEALLOCATE ( rrtm_swuflxc )
          DEALLOCATE ( rrtm_swhr  )
          DEALLOCATE ( rrtm_swhrc )

       ENDIF

!
!--    Open file for reading
       nc_stat = NF90_OPEN( rrtm_input_file, NF90_NOWRITE, id )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 549 )

!
!--    Inquire dimension of z axis and save in nz_snd
       nc_stat = NF90_INQ_DIMID( id, "Pressure", id_dim_zrad )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim_zrad, len = nz_snd )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 551 )

!
! !--    Allocate temporary array for storing pressure data
       ALLOCATE( hyp_snd_tmp(1:nz_snd) )
       hyp_snd_tmp = 0.0_wp


!--    Read pressure from file
       nc_stat = NF90_INQ_VARID( id, "Pressure", id_var )
       nc_stat = NF90_GET_VAR( id, id_var, hyp_snd_tmp(:), start = (/1/),      &
                               count = (/nz_snd/) )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 552 )

!
!--    Allocate temporary array for storing temperature data
       ALLOCATE( t_snd_tmp(1:nz_snd) )
       t_snd_tmp = 0.0_wp

!
!--    Read temperature from file
       nc_stat = NF90_INQ_VARID( id, "ReferenceTemperature", id_var )
       nc_stat = NF90_GET_VAR( id, id_var, t_snd_tmp(:), start = (/1/),        &
                               count = (/nz_snd/) )
       CALL netcdf_handle_error_rad( 'read_sounding_data', 553 )

!
!--    Calculate start of sounding data
       nz_snd_start = nz_snd + 1
       nz_snd_end   = nz_snd + 1

!
!--    Start filling vertical dimension at 10hPa above the model domain (hyp is
!--    in Pa, hyp_snd in hPa).
       DO  k = 1, nz_snd
          IF ( hyp_snd_tmp(k) < ( hyp(nzt+1) - 1000.0_wp) * 0.01_wp )  THEN
             nz_snd_start = k
             EXIT
          END IF
       END DO

       IF ( nz_snd_start <= nz_snd )  THEN
          nz_snd_end = nz_snd
       END IF


!
!--    Calculate of total grid points for RRTMG calculations
       nzt_rad = nzt + nz_snd_end - nz_snd_start + 1

!
!--    Save data above LES domain in hyp_snd, t_snd and q_snd
!--    Note: q_snd_tmp is not calculated at the moment (dry residual atmosphere)
       ALLOCATE( hyp_snd(nzb+1:nzt_rad) )
       ALLOCATE( t_snd(nzb+1:nzt_rad)   )
       ALLOCATE( q_snd(nzb+1:nzt_rad)   )
       hyp_snd = 0.0_wp
       t_snd = 0.0_wp
       q_snd = 0.0_wp

       hyp_snd(nzt+2:nzt_rad) = hyp_snd_tmp(nz_snd_start+1:nz_snd_end)
       t_snd(nzt+2:nzt_rad)   = t_snd_tmp(nz_snd_start+1:nz_snd_end)

       nc_stat = NF90_CLOSE( id )

!
!--    Calculate pressure levels on zu and zw grid. Sounding data is added at 
!--    top of the LES domain. This routine does not consider horizontal or 
!--    vertical variability of pressure and temperature
       ALLOCATE ( rrtm_play(0:0,nzb+1:nzt_rad+1)   )
       ALLOCATE ( rrtm_plev(0:0,nzb+1:nzt_rad+2)   )

       t_surface = pt_surface * ( surface_pressure / 1000.0_wp )**0.286_wp
       DO k = nzb+1, nzt+1
          rrtm_play(0,k) = hyp(k) * 0.01_wp
          rrtm_plev(0,k) = surface_pressure * ( (t_surface - g/cp * zw(k-1)) / &
                         t_surface )**(1.0_wp/0.286_wp)
       ENDDO

       DO k = nzt+2, nzt_rad
          rrtm_play(0,k) = hyp_snd(k)
          rrtm_plev(0,k) = 0.5_wp * ( rrtm_play(0,k) + rrtm_play(0,k-1) )
       ENDDO
       rrtm_plev(0,nzt_rad+1) = MAX( 0.5 * hyp_snd(nzt_rad),                   &
                                   1.5 * hyp_snd(nzt_rad)                      &
                                 - 0.5 * hyp_snd(nzt_rad-1) )
       rrtm_plev(0,nzt_rad+2)  = MIN( 1.0E-4_wp,                               &
                                      0.25_wp * rrtm_plev(0,nzt_rad+1) )

       rrtm_play(0,nzt_rad+1) = 0.5 * rrtm_plev(0,nzt_rad+1)

!
!--    Calculate temperature/humidity levels at top of the LES domain. 
!--    Currently, the temperature is taken from sounding data (might lead to a 
!--    temperature jump at interface. To do: Humidity is currently not 
!--    calculated above the LES domain.
       ALLOCATE ( rrtm_tlay(0:0,nzb+1:nzt_rad+1)   )
       ALLOCATE ( rrtm_tlev(0:0,nzb+1:nzt_rad+2)   )
       ALLOCATE ( rrtm_h2ovmr(0:0,nzb+1:nzt_rad+1) )

       DO k = nzt+8, nzt_rad
          rrtm_tlay(0,k)   = t_snd(k)
          rrtm_h2ovmr(0,k) = q_snd(k)
       ENDDO
       rrtm_tlay(0,nzt_rad+1) = 2.0_wp * rrtm_tlay(0,nzt_rad)                 &
                                - rrtm_tlay(0,nzt_rad-1)
       DO k = nzt+9, nzt_rad+1
          rrtm_tlev(0,k)   = rrtm_tlay(0,k-1) + (rrtm_tlay(0,k)                &
                             - rrtm_tlay(0,k-1))                               &
                             / ( rrtm_play(0,k) - rrtm_play(0,k-1) )           &
                             * ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
       ENDDO
       rrtm_h2ovmr(0,nzt_rad+1) = rrtm_h2ovmr(0,nzt_rad)

       rrtm_tlev(0,nzt_rad+2)   = 2.0_wp * rrtm_tlay(0,nzt_rad+1)              &
                                  - rrtm_tlev(0,nzt_rad)
!
!--    Allocate remaining RRTMG arrays
       ALLOCATE ( rrtm_cicewp(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_cldfr(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_cliqwp(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_reice(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_reliq(0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_lw_taucld(1:nbndlw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_lw_tauaer(0:0,nzb+1:nzt_rad+1,1:nbndlw+1) )
       ALLOCATE ( rrtm_sw_taucld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_ssacld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_asmcld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_fsfcld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
       ALLOCATE ( rrtm_sw_tauaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
       ALLOCATE ( rrtm_sw_ssaaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
       ALLOCATE ( rrtm_sw_asmaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) ) 
       ALLOCATE ( rrtm_sw_ecaer(0:0,nzb+1:nzt_rad+1,1:naerec+1) )    

!
!--    The ice phase is currently not considered in PALM
       rrtm_cicewp = 0.0_wp
       rrtm_reice  = 0.0_wp

!
!--    Set other parameters (move to NAMELIST parameters in the future)
       rrtm_lw_tauaer = 0.0_wp
       rrtm_lw_taucld = 0.0_wp
       rrtm_sw_taucld = 0.0_wp
       rrtm_sw_ssacld = 0.0_wp
       rrtm_sw_asmcld = 0.0_wp
       rrtm_sw_fsfcld = 0.0_wp
       rrtm_sw_tauaer = 0.0_wp
       rrtm_sw_ssaaer = 0.0_wp
       rrtm_sw_asmaer = 0.0_wp
       rrtm_sw_ecaer  = 0.0_wp


       ALLOCATE ( rrtm_swdflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_swuflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_swhr(0:0,nzb+1:nzt_rad+1)  )
       ALLOCATE ( rrtm_swuflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_swdflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_swhrc(0:0,nzb+1:nzt_rad+1) )

       rrtm_swdflx  = 0.0_wp
       rrtm_swuflx  = 0.0_wp
       rrtm_swhr    = 0.0_wp  
       rrtm_swuflxc = 0.0_wp
       rrtm_swdflxc = 0.0_wp
       rrtm_swhrc   = 0.0_wp

       ALLOCATE ( rrtm_lwdflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_lwuflx(0:0,nzb:nzt_rad+1)  )
       ALLOCATE ( rrtm_lwhr(0:0,nzb+1:nzt_rad+1)  )
       ALLOCATE ( rrtm_lwuflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_lwdflxc(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_lwhrc(0:0,nzb+1:nzt_rad+1) )

       rrtm_lwdflx  = 0.0_wp
       rrtm_lwuflx  = 0.0_wp
       rrtm_lwhr    = 0.0_wp  
       rrtm_lwuflxc = 0.0_wp
       rrtm_lwdflxc = 0.0_wp
       rrtm_lwhrc   = 0.0_wp

       ALLOCATE ( rrtm_lwuflx_dt(0:0,nzb:nzt_rad+1) )
       ALLOCATE ( rrtm_lwuflxc_dt(0:0,nzb:nzt_rad+1) )

       rrtm_lwuflx_dt = 0.0_wp
       rrtm_lwuflxc_dt = 0.0_wp

    END SUBROUTINE read_sounding_data


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read trace gas data from file
!------------------------------------------------------------------------------!
    SUBROUTINE read_trace_gas_data

       USE rrsw_ncpar

       IMPLICIT NONE

       INTEGER(iwp), PARAMETER :: num_trace_gases = 9 !< number of trace gases (absorbers)

       CHARACTER(LEN=5), DIMENSION(num_trace_gases), PARAMETER ::              & !< trace gas names
           trace_names = (/'O3   ', 'CO2  ', 'CH4  ', 'N2O  ', 'O2   ',        &
                           'CFC11', 'CFC12', 'CFC22', 'CCL4 '/)

       INTEGER(iwp) :: id,     & !< NetCDF id
                       k,      & !< loop index
                       m,      & !< loop index
                       n,      & !< loop index
                       nabs,   & !< number of absorbers
                       np,     & !< number of pressure levels
                       id_abs, & !< NetCDF id of the respective absorber
                       id_dim, & !< NetCDF id of asborber's dimension
                       id_var    !< NetCDf id ot the absorber

       REAL(wp) :: p_mls_l, p_mls_u, p_wgt_l, p_wgt_u, p_mls_m


       REAL(wp), DIMENSION(:), ALLOCATABLE   ::  p_mls,         & !< pressure levels for the absorbers
                                                 rrtm_play_tmp, & !< temporary array for pressure zu-levels
                                                 rrtm_plev_tmp, & !< temporary array for pressure zw-levels
                                                 trace_path_tmp   !< temporary array for storing trace gas path data

       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  trace_mls,      & !< array for storing the absorber amounts
                                                 trace_mls_path, & !< array for storing trace gas path data
                                                 trace_mls_tmp     !< temporary array for storing trace gas data


!
!--    In case of updates, deallocate arrays first (sufficient to check one
!--    array as the others are automatically allocated)
       IF ( ALLOCATED ( rrtm_o3vmr ) )  THEN
          DEALLOCATE ( rrtm_o3vmr  )
          DEALLOCATE ( rrtm_co2vmr )
          DEALLOCATE ( rrtm_ch4vmr )
          DEALLOCATE ( rrtm_n2ovmr )
          DEALLOCATE ( rrtm_o2vmr  )
          DEALLOCATE ( rrtm_cfc11vmr )
          DEALLOCATE ( rrtm_cfc12vmr )
          DEALLOCATE ( rrtm_cfc22vmr )
          DEALLOCATE ( rrtm_ccl4vmr  )
       ENDIF

!
!--    Allocate trace gas profiles
       ALLOCATE ( rrtm_o3vmr(0:0,1:nzt_rad+1)  )
       ALLOCATE ( rrtm_co2vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_ch4vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_n2ovmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_o2vmr(0:0,1:nzt_rad+1)  )
       ALLOCATE ( rrtm_cfc11vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_cfc12vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_cfc22vmr(0:0,1:nzt_rad+1) )
       ALLOCATE ( rrtm_ccl4vmr(0:0,1:nzt_rad+1)  )

!
!--    Open file for reading
       nc_stat = NF90_OPEN( rrtm_input_file, NF90_NOWRITE, id )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 549 )
!
!--    Inquire dimension ids and dimensions
       nc_stat = NF90_INQ_DIMID( id, "Pressure", id_dim )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, len = np) 
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )

       nc_stat = NF90_INQ_DIMID( id, "Absorber", id_dim )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, len = nabs ) 
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
   

!
!--    Allocate pressure, and trace gas arrays     
       ALLOCATE( p_mls(1:np) )
       ALLOCATE( trace_mls(1:num_trace_gases,1:np) ) 
       ALLOCATE( trace_mls_tmp(1:nabs,1:np) ) 


       nc_stat = NF90_INQ_VARID( id, "Pressure", id_var )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_GET_VAR( id, id_var, p_mls )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )

       nc_stat = NF90_INQ_VARID( id, "AbsorberAmountMLS", id_var )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
       nc_stat = NF90_GET_VAR( id, id_var, trace_mls_tmp )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )


!
!--    Write absorber amounts (mls) to trace_mls
       DO n = 1, num_trace_gases
          CALL getAbsorberIndex( TRIM( trace_names(n) ), id_abs )

          trace_mls(n,1:np) = trace_mls_tmp(id_abs,1:np)

!
!--       Replace missing values by zero
          WHERE ( trace_mls(n,:) > 2.0_wp )  
             trace_mls(n,:) = 0.0_wp
          END WHERE
       END DO

       DEALLOCATE ( trace_mls_tmp )

       nc_stat = NF90_CLOSE( id )
       CALL netcdf_handle_error_rad( 'read_trace_gas_data', 551 )

!
!--    Add extra pressure level for calculations of the trace gas paths
       ALLOCATE ( rrtm_play_tmp(1:nzt_rad+1) )
       ALLOCATE ( rrtm_plev_tmp(1:nzt_rad+2) )

       rrtm_play_tmp(1:nzt_rad)   = rrtm_play(0,1:nzt_rad) 
       rrtm_plev_tmp(1:nzt_rad+1) = rrtm_plev(0,1:nzt_rad+1)
       rrtm_play_tmp(nzt_rad+1)   = rrtm_plev(0,nzt_rad+1) * 0.5_wp
       rrtm_plev_tmp(nzt_rad+2)   = MIN( 1.0E-4_wp, 0.25_wp                    &
                                         * rrtm_plev(0,nzt_rad+1) )
 
!
!--    Calculate trace gas path (zero at surface) with interpolation to the
!--    sounding levels
       ALLOCATE ( trace_mls_path(1:nzt_rad+2,1:num_trace_gases) )

       trace_mls_path(nzb+1,:) = 0.0_wp
       
       DO k = nzb+2, nzt_rad+2
          DO m = 1, num_trace_gases
             trace_mls_path(k,m) = trace_mls_path(k-1,m)

!
!--          When the pressure level is higher than the trace gas pressure
!--          level, assume that 
             IF ( rrtm_plev_tmp(k-1) > p_mls(1) )  THEN             
                
                trace_mls_path(k,m) = trace_mls_path(k,m) + trace_mls(m,1)     &
                                      * ( rrtm_plev_tmp(k-1)                   &
                                          - MAX( p_mls(1), rrtm_plev_tmp(k) )  &
                                        ) / g
             ENDIF

!
!--          Integrate for each sounding level from the contributing p_mls 
!--          levels
             DO n = 2, np
!
!--             Limit p_mls so that it is within the model level
                p_mls_u = MIN( rrtm_plev_tmp(k-1),                             &
                          MAX( rrtm_plev_tmp(k), p_mls(n) ) )
                p_mls_l = MIN( rrtm_plev_tmp(k-1),                             &
                          MAX( rrtm_plev_tmp(k), p_mls(n-1) ) )

                IF ( p_mls_l > p_mls_u )  THEN

!
!--                Calculate weights for interpolation
                   p_mls_m = 0.5_wp * (p_mls_l + p_mls_u)
                   p_wgt_u = (p_mls(n-1) - p_mls_m) / (p_mls(n-1) - p_mls(n))
                   p_wgt_l = (p_mls_m - p_mls(n))   / (p_mls(n-1) - p_mls(n))

!
!--                Add level to trace gas path
                   trace_mls_path(k,m) = trace_mls_path(k,m)                   &
                                         +  ( p_wgt_u * trace_mls(m,n)         &
                                            + p_wgt_l * trace_mls(m,n-1) )     &
                                         * (p_mls_l - p_mls_u) / g
                ENDIF
             ENDDO

             IF ( rrtm_plev_tmp(k) < p_mls(np) )  THEN 
                trace_mls_path(k,m) = trace_mls_path(k,m) + trace_mls(m,np)    &
                                      * ( MIN( rrtm_plev_tmp(k-1), p_mls(np) ) &
                                          - rrtm_plev_tmp(k)                   &
                                        ) / g 
             ENDIF  
          ENDDO
       ENDDO


!
!--    Prepare trace gas path profiles
       ALLOCATE ( trace_path_tmp(1:nzt_rad+1) )

       DO m = 1, num_trace_gases

          trace_path_tmp(1:nzt_rad+1) = ( trace_mls_path(2:nzt_rad+2,m)        &
                                       - trace_mls_path(1:nzt_rad+1,m) ) * g   &
                                       / ( rrtm_plev_tmp(1:nzt_rad+1)          &
                                       - rrtm_plev_tmp(2:nzt_rad+2) )

!
!--       Save trace gas paths to the respective arrays
          SELECT CASE ( TRIM( trace_names(m) ) )

             CASE ( 'O3' )

                rrtm_o3vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CO2' )

                rrtm_co2vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CH4' )

                rrtm_ch4vmr(0,:) = trace_path_tmp(:)

             CASE ( 'N2O' )

                rrtm_n2ovmr(0,:) = trace_path_tmp(:)

             CASE ( 'O2' )

                rrtm_o2vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CFC11' )

                rrtm_cfc11vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CFC12' )

                rrtm_cfc12vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CFC22' )

                rrtm_cfc22vmr(0,:) = trace_path_tmp(:)

             CASE ( 'CCL4' )

                rrtm_ccl4vmr(0,:) = trace_path_tmp(:)

             CASE DEFAULT

          END SELECT

       ENDDO

       DEALLOCATE ( trace_path_tmp )
       DEALLOCATE ( trace_mls_path )
       DEALLOCATE ( rrtm_play_tmp )
       DEALLOCATE ( rrtm_plev_tmp )
       DEALLOCATE ( trace_mls )
       DEALLOCATE ( p_mls )

    END SUBROUTINE read_trace_gas_data


    SUBROUTINE netcdf_handle_error_rad( routine_name, errno )

       USE control_parameters,                                                 &
           ONLY:  message_string

       USE NETCDF

       USE pegrid

       IMPLICIT NONE

       CHARACTER(LEN=6) ::  message_identifier
       CHARACTER(LEN=*) ::  routine_name

       INTEGER(iwp) ::  errno

       IF ( nc_stat /= NF90_NOERR )  THEN

          WRITE( message_identifier, '(''NC'',I4.4)' )  errno
          message_string = TRIM( NF90_STRERROR( nc_stat ) )

          CALL message( routine_name, message_identifier, 2, 2, 0, 6, 1 )

       ENDIF

    END SUBROUTINE netcdf_handle_error_rad
#endif


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate temperature tendency due to radiative cooling/heating.
!> Cache-optimized version.
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_tendency_ij ( i, j, tend )

    USE cloud_parameters,                                                      &
        ONLY:  pt_d_t

    IMPLICIT NONE

    INTEGER(iwp) :: i, j, k !< loop indices

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) :: tend !< pt tendency term

    IF ( radiation_scheme == 'rrtmg' )  THEN
#if defined  ( __rrtmg )
!
!--    Calculate tendency based on heating rate
       DO k = nzb+1, nzt+1
          tend(k,j,i) = tend(k,j,i) + (rad_lw_hr(k,j,i) + rad_sw_hr(k,j,i))    &
                                         * pt_d_t(k) * d_seconds_hour
       ENDDO
#endif
    ENDIF

    END SUBROUTINE radiation_tendency_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate temperature tendency due to radiative cooling/heating.
!> Vector-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_tendency ( tend )

    USE cloud_parameters,                                                      &
        ONLY:  pt_d_t

    USE indices,                                                               &
        ONLY:  nxl, nxr, nyn, nys

    IMPLICIT NONE

    INTEGER(iwp) :: i, j, k !< loop indices

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) :: tend !< pt tendency term

    IF ( radiation_scheme == 'rrtmg' )  THEN
#if defined  ( __rrtmg )
!
!--    Calculate tendency based on heating rate
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO k = nzb+1, nzt+1
                tend(k,j,i) = tend(k,j,i) + ( rad_lw_hr(k,j,i)                 &
                                          +  rad_sw_hr(k,j,i) ) * pt_d_t(k)    &
                                          * d_seconds_hour
             ENDDO
          ENDDO
       ENDDO
#endif
    ENDIF


 END SUBROUTINE radiation_tendency

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!------------------------------------------------------------------------------!
SUBROUTINE radiation_3d_data_averaging( mode, variable )
 

    USE control_parameters

    USE indices

    USE kinds

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode    !< 
    CHARACTER (LEN=*) :: variable !< 

    INTEGER(iwp) ::  i !< 
    INTEGER(iwp) ::  j !< 
    INTEGER(iwp) ::  k !< 

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

             CASE ( 'rad_net*' )
                IF ( .NOT. ALLOCATED( rad_net_av ) )  THEN
                   ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_net_av = 0.0_wp

             CASE ( 'rad_lw_in' )
                IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
                   ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_in_av = 0.0_wp

             CASE ( 'rad_lw_out' )
                IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
                   ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_out_av = 0.0_wp

             CASE ( 'rad_lw_cs_hr' )
                IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
                   ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_cs_hr_av = 0.0_wp

             CASE ( 'rad_lw_hr' )
                IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
                   ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_hr_av = 0.0_wp

             CASE ( 'rad_sw_in' )
                IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
                   ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_in_av = 0.0_wp

             CASE ( 'rad_sw_out' )
                IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
                   ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_out_av = 0.0_wp

             CASE ( 'rad_sw_cs_hr' )
                IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
                   ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_cs_hr_av = 0.0_wp

             CASE ( 'rad_sw_hr' )
                IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
                   ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_hr_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'rad_net*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   rad_net_av(j,i) = rad_net_av(j,i) + rad_net(j,i)
                ENDDO
             ENDDO

          CASE ( 'rad_lw_in' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_lw_in_av(k,j,i) = rad_lw_in_av(k,j,i) + rad_lw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_lw_out' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_lw_out_av(k,j,i) = rad_lw_out_av(k,j,i) + rad_lw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_lw_cs_hr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_lw_cs_hr_av(k,j,i) = rad_lw_cs_hr_av(k,j,i) + rad_lw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_lw_hr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_lw_hr_av(k,j,i) = rad_lw_hr_av(k,j,i) + rad_lw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_sw_in' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_sw_in_av(k,j,i) = rad_sw_in_av(k,j,i) + rad_sw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_sw_out' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_sw_out_av(k,j,i) = rad_sw_out_av(k,j,i) + rad_sw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_sw_cs_hr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_sw_cs_hr_av(k,j,i) = rad_sw_cs_hr_av(k,j,i) + rad_sw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_sw_hr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_sw_hr_av(k,j,i) = rad_sw_hr_av(k,j,i) + rad_sw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

         CASE ( 'rad_net*' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   rad_net_av(j,i) = rad_net_av(j,i) / REAL( average_count_3d, KIND=wp )
                ENDDO
             ENDDO

          CASE ( 'rad_lw_in' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_lw_in_av(k,j,i) = rad_lw_in_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_lw_out' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_lw_out_av(k,j,i) = rad_lw_out_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_lw_cs_hr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_lw_cs_hr_av(k,j,i) = rad_lw_cs_hr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_lw_hr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_lw_hr_av(k,j,i) = rad_lw_hr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_sw_in' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_sw_in_av(k,j,i) = rad_sw_in_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_sw_out' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_sw_out_av(k,j,i) = rad_sw_out_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_sw_cs_hr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_sw_cs_hr_av(k,j,i) = rad_sw_cs_hr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'rad_sw_hr' )
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      rad_sw_hr_av(k,j,i) = rad_sw_hr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDDO

       END SELECT

    ENDIF

END SUBROUTINE radiation_3d_data_averaging


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
SUBROUTINE radiation_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )
    
    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)  ::  var         !< 
    LOGICAL, INTENT(OUT)           ::  found       !< 
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !< 
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !< 
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !< 

    found  = .TRUE.


!
!-- Check for the grid
    SELECT CASE ( TRIM( var ) )

       CASE ( 'rad_lw_cs_hr', 'rad_lw_hr', 'rad_sw_cs_hr', 'rad_sw_hr',        &
              'rad_lw_cs_hr_xy', 'rad_lw_hr_xy', 'rad_sw_cs_hr_xy',            &
              'rad_sw_hr_xy', 'rad_lw_cs_hr_xz', 'rad_lw_hr_xz',               &
              'rad_sw_cs_hr_xz', 'rad_sw_hr_xz', 'rad_lw_cs_hr_yz',            &
              'rad_lw_hr_yz', 'rad_sw_cs_hr_yz', 'rad_sw_hr_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'rad_lw_in', 'rad_lw_out', 'rad_sw_in', 'rad_sw_out',            &
              'rad_lw_in_xy', 'rad_lw_out_xy', 'rad_sw_in_xy','rad_sw_out_xy', &
              'rad_lw_in_xz', 'rad_lw_out_xz', 'rad_sw_in_xz','rad_sw_out_xz', &
              'rad_lw_in_yz', 'rad_lw_out_yz', 'rad_sw_in_yz','rad_sw_out_yz' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'


       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

        END SELECT

    END SUBROUTINE radiation_define_netcdf_grid

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_2d( av, variable, found, grid, mode,         &
                                      local_pf, two_d )
 
    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !< 
    CHARACTER (LEN=*) ::  mode     !< 
    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av !< 
    INTEGER(iwp) ::  i  !< 
    INTEGER(iwp) ::  j  !< 
    INTEGER(iwp) ::  k  !< 

    LOGICAL      ::  found !< 
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxlg:nxrg,nysg:nyng,nzb:nzt+1) ::  local_pf !< 

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'rad_net*_xy' )        ! 2d-array
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   local_pf(i,j,nzb+1) = rad_net(j,i)
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   local_pf(i,j,nzb+1) = rad_net_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

 
       CASE ( 'rad_lw_in_xy', 'rad_lw_in_xz', 'rad_lw_in_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_lw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_lw_in_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_lw_out_xy', 'rad_lw_out_xz', 'rad_lw_out_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_lw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_lw_out_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF   
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_lw_cs_hr_xy', 'rad_lw_cs_hr_xz', 'rad_lw_cs_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_lw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_lw_cs_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_lw_hr_xy', 'rad_lw_hr_xz', 'rad_lw_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_lw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_lw_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_sw_in_xy', 'rad_sw_in_xz', 'rad_sw_in_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_in_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_sw_out_xy', 'rad_sw_out_xz', 'rad_sw_out_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_out_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_sw_cs_hr_xy', 'rad_sw_cs_hr_xz', 'rad_sw_cs_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_cs_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_sw_hr_xy', 'rad_sw_hr_xz', 'rad_sw_hr_yz' )
          IF ( av == 0 ) THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng 
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT
 
 END SUBROUTINE radiation_data_output_2d


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_3d( av, variable, found, local_pf )
 

    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable !< 

    INTEGER(iwp) ::  av    !< 
    INTEGER(iwp) ::  i     !< 
    INTEGER(iwp) ::  j     !< 
    INTEGER(iwp) ::  k     !< 

    LOGICAL      ::  found !< 

    REAL(sp), DIMENSION(nxlg:nxrg,nysg:nyng,nzb:nzt+1) ::  local_pf !< 


    found = .TRUE.


    SELECT CASE ( TRIM( variable ) )

      CASE ( 'rad_sw_in' )
         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_sw_in(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_sw_in_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_out' )
         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_sw_out(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_sw_out_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_cs_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_sw_cs_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_sw_cs_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_sw_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_sw_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_in' )
         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_lw_in(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_lw_in_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_out' )
         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_lw_out(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_lw_out_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_cs_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_lw_cs_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_lw_cs_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_lw_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = rad_lw_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE radiation_data_output_3d

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining masked data output
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_mask( av, variable, found, local_pf )
 
    USE control_parameters
        
    USE indices
    
    USE kinds
    

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable   !< 

    INTEGER(iwp) ::  av   !< 
    INTEGER(iwp) ::  i    !< 
    INTEGER(iwp) ::  j    !< 
    INTEGER(iwp) ::  k    !< 

    LOGICAL ::  found     !< 

    REAL(wp),                                                                  &
       DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  &
          local_pf   !< 


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )


       CASE ( 'rad_lw_in' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_in(mask_k(mid,k),              &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_in_av(mask_k(mid,k),           &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_lw_out' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_out(mask_k(mid,k),             &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_out_av(mask_k(mid,k),          &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_lw_cs_hr' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_cs_hr(mask_k(mid,k),           &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_cs_hr_av(mask_k(mid,k),        &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_lw_hr' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_hr(mask_k(mid,k),              &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_lw_hr_av(mask_k(mid,k),           &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_sw_in' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_in(mask_k(mid,k),              &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_in_av(mask_k(mid,k),           &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_sw_out' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_out(mask_k(mid,k),             &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_out_av(mask_k(mid,k),          &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_sw_cs_hr' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_cs_hr(mask_k(mid,k),           &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_cs_hr_av(mask_k(mid,k),        &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'rad_sw_hr' )
          IF ( av == 0 )  THEN
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_hr(mask_k(mid,k),              &
                                            mask_j(mid,j),mask_i(mid,i))
                    ENDDO
                 ENDDO
              ENDDO
          ELSE
             DO  i = 1, mask_size_l(mid,1)
                DO  j = 1, mask_size_l(mid,2)
                   DO  k = 1, mask_size_l(mid,3)
                       local_pf(i,j,k) = rad_sw_hr_av(mask_k(mid,k),           &
                                               mask_j(mid,j),mask_i(mid,i))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE radiation_data_output_mask


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defines masked output variables
!------------------------------------------------------------------------------!
 SUBROUTINE radiation_last_actions
 

    USE control_parameters
        
    USE kinds

    IMPLICIT NONE

    IF ( write_binary(1:4) == 'true' )  THEN
       IF ( ALLOCATED( rad_net ) )  THEN
          WRITE ( 14 )  'rad_net             ';  WRITE ( 14 )  rad_net  
       ENDIF
       IF ( ALLOCATED( rad_net_av ) )  THEN
          WRITE ( 14 )  'rad_net_av          ';  WRITE ( 14 )  rad_net_av  
       ENDIF  
       IF ( ALLOCATED( rad_lw_in ) )  THEN
          WRITE ( 14 )  'rad_lw_in           ';  WRITE ( 14 )  rad_lw_in  
       ENDIF
       IF ( ALLOCATED( rad_lw_in_av ) )  THEN
          WRITE ( 14 )  'rad_lw_in_av        ';  WRITE ( 14 )  rad_lw_in_av  
       ENDIF 
       IF ( ALLOCATED( rad_lw_out ) )  THEN
          WRITE ( 14 )  'rad_lw_out          ';  WRITE ( 14 )  rad_lw_out 
       ENDIF 
       IF ( ALLOCATED( rad_lw_out_av ) )  THEN
          WRITE ( 14 )  'rad_lw_out_av       ';  WRITE ( 14 )  rad_lw_out_av  
       ENDIF 
       IF ( ALLOCATED( rad_lw_out_change_0 ) )  THEN
          WRITE ( 14 )  'rad_lw_out_change_0 '
          WRITE ( 14 )  rad_lw_out_change_0
       ENDIF
       IF ( ALLOCATED( rad_lw_cs_hr ) )  THEN
          WRITE ( 14 )  'rad_lw_cs_hr        ';  WRITE ( 14 )  rad_lw_cs_hr
       ENDIF
       IF ( ALLOCATED( rad_lw_cs_hr_av ) )  THEN
          WRITE ( 14 )  'rad_lw_cs_hr_av     ';  WRITE ( 14 )  rad_lw_cs_hr_av
       ENDIF
       IF ( ALLOCATED( rad_lw_hr ) )  THEN
          WRITE ( 14 )  'rad_lw_hr           ';  WRITE ( 14 )  rad_lw_hr
       ENDIF
       IF ( ALLOCATED( rad_lw_hr_av ) )  THEN
          WRITE ( 14 )  'rad_lw_hr_av        ';  WRITE ( 14 )  rad_lw_hr_av
       ENDIF
       IF ( ALLOCATED( rad_sw_in ) )  THEN
          WRITE ( 14 )  'rad_sw_in           ';  WRITE ( 14 )  rad_sw_in  
       ENDIF 
       IF ( ALLOCATED( rad_sw_in_av ) )  THEN
          WRITE ( 14 )  'rad_sw_in_av        ';  WRITE ( 14 )  rad_sw_in_av  
       ENDIF 
       IF ( ALLOCATED( rad_sw_out ) )  THEN
          WRITE ( 14 )  'rad_sw_out          ';  WRITE ( 14 )  rad_sw_out  
       ENDIF 
       IF ( ALLOCATED( rad_sw_out_av ) )  THEN
          WRITE ( 14 )  'rad_sw_out_av       ';  WRITE ( 14 )  rad_sw_out_av  
       ENDIF 
       IF ( ALLOCATED( rad_sw_cs_hr ) )  THEN
          WRITE ( 14 )  'rad_sw_cs_hr        ';  WRITE ( 14 )  rad_sw_cs_hr
       ENDIF
       IF ( ALLOCATED( rad_sw_cs_hr_av ) )  THEN
          WRITE ( 14 )  'rad_sw_cs_hr_av     ';  WRITE ( 14 )  rad_sw_cs_hr_av
       ENDIF
       IF ( ALLOCATED( rad_sw_hr ) )  THEN
          WRITE ( 14 )  'rad_sw_hr           ';  WRITE ( 14 )  rad_sw_hr
       ENDIF
       IF ( ALLOCATED( rad_sw_hr_av ) )  THEN
          WRITE ( 14 )  'rad_sw_hr_av        ';  WRITE ( 14 )  rad_sw_hr_av
       ENDIF

       WRITE ( 14 )  '*** end rad ***     '

    ENDIF

 END SUBROUTINE radiation_last_actions


SUBROUTINE radiation_read_restart_data( i, nxlfa, nxl_on_file, nxrfa, nxr_on_file,   &
                                     nynfa, nyn_on_file, nysfa, nys_on_file,   &
                                     offset_xa, offset_ya, overlap_count,      &
                                     tmp_2d, tmp_3d )
 

    USE control_parameters
        
    USE indices
    
    USE kinds
    
    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=20) :: field_char   !< 

    INTEGER(iwp) ::  i               !< 
    INTEGER(iwp) ::  k               !< 
    INTEGER(iwp) ::  nxlc            !< 
    INTEGER(iwp) ::  nxlf            !< 
    INTEGER(iwp) ::  nxl_on_file     !< 
    INTEGER(iwp) ::  nxrc            !< 
    INTEGER(iwp) ::  nxrf            !< 
    INTEGER(iwp) ::  nxr_on_file     !< 
    INTEGER(iwp) ::  nync            !< 
    INTEGER(iwp) ::  nynf            !< 
    INTEGER(iwp) ::  nyn_on_file     !< 
    INTEGER(iwp) ::  nysc            !< 
    INTEGER(iwp) ::  nysf            !< 
    INTEGER(iwp) ::  nys_on_file     !< 
    INTEGER(iwp) ::  overlap_count   !< 

    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxlfa       !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxrfa       !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nynfa       !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nysfa       !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_xa   !< 
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_ya   !< 

    REAL(wp),                                                                  &
       DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::&
          tmp_2d   !< 

    REAL(wp),                                                                  &
       DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::&
          tmp_3d   !< 



   IF ( initializing_actions == 'read_restart_data' )  THEN
      READ ( 13 )  field_char

      DO  WHILE ( TRIM( field_char ) /= '*** end rad ***' )

         DO  k = 1, overlap_count

            nxlf = nxlfa(i,k)
            nxlc = nxlfa(i,k) + offset_xa(i,k)
            nxrf = nxrfa(i,k)
            nxrc = nxrfa(i,k) + offset_xa(i,k)
            nysf = nysfa(i,k)
            nysc = nysfa(i,k) + offset_ya(i,k)
            nynf = nynfa(i,k)
            nync = nynfa(i,k) + offset_ya(i,k)


            SELECT CASE ( TRIM( field_char ) )

                CASE ( 'rad_net' )
                   IF ( .NOT. ALLOCATED( rad_net ) )  THEN
                      ALLOCATE( rad_net(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   rad_net(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_net_av' )
                   IF ( .NOT. ALLOCATED( rad_net_av ) )  THEN
                      ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   rad_net_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  = &
                                          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
                CASE ( 'rad_lw_in' )
                   IF ( .NOT. ALLOCATED( rad_lw_in ) )  THEN
                      ALLOCATE( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_lw_in(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_lw_in_av' )
                   IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
                      ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_lw_in_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_lw_out' )
                   IF ( .NOT. ALLOCATED( rad_lw_out ) )  THEN
                      ALLOCATE( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_lw_out(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_lw_out_av' )
                   IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
                      ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_lw_out_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_lw_out_change_0' )
                   IF ( .NOT. ALLOCATED( rad_lw_out_change_0 ) )  THEN
                      ALLOCATE( rad_lw_out_change_0(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   rad_lw_out_change_0(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)&
                              = tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_lw_cs_hr' )
                   IF ( .NOT. ALLOCATED( rad_lw_cs_hr ) )  THEN
                      ALLOCATE( rad_lw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_lw_cs_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                           tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_lw_cs_hr_av' )
                   IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
                      ALLOCATE( rad_lw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_lw_cs_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                           tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_lw_hr' )
                   IF ( .NOT. ALLOCATED( rad_lw_hr ) )  THEN
                      ALLOCATE( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_lw_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                           tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_lw_hr_av' )
                   IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
                      ALLOCATE( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_lw_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                           tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_sw_in' )
                   IF ( .NOT. ALLOCATED( rad_sw_in ) )  THEN
                      ALLOCATE( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_sw_in(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_sw_in_av' )
                   IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
                      ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_sw_in_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_sw_out' )
                   IF ( .NOT. ALLOCATED( rad_sw_out ) )  THEN
                      ALLOCATE( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_sw_out(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_sw_out_av' )
                   IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
                      ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF  
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_sw_out_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_sw_cs_hr' )
                   IF ( .NOT. ALLOCATED( rad_sw_cs_hr ) )  THEN
                      ALLOCATE( rad_sw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_sw_cs_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                           tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_sw_cs_hr_av' )
                   IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
                      ALLOCATE( rad_sw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_sw_cs_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                           tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_sw_hr' )
                   IF ( .NOT. ALLOCATED( rad_sw_hr ) )  THEN
                      ALLOCATE( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_sw_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                           tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'rad_sw_hr_av' )
                   IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
                      ALLOCATE( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rad_lw_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) = &
                           tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

               CASE DEFAULT
                  WRITE( message_string, * ) 'unknown variable named "',       &
                                        TRIM( field_char ), '" found in',      &
                                        '&data from prior run on PE ', myid
                  CALL message( 'radiation_read_restart_data', 'PA0441', 1, 2, 0, 6, &
                                 0 )

            END SELECT

         ENDDO

         READ ( 13 )  field_char

      ENDDO
   ENDIF

 END SUBROUTINE radiation_read_restart_data


 END MODULE radiation_model_mod
