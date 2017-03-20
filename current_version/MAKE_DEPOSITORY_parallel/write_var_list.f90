!> @file write_var_list.f90
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
! $Id: write_var_list.f90 2043 2016-11-02 13:48:37Z suehring $
!
! 2042 2016-11-02 13:47:31Z suehring
! Bugfix, write restart data for wall_heatflux, wall_qflux and wall_sflux
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1992 2016-08-12 15:14:59Z suehring
! top scalarflux added 
!
! 1957 2016-07-07 10:43:48Z suehring
! flight module added
! 
! 1849 2016-04-08 11:33:18Z hoffmann
! Adapted for modularization of microphysics
!
! 1833 2016-04-07 14:23:03Z raasch
! spectra_mod added
!
! 1831 2016-04-07 13:15:51Z hoffmann 
! turbulence renamed collision_turbulence, drizzle renamed 
! cloud_water_sedimentation
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module name changed + related changes
!
! 1705 2015-11-02 14:28:56Z maronga
! Bugfix: two lines required swapping
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added output of most_method, constant_flux_layer, zeta_min, zeta_max. Removed 
! output of prandtl_layer and rif_min, rif_max.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Adapted for RRTMG
! 
! 1551 2015-03-03 14:18:16Z maronga
! Typo removed
! 
! 1502 2014-12-03 18:22:31Z kanani
! Canopy module and parameters removed (parameters are always read from 
! canopy_par NAMELIST for initial and restart runs),
! Bugfix: added blanks in "cloud_top_radiation"-string to a total of 30 
! characters
! 
! 1496 2014-12-02 17:25:50Z maronga
! Renamed "radiation" -> "cloud_top_radiation"
! 
! 1484 2014-10-21 10:53:05Z kanani
! Changes in the course of the canopy-model modularization:
!   parameters alpha_lad, beta_lad, lai_beta added,
!   module plant_canopy_model_mod added,
!   drag_coefficient, leaf_surface_concentration and scalar_exchange_coefficient
!   renamed to canopy_drag_coeff, leaf_surface_conc and leaf_scalar_exch_coeff
!
! 1324 2014-03-21 09:13:16Z suehring
! Bugfix: ONLY statement for module netcdf_control removed 
!
! 1320 2014-03-20 08:40:49Z raasch
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
!
! 1308 2014-03-13 14:58:42Z fricke
! +do2d_xy_time_count, do2d_xz_time_count, do2d_yz_time_count,
! +do3d_time_count
!
! 1241 2013-10-30 11:36:58Z heinze
! +nudging
! +large_scale_forcing
!
! 1179 2013-06-14 05:57:58Z raasch
! +reference_state, ref_state
!
! 1115 2013-03-26 18:16:16Z hoffmann
! unused variables removed
!
! 1065 2012-11-22 17:42:36Z hoffmann
! +nc, c_sedimentation, turbulence, limiter_sedimentation
! -mu_constant, mu_constant_value
!
! 1053 2012-11-13 17:11:03Z hoffmann
! necessary expansions according to the two new prognostic equations (nr, qr) 
! of the two-moment cloud physics scheme:
! +bc_*_b, bc_*_t, bc_*_t_val, *_init, *_surface, *_surface_initial_change,
! +*_vertical_gradient, *_vertical_gradient_level, *_vertical_gradient_level_ind,
! +surface_waterflux_*
!
! in addition, steering parameters parameters of the two-moment cloud physics 
! scheme:   
! +cloud_scheme, +drizzle, +mu_constant, +mu_constant_value, +ventilation_effect
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! -adjust_mixing_length
!
! 1003 2012-09-14 14:35:53Z raasch
! -grid_matching
!
! 1001 2012-09-13 14:08:46Z raasch
! -cut_spline_overshoot, last_dt_change, long_filter_factor, overshoot_limit_*,
! ups_limit_*
!
! 978 2012-08-09 08:28:32Z fricke
! -km_damp_max, outflow_damping_width
! +pt_damping_factor, pt_damping_width
! +z0h_factor
!
! 940 2012-07-09 14:31:00Z raasch
! +neutral
!
! 927 2012-06-06 19:15:04Z raasch
! +masking_method
!
! 849 2012-03-15 10:35:09Z raasch
! first_call_advec_particles renamed first_call_lpm
!
! 824 2012-02-17 09:09:57Z raasch
! +curvature_solution_effects
!
! Revision 1.1  1998/03/18 20:20:38  raasch
! Initial revision
!
!
! Description:
! ------------
!> Writing values of control variables to restart-file (binary format).
!> This information are only written to the file opened by PE0.
!------------------------------------------------------------------------------!
 SUBROUTINE write_var_list
 

    USE arrays_3d,                                                             &
        ONLY:  inflow_damping_factor, mean_inflow_profiles, pt_init,           &
               q_init, ref_state, s_init, sa_init, u_init, ug, v_init, vg

    USE control_parameters
    
    USE flight_mod,                                                            &
        ONLY:  flight_write_restart_data
    
    USE grid_variables,                                                        &
        ONLY:  dx, dy
    
    USE indices,                                                               &
        ONLY:  nz, nx, ny

    USE microphysics_mod,                                                      &
        ONLY:  c_sedimentation, cloud_water_sedimentation,                     &
               collision_turbulence, limiter_sedimentation, nc_const,          &
               ventilation_effect

    USE model_1d,                                                              &
        ONLY:  damp_level_1d, dt_pr_1d, dt_run_control_1d, end_time_1d
    
    USE netcdf_interface,                                                      &
        ONLY:  netcdf_precision, output_for_t0
    
    USE particle_attributes,                                                   &
        ONLY:  curvature_solution_effects, time_sort_particles
    
    USE pegrid

    USE radiation_model_mod,                                                   &
        ONLY:  time_radiation

    USE spectra_mod,                                                           &
        ONLY:  average_count_sp

    USE statistics,                                                            &
        ONLY:  statistic_regions, hom, hom_sum, u_max, u_max_ijk, v_max,       &
               v_max_ijk, w_max, w_max_ijk

    
    IMPLICIT NONE

    CHARACTER (LEN=10) ::  binary_version   !< 


    binary_version = '4.1'

    WRITE ( 14 )  binary_version

    WRITE ( 14 )  'numprocs                      '
    WRITE ( 14 )  numprocs
    WRITE ( 14 )  'hor_index_bounds              '
    WRITE ( 14 )  hor_index_bounds
    WRITE ( 14 )  'nz                            '
    WRITE ( 14 )  nz
    WRITE ( 14 )  'max_pr_user                   '
    WRITE ( 14 )  max_pr_user
    WRITE ( 14 )  'statistic_regions             '
    WRITE ( 14 )  statistic_regions

!
!-- Caution: After changes in the following parameter-list, the
!-- -------  version number stored in the variable binary_version has to be
!--          increased. The same changes must also be done in the parameter-
!--          list in read_var_list.

    WRITE ( 14 )  'advected_distance_x           '
    WRITE ( 14 )  advected_distance_x
    WRITE ( 14 )  'advected_distance_y           '
    WRITE ( 14 )  advected_distance_y
    WRITE ( 14 )  'alpha_surface                 '
    WRITE ( 14 )  alpha_surface
    WRITE ( 14 )  'average_count_pr              '
    WRITE ( 14 )  average_count_pr
    WRITE ( 14 )  'average_count_sp              '
    WRITE ( 14 )  average_count_sp
    WRITE ( 14 )  'average_count_3d              '
    WRITE ( 14 )  average_count_3d
    WRITE ( 14 )  'bc_e_b                        '
    WRITE ( 14 )  bc_e_b
    WRITE ( 14 )  'bc_lr                         '
    WRITE ( 14 )  bc_lr
    WRITE ( 14 )  'bc_ns                         '
    WRITE ( 14 )  bc_ns
    WRITE ( 14 )  'bc_p_b                        '
    WRITE ( 14 )  bc_p_b
    WRITE ( 14 )  'bc_p_t                        '
    WRITE ( 14 )  bc_p_t
    WRITE ( 14 )  'bc_pt_b                       '
    WRITE ( 14 )  bc_pt_b
    WRITE ( 14 )  'bc_pt_t                       '
    WRITE ( 14 )  bc_pt_t
    WRITE ( 14 )  'bc_pt_t_val                   '
    WRITE ( 14 )  bc_pt_t_val
    WRITE ( 14 )  'bc_q_b                        '
    WRITE ( 14 )  bc_q_b
    WRITE ( 14 )  'bc_q_t                        '
    WRITE ( 14 )  bc_q_t
    WRITE ( 14 )  'bc_q_t_val                    '
    WRITE ( 14 )  bc_q_t_val
    WRITE ( 14 )  'bc_s_b                        '
    WRITE ( 14 )  bc_s_b
    WRITE ( 14 )  'bc_s_t                        '
    WRITE ( 14 )  bc_s_t
    WRITE ( 14 )  'bc_sa_t                       '
    WRITE ( 14 )  bc_sa_t
    WRITE ( 14 )  'bc_uv_b                       '
    WRITE ( 14 )  bc_uv_b
    WRITE ( 14 )  'bc_uv_t                       '
    WRITE ( 14 )  bc_uv_t
    WRITE ( 14 )  'bottom_salinityflux           '
    WRITE ( 14 )  bottom_salinityflux
    WRITE ( 14 )  'building_height               '
    WRITE ( 14 )  building_height
    WRITE ( 14 )  'building_length_x             '
    WRITE ( 14 )  building_length_x
    WRITE ( 14 )  'building_length_y             '
    WRITE ( 14 )  building_length_y
    WRITE ( 14 )  'building_wall_left            '
    WRITE ( 14 )  building_wall_left
    WRITE ( 14 )  'building_wall_south           '
    WRITE ( 14 )  building_wall_south
    WRITE ( 14 )  'call_psolver_at_all_substeps  '
    WRITE ( 14 )  call_psolver_at_all_substeps
    WRITE ( 14 )  'canyon_height                 '
    WRITE ( 14 )  canyon_height
    WRITE ( 14 )  'canyon_width_x                '
    WRITE ( 14 )  canyon_width_x
    WRITE ( 14 )  'canyon_width_y                '
    WRITE ( 14 )  canyon_width_y
    WRITE ( 14 )  'canyon_wall_left              '
    WRITE ( 14 )  canyon_wall_left
    WRITE ( 14 )  'canyon_wall_south             '
    WRITE ( 14 )  canyon_wall_south
    WRITE ( 14 )  'c_sedimentation               '
    WRITE ( 14 )  c_sedimentation
    WRITE ( 14 )  'cfl_factor                    '
    WRITE ( 14 )  cfl_factor
    WRITE ( 14 )  'cloud_droplets                '
    WRITE ( 14 )  cloud_droplets
    WRITE ( 14 )  'cloud_physics                 '
    WRITE ( 14 )  cloud_physics
    WRITE ( 14 )  'cloud_scheme                  '
    WRITE ( 14 )  cloud_scheme
    WRITE ( 14 )  'collective_wait               '
    WRITE ( 14 )  collective_wait
    WRITE ( 14 )  'conserve_volume_flow          '
    WRITE ( 14 )  conserve_volume_flow
    WRITE ( 14 )  'conserve_volume_flow_mode     '
    WRITE ( 14 )  conserve_volume_flow_mode
    WRITE ( 14 )  'coupling_start_time           '
    WRITE ( 14 )  coupling_start_time
    WRITE ( 14 )  'constant_flux_layer           '
    WRITE ( 14 )  constant_flux_layer
    WRITE ( 14 )  'current_timestep_number       '
    WRITE ( 14 )  current_timestep_number
    WRITE ( 14 )  'curvature_solution_effects    '
    WRITE ( 14 )  curvature_solution_effects
    WRITE ( 14 )  'cycle_mg                      '
    WRITE ( 14 )  cycle_mg
    WRITE ( 14 )  'damp_level_1d                 '
    WRITE ( 14 )  damp_level_1d
    WRITE ( 14 )  'dissipation_1d                '
    WRITE ( 14 )  dissipation_1d
    WRITE ( 14 )  'do2d_xy_time_count            '
    WRITE ( 14 )  do2d_xy_time_count
    WRITE ( 14 )  'do2d_xz_time_count            '
    WRITE ( 14 )  do2d_xz_time_count
    WRITE ( 14 )  'do2d_yz_time_count            '
    WRITE ( 14 )  do2d_yz_time_count
    WRITE ( 14 )  'do3d_time_count               '
    WRITE ( 14 )  do3d_time_count
    WRITE ( 14 )  'dp_external                   '
    WRITE ( 14 )  dp_external
    WRITE ( 14 )  'dp_level_b                    '
    WRITE ( 14 )  dp_level_b
    WRITE ( 14 )  'dp_smooth                     '
    WRITE ( 14 )  dp_smooth
    WRITE ( 14 )  'dpdxy                         '
    WRITE ( 14 )  dpdxy
    WRITE ( 14 )  'cloud_water_sedimentation     '
    WRITE ( 14 )  cloud_water_sedimentation
    WRITE ( 14 )  'dt_pr_1d                      '
    WRITE ( 14 )  dt_pr_1d
    WRITE ( 14 )  'dt_run_control_1d             '
    WRITE ( 14 )  dt_run_control_1d
    WRITE ( 14 )  'dt_3d                         '
    WRITE ( 14 )  dt_3d
    WRITE ( 14 )  'dvrp_filecount                '
    WRITE ( 14 )  dvrp_filecount
    WRITE ( 14 )  'dx                            '
    WRITE ( 14 )  dx
    WRITE ( 14 )  'dy                            '
    WRITE ( 14 )  dy
    WRITE ( 14 )  'dz                            '
    WRITE ( 14 )  dz
    WRITE ( 14 )  'dz_max                        '
    WRITE ( 14 )  dz_max
    WRITE ( 14 )  'dz_stretch_factor             '
    WRITE ( 14 )  dz_stretch_factor
    WRITE ( 14 )  'dz_stretch_level              '
    WRITE ( 14 )  dz_stretch_level
    WRITE ( 14 )  'e_min                         '
    WRITE ( 14 )  e_min
    WRITE ( 14 )  'end_time_1d                   '
    WRITE ( 14 )  end_time_1d
    WRITE ( 14 )  'fft_method                    '
    WRITE ( 14 )  fft_method
    WRITE ( 14 )  'first_call_lpm                '
    WRITE ( 14 )  first_call_lpm
    WRITE ( 14 )  'galilei_transformation        '
    WRITE ( 14 )  galilei_transformation
    WRITE ( 14 )  'hom                           '
    WRITE ( 14 )  hom
    WRITE ( 14 )  'hom_sum                       '
    WRITE ( 14 )  hom_sum
    WRITE ( 14 )  'humidity                      '
    WRITE ( 14 )  humidity
    IF ( ALLOCATED( inflow_damping_factor ) )  THEN
       WRITE ( 14 )  'inflow_damping_factor         '
       WRITE ( 14 )  inflow_damping_factor
    ENDIF
    WRITE ( 14 )  'inflow_damping_height         '
    WRITE ( 14 )  inflow_damping_height
    WRITE ( 14 )  'inflow_damping_width          '
    WRITE ( 14 )  inflow_damping_width
    WRITE ( 14 )  'inflow_disturbance_begin      '
    WRITE ( 14 )  inflow_disturbance_begin
    WRITE ( 14 )  'inflow_disturbance_end        '
    WRITE ( 14 )  inflow_disturbance_end
    WRITE ( 14 )  'km_constant                   '
    WRITE ( 14 )  km_constant
    WRITE ( 14 )  'large_scale_forcing           '
    WRITE ( 14 )  large_scale_forcing
    WRITE ( 14 )  'large_scale_subsidence        '
    WRITE ( 14 )  large_scale_subsidence
    WRITE ( 14 )  'limiter_sedimentation         '
    WRITE ( 14 )  limiter_sedimentation
    WRITE ( 14 )  'loop_optimization             '
    WRITE ( 14 )  loop_optimization
    WRITE ( 14 )  'masking_method                '
    WRITE ( 14 )  masking_method
    IF ( ALLOCATED( mean_inflow_profiles ) )  THEN
       WRITE ( 14 )  'mean_inflow_profiles          '
       WRITE ( 14 )  mean_inflow_profiles
    ENDIF
    WRITE ( 14 )  'mg_cycles                     '
    WRITE ( 14 )  mg_cycles
    WRITE ( 14 )  'mg_switch_to_pe0_level        '
    WRITE ( 14 )  mg_switch_to_pe0_level
    WRITE ( 14 )  'mixing_length_1d              '
    WRITE ( 14 )  mixing_length_1d
    WRITE ( 14 )  'momentum_advec                '
    WRITE ( 14 )  momentum_advec
    WRITE ( 14 )  'most_method                   '
    WRITE ( 14 )  most_method
    WRITE ( 14 )  'nc_const                      '
    WRITE ( 14 )  nc_const
    WRITE ( 14 )  'netcdf_precision              '
    WRITE ( 14 )  netcdf_precision
    WRITE ( 14 )  'neutral                       '
    WRITE ( 14 )  neutral
    WRITE ( 14 )  'ngsrb                         '
    WRITE ( 14 )  ngsrb
    WRITE ( 14 )  'nsor                          '
    WRITE ( 14 )  nsor
    WRITE ( 14 )  'nsor_ini                      '
    WRITE ( 14 )  nsor_ini
    WRITE ( 14 )  'nudging                       '
    WRITE ( 14 )  nudging
    WRITE ( 14 )  'num_leg                       '
    WRITE ( 14 )  num_leg
    WRITE ( 14 )  'nx                            '
    WRITE ( 14 )  nx
    WRITE ( 14 )  'ny                            '
    WRITE ( 14 )  ny
    WRITE ( 14 )  'ocean                         '
    WRITE ( 14 )  ocean
    WRITE ( 14 )  'old_dt                        '
    WRITE ( 14 )  old_dt
    WRITE ( 14 )  'omega                         '
    WRITE ( 14 )  omega
    WRITE ( 14 )  'omega_sor                     '
    WRITE ( 14 )  omega_sor
    WRITE ( 14 )  'output_for_t0                 '
    WRITE ( 14 )  output_for_t0
    WRITE ( 14 )  'passive_scalar                '
    WRITE ( 14 )  passive_scalar
    WRITE ( 14 )  'phi                           '
    WRITE ( 14 )  phi
    WRITE ( 14 )  'prandtl_number                '
    WRITE ( 14 )  prandtl_number
    WRITE ( 14 )  'precipitation                 '
    WRITE ( 14 )  precipitation
    WRITE ( 14 )  'psolver                       '
    WRITE ( 14 )  psolver
    WRITE ( 14 )  'pt_damping_factor             '
    WRITE ( 14 )  pt_damping_factor
    WRITE ( 14 )  'pt_damping_width              '
    WRITE ( 14 )  pt_damping_width
    WRITE ( 14 )  'pt_init                       '
    WRITE ( 14 )  pt_init
    WRITE ( 14 )  'pt_reference                  '
    WRITE ( 14 )  pt_reference
    WRITE ( 14 )  'pt_surface                    '
    WRITE ( 14 )  pt_surface
    WRITE ( 14 )  'pt_surface_initial_change     '
    WRITE ( 14 )  pt_surface_initial_change
    WRITE ( 14 )  'pt_vertical_gradient          '
    WRITE ( 14 )  pt_vertical_gradient
    WRITE ( 14 )  'pt_vertical_gradient_level    '
    WRITE ( 14 )  pt_vertical_gradient_level
    WRITE ( 14 )  'pt_vertical_gradient_level_ind'
    WRITE ( 14 )  pt_vertical_gradient_level_ind
    WRITE ( 14 )  'q_init                        '
    WRITE ( 14 )  q_init
    WRITE ( 14 )  'q_surface                     '
    WRITE ( 14 )  q_surface
    WRITE ( 14 )  'q_surface_initial_change      '
    WRITE ( 14 )  q_surface_initial_change
    WRITE ( 14 )  'q_vertical_gradient           '
    WRITE ( 14 )  q_vertical_gradient
    WRITE ( 14 )  'q_vertical_gradient_level     '
    WRITE ( 14 )  q_vertical_gradient_level
    WRITE ( 14 )  'q_vertical_gradient_level_ind '
    WRITE ( 14 )  q_vertical_gradient_level_ind
    WRITE ( 14 )  'cloud_top_radiation           '
    WRITE ( 14 )  cloud_top_radiation
    WRITE ( 14 )  'random_generator              '
    WRITE ( 14 )  random_generator
    WRITE ( 14 )  'random_heatflux               '
    WRITE ( 14 )  random_heatflux
    WRITE ( 14 )  'rayleigh_damping_factor       '
    WRITE ( 14 )  rayleigh_damping_factor
    WRITE ( 14 )  'rayleigh_damping_height       '
    WRITE ( 14 )  rayleigh_damping_height
    WRITE ( 14 )  'recycling_width               '
    WRITE ( 14 )  recycling_width
    WRITE ( 14 )  'reference_state               '
    WRITE ( 14 )  reference_state
    WRITE ( 14 )  'ref_state                     '
    WRITE ( 14 )  ref_state
    WRITE ( 14 )  'residual_limit                '
    WRITE ( 14 )  residual_limit
    WRITE ( 14 )  'roughness_length              '
    WRITE ( 14 )  roughness_length
    WRITE ( 14 )  'runnr                         '
    WRITE ( 14 )  runnr
    WRITE ( 14 )  'run_coupled                   '
    WRITE ( 14 )  run_coupled
    WRITE ( 14 )  's_init                        '
    WRITE ( 14 )  s_init
    WRITE ( 14 )  's_surface                     '
    WRITE ( 14 )  s_surface
    WRITE ( 14 )  's_surface_initial_change      '
    WRITE ( 14 )  s_surface_initial_change
    WRITE ( 14 )  's_vertical_gradient           '
    WRITE ( 14 )  s_vertical_gradient
    WRITE ( 14 )  's_vertical_gradient_level     '
    WRITE ( 14 )  s_vertical_gradient_level
    WRITE ( 14 )  's_vertical_gradient_level_ind '
    WRITE ( 14 )  s_vertical_gradient_level_ind 
    WRITE ( 14 )  'sa_init                       '
    WRITE ( 14 )  sa_init
    WRITE ( 14 )  'sa_surface                    '
    WRITE ( 14 )  sa_surface
    WRITE ( 14 )  'sa_vertical_gradient          '
    WRITE ( 14 )  sa_vertical_gradient
    WRITE ( 14 )  'sa_vertical_gradient_level    '
    WRITE ( 14 )  sa_vertical_gradient_level
    WRITE ( 14 )  'scalar_advec                  '
    WRITE ( 14 )  scalar_advec
    WRITE ( 14 )  'simulated_time                '
    WRITE ( 14 )  simulated_time
    WRITE ( 14 )  'surface_heatflux              '
    WRITE ( 14 )  surface_heatflux
    WRITE ( 14 )  'surface_pressure              '
    WRITE ( 14 )  surface_pressure
    WRITE ( 14 )  'surface_scalarflux            '
    WRITE ( 14 )  surface_scalarflux    
    WRITE ( 14 )  'surface_waterflux             '
    WRITE ( 14 )  surface_waterflux    
    WRITE ( 14 )  's_surface                     '
    WRITE ( 14 )  s_surface
    WRITE ( 14 )  's_surface_initial_change      '
    WRITE ( 14 )  s_surface_initial_change
    WRITE ( 14 )  's_vertical_gradient           '
    WRITE ( 14 )  s_vertical_gradient
    WRITE ( 14 )  's_vertical_gradient_level     '
    WRITE ( 14 )  s_vertical_gradient_level
    WRITE ( 14 )  'time_coupling                 '
    WRITE ( 14 )  time_coupling
    WRITE ( 14 )  'time_disturb                  '
    WRITE ( 14 )  time_disturb
    WRITE ( 14 )  'time_domask                   '
    WRITE ( 14 )  time_domask
    WRITE ( 14 )  'time_dopr                     '
    WRITE ( 14 )  time_dopr
    WRITE ( 14 )  'time_dopr_av                  '
    WRITE ( 14 )  time_dopr_av
    WRITE ( 14 )  'time_dopr_listing             '
    WRITE ( 14 )  time_dopr_listing
    WRITE ( 14 )  'time_dopts                    '
    WRITE ( 14 )  time_dopts
    WRITE ( 14 )  'time_dosp                     '
    WRITE ( 14 )  time_dosp
    WRITE ( 14 )  'time_dots                     '
    WRITE ( 14 )  time_dots
    WRITE ( 14 )  'time_do2d_xy                  '
    WRITE ( 14 )  time_do2d_xy
    WRITE ( 14 )  'time_do2d_xz                  '
    WRITE ( 14 )  time_do2d_xz
    WRITE ( 14 )  'time_do2d_yz                  '
    WRITE ( 14 )  time_do2d_yz
    WRITE ( 14 )  'time_do3d                     '
    WRITE ( 14 )  time_do3d
    WRITE ( 14 )  'time_do_av                    '
    WRITE ( 14 )  time_do_av
    WRITE ( 14 )  'time_do_sla                   '
    WRITE ( 14 )  time_do_sla
    WRITE ( 14 )  'time_dvrp                     '
    WRITE ( 14 )  time_dvrp
    WRITE ( 14 )  'time_radiation                '
    WRITE ( 14 )  time_radiation
    WRITE ( 14 )  'time_restart                  '
    WRITE ( 14 )  time_restart
    WRITE ( 14 )  'time_run_control              '
    WRITE ( 14 )  time_run_control
    WRITE ( 14 )  'time_since_reference_point    '
    WRITE ( 14 )  time_since_reference_point
    WRITE ( 14 )  'time_sort_particles           '
    WRITE ( 14 )  time_sort_particles
    WRITE ( 14 )  'timestep_scheme               '
    WRITE ( 14 )  timestep_scheme
    WRITE ( 14 )  'topography                    '
    WRITE ( 14 )  topography
    WRITE ( 14 )  'topography_grid_convention    '
    WRITE ( 14 )  topography_grid_convention
    WRITE ( 14 )  'top_heatflux                  '
    WRITE ( 14 )  top_heatflux
    WRITE ( 14 )  'top_momentumflux_u            '
    WRITE ( 14 )  top_momentumflux_u
    WRITE ( 14 )  'top_momentumflux_v            '
    WRITE ( 14 )  top_momentumflux_v
    WRITE ( 14 )  'top_salinityflux              '
    WRITE ( 14 )  top_salinityflux
    WRITE ( 14 )  'top_scalarflux                '
    WRITE ( 14 )  top_scalarflux
    WRITE ( 14 )  'tsc                           '
    WRITE ( 14 )  tsc
    WRITE ( 14 )  'collision_turbulence          '
    WRITE ( 14 )  collision_turbulence
    WRITE ( 14 )  'turbulent_inflow              '
    WRITE ( 14 )  turbulent_inflow
    WRITE ( 14 )  'u_bulk                        '
    WRITE ( 14 )  u_bulk
    WRITE ( 14 )  'u_init                        '
    WRITE ( 14 )  u_init
    WRITE ( 14 )  'u_max                         '
    WRITE ( 14 )  u_max
    WRITE ( 14 )  'u_max_ijk                     '
    WRITE ( 14 )  u_max_ijk
    WRITE ( 14 )  'ug                            '
    WRITE ( 14 )  ug
    WRITE ( 14 )  'ug_surface                    '
    WRITE ( 14 )  ug_surface
    WRITE ( 14 )  'ug_vertical_gradient          '
    WRITE ( 14 )  ug_vertical_gradient
    WRITE ( 14 )  'ug_vertical_gradient_level    '
    WRITE ( 14 )  ug_vertical_gradient_level
    WRITE ( 14 )  'ug_vertical_gradient_level_ind'
    WRITE ( 14 )  ug_vertical_gradient_level_ind
    WRITE ( 14 )  'use_surface_fluxes            '
    WRITE ( 14 )  use_surface_fluxes
    WRITE ( 14 )  'use_top_fluxes                '
    WRITE ( 14 )  use_top_fluxes
    WRITE ( 14 )  'use_ug_for_galilei_tr         '
    WRITE ( 14 )  use_ug_for_galilei_tr
    WRITE ( 14 )  'use_upstream_for_tke          '
    WRITE ( 14 )  use_upstream_for_tke
    WRITE ( 14 )  'v_bulk                        '
    WRITE ( 14 )  v_bulk
    WRITE ( 14 )  'v_init                        '
    WRITE ( 14 )  v_init
    WRITE ( 14 )  'v_max                         '
    WRITE ( 14 )  v_max
    WRITE ( 14 )  'v_max_ijk                     '
    WRITE ( 14 )  v_max_ijk
    WRITE ( 14 )  'ventilation_effect            '
    WRITE ( 14 )  ventilation_effect
    WRITE ( 14 )  'vg                            '
    WRITE ( 14 )  vg
    WRITE ( 14 )  'vg_surface                    '
    WRITE ( 14 )  vg_surface
    WRITE ( 14 )  'vg_vertical_gradient          '
    WRITE ( 14 )  vg_vertical_gradient
    WRITE ( 14 )  'vg_vertical_gradient_level    '
    WRITE ( 14 )  vg_vertical_gradient_level
    WRITE ( 14 )  'vg_vertical_gradient_level_ind'
    WRITE ( 14 )  vg_vertical_gradient_level_ind
    WRITE ( 14 )  'virtual_flight                '
    WRITE ( 14 )  virtual_flight
    WRITE ( 14 )  'volume_flow_area              '
    WRITE ( 14 )  volume_flow_area
    WRITE ( 14 )  'volume_flow_initial           '
    WRITE ( 14 )  volume_flow_initial
    WRITE ( 14 )  'wall_adjustment               '
    WRITE ( 14 )  wall_adjustment
    WRITE ( 14 )  'subs_vertical_gradient        '
    WRITE ( 14 )  subs_vertical_gradient
    WRITE ( 14 )  'subs_vertical_gradient_level  '
    WRITE ( 14 )  subs_vertical_gradient_level
    WRITE ( 14 )  'subs_vertical_gradient_level_i'
    WRITE ( 14 )  subs_vertical_gradient_level_i
    WRITE ( 14 )  'wall_heatflux                 '
    WRITE ( 14 )  wall_heatflux
    WRITE ( 14 )  'wall_qflux                    '
    WRITE ( 14 )  wall_qflux
    WRITE ( 14 )  'wall_sflux                    '
    WRITE ( 14 )  wall_sflux
    WRITE ( 14 )  'w_max                         '
    WRITE ( 14 )  w_max
    WRITE ( 14 )  'w_max_ijk                     '
    WRITE ( 14 )  w_max_ijk
    WRITE ( 14 )  'zeta_max                      '
    WRITE ( 14 )  zeta_max
    WRITE ( 14 )  'zeta_min                      '
    WRITE ( 14 )  zeta_min
    WRITE ( 14 )  'z0h_factor                    '
    WRITE ( 14 )  z0h_factor  
!
!-- Set the end-of-file mark
    WRITE ( 14 )  '*** end ***                   '
!
!-- If required, write restart data for virtual measurements.
    IF ( virtual_flight )  CALL flight_write_restart_data  
    
    
 END SUBROUTINE write_var_list
