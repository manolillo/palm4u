!> @file time_integration.f90
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
! ------------------
! kk: Added exchange_horiz for chemical species
! FKa: Minor formatting
! 
! Former revisions:
! -----------------
! $Id: time_integration.f90 2159 2017-02-22 18:01:07Z kanani $
!
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
! 
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
! 
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters,
! removed commented CALLs of global_min_max.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added CALLs for new urban surface model
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1976 2016-07-27 13:28:04Z maronga
! Simplified calls to radiation model
! 
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! 
! 1957 2016-07-07 10:43:48Z suehring
! flight module added
! 
! 1919 2016-05-27 14:51:23Z raasch
! Initial version of purely vertical nesting introduced. 
!
! 1918 2016-05-27 14:35:57Z raasch
! determination of time step moved to the end of the time step loop,
! the first time step is now always calculated before the time step loop (i.e.
! also in case of restart runs)
!
! 1914 2016-05-26 14:44:07Z witha
! Added call for wind turbine model
!
! 1878 2016-04-19 12:30:36Z hellstea
! Synchronization for nested runs rewritten
! 
! 1853 2016-04-11 09:00:35Z maronga
! Adjusted for use with radiation_scheme = constant
!
! 1849 2016-04-08 11:33:18Z hoffmann
! Adapted for modularization of microphysics
!
! 1833 2016-04-07 14:23:03Z raasch
! spectrum renamed spectra_mod, spectra related variables moved to spectra_mod
!
! 1831 2016-04-07 13:15:51Z hoffmann
! turbulence renamed collision_turbulence
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme replaced by microphysics_*
!
! 1808 2016-04-05 19:44:00Z raasch
! output message in case unscheduled radiation calls removed
!
! 1797 2016-03-21 16:50:28Z raasch
! introduction of different datatransfer modes
!
! 1791 2016-03-11 10:41:25Z raasch
! call of pmci_update_new removed
!
! 1786 2016-03-08 05:49:27Z raasch
! +module spectrum
!
! 1783 2016-03-06 18:36:17Z raasch
! switch back of netcdf data format for mask output moved to the mask output
! routine
!
! 1781 2016-03-03 15:12:23Z raasch
! some pmc calls removed at the beginning (before timeloop),
! pmc initialization moved to the main program
!
! 1764 2016-02-28 12:45:19Z raasch
! PMC_ACTIVE flags removed,
! bugfix: nest synchronization after first call of timestep
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1736 2015-12-04 08:56:33Z raasch
! no perturbations added to total domain if energy limit has been set zero
!
! 1691 2015-10-26 16:17:44Z maronga
! Added option for spin-ups without land surface and radiation models. Moved calls
! for radiation and lan surface schemes.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1671 2015-09-25 03:29:37Z raasch
! bugfix: ghostpoint exchange for array diss in case that sgs velocities are used
! for particles
!
! 1585 2015-04-30 07:05:52Z maronga
! Moved call of radiation scheme. Added support for RRTM
! 
! 1551 2015-03-03 14:18:16Z maronga
! Added interface for different radiation schemes.
! 
! 1496 2014-12-02 17:25:50Z maronga
! Added calls for the land surface model and radiation scheme
! 
! 1402 2014-05-09 14:25:13Z raasch
! location messages modified
! 
! 1384 2014-05-02 14:31:06Z raasch
! location messages added
! 
! 1380 2014-04-28 12:40:45Z heinze
! CALL of nudge_ref added 
! bc_pt_t_val and bc_q_t_val are updated in case nudging is used
! 
! 1365 2014-04-22 15:03:56Z boeske
! Reset sums_ls_l to zero at each timestep
! +sums_ls_l
! Calculation of reference state (previously in subroutine calc_mean_profile)

! 1342 2014-03-26 17:04:47Z kanani
! REAL constants defined as wp-kind
! 
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements 
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1308 2014-03-13 14:58:42Z fricke
! +netcdf_data_format_save
! For masked data, parallel netcdf output is not tested so far, hence 
! netcdf_data_format is switched back to non-paralell output.
!
! 1276 2014-01-15 13:40:41Z heinze
! Use LSF_DATA also in case of Dirichlet bottom boundary condition for scalars
!
! 1257 2013-11-08 15:18:40Z raasch
! acc-update-host directive for timestep removed
!
! 1241 2013-10-30 11:36:58Z heinze
! Generalize calc_mean_profile for wider use
! Determine shf and qsws in dependence on data from LSF_DATA
! Determine ug and vg in dependence on data from LSF_DATA
! 1221 2013-09-10 08:59:13Z raasch
! host update of arrays before timestep is called
!
! 1179 2013-06-14 05:57:58Z raasch
! mean profiles for reference state are only calculated if required,
! small bugfix for background communication
!
! 1171 2013-05-30 11:27:45Z raasch
! split of prognostic_equations deactivated (comment lines), for the time being
!
! 1128 2013-04-12 06:19:32Z raasch
! asynchronous transfer of ghost point data realized for acc-optimized version:
! prognostic_equations are first called two times for those points required for
! the left-right and north-south exchange, respectively, and then for the
! remaining points,
! those parts requiring global communication moved from prognostic_equations to
! here
!
! 1115 2013-03-26 18:16:16Z hoffmann
! calculation of qr and nr is restricted to precipitation
!
! 1113 2013-03-10 02:48:14Z raasch
! GPU-porting of boundary conditions,
! openACC directives updated
! formal parameter removed from routine boundary_conds
!
! 1111 2013-03-08 23:54:10Z raasch
! +internal timestep counter for cpu statistics added,
! openACC directives updated
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1065 2012-11-22 17:42:36Z hoffmann
! exchange of diss (dissipation rate) in case of turbulence = .TRUE. added
!
! 1053 2012-11-13 17:11:03Z hoffmann
! exchange of ghost points for nr, qr added
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1019 2012-09-28 06:46:45Z raasch
! non-optimized version of prognostic_equations removed
!
! 1015 2012-09-27 09:23:24Z raasch
! +call of prognostic_equations_acc
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog- and upstream-spline-scheme removed
!
! 849 2012-03-15 10:35:09Z raasch
! advec_particles renamed lpm, first_call_advec_particles renamed first_call_lpm
!
! 825 2012-02-19 03:03:44Z raasch
! wang_collision_kernel renamed wang_kernel
!
! Revision 1.1  1997/08/11 06:19:04  raasch
! Initial revision
!
!
! Description:
! ------------
!> Integration in time of the model equations, statistical analysis and graphic
!> output
!------------------------------------------------------------------------------!
 SUBROUTINE time_integration
 

    USE advec_ws,                                                              &
        ONLY:  ws_statistics

    USE arrays_3d,                                                             &
        ONLY:  diss, dzu, e, e_p, nr_p, prho, pt, pt_p, pt_init, q_init, q,    &
               ql, ql_c, ql_v, ql_vp, qr_p, q_p, ref_state, rho_ocean, s, s_p, sa_p, &
               tend, u, u_p, v, vpt, v_p, w, w_p

    USE calc_mean_profile_mod,                                                 &
        ONLY:  calc_mean_profile

    USE control_parameters,                                                    &
        ONLY:  advected_distance_x, advected_distance_y, average_count_3d,     &
               averaging_interval, averaging_interval_pr,                      &
               bc_lr_cyc, bc_ns_cyc, bc_pt_t_val,                              &
               bc_q_t_val, call_psolver_at_all_substeps, cloud_droplets,       &
               cloud_physics, constant_flux_layer, constant_heatflux,          &
               create_disturbances, dopr_n, constant_diffusion, coupling_mode, &
               coupling_start_time, current_timestep_number,                   &
               disturbance_created, disturbance_energy_limit, dist_range,      &
               do_sum, dt_3d, dt_averaging_input, dt_averaging_input_pr,       &
               dt_coupling, dt_data_output_av, dt_disturb, dt_do2d_xy,         &
               dt_do2d_xz, dt_do2d_yz, dt_do3d, dt_domask,dt_dopts, dt_dopr,   &
               dt_dopr_listing, dt_dots, dt_dvrp, dt_run_control,              &
               end_time, first_call_lpm, galilei_transformation, humidity,     &
               intermediate_timestep_count,                                    &
               intermediate_timestep_count_max, large_scale_forcing,           &
               loop_optimization, lsf_surf, lsf_vert, masks,                   &
               microphysics_seifert, mid, nest_domain,                         &
               neutral, nr_timesteps_this_run, nudging,                        &
               ocean, on_device, passive_scalar,                               &
               prho_reference, pt_reference, pt_slope_offset, random_heatflux, &
               run_coupled, simulated_time, simulated_time_chr,                &
               skip_time_do2d_xy, skip_time_do2d_xz, skip_time_do2d_yz,        &
               skip_time_do3d, skip_time_domask, skip_time_dopr,               &
               skip_time_data_output_av, sloping_surface,                      &
               stop_dt, terminate_coupled, terminate_run, timestep_scheme,     &
               time_coupling, time_do2d_xy, time_do2d_xz, time_do2d_yz,        &
               time_do3d, time_domask, time_dopr, time_dopr_av,                &
               time_dopr_listing, time_dopts, time_dosp, time_dosp_av,         &
               time_dots, time_do_av, time_do_sla, time_disturb, time_dvrp,    &
               time_run_control, time_since_reference_point,                   &
               turbulent_inflow, turbulent_outflow, urban_surface,             &
               use_initial_profile_as_reference,                               &
               use_single_reference_value, u_gtrans, v_gtrans, virtual_flight, &
               ws_scheme_mom, ws_scheme_sca

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE flight_mod,                                                            &
        ONLY:  flight_measurement


    USE indices,                                                               &
        ONLY:  i_left, i_right, j_north, j_south, nbgp, nx, nxl, nxlg, nxr,    &
               nxrg, nyn, nyng, nys, nysg, nzb, nzt, nzb_u_inner, nzb_v_inner

    USE interaction_droplets_ptq_mod,                                          &
        ONLY:  interaction_droplets_ptq

    USE interfaces

#ifdef KPP_CHEM
    USE kchem_driver,                                                          &
        ONLY:  chem_species, NSPEC, use_kpp_chemistry
#endif    
    
    USE kinds

    USE land_surface_model_mod,                                                &
        ONLY:  land_surface, lsm_energy_balance, lsm_soil_model,               &
               skip_time_do_lsm

    USE ls_forcing_mod,                                                        &
        ONLY:  ls_forcing_surf, ls_forcing_vert

    USE microphysics_mod,                                                      &
        ONLY: collision_turbulence

    USE nudge_mod,                                                             &
        ONLY:  calc_tnudge, nudge_ref

    USE particle_attributes,                                                   &
        ONLY:  particle_advection, particle_advection_start,                   &
               use_sgs_for_particles, wang_kernel

    USE pegrid

    USE pmc_interface,                                                         &
        ONLY:  nested_run, nesting_mode, pmci_datatrans,                       &
               pmci_ensure_nest_mass_conservation, pmci_synchronize

    USE production_e_mod,                                                      &
        ONLY:  production_e_init

    USE progress_bar,                                                          &
        ONLY:  finish_progress_bar, output_progress_bar

    USE prognostic_equations_mod,                                              &
        ONLY:  prognostic_equations_acc, prognostic_equations_cache,           &
               prognostic_equations_vector

    USE radiation_model_mod,                                                   &
        ONLY: dt_radiation, force_radiation_call, radiation, radiation_control,&
              skip_time_do_radiation, time_radiation

    USE spectra_mod,                                                           &
        ONLY: average_count_sp, averaging_interval_sp, calc_spectra, dt_dosp,  &
              skip_time_dosp

    USE statistics,                                                            &
        ONLY:  flow_statistics_called, hom, pr_palm, sums_ls_l, u_max,         &
               u_max_ijk, v_max, v_max_ijk, w_max, w_max_ijk

    USE surface_layer_fluxes_mod,                                              &
        ONLY:  surface_layer_fluxes

    USE urban_surface_mod,                                                     &
        ONLY:  usm_material_heat_model, usm_material_model,     &
               usm_radiation, usm_surface_energy_balance              

    USE user_actions_mod,                                                      &
        ONLY:  user_actions

    USE wind_turbine_model_mod,                                                &
        ONLY:  wind_turbine, wtm_forces



    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string          !< 
    INTEGER           ::  n
    
    REAL(wp) ::  dt_3d_old  !< temporary storage of timestep to be used for
                            !< steering of run control output interval


!
!-- At beginning determine the first time step
    CALL timestep

!
!-- Synchronize the timestep in case of nested run.
    IF ( nested_run )  THEN
!
!--    Synchronization by unifying the time step.
!--    Global minimum of all time-steps is used for all.
       CALL pmci_synchronize
    ENDIF

!
!-- Determine and print out the run control quantities before the first time
!-- step of this run. For the initial run, some statistics (e.g. divergence)
!-- need to be determined first.
    IF ( simulated_time == 0.0_wp )  CALL flow_statistics
    CALL run_control

!
!-- Data exchange between coupled models in case that a call has been omitted 
!-- at the end of the previous run of a job chain.
    IF ( coupling_mode /= 'uncoupled'  .AND.  run_coupled )  THEN
!
!--    In case of model termination initiated by the local model the coupler 
!--    must not be called because this would again cause an MPI hang. 
       DO WHILE ( time_coupling >= dt_coupling  .AND.  terminate_coupled == 0 )
          CALL surface_coupler
          time_coupling = time_coupling - dt_coupling
       ENDDO
       IF (time_coupling == 0.0_wp  .AND.                                      &
           time_since_reference_point < dt_coupling )                          &
       THEN
          time_coupling = time_since_reference_point
       ENDIF
    ENDIF

#if defined( __dvrp_graphics )
!
!-- Time measurement with dvrp software  
    CALL DVRP_LOG_EVENT( 2, current_timestep_number )
#endif

    CALL location_message( 'start with time-stepping', .TRUE. )
!
!-- Start of the time loop
    DO  WHILE ( simulated_time < end_time  .AND.  .NOT. stop_dt  .AND. &
                .NOT. terminate_run )

       CALL cpu_log( log_point_s(10), 'timesteps', 'start' )

!
!--    Determine ug, vg and w_subs in dependence on data from external file 
!--    LSF_DATA
       IF ( large_scale_forcing .AND. lsf_vert )  THEN
           CALL ls_forcing_vert ( simulated_time )
           sums_ls_l = 0.0_wp
       ENDIF

!
!--    Set pt_init and q_init to the current profiles taken from 
!--    NUDGING_DATA 
       IF ( nudging )  THEN
           CALL nudge_ref ( simulated_time )
!
!--        Store temperature gradient at the top boundary for possible Neumann
!--        boundary condition
           bc_pt_t_val = ( pt_init(nzt+1) - pt_init(nzt) ) / dzu(nzt+1)
           bc_q_t_val  = ( q_init(nzt+1) - q_init(nzt) ) / dzu(nzt+1)
       ENDIF

!
!--    Execute the user-defined actions
       CALL user_actions( 'before_timestep' )

!
!--    Calculate forces by wind turbines
       IF ( wind_turbine )  THEN

          CALL cpu_log( log_point(55), 'wind_turbine', 'start' )

          CALL wtm_forces

          CALL cpu_log( log_point(55), 'wind_turbine', 'stop' )

       ENDIF       
       
!
!--    Start of intermediate step loop
       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1

!
!--       Set the steering factors for the prognostic equations which depend
!--       on the timestep scheme
          CALL timestep_scheme_steering

!
!--       Calculate those variables needed in the tendency terms which need
!--       global communication
          IF ( .NOT. use_single_reference_value  .AND. &
               .NOT. use_initial_profile_as_reference )  THEN
!
!--          Horizontally averaged profiles to be used as reference state in
!--          buoyancy terms (WARNING: only the respective last call of
!--          calc_mean_profile defines the reference state!)
             IF ( .NOT. neutral )  THEN
                CALL calc_mean_profile( pt, 4 )
                ref_state(:)  = hom(:,1,4,0) ! this is used in the buoyancy term
             ENDIF
             IF ( ocean )  THEN
                CALL calc_mean_profile( rho_ocean, 64 )
                ref_state(:)  = hom(:,1,64,0)
             ENDIF
             IF ( humidity )  THEN
                CALL calc_mean_profile( vpt, 44 )
                ref_state(:)  = hom(:,1,44,0)
             ENDIF

          ENDIF

          IF ( .NOT. constant_diffusion )  CALL production_e_init
          IF ( ( ws_scheme_mom .OR. ws_scheme_sca )  .AND.  &
               intermediate_timestep_count == 1 )  CALL ws_statistics
!
!--       In case of nudging calculate current nudging time scale and horizontal
!--       means of u, v, pt and q
          IF ( nudging )  THEN
             CALL calc_tnudge( simulated_time )
             CALL calc_mean_profile( u, 1 )
             CALL calc_mean_profile( v, 2 )
             CALL calc_mean_profile( pt, 4 )
             CALL calc_mean_profile( q, 41 )
          ENDIF

!
!--       Solve the prognostic equations. A fast cache optimized version with
!--       only one single loop is used in case of Piascek-Williams advection
!--       scheme. NEC vector machines use a different version, because
!--       in the other versions a good vectorization is prohibited due to
!--       inlining problems.
          IF ( loop_optimization == 'cache' )  THEN
             CALL prognostic_equations_cache
          ELSEIF ( loop_optimization == 'vector' )  THEN
             CALL prognostic_equations_vector
          ELSEIF ( loop_optimization == 'acc' )  THEN
             i_left  = nxl;         i_right = nxr
             j_south = nys;         j_north = nyn
             CALL prognostic_equations_acc

!             i_left  = nxl;         i_right = nxl+nbgp-1
!             j_south = nys;         j_north = nyn
!             CALL prognostic_equations_acc
!             i_left  = nxr-nbgp+1;  i_right = nxr
!             j_south = nys;         j_north = nyn
!             CALL prognostic_equations_acc

!
!--          Exchange of ghost points (lateral boundary conditions)
             IF ( background_communication )  THEN

                CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'start' )
                
                send_receive = 'lr'
                sendrecv_in_background = .TRUE.
                req          = 0
                req_count    = 0

                IF ( numprocs == 1 )  THEN    ! workaround for single-core GPU runs
                   on_device = .TRUE.         ! to be removed after complete porting
                ELSE                          ! of ghost point exchange
                   !$acc update host( e_p, pt_p, u_p, v_p, w_p )
                ENDIF

                CALL exchange_horiz( u_p, nbgp )
                CALL exchange_horiz( v_p, nbgp )
                CALL exchange_horiz( w_p, nbgp )
                CALL exchange_horiz( pt_p, nbgp )
                IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e_p, nbgp )
                IF ( ocean )  THEN
                   CALL exchange_horiz( sa_p, nbgp )
                   CALL exchange_horiz( rho_ocean, nbgp )
                   CALL exchange_horiz( prho, nbgp )
                ENDIF
                IF ( humidity )  THEN
                   CALL exchange_horiz( q_p, nbgp )
                   IF ( cloud_physics .AND. microphysics_seifert )  THEN
                      CALL exchange_horiz( qr_p, nbgp )
                      CALL exchange_horiz( nr_p, nbgp )
                   ENDIF
                ENDIF
                IF ( cloud_droplets )  THEN
                   CALL exchange_horiz( ql, nbgp )
                   CALL exchange_horiz( ql_c, nbgp )
                   CALL exchange_horiz( ql_v, nbgp )
                   CALL exchange_horiz( ql_vp, nbgp )
                ENDIF
                IF ( wang_kernel  .OR.  collision_turbulence  .OR.             &
                     use_sgs_for_particles )  THEN
                   CALL exchange_horiz( diss, nbgp )
                ENDIF
                IF ( passive_scalar )  CALL exchange_horiz( s_p, nbgp )

                IF ( numprocs == 1 )  THEN    ! workaround for single-core GPU runs
                   on_device = .FALSE.        ! to be removed after complete porting
                ELSE                          ! of ghost point exchange
                   !$acc update device( e_p, pt_p, u_p, v_p, w_p )
                ENDIF

                sendrecv_in_background = .FALSE.

                CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'pause' )

             ENDIF

!             i_left  = nxl+nbgp;    i_right = nxr-nbgp
!             j_south = nys;         j_north = nys+nbgp-1
!             CALL prognostic_equations_acc
!             i_left  = nxl+nbgp;    i_right = nxr-nbgp
!             j_south = nyn-nbgp+1;  j_north = nyn
!             CALL prognostic_equations_acc

             IF ( background_communication )  THEN
                CALL cpu_log( log_point(41), 'exchange-horiz-wait', 'start' )
#if defined( __parallel )
                CALL MPI_WAITALL( req_count, req, wait_stat, ierr )
#endif
                CALL cpu_log( log_point(41), 'exchange-horiz-wait', 'pause' )

                CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'continue' )

                send_receive = 'ns'
                sendrecv_in_background = .TRUE.
                req          = 0
                req_count    = 0

                IF ( numprocs == 1 )  THEN    ! workaround for single-core GPU runs
                   on_device = .TRUE.         ! to be removed after complete porting
                ELSE                          ! of ghost point exchange
                   !$acc update host( e_p, pt_p, u_p, v_p, w_p )
                ENDIF

                CALL exchange_horiz( u_p, nbgp )
                CALL exchange_horiz( v_p, nbgp )
                CALL exchange_horiz( w_p, nbgp )
                CALL exchange_horiz( pt_p, nbgp )
                IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e_p, nbgp )
                IF ( ocean )  THEN
                   CALL exchange_horiz( sa_p, nbgp )
                   CALL exchange_horiz( rho_ocean, nbgp )
                  CALL exchange_horiz( prho, nbgp )
                ENDIF
                IF ( humidity )  THEN
                   CALL exchange_horiz( q_p, nbgp )
                   IF ( cloud_physics .AND. microphysics_seifert )  THEN
                      CALL exchange_horiz( qr_p, nbgp )
                      CALL exchange_horiz( nr_p, nbgp )
                   ENDIF
                ENDIF
                IF ( cloud_droplets )  THEN
                   CALL exchange_horiz( ql, nbgp )
                   CALL exchange_horiz( ql_c, nbgp )
                   CALL exchange_horiz( ql_v, nbgp )
                   CALL exchange_horiz( ql_vp, nbgp )
                ENDIF
                IF ( wang_kernel  .OR.  collision_turbulence  .OR.             &
                     use_sgs_for_particles )  THEN
                   CALL exchange_horiz( diss, nbgp )
                ENDIF
                IF ( passive_scalar )  CALL exchange_horiz( s_p, nbgp )

                IF ( numprocs == 1 )  THEN    ! workaround for single-core GPU runs
                   on_device = .FALSE.        ! to be removed after complete porting
                ELSE                          ! of ghost point exchange
                   !$acc update device( e_p, pt_p, u_p, v_p, w_p )
                ENDIF

                sendrecv_in_background = .FALSE.

                CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'stop' )

             ENDIF

!             i_left  = nxl+nbgp;    i_right = nxr-nbgp
!             j_south = nys+nbgp;    j_north = nyn-nbgp
!             CALL prognostic_equations_acc

             IF ( background_communication )  THEN
                CALL cpu_log( log_point(41), 'exchange-horiz-wait', 'continue' )
#if defined( __parallel )
                CALL MPI_WAITALL( req_count, req, wait_stat, ierr )
#endif
                send_receive = 'al'
                CALL cpu_log( log_point(41), 'exchange-horiz-wait', 'stop' )
             ENDIF

          ENDIF

!
!--       Particle transport/physics with the Lagrangian particle model
!--       (only once during intermediate steps, because it uses an Euler-step)
!--       ### particle model should be moved before prognostic_equations, in order
!--       to regard droplet interactions directly
          IF ( particle_advection  .AND.                         &
               simulated_time >= particle_advection_start  .AND. &
               intermediate_timestep_count == 1 )  THEN
             CALL lpm
             first_call_lpm = .FALSE.
          ENDIF

!
!--       Interaction of droplets with temperature and specific humidity.
!--       Droplet condensation and evaporation is calculated within
!--       advec_particles.
          IF ( cloud_droplets  .AND.  &
               intermediate_timestep_count == intermediate_timestep_count_max )&
          THEN
             CALL interaction_droplets_ptq
          ENDIF

!
!--       Exchange of ghost points (lateral boundary conditions)
          IF ( .NOT. background_communication )  THEN

             CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'start' )

             IF ( numprocs == 1 )  THEN    ! workaround for single-core GPU runs
                on_device = .TRUE.         ! to be removed after complete porting
             ELSE                          ! of ghost point exchange
                !$acc update host( e_p, pt_p, u_p, v_p, w_p )
             ENDIF

             CALL exchange_horiz( u_p, nbgp )
             CALL exchange_horiz( v_p, nbgp )
             CALL exchange_horiz( w_p, nbgp )
             CALL exchange_horiz( pt_p, nbgp )
             IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e_p, nbgp )
             IF ( ocean )  THEN
                CALL exchange_horiz( sa_p, nbgp )
                CALL exchange_horiz( rho_ocean, nbgp )
                CALL exchange_horiz( prho, nbgp )
             ENDIF
             IF ( humidity )  THEN
                CALL exchange_horiz( q_p, nbgp )
                IF ( cloud_physics .AND. microphysics_seifert )  THEN
                   CALL exchange_horiz( qr_p, nbgp )
                   CALL exchange_horiz( nr_p, nbgp )
                ENDIF
             ENDIF
             IF ( cloud_droplets )  THEN
                CALL exchange_horiz( ql, nbgp )
                CALL exchange_horiz( ql_c, nbgp )
                CALL exchange_horiz( ql_v, nbgp )
                CALL exchange_horiz( ql_vp, nbgp )
             ENDIF
             IF ( wang_kernel  .OR.  collision_turbulence  .OR.                &
                  use_sgs_for_particles )  THEN
                CALL exchange_horiz( diss, nbgp )
             ENDIF
             IF ( passive_scalar )  CALL exchange_horiz( s_p, nbgp )      

             IF ( numprocs == 1 )  THEN    ! workaround for single-core GPU runs
                on_device = .FALSE.        ! to be removed after complete porting
             ELSE                          ! of ghost point exchange
                !$acc update device( e_p, pt_p, u_p, v_p, w_p )
             ENDIF

!
!--          Ghost point exchange for chemical species
#ifdef KPP_CHEM
             IF ( use_kpp_chemistry )  THEN
                DO  n = 1, NSPEC
                   CALL exchange_horiz( chem_species(n)%conc_p, nbgp ) 
                ENDDO
             ENDIF
#endif
             CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'stop' )

          ENDIF

!
!--       Boundary conditions for the prognostic quantities (except of the
!--       velocities at the outflow in case of a non-cyclic lateral wall)
          CALL boundary_conds

!
!--       Swap the time levels in preparation for the next time step.
          CALL swap_timelevel

          IF ( nested_run )  THEN

             CALL cpu_log( log_point(60), 'nesting', 'start' )

!--          Domain nesting. The data transfer subroutines pmci_parent_datatrans
!--          and pmci_child_datatrans are called inside the wrapper
!--          subroutine pmci_datatrans according to the control parameters
!--          nesting_mode and nesting_datatransfer_mode.
!--          TO_DO: why is nesting_mode given as a parameter here?
             CALL pmci_datatrans( nesting_mode )

             IF ( TRIM( nesting_mode ) == 'two-way' .OR.                               &
                  nesting_mode == 'vertical' )  THEN
!
!--             Exchange_horiz is needed for all parent-domains after the
!--             anterpolation
                CALL exchange_horiz( u, nbgp )
                CALL exchange_horiz( v, nbgp )
                CALL exchange_horiz( w, nbgp )
                IF ( .NOT. neutral            )  CALL exchange_horiz( pt, nbgp )
                IF ( humidity                 )  CALL exchange_horiz( q, nbgp  )
                IF ( passive_scalar           )  CALL exchange_horiz( s, nbgp  )
                IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e, nbgp  )
             ENDIF
!
!--          Correct the w top-BC in nest domains to ensure mass conservation.
!--          This action must never be done for the root domain. Vertical 
!--          nesting implies mass conservation.
             IF ( nest_domain )  THEN
                CALL pmci_ensure_nest_mass_conservation
             ENDIF

             CALL cpu_log( log_point(60), 'nesting', 'stop' )

          ENDIF

!
!--       Temperature offset must be imposed at cyclic boundaries in x-direction
!--       when a sloping surface is used
          IF ( sloping_surface )  THEN
             IF ( nxl ==  0 )  pt(:,:,nxlg:nxl-1) = pt(:,:,nxlg:nxl-1) - &
                                                    pt_slope_offset
             IF ( nxr == nx )  pt(:,:,nxr+1:nxrg) = pt(:,:,nxr+1:nxrg) + &
                                                    pt_slope_offset
          ENDIF

!
!--       Impose a turbulent inflow using the recycling method
          IF ( turbulent_inflow )  CALL  inflow_turbulence

!
!--       Set values at outflow boundary using the special outflow condition
          IF ( turbulent_outflow )  CALL  outflow_turbulence

!
!--       Impose a random perturbation on the horizontal velocity field
          IF ( create_disturbances  .AND.                                      &
               ( call_psolver_at_all_substeps  .AND.                           &
               intermediate_timestep_count == intermediate_timestep_count_max )&
          .OR. ( .NOT. call_psolver_at_all_substeps  .AND.                     &
               intermediate_timestep_count == 1 ) )                            &
          THEN
             time_disturb = time_disturb + dt_3d
             IF ( time_disturb >= dt_disturb )  THEN
                !$acc update host( u, v )
                IF ( numprocs == 1 )  on_device = .FALSE.  ! workaround, remove later
                IF ( disturbance_energy_limit /= 0.0_wp  .AND.                 &
                     hom(nzb+5,1,pr_palm,0) < disturbance_energy_limit )  THEN
                   CALL disturb_field( nzb_u_inner, tend, u )
                   CALL disturb_field( nzb_v_inner, tend, v )
                ELSEIF ( .NOT. bc_lr_cyc  .OR.  .NOT. bc_ns_cyc )  THEN
!
!--                Runs with a non-cyclic lateral wall need perturbations
!--                near the inflow throughout the whole simulation
                   dist_range = 1
                   CALL disturb_field( nzb_u_inner, tend, u )
                   CALL disturb_field( nzb_v_inner, tend, v )
                   dist_range = 0
                ENDIF
                IF ( numprocs == 1 )  on_device = .TRUE.  ! workaround, remove later
                !$acc update device( u, v )
                time_disturb = time_disturb - dt_disturb
             ENDIF
          ENDIF

!
!--       Reduce the velocity divergence via the equation for perturbation
!--       pressure.
          IF ( intermediate_timestep_count == 1  .OR. &
                call_psolver_at_all_substeps )  THEN
             CALL pres
          ENDIF

!
!--       If required, compute liquid water content
          IF ( cloud_physics )  THEN
             CALL calc_liquid_water_content
             !$acc update device( ql )
          ENDIF
! 
!--       If required, compute virtual potential temperature 
          IF ( humidity )  THEN 
             CALL compute_vpt 
             !$acc update device( vpt ) 
          ENDIF 

!
!--       Compute the diffusion quantities
          IF ( .NOT. constant_diffusion )  THEN

!
!--          Determine surface fluxes shf and qsws and surface values 
!--          pt_surface and q_surface in dependence on data from external 
!--          file LSF_DATA respectively
             IF ( ( large_scale_forcing .AND. lsf_surf ) .AND. &
                 intermediate_timestep_count == intermediate_timestep_count_max )&
             THEN
                CALL ls_forcing_surf ( simulated_time )
             ENDIF

!
!--          First the vertical fluxes in the surface (constant flux) layer are computed
             IF ( constant_flux_layer )  THEN
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'start' )
                CALL surface_layer_fluxes
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'stop' )
             ENDIF

!
!--          If required, solve the energy balance for the surface and run soil 
!--          model
             IF ( land_surface .AND. simulated_time > skip_time_do_lsm)  THEN

                CALL cpu_log( log_point(54), 'land_surface', 'start' )
                CALL lsm_energy_balance
                CALL lsm_soil_model
                CALL cpu_log( log_point(54), 'land_surface', 'stop' )
             ENDIF

!
!--          If required, solve the energy balance for urban surfaces and run 
!--          the material heat model
             IF (urban_surface) THEN
                CALL cpu_log( log_point(74), 'urban_surface', 'start' )
                CALL usm_surface_energy_balance
                IF ( usm_material_model )  THEN
                   CALL usm_material_heat_model
                ENDIF
                CALL cpu_log( log_point(74), 'urban_surface', 'stop' )
             ENDIF

!
!--          Compute the diffusion coefficients
             CALL cpu_log( log_point(17), 'diffusivities', 'start' )
             IF ( .NOT. humidity ) THEN
                IF ( ocean )  THEN
                   CALL diffusivities( prho, prho_reference )
                ELSE
                   CALL diffusivities( pt, pt_reference )
                ENDIF
             ELSE
                CALL diffusivities( vpt, pt_reference )
             ENDIF
             CALL cpu_log( log_point(17), 'diffusivities', 'stop' )

          ENDIF

!
!--       If required, calculate radiative fluxes and heating rates
          IF ( radiation .AND. intermediate_timestep_count                     &
               == intermediate_timestep_count_max .AND. simulated_time >    &
               skip_time_do_radiation )  THEN

               time_radiation = time_radiation + dt_3d

             IF ( time_radiation >= dt_radiation .OR. force_radiation_call )   &
             THEN

                CALL cpu_log( log_point(50), 'radiation', 'start' )

                IF ( .NOT. force_radiation_call )  THEN
                   time_radiation = time_radiation - dt_radiation
                ENDIF

                CALL radiation_control

                CALL cpu_log( log_point(50), 'radiation', 'stop' )

                IF (urban_surface)  THEN
                   CALL cpu_log( log_point(75), 'usm_radiation', 'start' )
                   CALL usm_radiation
                   CALL cpu_log( log_point(75), 'usm_radiation', 'stop' )
                ENDIF

             ENDIF
          ENDIF

       ENDDO   ! Intermediate step loop

!
!--    Increase simulation time and output times
       nr_timesteps_this_run      = nr_timesteps_this_run + 1
       current_timestep_number    = current_timestep_number + 1
       simulated_time             = simulated_time   + dt_3d
       simulated_time_chr         = time_to_string( simulated_time )
       time_since_reference_point = simulated_time - coupling_start_time



       IF ( simulated_time >= skip_time_data_output_av )  THEN
          time_do_av         = time_do_av       + dt_3d
       ENDIF
       IF ( simulated_time >= skip_time_do2d_xy )  THEN
          time_do2d_xy       = time_do2d_xy     + dt_3d
       ENDIF
       IF ( simulated_time >= skip_time_do2d_xz )  THEN
          time_do2d_xz       = time_do2d_xz     + dt_3d
       ENDIF
       IF ( simulated_time >= skip_time_do2d_yz )  THEN
          time_do2d_yz       = time_do2d_yz     + dt_3d
       ENDIF
       IF ( simulated_time >= skip_time_do3d    )  THEN
          time_do3d          = time_do3d        + dt_3d
       ENDIF
       DO  mid = 1, masks
          IF ( simulated_time >= skip_time_domask(mid) )  THEN
             time_domask(mid)= time_domask(mid) + dt_3d
          ENDIF
       ENDDO
       time_dvrp          = time_dvrp        + dt_3d
       IF ( simulated_time >= skip_time_dosp )  THEN
          time_dosp       = time_dosp        + dt_3d
       ENDIF
       time_dots          = time_dots        + dt_3d
       IF ( .NOT. first_call_lpm )  THEN
          time_dopts      = time_dopts       + dt_3d
       ENDIF
       IF ( simulated_time >= skip_time_dopr )  THEN
          time_dopr       = time_dopr        + dt_3d
       ENDIF
       time_dopr_listing          = time_dopr_listing        + dt_3d
       time_run_control   = time_run_control + dt_3d

!
!--    Data exchange between coupled models
       IF ( coupling_mode /= 'uncoupled'  .AND.  run_coupled )  THEN
          time_coupling = time_coupling + dt_3d

!
!--       In case of model termination initiated by the local model 
!--       (terminate_coupled > 0), the coupler must be skipped because it would 
!--       cause an MPI intercomminucation hang.
!--       If necessary, the coupler will be called at the beginning of the 
!--       next restart run.
          DO WHILE ( time_coupling >= dt_coupling .AND. terminate_coupled == 0 )
             CALL surface_coupler
             time_coupling = time_coupling - dt_coupling
          ENDDO
       ENDIF

!
!--    Execute user-defined actions
       CALL user_actions( 'after_integration' )

!
!--    If Galilei transformation is used, determine the distance that the
!--    model has moved so far
       IF ( galilei_transformation )  THEN
          advected_distance_x = advected_distance_x + u_gtrans * dt_3d
          advected_distance_y = advected_distance_y + v_gtrans * dt_3d
       ENDIF

!
!--    Check, if restart is necessary (because cpu-time is expiring or
!--    because it is forced by user) and set stop flag
!--    This call is skipped if the remote model has already initiated a restart.
       IF ( .NOT. terminate_run )  CALL check_for_restart

!
!--    Carry out statistical analysis and output at the requested output times.
!--    The MOD function is used for calculating the output time counters (like
!--    time_dopr) in order to regard a possible decrease of the output time
!--    interval in case of restart runs

!
!--    Set a flag indicating that so far no statistics have been created
!--    for this time step
       flow_statistics_called = .FALSE.

!
!--    If required, call flow_statistics for averaging in time
       IF ( averaging_interval_pr /= 0.0_wp  .AND.  &
            ( dt_dopr - time_dopr ) <= averaging_interval_pr  .AND.  &
            simulated_time >= skip_time_dopr )  THEN
          time_dopr_av = time_dopr_av + dt_3d
          IF ( time_dopr_av >= dt_averaging_input_pr )  THEN
             do_sum = .TRUE.
             time_dopr_av = MOD( time_dopr_av, &
                                    MAX( dt_averaging_input_pr, dt_3d ) )
          ENDIF
       ENDIF
       IF ( do_sum )  CALL flow_statistics

!
!--    Sum-up 3d-arrays for later output of time-averaged 2d/3d/masked data
       IF ( averaging_interval /= 0.0_wp  .AND.                                &
            ( dt_data_output_av - time_do_av ) <= averaging_interval  .AND. &
            simulated_time >= skip_time_data_output_av )                    &
       THEN
          time_do_sla = time_do_sla + dt_3d
          IF ( time_do_sla >= dt_averaging_input )  THEN
             CALL sum_up_3d_data
             average_count_3d = average_count_3d + 1
             time_do_sla = MOD( time_do_sla, MAX( dt_averaging_input, dt_3d ) )
          ENDIF
       ENDIF

!
!--    Calculate spectra for time averaging
       IF ( averaging_interval_sp /= 0.0_wp  .AND.  &
            ( dt_dosp - time_dosp ) <= averaging_interval_sp  .AND.  &
            simulated_time >= skip_time_dosp )  THEN
          time_dosp_av = time_dosp_av + dt_3d
          IF ( time_dosp_av >= dt_averaging_input_pr )  THEN
             CALL calc_spectra
             time_dosp_av = MOD( time_dosp_av, &
                                 MAX( dt_averaging_input_pr, dt_3d ) )
          ENDIF
       ENDIF

!
!--    Call flight module and output data
       IF ( virtual_flight )  THEN
          CALL flight_measurement
          CALL data_output_flight
       ENDIF

!
!--    Profile output (ASCII) on file
       IF ( time_dopr_listing >= dt_dopr_listing )  THEN
          CALL print_1d
          time_dopr_listing = MOD( time_dopr_listing, MAX( dt_dopr_listing, &
                                                           dt_3d ) )
       ENDIF

!
!--    Graphic output for PROFIL
       IF ( time_dopr >= dt_dopr )  THEN
          IF ( dopr_n /= 0 )  CALL data_output_profiles
          time_dopr = MOD( time_dopr, MAX( dt_dopr, dt_3d ) )
          time_dopr_av = 0.0_wp    ! due to averaging (see above)
       ENDIF

!
!--    Graphic output for time series
       IF ( time_dots >= dt_dots )  THEN
          CALL data_output_tseries
          time_dots = MOD( time_dots, MAX( dt_dots, dt_3d ) )
       ENDIF

!
!--    Output of spectra (formatted for use with PROFIL), in case of no
!--    time averaging, spectra has to be calculated before
       IF ( time_dosp >= dt_dosp )  THEN
          IF ( average_count_sp == 0 )  CALL calc_spectra
          CALL data_output_spectra
          time_dosp = MOD( time_dosp, MAX( dt_dosp, dt_3d ) )
       ENDIF

!
!--    2d-data output (cross-sections)
       IF ( time_do2d_xy >= dt_do2d_xy )  THEN
          CALL data_output_2d( 'xy', 0 )
          time_do2d_xy = MOD( time_do2d_xy, MAX( dt_do2d_xy, dt_3d ) )
       ENDIF
       IF ( time_do2d_xz >= dt_do2d_xz )  THEN
          CALL data_output_2d( 'xz', 0 )
          time_do2d_xz = MOD( time_do2d_xz, MAX( dt_do2d_xz, dt_3d ) )
       ENDIF
       IF ( time_do2d_yz >= dt_do2d_yz )  THEN
          CALL data_output_2d( 'yz', 0 )
          time_do2d_yz = MOD( time_do2d_yz, MAX( dt_do2d_yz, dt_3d ) )
       ENDIF

!
!--    3d-data output (volume data)
       IF ( time_do3d >= dt_do3d )  THEN
          CALL data_output_3d( 0 )
          time_do3d = MOD( time_do3d, MAX( dt_do3d, dt_3d ) )
       ENDIF

!
!--    Masked data output
       DO  mid = 1, masks
          IF ( time_domask(mid) >= dt_domask(mid) )  THEN
             CALL data_output_mask( 0 )
             time_domask(mid) = MOD( time_domask(mid),  &
                                     MAX( dt_domask(mid), dt_3d ) )
          ENDIF
       ENDDO

!
!--    Output of time-averaged 2d/3d/masked data
       IF ( time_do_av >= dt_data_output_av )  THEN
          CALL average_3d_data
          CALL data_output_2d( 'xy', 1 )
          CALL data_output_2d( 'xz', 1 )
          CALL data_output_2d( 'yz', 1 )
          CALL data_output_3d( 1 )
          DO  mid = 1, masks
             CALL data_output_mask( 1 )
          ENDDO
          time_do_av = MOD( time_do_av, MAX( dt_data_output_av, dt_3d ) )
       ENDIF

!
!--    Output of particle time series
       IF ( particle_advection )  THEN
          IF ( time_dopts >= dt_dopts  .OR. &
               ( simulated_time >= particle_advection_start  .AND. &
                 first_call_lpm ) )  THEN
             CALL data_output_ptseries
             time_dopts = MOD( time_dopts, MAX( dt_dopts, dt_3d ) )
          ENDIF
       ENDIF

!
!--    Output of dvrp-graphics (isosurface, particles, slicer)
#if defined( __dvrp_graphics )
       CALL DVRP_LOG_EVENT( -2, current_timestep_number-1 )
#endif
       IF ( time_dvrp >= dt_dvrp )  THEN
          CALL data_output_dvrp
          time_dvrp = MOD( time_dvrp, MAX( dt_dvrp, dt_3d ) )
       ENDIF
#if defined( __dvrp_graphics )
       CALL DVRP_LOG_EVENT( 2, current_timestep_number )
#endif

!
!--    If required, set the heat flux for the next time step at a random value
       IF ( constant_heatflux  .AND.  random_heatflux )  CALL disturb_heatflux

!
!--    Execute user-defined actions
       CALL user_actions( 'after_timestep' )

!
!--    Determine size of next time step. Save timestep dt_3d because it is
!--    newly calculated in routine timestep, but required further below for
!--    steering the run control output interval
       dt_3d_old = dt_3d
       CALL timestep

!
!--    Synchronize the timestep in case of nested run.
       IF ( nested_run )  THEN
!
!--       Synchronize by unifying the time step.
!--       Global minimum of all time-steps is used for all.
          CALL pmci_synchronize
       ENDIF

!
!--    Computation and output of run control parameters.
!--    This is also done whenever perturbations have been imposed
       IF ( time_run_control >= dt_run_control  .OR.                     &
            timestep_scheme(1:5) /= 'runge'  .OR.  disturbance_created ) &
       THEN
          CALL run_control
          IF ( time_run_control >= dt_run_control )  THEN
             time_run_control = MOD( time_run_control, &
                                     MAX( dt_run_control, dt_3d_old ) )
          ENDIF
       ENDIF

!
!--    Output elapsed simulated time in form of a progress bar on stdout
       IF ( myid == 0 )  CALL output_progress_bar

       CALL cpu_log( log_point_s(10), 'timesteps', 'stop' )


    ENDDO   ! time loop

    IF ( myid == 0 )  CALL finish_progress_bar

#if defined( __dvrp_graphics )
    CALL DVRP_LOG_EVENT( -2, current_timestep_number )
#endif

    CALL location_message( 'finished time-stepping', .TRUE. )

 END SUBROUTINE time_integration
