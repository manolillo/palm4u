!> @file prognostic_equations.f90
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
! kk: Added chemistry calculations
! bK: Added missing initialization of s_init (tendency_caller)
! FKa: Bugfix concerning exchange of ghost points after chemical reactions
!      (separate loop now)
! FKa: nbgp added to ONLY list, minor formatting, added some cpu log_points
! 
! Former revisions:
! -----------------
! $Id: prognostic_equations.f90 2159 2017-02-22 18:01:07Z kanani $
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
! 
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added pt tendency calculation based on energy balance at urban surfaces
! (new urban surface model)
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1976 2016-07-27 13:28:04Z maronga
! Simplied calls to radiation model
! 
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! 
! 1914 2016-05-26 14:44:07Z witha
! Added calls for wind turbine model
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 1826 2016-04-07 12:01:39Z maronga
! Renamed canopy model calls.
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! Kessler microphysics scheme moved to microphysics.
!
! 1757 2016-02-22 15:49:32Z maronga
! 
! 1691 2015-10-26 16:17:44Z maronga
! Added optional model spin-up without radiation / land surface model calls.
! Formatting corrections.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1585 2015-04-30 07:05:52Z maronga
! Added call for temperature tendency calculation due to radiative flux divergence
! 
! 1517 2015-01-07 19:12:25Z hoffmann
! advec_s_bc_mod addded, since advec_s_bc is now a module
!
! 1496 2014-12-02 17:25:50Z maronga
! Renamed "radiation" -> "cloud_top_radiation"
! 
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   parameters cthf and plant_canopy moved to module plant_canopy_model_mod.
! Removed double-listing of use_upstream_for_tke in ONLY-list of module
! control_parameters
! 
! 1409 2014-05-23 12:11:32Z suehring
! Bugfix: i_omp_start changed for advec_u_ws at left inflow and outflow boundary. 
! This ensures that left-hand side fluxes are also calculated for nxl in that
! case, even though the solution at nxl is overwritten in boundary_conds() 
! 
! 1398 2014-05-07 11:15:00Z heinze
! Rayleigh-damping for horizontal velocity components changed: instead of damping
! against ug and vg, damping against u_init and v_init is used to allow for a 
! homogenized treatment in case of nudging
! 
! 1380 2014-04-28 12:40:45Z heinze
! Change order of calls for scalar prognostic quantities: 
! ls_advec -> nudging -> subsidence since initial profiles 
! 
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY lists
! 
! 1365 2014-04-22 15:03:56Z boeske
! Calls of ls_advec for large scale advection added, 
! subroutine subsidence is only called if use_subsidence_tendencies = .F.,
! new argument ls_index added to the calls of subsidence
! +ls_index
! 
! 1361 2014-04-16 15:17:48Z hoffmann
! Two-moment microphysics moved to the start of prognostic equations. This makes
! the 3d arrays for tend_q, tend_qr, tend_pt and tend_pt redundant.
! Additionally, it is allowed to call the microphysics just once during the time
! step (not at each sub-time step).
!
! Two-moment cloud physics added for vector and accelerator optimization.
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
! 
! 1337 2014-03-25 15:11:48Z heinze
! Bugfix: REAL constants provided with KIND-attribute
! 
! 1332 2014-03-25 11:59:43Z suehring
! Bugfix: call advec_ws or advec_pw for TKE only if NOT use_upstream_for_tke 
! 
! 1330 2014-03-24 17:29:32Z suehring
! In case of SGS-particle velocity advection of TKE is also allowed with 
! dissipative 5th-order scheme.
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop vector clauses removed, independent clauses added
!
! 1246 2013-11-01 08:59:45Z heinze
! enable nudging also for accelerator version
!
! 1241 2013-10-30 11:36:58Z heinze
! usage of nudging enabled (so far not implemented for accelerator version)
!
! 1179 2013-06-14 05:57:58Z raasch
! two arguments removed from routine buoyancy, ref_state updated on device
!
! 1128 2013-04-12 06:19:32Z raasch
! those parts requiring global communication moved to time_integration,
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1115 2013-03-26 18:16:16Z hoffmann
! optimized cloud physics: calculation of microphysical tendencies transfered 
! to microphysics.f90; qr and nr are only calculated if precipitation is required
!
! 1111 2013-03-08 23:54:10Z raasch
! update directives for prognostic quantities removed
!
! 1106 2013-03-04 05:31:38Z raasch
! small changes in code formatting
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1053 2012-11-13 17:11:03Z hoffmann
! implementation of two new prognostic equations for rain drop concentration (nr)
! and rain water content (qr)
!
! currently, only available for cache loop optimization
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1019 2012-09-28 06:46:45Z raasch
! non-optimized version of prognostic_equations removed
!
! 1015 2012-09-27 09:23:24Z raasch
! new branch prognostic_equations_acc
! OpenACC statements added + code changes required for GPU optimization
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog- and upstream-spline-scheme removed
!
! 978 2012-08-09 08:28:32Z fricke
! km_damp_x and km_damp_y removed in calls of diffusion_u and diffusion_v
! add ptdf_x, ptdf_y for damping the potential temperature at the inflow
! boundary in case of non-cyclic lateral boundaries
! Bugfix: first thread index changes for WS-scheme at the inflow
!
! 940 2012-07-09 14:31:00Z raasch
! temperature equation can be switched off
!
! Revision 1.1  2000/04/13 14:56:27  schroeter
! Initial revision
!
!
! Description:
! ------------
!> Solving the prognostic equations.
!>
!> @todo tendency_caller needs to be a subroutine in kchem_driver
!>       (SUBROUTINE kchem_tendency)
!> @todo CALL tendency_caller-->kchem_tendency has to be adjusted in case of
!>       other advection scheme is used. With the provided arguments right now
!>       it works only for ws-scheme
!------------------------------------------------------------------------------!
 MODULE prognostic_equations_mod

 

    USE arrays_3d,                                                             &
        ONLY:  diss_l_e, diss_l_nr, diss_l_pt, diss_l_q, diss_l_qr,            &
               diss_l_s, diss_l_sa, diss_s_e, diss_s_nr, diss_s_pt, diss_s_q,  &
               diss_s_qr, diss_s_s, diss_s_sa, e, e_p, flux_s_e, flux_s_nr,    &
               flux_s_pt, flux_s_q, flux_s_qr, flux_s_s, flux_s_sa, flux_l_e,  &
               flux_l_nr, flux_l_pt, flux_l_q, flux_l_qr, flux_l_s, flux_l_sa, &
               nr, nr_p, nrsws, nrswst, pt, ptdf_x, ptdf_y, pt_init, pt_p,     &
               prho, q, q_init, q_p, qsws, qswst, qr, qr_p, qrsws, qrswst, rdf,&
               rdf_sc, ref_state, rho_ocean, s, s_init, s_p, sa, sa_init, sa_p,      &
               saswsb, saswst, shf, ssws, sswst, tend,                         &
               te_m, tnr_m, tpt_m, tq_m, tqr_m, ts_m, tsa_m, tswst, tu_m, tv_m,&
               tw_m, u, ug, u_init, u_p, v, vg, vpt, v_init, v_p, w, w_p
        
    USE control_parameters,                                                    &
        ONLY:  call_microphysics_at_all_substeps, cloud_physics,               &
               cloud_top_radiation, constant_diffusion, dp_external,           &
               dp_level_ind_b, dp_smooth_factor, dpdxy, dt_3d, humidity,       &
               inflow_l, intermediate_timestep_count,                          &
               intermediate_timestep_count_max, large_scale_forcing,           &
               large_scale_subsidence, microphysics_seifert,                   &
               microphysics_sat_adjust, neutral, nudging, ocean, outflow_l,    &
               outflow_s, passive_scalar, prho_reference, prho_reference,      &
               prho_reference, pt_reference, pt_reference, pt_reference,       &
               scalar_advec, scalar_advec, simulated_time, sloping_surface,    &
               timestep_scheme, tsc, urban_surface, use_subsidence_tendencies, &
               use_upstream_for_tke, wall_heatflux,                            &
               wall_nrflux, wall_qflux, wall_qflux, wall_qflux, wall_qrflux,   &
               wall_salinityflux, wall_sflux, ws_scheme_mom, ws_scheme_sca

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE eqn_state_seawater_mod,                                                &
        ONLY:  eqn_state_seawater

    USE indices,                                                               &
        ONLY:  i_left, i_right, j_north, j_south, nxl, nxlu, nxr, nyn, nys,    &
               nysv, nzb_s_inner, nzb_u_inner, nzb_v_inner, nzb_w_inner, nzb, nzt, &
               nbgp

    USE advec_ws,                                                              &
        ONLY:  advec_s_ws, advec_s_ws_acc, advec_u_ws, advec_u_ws_acc,         &
               advec_v_ws, advec_v_ws_acc, advec_w_ws, advec_w_ws_acc

    USE advec_s_bc_mod,                                                        &
        ONLY:  advec_s_bc

    USE advec_s_pw_mod,                                                        &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                        &
        ONLY:  advec_s_up

    USE advec_u_pw_mod,                                                        &
        ONLY:  advec_u_pw

    USE advec_u_up_mod,                                                        &
        ONLY:  advec_u_up

    USE advec_v_pw_mod,                                                        &
        ONLY:  advec_v_pw

    USE advec_v_up_mod,                                                        &
        ONLY:  advec_v_up

    USE advec_w_pw_mod,                                                        &
        ONLY:  advec_w_pw

    USE advec_w_up_mod,                                                        &
        ONLY:  advec_w_up

    USE buoyancy_mod,                                                          &
        ONLY:  buoyancy, buoyancy_acc

    USE calc_radiation_mod,                                                    &
        ONLY:  calc_radiation
  
    USE coriolis_mod,                                                          &
        ONLY:  coriolis, coriolis_acc

    USE diffusion_e_mod,                                                       &
        ONLY:  diffusion_e, diffusion_e_acc

    USE diffusion_s_mod,                                                       &
        ONLY:  diffusion_s, diffusion_s_acc

    USE diffusion_u_mod,                                                       &
        ONLY:  diffusion_u, diffusion_u_acc

    USE diffusion_v_mod,                                                       &
        ONLY:  diffusion_v, diffusion_v_acc

    USE diffusion_w_mod,                                                       &
        ONLY:  diffusion_w, diffusion_w_acc

    USE kinds

    USE ls_forcing_mod,                                                        &
        ONLY:  ls_advec

    USE microphysics_mod,                                                      &
        ONLY:  microphysics_control

    USE nudge_mod,                                                             &
        ONLY:  nudge

    USE plant_canopy_model_mod,                                                &
        ONLY:  cthf, plant_canopy, pcm_tendency

    USE production_e_mod,                                                      &
        ONLY:  production_e, production_e_acc

    USE radiation_model_mod,                                                   &
        ONLY:  radiation, radiation_tendency,                                  &
               skip_time_do_radiation

    USE statistics,                                                            &
        ONLY:  hom

    USE subsidence_mod,                                                        &
        ONLY:  subsidence

    USE urban_surface_mod,                                                     &
        ONLY:  usm_wall_heat_flux

    USE user_actions_mod,                                                      &
        ONLY:  user_actions

    USE wind_turbine_model_mod,                                                &
        ONLY:  wind_turbine, wtm_tendencies

#ifdef KPP_CHEM
    USE kchem_driver,                                                          &
        ONLY: chem_species, kchem_integrate, NSPEC, use_kpp_chemistry
#endif

    PRIVATE
    PUBLIC prognostic_equations_cache, prognostic_equations_vector, &
           prognostic_equations_acc

    INTERFACE prognostic_equations_cache
       MODULE PROCEDURE prognostic_equations_cache
    END INTERFACE prognostic_equations_cache

    INTERFACE prognostic_equations_vector
       MODULE PROCEDURE prognostic_equations_vector
    END INTERFACE prognostic_equations_vector

    INTERFACE prognostic_equations_acc
       MODULE PROCEDURE prognostic_equations_acc
    END INTERFACE prognostic_equations_acc


 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Version with one optimized loop over all equations. It is only allowed to
!> be called for the Wicker and Skamarock or Piascek-Williams advection scheme.
!>
!> Here the calls of most subroutines are embedded in two DO loops over i and j,
!> so communication between CPUs is not allowed (does not make sense) within
!> these loops.
!>
!> (Optimized to avoid cache missings, i.e. for Power4/5-architectures.)
!------------------------------------------------------------------------------!
 
 SUBROUTINE prognostic_equations_cache


    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  i_omp_start         !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  omp_get_thread_num  !<
    INTEGER(iwp) ::  tn = 0              !<
    
    LOGICAL      ::  loop_start          !<
    INTEGER      ::  n


!
!-- Time measurement can only be performed for the whole set of equations
    CALL cpu_log( log_point(32), 'all progn.equations', 'start' )

!
!-- Calculation of chemical reactions. This is done outside of main loop,
!-- since exchange of ghost points is required after this update of the
!-- concentrations of chemical species                                    
#ifdef KPP_CHEM
    IF ( use_kpp_chemistry )  THEN
       CALL cpu_log( log_point(32), 'all progn.equations', 'pause' )
       CALL cpu_log( log_point_s(82), 'chemistry reactions ', 'start' )
       
       DO  i = nxl, nxr
          DO  j = nys, nyn
             CALL kchem_integrate (i,j)                                                
          ENDDO
       ENDDO

       CALL cpu_log( log_point_s(84), 'chemistry exch-horiz ', 'start' )
!--    Loop over chemical species       
       DO  n = 1, NSPEC
          CALL exchange_horiz( chem_species(n)%conc, nbgp )                        
       ENDDO
       CALL cpu_log( log_point_s(84), 'chemistry exch-horiz ', 'stop' )
       
       CALL cpu_log( log_point_s(82), 'chemistry reactions ', 'stop' )
       CALL cpu_log( log_point(32), 'all progn.equations', 'continue' )
    ENDIF       
#endif     

!-- Loop over all prognostic equations
!$OMP PARALLEL private (i,i_omp_start,j,k,loop_start,tn)

!$  tn = omp_get_thread_num() 
    loop_start = .TRUE.
!$OMP DO
    DO  i = nxl, nxr

!
!--    Store the first loop index. It differs for each thread and is required
!--    later in advec_ws
       IF ( loop_start )  THEN
          loop_start  = .FALSE.
          i_omp_start = i 
       ENDIF

       DO  j = nys, nyn
!
!--       If required, calculate cloud microphysics
          IF ( cloud_physics  .AND.  .NOT. microphysics_sat_adjust  .AND.      &
               ( intermediate_timestep_count == 1  .OR.                        &
                 call_microphysics_at_all_substeps )                           &
             )  THEN
             CALL microphysics_control( i, j )
          ENDIF
!
!--       Tendency terms for u-velocity component
          IF ( .NOT. outflow_l  .OR.  i > nxl )  THEN

             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_mom )  THEN
                   CALL advec_u_ws( i, j, i_omp_start, tn )
                ELSE 
                   CALL advec_u_pw( i, j )
                ENDIF 
             ELSE
                CALL advec_u_up( i, j )
             ENDIF
             CALL diffusion_u( i, j )
             CALL coriolis( i, j, 1 )
             IF ( sloping_surface  .AND.  .NOT. neutral )  THEN
                CALL buoyancy( i, j, pt, 1 )
             ENDIF

!
!--          Drag by plant canopy
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 1 )

!
!--          External pressure gradient
             IF ( dp_external )  THEN
                DO  k = dp_level_ind_b+1, nzt
                   tend(k,j,i) = tend(k,j,i) - dpdxy(1) * dp_smooth_factor(k)
                ENDDO
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'u' )

!
!--          Forces by wind turbines
             IF ( wind_turbine )  CALL wtm_tendencies( i, j, 1 )

             CALL user_actions( i, j, 'u-tendency' )
!
!--          Prognostic equation for u-velocity component
             DO  k = nzb_u_inner(j,i)+1, nzt
                u_p(k,j,i) = u(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * tu_m(k,j,i) )       &
                                      - tsc(5) * rdf(k) * ( u(k,j,i) - u_init(k) )
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb_u_inner(j,i)+1, nzt
                      tu_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb_u_inner(j,i)+1, nzt
                      tu_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tu_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF

!
!--       Tendency terms for v-velocity component
          IF ( .NOT. outflow_s  .OR.  j > nys )  THEN

             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_mom )  THEN
                    CALL advec_v_ws( i, j, i_omp_start, tn )
                ELSE 
                    CALL advec_v_pw( i, j )
                ENDIF
             ELSE
                CALL advec_v_up( i, j )
             ENDIF
             CALL diffusion_v( i, j )
             CALL coriolis( i, j, 2 )

!
!--          Drag by plant canopy
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 2 )        

!
!--          External pressure gradient
             IF ( dp_external )  THEN
                DO  k = dp_level_ind_b+1, nzt
                   tend(k,j,i) = tend(k,j,i) - dpdxy(2) * dp_smooth_factor(k)
                ENDDO
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'v' )

!
!--          Forces by wind turbines
             IF ( wind_turbine )  CALL wtm_tendencies( i, j, 2 )

             CALL user_actions( i, j, 'v-tendency' )
!
!--          Prognostic equation for v-velocity component
             DO  k = nzb_v_inner(j,i)+1, nzt
                v_p(k,j,i) = v(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * tv_m(k,j,i) )       &
                                      - tsc(5) * rdf(k) * ( v(k,j,i) - v_init(k) )
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb_v_inner(j,i)+1, nzt
                      tv_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb_v_inner(j,i)+1, nzt
                      tv_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tv_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF

!
!--       Tendency terms for w-velocity component
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_mom )  THEN
                CALL advec_w_ws( i, j, i_omp_start, tn )
             ELSE 
                CALL advec_w_pw( i, j )
             END IF
          ELSE
             CALL advec_w_up( i, j )
          ENDIF
          CALL diffusion_w( i, j )
          CALL coriolis( i, j, 3 )

          IF ( .NOT. neutral )  THEN
             IF ( ocean )  THEN
                CALL buoyancy( i, j, rho_ocean, 3 )
             ELSE
                IF ( .NOT. humidity )  THEN
                   CALL buoyancy( i, j, pt, 3 )
                ELSE
                   CALL buoyancy( i, j, vpt, 3 )
                ENDIF
             ENDIF
          ENDIF

!
!--       Drag by plant canopy
          IF ( plant_canopy )  CALL pcm_tendency( i, j, 3 )

!
!--       Forces by wind turbines
          IF ( wind_turbine )  CALL wtm_tendencies( i, j, 3 )

          CALL user_actions( i, j, 'w-tendency' )

!
!--       Prognostic equation for w-velocity component
          DO  k = nzb_w_inner(j,i)+1, nzt-1
             w_p(k,j,i) = w(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +          &
                                               tsc(3) * tw_m(k,j,i) )          &
                                   - tsc(5) * rdf(k) * w(k,j,i)
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb_w_inner(j,i)+1, nzt-1
                   tw_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb_w_inner(j,i)+1, nzt-1
                   tw_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tw_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

!
!--       If required, compute prognostic equation for potential temperature
          IF ( .NOT. neutral )  THEN
!
!--          Tendency terms for potential temperature
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( ws_scheme_sca )  THEN
                      CALL advec_s_ws( i, j, pt, 'pt', flux_s_pt, diss_s_pt, &
                                       flux_l_pt, diss_l_pt, i_omp_start, tn )
                   ELSE
                      CALL advec_s_pw( i, j, pt )
                   ENDIF
             ELSE
                CALL advec_s_up( i, j, pt )
             ENDIF
             CALL diffusion_s( i, j, pt, shf, tswst, wall_heatflux )

!
!--          Tendency pt from wall heat flux from urban surface
             IF ( urban_surface )  THEN
                CALL usm_wall_heat_flux( i, j )
             ENDIF

!
!--          If required compute heating/cooling due to long wave radiation
!--          processes
             IF ( cloud_top_radiation )  THEN
                CALL calc_radiation( i, j )
             ENDIF

!
!--          Consideration of heat sources within the plant canopy
             IF ( plant_canopy  .AND.  cthf /= 0.0_wp )  THEN
                CALL pcm_tendency( i, j, 4 )
             ENDIF

!
!--          Large scale advection
             IF ( large_scale_forcing )  THEN
                CALL ls_advec( i, j, simulated_time, 'pt' )
             ENDIF     

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'pt' ) 

!
!--          If required, compute effect of large-scale subsidence/ascent
             IF ( large_scale_subsidence  .AND.                                &
                  .NOT. use_subsidence_tendencies )  THEN
                CALL subsidence( i, j, tend, pt, pt_init, 2 )
             ENDIF

!
!--          If required, add tendency due to radiative heating/cooling
             IF ( radiation  .AND.                                             &
                  simulated_time > skip_time_do_radiation )  THEN
                CALL radiation_tendency ( i, j, tend )
             ENDIF


             CALL user_actions( i, j, 'pt-tendency' )

!
!--          Prognostic equation for potential temperature
             DO  k = nzb_s_inner(j,i)+1, nzt
                pt_p(k,j,i) = pt(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +     &
                                                    tsc(3) * tpt_m(k,j,i) )    &
                                        - tsc(5) * ( pt(k,j,i) - pt_init(k) ) *&
                                          ( rdf_sc(k) + ptdf_x(i) + ptdf_y(j) )
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tpt_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tpt_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                                      5.3125_wp * tpt_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF

!
!--       If required, compute prognostic equation for salinity
          IF ( ocean )  THEN

!
!--          Tendency-terms for salinity
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' ) &
             THEN
                IF ( ws_scheme_sca )  THEN
                    CALL advec_s_ws( i, j, sa, 'sa', flux_s_sa,  &
                                diss_s_sa, flux_l_sa, diss_l_sa, i_omp_start, tn  )
                ELSE 
                    CALL advec_s_pw( i, j, sa )
                ENDIF
             ELSE
                CALL advec_s_up( i, j, sa )
             ENDIF
             CALL diffusion_s( i, j, sa, saswsb, saswst, wall_salinityflux )

             CALL user_actions( i, j, 'sa-tendency' )

!
!--          Prognostic equation for salinity
             DO  k = nzb_s_inner(j,i)+1, nzt
                sa_p(k,j,i) = sa(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +     &
                                                    tsc(3) * tsa_m(k,j,i) )    &
                                        - tsc(5) * rdf_sc(k) *                 &
                                          ( sa(k,j,i) - sa_init(k) )
                IF ( sa_p(k,j,i) < 0.0_wp )  sa_p(k,j,i) = 0.1_wp * sa(k,j,i)
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tsa_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tsa_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                                      5.3125_wp * tsa_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

!
!--          Calculate density by the equation of state for seawater
             CALL eqn_state_seawater( i, j )

          ENDIF

!
!--       If required, compute prognostic equation for total water content
          IF ( humidity )  THEN

!
!--          Tendency-terms for total water content / scalar
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' ) &
             THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( i, j, q, 'q', flux_s_q, & 
                                diss_s_q, flux_l_q, diss_l_q, i_omp_start, tn )
                ELSE 
                   CALL advec_s_pw( i, j, q )
                ENDIF
             ELSE
                CALL advec_s_up( i, j, q )
             ENDIF
             CALL diffusion_s( i, j, q, qsws, qswst, wall_qflux )

!
!--          Sink or source of humidity due to canopy elements
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 5 )

!
!--          Large scale advection
             IF ( large_scale_forcing )  THEN
                CALL ls_advec( i, j, simulated_time, 'q' )
             ENDIF

!
!--          Nudging
             IF ( nudging )  CALL nudge( i, j, simulated_time, 'q' ) 

!
!--          If required compute influence of large-scale subsidence/ascent
             IF ( large_scale_subsidence  .AND.                                &
                  .NOT. use_subsidence_tendencies )  THEN
                CALL subsidence( i, j, tend, q, q_init, 3 )
             ENDIF

             CALL user_actions( i, j, 'q-tendency' )

!
!--          Prognostic equation for total water content / scalar
             DO  k = nzb_s_inner(j,i)+1, nzt
                q_p(k,j,i) = q(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * tq_m(k,j,i) )       &
                                      - tsc(5) * rdf_sc(k) *                   &
                                        ( q(k,j,i) - q_init(k) )
                IF ( q_p(k,j,i) < 0.0_wp )  q_p(k,j,i) = 0.1_wp * q(k,j,i)
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tq_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tq_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                                     5.3125_wp * tq_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

!
!--          If required, calculate prognostic equations for rain water content 
!--          and rain drop concentration
             IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
!
!--             Calculate prognostic equation for rain water content
                tend(:,j,i) = 0.0_wp
                IF ( timestep_scheme(1:5) == 'runge' ) &
                THEN
                   IF ( ws_scheme_sca )  THEN
                      CALL advec_s_ws( i, j, qr, 'qr', flux_s_qr,       & 
                                       diss_s_qr, flux_l_qr, diss_l_qr, &
                                       i_omp_start, tn )
                   ELSE 
                      CALL advec_s_pw( i, j, qr )
                   ENDIF
                ELSE
                   CALL advec_s_up( i, j, qr )
                ENDIF
                CALL diffusion_s( i, j, qr, qrsws, qrswst, wall_qrflux )

!
!--             Prognostic equation for rain water content
                DO  k = nzb_s_inner(j,i)+1, nzt
                   qr_p(k,j,i) = qr(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +  &
                                                       tsc(3) * tqr_m(k,j,i) ) &
                                           - tsc(5) * rdf_sc(k) * qr(k,j,i)
                   IF ( qr_p(k,j,i) < 0.0_wp )  qr_p(k,j,i) = 0.0_wp
                ENDDO
!
!--             Calculate tendencies for the next Runge-Kutta step
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( intermediate_timestep_count == 1 )  THEN
                      DO  k = nzb_s_inner(j,i)+1, nzt
                         tqr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ELSEIF ( intermediate_timestep_count < &
                            intermediate_timestep_count_max )  THEN
                      DO  k = nzb_s_inner(j,i)+1, nzt
                        tqr_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                                        5.3125_wp * tqr_m(k,j,i)
                      ENDDO
                   ENDIF
                ENDIF

!
!--             Calculate prognostic equation for rain drop concentration.
                tend(:,j,i) = 0.0_wp
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( ws_scheme_sca )  THEN
                      CALL advec_s_ws( i, j, nr, 'nr', flux_s_nr,    & 
                                    diss_s_nr, flux_l_nr, diss_l_nr, &
                                    i_omp_start, tn )
                   ELSE 
                      CALL advec_s_pw( i, j, nr )
                   ENDIF
                ELSE
                   CALL advec_s_up( i, j, nr )
                ENDIF
                CALL diffusion_s( i, j, nr, nrsws, nrswst, wall_nrflux )

!
!--             Prognostic equation for rain drop concentration
                DO  k = nzb_s_inner(j,i)+1, nzt
                   nr_p(k,j,i) = nr(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +  &
                                                       tsc(3) * tnr_m(k,j,i) ) &
                                           - tsc(5) * rdf_sc(k) * nr(k,j,i)
                   IF ( nr_p(k,j,i) < 0.0_wp )  nr_p(k,j,i) = 0.0_wp
                ENDDO
!
!--             Calculate tendencies for the next Runge-Kutta step
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( intermediate_timestep_count == 1 )  THEN
                      DO  k = nzb_s_inner(j,i)+1, nzt
                         tnr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ELSEIF ( intermediate_timestep_count < &
                            intermediate_timestep_count_max )  THEN
                      DO  k = nzb_s_inner(j,i)+1, nzt
                         tnr_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                                        5.3125_wp * tnr_m(k,j,i)
                      ENDDO
                   ENDIF
                ENDIF

             ENDIF

          ENDIF
          
!
!--       If required, compute prognostic equation for scalar
          IF ( passive_scalar )  THEN
!
!--          Tendency-terms for total water content / scalar
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' ) &
             THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( i, j, s, 's', flux_s_s, & 
                                diss_s_s, flux_l_s, diss_l_s, i_omp_start, tn )
                ELSE 
                   CALL advec_s_pw( i, j, s )
                ENDIF
             ELSE
                CALL advec_s_up( i, j, s )
             ENDIF
             CALL diffusion_s( i, j, s, ssws, sswst, wall_sflux )

!
!--          Sink or source of scalar concentration due to canopy elements
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 7 )

!
!--          Large scale advection, still need to be extended for scalars
!              IF ( large_scale_forcing )  THEN
!                 CALL ls_advec( i, j, simulated_time, 's' )
!              ENDIF

!
!--          Nudging, still need to be extended for scalars
!              IF ( nudging )  CALL nudge( i, j, simulated_time, 's' ) 

!
!--          If required compute influence of large-scale subsidence/ascent.
!--          Note, the last argument is of no meaning in this case, as it is 
!--          only used in conjunction with large_scale_forcing, which is to 
!--          date not implemented for scalars.
             IF ( large_scale_subsidence  .AND.                                &
                  .NOT. use_subsidence_tendencies  .AND.                       &
                  .NOT. large_scale_forcing )  THEN
                CALL subsidence( i, j, tend, s, s_init, 3 )
             ENDIF

             CALL user_actions( i, j, 's-tendency' )

!
!--          Prognostic equation for scalar
             DO  k = nzb_s_inner(j,i)+1, nzt
                s_p(k,j,i) = s(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * ts_m(k,j,i) )       &
                                      - tsc(5) * rdf_sc(k) *                   &
                                        ( s(k,j,i) - s_init(k) )
                IF ( s_p(k,j,i) < 0.0_wp )  s_p(k,j,i) = 0.1_wp * s(k,j,i)
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      ts_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      ts_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                                     5.3125_wp * ts_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF          
!
!--       If required, compute prognostic equation for turbulent kinetic 
!--       energy (TKE)
          IF ( .NOT. constant_diffusion )  THEN

!
!--          Tendency-terms for TKE
             tend(:,j,i) = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge'  &
                 .AND.  .NOT. use_upstream_for_tke )  THEN 
                 IF ( ws_scheme_sca )  THEN
                     CALL advec_s_ws( i, j, e, 'e', flux_s_e, diss_s_e, &
                                      flux_l_e, diss_l_e , i_omp_start, tn )
                 ELSE
                     CALL advec_s_pw( i, j, e )
                 ENDIF
             ELSE
                CALL advec_s_up( i, j, e )
             ENDIF
             IF ( .NOT. humidity )  THEN
                IF ( ocean )  THEN
                   CALL diffusion_e( i, j, prho, prho_reference )
                ELSE
                   CALL diffusion_e( i, j, pt, pt_reference )
                ENDIF
             ELSE
                CALL diffusion_e( i, j, vpt, pt_reference )
             ENDIF
             CALL production_e( i, j )

!
!--          Additional sink term for flows through plant canopies
             IF ( plant_canopy )  CALL pcm_tendency( i, j, 6 ) 

             CALL user_actions( i, j, 'e-tendency' )

!
!--          Prognostic equation for TKE.
!--          Eliminate negative TKE values, which can occur due to numerical
!--          reasons in the course of the integration. In such cases the old
!--          TKE value is reduced by 90%.
             DO  k = nzb_s_inner(j,i)+1, nzt
                e_p(k,j,i) = e(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * te_m(k,j,i) )
                IF ( e_p(k,j,i) < 0.0_wp )  e_p(k,j,i) = 0.1_wp * e(k,j,i)
             ENDDO

!
!--          Calculate tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      te_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ELSEIF ( intermediate_timestep_count < &
                         intermediate_timestep_count_max )  THEN
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      te_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                                     5.3125_wp * te_m(k,j,i)
                   ENDDO
                ENDIF
             ENDIF

          ENDIF   ! TKE equation

#ifdef KPP_CHEM
          IF ( use_kpp_chemistry )  THEN
             CALL cpu_log( log_point(32), 'all progn.equations', 'pause' )  
             CALL cpu_log( log_point_s(83), 'chemistry advec+diff+prog ', 'start' )
             
!--          Loop over chemical species
             DO  n = 1, NSPEC  
                CALL tendency_caller ( chem_species(n)%conc_p, chem_species(n)%conc, &
                                       chem_species(n)%tconc_m,                      &
                                       i, j, i_omp_start, tn,                        &
                                       chem_species(n)%ssws, chem_species(n)%sswst,  &
                                       chem_species(n)%flux_s, chem_species(n)%diss_s, &
                                       chem_species(n)%flux_l, chem_species(n)%diss_l)       
             ENDDO

             CALL cpu_log( log_point_s(83), 'chemistry advec+diff+prog ', 'stop' )
             CALL cpu_log( log_point(32), 'all progn.equations', 'continue' )
          ENDIF                
#endif

       ENDDO
    ENDDO
!$OMP END PARALLEL


    CALL cpu_log( log_point(32), 'all progn.equations', 'stop' )


 END SUBROUTINE prognostic_equations_cache


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Version for vector machines
!------------------------------------------------------------------------------!
 
 SUBROUTINE prognostic_equations_vector


    IMPLICIT NONE

    INTEGER(iwp) ::  i    !<
    INTEGER(iwp) ::  j    !<
    INTEGER(iwp) ::  k    !<

    REAL(wp)     ::  sbt  !<


!
!-- If required, calculate cloud microphysical impacts
    IF ( cloud_physics  .AND.  .NOT. microphysics_sat_adjust  .AND.            &
         ( intermediate_timestep_count == 1  .OR.                              &
           call_microphysics_at_all_substeps )                                 &
       )  THEN
       CALL cpu_log( log_point(51), 'microphysics', 'start' )
       CALL microphysics_control
       CALL cpu_log( log_point(51), 'microphysics', 'stop' )
    ENDIF

!
!-- u-velocity component
    CALL cpu_log( log_point(5), 'u-equation', 'start' )

    tend = 0.0_wp
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_u_ws
       ELSE
          CALL advec_u_pw
       ENDIF
    ELSE
       CALL advec_u_up
    ENDIF
    CALL diffusion_u
    CALL coriolis( 1 )
    IF ( sloping_surface  .AND.  .NOT. neutral )  THEN
       CALL buoyancy( pt, 1 )
    ENDIF

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 1 )

!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = nxlu, nxr
          DO  j = nys, nyn
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(1) * dp_smooth_factor(k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Nudging
    IF ( nudging )  CALL nudge( simulated_time, 'u' )

!
!-- Forces by wind turbines
    IF ( wind_turbine )  CALL wtm_tendencies( 1 )

    CALL user_actions( 'u-tendency' )

!
!-- Prognostic equation for u-velocity component
    DO  i = nxlu, nxr
       DO  j = nys, nyn
          DO  k = nzb_u_inner(j,i)+1, nzt
             u_p(k,j,i) = u(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +          &
                                               tsc(3) * tu_m(k,j,i) )          &
                                   - tsc(5) * rdf(k) * ( u(k,j,i) - u_init(k) )
          ENDDO
       ENDDO
    ENDDO

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  i = nxlu, nxr
             DO  j = nys, nyn
                DO  k = nzb_u_inner(j,i)+1, nzt
                   tu_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          DO  i = nxlu, nxr
             DO  j = nys, nyn
                DO  k = nzb_u_inner(j,i)+1, nzt
                   tu_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tu_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(5), 'u-equation', 'stop' )

!
!-- v-velocity component
    CALL cpu_log( log_point(6), 'v-equation', 'start' )

    tend = 0.0_wp
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_v_ws
       ELSE 
          CALL advec_v_pw
       END IF
    ELSE
       CALL advec_v_up
    ENDIF
    CALL diffusion_v
    CALL coriolis( 2 )

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 2 )

!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = nxl, nxr
          DO  j = nysv, nyn
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(2) * dp_smooth_factor(k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Nudging
    IF ( nudging )  CALL nudge( simulated_time, 'v' )

!
!-- Forces by wind turbines
    IF ( wind_turbine )  CALL wtm_tendencies( 2 )

    CALL user_actions( 'v-tendency' )

!
!-- Prognostic equation for v-velocity component
    DO  i = nxl, nxr
       DO  j = nysv, nyn
          DO  k = nzb_v_inner(j,i)+1, nzt
             v_p(k,j,i) = v(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +          &
                                               tsc(3) * tv_m(k,j,i) )          &
                                   - tsc(5) * rdf(k) * ( v(k,j,i) - v_init(k) )
          ENDDO
       ENDDO
    ENDDO

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  i = nxl, nxr
             DO  j = nysv, nyn
                DO  k = nzb_v_inner(j,i)+1, nzt
                   tv_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          DO  i = nxl, nxr
             DO  j = nysv, nyn
                DO  k = nzb_v_inner(j,i)+1, nzt
                   tv_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tv_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(6), 'v-equation', 'stop' )

!
!-- w-velocity component
    CALL cpu_log( log_point(7), 'w-equation', 'start' )

    tend = 0.0_wp
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_w_ws
       ELSE
          CALL advec_w_pw
       ENDIF
    ELSE
       CALL advec_w_up
    ENDIF
    CALL diffusion_w
    CALL coriolis( 3 )

    IF ( .NOT. neutral )  THEN
       IF ( ocean )  THEN
          CALL buoyancy( rho_ocean, 3 )
       ELSE
          IF ( .NOT. humidity )  THEN
             CALL buoyancy( pt, 3 )
          ELSE
             CALL buoyancy( vpt, 3 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 3 )

!
!-- Forces by wind turbines
    IF ( wind_turbine )  CALL wtm_tendencies( 3 )

    CALL user_actions( 'w-tendency' )

!
!-- Prognostic equation for w-velocity component
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb_w_inner(j,i)+1, nzt-1
             w_p(k,j,i) = w(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +          &
                                               tsc(3) * tw_m(k,j,i) )          &
                                   - tsc(5) * rdf(k) * w(k,j,i)
          ENDDO
       ENDDO
    ENDDO

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_w_inner(j,i)+1, nzt-1
                   tw_m(k,j,i) = tend(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_w_inner(j,i)+1, nzt-1
                   tw_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tw_m(k,j,i)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(7), 'w-equation', 'stop' )


!
!-- If required, compute prognostic equation for potential temperature
    IF ( .NOT. neutral )  THEN

       CALL cpu_log( log_point(13), 'pt-equation', 'start' )

!
!--    pt-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( pt, 'pt' )

       ENDIF

!
!--    pt-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( pt, 'pt' )
             ELSE
                CALL advec_s_pw( pt )
             ENDIF
          ELSE
             CALL advec_s_up( pt )
          ENDIF
       ENDIF

       CALL diffusion_s( pt, shf, tswst, wall_heatflux )

!
!--    Tendency pt from wall heat flux from urban surface
       IF ( urban_surface )  THEN
          CALL usm_wall_heat_flux
       ENDIF

!
!--    If required compute heating/cooling due to long wave radiation processes
       IF ( cloud_top_radiation )  THEN
          CALL calc_radiation
       ENDIF

!
!--    Consideration of heat sources within the plant canopy
       IF ( plant_canopy .AND. ( cthf /= 0.0_wp ) ) THEN
          CALL pcm_tendency( 4 )
       ENDIF

!
!--    Large scale advection
       IF ( large_scale_forcing )  THEN
          CALL ls_advec( simulated_time, 'pt' )
       ENDIF

!
!--    Nudging
       IF ( nudging )  CALL nudge( simulated_time, 'pt' ) 

!
!--    If required compute influence of large-scale subsidence/ascent
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies )  THEN
          CALL subsidence( tend, pt, pt_init, 2 )
       ENDIF

!
!--    If required, add tendency due to radiative heating/cooling
       IF ( radiation  .AND.                                                   &
            simulated_time > skip_time_do_radiation )  THEN
            CALL radiation_tendency ( tend )
       ENDIF

       CALL user_actions( 'pt-tendency' )

!
!--    Prognostic equation for potential temperature
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_s_inner(j,i)+1, nzt
                pt_p(k,j,i) = pt(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +        &
                                                    tsc(3) * tpt_m(k,j,i) )    &
                                        - tsc(5) * ( pt(k,j,i) - pt_init(k) ) *&
                                          ( rdf_sc(k) + ptdf_x(i) + ptdf_y(j) )
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tpt_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tpt_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                                      5.3125_wp * tpt_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(13), 'pt-equation', 'stop' )

    ENDIF

!
!-- If required, compute prognostic equation for salinity
    IF ( ocean )  THEN

       CALL cpu_log( log_point(37), 'sa-equation', 'start' )

!
!--    sa-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( sa, 'sa' )

       ENDIF

!
!--    sa-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                 CALL advec_s_ws( sa, 'sa' )
             ELSE
                 CALL advec_s_pw( sa )
             ENDIF
          ELSE
             CALL advec_s_up( sa )
          ENDIF
       ENDIF

       CALL diffusion_s( sa, saswsb, saswst, wall_salinityflux )
       
       CALL user_actions( 'sa-tendency' )

!
!--    Prognostic equation for salinity
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_s_inner(j,i)+1, nzt
                sa_p(k,j,i) = sa(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +        &
                                                    tsc(3) * tsa_m(k,j,i) )    &
                                        - tsc(5) * rdf_sc(k) *                 &
                                          ( sa(k,j,i) - sa_init(k) )
                IF ( sa_p(k,j,i) < 0.0_wp )  sa_p(k,j,i) = 0.1_wp * sa(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tsa_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tsa_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                                      5.3125_wp * tsa_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(37), 'sa-equation', 'stop' )

!
!--    Calculate density by the equation of state for seawater
       CALL cpu_log( log_point(38), 'eqns-seawater', 'start' )
       CALL eqn_state_seawater
       CALL cpu_log( log_point(38), 'eqns-seawater', 'stop' )

    ENDIF

!
!-- If required, compute prognostic equation for total water content
    IF ( humidity )  THEN

       CALL cpu_log( log_point(29), 'q-equation', 'start' )

!
!--    Scalar/q-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( q, 'q' )

       ENDIF

!
!--    Scalar/q-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( q, 'q' )
             ELSE
                CALL advec_s_pw( q )
             ENDIF
          ELSE
             CALL advec_s_up( q )
          ENDIF
       ENDIF

       CALL diffusion_s( q, qsws, qswst, wall_qflux )
       
!
!--    Sink or source of humidity due to canopy elements
       IF ( plant_canopy ) CALL pcm_tendency( 5 )

!
!--    Large scale advection
       IF ( large_scale_forcing )  THEN
          CALL ls_advec( simulated_time, 'q' )
       ENDIF

!
!--    Nudging
       IF ( nudging )  CALL nudge( simulated_time, 'q' ) 

!
!--    If required compute influence of large-scale subsidence/ascent
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies )  THEN
         CALL subsidence( tend, q, q_init, 3 )
       ENDIF

       CALL user_actions( 'q-tendency' )

!
!--    Prognostic equation for total water content
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_s_inner(j,i)+1, nzt
                q_p(k,j,i) = q(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +          &
                                                  tsc(3) * tq_m(k,j,i) )       &
                                      - tsc(5) * rdf_sc(k) *                   &
                                        ( q(k,j,i) - q_init(k) )
                IF ( q_p(k,j,i) < 0.0_wp )  q_p(k,j,i) = 0.1_wp * q(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tq_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      tq_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tq_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(29), 'q-equation', 'stop' )

!
!--    If required, calculate prognostic equations for rain water content 
!--    and rain drop concentration
       IF ( cloud_physics  .AND.  microphysics_seifert )  THEN

          CALL cpu_log( log_point(52), 'qr-equation', 'start' )

!
!--       Calculate prognostic equation for rain water content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qr, 'qr' )

          ENDIF

!
!--       qr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( qr, 'qr' )
                ELSE
                   CALL advec_s_pw( qr )
                ENDIF
             ELSE
                CALL advec_s_up( qr )
             ENDIF
          ENDIF

          CALL diffusion_s( qr, qrsws, qrswst, wall_qrflux )

          CALL user_actions( 'qr-tendency' )

!
!--       Prognostic equation for rain water content
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_s_inner(j,i)+1, nzt
                   qr_p(k,j,i) = qr(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +     &
                                                     tsc(3) * tqr_m(k,j,i) )   &
                                           - tsc(5) * rdf_sc(k) * qr(k,j,i)
                   IF ( qr_p(k,j,i) < 0.0_wp )  qr_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_s_inner(j,i)+1, nzt
                         tqr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_s_inner(j,i)+1, nzt
                         tqr_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * &
                                                                   tqr_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(52), 'qr-equation', 'stop' )
          CALL cpu_log( log_point(53), 'nr-equation', 'start' )

!
!--       Calculate prognostic equation for rain drop concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( nr, 'nr' )

          ENDIF

!
!--       nr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( nr, 'nr' )
                ELSE
                   CALL advec_s_pw( nr )
                ENDIF
             ELSE
                CALL advec_s_up( nr )
             ENDIF
          ENDIF

          CALL diffusion_s( nr, nrsws, nrswst, wall_nrflux )

!
!--       Prognostic equation for rain drop concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_s_inner(j,i)+1, nzt
                   nr_p(k,j,i) = nr(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +     &
                                                     tsc(3) * tnr_m(k,j,i) )   &
                                           - tsc(5) * rdf_sc(k) * nr(k,j,i)
                   IF ( nr_p(k,j,i) < 0.0_wp )  nr_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_s_inner(j,i)+1, nzt
                         tnr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_s_inner(j,i)+1, nzt
                         tnr_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * &
                                                                   tnr_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

          CALL cpu_log( log_point(53), 'nr-equation', 'stop' )

       ENDIF

    ENDIF
!
!-- If required, compute prognostic equation for scalar
    IF ( passive_scalar )  THEN

       CALL cpu_log( log_point(66), 's-equation', 'start' )

!
!--    Scalar/q-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( s, 's' )

       ENDIF

!
!--    Scalar/q-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( s, 's' )
             ELSE
                CALL advec_s_pw( s )
             ENDIF
          ELSE
             CALL advec_s_up( s )
          ENDIF
       ENDIF

       CALL diffusion_s( s, ssws, sswst, wall_sflux )
       
!
!--    Sink or source of humidity due to canopy elements
       IF ( plant_canopy ) CALL pcm_tendency( 7 )

!
!--    Large scale advection. Not implemented for scalars so far.
!        IF ( large_scale_forcing )  THEN
!           CALL ls_advec( simulated_time, 'q' )
!        ENDIF

!
!--    Nudging. Not implemented for scalars so far.
!        IF ( nudging )  CALL nudge( simulated_time, 'q' ) 

!
!--    If required compute influence of large-scale subsidence/ascent.
!--    Not implemented for scalars so far.
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies  .AND.                             &
            .NOT. large_scale_forcing )  THEN
         CALL subsidence( tend, s, s_init, 3 )
       ENDIF

       CALL user_actions( 's-tendency' )

!
!--    Prognostic equation for total water content
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_s_inner(j,i)+1, nzt
                s_p(k,j,i) = s(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +          &
                                                  tsc(3) * ts_m(k,j,i) )       &
                                      - tsc(5) * rdf_sc(k) *                   &
                                        ( s(k,j,i) - s_init(k) )
                IF ( s_p(k,j,i) < 0.0_wp )  s_p(k,j,i) = 0.1_wp * s(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      ts_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      ts_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * ts_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(66), 's-equation', 'stop' )

    ENDIF
!
!-- If required, compute prognostic equation for turbulent kinetic 
!-- energy (TKE)
    IF ( .NOT. constant_diffusion )  THEN

       CALL cpu_log( log_point(16), 'tke-equation', 'start' )

       sbt = tsc(2)
       IF ( .NOT. use_upstream_for_tke )  THEN
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( e, 'e' )

          ENDIF
       ENDIF

!
!--    TKE-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme'  .OR.  use_upstream_for_tke )  THEN
          IF ( use_upstream_for_tke )  THEN
             tend = 0.0_wp
             CALL advec_s_up( e )
          ELSE
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( e, 'e' )
                ELSE 
                   CALL advec_s_pw( e )
                ENDIF
             ELSE
                CALL advec_s_up( e )
             ENDIF
          ENDIF
       ENDIF

       IF ( .NOT. humidity )  THEN
          IF ( ocean )  THEN
             CALL diffusion_e( prho, prho_reference )
          ELSE
             CALL diffusion_e( pt, pt_reference )
          ENDIF
       ELSE
          CALL diffusion_e( vpt, pt_reference )
       ENDIF

       CALL production_e

!
!--    Additional sink term for flows through plant canopies
       IF ( plant_canopy )  CALL pcm_tendency( 6 )
       CALL user_actions( 'e-tendency' )

!
!--    Prognostic equation for TKE.
!--    Eliminate negative TKE values, which can occur due to numerical
!--    reasons in the course of the integration. In such cases the old TKE
!--    value is reduced by 90%.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_s_inner(j,i)+1, nzt
                e_p(k,j,i) = e(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +          &
                                                  tsc(3) * te_m(k,j,i) )
                IF ( e_p(k,j,i) < 0.0_wp )  e_p(k,j,i) = 0.1_wp * e(k,j,i)
             ENDDO
          ENDDO
       ENDDO

!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      te_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      te_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * te_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       CALL cpu_log( log_point(16), 'tke-equation', 'stop' )

    ENDIF

 END SUBROUTINE prognostic_equations_vector


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Version for accelerator boards
!------------------------------------------------------------------------------!
 
 SUBROUTINE prognostic_equations_acc


    IMPLICIT NONE

    INTEGER(iwp) ::  i           !<
    INTEGER(iwp) ::  j           !<
    INTEGER(iwp) ::  k           !<
    INTEGER(iwp) ::  runge_step  !<

    REAL(wp)     ::  sbt         !<

!
!-- Set switch for intermediate Runge-Kutta step
    runge_step = 0
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          runge_step = 1
       ELSEIF ( intermediate_timestep_count < &
                intermediate_timestep_count_max )  THEN
          runge_step = 2
       ENDIF
    ENDIF

!
!-- If required, calculate cloud microphysical impacts (two-moment scheme)
    IF ( cloud_physics  .AND.  .NOT. microphysics_sat_adjust  .AND.            &
         ( intermediate_timestep_count == 1  .OR.                              &
           call_microphysics_at_all_substeps )                                 &
       )  THEN
       CALL cpu_log( log_point(51), 'microphysics', 'start' )
       CALL microphysics_control
       CALL cpu_log( log_point(51), 'microphysics', 'stop' )
    ENDIF

!
!-- u-velocity component
!++ Statistics still not completely ported to accelerators
    !$acc update device( hom, ref_state )
    CALL cpu_log( log_point(5), 'u-equation', 'start' )

    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_u_ws_acc
       ELSE
          tend = 0.0_wp   ! to be removed later??
          CALL advec_u_pw
       ENDIF
    ELSE
       CALL advec_u_up
    ENDIF
    CALL diffusion_u_acc
    CALL coriolis_acc( 1 )
    IF ( sloping_surface  .AND.  .NOT. neutral )  THEN
       CALL buoyancy( pt, 1 )
    ENDIF

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 1 )

!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(1) * dp_smooth_factor(k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Nudging
    IF ( nudging )  CALL nudge( simulated_time, 'u' )

!
!-- Forces by wind turbines
    IF ( wind_turbine )  CALL wtm_tendencies( 1 )

    CALL user_actions( 'u-tendency' )

!
!-- Prognostic equation for u-velocity component
    !$acc kernels present( nzb_u_inner, rdf, tend, tu_m, u, u_init, u_p )
    !$acc loop independent
    DO  i = i_left, i_right
       !$acc loop independent
       DO  j = j_south, j_north
          !$acc loop independent
          DO  k = 1, nzt
             IF ( k > nzb_u_inner(j,i) )  THEN
                u_p(k,j,i) = u(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * tu_m(k,j,i) )       &
                                      - tsc(5) * rdf(k) * ( u(k,j,i) - u_init(k) )
!
!--             Tendencies for the next Runge-Kutta step
                IF ( runge_step == 1 )  THEN
                   tu_m(k,j,i) = tend(k,j,i)
                ELSEIF ( runge_step == 2 )  THEN
                   tu_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tu_m(k,j,i)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !$acc end kernels

    CALL cpu_log( log_point(5), 'u-equation', 'stop' )

!
!-- v-velocity component
    CALL cpu_log( log_point(6), 'v-equation', 'start' )

    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_v_ws_acc
       ELSE
          tend = 0.0_wp    ! to be removed later??
          CALL advec_v_pw
       END IF
    ELSE
       CALL advec_v_up
    ENDIF
    CALL diffusion_v_acc
    CALL coriolis_acc( 2 )

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 2 )

!
!-- External pressure gradient
    IF ( dp_external )  THEN
       DO  i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = dp_level_ind_b+1, nzt
                tend(k,j,i) = tend(k,j,i) - dpdxy(2) * dp_smooth_factor(k)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

!
!-- Nudging
    IF ( nudging )  CALL nudge( simulated_time, 'v' )

!
!-- Forces by wind turbines
    IF ( wind_turbine )  CALL wtm_tendencies( 2 )

    CALL user_actions( 'v-tendency' )

!
!-- Prognostic equation for v-velocity component
    !$acc kernels present( nzb_v_inner, rdf, tend, tv_m, v, v_init, v_p )
    !$acc loop independent
    DO  i = i_left, i_right
       !$acc loop independent
       DO  j = j_south, j_north
          !$acc loop independent
          DO  k = 1, nzt
             IF ( k > nzb_v_inner(j,i) )  THEN
                v_p(k,j,i) = v(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * tv_m(k,j,i) )       &
                                      - tsc(5) * rdf(k) * ( v(k,j,i) - v_init(k) )
!
!--             Tendencies for the next Runge-Kutta step
                IF ( runge_step == 1 )  THEN
                   tv_m(k,j,i) = tend(k,j,i)
                ELSEIF ( runge_step == 2 )  THEN
                   tv_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tv_m(k,j,i)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !$acc end kernels

    CALL cpu_log( log_point(6), 'v-equation', 'stop' )

!
!-- w-velocity component
    CALL cpu_log( log_point(7), 'w-equation', 'start' )

    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_mom )  THEN
          CALL advec_w_ws_acc
       ELSE
          tend = 0.0_wp    ! to be removed later??
          CALL advec_w_pw
       ENDIF
    ELSE
       CALL advec_w_up
    ENDIF
    CALL diffusion_w_acc
    CALL coriolis_acc( 3 )

    IF ( .NOT. neutral )  THEN
       IF ( ocean )  THEN
          CALL buoyancy( rho_ocean, 3 )
       ELSE
          IF ( .NOT. humidity )  THEN
             CALL buoyancy_acc( pt, 3 )
          ELSE
             CALL buoyancy( vpt, 3 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Drag by plant canopy
    IF ( plant_canopy )  CALL pcm_tendency( 3 )

!
!-- Forces by wind turbines
    IF ( wind_turbine )  CALL wtm_tendencies( 3 )

    CALL user_actions( 'w-tendency' )

!
!-- Prognostic equation for w-velocity component
    !$acc kernels present( nzb_w_inner, rdf, tend, tw_m, w, w_p )
    !$acc loop independent
    DO  i = i_left, i_right
       !$acc loop independent
       DO  j = j_south, j_north
          !$acc loop independent
          DO  k = 1, nzt-1
             IF ( k > nzb_w_inner(j,i) )  THEN
                w_p(k,j,i) = w(k,j,i) + dt_3d * ( tsc(2) * tend(k,j,i) +       &
                                                  tsc(3) * tw_m(k,j,i) )       &
                                      - tsc(5) * rdf(k) * w(k,j,i)
   !
   !--          Tendencies for the next Runge-Kutta step
                IF ( runge_step == 1 )  THEN
                   tw_m(k,j,i) = tend(k,j,i)
                ELSEIF ( runge_step == 2 )  THEN
                   tw_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tw_m(k,j,i)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !$acc end kernels

    CALL cpu_log( log_point(7), 'w-equation', 'stop' )


!
!-- If required, compute prognostic equation for potential temperature
    IF ( .NOT. neutral )  THEN

       CALL cpu_log( log_point(13), 'pt-equation', 'start' )

!
!--    pt-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( pt, 'pt' )

       ENDIF

!
!--    pt-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws_acc( pt, 'pt' )
             ELSE
                tend = 0.0_wp    ! to be removed later??
                CALL advec_s_pw( pt )
             ENDIF
          ELSE
             CALL advec_s_up( pt )
          ENDIF
       ENDIF

       CALL diffusion_s_acc( pt, shf, tswst, wall_heatflux )

!
!--    Tendency pt from wall heat flux from urban surface
       IF ( urban_surface )  THEN
          CALL usm_wall_heat_flux
       ENDIF

!
!--    If required compute heating/cooling due to long wave radiation processes
       IF ( cloud_top_radiation )  THEN
          CALL calc_radiation
       ENDIF

!
!--    Consideration of heat sources within the plant canopy
       IF ( plant_canopy .AND. ( cthf /= 0.0_wp ) ) THEN
          CALL pcm_tendency( 4 )
       ENDIF

!
!--    Large scale advection
       IF ( large_scale_forcing )  THEN
          CALL ls_advec( simulated_time, 'pt' )
       ENDIF

!
!--    Nudging
       IF ( nudging )  CALL nudge( simulated_time, 'pt' ) 

!
!--    If required compute influence of large-scale subsidence/ascent
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies )  THEN
          CALL subsidence( tend, pt, pt_init, 2 )
       ENDIF

       IF ( radiation .AND.                                                    &
            simulated_time > skip_time_do_radiation )  THEN
            CALL radiation_tendency ( tend )
       ENDIF

       CALL user_actions( 'pt-tendency' )

!
!--    Prognostic equation for potential temperature
       !$acc kernels present( nzb_s_inner, rdf_sc, ptdf_x, ptdf_y, pt_init ) &
       !$acc         present( tend, tpt_m, pt, pt_p )
       !$acc loop independent
       DO  i = i_left, i_right
          !$acc loop independent
          DO  j = j_south, j_north
             !$acc loop independent
             DO  k = 1, nzt
                IF ( k > nzb_s_inner(j,i) )  THEN
                   pt_p(k,j,i) = pt(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +        &
                                                       tsc(3) * tpt_m(k,j,i) )    &
                                           - tsc(5) * ( pt(k,j,i) - pt_init(k) ) *&
                                             ( rdf_sc(k) + ptdf_x(i) + ptdf_y(j) )
!
!--                Tendencies for the next Runge-Kutta step
                   IF ( runge_step == 1 )  THEN
                      tpt_m(k,j,i) = tend(k,j,i)
                   ELSEIF ( runge_step == 2 )  THEN
                      tpt_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tpt_m(k,j,i)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       !$acc end kernels

       CALL cpu_log( log_point(13), 'pt-equation', 'stop' )

    ENDIF

!
!-- If required, compute prognostic equation for salinity
    IF ( ocean )  THEN

       CALL cpu_log( log_point(37), 'sa-equation', 'start' )

!
!--    sa-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( sa, 'sa' )

       ENDIF

!
!--    sa-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                 CALL advec_s_ws( sa, 'sa' )
             ELSE
                 CALL advec_s_pw( sa )
             ENDIF
          ELSE
             CALL advec_s_up( sa )
          ENDIF
       ENDIF

       CALL diffusion_s( sa, saswsb, saswst, wall_salinityflux )

       CALL user_actions( 'sa-tendency' )

!
!--    Prognostic equation for salinity
       DO  i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = nzb_s_inner(j,i)+1, nzt
                sa_p(k,j,i) = sa(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +        &
                                                    tsc(3) * tsa_m(k,j,i) )    &
                                        - tsc(5) * rdf_sc(k) *                 &
                                          ( sa(k,j,i) - sa_init(k) )
                IF ( sa_p(k,j,i) < 0.0_wp )  sa_p(k,j,i) = 0.1_wp * sa(k,j,i)
!
!--             Tendencies for the next Runge-Kutta step
                IF ( runge_step == 1 )  THEN
                   tsa_m(k,j,i) = tend(k,j,i)
                ELSEIF ( runge_step == 2 )  THEN
                   tsa_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tsa_m(k,j,i)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point(37), 'sa-equation', 'stop' )

!
!--    Calculate density by the equation of state for seawater
       CALL cpu_log( log_point(38), 'eqns-seawater', 'start' )
       CALL eqn_state_seawater
       CALL cpu_log( log_point(38), 'eqns-seawater', 'stop' )

    ENDIF

!
!-- If required, compute prognostic equation for total water content
    IF ( humidity )  THEN

       CALL cpu_log( log_point(29), 'q-equation', 'start' )

!
!--    Scalar/q-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( q, 'q' )

       ENDIF

!
!--    Scalar/q-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( q, 'q' )
             ELSE
                CALL advec_s_pw( q )
             ENDIF
          ELSE
             CALL advec_s_up( q )
          ENDIF
       ENDIF

       CALL diffusion_s( q, qsws, qswst, wall_qflux )

!
!--    Sink or source of scalar concentration due to canopy elements
       IF ( plant_canopy ) CALL pcm_tendency( 5 )

!
!--    Large scale advection
       IF ( large_scale_forcing )  THEN
          CALL ls_advec( simulated_time, 'q' )
       ENDIF

!
!--    Nudging
       IF ( nudging )  CALL nudge( simulated_time, 'q' ) 

!
!--    If required compute influence of large-scale subsidence/ascent
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies )  THEN
         CALL subsidence( tend, q, q_init, 3 )
       ENDIF

       CALL user_actions( 'q-tendency' )

!
!--    Prognostic equation for total water content / scalar
       DO  i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = nzb_s_inner(j,i)+1, nzt
                q_p(k,j,i) = q(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +          &
                                                  tsc(3) * tq_m(k,j,i) )       &
                                      - tsc(5) * rdf_sc(k) *                   &
                                        ( q(k,j,i) - q_init(k) )
                IF ( q_p(k,j,i) < 0.0_wp )  q_p(k,j,i) = 0.1_wp * q(k,j,i)
!
!--             Tendencies for the next Runge-Kutta step
                IF ( runge_step == 1 )  THEN
                   tq_m(k,j,i) = tend(k,j,i)
                ELSEIF ( runge_step == 2 )  THEN
                   tq_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * tq_m(k,j,i)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point(29), 'q-equation', 'stop' )

!
!--    If required, calculate prognostic equations for rain water content 
!--    and rain drop concentration
       IF ( cloud_physics  .AND.  microphysics_seifert )  THEN

          CALL cpu_log( log_point(52), 'qr-equation', 'start' )
!
!--       qr-tendency terms with communication
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qr, 'qr' )

          ENDIF

!
!--       qr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( qr, 'qr' )
                ELSE
                   CALL advec_s_pw( qr )
                ENDIF
             ELSE
                CALL advec_s_up( qr )
             ENDIF
          ENDIF

          CALL diffusion_s( qr, qrsws, qrswst, wall_qrflux )

!
!--       Prognostic equation for rain water content
          DO  i = i_left, i_right
             DO  j = j_south, j_north
                DO  k = nzb_s_inner(j,i)+1, nzt
                   qr_p(k,j,i) = qr(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +     &
                                                       tsc(3) * tqr_m(k,j,i) ) &
                                           - tsc(5) * rdf_sc(k) * qr(k,j,i)
                   IF ( qr_p(k,j,i) < 0.0_wp )  qr_p(k,j,i) = 0.0_wp
!
!--                Tendencies for the next Runge-Kutta step
                   IF ( runge_step == 1 )  THEN
                      tqr_m(k,j,i) = tend(k,j,i)
                   ELSEIF ( runge_step == 2 )  THEN
                      tqr_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp *    &
                                                                tqr_m(k,j,i)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          CALL cpu_log( log_point(52), 'qr-equation', 'stop' )
          CALL cpu_log( log_point(53), 'nr-equation', 'start' )

!
!--       nr-tendency terms with communication
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( nr, 'nr' )

          ENDIF

!
!--       nr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws( nr, 'nr' )
                ELSE
                   CALL advec_s_pw( nr )
                ENDIF
             ELSE
                CALL advec_s_up( nr )
             ENDIF
          ENDIF

          CALL diffusion_s( nr, nrsws, nrswst, wall_nrflux )

!
!--       Prognostic equation for rain drop concentration
          DO  i = i_left, i_right
             DO  j = j_south, j_north
                DO  k = nzb_s_inner(j,i)+1, nzt
                   nr_p(k,j,i) = nr(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +     &
                                                       tsc(3) * tnr_m(k,j,i) ) &
                                           - tsc(5) * rdf_sc(k) * nr(k,j,i)
                   IF ( nr_p(k,j,i) < 0.0_wp )  nr_p(k,j,i) = 0.0_wp
!
!--                Tendencies for the next Runge-Kutta step
                   IF ( runge_step == 1 )  THEN
                      tnr_m(k,j,i) = tend(k,j,i)
                   ELSEIF ( runge_step == 2 )  THEN
                      tnr_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp *    &
                                                                tnr_m(k,j,i)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          CALL cpu_log( log_point(53), 'nr-equation', 'stop' )

       ENDIF

    ENDIF

!
!-- If required, compute prognostic equation for scalar
    IF ( passive_scalar )  THEN

       CALL cpu_log( log_point(66), 's-equation', 'start' )

!
!--    Scalar/q-tendency terms with communication
       sbt = tsc(2)
       IF ( scalar_advec == 'bc-scheme' )  THEN

          IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--          Bott-Chlond scheme always uses Euler time step. Thus:
             sbt = 1.0_wp
          ENDIF
          tend = 0.0_wp
          CALL advec_s_bc( s, 's' )

       ENDIF

!
!--    Scalar/q-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme' )  THEN
          tend = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( s, 's' )
             ELSE
                CALL advec_s_pw( s )
             ENDIF
          ELSE
             CALL advec_s_up( s )
          ENDIF
       ENDIF

       CALL diffusion_s( s, ssws, sswst, wall_sflux )

!
!--    Sink or source of scalar concentration due to canopy elements
       IF ( plant_canopy ) CALL pcm_tendency( 7 )

!
!--    Large scale advection. Not implemented so far.
!        IF ( large_scale_forcing )  THEN
!           CALL ls_advec( simulated_time, 's' )
!        ENDIF

!
!--    Nudging. Not implemented so far.
!        IF ( nudging )  CALL nudge( simulated_time, 's' ) 

!
!--    If required compute influence of large-scale subsidence/ascent.
!--    Not implemented so far.
       IF ( large_scale_subsidence  .AND.                                      &
            .NOT. use_subsidence_tendencies  .AND.                             &
            .NOT. large_scale_forcing )  THEN
         CALL subsidence( tend, s, s_init, 3 )
       ENDIF

       CALL user_actions( 's-tendency' )

!
!--    Prognostic equation for total water content / scalar
       DO  i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = nzb_s_inner(j,i)+1, nzt
                s_p(k,j,i) = s(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +          &
                                                  tsc(3) * ts_m(k,j,i) )       &
                                      - tsc(5) * rdf_sc(k) *                   &
                                        ( s(k,j,i) - s_init(k) )
                IF ( s_p(k,j,i) < 0.0_wp )  s_p(k,j,i) = 0.1_wp * s(k,j,i)
!
!--             Tendencies for the next Runge-Kutta step
                IF ( runge_step == 1 )  THEN
                   ts_m(k,j,i) = tend(k,j,i)
                ELSEIF ( runge_step == 2 )  THEN
                   ts_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * ts_m(k,j,i)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       CALL cpu_log( log_point(66), 's-equation', 'stop' )

    ENDIF
!
!-- If required, compute prognostic equation for turbulent kinetic
!-- energy (TKE)
    IF ( .NOT. constant_diffusion )  THEN

       CALL cpu_log( log_point(16), 'tke-equation', 'start' )

       sbt = tsc(2)
       IF ( .NOT. use_upstream_for_tke )  THEN
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( e, 'e' )

          ENDIF
       ENDIF

!
!--    TKE-tendency terms with no communication
       IF ( scalar_advec /= 'bc-scheme'  .OR.  use_upstream_for_tke )  THEN
          IF ( use_upstream_for_tke )  THEN
             tend = 0.0_wp
             CALL advec_s_up( e )
          ELSE
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   CALL advec_s_ws_acc( e, 'e' )
                ELSE
                   tend = 0.0_wp    ! to be removed later??
                   CALL advec_s_pw( e )
                ENDIF
             ELSE
                tend = 0.0_wp    ! to be removed later??
                CALL advec_s_up( e )
             ENDIF
          ENDIF
       ENDIF
!!
       IF ( .NOT. humidity )  THEN
          IF ( ocean )  THEN
             CALL diffusion_e( prho, prho_reference )
          ELSE
             CALL diffusion_e_acc( pt, pt_reference )
          ENDIF
       ELSE
          CALL diffusion_e( vpt, pt_reference )
       ENDIF

       CALL production_e_acc

!
!--    Additional sink term for flows through plant canopies
       IF ( plant_canopy )  CALL pcm_tendency( 6 )
       CALL user_actions( 'e-tendency' )

!
!--    Prognostic equation for TKE.
!--    Eliminate negative TKE values, which can occur due to numerical
!--    reasons in the course of the integration. In such cases the old TKE
!--    value is reduced by 90%.
       !$acc kernels present( e, e_p, nzb_s_inner, tend, te_m )
       !$acc loop independent
       DO  i = i_left, i_right
          !$acc loop independent
          DO  j = j_south, j_north
             !$acc loop independent
             DO  k = 1, nzt
                IF ( k > nzb_s_inner(j,i) )  THEN
                   e_p(k,j,i) = e(k,j,i) + dt_3d * ( sbt * tend(k,j,i) +          &
                                                     tsc(3) * te_m(k,j,i) )
                   IF ( e_p(k,j,i) < 0.0_wp )  e_p(k,j,i) = 0.1_wp * e(k,j,i)
!
!--                Tendencies for the next Runge-Kutta step
                   IF ( runge_step == 1 )  THEN
                      te_m(k,j,i) = tend(k,j,i)
                   ELSEIF ( runge_step == 2 )  THEN
                      te_m(k,j,i) = -9.5625_wp * tend(k,j,i) + 5.3125_wp * te_m(k,j,i)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       !$acc end kernels

       CALL cpu_log( log_point(16), 'tke-equation', 'stop' )

    ENDIF

 END SUBROUTINE prognostic_equations_acc

 ! Module private subroutines

 SUBROUTINE tendency_caller ( rs_p, rs, trs_m, i, j, i_omp_start, tn, ssws_kc, sswst_kc, flux_s, diss_s, flux_l, diss_l )
    USE pegrid,   ONLY: threads_per_task
    USE indices,  ONLY: nysg,nyng,nxlg,nxrg
    IMPLICIT NONE

    REAL(wp), DIMENSION(:,:,:), POINTER                         :: rs_p, rs, trs_m

    INTEGER(iwp),INTENT(IN) :: i, j, i_omp_start, tn
    REAL(wp),INTENT(IN),DIMENSION(nysg:nyng,nxlg:nxrg)          :: ssws_kc  !<
    REAL(wp),INTENT(IN),DIMENSION(nysg:nyng,nxlg:nxrg)          :: sswst_kc !<
    REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1)         :: flux_s   !<
    REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1)         :: diss_s   !<

    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) :: flux_l   !<
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) :: diss_l   !<

!-- local variables

    INTEGER :: k

    LOGICAL :: test
    s_init = 0.1                                    !bK 
    !
    !--    Tendency-terms for reactive scalar
    tend(:,j,i) = 0.0_wp

!    
!-- Advection terms
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( ws_scheme_sca )  THEN
          CALL advec_s_ws( i, j, rs, 'kc', flux_s, diss_s,           &
             flux_l, diss_l, i_omp_start, tn )
       ELSE
          CALL advec_s_pw( i, j, rs )
       ENDIF
    ELSE
       CALL advec_s_up( i, j, rs )
    ENDIF

!
!-- Diffusion terms ( ie hinteren 3 sind 0 )
!!!     CALL diffusion_s( i, j, rs, sswrc, sswrct, wall_rcflux )     
    CALL diffusion_s( i, j, rs, ssws_kc, sswst_kc, wall_sflux )   !kk use ssws and sswst instead of sswrc and sswrct   OK?

!    
!-- Prognostic equation for scalar
    DO  k = nzb_s_inner(j,i)+1, nzt
       rs_p(k,j,i) = rs(k,j,i) + dt_3d  * ( tsc(2) * tend(k,j,i) + tsc(3) * trs_m(k,j,i) )            &
                               - tsc(5) * rdf_sc(k) * ( rs(k,j,i) - s_init(k) )  !kk rs_init wird s_init OK?
       IF ( rs_p(k,j,i) < 0.0_wp )  rs_p(k,j,i) = 0.1_wp * rs(k,j,i)    !FKS6
    ENDDO

!
!-- Calculate tendencies for the next Runge-Kutta step
    IF ( timestep_scheme(1:5) == 'runge' )  THEN
       IF ( intermediate_timestep_count == 1 )  THEN
          DO  k = nzb_s_inner(j,i)+1, nzt
             trs_m(k,j,i) = tend(k,j,i)
          ENDDO
       ELSEIF ( intermediate_timestep_count < &
          intermediate_timestep_count_max )  THEN
          DO  k = nzb_s_inner(j,i)+1, nzt
             trs_m(k,j,i) = -9.5625_wp * tend(k,j,i) + &
                5.3125_wp * trs_m(k,j,i)
          ENDDO
       ENDIF
    ENDIF
 
 END SUBROUTINE


 END MODULE prognostic_equations_mod
