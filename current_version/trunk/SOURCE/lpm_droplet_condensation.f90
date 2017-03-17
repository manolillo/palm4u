!> @file lpm_droplet_condensation.f90
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
! 
! 
! Former revisions:
! -----------------
! $Id: lpm_droplet_condensation.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1890 2016-04-22 08:52:11Z hoffmann
! Some improvements of the Rosenbrock method. If the Rosenbrock method needs more
! than 40 iterations to find a sufficient time setp, the model is not aborted.
! This might lead to small erros. But they can be assumend as negligible, since 
! the maximum timesetp of the Rosenbrock method after 40 iterations will be 
! smaller than 10^-16 s.  
! 
! 1871 2016-04-15 11:46:09Z hoffmann
! Initialization of aerosols added.
!
! 1849 2016-04-08 11:33:18Z hoffmann
! Interpolation of supersaturation has been removed because it is not in
! accordance with the release/depletion of latent heat/water vapor in
! interaction_droplets_ptq.
! Calculation of particle Reynolds number has been corrected.
! eps_ros added from modules.
!
! 1831 2016-04-07 13:15:51Z hoffmann
! curvature_solution_effects moved to particle_attributes
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
! Kind definition added to all floating point numbers.
!
! 1346 2014-03-27 13:18:20Z heinze
! Bugfix: REAL constants provided with KIND-attribute especially in call of 
! intrinsic function like MAX, MIN, SIGN
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1318 2014-03-17 13:35:16Z raasch
! module interfaces removed
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1071 2012-11-29 16:54:55Z franke
! Ventilation effect for evaporation of large droplets included
! Check for unreasonable results included in calculation of Rosenbrock method
! since physically unlikely results were observed and for the same
! reason the first internal time step in Rosenbrock method should be < 1.0E02 in
! case of evaporation
! Unnecessary calculation of ql_int removed
! Unnecessary calculations in Rosenbrock method (d2rdt2, drdt_m, dt_ros_last)
! removed
! Bugfix: factor in calculation of surface tension changed from 0.00155 to
! 0.000155
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! initial revision (former part of advec_particles)
!
!
! Description:
! ------------
!> Calculates change in droplet radius by condensation/evaporation, using
!> either an analytic formula or by numerically integrating the radius growth
!> equation including curvature and solution effects using Rosenbrocks method
!> (see Numerical recipes in FORTRAN, 2nd edition, p. 731).
!> The analytical formula and growth equation follow those given in
!> Rogers and Yau (A short course in cloud physics, 3rd edition, p. 102/103).
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_droplet_condensation (ip,jp,kp)
 

    USE arrays_3d,                                                             &
        ONLY:  hyp, pt, q,  ql_c, ql_v

    USE cloud_parameters,                                                      &
        ONLY:  l_d_rv, l_v, rho_l, r_v

    USE constants,                                                             &
        ONLY:  pi

    USE control_parameters,                                                    &
        ONLY:  dt_3d, dz, message_string, molecular_viscosity, rho_surface

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE lpm_collision_kernels_mod,                                             &
        ONLY:  rclass_lbound, rclass_ubound

    USE kinds

    USE particle_attributes,                                                   &
        ONLY:  curvature_solution_effects, hall_kernel,                        &
               molecular_weight_of_solute, molecular_weight_of_water,          &
               number_of_particles, particles, radius_classes, rho_s,          &
               use_kernel_tables, vanthoff, wang_kernel


    IMPLICIT NONE

    INTEGER(iwp) :: ip                         !<
    INTEGER(iwp) :: internal_timestep_count    !<
    INTEGER(iwp) :: jp                         !<
    INTEGER(iwp) :: jtry                       !<
    INTEGER(iwp) :: kp                         !<
    INTEGER(iwp) :: n                          !<
    INTEGER(iwp) :: ros_count                  !<
 
    INTEGER(iwp), PARAMETER ::  maxtry = 40    !<

    LOGICAL ::  repeat                         !<

    REAL(wp) ::  aa                            !<
    REAL(wp) ::  afactor                       !< curvature effects
    REAL(wp) ::  arg                           !<
    REAL(wp) ::  bfactor                       !< solute effects
    REAL(wp) ::  ddenom                        !<
    REAL(wp) ::  delta_r                       !<
    REAL(wp) ::  diameter                      !< diameter of cloud droplets
    REAL(wp) ::  diff_coeff_v                  !< diffusivity for water vapor
    REAL(wp) ::  drdt                          !<
    REAL(wp) ::  drdt_ini                      !<
    REAL(wp) ::  dt_ros                        !<
    REAL(wp) ::  dt_ros_next                   !<
    REAL(wp) ::  dt_ros_sum                    !<
    REAL(wp) ::  dt_ros_sum_ini                !<
    REAL(wp) ::  d2rdtdr                       !<
    REAL(wp) ::  errmax                        !<
    REAL(wp) ::  e_a                           !< current vapor pressure
    REAL(wp) ::  e_s                           !< current saturation vapor pressure
    REAL(wp) ::  err_ros                       !<
    REAL(wp) ::  g1                            !<
    REAL(wp) ::  g2                            !<
    REAL(wp) ::  g3                            !<
    REAL(wp) ::  g4                            !<
    REAL(wp) ::  r_ros                         !<
    REAL(wp) ::  r_ros_ini                     !<
    REAL(wp) ::  sigma                         !<
    REAL(wp) ::  thermal_conductivity_v        !< thermal conductivity for water
    REAL(wp) ::  t_int                         !< temperature
    REAL(wp) ::  w_s                           !< terminal velocity of droplets
    REAL(wp) ::  re_p                          !<

!
!-- Parameters for Rosenbrock method
    REAL(wp), PARAMETER ::  a21 = 2.0_wp               !<
    REAL(wp), PARAMETER ::  a31 = 48.0_wp / 25.0_wp    !<
    REAL(wp), PARAMETER ::  a32 = 6.0_wp / 25.0_wp     !<
    REAL(wp), PARAMETER ::  b1 = 19.0_wp / 9.0_wp      !<
    REAL(wp), PARAMETER ::  b2 = 0.5_wp                !<
    REAL(wp), PARAMETER ::  b3 = 25.0_wp / 108.0_wp    !<
    REAL(wp), PARAMETER ::  b4 = 125.0_wp / 108.0_wp   !<
    REAL(wp), PARAMETER ::  c21 = -8.0_wp              !<
    REAL(wp), PARAMETER ::  c31 = 372.0_wp / 25.0_wp   !<
    REAL(wp), PARAMETER ::  c32 = 12.0_wp / 5.0_wp     !<
    REAL(wp), PARAMETER ::  c41 = -112.0_wp / 125.0_wp !<
    REAL(wp), PARAMETER ::  c42 = -54.0_wp / 125.0_wp  !<
    REAL(wp), PARAMETER ::  c43 = -2.0_wp / 5.0_wp     !<
    REAL(wp), PARAMETER ::  errcon = 0.1296_wp         !<
    REAL(wp), PARAMETER ::  e1 = 17.0_wp / 54.0_wp     !<
    REAL(wp), PARAMETER ::  e2 = 7.0_wp / 36.0_wp      !<
    REAL(wp), PARAMETER ::  e3 = 0.0_wp                !<
    REAL(wp), PARAMETER ::  e4 = 125.0_wp / 108.0_wp   !<
    REAL(wp), PARAMETER ::  eps_ros = 1.0E-3_wp        !< accuracy of Rosenbrock method
    REAL(wp), PARAMETER ::  gam = 0.5_wp               !<
    REAL(wp), PARAMETER ::  grow = 1.5_wp              !<
    REAL(wp), PARAMETER ::  pgrow = -0.25_wp           !<
    REAL(wp), PARAMETER ::  pshrnk = -1.0_wp /3.0_wp   !<
    REAL(wp), PARAMETER ::  shrnk = 0.5_wp             !<

!
!-- Parameters for terminal velocity
    REAL(wp), PARAMETER ::  a_rog = 9.65_wp      !< parameter for fall velocity
    REAL(wp), PARAMETER ::  b_rog = 10.43_wp     !< parameter for fall velocity
    REAL(wp), PARAMETER ::  c_rog = 0.6_wp       !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_cap_rog = 4.0_wp   !< parameter for fall velocity
    REAL(wp), PARAMETER ::  k_low_rog = 12.0_wp  !< parameter for fall velocity
    REAL(wp), PARAMETER ::  d0_rog = 0.745_wp    !< separation diameter

    REAL(wp), DIMENSION(number_of_particles) ::  ventilation_effect     !<
    REAL(wp), DIMENSION(number_of_particles) ::  new_r                  !<



    CALL cpu_log( log_point_s(42), 'lpm_droplet_condens', 'start' )

!
!-- Calculate temperature, saturation vapor pressure and current vapor pressure
    t_int = pt(kp,jp,ip) * ( hyp(kp) / 100000.0_wp )**0.286_wp
    e_s   = 611.0_wp * EXP( l_d_rv * ( 3.6609E-3_wp - 1.0_wp / t_int ) )
    e_a   = q(kp,jp,ip) * hyp(kp) / ( 0.378_wp * q(kp,jp,ip) + 0.622_wp )
!
!-- Thermal conductivity for water (from Rogers and Yau, Table 7.1),
!-- diffusivity for water vapor (after Hall und Pruppacher, 1976)
    thermal_conductivity_v = 7.94048E-05_wp * t_int + 0.00227011_wp
    diff_coeff_v           = 0.211E-4_wp * ( t_int / 273.15_wp )**1.94_wp * &
                             ( 101325.0_wp / hyp(kp) )
!
!-- Calculate effects of heat conductivity and diffusion of water vapor on the 
!-- condensation/evaporation process (typically known as 1.0 / (F_k + F_d) )
    ddenom  = 1.0_wp / ( rho_l * r_v * t_int / ( e_s * diff_coeff_v ) +        &
                         ( l_v / ( r_v * t_int ) - 1.0_wp ) * rho_l *          &
                         l_v / ( thermal_conductivity_v * t_int )              &
                       )

    new_r = 0.0_wp

!
!-- Determine ventilation effect on evaporation of large drops
    DO  n = 1, number_of_particles

       IF ( particles(n)%radius >= 4.0E-5_wp  .AND.  e_a / e_s < 1.0_wp )  THEN
!
!--       Terminal velocity is computed for vertical direction (Rogers et al.,
!--       1993, J. Appl. Meteorol.)
          diameter = particles(n)%radius * 2000.0_wp !diameter in mm
          IF ( diameter <= d0_rog )  THEN
             w_s = k_cap_rog * diameter * ( 1.0_wp - EXP( -k_low_rog * diameter ) )
          ELSE
             w_s = a_rog - b_rog * EXP( -c_rog * diameter )
          ENDIF
!
!--       First calculate droplet's Reynolds number
          re_p = 2.0_wp * particles(n)%radius * w_s / molecular_viscosity
!
!--       Ventilation coefficient (Rogers and Yau, 1989):
          IF ( re_p > 2.5_wp )  THEN
             ventilation_effect(n) = 0.78_wp + 0.28_wp * SQRT( re_p )
          ELSE
             ventilation_effect(n) = 1.0_wp + 0.09_wp * re_p
          ENDIF
       ELSE
!
!--       For small droplets or in supersaturated environments, the ventilation 
!--       effect does not play a role
          ventilation_effect(n) = 1.0_wp
       ENDIF
    ENDDO

!
!-- Use analytic model for condensational growth
    IF( .NOT. curvature_solution_effects ) then
       DO  n = 1, number_of_particles
          arg = particles(n)%radius**2 + 2.0_wp * dt_3d * ddenom *             &
                                         ventilation_effect(n) *               &
                                         ( e_a / e_s - 1.0_wp )
          arg = MAX( arg, 1.0E-16_wp )
          new_r(n) = SQRT( arg )
       ENDDO
    ENDIF

!
!-- If selected, use numerical solution of the condensational growth
!-- equation (e.g., for studying the activation of aerosols).
!-- Curvature and solutions effects are included in growth equation.
!-- Change in Radius is calculated with a 4th-order Rosenbrock method
!-- for stiff o.d.e's with monitoring local truncation error to adjust
!-- stepsize (see Numerical recipes in FORTRAN, 2nd edition, p. 731).
    DO  n = 1, number_of_particles
       IF ( curvature_solution_effects )  THEN

          ros_count = 0
          repeat = .TRUE.
!
!--       Carry out the Rosenbrock algorithm. In case of unreasonable results
!--       the switch "repeat" will be set true and the algorithm will be carried
!--       out again with the internal time step set to its initial (small) value.
!--       Unreasonable results may occur if the external conditions, especially 
!--       the supersaturation, has significantly changed compared to the last 
!--       PALM timestep.
          DO WHILE ( repeat )

             repeat = .FALSE.
!
!--          Curvature effect (afactor) with surface tension parameterization
!--          by Straka (2009)
             sigma = 0.0761_wp - 0.000155_wp * ( t_int - 273.15_wp )
             afactor = 2.0_wp * sigma / ( rho_l * r_v * t_int )
!
!--          Solute effect (bfactor)
             bfactor = vanthoff * rho_s * particles(n)%rvar2**3 *              &
                       molecular_weight_of_water /                             &
                       ( rho_l * molecular_weight_of_solute )

             r_ros = particles(n)%radius
             dt_ros_sum  = 0.0_wp      ! internal integrated time (s)
             internal_timestep_count = 0
!
!--          Take internal time step values from the end of last PALM time step
             dt_ros_next = particles(n)%rvar1

!
!--          Internal time step should not be > 1.0E-2 and < 1.0E-6 
             IF ( dt_ros_next > 1.0E-2_wp )  THEN
                dt_ros_next = 1.0E-2_wp
             ELSEIF ( dt_ros_next < 1.0E-6_wp )  THEN
                dt_ros_next = 1.0E-6_wp
             ENDIF

!
!--          If calculation of Rosenbrock method is repeated due to unreasonalble
!--          results during previous try the initial internal time step has to be
!--          reduced
             IF ( ros_count > 1 )  THEN
                dt_ros_next = dt_ros_next * 0.1_wp
             ELSEIF ( ros_count > 5 )  THEN
!
!--             Prevent creation of infinite loop
                message_string = 'ros_count > 5 in Rosenbrock method'
                CALL message( 'lpm_droplet_condensation', 'PA0018', 2, 2, &
                               0, 6, 0 )
             ENDIF

!
!--          Internal time step must not be larger than PALM time step
             dt_ros = MIN( dt_ros_next, dt_3d )

!
!--          Integrate growth equation in time unless PALM time step is reached
             DO WHILE ( dt_ros_sum < dt_3d )

                internal_timestep_count = internal_timestep_count + 1

!
!--             Derivative at starting value
                drdt = ddenom * ventilation_effect(n) * ( e_a / e_s - 1.0_wp - &
                                                          afactor / r_ros +    &
                                                          bfactor / r_ros**3   &
                                                        ) / r_ros

                drdt_ini       = drdt
                dt_ros_sum_ini = dt_ros_sum
                r_ros_ini      = r_ros

!
!--             Calculate radial derivative of dr/dt
                d2rdtdr = ddenom * ventilation_effect(n) *                     &
                                       ( ( 1.0_wp - e_a / e_s ) / r_ros**2 +   &
                                         2.0_wp * afactor / r_ros**3 -         &
                                         4.0_wp * bfactor / r_ros**5           &
                                       )
!
!--             Adjust stepsize unless required accuracy is reached
                DO  jtry = 1, maxtry+1

                   IF ( jtry == maxtry+1 )  THEN
                      message_string = 'maxtry > 40 in Rosenbrock method'
                      CALL message( 'lpm_droplet_condensation', 'PA0347', 0,   &
                                    1, 0, 6, 0 )
                   ENDIF

                   aa    = 1.0_wp / ( gam * dt_ros ) - d2rdtdr
                   g1    = drdt_ini / aa
                   r_ros = r_ros_ini + a21 * g1
                   drdt  = ddenom * ventilation_effect(n) * ( e_a / e_s - 1.0_wp - &
                                                              afactor / r_ros +    &
                                                              bfactor / r_ros**3   &
                                                            ) / r_ros

                   g2    = ( drdt + c21 * g1 / dt_ros )&
                           / aa
                   r_ros = r_ros_ini + a31 * g1 + a32 * g2
                   drdt  = ddenom * ventilation_effect(n) * ( e_a / e_s - 1.0_wp - &
                                                              afactor / r_ros +    &
                                                              bfactor / r_ros**3   &
                                                            ) / r_ros

                   g3    = ( drdt +  &
                             ( c31 * g1 + c32 * g2 ) / dt_ros ) / aa
                   g4    = ( drdt +  &
                             ( c41 * g1 + c42 * g2 + c43 * g3 ) / dt_ros ) / aa
                   r_ros = r_ros_ini + b1 * g1 + b2 * g2 + b3 * g3 + b4 * g4

                   dt_ros_sum = dt_ros_sum_ini + dt_ros

                   IF ( dt_ros_sum == dt_ros_sum_ini )  THEN
                      message_string = 'zero stepsize in Rosenbrock method'
                      CALL message( 'lpm_droplet_condensation', 'PA0348', 2,   &
                                    2, 0, 6, 0 )
                   ENDIF
!
!--                Calculate error
                   err_ros = e1 * g1 + e2 * g2 + e3 * g3 + e4 * g4
                   errmax  = 0.0_wp
                   errmax  = MAX( errmax, ABS( err_ros / r_ros_ini ) ) / eps_ros
!
!--                Leave loop if accuracy is sufficient, otherwise try again
!--                with a reduced stepsize
                   IF ( errmax <= 1.0_wp )  THEN
                      EXIT
                   ELSE
                      dt_ros = MAX( ABS( 0.9_wp * dt_ros * errmax**pshrnk ),   &
                                    shrnk * ABS( dt_ros ) )
                   ENDIF

                ENDDO  ! loop for stepsize adjustment

!
!--             Calculate next internal time step
                IF ( errmax > errcon )  THEN
                   dt_ros_next = 0.9_wp * dt_ros * errmax**pgrow
                ELSE
                   dt_ros_next = grow * dt_ros
                ENDIF

!
!--             Estimated time step is reduced if the PALM time step is exceeded
                IF ( ( dt_ros_next + dt_ros_sum ) >= dt_3d )  THEN
                   dt_ros = dt_3d - dt_ros_sum
                ELSE
                   dt_ros = dt_ros_next
                ENDIF

             ENDDO
!
!--          Store internal time step value for next PALM step
             particles(n)%rvar1 = dt_ros_next

!
!--          Radius should not fall below 1E-8 because Rosenbrock method may
!--          lead to errors otherwise
             new_r(n) = MAX( r_ros, particles(n)%rvar2 )
!
!--          Check if calculated droplet radius change is reasonable since in
!--          case of droplet evaporation the Rosenbrock method may lead to
!--          secondary solutions which are physically unlikely.
!--          Due to the solution effect the droplets may grow for relative
!--          humidities below 100%, but change of radius should not be too 
!--          large. In case of unreasonable droplet growth the Rosenbrock 
!--          method is recalculated using a smaller initial time step.
!--          Limiting values are tested for droplets down to 1.0E-7
             IF ( new_r(n) - particles(n)%radius >= 3.0E-7_wp  .AND.  &
                  e_a / e_s < 0.97_wp )  THEN
                ros_count = ros_count + 1
                repeat = .TRUE.
             ENDIF

          ENDDO    ! Rosenbrock method

       ENDIF

       delta_r = new_r(n) - particles(n)%radius

!
!--    Sum up the change in volume of liquid water for the respective grid
!--    volume (this is needed later in lpm_calc_liquid_water_content for
!--    calculating the release of latent heat)
       ql_c(kp,jp,ip) = ql_c(kp,jp,ip) + particles(n)%weight_factor *          &
                                   rho_l * 1.33333333_wp * pi *                &
                                   ( new_r(n)**3 - particles(n)%radius**3 ) /  &
                                   ( rho_surface * dx * dy * dz )
       IF ( ql_c(kp,jp,ip) > 100.0_wp )  THEN
          WRITE( message_string, * ) 'k=',kp,' j=',jp,' i=',ip,      &
                       ' ql_c=',ql_c(kp,jp,ip), ' &part(',n,')%wf=', &
                       particles(n)%weight_factor,' delta_r=',delta_r
          CALL message( 'lpm_droplet_condensation', 'PA0143', 2, 2, -1, 6, 1 )
       ENDIF

!
!--    Change the droplet radius
       IF ( ( new_r(n) - particles(n)%radius ) < 0.0_wp  .AND.        &
            new_r(n) < 0.0_wp )  THEN
          WRITE( message_string, * ) '#1 k=',kp,' j=',jp,' i=',ip,    &
                       ' e_s=',e_s, ' e_a=',e_a,' t_int=',t_int,      &
                       ' &delta_r=',delta_r,                          &
                       ' particle_radius=',particles(n)%radius
          CALL message( 'lpm_droplet_condensation', 'PA0144', 2, 2, -1, 6, 1 )
       ENDIF

!
!--    Sum up the total volume of liquid water (needed below for
!--    re-calculating the weighting factors)
       ql_v(kp,jp,ip) = ql_v(kp,jp,ip) + particles(n)%weight_factor * new_r(n)**3

       particles(n)%radius = new_r(n)

!
!--    Determine radius class of the particle needed for collision
       IF ( ( hall_kernel  .OR.  wang_kernel )  .AND.  use_kernel_tables )     &
       THEN
          particles(n)%class = ( LOG( new_r(n) ) - rclass_lbound ) /           &
                               ( rclass_ubound - rclass_lbound ) *             &
                               radius_classes
          particles(n)%class = MIN( particles(n)%class, radius_classes )
          particles(n)%class = MAX( particles(n)%class, 1 )
       ENDIF

    ENDDO

    CALL cpu_log( log_point_s(42), 'lpm_droplet_condens', 'stop' )


 END SUBROUTINE lpm_droplet_condensation
