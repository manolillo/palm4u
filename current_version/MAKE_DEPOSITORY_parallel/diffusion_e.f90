!> @file diffusion_e.f90
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
! $Id: diffusion_e.f90 2038 2016-10-26 11:16:56Z knoop $
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! Adapted for modularization of microphysics
! 
! 1831 2016-04-07 13:15:51Z hoffmann 
! turbulence renamed collision_turbulence
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY list
! rif removed from acc-present-list
! 
! 1340 2014-03-25 19:45:13Z kanani
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1257 2013-11-08 15:18:40Z raasch
! openacc loop and loop vector clauses removed
!
! 1179 2013-06-14 05:57:58Z raasch
! use_reference renamed use_single_reference_value
!
! 1171 2013-05-30 11:27:45Z raasch
! use_reference-case activated in accelerator version
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1065 2012-11-22 17:42:36Z hoffmann
! Enabled the claculation of diss in case of turbulence = .TRUE. (parameterized 
! effects of turbulence on autoconversion and accretion in two-moments cloud 
! physics scheme). 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator version (*_acc) added,
! adjustment of mixing length to the Prandtl mixing length at first grid point
! above ground removed
!
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! 1001 2012-09-13 14:08:46Z raasch
! most arrays comunicated by module instead of parameter list
!
! 825 2012-02-19 03:03:44Z raasch
! wang_collision_kernel renamed wang_kernel
!
! Revision 1.1  1997/09/19 07:40:24  raasch
! Initial revision
!
!
! Description:
! ------------
!> Diffusion- and dissipation terms for the TKE
!------------------------------------------------------------------------------!
 MODULE diffusion_e_mod
 

    PRIVATE
    PUBLIC diffusion_e, diffusion_e_acc
    

    INTERFACE diffusion_e
       MODULE PROCEDURE diffusion_e
       MODULE PROCEDURE diffusion_e_ij
    END INTERFACE diffusion_e
 
    INTERFACE diffusion_e_acc
       MODULE PROCEDURE diffusion_e_acc
    END INTERFACE diffusion_e_acc

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_e( var, var_reference )

       USE arrays_3d,                                                          &
           ONLY:  dd2zu, ddzu, ddzw, diss, e, km, l_grid, tend, zu, zw,        &
                  drho_air, rho_air_zw
           
       USE control_parameters,                                                 &
           ONLY:  atmos_ocean_sign, g, use_single_reference_value, &
                  wall_adjustment, wall_adjustment_factor

       USE grid_variables,                                                     &
           ONLY:  ddx2, ddy2
           
       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzb_s_inner,&
                  nzt
           
       USE kinds

       USE microphysics_mod,                                                   &
           ONLY:  collision_turbulence

       USE particle_attributes,                                                &
           ONLY:  use_sgs_for_particles, wang_kernel

       IMPLICIT NONE

       INTEGER(iwp) ::  i              !< 
       INTEGER(iwp) ::  j              !< 
       INTEGER(iwp) ::  k              !< 
       REAL(wp)     ::  dvar_dz        !< 
       REAL(wp)     ::  l_stable       !< 
       REAL(wp)     ::  var_reference  !< 

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !< 
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  var  !< 
#endif
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  dissipation  !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  l            !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  ll           !< 
 

!
!--    This if clause must be outside the k-loop because otherwise
!--    runtime errors occur with -C hopt on NEC
       IF ( use_single_reference_value )  THEN

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_s_inner(j,i)+1, nzt
!
!--                Calculate the mixing length (for dissipation)
                   dvar_dz = atmos_ocean_sign * &
                             ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)
                   IF ( dvar_dz > 0.0_wp ) THEN
                      l_stable = 0.76_wp * SQRT( e(k,j,i) ) / &
                                    SQRT( g / var_reference * dvar_dz ) + 1E-5_wp
                   ELSE
                      l_stable = l_grid(k)
                   ENDIF
!
!--                Adjustment of the mixing length
                   IF ( wall_adjustment )  THEN
                      l(k,j)  = MIN( wall_adjustment_factor *          &
                                     ( zu(k) - zw(nzb_s_inner(j,i)) ), &
                                     l_grid(k), l_stable )
                      ll(k,j) = MIN( wall_adjustment_factor *          &
                                     ( zu(k) - zw(nzb_s_inner(j,i)) ), &
                                     l_grid(k) )
                   ELSE
                      l(k,j)  = MIN( l_grid(k), l_stable )
                      ll(k,j) = l_grid(k)
                   ENDIF

                ENDDO
             ENDDO

!
!--          Calculate the tendency terms
             DO  j = nys, nyn
                DO  k = nzb_s_inner(j,i)+1, nzt

                    dissipation(k,j) = ( 0.19_wp + 0.74_wp * l(k,j) / ll(k,j) ) * &
                                       e(k,j,i) * SQRT( e(k,j,i) ) / l(k,j)

                    tend(k,j,i) = tend(k,j,i)                                  &
                                        + (                                    &
                          ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )  &
                        - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )  &
                                          ) * ddx2                             &
                                        + (                                    &
                          ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )  &
                        - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )  &
                                          ) * ddy2                             &
                                        + (                                    &
               ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1) &
                                                             * rho_air_zw(k)   &
             - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)   &
                                                             * rho_air_zw(k-1) &
                                          ) * ddzw(k) * drho_air(k)            &
                             - dissipation(k,j)

                ENDDO
             ENDDO

!
!--          Store dissipation if needed for calculating the sgs particle
!--          velocities
             IF ( use_sgs_for_particles  .OR.  wang_kernel  .OR.               &
                  collision_turbulence )  THEN
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      diss(k,j,i) = dissipation(k,j)
                   ENDDO
                ENDDO
             ENDIF

          ENDDO

       ELSE

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_s_inner(j,i)+1, nzt
!
!--                Calculate the mixing length (for dissipation)
                   dvar_dz = atmos_ocean_sign * &
                             ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)
                   IF ( dvar_dz > 0.0_wp ) THEN
                      l_stable = 0.76_wp * SQRT( e(k,j,i) ) / &
                                           SQRT( g / var(k,j,i) * dvar_dz ) + 1E-5_wp
                   ELSE
                      l_stable = l_grid(k)
                   ENDIF
!
!--                Adjustment of the mixing length
                   IF ( wall_adjustment )  THEN
                      l(k,j)  = MIN( wall_adjustment_factor *          &
                                     ( zu(k) - zw(nzb_s_inner(j,i)) ), &
                                     l_grid(k), l_stable )
                      ll(k,j) = MIN( wall_adjustment_factor *          &
                                     ( zu(k) - zw(nzb_s_inner(j,i)) ), &
                                     l_grid(k) )
                   ELSE
                      l(k,j)  = MIN( l_grid(k), l_stable )
                      ll(k,j) = l_grid(k)
                   ENDIF

                ENDDO
             ENDDO

!
!--          Calculate the tendency terms
             DO  j = nys, nyn
                DO  k = nzb_s_inner(j,i)+1, nzt

                    dissipation(k,j) = ( 0.19_wp + 0.74_wp * l(k,j) / ll(k,j) ) * &
                                             e(k,j,i) * SQRT( e(k,j,i) ) / l(k,j)

                    tend(k,j,i) = tend(k,j,i)                                  &
                                        + (                                    &
                          ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )  &
                        - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )  &
                                          ) * ddx2                             &
                                        + (                                    &
                          ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )  &
                        - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )  &
                                          ) * ddy2                             &
                                        + (                                    &
               ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1) &
                                                             * rho_air_zw(k)   &
             - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)   &
                                                             * rho_air_zw(k-1) &
                                          ) * ddzw(k) * drho_air(k)            &
                             - dissipation(k,j)

                ENDDO
             ENDDO

!
!--          Store dissipation if needed for calculating the sgs particle
!--          velocities
             IF ( use_sgs_for_particles  .OR.  wang_kernel  .OR.               &
                  collision_turbulence )  THEN
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzt
                      diss(k,j,i) = dissipation(k,j)
                   ENDDO
                ENDDO
             ENDIF

          ENDDO

       ENDIF

!
!--    Boundary condition for dissipation
       IF ( use_sgs_for_particles  .OR.  wang_kernel  .OR.  &
            collision_turbulence )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                diss(nzb_s_inner(j,i),j,i) = diss(nzb_s_inner(j,i)+1,j,i)
             ENDDO
          ENDDO
       ENDIF

    END SUBROUTINE diffusion_e


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points - accelerator version
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_e_acc( var, var_reference )

       USE arrays_3d,                                                          &
           ONLY:  dd2zu, ddzu, ddzw, diss, e, km, l_grid, tend, zu, zw,        &
                  drho_air, rho_air_zw
          
       USE control_parameters,                                                 &
           ONLY:  atmos_ocean_sign, g, use_single_reference_value,             &
                  wall_adjustment, wall_adjustment_factor

       USE grid_variables,                                                     &
           ONLY:  ddx2, ddy2
           
       USE indices,                                                            &
           ONLY:  i_left, i_right, j_north, j_south, nxlg, nxrg, nyng, nysg,   &
                  nzb, nzb_s_inner, nzt
           
       USE kinds

       USE microphysics_mod,                                                   &
           ONLY:  collision_turbulence

       USE particle_attributes,                                                &
           ONLY:  use_sgs_for_particles, wang_kernel

       IMPLICIT NONE

       INTEGER(iwp) ::  i              !< 
       INTEGER(iwp) ::  j              !< 
       INTEGER(iwp) ::  k              !< 
       REAL(wp)     ::  dissipation    !< 
       REAL(wp)     ::  dvar_dz        !< 
       REAL(wp)     ::  l              !< 
       REAL(wp)     ::  ll             !< 
       REAL(wp)     ::  l_stable       !< 
       REAL(wp)     ::  var_reference  !< 

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !< 
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  var  !< 
#endif


!
!--    This if clause must be outside the k-loop because otherwise
!--    runtime errors occur with -C hopt on NEC
       IF ( use_single_reference_value )  THEN

          !$acc kernels present( ddzu, ddzw, dd2zu, diss, e, km, l_grid ) &
          !$acc         present( nzb_s_inner, tend, var, zu, zw )
          DO  i = i_left, i_right
             DO  j = j_south, j_north
                DO  k = 1, nzt

                   IF ( k > nzb_s_inner(j,i) )  THEN
!
!--                   Calculate the mixing length (for dissipation)
                      dvar_dz = atmos_ocean_sign * &
                                ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)
                      IF ( dvar_dz > 0.0_wp ) THEN
                         l_stable = 0.76_wp * SQRT( e(k,j,i) ) / &
                                       SQRT( g / var_reference * dvar_dz ) + 1E-5_wp
                      ELSE
                         l_stable = l_grid(k)
                      ENDIF
!
!--                   Adjustment of the mixing length
                      IF ( wall_adjustment )  THEN
                         l  = MIN( wall_adjustment_factor *          &
                                   ( zu(k) - zw(nzb_s_inner(j,i)) ), &
                                   l_grid(k), l_stable )
                         ll = MIN( wall_adjustment_factor *          &
                                   ( zu(k) - zw(nzb_s_inner(j,i)) ), &
                                   l_grid(k) )
                      ELSE
                         l  = MIN( l_grid(k), l_stable )
                         ll = l_grid(k)
                      ENDIF
!
!--                   Calculate the tendency terms
                      dissipation = ( 0.19_wp + 0.74_wp * l / ll ) * &
                                          e(k,j,i) * SQRT( e(k,j,i) ) / l

                      tend(k,j,i) = tend(k,j,i)                                  &
                                        + (                                    &
                          ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )  &
                        - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )  &
                                          ) * ddx2                             &
                                        + (                                    &
                          ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )  &
                        - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )  &
                                          ) * ddy2                             &
                                        + (                                    &
               ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1) &
                                                             * rho_air_zw(k)   &
             - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)   &
                                                             * rho_air_zw(k-1) &
                                          ) * ddzw(k) * drho_air(k)            &
                                  - dissipation

!
!--                   Store dissipation if needed for calculating the sgs particle
!--                   velocities
                      IF ( use_sgs_for_particles  .OR.  wang_kernel  .OR.      &
                           collision_turbulence )  THEN
                         diss(k,j,i) = dissipation
                      ENDIF

                   ENDIF

                ENDDO
             ENDDO
          ENDDO
          !$acc end kernels

       ELSE

          !$acc kernels present( ddzu, ddzw, dd2zu, diss, e, km, l_grid ) &
          !$acc         present( nzb_s_inner, tend, var, zu, zw )
          DO  i = i_left, i_right
             DO  j = j_south, j_north
                DO  k = 1, nzt

                   IF ( k > nzb_s_inner(j,i) )  THEN
!
!--                   Calculate the mixing length (for dissipation)
                      dvar_dz = atmos_ocean_sign * &
                                ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)
                      IF ( dvar_dz > 0.0_wp ) THEN
                         l_stable = 0.76_wp * SQRT( e(k,j,i) ) / &
                                              SQRT( g / var(k,j,i) * dvar_dz ) + 1E-5_wp
                      ELSE
                         l_stable = l_grid(k)
                      ENDIF
!
!--                   Adjustment of the mixing length
                      IF ( wall_adjustment )  THEN
                         l  = MIN( wall_adjustment_factor *          &
                                   ( zu(k) - zw(nzb_s_inner(j,i)) ), &
                                     l_grid(k), l_stable )
                         ll = MIN( wall_adjustment_factor *          &
                                   ( zu(k) - zw(nzb_s_inner(j,i)) ), &
                                   l_grid(k) )
                      ELSE
                         l  = MIN( l_grid(k), l_stable )
                         ll = l_grid(k)
                      ENDIF
!
!--                   Calculate the tendency terms
                      dissipation = ( 0.19_wp + 0.74_wp * l / ll ) * &
                                          e(k,j,i) * SQRT( e(k,j,i) ) / l

                      tend(k,j,i) = tend(k,j,i)                                &
                                        + (                                    &
                          ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )  &
                        - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )  &
                                          ) * ddx2                             &
                                        + (                                    &
                          ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )  &
                        - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )  &
                                          ) * ddy2                             &
                                        + (                                    &
               ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1) &
                                                             * rho_air_zw(k)   &
             - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)   &
                                                             * rho_air_zw(k-1) &
                                          ) * ddzw(k) * drho_air(k)            &
                                  - dissipation

!
!--                   Store dissipation if needed for calculating the sgs
!--                   particle  velocities
                      IF ( use_sgs_for_particles  .OR.  wang_kernel  .OR.      &
                           collision_turbulence )  THEN
                         diss(k,j,i) = dissipation
                      ENDIF

                   ENDIF

                ENDDO
             ENDDO
          ENDDO
          !$acc end kernels

       ENDIF

!
!--    Boundary condition for dissipation
       IF ( use_sgs_for_particles  .OR.  wang_kernel  .OR.                     &
            collision_turbulence )  THEN
          !$acc kernels present( diss, nzb_s_inner )
          DO  i = i_left, i_right
             DO  j = j_south, j_north
                diss(nzb_s_inner(j,i),j,i) = diss(nzb_s_inner(j,i)+1,j,i)
             ENDDO
          ENDDO
          !$acc end kernels
       ENDIF

    END SUBROUTINE diffusion_e_acc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_e_ij( i, j, var, var_reference )

       USE arrays_3d,                                                          &
           ONLY:  dd2zu, ddzu, ddzw, diss, e, km, l_grid, tend, zu, zw,        &
                  drho_air, rho_air_zw
          
       USE control_parameters,                                                 &
           ONLY:  atmos_ocean_sign, g, use_single_reference_value,             &
                  wall_adjustment, wall_adjustment_factor

       USE grid_variables,                                                     &
           ONLY:  ddx2, ddy2
           
       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzb_s_inner, nzt
           
       USE kinds

       USE microphysics_mod,                                                   &
           ONLY:  collision_turbulence

       USE particle_attributes,                                                &
           ONLY:  use_sgs_for_particles, wang_kernel

       IMPLICIT NONE

       INTEGER(iwp) ::  i              !< 
       INTEGER(iwp) ::  j              !< 
       INTEGER(iwp) ::  k              !< 
       REAL(wp)     ::  dvar_dz        !< 
       REAL(wp)     ::  l_stable       !<
       REAL(wp)     ::  var_reference  !<

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  var     !<
#endif
       REAL(wp), DIMENSION(nzb+1:nzt) ::  dissipation  !<
       REAL(wp), DIMENSION(nzb+1:nzt) ::  l            !<
       REAL(wp), DIMENSION(nzb+1:nzt) ::  ll           !<


!
!--    Calculate the mixing length (for dissipation)
       DO  k = nzb_s_inner(j,i)+1, nzt
          dvar_dz = atmos_ocean_sign * &
                    ( var(k+1,j,i) - var(k-1,j,i) ) * dd2zu(k)
          IF ( dvar_dz > 0.0_wp ) THEN
             IF ( use_single_reference_value )  THEN
                l_stable = 0.76_wp * SQRT( e(k,j,i) ) / &
                                     SQRT( g / var_reference * dvar_dz ) + 1E-5_wp
             ELSE
                l_stable = 0.76_wp * SQRT( e(k,j,i) ) / &
                                     SQRT( g / var(k,j,i) * dvar_dz ) + 1E-5_wp
             ENDIF
          ELSE
             l_stable = l_grid(k)
          ENDIF
!
!--       Adjustment of the mixing length
          IF ( wall_adjustment )  THEN
             l(k)  = MIN( wall_adjustment_factor *                     &
                          ( zu(k) - zw(nzb_s_inner(j,i)) ), l_grid(k), &
                          l_stable )
             ll(k) = MIN( wall_adjustment_factor *                     &
                          ( zu(k) - zw(nzb_s_inner(j,i)) ), l_grid(k) )
          ELSE
             l(k)  = MIN( l_grid(k), l_stable )
             ll(k) = l_grid(k)
          ENDIF
!
!--       Calculate the tendency term
          dissipation(k) = ( 0.19_wp + 0.74_wp * l(k) / ll(k) ) * e(k,j,i) * &
                                 SQRT( e(k,j,i) ) / l(k)

          tend(k,j,i) = tend(k,j,i)                                           &
                                       + (                                    &
                         ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )  &
                       - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )  &
                                         ) * ddx2                             &
                                       + (                                    &
                         ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )  &
                       - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )  &
                                         ) * ddy2                             &
                                       + (                                    &
              ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1) &
                                                            * rho_air_zw(k)   &
            - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)   &
                                                            * rho_air_zw(k-1) &
                                         ) * ddzw(k) * drho_air(k)            &
                                       - dissipation(k)

       ENDDO

!
!--    Store dissipation if needed for calculating the sgs particle velocities
       IF ( use_sgs_for_particles  .OR.  wang_kernel  .OR.                     &
            collision_turbulence )  THEN
          DO  k = nzb_s_inner(j,i)+1, nzt
             diss(k,j,i) = dissipation(k)
          ENDDO
!
!--       Boundary condition for dissipation
          diss(nzb_s_inner(j,i),j,i) = diss(nzb_s_inner(j,i)+1,j,i)
       ENDIF

    END SUBROUTINE diffusion_e_ij

 END MODULE diffusion_e_mod
