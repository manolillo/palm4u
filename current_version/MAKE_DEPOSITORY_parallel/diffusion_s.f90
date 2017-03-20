!> @file diffusion_s.f90
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
! $Id: diffusion_s.f90 2038 2016-10-26 11:16:56Z knoop $
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
! 
! 
! 1691 2015-10-26 16:17:44Z maronga
! Formatting corrections.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY list
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
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator version (*_acc) added
!
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! 1001 2012-09-13 14:08:46Z raasch
! some arrays comunicated by module instead of parameter list
!
! Revision 1.1  2000/04/13 14:54:02  schroeter
! Initial revision
!
!
! Description:
! ------------
!> Diffusion term of scalar quantities (temperature and water content)
!------------------------------------------------------------------------------!
 MODULE diffusion_s_mod
 

    PRIVATE
    PUBLIC diffusion_s, diffusion_s_acc

    INTERFACE diffusion_s
       MODULE PROCEDURE diffusion_s
       MODULE PROCEDURE diffusion_s_ij
    END INTERFACE diffusion_s

    INTERFACE diffusion_s_acc
       MODULE PROCEDURE diffusion_s_acc
    END INTERFACE diffusion_s_acc

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_s( s, s_flux_b, s_flux_t, wall_s_flux )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, kh, tend, drho_air, rho_air_zw
       
       USE control_parameters,                                                 & 
           ONLY: use_surface_fluxes, use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx2, ddy2, fwxm, fwxp, fwym, fwyp, wall_w_x, wall_w_y
       
       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb,             &
                  nzb_diff_s_inner, nzb_s_inner, nzb_s_outer, nzt, nzt_diff
       
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< 
       INTEGER(iwp) ::  j                 !< 
       INTEGER(iwp) ::  k                 !< 
       REAL(wp)     ::  wall_s_flux(0:4)  !< 
       REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  s_flux_b, s_flux_t !<
#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s  !< 
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  s  !< 
#endif

       DO  i = nxl, nxr
          DO  j = nys,nyn
!
!--          Compute horizontal diffusion
             DO  k = nzb_s_outer(j,i)+1, nzt

                tend(k,j,i) = tend(k,j,i)                                      &
                                          + 0.5_wp * (                         &
                        ( kh(k,j,i) + kh(k,j,i+1) ) * ( s(k,j,i+1)-s(k,j,i) )  &
                      - ( kh(k,j,i) + kh(k,j,i-1) ) * ( s(k,j,i)-s(k,j,i-1) )  &
                                                     ) * ddx2                  &
                                          + 0.5_wp * (                         &
                        ( kh(k,j,i) + kh(k,j+1,i) ) * ( s(k,j+1,i)-s(k,j,i) )  &
                      - ( kh(k,j,i) + kh(k,j-1,i) ) * ( s(k,j,i)-s(k,j-1,i) )  &
                                                     ) * ddy2
             ENDDO

!
!--          Apply prescribed horizontal wall heatflux where necessary
             IF ( ( wall_w_x(j,i) /= 0.0_wp ) .OR. ( wall_w_y(j,i) /= 0.0_wp ) &
                )  THEN
                DO  k = nzb_s_inner(j,i)+1, nzb_s_outer(j,i)

                   tend(k,j,i) = tend(k,j,i)                                   &
                                                + ( fwxp(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j,i+1) ) * ( s(k,j,i+1)-s(k,j,i) )  &
                        + ( 1.0_wp - fwxp(j,i) ) * wall_s_flux(1)              &
                                                   -fwxm(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j,i-1) ) * ( s(k,j,i)-s(k,j,i-1) )  &
                        + ( 1.0_wp - fwxm(j,i) ) * wall_s_flux(2)              &
                                                  ) * ddx2                     &
                                                + ( fwyp(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j+1,i) ) * ( s(k,j+1,i)-s(k,j,i) )  &
                        + ( 1.0_wp - fwyp(j,i) ) * wall_s_flux(3)              &
                                                   -fwym(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j-1,i) ) * ( s(k,j,i)-s(k,j-1,i) )  &
                        + ( 1.0_wp - fwym(j,i) ) * wall_s_flux(4)              &
                                                  ) * ddy2
                ENDDO
             ENDIF

!
!--          Compute vertical diffusion. In case that surface fluxes have been
!--          prescribed or computed at bottom and/or top, index k starts/ends at
!--          nzb+2 or nzt-1, respectively.
             DO  k = nzb_diff_s_inner(j,i), nzt_diff

                tend(k,j,i) = tend(k,j,i)                                      &
                                       + 0.5_wp * (                            &
            ( kh(k,j,i) + kh(k+1,j,i) ) * ( s(k+1,j,i)-s(k,j,i) ) * ddzu(k+1)  &
                                                            * rho_air_zw(k)    &
          - ( kh(k,j,i) + kh(k-1,j,i) ) * ( s(k,j,i)-s(k-1,j,i) ) * ddzu(k)    &
                                                            * rho_air_zw(k-1)  &
                                                  ) * ddzw(k) * drho_air(k)
             ENDDO

!
!--          Vertical diffusion at the first computational gridpoint along
!--          z-direction
             IF ( use_surface_fluxes )  THEN

                k = nzb_s_inner(j,i)+1

                tend(k,j,i) = tend(k,j,i)                                      &
                                       + ( 0.5_wp * ( kh(k,j,i)+kh(k+1,j,i) )  &
                                                  * ( s(k+1,j,i)-s(k,j,i) )    &
                                                  * ddzu(k+1)                  &
                                                  * rho_air_zw(k)              &
                                           + s_flux_b(j,i)                     &
                                         ) * ddzw(k) * drho_air(k)

             ENDIF

!
!--          Vertical diffusion at the last computational gridpoint along
!--          z-direction
             IF ( use_top_fluxes )  THEN

                k = nzt

                tend(k,j,i) = tend(k,j,i)                                      &
                                       + ( - s_flux_t(j,i)                     &
                                           - 0.5_wp * ( kh(k-1,j,i)+kh(k,j,i) )&
                                                    * ( s(k,j,i)-s(k-1,j,i) )  &
                                                    * ddzu(k)                  &
                                                    * rho_air_zw(k-1)          &
                                         ) * ddzw(k) * drho_air(k)

             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_s


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points - accelerator version
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_s_acc( s, s_flux_b, s_flux_t, wall_s_flux )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, kh, tend, drho_air, rho_air_zw
           
       USE control_parameters,                                                 & 
           ONLY: use_surface_fluxes, use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx2, ddy2, fwxm, fwxp, fwym, fwyp, wall_w_x, wall_w_y
       
       USE indices, &
           ONLY: i_left, i_right, j_north, j_south, nxlg, nxrg, nyng, nysg,    &
                 nzb, nzb_diff_s_inner, nzb_s_inner, nzb_s_outer, nzt, nzt_diff
           
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< 
       INTEGER(iwp) ::  j                 !< 
       INTEGER(iwp) ::  k                 !< 
       REAL(wp)     ::  wall_s_flux(0:4)  !< 
       REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  s_flux_b !< 
       REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  s_flux_t !< 
#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s  !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  s  !<
#endif

       !$acc kernels present( ddzu, ddzw, fwxm, fwxp, fwym, fwyp, kh )        &
       !$acc         present( nzb_diff_s_inner, nzb_s_inner, nzb_s_outer, s ) &
       !$acc         present( s_flux_b, s_flux_t, tend, wall_s_flux )         &
       !$acc         present( wall_w_x, wall_w_y )
       DO  i = i_left, i_right
          DO  j = j_south, j_north
!
!--          Compute horizontal diffusion
             DO  k = 1, nzt
                IF ( k > nzb_s_outer(j,i) )  THEN

                   tend(k,j,i) = tend(k,j,i)                                   &
                                          + 0.5_wp * (                         &
                        ( kh(k,j,i) + kh(k,j,i+1) ) * ( s(k,j,i+1)-s(k,j,i) )  &
                      - ( kh(k,j,i) + kh(k,j,i-1) ) * ( s(k,j,i)-s(k,j,i-1) )  &
                                                     ) * ddx2                  &
                                          + 0.5_wp * (                         &
                        ( kh(k,j,i) + kh(k,j+1,i) ) * ( s(k,j+1,i)-s(k,j,i) )  &
                      - ( kh(k,j,i) + kh(k,j-1,i) ) * ( s(k,j,i)-s(k,j-1,i) )  &
                                                     ) * ddy2
                ENDIF
             ENDDO

!
!--          Apply prescribed horizontal wall heatflux where necessary
             DO  k = 1, nzt
                IF ( k > nzb_s_inner(j,i)  .AND.  k <= nzb_s_outer(j,i)  .AND. &
                     ( wall_w_x(j,i) /= 0.0_wp  .OR.  wall_w_y(j,i) /= 0.0_wp ) )    &
                THEN
                   tend(k,j,i) = tend(k,j,i)                                   &
                                                + ( fwxp(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j,i+1) ) * ( s(k,j,i+1)-s(k,j,i) )  &
                        + ( 1.0_wp - fwxp(j,i) ) * wall_s_flux(1)              &
                                                   -fwxm(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j,i-1) ) * ( s(k,j,i)-s(k,j,i-1) )  &
                        + ( 1.0_wp - fwxm(j,i) ) * wall_s_flux(2)              &
                                                  ) * ddx2                     &
                                                + ( fwyp(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j+1,i) ) * ( s(k,j+1,i)-s(k,j,i) )  &
                        + ( 1.0_wp - fwyp(j,i) ) * wall_s_flux(3)              &
                                                   -fwym(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j-1,i) ) * ( s(k,j,i)-s(k,j-1,i) )  &
                        + ( 1.0_wp - fwym(j,i) ) * wall_s_flux(4)              &
                                                  ) * ddy2
                ENDIF
             ENDDO

!
!--          Compute vertical diffusion. In case that surface fluxes have been
!--          prescribed or computed at bottom and/or top, index k starts/ends at
!--          nzb+2 or nzt-1, respectively.
             DO  k = 1, nzt_diff
                IF ( k >= nzb_diff_s_inner(j,i) )  THEN
                   tend(k,j,i) = tend(k,j,i)                                   &
                                       + 0.5_wp * (                            &
            ( kh(k,j,i) + kh(k+1,j,i) ) * ( s(k+1,j,i)-s(k,j,i) ) * ddzu(k+1)  &
                                                            * rho_air_zw(k)    &
          - ( kh(k,j,i) + kh(k-1,j,i) ) * ( s(k,j,i)-s(k-1,j,i) ) * ddzu(k)    &
                                                            * rho_air_zw(k-1)  &
                                                  ) * ddzw(k) * drho_air(k)
                ENDIF
             ENDDO

!
!--          Vertical diffusion at the first computational gridpoint along
!--          z-direction
             DO  k = 1, nzt
                IF ( use_surface_fluxes  .AND.  k == nzb_s_inner(j,i)+1 )  THEN
                   tend(k,j,i) = tend(k,j,i)                                   &
                                          + ( 0.5_wp * ( kh(k,j,i)+kh(k+1,j,i) )&
                                                     * ( s(k+1,j,i)-s(k,j,i) ) &
                                                     * ddzu(k+1)               &
                                                     * rho_air_zw(k)           &
                                              + s_flux_b(j,i)                  &
                                            ) * ddzw(k) * drho_air(k)
                ENDIF

!
!--             Vertical diffusion at the last computational gridpoint along
!--             z-direction
                IF ( use_top_fluxes  .AND.  k == nzt )  THEN
                   tend(k,j,i) = tend(k,j,i)                                   &
                                          + ( - s_flux_t(j,i)                  &
                                              - 0.5_wp * ( kh(k-1,j,i)+kh(k,j,i) )&
                                                       * ( s(k,j,i)-s(k-1,j,i) )  &
                                                       * ddzu(k)                  &
                                                       * rho_air_zw(k-1)          &
                                            ) * ddzw(k) * drho_air(k)
                ENDIF
             ENDDO

          ENDDO
       ENDDO
       !$acc end kernels

    END SUBROUTINE diffusion_s_acc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_s_ij( i, j, s, s_flux_b, s_flux_t, wall_s_flux )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, kh, tend, drho_air, rho_air_zw
           
       USE control_parameters,                                                 & 
           ONLY: use_surface_fluxes, use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx2, ddy2, fwxm, fwxp, fwym, fwyp, wall_w_x, wall_w_y
       
       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzb_diff_s_inner, nzb_s_inner,  &
                  nzb_s_outer, nzt, nzt_diff
       
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< 
       INTEGER(iwp) ::  j                 !< 
       INTEGER(iwp) ::  k                 !< 
       REAL(wp)     ::  wall_s_flux(0:4)  !< 
       REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  s_flux_b  !< 
       REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  s_flux_t  !< 
#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s !< 
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  s  !< 
#endif

!
!--    Compute horizontal diffusion
       DO  k = nzb_s_outer(j,i)+1, nzt

          tend(k,j,i) = tend(k,j,i)                                            &
                                          + 0.5_wp * (                         &
                        ( kh(k,j,i) + kh(k,j,i+1) ) * ( s(k,j,i+1)-s(k,j,i) )  &
                      - ( kh(k,j,i) + kh(k,j,i-1) ) * ( s(k,j,i)-s(k,j,i-1) )  &
                                                     ) * ddx2                  &
                                          + 0.5_wp * (                         &
                        ( kh(k,j,i) + kh(k,j+1,i) ) * ( s(k,j+1,i)-s(k,j,i) )  &
                      - ( kh(k,j,i) + kh(k,j-1,i) ) * ( s(k,j,i)-s(k,j-1,i) )  &
                                                     ) * ddy2
       ENDDO

!
!--    Apply prescribed horizontal wall heatflux where necessary
       IF ( ( wall_w_x(j,i) /= 0.0_wp ) .OR. ( wall_w_y(j,i) /= 0.0_wp ) )     &
       THEN
          DO  k = nzb_s_inner(j,i)+1, nzb_s_outer(j,i)

             tend(k,j,i) = tend(k,j,i)                                         &
                                                + ( fwxp(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j,i+1) ) * ( s(k,j,i+1)-s(k,j,i) )  &
                        + ( 1.0_wp - fwxp(j,i) ) * wall_s_flux(1)              &
                                                   -fwxm(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j,i-1) ) * ( s(k,j,i)-s(k,j,i-1) )  &
                        + ( 1.0_wp - fwxm(j,i) ) * wall_s_flux(2)              &
                                                  ) * ddx2                     &
                                                + ( fwyp(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j+1,i) ) * ( s(k,j+1,i)-s(k,j,i) )  &
                        + ( 1.0_wp - fwyp(j,i) ) * wall_s_flux(3)              &
                                                   -fwym(j,i) * 0.5_wp *       &
                        ( kh(k,j,i) + kh(k,j-1,i) ) * ( s(k,j,i)-s(k,j-1,i) )  &
                        + ( 1.0_wp - fwym(j,i) ) * wall_s_flux(4)              &
                                                  ) * ddy2
          ENDDO
       ENDIF

!
!--    Compute vertical diffusion. In case that surface fluxes have been
!--    prescribed or computed at bottom and/or top, index k starts/ends at
!--    nzb+2 or nzt-1, respectively.
       DO  k = nzb_diff_s_inner(j,i), nzt_diff

          tend(k,j,i) = tend(k,j,i)                                            &
                                       + 0.5_wp * (                            &
            ( kh(k,j,i) + kh(k+1,j,i) ) * ( s(k+1,j,i)-s(k,j,i) ) * ddzu(k+1)  &
                                                            * rho_air_zw(k)    &
          - ( kh(k,j,i) + kh(k-1,j,i) ) * ( s(k,j,i)-s(k-1,j,i) ) * ddzu(k)    &
                                                            * rho_air_zw(k-1)  &
                                                  ) * ddzw(k) * drho_air(k)
       ENDDO

!
!--    Vertical diffusion at the first computational gridpoint along z-direction
       IF ( use_surface_fluxes )  THEN

          k = nzb_s_inner(j,i)+1

          tend(k,j,i) = tend(k,j,i) + ( 0.5_wp * ( kh(k,j,i)+kh(k+1,j,i) )     &
                                               * ( s(k+1,j,i)-s(k,j,i) )       &
                                               * ddzu(k+1)                     &
                                               * rho_air_zw(k)                 &
                                        + s_flux_b(j,i)                        &
                                      ) * ddzw(k) * drho_air(k)

       ENDIF

!
!--    Vertical diffusion at the last computational gridpoint along z-direction
       IF ( use_top_fluxes )  THEN

          k = nzt

          tend(k,j,i) = tend(k,j,i) + ( - s_flux_t(j,i)                        &
                                      - 0.5_wp * ( kh(k-1,j,i)+kh(k,j,i) )     &
                                               * ( s(k,j,i)-s(k-1,j,i) )       &
                                               * ddzu(k)                       &
                                               * rho_air_zw(k-1)               &
                                      ) * ddzw(k) * drho_air(k)

       ENDIF

    END SUBROUTINE diffusion_s_ij

 END MODULE diffusion_s_mod
