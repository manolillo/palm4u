!> @file diffusion_u.f90
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
! $Id: diffusion_u.f90 2038 2016-10-26 11:16:56Z knoop $
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
! 1740 2016-01-13 08:19:40Z raasch
! unnecessary calculations of kmzm and kmzp in wall bounded parts removed
!
! 1691 2015-10-26 16:17:44Z maronga
! Formatting corrections.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
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
! openacc loop and loop vector clauses removed, declare create moved after
! the FORTRAN declaration statement
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator version (*_acc) added
!
! 1001 2012-09-13 14:08:46Z raasch
! arrays comunicated by module instead of parameter list
!
! 978 2012-08-09 08:28:32Z fricke
! outflow damping layer removed
! kmym_x/_y and kmyp_x/_y change to kmym and kmyp
!
! Revision 1.1  1997/09/12 06:23:51  raasch
! Initial revision
!
!
! Description:
! ------------
!> Diffusion term of the u-component
!> @todo additional damping (needed for non-cyclic bc) causes bad vectorization
!>       and slows down the speed on NEC about 5-10%
!------------------------------------------------------------------------------!
 MODULE diffusion_u_mod
 

    USE wall_fluxes_mod

    PRIVATE
    PUBLIC diffusion_u, diffusion_u_acc

    INTERFACE diffusion_u
       MODULE PROCEDURE diffusion_u
       MODULE PROCEDURE diffusion_u_ij
    END INTERFACE diffusion_u

    INTERFACE diffusion_u_acc
       MODULE PROCEDURE diffusion_u_acc
    END INTERFACE diffusion_u_acc

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_u

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, km, tend, u, usws, uswst, v, w,                  &
                  drho_air, rho_air_zw
       
       USE control_parameters,                                                 &
           ONLY:  constant_top_momentumflux, topography, use_surface_fluxes,   &
                  use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx, ddx2, ddy, fym, fyp, wall_u
       
       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nzb, nzb_diff_u, nzb_u_inner,      &
                  nzb_u_outer, nzt, nzt_diff
       
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< 
       INTEGER(iwp) ::  j     !< 
       INTEGER(iwp) ::  k     !< 
       REAL(wp)     ::  kmym  !<
       REAL(wp)     ::  kmyp  !<
       REAL(wp)     ::  kmzm  !<
       REAL(wp)     ::  kmzp  !<

       REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  usvs  !< 

!
!--    First calculate horizontal momentum flux u'v' at vertical walls,
!--    if neccessary
       IF ( topography /= 'flat' )  THEN
          CALL wall_fluxes( usvs, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, nzb_u_inner, &
                            nzb_u_outer, wall_u )
       ENDIF

       DO  i = nxlu, nxr
          DO  j = nys, nyn
!
!--          Compute horizontal diffusion
             DO  k = nzb_u_outer(j,i)+1, nzt
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmyp = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k,j+1,i)+km(k,j,i-1)+km(k,j+1,i-1) )
                kmym = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k,j-1,i)+km(k,j,i-1)+km(k,j-1,i-1) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + 2.0_wp * (                                           &
                      &           km(k,j,i)   * ( u(k,j,i+1) - u(k,j,i)   )    &
                      &         - km(k,j,i-1) * ( u(k,j,i)   - u(k,j,i-1) )    &
                      &            ) * ddx2                                    &
                      & + ( kmyp * ( u(k,j+1,i) - u(k,j,i)     ) * ddy         &
                      &   + kmyp * ( v(k,j+1,i) - v(k,j+1,i-1) ) * ddx         &
                      &   - kmym * ( u(k,j,i) - u(k,j-1,i) ) * ddy             &
                      &   - kmym * ( v(k,j,i) - v(k,j,i-1) ) * ddx             &
                      &   ) * ddy
             ENDDO

!
!--          Wall functions at the north and south walls, respectively
             IF ( wall_u(j,i) /= 0.0_wp )  THEN

                DO  k = nzb_u_inner(j,i)+1, nzb_u_outer(j,i)
                   kmyp = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j+1,i)+km(k,j,i-1)+km(k,j+1,i-1) )
                   kmym = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j-1,i)+km(k,j,i-1)+km(k,j-1,i-1) )

                   tend(k,j,i) = tend(k,j,i)                                   &
                                 + 2.0_wp * (                                  &
                                       km(k,j,i)   * ( u(k,j,i+1) - u(k,j,i) ) &
                                     - km(k,j,i-1) * ( u(k,j,i) - u(k,j,i-1) ) &
                                            ) * ddx2                           &
                                 + (   fyp(j,i) * (                            &
                                  kmyp * ( u(k,j+1,i) - u(k,j,i)     ) * ddy   &
                                + kmyp * ( v(k,j+1,i) - v(k,j+1,i-1) ) * ddx   &
                                                  )                            &
                                     - fym(j,i) * (                            &
                                  kmym * ( u(k,j,i) - u(k,j-1,i) ) * ddy       &
                                + kmym * ( v(k,j,i) - v(k,j,i-1) ) * ddx       &
                                                  )                            &
                                     + wall_u(j,i) * usvs(k,j,i)               &
                                   ) * ddy
                ENDDO
             ENDIF

!
!--          Compute vertical diffusion. In case of simulating a Prandtl layer,
!--          index k starts at nzb_u_inner+2.
             DO  k = nzb_diff_u(j,i), nzt_diff
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzp = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) )
                kmzm = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k-1,j,i)+km(k,j,i-1)+km(k-1,j,i-1) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( kmzp * ( ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)   &
                      &            + ( w(k,j,i)   - w(k,j,i-1) ) * ddx         &
                      &            ) * rho_air_zw(k)                           &
                      &   - kmzm * ( ( u(k,j,i)   - u(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j,i-1) ) * ddx       &
                      &            ) * rho_air_zw(k-1)                         &
                      &   ) * ddzw(k) * drho_air(k)
             ENDDO

!
!--          Vertical diffusion at the first grid point above the surface,
!--          if the momentum flux at the bottom is given by the Prandtl law or
!--          if it is prescribed by the user.
!--          Difference quotient of the momentum flux is not formed over half
!--          of the grid spacing (2.0*ddzw(k)) any more, since the comparison
!--          with other (LES) models showed that the values of the momentum
!--          flux becomes too large in this case.
!--          The term containing w(k-1,..) (see above equation) is removed here
!--          because the vertical velocity is assumed to be zero at the surface.
             IF ( use_surface_fluxes )  THEN
                k = nzb_u_inner(j,i)+1
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzp = 0.25_wp *                                               &
                      ( km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( kmzp * ( ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)   &
                      &            + ( w(k,j,i)   - w(k,j,i-1) ) * ddx         &
                      &            ) * rho_air_zw(k)                           &
                      &   - ( -usws(j,i) )                                     &
                      &   ) * ddzw(k) * drho_air(k)
             ENDIF

!
!--          Vertical diffusion at the first gridpoint below the top boundary,
!--          if the momentum flux at the top is prescribed by the user
             IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
                k = nzt
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzm = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k-1,j,i)+km(k,j,i-1)+km(k-1,j,i-1) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( ( -uswst(j,i) )                                    &
                      &   - kmzm * ( ( u(k,j,i)   - u(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j,i-1) ) * ddx       &
                      &            ) * rho_air_zw(k-1)                         &
                      &   ) * ddzw(k) * drho_air(k)
             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_u


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points - accelerator version
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_u_acc

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, km, tend, u, usws, uswst, v, w,                  &
                  drho_air, rho_air_zw
       
       USE control_parameters,                                                 &
           ONLY:  constant_top_momentumflux, topography, use_surface_fluxes,   &
                  use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx, ddx2, ddy, fym, fyp, wall_u
       
       USE indices,                                                            &
           ONLY:  i_left, i_right, j_north, j_south, nxl, nxr, nyn, nys, nzb,  &
                  nzb_diff_u, nzb_u_inner, nzb_u_outer, nzt, nzt_diff
       
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< 
       INTEGER(iwp) ::  j     !< 
       INTEGER(iwp) ::  k     !< 
       REAL(wp)     ::  kmym  !<
       REAL(wp)     ::  kmyp  !<
       REAL(wp)     ::  kmzm  !<
       REAL(wp)     ::  kmzp  !<

       REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  usvs  !< 
       !$acc declare create ( usvs )

!
!--    First calculate horizontal momentum flux u'v' at vertical walls,
!--    if neccessary
       IF ( topography /= 'flat' )  THEN
          CALL wall_fluxes_acc( usvs, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,          &
                                nzb_u_inner, nzb_u_outer, wall_u )
       ENDIF

       !$acc kernels present ( u, v, w, km, tend, usws, uswst )                &
       !$acc         present ( ddzu, ddzw, fym, fyp, wall_u )                  &
       !$acc         present ( nzb_u_inner, nzb_u_outer, nzb_diff_u )
       DO  i = i_left, i_right
          DO  j = j_south, j_north
!
!--          Compute horizontal diffusion
             DO  k = 1, nzt
                IF ( k > nzb_u_outer(j,i) )  THEN
!
!--                Interpolate eddy diffusivities on staggered gridpoints
                   kmyp = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j+1,i)+km(k,j,i-1)+km(k,j+1,i-1) )
                   kmym = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j-1,i)+km(k,j,i-1)+km(k,j-1,i-1) )

                   tend(k,j,i) = tend(k,j,i)                                   &
                         & + 2.0_wp * (                                        &
                         &           km(k,j,i)   * ( u(k,j,i+1) - u(k,j,i)   ) &
                         &         - km(k,j,i-1) * ( u(k,j,i)   - u(k,j,i-1) ) &
                         &            ) * ddx2                                 &
                         & + ( kmyp * ( u(k,j+1,i) - u(k,j,i)     ) * ddy      &
                         &   + kmyp * ( v(k,j+1,i) - v(k,j+1,i-1) ) * ddx      &
                         &   - kmym * ( u(k,j,i) - u(k,j-1,i) ) * ddy          &
                         &   - kmym * ( v(k,j,i) - v(k,j,i-1) ) * ddx          &
                         &   ) * ddy
                ENDIF
             ENDDO

!
!--          Wall functions at the north and south walls, respectively
             DO  k = 1, nzt
                IF( k > nzb_u_inner(j,i)  .AND.  k <= nzb_u_outer(j,i)  .AND.  &
                    wall_u(j,i) /= 0.0_wp )  THEN

                   kmyp = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j+1,i)+km(k,j,i-1)+km(k,j+1,i-1) )
                   kmym = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j-1,i)+km(k,j,i-1)+km(k,j-1,i-1) )

                   tend(k,j,i) = tend(k,j,i)                                   &
                                 + 2.0_wp * (                                  &
                                       km(k,j,i)   * ( u(k,j,i+1) - u(k,j,i) ) &
                                     - km(k,j,i-1) * ( u(k,j,i) - u(k,j,i-1) ) &
                                            ) * ddx2                           &
                                 + (   fyp(j,i) * (                            &
                                  kmyp * ( u(k,j+1,i) - u(k,j,i)     ) * ddy   &
                                + kmyp * ( v(k,j+1,i) - v(k,j+1,i-1) ) * ddx   &
                                                  )                            &
                                     - fym(j,i) * (                            &
                                  kmym * ( u(k,j,i) - u(k,j-1,i) ) * ddy       &
                                + kmym * ( v(k,j,i) - v(k,j,i-1) ) * ddx       &
                                                  )                            &
                                     + wall_u(j,i) * usvs(k,j,i)               &
                                   ) * ddy
                ENDIF
             ENDDO

!
!--          Compute vertical diffusion. In case of simulating a Prandtl layer,
!--          index k starts at nzb_u_inner+2.
             DO  k = 1, nzt_diff
                IF ( k >= nzb_diff_u(j,i) )  THEN
!
!--                Interpolate eddy diffusivities on staggered gridpoints
                   kmzp = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) )
                   kmzm = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k-1,j,i)+km(k,j,i-1)+km(k-1,j,i-1) )

                   tend(k,j,i) = tend(k,j,i)                                   &
                         & + ( kmzp * ( ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)&
                         &            + ( w(k,j,i)   - w(k,j,i-1) ) * ddx      &
                         &            ) * rho_air_zw(k)                        &
                         &   - kmzm * ( ( u(k,j,i)   - u(k-1,j,i)   ) * ddzu(k)&
                         &            + ( w(k-1,j,i) - w(k-1,j,i-1) ) * ddx    &
                         &            ) * rho_air_zw(k-1)                      &
                         &   ) * ddzw(k) * drho_air(k)
                ENDIF
             ENDDO

          ENDDO
       ENDDO

!
!--    Vertical diffusion at the first grid point above the surface,
!--    if the momentum flux at the bottom is given by the Prandtl law or
!--    if it is prescribed by the user.
!--    Difference quotient of the momentum flux is not formed over half
!--    of the grid spacing (2.0*ddzw(k)) any more, since the comparison
!--    with other (LES) models showed that the values of the momentum
!--    flux becomes too large in this case.
!--    The term containing w(k-1,..) (see above equation) is removed here
!--    because the vertical velocity is assumed to be zero at the surface.
       IF ( use_surface_fluxes )  THEN

          DO  i = i_left, i_right
             DO  j = j_south, j_north
          
                k = nzb_u_inner(j,i)+1
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzp = 0.25_wp *                                               &
                      ( km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( kmzp * ( ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)   &
                      &            + ( w(k,j,i)   - w(k,j,i-1) ) * ddx         &
                      &            ) * rho_air_zw(k)                           &
                      &   - ( -usws(j,i) )                                     &
                      &   ) * ddzw(k) * drho_air(k)
             ENDDO
          ENDDO

       ENDIF

!
!--    Vertical diffusion at the first gridpoint below the top boundary,
!--    if the momentum flux at the top is prescribed by the user
       IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN

          k = nzt

          DO  i = i_left, i_right
             DO  j = j_south, j_north

!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzm = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k-1,j,i)+km(k,j,i-1)+km(k-1,j,i-1) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( ( -uswst(j,i) )                                    &
                      &   - kmzm * ( ( u(k,j,i)   - u(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j,i-1) ) * ddx       &
                      &            ) * rho_air_zw(k-1)                         &
                      &   ) * ddzw(k) * drho_air(k)
             ENDDO
          ENDDO

       ENDIF
       !$acc end kernels

    END SUBROUTINE diffusion_u_acc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_u_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, km, tend, u, usws, uswst, v, w,                  &
                  drho_air, rho_air_zw
       
       USE control_parameters,                                                 &
           ONLY:  constant_top_momentumflux, use_surface_fluxes, use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx, ddx2, ddy, fym, fyp, wall_u
       
       USE indices,                                                            &
           ONLY:  nzb, nzb_diff_u, nzb_u_inner, nzb_u_outer, nzt, nzt_diff
       
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< 
       INTEGER(iwp) ::  j     !< 
       INTEGER(iwp) ::  k     !< 
       REAL(wp)     ::  kmym  !<
       REAL(wp)     ::  kmyp  !<
       REAL(wp)     ::  kmzm  !<
       REAL(wp)     ::  kmzp  !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  usvs  !< 

!
!--    Compute horizontal diffusion
       DO  k = nzb_u_outer(j,i)+1, nzt
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmyp = 0.25_wp * ( km(k,j,i)+km(k,j+1,i)+km(k,j,i-1)+km(k,j+1,i-1) )
          kmym = 0.25_wp * ( km(k,j,i)+km(k,j-1,i)+km(k,j,i-1)+km(k,j-1,i-1) )

          tend(k,j,i) = tend(k,j,i)                                            &
                      & + 2.0_wp * (                                           &
                      &           km(k,j,i)   * ( u(k,j,i+1) - u(k,j,i)   )    &
                      &         - km(k,j,i-1) * ( u(k,j,i)   - u(k,j,i-1) )    &
                      &            ) * ddx2                                    &
                      & + ( kmyp * ( u(k,j+1,i) - u(k,j,i)     ) * ddy         &
                      &   + kmyp * ( v(k,j+1,i) - v(k,j+1,i-1) ) * ddx         &
                      &   - kmym * ( u(k,j,i) - u(k,j-1,i) ) * ddy             &
                      &   - kmym * ( v(k,j,i) - v(k,j,i-1) ) * ddx             &
                      &   ) * ddy
       ENDDO

!
!--    Wall functions at the north and south walls, respectively
       IF ( wall_u(j,i) /= 0.0_wp )  THEN

!
!--       Calculate the horizontal momentum flux u'v'
          CALL wall_fluxes( i, j, nzb_u_inner(j,i)+1, nzb_u_outer(j,i),        &
                            usvs, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )

          DO  k = nzb_u_inner(j,i)+1, nzb_u_outer(j,i)
             kmyp = 0.25_wp * ( km(k,j,i)+km(k,j+1,i)+km(k,j,i-1)+km(k,j+1,i-1) )
             kmym = 0.25_wp * ( km(k,j,i)+km(k,j-1,i)+km(k,j,i-1)+km(k,j-1,i-1) )

             tend(k,j,i) = tend(k,j,i)                                         &
                                 + 2.0_wp * (                                  &
                                       km(k,j,i)   * ( u(k,j,i+1) - u(k,j,i) ) &
                                     - km(k,j,i-1) * ( u(k,j,i) - u(k,j,i-1) ) &
                                            ) * ddx2                           &
                                 + (   fyp(j,i) * (                            &
                                  kmyp * ( u(k,j+1,i) - u(k,j,i)     ) * ddy   &
                                + kmyp * ( v(k,j+1,i) - v(k,j+1,i-1) ) * ddx   &
                                                  )                            &
                                     - fym(j,i) * (                            &
                                  kmym * ( u(k,j,i) - u(k,j-1,i) ) * ddy       &
                                + kmym * ( v(k,j,i) - v(k,j,i-1) ) * ddx       &
                                                  )                            &
                                     + wall_u(j,i) * usvs(k)                   &
                                   ) * ddy
          ENDDO
       ENDIF

!
!--    Compute vertical diffusion. In case of simulating a Prandtl layer,
!--    index k starts at nzb_u_inner+2.
       DO  k = nzb_diff_u(j,i), nzt_diff
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmzp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) )
          kmzm = 0.25_wp * ( km(k,j,i)+km(k-1,j,i)+km(k,j,i-1)+km(k-1,j,i-1) )

          tend(k,j,i) = tend(k,j,i)                                            &
                      & + ( kmzp * ( ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)   &
                      &            + ( w(k,j,i)   - w(k,j,i-1) ) * ddx         &
                      &            ) * rho_air_zw(k)                           &
                      &   - kmzm * ( ( u(k,j,i)   - u(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j,i-1) ) * ddx       &
                      &            ) * rho_air_zw(k-1)                         &
                      &   ) * ddzw(k) * drho_air(k)
       ENDDO

!
!--    Vertical diffusion at the first grid point above the surface, if the
!--    momentum flux at the bottom is given by the Prandtl law or if it is
!--    prescribed by the user.
!--    Difference quotient of the momentum flux is not formed over half of
!--    the grid spacing (2.0*ddzw(k)) any more, since the comparison with 
!--    other (LES) models showed that the values of the momentum flux becomes
!--    too large in this case.
!--    The term containing w(k-1,..) (see above equation) is removed here
!--    because the vertical velocity is assumed to be zero at the surface.
       IF ( use_surface_fluxes )  THEN
          k = nzb_u_inner(j,i)+1
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmzp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) )

          tend(k,j,i) = tend(k,j,i)                                            &
                      & + ( kmzp * ( ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)   &
                      &            + ( w(k,j,i)   - w(k,j,i-1) ) * ddx         &
                      &            ) * rho_air_zw(k)                           &
                      &   - ( -usws(j,i) )                                     &
                      &   ) * ddzw(k) * drho_air(k)
       ENDIF

!
!--    Vertical diffusion at the first gridpoint below the top boundary,
!--    if the momentum flux at the top is prescribed by the user
       IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
          k = nzt
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmzm = 0.25_wp * ( km(k,j,i)+km(k-1,j,i)+km(k,j,i-1)+km(k-1,j,i-1) )

          tend(k,j,i) = tend(k,j,i)                                            &
                      & + ( ( -uswst(j,i) )                                    &
                      &   - kmzm * ( ( u(k,j,i)   - u(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j,i-1) ) * ddx       &
                      &            ) * rho_air_zw(k-1)                         &
                      &   ) * ddzw(k) * drho_air(k)
       ENDIF

    END SUBROUTINE diffusion_u_ij

 END MODULE diffusion_u_mod
