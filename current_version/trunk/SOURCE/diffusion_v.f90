!> @file diffusion_v.f90
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
! $Id: diffusion_v.f90 2038 2016-10-26 11:16:56Z knoop $
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
! kmxm_x/_y and kmxp_x/_y change to kmxm and kmxp
!
! Revision 1.1  1997/09/12 06:24:01  raasch
! Initial revision
!
!
! Description:
! ------------
!> Diffusion term of the v-component
!------------------------------------------------------------------------------!
 MODULE diffusion_v_mod
 

    USE wall_fluxes_mod

    PRIVATE
    PUBLIC diffusion_v, diffusion_v_acc

    INTERFACE diffusion_v
       MODULE PROCEDURE diffusion_v
       MODULE PROCEDURE diffusion_v_ij
    END INTERFACE diffusion_v

    INTERFACE diffusion_v_acc
       MODULE PROCEDURE diffusion_v_acc
    END INTERFACE diffusion_v_acc

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_v

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, km, tend, u, v, vsws, vswst, w,                  &
                  drho_air, rho_air_zw
       
       USE control_parameters,                                                 &
           ONLY:  constant_top_momentumflux, topography, use_surface_fluxes,   &
                  use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx, ddy, ddy2, fxm, fxp, wall_v
       
       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nysv, nzb, nzb_diff_v, nzb_v_inner,      &
                  nzb_v_outer, nzt, nzt_diff
       
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< 
       INTEGER(iwp) ::  j     !< 
       INTEGER(iwp) ::  k     !< 
       REAL(wp)     ::  kmxm  !< 
       REAL(wp)     ::  kmxp  !< 
       REAL(wp)     ::  kmzm  !< 
       REAL(wp)     ::  kmzp  !< 

       REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  vsus  !< 

!
!--    First calculate horizontal momentum flux v'u' at vertical walls,
!--    if neccessary
       IF ( topography /= 'flat' )  THEN
          CALL wall_fluxes( vsus, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, nzb_v_inner, &
                            nzb_v_outer, wall_v )
       ENDIF

       DO  i = nxl, nxr
          DO  j = nysv, nyn
!
!--          Compute horizontal diffusion
             DO  k = nzb_v_outer(j,i)+1, nzt
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmxp = 0.25_wp * &
                       ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
                kmxm = 0.25_wp * &
                       ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( kmxp * ( v(k,j,i+1) - v(k,j,i)     ) * ddx         &
                      &   + kmxp * ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy         &
                      &   - kmxm * ( v(k,j,i) - v(k,j,i-1) ) * ddx             &
                      &   - kmxm * ( u(k,j,i) - u(k,j-1,i) ) * ddy             &
                      &   ) * ddx                                              &
                      & + 2.0_wp * (                                           &
                      &           km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i) )      &
                      &         - km(k,j-1,i) * ( v(k,j,i) - v(k,j-1,i) )      &
                      &            ) * ddy2
             ENDDO

!
!--          Wall functions at the left and right walls, respectively
             IF ( wall_v(j,i) /= 0.0_wp )  THEN

                DO  k = nzb_v_inner(j,i)+1, nzb_v_outer(j,i)
                   kmxp = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
                   kmxm = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )
                   
                   tend(k,j,i) = tend(k,j,i)                                   &
                                 + 2.0_wp * (                                  &
                                       km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i) ) &
                                     - km(k,j-1,i) * ( v(k,j,i) - v(k,j-1,i) ) &
                                            ) * ddy2                           &
                                 + (   fxp(j,i) * (                            &
                                  kmxp * ( v(k,j,i+1) - v(k,j,i)     ) * ddx   &
                                + kmxp * ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy   &
                                                  )                            &
                                     - fxm(j,i) * (                            &
                                  kmxm * ( v(k,j,i) - v(k,j,i-1) ) * ddx       &
                                + kmxm * ( u(k,j,i) - u(k,j-1,i) ) * ddy       &
                                                  )                            &
                                     + wall_v(j,i) * vsus(k,j,i)               &
                                   ) * ddx
                ENDDO
             ENDIF

!
!--          Compute vertical diffusion. In case of simulating a Prandtl
!--          layer, index k starts at nzb_v_inner+2.
             DO  k = nzb_diff_v(j,i), nzt_diff
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzp = 0.25_wp * &
                       ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )
                kmzm = 0.25_wp * &
                       ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i) ) * ddzu(k+1)     &
                      &            + ( w(k,j,i) - w(k,j-1,i) ) * ddy           &
                      &            ) * rho_air_zw(k)                           &
                      &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy       &
                      &            ) * rho_air_zw(k-1)                         &
                      &   ) * ddzw(k) * drho_air(k)
             ENDDO

!
!--          Vertical diffusion at the first grid point above the surface,
!--          if the momentum flux at the bottom is given by the Prandtl law
!--          or if it is prescribed by the user.
!--          Difference quotient of the momentum flux is not formed over
!--          half of the grid spacing (2.0*ddzw(k)) any more, since the
!--          comparison with other (LES) models showed that the values of
!--          the momentum flux becomes too large in this case.
!--          The term containing w(k-1,..) (see above equation) is removed here
!--          because the vertical velocity is assumed to be zero at the surface.
             IF ( use_surface_fluxes )  THEN
                k = nzb_v_inner(j,i)+1
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzp = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)   &
                      &            + ( w(k,j,i)   - w(k,j-1,i) ) * ddy         &
                      &            ) * rho_air_zw(k)                           &
                      &   - ( -vsws(j,i) )                                     &
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
                       ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( ( -vswst(j,i) )                                    &
                      &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy       &
                      &            ) * rho_air_zw(k-1)                         &
                      &   ) * ddzw(k) * drho_air(k)
             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_v


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points - accelerator version
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_v_acc

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, km, tend, u, v, vsws, vswst, w,                  &
                  drho_air, rho_air_zw
       
       USE control_parameters,                                                 &
           ONLY:  constant_top_momentumflux, topography, use_surface_fluxes,   &
                  use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx, ddy, ddy2, fxm, fxp, wall_v
       
       USE indices,                                                            &
           ONLY:  i_left, i_right, j_north, j_south, nxl, nxr, nyn, nys, nzb,  &
                  nzb_diff_v, nzb_v_inner, nzb_v_outer, nzt, nzt_diff
       
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< 
       INTEGER(iwp) ::  j     !< 
       INTEGER(iwp) ::  k     !< 
       REAL(wp)     ::  kmxm  !< 
       REAL(wp)     ::  kmxp  !< 
       REAL(wp)     ::  kmzm  !< 
       REAL(wp)     ::  kmzp  !< 

       REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr) ::  vsus  !< 
       !$acc declare create ( vsus )

!
!--    First calculate horizontal momentum flux v'u' at vertical walls,
!--    if neccessary
       IF ( topography /= 'flat' )  THEN
          CALL wall_fluxes_acc( vsus, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp,          &
                                nzb_v_inner, nzb_v_outer, wall_v )
       ENDIF

       !$acc kernels present ( u, v, w, km, tend, vsws, vswst )                &
       !$acc         present ( ddzu, ddzw, fxm, fxp, wall_v )                  &
       !$acc         present ( nzb_v_inner, nzb_v_outer, nzb_diff_v )
       DO  i = i_left, i_right
          DO  j = j_south, j_north
!
!--          Compute horizontal diffusion
             DO  k = 1, nzt
                IF ( k > nzb_v_outer(j,i) )  THEN
!
!--                Interpolate eddy diffusivities on staggered gridpoints
                   kmxp = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
                   kmxm = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )

                   tend(k,j,i) = tend(k,j,i)                                   &
                         & + ( kmxp * ( v(k,j,i+1) - v(k,j,i)     ) * ddx      &
                         &   + kmxp * ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy      &
                         &   - kmxm * ( v(k,j,i) - v(k,j,i-1) ) * ddx          &
                         &   - kmxm * ( u(k,j,i) - u(k,j-1,i) ) * ddy          &
                         &   ) * ddx                                           &
                         & + 2.0_wp * (                                        &
                         &           km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i) )   &
                         &         - km(k,j-1,i) * ( v(k,j,i) - v(k,j-1,i) )   &
                         &            ) * ddy2
                ENDIF
             ENDDO

!
!--          Wall functions at the left and right walls, respectively
             DO  k = 1, nzt
                IF( k > nzb_v_inner(j,i)  .AND.  k <= nzb_v_outer(j,i)  .AND.  &
                    wall_v(j,i) /= 0.0_wp )  THEN

                   kmxp = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
                   kmxm = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )
                   
                   tend(k,j,i) = tend(k,j,i)                                   &
                                 + 2.0_wp * (                                  &
                                       km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i) ) &
                                     - km(k,j-1,i) * ( v(k,j,i) - v(k,j-1,i) ) &
                                            ) * ddy2                           &
                                 + (   fxp(j,i) * (                            &
                                  kmxp * ( v(k,j,i+1) - v(k,j,i)     ) * ddx   &
                                + kmxp * ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy   &
                                                  )                            &
                                     - fxm(j,i) * (                            &
                                  kmxm * ( v(k,j,i) - v(k,j,i-1) ) * ddx       &
                                + kmxm * ( u(k,j,i) - u(k,j-1,i) ) * ddy       &
                                                  )                            &
                                     + wall_v(j,i) * vsus(k,j,i)               &
                                   ) * ddx
                ENDIF
             ENDDO

!
!--          Compute vertical diffusion. In case of simulating a Prandtl
!--          layer, index k starts at nzb_v_inner+2.
             DO  k = 1, nzt_diff
                IF ( k >= nzb_diff_v(j,i) )  THEN
!
!--                Interpolate eddy diffusivities on staggered gridpoints
                   kmzp = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )
                   kmzm = 0.25_wp *                                            &
                          ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

                   tend(k,j,i) = tend(k,j,i)                                   &
                         & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)&
                         &            + ( w(k,j,i)   - w(k,j-1,i) ) * ddy      &
                         &            ) * rho_air_zw(k)                        &
                         &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)&
                         &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy    &
                         &            ) * rho_air_zw(k-1)                      &
                         &   ) * ddzw(k) * drho_air(k)
                ENDIF
             ENDDO

          ENDDO
       ENDDO

!
!--    Vertical diffusion at the first grid point above the surface,
!--    if the momentum flux at the bottom is given by the Prandtl law
!--    or if it is prescribed by the user.
!--    Difference quotient of the momentum flux is not formed over
!--    half of the grid spacing (2.0*ddzw(k)) any more, since the
!--    comparison with other (LES) models showed that the values of
!--    the momentum flux becomes too large in this case.
!--    The term containing w(k-1,..) (see above equation) is removed here
!--    because the vertical velocity is assumed to be zero at the surface.
       IF ( use_surface_fluxes )  THEN

          DO  i = i_left, i_right
             DO  j = j_south, j_north
          
                k = nzb_v_inner(j,i)+1
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzp = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)   &
                      &            + ( w(k,j,i)   - w(k,j-1,i) ) * ddy         &
                      &            ) * rho_air_zw(k)                           &
                      &   - ( -vsws(j,i) )                                     &
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
                       ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                      &
                      & + ( ( -vswst(j,i) )                                    &
                      &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy       &
                      &            ) * rho_air_zw(k-1)                         &
                      &   ) * ddzw(k) * drho_air(k)
             ENDDO
          ENDDO

       ENDIF
       !$acc end kernels

    END SUBROUTINE diffusion_v_acc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_v_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, km, tend, u, v, vsws, vswst, w,                  &
                  drho_air, rho_air_zw
       
       USE control_parameters,                                                 &
           ONLY:  constant_top_momentumflux, use_surface_fluxes, use_top_fluxes
       
       USE grid_variables,                                                     &
           ONLY:  ddx, ddy, ddy2, fxm, fxp, wall_v
       
       USE indices,                                                            &
           ONLY:  nzb, nzb_diff_v, nzb_v_inner, nzb_v_outer, nzt, nzt_diff
       
       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< 
       INTEGER(iwp) ::  j     !< 
       INTEGER(iwp) ::  k     !< 
       REAL(wp)     ::  kmxm  !< 
       REAL(wp)     ::  kmxp  !< 
       REAL(wp)     ::  kmzm  !< 
       REAL(wp)     ::  kmzp  !< 

       REAL(wp), DIMENSION(nzb:nzt+1) ::  vsus  !< 

!
!--    Compute horizontal diffusion
       DO  k = nzb_v_outer(j,i)+1, nzt
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmxp = 0.25_wp * ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
          kmxm = 0.25_wp * ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )

          tend(k,j,i) = tend(k,j,i)                                            &
                      & + ( kmxp * ( v(k,j,i+1) - v(k,j,i)     ) * ddx         &
                      &   + kmxp * ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy         &
                      &   - kmxm * ( v(k,j,i) - v(k,j,i-1) ) * ddx             &
                      &   - kmxm * ( u(k,j,i) - u(k,j-1,i) ) * ddy             &
                      &   ) * ddx                                              &
                      & + 2.0_wp * (                                           &
                      &           km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i) )      &
                      &         - km(k,j-1,i) * ( v(k,j,i) - v(k,j-1,i) )      &
                      &            ) * ddy2
       ENDDO

!
!--    Wall functions at the left and right walls, respectively
       IF ( wall_v(j,i) /= 0.0_wp )  THEN

!
!--       Calculate the horizontal momentum flux v'u'
          CALL wall_fluxes( i, j, nzb_v_inner(j,i)+1, nzb_v_outer(j,i),        &
                            vsus, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp )

          DO  k = nzb_v_inner(j,i)+1, nzb_v_outer(j,i)
             kmxp = 0.25_wp *                                                  &
                    ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
             kmxm = 0.25_wp *                                                  &
                    ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )

             tend(k,j,i) = tend(k,j,i)                                         &
                                 + 2.0_wp * (                                  &
                                       km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i) ) &
                                     - km(k,j-1,i) * ( v(k,j,i) - v(k,j-1,i) ) &
                                            ) * ddy2                           &
                                 + (   fxp(j,i) * (                            &
                                  kmxp * ( v(k,j,i+1) - v(k,j,i)     ) * ddx   &
                                + kmxp * ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy   &
                                                  )                            &
                                     - fxm(j,i) * (                            &
                                  kmxm * ( v(k,j,i) - v(k,j,i-1) ) * ddx       &
                                + kmxm * ( u(k,j,i) - u(k,j-1,i) ) * ddy       &
                                                  )                            &
                                     + wall_v(j,i) * vsus(k)                   &
                                   ) * ddx
          ENDDO
       ENDIF

!
!--    Compute vertical diffusion. In case of simulating a Prandtl layer,
!--    index k starts at nzb_v_inner+2.
       DO  k = nzb_diff_v(j,i), nzt_diff
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmzp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )
          kmzm = 0.25_wp * ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

          tend(k,j,i) = tend(k,j,i)                                            &
                      & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i) ) * ddzu(k+1)     &
                      &            + ( w(k,j,i) - w(k,j-1,i) ) * ddy           &
                      &            ) * rho_air_zw(k)                           &
                      &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy       &
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
          k = nzb_v_inner(j,i)+1
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmzp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )

          tend(k,j,i) = tend(k,j,i)                                            &
                      & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)   &
                      &            + ( w(k,j,i)   - w(k,j-1,i) ) * ddy         &
                      &            ) * rho_air_zw(k)                           &
                      &   - ( -vsws(j,i) )                                     &
                      &   ) * ddzw(k) * drho_air(k)
       ENDIF

!
!--    Vertical diffusion at the first gridpoint below the top boundary,
!--    if the momentum flux at the top is prescribed by the user
       IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
          k = nzt
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmzm = 0.25_wp * ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

          tend(k,j,i) = tend(k,j,i)                                            &
                      & + ( ( -vswst(j,i) )                                    &
                      &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)   &
                      &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy       &
                      &            ) * rho_air_zw(k-1)                         &
                      &   ) * ddzw(k) * drho_air(k)
       ENDIF

    END SUBROUTINE diffusion_v_ij

 END MODULE diffusion_v_mod
