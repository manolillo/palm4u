!> @file nudging_mod.f90
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
! $Id: nudging_mod.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 
! 1757 2016-02-22 15:49:32Z maronga
! Bugfix: allow for using higher vertical resolution in nudging file than grid
! spacing in the LES model
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1398 2014-05-07 11:15:00Z heinze
! Subroutine nudge_ref is extended to set u_init and v_init to the current 
! nudging profiles
! 
! 1382 2014-04-30 12:15:41Z boeske
! Changed the weighting factor that is used in the summation of nudging 
! tendencies for profile data output from weight_pres to weight_substep,
! added Neumann boundary conditions for profile data output of nudging terms at
! nzt+1
! 
! 1380 2014-04-28 12:40:45Z heinze
! Subroutine nudge_ref added to account for proper upper scalar boundary 
! conditions in case of nudging
! 
! 1365 2014-04-22 15:03:56Z boeske
! Variable t renamed nt, variable currtnudge renamed tmp_tnudge,
! summation of nudging tendencies for data output added
! +sums_ls_l, tmp_tend
! Added new subroutine calc_tnudge, which calculates the current nudging time 
! scale at each time step
! 
! 1355 2014-04-10 10:21:29Z heinze
! Error message specified.
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
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
! 1268 2013-12-12 09:47:53Z heinze
! bugfix: argument of calc_mean_profile corrected
!
! 1251 2013-11-07 08:14:30Z heinze
! bugfix: calculate dtm and dtp also in vector version
!
! 1249 2013-11-06 10:45:47Z heinze
! remove call of user module
! reformatting
!
! 1241 2013-10-30 11:36:58Z heinze
! Initial revision
!
! Description:
! ------------
!> Nudges u, v, pt and q to given profiles on a relaxation timescale tnudge. 
!> Profiles are read in from NUDGIN_DATA. Code is based on Neggers et al. (2012) 
!> and also part of DALES and UCLA-LES.
!--------------------------------------------------------------------------------!
 MODULE nudge_mod
 

    PRIVATE
    PUBLIC init_nudge, calc_tnudge, nudge, nudge_ref
    SAVE

    INTERFACE nudge
       MODULE PROCEDURE nudge
       MODULE PROCEDURE nudge_ij
    END INTERFACE nudge

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE init_nudge

       USE arrays_3d,                                                          &
           ONLY:  ptnudge, qnudge, timenudge, tmp_tnudge, tnudge, unudge,      &
                  vnudge, wnudge, zu

       USE control_parameters,                                                 &
           ONLY:  dt_3d, lptnudge, lqnudge, lunudge, lvnudge, lwnudge,         &
                   message_string, ntnudge

       USE indices,                                                            &
           ONLY:  nzb, nzt

       USE kinds

       IMPLICIT NONE


       INTEGER(iwp) ::  finput = 90  !<
       INTEGER(iwp) ::  ierrn        !<
       INTEGER(iwp) ::  k            !<
       INTEGER(iwp) ::  nt            !<

       CHARACTER(1) ::  hash     !<

       REAL(wp) ::  highheight   !<
       REAL(wp) ::  highqnudge   !<
       REAL(wp) ::  highptnudge  !<
       REAL(wp) ::  highunudge   !<
       REAL(wp) ::  highvnudge   !<
       REAL(wp) ::  highwnudge   !<
       REAL(wp) ::  hightnudge   !<

       REAL(wp) ::  lowheight    !<
       REAL(wp) ::  lowqnudge    !<
       REAL(wp) ::  lowptnudge   !<
       REAL(wp) ::  lowunudge    !<
       REAL(wp) ::  lowvnudge    !<
       REAL(wp) ::  lowwnudge    !<
       REAL(wp) ::  lowtnudge    !<

       REAL(wp) ::  fac          !<

       ALLOCATE( ptnudge(nzb:nzt+1,1:ntnudge), qnudge(nzb:nzt+1,1:ntnudge), &
                 tnudge(nzb:nzt+1,1:ntnudge), unudge(nzb:nzt+1,1:ntnudge),  &
                 vnudge(nzb:nzt+1,1:ntnudge), wnudge(nzb:nzt+1,1:ntnudge)  )

       ALLOCATE( tmp_tnudge(nzb:nzt) )

       ALLOCATE( timenudge(0:ntnudge) )

       ptnudge = 0.0_wp; qnudge = 0.0_wp; tnudge = 0.0_wp; unudge = 0.0_wp
       vnudge = 0.0_wp; wnudge = 0.0_wp; timenudge = 0.0_wp
!
!--    Initialize array tmp_nudge with a current nudging time scale of 6 hours
       tmp_tnudge = 21600.0_wp

       nt = 0
       OPEN ( finput, FILE='NUDGING_DATA', STATUS='OLD', &
              FORM='FORMATTED', IOSTAT=ierrn )

       IF ( ierrn /= 0 )  THEN
          message_string = 'file NUDGING_DATA does not exist'
          CALL message( 'nudging', 'PA0365', 1, 2, 0, 6, 0 )
       ENDIF

       ierrn = 0

 rloop:DO
          nt = nt + 1
          hash = "#"
          ierrn = 1 ! not zero
!
!--       Search for the next line consisting of "# time", 
!--       from there onwards the profiles will be read
          DO WHILE ( .NOT. ( hash == "#" .AND. ierrn == 0 ) ) 
          
            READ ( finput, *, IOSTAT=ierrn ) hash, timenudge(nt)
            IF ( ierrn < 0 )  EXIT rloop

          ENDDO

          ierrn = 0
          READ ( finput, *, IOSTAT=ierrn ) lowheight, lowtnudge, lowunudge,   &
                                           lowvnudge, lowwnudge , lowptnudge, &
                                           lowqnudge

          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file NUDGING_DATA'
             CALL message( 'nudging', 'PA0366', 1, 2, 0, 6, 0 )
          ENDIF

          ierrn = 0
          READ ( finput, *, IOSTAT=ierrn ) highheight, hightnudge, highunudge,   &
                                           highvnudge, highwnudge , highptnudge, &
                                           highqnudge

          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file NUDGING_DATA'
             CALL message( 'nudging', 'PA0366', 1, 2, 0, 6, 0 )
          ENDIF

          DO  k = nzb, nzt+1
             DO WHILE ( highheight < zu(k) )
                lowheight  = highheight
                lowtnudge  = hightnudge
                lowunudge  = highunudge
                lowvnudge  = highvnudge
                lowwnudge  = highwnudge
                lowptnudge = highptnudge
                lowqnudge  = highqnudge
 
                ierrn = 0
                READ ( finput, *, IOSTAT=ierrn )  highheight , hightnudge , &
                                                  highunudge , highvnudge , &
                                                  highwnudge , highptnudge, &
                                                  highqnudge
                IF (ierrn /= 0 )  THEN
                   WRITE( message_string, * ) 'zu(nzt+1) = ', zu(nzt+1), 'm is ',&
                        'higher than the maximum height in NUDING_DATA which ',  &
                        'is ', lowheight, 'm. Interpolation on PALM ',           &
                        'grid is not possible.'
                   CALL message( 'nudging', 'PA0364', 1, 2, 0, 6, 0 )
                ENDIF
             ENDDO

!
!--          Interpolation of prescribed profiles in space 

             fac = ( highheight - zu(k) ) / ( highheight - lowheight )

             tnudge(k,nt)  = fac * lowtnudge  + ( 1.0_wp - fac ) * hightnudge
             unudge(k,nt)  = fac * lowunudge  + ( 1.0_wp - fac ) * highunudge
             vnudge(k,nt)  = fac * lowvnudge  + ( 1.0_wp - fac ) * highvnudge
             wnudge(k,nt)  = fac * lowwnudge  + ( 1.0_wp - fac ) * highwnudge
             ptnudge(k,nt) = fac * lowptnudge + ( 1.0_wp - fac ) * highptnudge
             qnudge(k,nt)  = fac * lowqnudge  + ( 1.0_wp - fac ) * highqnudge
          ENDDO

       ENDDO rloop

       CLOSE ( finput )

!
!--    Prevent nudging if nudging profiles exhibt too small values
!--    not used so far
       lptnudge  = ANY( ABS( ptnudge ) > 1.0e-8_wp )
       lqnudge   = ANY( ABS( qnudge )  > 1.0e-8_wp )
       lunudge   = ANY( ABS( unudge )  > 1.0e-8_wp )
       lvnudge   = ANY( ABS( vnudge )  > 1.0e-8_wp )
       lwnudge   = ANY( ABS( wnudge )  > 1.0e-8_wp )

    END SUBROUTINE init_nudge


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE calc_tnudge ( time )

       USE arrays_3d,                                                          &
           ONLY:  timenudge, tmp_tnudge, tnudge

       USE control_parameters,                                                 &
           ONLY:  dt_3d 

       USE indices,                                                            &
           ONLY:  nzb, nzt

       USE kinds

       IMPLICIT NONE


       REAL(wp) ::  dtm         !<
       REAL(wp) ::  dtp         !<
       REAL(wp) ::  time        !<

       INTEGER(iwp) ::  k   !<
       INTEGER(iwp) ::  nt  !<

       nt = 1
       DO WHILE ( time > timenudge(nt) )
         nt = nt+1
       ENDDO
       IF ( time /= timenudge(1) ) THEN
         nt = nt-1
       ENDIF

       dtm = ( time - timenudge(nt) ) / ( timenudge(nt+1) - timenudge(nt) )
       dtp = ( timenudge(nt+1) - time ) / ( timenudge(nt+1) - timenudge(nt) )

       DO  k = nzb, nzt
          tmp_tnudge(k) = MAX( dt_3d, tnudge(k,nt) * dtp + tnudge(k,nt+1) * dtm )
       ENDDO

    END SUBROUTINE calc_tnudge

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE nudge ( time, prog_var )

       USE arrays_3d,                                                          &
           ONLY:  pt, ptnudge, q, qnudge, tend, timenudge, tmp_tnudge, tnudge, &
                  u, unudge, v, vnudge

       USE control_parameters,                                                 &
           ONLY:  dt_3d, intermediate_timestep_count, message_string

       USE indices,                                                            &
           ONLY:  nxl, nxr, nys, nyn, nzb, nzb_u_inner, nzt

       USE kinds

       USE statistics,                                                         &
           ONLY:  hom, sums_ls_l, weight_substep

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  prog_var  !<

       REAL(wp) ::  tmp_tend    !<
       REAL(wp) ::  dtm         !<
       REAL(wp) ::  dtp         !<
       REAL(wp) ::  time        !<

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<
       INTEGER(iwp) ::  nt  !<


       nt = 1
       DO WHILE ( time > timenudge(nt) )
         nt = nt+1
       ENDDO
       IF ( time /= timenudge(1) ) THEN
         nt = nt-1
       ENDIF

       dtm = ( time - timenudge(nt) ) / ( timenudge(nt+1) - timenudge(nt) )
       dtp = ( timenudge(nt+1) - time ) / ( timenudge(nt+1) - timenudge(nt) )

       SELECT CASE ( prog_var )

          CASE ( 'u' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb_u_inner(j,i)+1, nzt

                      tmp_tend = - ( hom(k,1,1,0) - ( unudge(k,nt) * dtp +     &
                                     unudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend

                      sums_ls_l(k,6) = sums_ls_l(k,6) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO
                  
                   sums_ls_l(nzt+1,6) = sums_ls_l(nzt,6)
 
                ENDDO
            ENDDO

          CASE ( 'v' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb_u_inner(j,i)+1, nzt

                      tmp_tend = - ( hom(k,1,2,0) - ( vnudge(k,nt) * dtp +     &
                                     vnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend

                      sums_ls_l(k,7) = sums_ls_l(k,7) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO
                  
                   sums_ls_l(nzt+1,7) = sums_ls_l(nzt,7)

                ENDDO
            ENDDO

          CASE ( 'pt' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb_u_inner(j,i)+1, nzt

                      tmp_tend = - ( hom(k,1,4,0) - ( ptnudge(k,nt) * dtp +    &
                                     ptnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend

                      sums_ls_l(k,4) = sums_ls_l(k,4) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO

                   sums_ls_l(nzt+1,4) = sums_ls_l(nzt,4)

                ENDDO
            ENDDO

          CASE ( 'q' )

             DO  i = nxl, nxr
                DO  j = nys, nyn

                   DO  k = nzb_u_inner(j,i)+1, nzt

                      tmp_tend = - ( hom(k,1,41,0) - ( qnudge(k,nt) * dtp +    &
                                     qnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                      tend(k,j,i) = tend(k,j,i) + tmp_tend

                      sums_ls_l(k,5) = sums_ls_l(k,5) + tmp_tend *             &
                                     weight_substep(intermediate_timestep_count)
                   ENDDO
                  
                   sums_ls_l(nzt+1,5) = sums_ls_l(nzt,5)

                ENDDO
            ENDDO

          CASE DEFAULT
             message_string = 'unknown prognostic variable "' // prog_var // '"'
             CALL message( 'nudge', 'PA0367', 1, 2, 0, 6, 0 )

       END SELECT

    END SUBROUTINE nudge


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!

    SUBROUTINE nudge_ij( i, j, time, prog_var )

       USE arrays_3d,                                                          &
           ONLY:  pt, ptnudge, q, qnudge, tend, timenudge, tmp_tnudge, tnudge, &
                  u, unudge, v, vnudge

       USE control_parameters,                                                 &
           ONLY:  dt_3d, intermediate_timestep_count, message_string

       USE indices,                                                            &
           ONLY:  nxl, nxr, nys, nyn, nzb, nzb_u_inner, nzt

       USE kinds

       USE statistics,                                                         &
           ONLY:  hom, sums_ls_l, weight_substep

       IMPLICIT NONE


       CHARACTER (LEN=*) ::  prog_var  !<

       REAL(wp) ::  tmp_tend    !<
       REAL(wp) ::  dtm         !<
       REAL(wp) ::  dtp         !<
       REAL(wp) ::  time        !<

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<
       INTEGER(iwp) ::  nt  !<


       nt = 1
       DO WHILE ( time > timenudge(nt) )
         nt = nt+1
       ENDDO
       IF ( time /= timenudge(1) )  THEN
         nt = nt-1
       ENDIF

       dtm = ( time - timenudge(nt) ) / ( timenudge(nt+1) - timenudge(nt) )
       dtp = ( timenudge(nt+1) - time ) / ( timenudge(nt+1) - timenudge(nt) )

       SELECT CASE ( prog_var )

          CASE ( 'u' )

             DO  k = nzb_u_inner(j,i)+1, nzt

                tmp_tend = - ( hom(k,1,1,0) - ( unudge(k,nt) * dtp +           &
                               unudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend

                sums_ls_l(k,6) = sums_ls_l(k,6) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,6) = sums_ls_l(nzt,6)

          CASE ( 'v' )

             DO  k = nzb_u_inner(j,i)+1, nzt

                tmp_tend = - ( hom(k,1,2,0) - ( vnudge(k,nt) * dtp +           &
                               vnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend

                sums_ls_l(k,7) = sums_ls_l(k,7) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,7) = sums_ls_l(nzt,7)

          CASE ( 'pt' )

             DO  k = nzb_u_inner(j,i)+1, nzt

                tmp_tend = - ( hom(k,1,4,0) - ( ptnudge(k,nt) * dtp +          &
                               ptnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend

                sums_ls_l(k,4) = sums_ls_l(k,4) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,4) = sums_ls_l(nzt,4)


          CASE ( 'q' )

             DO  k = nzb_u_inner(j,i)+1, nzt

                tmp_tend = - ( hom(k,1,41,0) - ( qnudge(k,nt) * dtp +          &
                               qnudge(k,nt+1) * dtm ) ) / tmp_tnudge(k)

                tend(k,j,i) = tend(k,j,i) + tmp_tend

                sums_ls_l(k,5) = sums_ls_l(k,5) + tmp_tend                     &
                                 * weight_substep(intermediate_timestep_count)
             ENDDO

             sums_ls_l(nzt+1,5) = sums_ls_l(nzt,5)

          CASE DEFAULT
             message_string = 'unknown prognostic variable "' // prog_var // '"'
             CALL message( 'nudge', 'PA0367', 1, 2, 0, 6, 0 )

       END SELECT


    END SUBROUTINE nudge_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE nudge_ref ( time )

       USE arrays_3d,                                                          &
           ONLY:  time_vert, ptnudge, pt_init, qnudge, q_init, unudge, u_init, &
                  vnudge, v_init

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  nt                    !<

       REAL(wp)             ::  fac           !<
       REAL(wp), INTENT(in) ::  time          !<

!
!--    Interpolation in time of NUDGING_DATA for pt_init and q_init. This is 
!--    needed for correct upper boundary conditions for pt and q and in case that 
!      large scale subsidence as well as scalar Rayleigh-damping are used
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
        nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

       pt_init = ptnudge(:,nt) + fac * ( ptnudge(:,nt+1) - ptnudge(:,nt) )
       q_init  = qnudge(:,nt) + fac * ( qnudge(:,nt+1) - qnudge(:,nt) )
       u_init  = unudge(:,nt) + fac * ( unudge(:,nt+1) - unudge(:,nt) )
       v_init  = vnudge(:,nt) + fac * ( vnudge(:,nt+1) - vnudge(:,nt) )

    END SUBROUTINE nudge_ref

 END MODULE nudge_mod
