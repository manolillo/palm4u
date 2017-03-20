!> @file ls_forcing_mod.f90
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
! $Id: ls_forcing_mod.f90 2038 2016-10-26 11:16:56Z knoop $
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
! 
! 1602 2015-06-22 07:50:56Z heinze 
! PA0370 changed to PA0363
!
! 1382 2014-04-30 12:15:41Z boeske
! Renamed variables which store large scale forcing tendencies
! pt_lsa -> td_lsa_lpt, pt_subs -> td_sub_lpt, 
! q_lsa  -> td_lsa_q,   q_subs  -> td_sub_q,
! high|lowpt_lsa -> high|low_td_lsa_lpt, ...
! 
! 1365 2014-04-22 15:03:56Z boeske
! Usage of large scale forcing for pt and q enabled:
! Added new subroutine ls_advec for horizontal large scale advection and large 
! scale subsidence, 
! error message in init_ls_forcing specified,
! variable t renamed nt
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
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
! 1299 2014-03-06 13:15:21Z heinze
! Ensure a zero large scale vertical velocity at the surface
! Bugfix: typo in case of boundary condition in if-clause
!
! 1276 2014-01-15 13:40:41Z heinze
! Use LSF_DATA also in case of Dirichlet bottom boundary condition for scalars
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
!> Calculates large scale forcings (geostrophic wind and subsidence velocity) as 
!> well as surfaces fluxes dependent on time given in an external file (LSF_DATA).
!> Code is based in parts on DALES and UCLA-LES.
!--------------------------------------------------------------------------------!
 MODULE ls_forcing_mod
 

    PRIVATE
    PUBLIC init_ls_forcing, ls_forcing_surf, ls_forcing_vert, ls_advec
    SAVE

    INTERFACE ls_advec
       MODULE PROCEDURE ls_advec
       MODULE PROCEDURE ls_advec_ij
    END INTERFACE ls_advec

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE init_ls_forcing

       USE arrays_3d,                                                          &
           ONLY:  p_surf, pt_surf, q_surf, qsws_surf, shf_surf, td_lsa_lpt,    &
                  td_lsa_q, td_sub_lpt, td_sub_q, time_surf, time_vert,        &
                  heatflux_input_conversion, waterflux_input_conversion,       &
                  ug_vert, vg_vert, wsubs_vert, zu

       USE control_parameters,                                                 &
           ONLY:  end_time, lsf_surf, lsf_vert, message_string, nlsf

       USE indices,                                                            &
           ONLY:  ngp_sums_ls, nzb, nz, nzt

       USE kinds

       USE statistics,                                                         &
           ONLY:  sums_ls_l


       IMPLICIT NONE

       CHARACTER(100) ::  chmess      !<
       CHARACTER(1)   ::  hash        !<

       INTEGER(iwp) ::  ierrn         !<
       INTEGER(iwp) ::  finput = 90   !<
       INTEGER(iwp) ::  k             !<
       INTEGER(iwp) ::  nt             !<

       REAL(wp) ::  fac               !<
       REAL(wp) ::  highheight        !<
       REAL(wp) ::  highug_vert       !<
       REAL(wp) ::  highvg_vert       !<
       REAL(wp) ::  highwsubs_vert    !<
       REAL(wp) ::  lowheight         !<
       REAL(wp) ::  lowug_vert        !<
       REAL(wp) ::  lowvg_vert        !<
       REAL(wp) ::  lowwsubs_vert     !<
       REAL(wp) ::  high_td_lsa_lpt   !<
       REAL(wp) ::  low_td_lsa_lpt    !<
       REAL(wp) ::  high_td_lsa_q     !<
       REAL(wp) ::  low_td_lsa_q      !<
       REAL(wp) ::  high_td_sub_lpt   !<
       REAL(wp) ::  low_td_sub_lpt    !<
       REAL(wp) ::  high_td_sub_q     !<
       REAL(wp) ::  low_td_sub_q      !<
       REAL(wp) ::  r_dummy           !<

       ALLOCATE( p_surf(0:nlsf), pt_surf(0:nlsf), q_surf(0:nlsf),              &
                 qsws_surf(0:nlsf), shf_surf(0:nlsf),                          &
                 td_lsa_lpt(nzb:nzt+1,0:nlsf), td_lsa_q(nzb:nzt+1,0:nlsf),     &
                 td_sub_lpt(nzb:nzt+1,0:nlsf), td_sub_q(nzb:nzt+1,0:nlsf),     &
                 time_vert(0:nlsf), time_surf(0:nlsf),                         &
                 ug_vert(nzb:nzt+1,0:nlsf), vg_vert(nzb:nzt+1,0:nlsf),         &
                 wsubs_vert(nzb:nzt+1,0:nlsf) )

       p_surf = 0.0_wp; pt_surf = 0.0_wp; q_surf = 0.0_wp; qsws_surf = 0.0_wp
       shf_surf = 0.0_wp; time_vert = 0.0_wp; td_lsa_lpt = 0.0_wp
       td_lsa_q = 0.0_wp; td_sub_lpt = 0.0_wp; td_sub_q = 0.0_wp
       time_surf = 0.0_wp; ug_vert = 0.0_wp; vg_vert = 0.0_wp
       wsubs_vert = 0.0_wp

!
!--    Array for storing large scale forcing and nudging tendencies at each 
!--    timestep for data output
       ALLOCATE( sums_ls_l(nzb:nzt+1,0:7) )
       sums_ls_l = 0.0_wp

       ngp_sums_ls = (nz+2)*6

       OPEN ( finput, FILE='LSF_DATA', STATUS='OLD', &
              FORM='FORMATTED', IOSTAT=ierrn )

       IF ( ierrn /= 0 )  THEN
          message_string = 'file LSF_DATA does not exist'
          CALL message( 'ls_forcing', 'PA0368', 1, 2, 0, 6, 0 )
       ENDIF

       ierrn = 0
!
!--    First three lines of LSF_DATA contain header
       READ ( finput, FMT='(a100)', IOSTAT=ierrn ) chmess
       READ ( finput, FMT='(a100)', IOSTAT=ierrn ) chmess
       READ ( finput, FMT='(a100)', IOSTAT=ierrn ) chmess

       IF ( ierrn /= 0 )  THEN
          message_string = 'errors in file LSF_DATA'
          CALL message( 'ls_forcing', 'PA0369', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Surface values are read in
       nt     = 0
       ierrn = 0

       DO WHILE ( time_surf(nt) < end_time )
          nt = nt + 1
          READ ( finput, *, IOSTAT = ierrn ) time_surf(nt), shf_surf(nt),      &
                                             qsws_surf(nt), pt_surf(nt),       &
                                             q_surf(nt), p_surf(nt)

          IF ( ierrn < 0 )  THEN
            WRITE ( message_string, * ) 'No time dependent surface variables ',&
                              'in&LSF_DATA for end of run found'

             CALL message( 'ls_forcing', 'PA0363', 1, 2, 0, 6, 0 )
          ENDIF
       ENDDO

       shf_surf  = shf_surf  * heatflux_input_conversion(nzb)
       qsws_surf = qsws_surf * waterflux_input_conversion(nzb)

       IF ( time_surf(1) > end_time )  THEN
          WRITE ( message_string, * ) 'No time dependent surface variables in ',&
                                     '&LSF_DATA for end of run found - ',      &
                                     'lsf_surf is set to FALSE'
          CALL message( 'ls_forcing', 'PA0371', 0, 0, 0, 6, 0 )
          lsf_surf = .FALSE.
       ENDIF

!
!--    Go to the end of the list with surface variables
       DO WHILE ( ierrn == 0 )
          READ ( finput, *, IOSTAT = ierrn ) r_dummy
       ENDDO

!
!--    Profiles of ug, vg and w_subs are read in (large scale forcing)

       nt = 0
       DO WHILE ( time_vert(nt) < end_time )
          nt = nt + 1
          hash = "#"
          ierrn = 1 ! not zero
!
!--       Search for the next line consisting of "# time", 
!--       from there onwards the profiles will be read
          DO WHILE ( .NOT. ( hash == "#" .AND. ierrn == 0 ) ) 
             READ ( finput, *, IOSTAT=ierrn ) hash, time_vert(nt)
             IF ( ierrn < 0 )  THEN 
                WRITE( message_string, * ) 'No time dependent vertical profiles',&
                                 ' in&LSF_DATA for end of run found'
                CALL message( 'ls_forcing', 'PA0372', 1, 2, 0, 6, 0 )
             ENDIF
          ENDDO

          IF ( nt == 1 .AND. time_vert(nt) > end_time ) EXIT

          READ ( finput, *, IOSTAT=ierrn ) lowheight, lowug_vert, lowvg_vert,  &
                                           lowwsubs_vert, low_td_lsa_lpt,      &
                                           low_td_lsa_q, low_td_sub_lpt,       &
                                           low_td_sub_q
          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file LSF_DATA'
             CALL message( 'ls_forcing', 'PA0369', 1, 2, 0, 6, 0 )
          ENDIF

          READ ( finput, *, IOSTAT=ierrn ) highheight, highug_vert,            &
                                           highvg_vert, highwsubs_vert,        &
                                           high_td_lsa_lpt, high_td_lsa_q,     &
                                           high_td_sub_lpt, high_td_sub_q
      
          IF ( ierrn /= 0 )  THEN
             message_string = 'errors in file LSF_DATA'
             CALL message( 'ls_forcing', 'PA0369', 1, 2, 0, 6, 0 )
          ENDIF


          DO  k = nzb, nzt+1
             IF ( highheight < zu(k) )  THEN
                lowheight      = highheight
                lowug_vert     = highug_vert
                lowvg_vert     = highvg_vert
                lowwsubs_vert  = highwsubs_vert
                low_td_lsa_lpt = high_td_lsa_lpt
                low_td_lsa_q   = high_td_lsa_q
                low_td_sub_lpt = high_td_sub_lpt
                low_td_sub_q   = high_td_sub_q

                ierrn = 0
                READ ( finput, *, IOSTAT=ierrn ) highheight, highug_vert,      &
                                                 highvg_vert, highwsubs_vert,  &
                                                 high_td_lsa_lpt,              &
                                                 high_td_lsa_q,                &
                                                 high_td_sub_lpt, high_td_sub_q

                IF ( ierrn /= 0 )  THEN
                   WRITE( message_string, * ) 'zu(nzt+1) = ', zu(nzt+1), 'm ', &
                        'is higher than the maximum height in LSF_DATA which ',&
                        'is ', lowheight, 'm. Interpolation on PALM ',         &
                        'grid is not possible.'
                   CALL message( 'ls_forcing', 'PA0395', 1, 2, 0, 6, 0 )
                ENDIF

             ENDIF

!
!--          Interpolation of prescribed profiles in space 
             fac = (highheight-zu(k))/(highheight - lowheight)

             ug_vert(k,nt)    = fac * lowug_vert                               &
                                + ( 1.0_wp - fac ) * highug_vert
             vg_vert(k,nt)    = fac * lowvg_vert                               &
                                + ( 1.0_wp - fac ) * highvg_vert
             wsubs_vert(k,nt) = fac * lowwsubs_vert                            &
                                + ( 1.0_wp - fac ) * highwsubs_vert

             td_lsa_lpt(k,nt) = fac * low_td_lsa_lpt                           &
                                + ( 1.0_wp - fac ) * high_td_lsa_lpt
             td_lsa_q(k,nt)   = fac * low_td_lsa_q                             &
                                + ( 1.0_wp - fac ) * high_td_lsa_q
             td_sub_lpt(k,nt) = fac * low_td_sub_lpt                           &
                                + ( 1.0_wp - fac ) * high_td_sub_lpt
             td_sub_q(k,nt)   = fac * low_td_sub_q                             &
                                + ( 1.0_wp - fac ) * high_td_sub_q

          ENDDO

       ENDDO 

!
!--    Large scale vertical velocity has to be zero at the surface
       wsubs_vert(nzb,:) = 0.0_wp
 
       IF ( time_vert(1) > end_time )  THEN
          WRITE ( message_string, * ) 'Time dependent large scale profile ',   &
                             'forcing from&LSF_DATA sets in after end of ' ,   &
                             'simulation - lsf_vert is set to FALSE'
          CALL message( 'ls_forcing', 'PA0373', 0, 0, 0, 6, 0 )
          lsf_vert = .FALSE.
       ENDIF

       CLOSE( finput )


    END SUBROUTINE init_ls_forcing 


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE ls_forcing_surf ( time )

       USE arrays_3d,                                                          &
           ONLY:  p_surf, pt_surf, q_surf, qsws, qsws_surf, shf, shf_surf,     &
                  time_surf, time_vert, ug, ug_vert, vg, vg_vert

       USE control_parameters,                                                 &
           ONLY:  bc_q_b, ibc_pt_b, ibc_q_b, pt_surface, q_surface,            &
                  surface_pressure

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  nt                     !<

       REAL(wp)             :: fac            !<
       REAL(wp), INTENT(in) :: time           !<

!
!--    Interpolation in time of LSF_DATA at the surface
       nt = 1
       DO WHILE ( time > time_surf(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_surf(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time -time_surf(nt) ) / ( time_surf(nt+1) - time_surf(nt) )

       IF ( ibc_pt_b == 0 )  THEN
!
!--       In case of Dirichlet boundary condition shf must not 
!--       be set - it is calculated via MOST in prandtl_fluxes
          pt_surface = pt_surf(nt) + fac * ( pt_surf(nt+1) - pt_surf(nt) )

       ELSEIF ( ibc_pt_b == 1 )  THEN
!
!--       In case of Neumann boundary condition pt_surface is needed for 
!--       calculation of reference density
          shf        = shf_surf(nt) + fac * ( shf_surf(nt+1) - shf_surf(nt) )
          pt_surface = pt_surf(nt) + fac * ( pt_surf(nt+1) - pt_surf(nt) )

       ENDIF

       IF ( ibc_q_b == 0 )  THEN
!
!--       In case of Dirichlet boundary condition qsws must not 
!--       be set - it is calculated via MOST in prandtl_fluxes
          q_surface = q_surf(nt) + fac * ( q_surf(nt+1) - q_surf(nt) )

       ELSEIF ( ibc_q_b == 1 )  THEN

          qsws = qsws_surf(nt) + fac * ( qsws_surf(nt+1) - qsws_surf(nt) )

       ENDIF

       surface_pressure = p_surf(nt) + fac * ( p_surf(nt+1) - p_surf(nt) )

    END SUBROUTINE ls_forcing_surf 


!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE ls_forcing_vert ( time )

       USE arrays_3d,                                                          &
           ONLY:  time_vert, ug, ug_vert, vg, vg_vert, w_subs, wsubs_vert

       USE control_parameters,                                                 &
           ONLY:  large_scale_subsidence

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  nt                     !<

       REAL(wp)             ::  fac           !<
       REAL(wp), INTENT(in) ::  time          !<

!
!--    Interpolation in time of LSF_DATA for ug, vg and w_subs
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

       ug     = ug_vert(:,nt) + fac * ( ug_vert(:,nt+1) - ug_vert(:,nt) )
       vg     = vg_vert(:,nt) + fac * ( vg_vert(:,nt+1) - vg_vert(:,nt) )

       IF ( large_scale_subsidence )  THEN
          w_subs = wsubs_vert(:,nt)                                            &
                   + fac * ( wsubs_vert(:,nt+1) - wsubs_vert(:,nt) )
       ENDIF

    END SUBROUTINE ls_forcing_vert


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE ls_advec ( time, prog_var )

       USE arrays_3d,                                                          &
           ONLY:  td_lsa_lpt, td_lsa_q, td_sub_lpt, td_sub_q, tend, time_vert        
       
       USE control_parameters,                                                 &
           ONLY:  large_scale_subsidence, use_subsidence_tendencies
       
       USE indices
       
       USE kinds

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  prog_var   !< 

       REAL(wp), INTENT(in)  :: time    !< 
       REAL(wp) :: fac                  !<  

       INTEGER(iwp) ::  i               !< 
       INTEGER(iwp) ::  j               !< 
       INTEGER(iwp) ::  k               !< 
       INTEGER(iwp) ::  nt               !< 

!
!--    Interpolation in time of LSF_DATA 
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

!
!--    Add horizontal large scale advection tendencies of pt and q 
       SELECT CASE ( prog_var )

          CASE ( 'pt' )

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_u_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) + td_lsa_lpt(k,nt) + fac *     &
                                    ( td_lsa_lpt(k,nt+1) - td_lsa_lpt(k,nt) )
                   ENDDO
                ENDDO
             ENDDO

          CASE ( 'q' )

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_u_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) + td_lsa_q(k,nt) + fac *       &
                                    ( td_lsa_q(k,nt+1) - td_lsa_q(k,nt) )
                   ENDDO
                ENDDO
             ENDDO

       END SELECT

!
!--    Subsidence of pt and q with prescribed subsidence tendencies
       IF ( large_scale_subsidence .AND. use_subsidence_tendencies )  THEN

          SELECT CASE ( prog_var )

             CASE ( 'pt' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_u_inner(j,i)+1, nzt
                         tend(k,j,i) = tend(k,j,i) + td_sub_lpt(k,nt) + fac *  &
                                       ( td_sub_lpt(k,nt+1) - td_sub_lpt(k,nt) )
                      ENDDO
                   ENDDO
                ENDDO
  
             CASE ( 'q' )

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_u_inner(j,i)+1, nzt
                         tend(k,j,i) = tend(k,j,i) + td_sub_q(k,nt) + fac *    &
                                       ( td_sub_q(k,nt+1) - td_sub_q(k,nt) )
                      ENDDO
                   ENDDO
                ENDDO

          END SELECT

       ENDIF

    END SUBROUTINE ls_advec


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE ls_advec_ij ( i, j, time, prog_var )

       USE arrays_3d,                                                          &
           ONLY:  td_lsa_lpt, td_lsa_q, td_sub_lpt, td_sub_q, tend, time_vert        
       
       USE control_parameters,                                                 &
           ONLY:  large_scale_subsidence, use_subsidence_tendencies
       
       USE indices
       
       USE kinds

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  prog_var   !< 

       REAL(wp), INTENT(in)  :: time    !< 
       REAL(wp) :: fac                  !< 

       INTEGER(iwp) ::  i               !< 
       INTEGER(iwp) ::  j               !< 
       INTEGER(iwp) ::  k               !< 
       INTEGER(iwp) ::  nt               !< 

!
!--    Interpolation in time of LSF_DATA 
       nt = 1
       DO WHILE ( time > time_vert(nt) )
          nt = nt + 1
       ENDDO
       IF ( time /= time_vert(nt) )  THEN
         nt = nt - 1
       ENDIF

       fac = ( time-time_vert(nt) ) / ( time_vert(nt+1)-time_vert(nt) )

!
!--    Add horizontal large scale advection tendencies of pt and q 
       SELECT CASE ( prog_var )

          CASE ( 'pt' )

             DO  k = nzb_u_inner(j,i)+1, nzt
                tend(k,j,i) = tend(k,j,i) + td_lsa_lpt(k,nt)                   &
                              + fac * ( td_lsa_lpt(k,nt+1) - td_lsa_lpt(k,nt) )
             ENDDO

          CASE ( 'q' )

             DO  k = nzb_u_inner(j,i)+1, nzt
                tend(k,j,i) = tend(k,j,i) + td_lsa_q(k,nt)                     &
                              + fac * ( td_lsa_q(k,nt+1) - td_lsa_q(k,nt) )
             ENDDO

       END SELECT

!
!--    Subsidence of pt and q with prescribed profiles
       IF ( large_scale_subsidence .AND. use_subsidence_tendencies )  THEN

          SELECT CASE ( prog_var )

             CASE ( 'pt' )

                DO  k = nzb_u_inner(j,i)+1, nzt
                   tend(k,j,i) = tend(k,j,i) + td_sub_lpt(k,nt) + fac *        &
                                 ( td_sub_lpt(k,nt+1) - td_sub_lpt(k,nt) )
                ENDDO
  
             CASE ( 'q' )

                DO  k = nzb_u_inner(j,i)+1, nzt
                   tend(k,j,i) = tend(k,j,i) + td_sub_q(k,nt) + fac *          &
                                 ( td_sub_q(k,nt+1) - td_sub_q(k,nt) )
                ENDDO

          END SELECT

       ENDIF

    END SUBROUTINE ls_advec_ij


 END MODULE ls_forcing_mod
