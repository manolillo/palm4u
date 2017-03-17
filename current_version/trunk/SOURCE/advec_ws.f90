!> @file advec_ws.f90
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
! $Id: advec_ws.f90 2048 2016-11-08 12:18:06Z ketelsen $ 
! 
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1996 2016-08-18 11:42:29Z suehring
! Bugfix concerning calculation of turbulent of turbulent fluxes
! 
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! 
! 1942 2016-06-14 12:18:18Z suehring
! Initialization of flags for ws-scheme moved from init_grid.
! 
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme removed, microphysics_seifert added
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1630 2015-08-26 16:57:23Z suehring
! 
! 
! 1629 2015-08-26 16:56:11Z suehring
! Bugfix concerning wall_flags at left and south PE boundaries
!
! 1581 2015-04-10 13:45:59Z suehring
! 
! 
! 1580 2015-04-10 13:43:49Z suehring
! Bugfix: statistical evaluation of scalar fluxes in case of monotonic limiter
!
! 1567 2015-03-10 17:57:55Z suehring
! Bugfixes in monotonic limiter.
!
! 2015-03-09 13:10:37Z heinze
! Bugfix: REAL constants provided with KIND-attribute in call of 
! intrinsic functions like MAX and MIN 
! 
! 1557 2015-03-05 16:43:04Z suehring
! Enable monotone advection for scalars using monotonic limiter
!
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY list
! 
! 1361 2014-04-16 15:17:48Z hoffmann
! accelerator and vector version for qr and nr added
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute,
! module kinds added
! some formatting adjustments
!
! 1322 2014-03-20 16:38:49Z raasch
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
!
! 1257 2013-11-08 15:18:40Z raasch
! accelerator loop directives removed
!
! 1221 2013-09-10 08:59:13Z raasch
! wall_flags_00 introduced, which holds bits 32-...
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1115 2013-03-26 18:16:16Z hoffmann
! calculation of qr and nr is restricted to precipitation
!
! 1053 2012-11-13 17:11:03Z hoffmann
! necessary expansions according to the two new prognostic equations (nr, qr) 
! of the two-moment cloud physics scheme:
! +flux_l_*, flux_s_*, diss_l_*, diss_s_*, sums_ws*s_ws_l 
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1027 2012-10-15 17:18:39Z suehring
! Bugfix in calculation indices k_mm, k_pp in accelerator version
!
! 1019 2012-09-28 06:46:45Z raasch
! small change in comment lines
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator versions (*_acc) added
!
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! 888 2012-04-20 15:03:46Z suehring
! Number of IBITS() calls with identical arguments is reduced.
!
! 862 2012-03-26 14:21:38Z suehring
! ws-scheme also work with topography in combination with vector version.
! ws-scheme also work with outflow boundaries in combination with
! vector version.
! Degradation of the applied order of scheme is now steered by multiplying with
! Integer wall_flags_0. 2nd order scheme, WS3 and WS5 are calculated on each
! grid point and mulitplied with the appropriate flag.
! 2nd order numerical dissipation term changed. Now the appropriate 2nd order
! term derived according to the 4th and 6th order terms is applied. It turns
! out that diss_2nd does not provide sufficient dissipation near walls.
! Therefore, the function diss_2nd is removed.
! Near walls a divergence correction is necessary to overcome numerical
! instabilities due to too less divergence reduction of the velocity field.
! boundary_flags and logicals steering the degradation are removed.
! Empty SUBROUTINE local_diss removed.
! Further formatting adjustments.
!
! 801 2012-01-10 17:30:36Z suehring
! Bugfix concerning OpenMP parallelization. Summation of sums_wsus_ws_l,
! sums_wsvs_ws_l, sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l, sums_wspts_ws_l,
! sums_wsqs_ws_l, sums_wssas_ws_l is now thread-safe by adding an additional
! dimension.
!
! Initial revision
!
! 411 2009-12-11 12:31:43 Z suehring
!
! Description:
! ------------
!> Advection scheme for scalars and momentum using the flux formulation of
!> Wicker and Skamarock 5th order. Additionally the module contains of a 
!> routine using for initialisation and steering of the statical evaluation. 
!> The computation of turbulent fluxes takes place inside the advection
!> routines.
!> Near non-cyclic boundaries the order of the applied advection scheme is
!> degraded.
!> A divergence correction is applied. It is necessary for topography, since
!> the divergence is not sufficiently reduced, resulting in erroneous fluxes and
!> partly numerical instabilities. 
!-----------------------------------------------------------------------------!
 MODULE advec_ws

 

    PRIVATE
    PUBLIC   advec_s_ws, advec_s_ws_acc, advec_u_ws, advec_u_ws_acc,          &
             advec_v_ws, advec_v_ws_acc, advec_w_ws, advec_w_ws_acc,          &
             ws_init, ws_init_flags, ws_statistics

    INTERFACE ws_init
       MODULE PROCEDURE ws_init
    END INTERFACE ws_init

    INTERFACE ws_init_flags
       MODULE PROCEDURE ws_init_flags
    END INTERFACE ws_init_flags

    INTERFACE ws_statistics
       MODULE PROCEDURE ws_statistics
    END INTERFACE ws_statistics

    INTERFACE advec_s_ws
       MODULE PROCEDURE advec_s_ws
       MODULE PROCEDURE advec_s_ws_ij
    END INTERFACE advec_s_ws

    INTERFACE advec_u_ws
       MODULE PROCEDURE advec_u_ws
       MODULE PROCEDURE advec_u_ws_ij
    END INTERFACE advec_u_ws

    INTERFACE advec_u_ws_acc
       MODULE PROCEDURE advec_u_ws_acc
    END INTERFACE advec_u_ws_acc

    INTERFACE advec_v_ws
       MODULE PROCEDURE advec_v_ws
       MODULE PROCEDURE advec_v_ws_ij
    END INTERFACE advec_v_ws

    INTERFACE advec_v_ws_acc
       MODULE PROCEDURE advec_v_ws_acc
    END INTERFACE advec_v_ws_acc

    INTERFACE advec_w_ws
       MODULE PROCEDURE advec_w_ws
       MODULE PROCEDURE advec_w_ws_ij
    END INTERFACE advec_w_ws

    INTERFACE advec_w_ws_acc
       MODULE PROCEDURE advec_w_ws_acc
    END INTERFACE advec_w_ws_acc

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of WS-scheme
!------------------------------------------------------------------------------!
    SUBROUTINE ws_init

       USE arrays_3d,                                                          &
           ONLY:  diss_l_e, diss_l_nr, diss_l_pt, diss_l_q, diss_l_qr,         &
                  diss_l_s, diss_l_sa, diss_l_u, diss_l_v, diss_l_w,  flux_l_e,&
                  flux_l_nr, flux_l_pt, flux_l_q, flux_l_qr, flux_l_s,         &
                  flux_l_sa, flux_l_u, flux_l_v, flux_l_w, diss_s_e, diss_s_nr,&
                  diss_s_pt, diss_s_q, diss_s_qr, diss_s_s, diss_s_sa,         &
                  diss_s_u, diss_s_v, diss_s_w, flux_s_e, flux_s_nr, flux_s_pt,&
                  flux_s_q, flux_s_qr, flux_s_s, flux_s_sa, flux_s_u, flux_s_v,&
                  flux_s_w

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5, adv_sca_1, adv_sca_3,       &
                  adv_sca_5

       USE control_parameters,                                                 &
           ONLY:  cloud_physics, humidity, loop_optimization,                  &
                  monotonic_adjustment, passive_scalar, microphysics_seifert,  &
                  ocean, ws_scheme_mom, ws_scheme_sca

       USE indices,                                                            &
           ONLY:  nyn, nys, nzb, nzt

       USE kinds
       
       USE pegrid

       USE statistics,                                                         &
           ONLY:  sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l, sums_wsnrs_ws_l,&
                  sums_wspts_ws_l, sums_wsqrs_ws_l, sums_wsqs_ws_l,            &
                  sums_wsss_ws_l, sums_wssas_ws_l,  sums_wsss_ws_l,               &
                  sums_wsus_ws_l, sums_wsvs_ws_l  

!
!--    Set the appropriate factors for scalar and momentum advection.
       adv_sca_5 = 1.0_wp /  60.0_wp
       adv_sca_3 = 1.0_wp /  12.0_wp
       adv_sca_1 = 1.0_wp /   2.0_wp
       adv_mom_5 = 1.0_wp / 120.0_wp
       adv_mom_3 = 1.0_wp /  24.0_wp
       adv_mom_1 = 1.0_wp /   4.0_wp
!         
!--    Arrays needed for statical evaluation of fluxes.
       IF ( ws_scheme_mom )  THEN

          ALLOCATE( sums_wsus_ws_l(nzb:nzt+1,0:threads_per_task-1),            &
                    sums_wsvs_ws_l(nzb:nzt+1,0:threads_per_task-1),            &
                    sums_us2_ws_l(nzb:nzt+1,0:threads_per_task-1),             &
                    sums_vs2_ws_l(nzb:nzt+1,0:threads_per_task-1),             &
                    sums_ws2_ws_l(nzb:nzt+1,0:threads_per_task-1) )

          sums_wsus_ws_l = 0.0_wp
          sums_wsvs_ws_l = 0.0_wp
          sums_us2_ws_l  = 0.0_wp
          sums_vs2_ws_l  = 0.0_wp
          sums_ws2_ws_l  = 0.0_wp

       ENDIF

       IF ( ws_scheme_sca )  THEN

          ALLOCATE( sums_wspts_ws_l(nzb:nzt+1,0:threads_per_task-1) )
          sums_wspts_ws_l = 0.0_wp

          IF ( humidity  )  THEN
             ALLOCATE( sums_wsqs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             sums_wsqs_ws_l = 0.0_wp
          ENDIF
          
          IF ( passive_scalar )  THEN
             ALLOCATE( sums_wsss_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             sums_wsss_ws_l = 0.0_wp
          ENDIF

          IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
             ALLOCATE( sums_wsqrs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             ALLOCATE( sums_wsnrs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             sums_wsqrs_ws_l = 0.0_wp
             sums_wsnrs_ws_l = 0.0_wp
          ENDIF

          IF ( ocean )  THEN
             ALLOCATE( sums_wssas_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             sums_wssas_ws_l = 0.0_wp
          ENDIF

       ENDIF

!
!--    Arrays needed for reasons of speed optimization for cache version.
!--    For the vector version the buffer arrays are not necessary,
!--    because the the fluxes can swapped directly inside the loops of the 
!--    advection routines.
       IF ( loop_optimization /= 'vector' )  THEN

          IF ( ws_scheme_mom )  THEN

             ALLOCATE( flux_s_u(nzb+1:nzt,0:threads_per_task-1),               &
                       flux_s_v(nzb+1:nzt,0:threads_per_task-1),               &
                       flux_s_w(nzb+1:nzt,0:threads_per_task-1),               &
                       diss_s_u(nzb+1:nzt,0:threads_per_task-1),               &
                       diss_s_v(nzb+1:nzt,0:threads_per_task-1),               &
                       diss_s_w(nzb+1:nzt,0:threads_per_task-1) )
             ALLOCATE( flux_l_u(nzb+1:nzt,nys:nyn,0:threads_per_task-1),       &
                       flux_l_v(nzb+1:nzt,nys:nyn,0:threads_per_task-1),       &
                       flux_l_w(nzb+1:nzt,nys:nyn,0:threads_per_task-1),       &
                       diss_l_u(nzb+1:nzt,nys:nyn,0:threads_per_task-1),       &
                       diss_l_v(nzb+1:nzt,nys:nyn,0:threads_per_task-1),       &
                       diss_l_w(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )

          ENDIF

          IF ( ws_scheme_sca )  THEN

             ALLOCATE( flux_s_pt(nzb+1:nzt,0:threads_per_task-1),              &
                       flux_s_e(nzb+1:nzt,0:threads_per_task-1),               &
                       diss_s_pt(nzb+1:nzt,0:threads_per_task-1),              &
                       diss_s_e(nzb+1:nzt,0:threads_per_task-1) ) 
             ALLOCATE( flux_l_pt(nzb+1:nzt,nys:nyn,0:threads_per_task-1),      &
                       flux_l_e(nzb+1:nzt,nys:nyn,0:threads_per_task-1),       &
                       diss_l_pt(nzb+1:nzt,nys:nyn,0:threads_per_task-1),      &
                       diss_l_e(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )

             IF ( humidity )  THEN
                ALLOCATE( flux_s_q(nzb+1:nzt,0:threads_per_task-1),            &
                          diss_s_q(nzb+1:nzt,0:threads_per_task-1) )
                ALLOCATE( flux_l_q(nzb+1:nzt,nys:nyn,0:threads_per_task-1),    &
                          diss_l_q(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
             ENDIF
             
             IF ( passive_scalar )  THEN
                ALLOCATE( flux_s_s(nzb+1:nzt,0:threads_per_task-1),            &
                          diss_s_s(nzb+1:nzt,0:threads_per_task-1) )
                ALLOCATE( flux_l_s(nzb+1:nzt,nys:nyn,0:threads_per_task-1),    &
                          diss_l_s(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
             ENDIF             

             IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
                ALLOCATE( flux_s_qr(nzb+1:nzt,0:threads_per_task-1),           &
                          diss_s_qr(nzb+1:nzt,0:threads_per_task-1),           &
                          flux_s_nr(nzb+1:nzt,0:threads_per_task-1),           &
                          diss_s_nr(nzb+1:nzt,0:threads_per_task-1) )
                ALLOCATE( flux_l_qr(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          diss_l_qr(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          flux_l_nr(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          diss_l_nr(nzb+1:nzt,nys:nyn,0:threads_per_task-1) ) 
             ENDIF

             IF ( ocean )  THEN
                ALLOCATE( flux_s_sa(nzb+1:nzt,0:threads_per_task-1),           &
                          diss_s_sa(nzb+1:nzt,0:threads_per_task-1) )
                ALLOCATE( flux_l_sa(nzb+1:nzt,nys:nyn,0:threads_per_task-1),   &
                          diss_l_sa(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
             ENDIF

          ENDIF

       ENDIF

    END SUBROUTINE ws_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of flags for WS-scheme used to degrade the order of the scheme
!> near walls.
!------------------------------------------------------------------------------!
    SUBROUTINE ws_init_flags

       USE control_parameters,                                                 &
           ONLY:  inflow_l, inflow_n, inflow_r, inflow_s, momentum_advec,      &
                  nest_bound_l, nest_bound_n, nest_bound_r, nest_bound_s,      &
                  outflow_l, outflow_n, outflow_r, outflow_s, scalar_advec

       USE indices,                                                            &
           ONLY:  nbgp, nxl, nxlu, nxr, nyn, nys, nysv, nzb, nzb_s_inner,      &
                  nzb_u_inner, nzb_v_inner, nzb_w_inner, nzt, wall_flags_0,    &
                  wall_flags_00

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i  !< index variable along x
       INTEGER(iwp) ::  j  !< index variable along y
       INTEGER(iwp) ::  k  !< index variable along z

       LOGICAL      ::  flag_set !< steering variable for advection flags
   

       IF ( scalar_advec == 'ws-scheme' .OR.                                   &
            scalar_advec == 'ws-scheme-mono' )  THEN
!
!--       Set flags to steer the degradation of the advection scheme in advec_ws
!--       near topography, inflow- and outflow boundaries as well as bottom and
!--       top of model domain. wall_flags_0 remains zero for all non-prognostic
!--       grid points.
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb_s_inner(j,i)+1, nzt
!
!--                scalar - x-direction
!--                WS1 (0), WS3 (1), WS5 (2)
                   IF ( ( k <= nzb_s_inner(j,i+1) .OR. k <= nzb_s_inner(j,i+2) &   
                     .OR. k <= nzb_s_inner(j,i-1) )                            &
                       .OR. ( ( inflow_l .OR. outflow_l .OR. nest_bound_l )    &
                              .AND. i == nxl   )    .OR.                       &
                            ( ( inflow_r .OR. outflow_r .OR. nest_bound_r )    &
                              .AND. i == nxr   ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 0 )
                   ELSEIF ( ( k <= nzb_s_inner(j,i+3) .AND. k > nzb_s_inner(j,i+1)&
                                                      .AND. k > nzb_s_inner(j,i+2)&
                                                      .AND. k > nzb_s_inner(j,i-1)&
                            )                       .OR.                          &
                            ( k <= nzb_s_inner(j,i-2) .AND. k > nzb_s_inner(j,i+1)&
                                                      .AND. k > nzb_s_inner(j,i+2)&
                                                      .AND. k > nzb_s_inner(j,i-1)&
                            )                                                     &
                                                    .OR.                          &
                            ( ( inflow_r .OR. outflow_r .OR. nest_bound_r )       &
                              .AND. i == nxr-1 )    .OR.                          &
                            ( ( inflow_l .OR. outflow_l .OR. nest_bound_l )       &
                              .AND. i == nxlu  ) )                                &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 1 )
                   ELSEIF ( k > nzb_s_inner(j,i+1) .AND. k > nzb_s_inner(j,i+2)&
                      .AND. k > nzb_s_inner(j,i+3) .AND. k > nzb_s_inner(j,i-1)&
                      .AND. k > nzb_s_inner(j,i-2) )                           &
                   THEN 
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 2 )
                   ENDIF
!
!--                scalar - y-direction
!--                WS1 (3), WS3 (4), WS5 (5)
                   IF ( ( k <= nzb_s_inner(j+1,i)  .OR. k <= nzb_s_inner(j+2,i)   &   
                                                   .OR. k <= nzb_s_inner(j-1,i) ) &
                                                    .OR.                          &
                            ( ( inflow_s .OR. outflow_s .OR. nest_bound_s )       &
                              .AND. j == nys   )    .OR.                          &
                            ( ( inflow_n .OR. outflow_n .OR. nest_bound_n )       &
                              .AND. j == nyn   ) )                                &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 3 )
!
!--                WS3
                   ELSEIF ( ( k <= nzb_s_inner(j+3,i) .AND. k > nzb_s_inner(j+1,i)&
                                                      .AND. k > nzb_s_inner(j+2,i)&
                                                      .AND. k > nzb_s_inner(j-1,i)&
                            )                           .OR.                      &
                            ( k <= nzb_s_inner(j-2,i) .AND. k > nzb_s_inner(j+1,i)&
                                                      .AND. k > nzb_s_inner(j+2,i)&
                                                      .AND. k > nzb_s_inner(j-1,i)&
                            )                                                     &
                                                        .OR.                      &
                            ( ( inflow_s .OR. outflow_s .OR. nest_bound_s )       &
                              .AND. j == nysv  )    .OR.                          &
                            ( ( inflow_n .OR. outflow_n .OR. nest_bound_n )       &
                              .AND. j == nyn-1 ) )                                &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 4 )
!
!--                WS5
                   ELSEIF ( k > nzb_s_inner(j+1,i) .AND. k > nzb_s_inner(j+2,i)&
                      .AND. k > nzb_s_inner(j+3,i) .AND. k > nzb_s_inner(j-1,i)&
                      .AND. k > nzb_s_inner(j-2,i) )                           &
                   THEN 
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 5 )
                   ENDIF
! 
!--                scalar - z-direction
!--                WS1 (6), WS3 (7), WS5 (8)
                   flag_set = .FALSE.
                   IF ( k == nzb_s_inner(j,i) + 1 .OR. k == nzt )  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 6 )
                      flag_set = .TRUE.
                   ELSEIF ( k == nzb_s_inner(j,i) + 2 .OR. k == nzt - 1 )  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 7 )
                      flag_set = .TRUE.
                   ELSEIF ( k > nzb_s_inner(j,i) .AND. .NOT. flag_set )  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 8 )
                   ENDIF

                ENDDO
             ENDDO
          ENDDO
       ENDIF

       IF ( momentum_advec == 'ws-scheme' )  THEN
!
!--       Set wall_flags_0 to steer the degradation of the advection scheme in advec_ws
!--       near topography, inflow- and outflow boundaries as well as bottom and
!--       top of model domain. wall_flags_0 remains zero for all non-prognostic
!--       grid points.
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
!  
!--                At first, set flags to WS1. 
!--                Since fluxes are swapped in advec_ws.f90, this is necessary to 
!--                in order to handle the left/south flux.
!--                near vertical walls. 
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 9 )
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 12 )
!
!--                u component - x-direction
!--                WS1 (9), WS3 (10), WS5 (11)
                   IF ( k <= nzb_u_inner(j,i+1)     .OR.                       &
                            ( ( inflow_l .OR. outflow_l .OR. nest_bound_l )    &
                              .AND. i <= nxlu  )    .OR.                       &
                            ( ( inflow_r .OR. outflow_r .OR. nest_bound_r )    &
                              .AND. i == nxr   ) )                             &
                   THEN
                       wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 9 )
                   ELSEIF ( ( k <= nzb_u_inner(j,i+2) .AND.                    &
                              k >  nzb_u_inner(j,i+1) ) .OR.                   &
                              k <= nzb_u_inner(j,i-1)                          &
                                                        .OR.                   &
                            ( ( inflow_r .OR. outflow_r .OR. nest_bound_r )    &
                              .AND. i == nxr-1 )    .OR.                       &
                            ( ( inflow_l .OR. outflow_l .OR. nest_bound_l )    &
                              .AND. i == nxlu+1) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 10 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 9 )
                   ELSEIF ( k > nzb_u_inner(j,i+1) .AND. k > nzb_u_inner(j,i+2)   &
                                                   .AND. k > nzb_u_inner(j,i-1) ) &
                   THEN    
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 11 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 9 )
                   ENDIF
!
!--                u component - y-direction
!--                WS1 (12), WS3 (13), WS5 (14)
                   IF ( k <= nzb_u_inner(j+1,i) .OR.                           &
                            ( ( inflow_s .OR. outflow_s .OR. nest_bound_s )    &
                              .AND. j == nys   )    .OR.                       &
                            ( ( inflow_n .OR. outflow_n .OR. nest_bound_n )    &
                              .AND. j == nyn   ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 12 )
                   ELSEIF ( ( k <= nzb_u_inner(j+2,i) .AND.                    &
                              k >  nzb_u_inner(j+1,i) ) .OR.                   &
                              k <= nzb_u_inner(j-1,i)                          &
                                                        .OR.                   &
                            ( ( inflow_s .OR. outflow_s .OR. nest_bound_s )    &
                              .AND. j == nysv  )    .OR.                       &
                            ( ( inflow_n .OR. outflow_n .OR. nest_bound_n )    &
                              .AND. j == nyn-1 ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 13 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 12 )
                   ELSEIF ( k > nzb_u_inner(j+1,i) .AND. k > nzb_u_inner(j+2,i)   &
                                                   .AND. k > nzb_u_inner(j-1,i) ) &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 14 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 12 )
                   ENDIF
!
!--                u component - z-direction
!--                WS1 (15), WS3 (16), WS5 (17)
                   flag_set = .FALSE.
                   IF ( k == nzb_u_inner(j,i) + 1 .OR. k == nzt )  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 15 )
                      flag_set = .TRUE.
                   ELSEIF ( k == nzb_u_inner(j,i) + 2 .OR. k == nzt - 1 )  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 16 )
                      flag_set = .TRUE.
                   ELSEIF ( k > nzb_u_inner(j,i) .AND. .NOT. flag_set )  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 17 )
                   ENDIF

                ENDDO
             ENDDO
          ENDDO

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
!
!--                At first, set flags to WS1. 
!--                Since fluxes are swapped in advec_ws.f90, this is necessary to 
!--                in order to handle the left/south flux.
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 18 )
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 21 )
!
!--                v component - x-direction
!--                WS1 (18), WS3 (19), WS5 (20)
                   IF ( k <= nzb_v_inner(j,i+1) .OR.                           &
                            ( ( inflow_l .OR. outflow_l .OR. nest_bound_l )    &
                              .AND. i == nxl   )    .OR.                       &
                            ( ( inflow_r .OR. outflow_r .OR. nest_bound_r )    &
                              .AND. i == nxr   ) )                             &
                  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 18 )
!
!--                WS3
                   ELSEIF ( ( k <= nzb_v_inner(j,i+2) .AND.                    &
                              k >  nzb_v_inner(j,i+1) ) .OR.                   &
                              k <= nzb_v_inner(j,i-1)                          &
                                                    .OR.                       &
                            ( ( inflow_r .OR. outflow_r .OR. nest_bound_r )    &
                              .AND. i == nxr-1 )    .OR.                       &
                            ( ( inflow_l .OR. outflow_l .OR. nest_bound_l )    &
                              .AND. i == nxlu  ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 19 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 18 )
                   ELSEIF ( k > nzb_v_inner(j,i+1) .AND. k > nzb_v_inner(j,i+2)   &
                                                   .AND. k > nzb_v_inner(j,i-1) ) &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 20 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 18 )
                   ENDIF
!
!--                v component - y-direction
!--                WS1 (21), WS3 (22), WS5 (23)
                   IF ( k <= nzb_v_inner(j+1,i) .OR.                           &
                            ( ( inflow_s .OR. outflow_s .OR. nest_bound_s )    &
                              .AND. j <= nysv  )    .OR.                       &
                            ( ( inflow_n .OR. outflow_n .OR. nest_bound_n )    &
                              .AND. j == nyn   ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 21 )
                   ELSEIF ( ( k <= nzb_v_inner(j+2,i) .AND.                    &
                              k >  nzb_v_inner(j+1,i) ) .OR.                   &
                              k <= nzb_v_inner(j-1,i)                          &
                                                        .OR.                   &
                            ( ( inflow_s .OR. outflow_s .OR. nest_bound_s )    &
                              .AND. j == nysv+1)    .OR.                       &
                            ( ( inflow_n .OR. outflow_n .OR. nest_bound_n )    &
                              .AND. j == nyn-1 ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 22 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 21 )
                   ELSEIF ( k > nzb_v_inner(j+1,i) .AND. k > nzb_v_inner(j+2,i)   &
                                                   .AND. k > nzb_v_inner(j-1,i) ) &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 23 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 21 )
                   ENDIF
!
!--                v component - z-direction
!--                WS1 (24), WS3 (25), WS5 (26)
                   flag_set = .FALSE.
                   IF ( k == nzb_v_inner(j,i) + 1 .OR. k == nzt )  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 24 )
                      flag_set = .TRUE.
                   ELSEIF ( k == nzb_v_inner(j,i) + 2 .OR. k == nzt - 1 )  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 25 )
                      flag_set = .TRUE.
                   ELSEIF ( k > nzb_v_inner(j,i) .AND. .NOT. flag_set )  THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 26 )
                   ENDIF

                ENDDO
             ENDDO
          ENDDO
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
!
!--                At first, set flags to WS1. 
!--                Since fluxes are swapped in advec_ws.f90, this is necessary to 
!--                in order to handle the left/south flux.
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 27 )
                   wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 30 )
!
!--                w component - x-direction
!--                WS1 (27), WS3 (28), WS5 (29)
                   IF ( k <= nzb_w_inner(j,i+1) .OR.                           &
                            ( ( inflow_l .OR. outflow_l .OR. nest_bound_l )    &
                              .AND. i == nxl   )    .OR.                       &
                            ( ( inflow_r .OR. outflow_r .OR. nest_bound_r )    &
                              .AND. i == nxr   ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 27 )
                   ELSEIF ( ( k <= nzb_w_inner(j,i+2) .AND.                    &
                              k >  nzb_w_inner(j,i+1) ) .OR.                   &
                              k <= nzb_w_inner(j,i-1)                          &
                                                        .OR.                   &
                            ( ( inflow_r .OR. outflow_r .OR. nest_bound_r )    &
                              .AND. i == nxr-1 )    .OR.                       &
                            ( ( inflow_l .OR. outflow_l .OR. nest_bound_l )    &
                              .AND. i == nxlu  ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 28 )
!   
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 27 )
                   ELSEIF ( k > nzb_w_inner(j,i+1) .AND. k > nzb_w_inner(j,i+2)   &
                                                   .AND. k > nzb_w_inner(j,i-1) ) &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i),29 )
!   
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 27 )
                   ENDIF
!
!--                w component - y-direction
!--                WS1 (30), WS3 (31), WS5 (32)
                   IF ( k <= nzb_w_inner(j+1,i) .OR.                           &
                            ( ( inflow_s .OR. outflow_s .OR. nest_bound_s )    &
                              .AND. j == nys   )    .OR.                       &
                            ( ( inflow_n .OR. outflow_n .OR. nest_bound_n )    &
                              .AND. j == nyn   ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 30 )
                   ELSEIF ( ( k <= nzb_w_inner(j+2,i) .AND.                    &
                              k >  nzb_w_inner(j+1,i) ) .OR.                   &
                               k <= nzb_w_inner(j-1,i)                         &
                                                        .OR.                   &
                            ( ( inflow_s .OR. outflow_s .OR. nest_bound_s )    &
                              .AND. j == nysv  )    .OR.                       &
                            ( ( inflow_n .OR. outflow_n .OR. nest_bound_n )    &
                              .AND. j == nyn-1 ) )                             &
                   THEN
                      wall_flags_0(k,j,i) = IBSET( wall_flags_0(k,j,i), 31 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 30 )
                   ELSEIF ( k > nzb_w_inner(j+1,i) .AND. k > nzb_w_inner(j+2,i)   &
                                                   .AND. k > nzb_w_inner(j-1,i) ) &
                   THEN
                      wall_flags_00(k,j,i) = IBSET( wall_flags_00(k,j,i), 0 )
!
!--                   Clear flag for WS1
                      wall_flags_0(k,j,i) = IBCLR( wall_flags_0(k,j,i), 30 )
                   ENDIF
!
!--                w component - z-direction
!--                WS1 (33), WS3 (34), WS5 (35)
                   flag_set = .FALSE.
                   IF ( k == nzb_w_inner(j,i) .OR. k == nzb_w_inner(j,i) + 1   &
                                              .OR. k == nzt )  THEN
!
!--                   Please note, at k == nzb_w_inner(j,i) a flag is explictely
!--                   set, although this is not a prognostic level. However, 
!--                   contrary to the advection of u,v and s this is necessary
!--                   because flux_t(nzb_w_inner(j,i)) is used for the tendency
!--                   at k == nzb_w_inner(j,i)+1.
                      wall_flags_00(k,j,i) = IBSET( wall_flags_00(k,j,i), 1 )
                      flag_set = .TRUE.
                   ELSEIF ( k == nzb_w_inner(j,i) + 2 .OR. k == nzt - 1 )  THEN
                      wall_flags_00(k,j,i) = IBSET( wall_flags_00(k,j,i), 2 )
                      flag_set = .TRUE.
                   ELSEIF ( k > nzb_w_inner(j,i) .AND. .NOT. flag_set )  THEN
                      wall_flags_00(k,j,i) = IBSET( wall_flags_00(k,j,i), 3 )
                   ENDIF

                ENDDO
             ENDDO
          ENDDO

       ENDIF


!
!--    Exchange 3D integer wall_flags.
       IF ( momentum_advec == 'ws-scheme' .OR. scalar_advec == 'ws-scheme'  &
       .OR. scalar_advec == 'ws-scheme-mono' )  THEN  
!
!--       Exchange ghost points for advection flags
          CALL exchange_horiz_int( wall_flags_0,  nbgp )
          CALL exchange_horiz_int( wall_flags_00, nbgp )
!
!--       Set boundary flags at inflow and outflow boundary in case of 
!--       non-cyclic boundary conditions.
         IF ( inflow_l .OR. outflow_l .OR. nest_bound_l )  THEN
             wall_flags_0(:,:,nxl-1)  = wall_flags_0(:,:,nxl)
             wall_flags_00(:,:,nxl-1) = wall_flags_00(:,:,nxl)
         ENDIF

         IF ( inflow_r .OR. outflow_r .OR. nest_bound_r )  THEN
            wall_flags_0(:,:,nxr+1)  = wall_flags_0(:,:,nxr)
            wall_flags_00(:,:,nxr+1) = wall_flags_00(:,:,nxr)
          ENDIF

          IF ( inflow_n .OR. outflow_n .OR. nest_bound_n )  THEN
             wall_flags_0(:,nyn+1,:)  = wall_flags_0(:,nyn,:)
             wall_flags_00(:,nyn+1,:) = wall_flags_00(:,nyn,:)
          ENDIF

          IF ( inflow_s .OR. outflow_s  .OR. nest_bound_s )  THEN
             wall_flags_0(:,nys-1,:)  = wall_flags_0(:,nys,:)
             wall_flags_00(:,nys-1,:) = wall_flags_00(:,nys,:)
          ENDIF
  
       ENDIF


    END SUBROUTINE ws_init_flags


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize variables used for storing statistic quantities (fluxes, variances)
!------------------------------------------------------------------------------!
    SUBROUTINE ws_statistics
    
       USE control_parameters,                                                 &
           ONLY:  cloud_physics, humidity, passive_scalar, ocean,              &
                  microphysics_seifert, ws_scheme_mom, ws_scheme_sca

       USE kinds

       USE statistics,                                                         &
           ONLY:  sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l, sums_wsnrs_ws_l,&
                  sums_wspts_ws_l, sums_wsqrs_ws_l, sums_wsqs_ws_l,            &
                  sums_wsss_ws_l, sums_wssas_ws_l,  sums_wsus_ws_l,            &
                  sums_wsvs_ws_l  

       IMPLICIT NONE

!       
!--    The arrays needed for statistical evaluation are set to to 0 at the 
!--    beginning of prognostic_equations.
       IF ( ws_scheme_mom )  THEN
          sums_wsus_ws_l = 0.0_wp
          sums_wsvs_ws_l = 0.0_wp
          sums_us2_ws_l  = 0.0_wp
          sums_vs2_ws_l  = 0.0_wp
          sums_ws2_ws_l  = 0.0_wp
       ENDIF

       IF ( ws_scheme_sca )  THEN
          sums_wspts_ws_l = 0.0_wp
          IF ( humidity       )  sums_wsqs_ws_l = 0.0_wp
          IF ( passive_scalar )  sums_wsss_ws_l = 0.0_wp
          IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
             sums_wsqrs_ws_l = 0.0_wp
             sums_wsnrs_ws_l = 0.0_wp
          ENDIF
          IF ( ocean )  sums_wssas_ws_l = 0.0_wp

       ENDIF

    END SUBROUTINE ws_statistics


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Scalar advection - Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_ws_ij( i, j, sk, sk_char, swap_flux_y_local,            &
                              swap_diss_y_local, swap_flux_x_local,            &
                              swap_diss_x_local, i_omp, tn )

       USE arrays_3d,                                                          &
           ONLY:  ddzw, drho_air, tend, u, v, w, rho_air, rho_air_zw

       USE constants,                                                          &
           ONLY:  adv_sca_1, adv_sca_3, adv_sca_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count, monotonic_adjustment, u_gtrans, &
                  v_gtrans 

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzb_max,    &
                  nzt, wall_flags_0

       USE kinds

       USE pegrid

       USE statistics,                                                         &
           ONLY:  hom, sums_wsnrs_ws_l, sums_wspts_ws_l, sums_wsqrs_ws_l,      &
                  sums_wsqs_ws_l, sums_wssas_ws_l, sums_wsss_ws_l,             &
                  weight_substep

       IMPLICIT NONE

       CHARACTER (LEN = *), INTENT(IN) ::  sk_char !<
       
       INTEGER(iwp) ::  i     !<
       INTEGER(iwp) ::  ibit0 !<
       INTEGER(iwp) ::  ibit1 !<
       INTEGER(iwp) ::  ibit2 !<
       INTEGER(iwp) ::  ibit3 !<
       INTEGER(iwp) ::  ibit4 !<
       INTEGER(iwp) ::  ibit5 !<
       INTEGER(iwp) ::  ibit6 !<
       INTEGER(iwp) ::  ibit7 !<
       INTEGER(iwp) ::  ibit8 !<
       INTEGER(iwp) ::  i_omp !<
       INTEGER(iwp) ::  j     !<
       INTEGER(iwp) ::  k     !<
       INTEGER(iwp) ::  k_mm  !<
       INTEGER(iwp) ::  k_mmm !<
       INTEGER(iwp) ::  k_pp  !<
       INTEGER(iwp) ::  k_ppp !<
       INTEGER(iwp) ::  tn    !<
       
       REAL(wp)     ::  diss_d !<
       REAL(wp)     ::  div    !<
       REAL(wp)     ::  flux_d !<
       REAL(wp)     ::  fd_1   !<
       REAL(wp)     ::  fl_1   !<
       REAL(wp)     ::  fn_1   !<
       REAL(wp)     ::  fr_1   !<
       REAL(wp)     ::  fs_1   !<
       REAL(wp)     ::  ft_1   !<
       REAL(wp)     ::  phi_d  !<
       REAL(wp)     ::  phi_l  !<
       REAL(wp)     ::  phi_n  !<
       REAL(wp)     ::  phi_r  !<
       REAL(wp)     ::  phi_s  !<
       REAL(wp)     ::  phi_t  !<
       REAL(wp)     ::  rd     !<
       REAL(wp)     ::  rl     !<
       REAL(wp)     ::  rn     !<
       REAL(wp)     ::  rr     !<
       REAL(wp)     ::  rs     !<
       REAL(wp)     ::  rt     !<
       REAL(wp)     ::  u_comp !<
       REAL(wp)     ::  v_comp !<
       
#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER    ::  sk     !<
#endif
       REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_n !<
       REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_r !<
       REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_t !<
       REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_n !<
       REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_r !<
       REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_t !<
       
       REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1) ::  swap_diss_y_local !<
       REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1) ::  swap_flux_y_local !<
       
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) ::  swap_diss_x_local !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) ::  swap_flux_x_local !<
       

!
!--    Compute southside fluxes of the respective PE bounds.
       IF ( j == nys )  THEN
!
!--       Up to the top of the highest topography.
          DO  k = nzb+1, nzb_max

             ibit5 = IBITS(wall_flags_0(k,j-1,i),5,1)
             ibit4 = IBITS(wall_flags_0(k,j-1,i),4,1)
             ibit3 = IBITS(wall_flags_0(k,j-1,i),3,1)

             v_comp                  = v(k,j,i) - v_gtrans
             swap_flux_y_local(k,tn) = v_comp *         (                     &
                                               ( 37.0_wp * ibit5 * adv_sca_5  &
                                            +     7.0_wp * ibit4 * adv_sca_3  &
                                            +              ibit3 * adv_sca_1  &
                                               ) *                            &
                                           ( sk(k,j,i)  + sk(k,j-1,i)     )   &
                                         -     (  8.0_wp * ibit5 * adv_sca_5  &
                                            +              ibit4 * adv_sca_3  &
                                                ) *                           &
                                           ( sk(k,j+1,i) + sk(k,j-2,i)    )   &
                                         +     (           ibit5 * adv_sca_5  &
                                               ) *                            &
                                           ( sk(k,j+2,i) + sk(k,j-3,i)    )   &
                                                        )

             swap_diss_y_local(k,tn) = -ABS( v_comp ) * (                     &
                                               ( 10.0_wp * ibit5 * adv_sca_5  &
                                            +     3.0_wp * ibit4 * adv_sca_3  &
                                            +              ibit3 * adv_sca_1  &
                                               ) *                            &
                                            ( sk(k,j,i)   - sk(k,j-1,i)  )    &
                                        -      (  5.0_wp * ibit5 * adv_sca_5  &
                                            +              ibit4 * adv_sca_3  &
                                            ) *                               &
                                            ( sk(k,j+1,i) - sk(k,j-2,i)  )    &
                                        +      (           ibit5 * adv_sca_5  &
                                               ) *                            &
                                            ( sk(k,j+2,i) - sk(k,j-3,i)  )    &
                                                        )

          ENDDO
!
!--       Above to the top of the highest topography. No degradation necessary.
          DO  k = nzb_max+1, nzt

             v_comp                  = v(k,j,i) - v_gtrans
             swap_flux_y_local(k,tn) = v_comp * (                             &
                                    37.0_wp * ( sk(k,j,i)   + sk(k,j-1,i) )   &
                                  -  8.0_wp * ( sk(k,j+1,i) + sk(k,j-2,i) )   &
                                  +           ( sk(k,j+2,i) + sk(k,j-3,i) )   &
                                                ) * adv_sca_5
              swap_diss_y_local(k,tn) = -ABS( v_comp ) * (                    &
                                    10.0_wp * ( sk(k,j,i)   - sk(k,j-1,i) )   &
                                  -  5.0_wp * ( sk(k,j+1,i) - sk(k,j-2,i) )   &
                                  +             sk(k,j+2,i) - sk(k,j-3,i)     &
                                                         ) * adv_sca_5

          ENDDO

       ENDIF
!
!--    Compute leftside fluxes of the respective PE bounds.
       IF ( i == i_omp )  THEN
        
          DO  k = nzb+1, nzb_max

             ibit2 = IBITS(wall_flags_0(k,j,i-1),2,1)
             ibit1 = IBITS(wall_flags_0(k,j,i-1),1,1)
             ibit0 = IBITS(wall_flags_0(k,j,i-1),0,1)

             u_comp                     = u(k,j,i) - u_gtrans
             swap_flux_x_local(k,j,tn) = u_comp * (                           &
                                               ( 37.0_wp * ibit2 * adv_sca_5  &
                                            +     7.0_wp * ibit1 * adv_sca_3  &
                                            +              ibit0 * adv_sca_1  &
                                               ) *                            &
                                            ( sk(k,j,i)   + sk(k,j,i-1)    )  &
                                        -      (  8.0_wp * ibit2 * adv_sca_5  &
                                            +              ibit1 * adv_sca_3  &
                                               ) *                            &
                                            ( sk(k,j,i+1) + sk(k,j,i-2)    )  &
                                        +      (           ibit2 * adv_sca_5  &
                                               ) *                            &
                                            ( sk(k,j,i+2) + sk(k,j,i-3)    )  &
                                                  )

              swap_diss_x_local(k,j,tn) = -ABS( u_comp ) * (                  &
                                               ( 10.0_wp * ibit2 * adv_sca_5  &
                                            +     3.0_wp * ibit1 * adv_sca_3  &
                                            +              ibit0 * adv_sca_1  &
                                               ) *                            &
                                            ( sk(k,j,i)   - sk(k,j,i-1)    )  &
                                        -      (  5.0_wp * ibit2 * adv_sca_5  &
                                            +              ibit1 * adv_sca_3  &
                                               ) *                            &
                                            ( sk(k,j,i+1) - sk(k,j,i-2)    )  &
                                        +      (           ibit2 * adv_sca_5  &
                                               ) *                            &
                                            ( sk(k,j,i+2) - sk(k,j,i-3)    )  &
                                                           )

          ENDDO

          DO  k = nzb_max+1, nzt

             u_comp                 = u(k,j,i) - u_gtrans
             swap_flux_x_local(k,j,tn) = u_comp * (                           &
                                      37.0_wp * ( sk(k,j,i)   + sk(k,j,i-1) ) &
                                    -  8.0_wp * ( sk(k,j,i+1) + sk(k,j,i-2) ) &
                                    +           ( sk(k,j,i+2) + sk(k,j,i-3) ) &
                                                  ) * adv_sca_5

             swap_diss_x_local(k,j,tn) = -ABS( u_comp ) * (                   &
                                      10.0_wp * ( sk(k,j,i)   - sk(k,j,i-1) ) &
                                    -  5.0_wp * ( sk(k,j,i+1) - sk(k,j,i-2) ) &
                                    +           ( sk(k,j,i+2) - sk(k,j,i-3) ) &
                                                          ) * adv_sca_5

          ENDDO
           
       ENDIF

       flux_t(0) = 0.0_wp
       diss_t(0) = 0.0_wp
       flux_d    = 0.0_wp
       diss_d    = 0.0_wp
!        
!--    Now compute the fluxes and tendency terms for the horizontal and
!--    vertical parts up to the top of the highest topography.
       DO  k = nzb+1, nzb_max
!
!--       Note: It is faster to conduct all multiplications explicitly, e.g.
!--       * adv_sca_5 ... than to determine a factor and multiplicate the
!--       flux at the end. 

          ibit2 = IBITS(wall_flags_0(k,j,i),2,1)
          ibit1 = IBITS(wall_flags_0(k,j,i),1,1)
          ibit0 = IBITS(wall_flags_0(k,j,i),0,1)

          u_comp    = u(k,j,i+1) - u_gtrans
          flux_r(k) = u_comp * (                                              &
                     ( 37.0_wp * ibit2 * adv_sca_5                            &
                  +     7.0_wp * ibit1 * adv_sca_3                            &
                  +              ibit0 * adv_sca_1                            &
                     ) *                                                      &
                             ( sk(k,j,i+1) + sk(k,j,i)   )                    &
              -      (  8.0_wp * ibit2 * adv_sca_5                            &
                  +              ibit1 * adv_sca_3                            &
                     ) *                                                      &
                             ( sk(k,j,i+2) + sk(k,j,i-1) )                    &
              +      (           ibit2 * adv_sca_5                            &
                     ) *                                                      &
                             ( sk(k,j,i+3) + sk(k,j,i-2) )                    &
                               )

          diss_r(k) = -ABS( u_comp ) * (                                      &
                     ( 10.0_wp * ibit2 * adv_sca_5                            &
                  +     3.0_wp * ibit1 * adv_sca_3                            &
                  +              ibit0 * adv_sca_1                            &
                     ) *                                                      &
                             ( sk(k,j,i+1) - sk(k,j,i)  )                     &
              -      (  5.0_wp * ibit2 * adv_sca_5                            &
                  +              ibit1 * adv_sca_3                            &
                     ) *                                                      &
                             ( sk(k,j,i+2) - sk(k,j,i-1) )                    &
              +      (           ibit2 * adv_sca_5                            &
                     ) *                                                      &
                             ( sk(k,j,i+3) - sk(k,j,i-2) )                    &
                                       )

          ibit5 = IBITS(wall_flags_0(k,j,i),5,1)
          ibit4 = IBITS(wall_flags_0(k,j,i),4,1)
          ibit3 = IBITS(wall_flags_0(k,j,i),3,1)

          v_comp    = v(k,j+1,i) - v_gtrans
          flux_n(k) = v_comp * (                                              &
                     ( 37.0_wp * ibit5 * adv_sca_5                            &
                  +     7.0_wp * ibit4 * adv_sca_3                            &
                  +              ibit3 * adv_sca_1                            &
                     ) *                                                      &
                             ( sk(k,j+1,i) + sk(k,j,i)   )                    &
              -      (  8.0_wp * ibit5 * adv_sca_5                            &
                  +              ibit4 * adv_sca_3                            &
                     ) *                                                      &
                             ( sk(k,j+2,i) + sk(k,j-1,i) )                    &
              +      (           ibit5 * adv_sca_5                            &
                     ) *                                                      &
                             ( sk(k,j+3,i) + sk(k,j-2,i) )                    &
                               )

          diss_n(k) = -ABS( v_comp ) * (                                      &
                     ( 10.0_wp * ibit5 * adv_sca_5                            &
                  +     3.0_wp * ibit4 * adv_sca_3                            &
                  +              ibit3 * adv_sca_1                            &
                     ) *                                                      &
                             ( sk(k,j+1,i) - sk(k,j,i)   )                    &
              -      (  5.0_wp * ibit5 * adv_sca_5                            &
                  +              ibit4 * adv_sca_3                            &
                     ) *                                                      &
                             ( sk(k,j+2,i) - sk(k,j-1,i) )                    &
              +      (           ibit5 * adv_sca_5                            &
                     ) *                                                      &
                             ( sk(k,j+3,i) - sk(k,j-2,i) )                    &
                                       )
!
!--       k index has to be modified near bottom and top, else array
!--       subscripts will be exceeded.
          ibit8 = IBITS(wall_flags_0(k,j,i),8,1)
          ibit7 = IBITS(wall_flags_0(k,j,i),7,1)
          ibit6 = IBITS(wall_flags_0(k,j,i),6,1)

          k_ppp = k + 3 * ibit8
          k_pp  = k + 2 * ( 1 - ibit6  )
          k_mm  = k - 2 * ibit8


          flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                            &
                     ( 37.0_wp * ibit8 * adv_sca_5                            &
                  +     7.0_wp * ibit7 * adv_sca_3                            &
                  +              ibit6 * adv_sca_1                            &
                     ) *                                                      &
                             ( sk(k+1,j,i)  + sk(k,j,i)    )                  &
              -      (  8.0_wp * ibit8 * adv_sca_5                            &
                  +              ibit7 * adv_sca_3                            &
                     ) *                                                      &
                             ( sk(k_pp,j,i) + sk(k-1,j,i)  )                  &
              +      (           ibit8 * adv_sca_5                            &
                     ) *     ( sk(k_ppp,j,i)+ sk(k_mm,j,i) )                  &
                                 )

          diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                    &
                     ( 10.0_wp * ibit8 * adv_sca_5                            &
                  +     3.0_wp * ibit7 * adv_sca_3                            &
                  +              ibit6 * adv_sca_1                            &
                     ) *                                                      &
                             ( sk(k+1,j,i)   - sk(k,j,i)    )                 &
              -      (  5.0_wp * ibit8 * adv_sca_5                            &
                  +              ibit7 * adv_sca_3                            &
                     ) *                                                      &
                             ( sk(k_pp,j,i)  - sk(k-1,j,i)  )                 &
              +      (           ibit8 * adv_sca_5                            &
                     ) *                                                      &
                             ( sk(k_ppp,j,i) - sk(k_mm,j,i) )                 &
                                         )
!
!--       Apply monotonic adjustment.
          IF ( monotonic_adjustment )  THEN
!
!--          At first, calculate first order fluxes.
             u_comp = u(k,j,i) - u_gtrans
             fl_1   =  ( u_comp   * ( sk(k,j,i) + sk(k,j,i-1) )               &
                   -ABS( u_comp ) * ( sk(k,j,i) - sk(k,j,i-1) )               &
                       ) * adv_sca_1 

             u_comp = u(k,j,i+1) - u_gtrans
             fr_1   =  ( u_comp   * ( sk(k,j,i+1) + sk(k,j,i) )               &
                   -ABS( u_comp ) * ( sk(k,j,i+1) - sk(k,j,i) )               &
                       ) * adv_sca_1 

             v_comp = v(k,j,i) - v_gtrans
             fs_1   =  ( v_comp   * ( sk(k,j,i) + sk(k,j-1,i) )               &
                   -ABS( v_comp ) * ( sk(k,j,i) - sk(k,j-1,i) )               &
                       ) * adv_sca_1 

             v_comp = v(k,j+1,i) - v_gtrans
             fn_1   =  ( v_comp   * ( sk(k,j+1,i) + sk(k,j,i) )               &
                   -ABS( v_comp ) * ( sk(k,j+1,i) - sk(k,j,i) )               &
                       ) * adv_sca_1 

             fd_1   =  ( w(k-1,j,i)   * ( sk(k,j,i) + sk(k-1,j,i) )           &
                   -ABS( w(k-1,j,i) ) * ( sk(k,j,i) - sk(k-1,j,i) )           &
                       ) * adv_sca_1 * rho_air_zw(k)

             ft_1   =  ( w(k,j,i)   * ( sk(k+1,j,i) + sk(k,j,i) )             &
                   -ABS( w(k,j,i) ) * ( sk(k+1,j,i) - sk(k,j,i) )             &
                      ) * adv_sca_1 * rho_air_zw(k)
!
!--          Calculate ratio of upwind gradients. Note, Min/Max is just to 
!--          avoid if statements.
             rl     = ( MAX( 0.0_wp, u(k,j,i) - u_gtrans ) *                  & 
                           ABS( ( sk(k,j,i-1) - sk(k,j,i-2)            ) /    &
                                ( sk(k,j,i)   - sk(k,j,i-1) + 1E-20_wp )      &
                              ) +                                             & 
                        MIN( 0.0_wp, u(k,j,i) - u_gtrans ) *                  &
                           ABS( ( sk(k,j,i)   - sk(k,j,i+1)            ) /    &
                                ( sk(k,j,i-1) - sk(k,j,i)   + 1E-20_wp )      &
                              )                                               &
                      ) / ABS( u(k,j,i) - u_gtrans + 1E-20_wp )

             rr     = ( MAX( 0.0_wp, u(k,j,i+1) - u_gtrans ) *                & 
                           ABS( ( sk(k,j,i)   - sk(k,j,i-1)            ) /    &
                                ( sk(k,j,i+1) - sk(k,j,i)   + 1E-20_wp )      &
                              ) +                                             & 
                        MIN( 0.0_wp, u(k,j,i+1) - u_gtrans ) *                &
                           ABS( ( sk(k,j,i+1) - sk(k,j,i+2)            ) /    &
                                ( sk(k,j,i)   - sk(k,j,i+1) + 1E-20_wp )      &
                              )                                               &
                      ) / ABS( u(k,j,i+1) - u_gtrans + 1E-20_wp )

             rs     = ( MAX( 0.0_wp, v(k,j,i) - v_gtrans ) *                  & 
                           ABS( ( sk(k,j-1,i) - sk(k,j-2,i)            ) /    &
                                ( sk(k,j,i)   - sk(k,j-1,i) + 1E-20_wp )      &
                              ) +                                             & 
                        MIN( 0.0_wp, v(k,j,i) - v_gtrans ) *                  &
                           ABS( ( sk(k,j,i)   - sk(k,j+1,i)            ) /    &
                                ( sk(k,j-1,i) - sk(k,j,i)   + 1E-20_wp )      &
                              )                                               &
                      ) / ABS( v(k,j,i) - v_gtrans + 1E-20_wp )

             rn     = ( MAX( 0.0_wp, v(k,j+1,i) - v_gtrans ) *                & 
                           ABS( ( sk(k,j,i)   - sk(k,j-1,i)            ) /    &
                                ( sk(k,j+1,i) - sk(k,j,i)   + 1E-20_wp )      &
                              ) +                                             & 
                       MIN( 0.0_wp, v(k,j+1,i) - v_gtrans ) *                 &
                           ABS( ( sk(k,j+1,i) - sk(k,j+2,i)            ) /    &
                                ( sk(k,j,i)   - sk(k,j+1,i) + 1E-20_wp )      &
                              )                                               &
                      ) / ABS( v(k,j+1,i) - v_gtrans + 1E-20_wp )     
!   
!--          Reuse k_mm and compute k_mmm for the vertical gradient ratios. 
!--          Note, for vertical advection below the third grid point above 
!--          surface ( or below the model top) rd and rt are set to 0, i.e. 
!--          use of first order scheme is enforced.
             k_mmm  = k - 3 * ibit8

             rd     = ( MAX( 0.0_wp, w(k-1,j,i) ) *                           & 
                           ABS( ( sk(k_mm,j,i) - sk(k_mmm,j,i)           ) /  &
                                ( sk(k-1,j,i)  - sk(k_mm,j,i) + 1E-20_wp )    &
                              ) +                                             & 
                        MIN( 0.0_wp, w(k-1,j,i) ) *                           &
                           ABS( ( sk(k-1,j,i) - sk(k,j,i)            ) /      &
                                ( sk(k_mm,j,i) - sk(k-1,j,i)   + 1E-20_wp )   &
                              )                                               &
                      ) * ibit8 / ABS( w(k-1,j,i) + 1E-20_wp ) 
 
             rt     = ( MAX( 0.0_wp, w(k,j,i) ) *                             & 
                           ABS( ( sk(k,j,i)   - sk(k-1,j,i)            ) /    &
                                ( sk(k+1,j,i) - sk(k,j,i)   + 1E-20_wp )      &
                              ) +                                             & 
                        MIN( 0.0_wp, w(k,j,i) ) *                             &
                           ABS( ( sk(k+1,j,i) - sk(k_pp,j,i)           ) /    &
                                ( sk(k,j,i)   - sk(k+1,j,i) + 1E-20_wp )      &
                              )                                               &
                      ) * ibit8 / ABS( w(k,j,i) + 1E-20_wp )
!
!--           Calculate empirical limiter function (van Albada2 limiter).
              phi_l = MIN( 1.0_wp, ( 2.0_wp * ABS( rl ) ) /                   &
                                         ( rl**2 + 1.0_wp ) ) 
              phi_r = MIN( 1.0_wp, ( 2.0_wp * ABS( rr ) ) /                   &
                                         ( rr**2 + 1.0_wp ) ) 
              phi_s = MIN( 1.0_wp, ( 2.0_wp * ABS( rs ) ) /                   &
                                         ( rs**2 + 1.0_wp ) ) 
              phi_n = MIN( 1.0_wp, ( 2.0_wp * ABS( rn ) ) /                   &
                                         ( rn**2 + 1.0_wp ) ) 
              phi_d = MIN( 1.0_wp, ( 2.0_wp * ABS( rd ) ) /                   &
                                         ( rd**2 + 1.0_wp ) ) 
              phi_t = MIN( 1.0_wp, ( 2.0_wp * ABS( rt ) ) /                   &
                                         ( rt**2 + 1.0_wp ) ) 
!
!--           Calculate the resulting monotone flux. 
              swap_flux_x_local(k,j,tn) = fl_1 - phi_l *                      &
                                        ( fl_1 - swap_flux_x_local(k,j,tn)  )
              flux_r(k)                 = fr_1 - phi_r *                      &
                                        ( fr_1 - flux_r(k)                  )
              swap_flux_y_local(k,tn)   = fs_1 - phi_s *                      &
                                        ( fs_1 - swap_flux_y_local(k,tn)    )
              flux_n(k)                 = fn_1 - phi_n *                      &
                                        ( fn_1 - flux_n(k)                  )
              flux_d                    = fd_1 - phi_d *                      &
                                        ( fd_1 - flux_d                     )
              flux_t(k)                 = ft_1 - phi_t *                      &
                                        ( ft_1 - flux_t(k)                  )
!
!--          Moreover, modify dissipation flux according to the limiter. 
             swap_diss_x_local(k,j,tn) = swap_diss_x_local(k,j,tn) * phi_l
             diss_r(k)                 = diss_r(k)                 * phi_r
             swap_diss_y_local(k,tn)   = swap_diss_y_local(k,tn)   * phi_s
             diss_n(k)                 = diss_n(k)                 * phi_n
             diss_d                    = diss_d                    * phi_d
             diss_t(k)                 = diss_t(k)                 * phi_t

          ENDIF
!
!--       Calculate the divergence of the velocity field. A respective
!--       correction is needed to overcome numerical instabilities caused
!--       by a not sufficient reduction of divergences near topography. 
          div         =   ( u(k,j,i+1) * ( ibit0 + ibit1 + ibit2 )             &
                          - u(k,j,i)   * ( IBITS(wall_flags_0(k,j,i-1),0,1)    &
                                         + IBITS(wall_flags_0(k,j,i-1),1,1)    &
                                         + IBITS(wall_flags_0(k,j,i-1),2,1)    &
                                         )                                     &
                          ) * rho_air(k) * ddx                                 &
                        + ( v(k,j+1,i) * ( ibit3 + ibit4 + ibit5 )             &
                          - v(k,j,i)   * ( IBITS(wall_flags_0(k,j-1,i),3,1)    &
                                         + IBITS(wall_flags_0(k,j-1,i),4,1)    &
                                         + IBITS(wall_flags_0(k,j-1,i),5,1)    &
                                         )                                     &
                          ) * rho_air(k) * ddy                                 &
                        + ( w(k,j,i) * rho_air_zw(k) *                         &
                                         ( ibit6 + ibit7 + ibit8 )             &
                          - w(k-1,j,i) * rho_air_zw(k-1) *                     &
                                         ( IBITS(wall_flags_0(k-1,j,i),6,1)    &
                                         + IBITS(wall_flags_0(k-1,j,i),7,1)    &
                                         + IBITS(wall_flags_0(k-1,j,i),8,1)    &
                                         )                                     &      
                          ) * ddzw(k)


          tend(k,j,i) = tend(k,j,i) - (                                       &
                        ( flux_r(k) + diss_r(k) - swap_flux_x_local(k,j,tn) - &
                          swap_diss_x_local(k,j,tn)            ) * ddx        &
                      + ( flux_n(k) + diss_n(k) - swap_flux_y_local(k,tn)   - &
                          swap_diss_y_local(k,tn)              ) * ddy        &
                      + ( ( flux_t(k) + diss_t(k) ) -                         &
                          ( flux_d    + diss_d    )                           &
                                                    ) * drho_air(k) * ddzw(k) &
                                      ) + sk(k,j,i) * div

          swap_flux_y_local(k,tn)   = flux_n(k)
          swap_diss_y_local(k,tn)   = diss_n(k)
          swap_flux_x_local(k,j,tn) = flux_r(k)
          swap_diss_x_local(k,j,tn) = diss_r(k)
          flux_d                    = flux_t(k)
          diss_d                    = diss_t(k)

       ENDDO
!
!--    Now compute the fluxes and tendency terms for the horizontal and
!--    vertical parts above the top of the highest topography. No degradation
!--    for the horizontal parts, but for the vertical it is stell needed.
       DO  k = nzb_max+1, nzt

          u_comp    = u(k,j,i+1) - u_gtrans
          flux_r(k) = u_comp * (                                              &
                      37.0_wp * ( sk(k,j,i+1) + sk(k,j,i)   )                 &
                    -  8.0_wp * ( sk(k,j,i+2) + sk(k,j,i-1) )                 &
                    +           ( sk(k,j,i+3) + sk(k,j,i-2) ) ) * adv_sca_5
          diss_r(k) = -ABS( u_comp ) * (                                      &
                      10.0_wp * ( sk(k,j,i+1) - sk(k,j,i)   )                 &
                    -  5.0_wp * ( sk(k,j,i+2) - sk(k,j,i-1) )                 &
                    +           ( sk(k,j,i+3) - sk(k,j,i-2) ) ) * adv_sca_5

          v_comp    = v(k,j+1,i) - v_gtrans
          flux_n(k) = v_comp * (                                              &
                      37.0_wp * ( sk(k,j+1,i) + sk(k,j,i)   )                 &
                    -  8.0_wp * ( sk(k,j+2,i) + sk(k,j-1,i) )                 &
                    +           ( sk(k,j+3,i) + sk(k,j-2,i) ) ) * adv_sca_5
          diss_n(k) = -ABS( v_comp ) * (                                      &
                      10.0_wp * ( sk(k,j+1,i) - sk(k,j,i)   )                 &
                    -  5.0_wp * ( sk(k,j+2,i) - sk(k,j-1,i) )                 &
                    +           ( sk(k,j+3,i) - sk(k,j-2,i) ) ) * adv_sca_5
!
!--       k index has to be modified near bottom and top, else array
!--       subscripts will be exceeded.
          ibit8 = IBITS(wall_flags_0(k,j,i),8,1)
          ibit7 = IBITS(wall_flags_0(k,j,i),7,1)
          ibit6 = IBITS(wall_flags_0(k,j,i),6,1)

          k_ppp = k + 3 * ibit8
          k_pp  = k + 2 * ( 1 - ibit6  )
          k_mm  = k - 2 * ibit8

          flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                            &
                    ( 37.0_wp * ibit8 * adv_sca_5                             &
                 +     7.0_wp * ibit7 * adv_sca_3                             &
                 +              ibit6 * adv_sca_1                             &
                    ) *                                                       &
                             ( sk(k+1,j,i)  + sk(k,j,i)   )                   &
              -     (  8.0_wp * ibit8 * adv_sca_5                             &
                  +              ibit7 * adv_sca_3                            &
                    ) *                                                       &
                             ( sk(k_pp,j,i) + sk(k-1,j,i) )                   &
              +     (           ibit8 * adv_sca_5                             &
                    ) *     ( sk(k_ppp,j,i)+ sk(k_mm,j,i) )                   &
                                 )

          diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                    &
                    ( 10.0_wp * ibit8 * adv_sca_5                             &
                 +     3.0_wp * ibit7 * adv_sca_3                             &
                 +              ibit6 * adv_sca_1                             &
                    ) *                                                       &
                             ( sk(k+1,j,i)   - sk(k,j,i)    )                 &
              -     (  5.0_wp * ibit8 * adv_sca_5                             &
                 +              ibit7 * adv_sca_3                             &
                    ) *                                                       &
                             ( sk(k_pp,j,i)  - sk(k-1,j,i)  )                 &
              +     (           ibit8 * adv_sca_5                             &
                    ) *                                                       &
                             ( sk(k_ppp,j,i) - sk(k_mm,j,i) )                 &
                                         )


!
!--       Apply monotonic adjustment.
          IF ( monotonic_adjustment )  THEN
!
!--          At first, calculate first order fluxes.
             u_comp = u(k,j,i) - u_gtrans
             fl_1   =  ( u_comp   * ( sk(k,j,i) + sk(k,j,i-1) )               &
                   -ABS( u_comp ) * ( sk(k,j,i) - sk(k,j,i-1) )               &
                       ) * adv_sca_1 

             u_comp = u(k,j,i+1) - u_gtrans
             fr_1   =  ( u_comp   * ( sk(k,j,i+1) + sk(k,j,i) )               &
                   -ABS( u_comp ) * ( sk(k,j,i+1) - sk(k,j,i) )               &
                       ) * adv_sca_1 

             v_comp = v(k,j,i) - v_gtrans
             fs_1   =  ( v_comp   * ( sk(k,j,i) + sk(k,j-1,i) )               &
                   -ABS( v_comp ) * ( sk(k,j,i) - sk(k,j-1,i) )               &
                       ) * adv_sca_1 

             v_comp = v(k,j+1,i) - v_gtrans
             fn_1   =  ( v_comp   * ( sk(k,j+1,i) + sk(k,j,i) )               &
                   -ABS( v_comp ) * ( sk(k,j+1,i) - sk(k,j,i) )               &
                       ) * adv_sca_1 

             fd_1   =  ( w(k-1,j,i)   * ( sk(k,j,i) + sk(k-1,j,i) )           &
                   -ABS( w(k-1,j,i) ) * ( sk(k,j,i) - sk(k-1,j,i) )           &
                       ) * adv_sca_1  * rho_air_zw(k)

             ft_1   =  ( w(k,j,i)   * ( sk(k+1,j,i) + sk(k,j,i) )             &
                   -ABS( w(k,j,i) ) * ( sk(k+1,j,i) - sk(k,j,i) )             &
                       ) * adv_sca_1  * rho_air_zw(k)
!
!--          Calculate ratio of upwind gradients. Note, Min/Max is just to 
!--          avoid if statements. 
             rl     = ( MAX( 0.0_wp, u(k,j,i) - u_gtrans ) *                  & 
                           ABS( ( sk(k,j,i-1) - sk(k,j,i-2)            ) /    &
                                ( sk(k,j,i)   - sk(k,j,i-1) + 1E-20_wp )      &
                              ) +                                             & 
                        MIN( 0.0_wp, u(k,j,i) - u_gtrans ) *                  &
                           ABS( ( sk(k,j,i)   - sk(k,j,i+1)            ) /    &
                                ( sk(k,j,i-1) - sk(k,j,i)   + 1E-20_wp )      &
                              )                                               &
                      ) / ABS( u(k,j,i) - u_gtrans + 1E-20_wp )

             rr     = ( MAX( 0.0_wp, u(k,j,i+1) - u_gtrans ) *                & 
                           ABS( ( sk(k,j,i)   - sk(k,j,i-1)            ) /    &
                                ( sk(k,j,i+1) - sk(k,j,i)   + 1E-20_wp )      &
                              ) +                                             & 
                        MIN( 0.0_wp, u(k,j,i+1) - u_gtrans ) *                &
                           ABS( ( sk(k,j,i+1) - sk(k,j,i+2)            ) /    &
                                ( sk(k,j,i)   - sk(k,j,i+1) + 1E-20_wp )      &
                              )                                               &
                      ) / ABS( u(k,j,i+1) - u_gtrans + 1E-20_wp )

             rs     = ( MAX( 0.0_wp, v(k,j,i) - v_gtrans ) *                  & 
                           ABS( ( sk(k,j-1,i) - sk(k,j-2,i)            ) /    &
                                ( sk(k,j,i)   - sk(k,j-1,i) + 1E-20_wp )      &
                              ) +                                             & 
                        MIN( 0.0_wp, v(k,j,i) - v_gtrans ) *                  &
                           ABS( ( sk(k,j,i)   - sk(k,j+1,i)            ) /    &
                                ( sk(k,j-1,i) - sk(k,j,i)   + 1E-20_wp )      &
                              )                                               &
                      ) / ABS( v(k,j,i) - v_gtrans + 1E-20_wp )

             rn     = ( MAX( 0.0_wp, v(k,j+1,i) - v_gtrans ) *                & 
                           ABS( ( sk(k,j,i)   - sk(k,j-1,i)            ) /    &
                                ( sk(k,j+1,i) - sk(k,j,i)   + 1E-20_wp )      &
                              ) +                                             & 
                       MIN( 0.0_wp, v(k,j+1,i) - v_gtrans ) *                 &
                           ABS( ( sk(k,j+1,i) - sk(k,j+2,i)            ) /    &
                                ( sk(k,j,i)   - sk(k,j+1,i) + 1E-20_wp )      &
                              )                                               &
                      ) / ABS( v(k,j+1,i) - v_gtrans + 1E-20_wp )     
!
!--          Reuse k_mm and compute k_mmm for the vertical gradient ratios. 
!--          Note, for vertical advection below the third grid point above 
!--          surface ( or below the model top) rd and rt are set to 0, i.e. 
!--          use of first order scheme is enforced.
             k_mmm  = k - 3 * ibit8

             rd     = ( MAX( 0.0_wp, w(k-1,j,i) ) *                           & 
                           ABS( ( sk(k_mm,j,i) - sk(k_mmm,j,i)           ) /  &
                                ( sk(k-1,j,i)  - sk(k_mm,j,i) + 1E-20_wp )    &
                              ) +                                             & 
                        MIN( 0.0_wp, w(k-1,j,i) ) *                           &
                           ABS( ( sk(k-1,j,i) - sk(k,j,i)            ) /      &
                                ( sk(k_mm,j,i) - sk(k-1,j,i)   + 1E-20_wp )   &
                              )                                               &
                      ) * ibit8 / ABS( w(k-1,j,i) + 1E-20_wp ) 
 
             rt     = ( MAX( 0.0_wp, w(k,j,i) ) *                             & 
                           ABS( ( sk(k,j,i)   - sk(k-1,j,i)            ) /    &
                                ( sk(k+1,j,i) - sk(k,j,i)   + 1E-20_wp )      &
                              ) +                                             & 
                        MIN( 0.0_wp, w(k,j,i) ) *                             &
                           ABS( ( sk(k+1,j,i) - sk(k_pp,j,i)           ) /    &
                                ( sk(k,j,i)   - sk(k+1,j,i) + 1E-20_wp )      &
                              )                                               &
                      ) * ibit8 / ABS( w(k,j,i) + 1E-20_wp ) 
!
!--           Calculate empirical limiter function (van Albada2 limiter).
              phi_l = MIN( 1.0_wp, ( 2.0_wp * ABS( rl ) ) /                   &
                                         ( rl**2 + 1.0_wp ) ) 
              phi_r = MIN( 1.0_wp, ( 2.0_wp * ABS( rr ) ) /                   &
                                         ( rr**2 + 1.0_wp ) ) 
              phi_s = MIN( 1.0_wp, ( 2.0_wp * ABS( rs ) ) /                   &
                                         ( rs**2 + 1.0_wp ) ) 
              phi_n = MIN( 1.0_wp, ( 2.0_wp * ABS( rn ) ) /                   &
                                         ( rn**2 + 1.0_wp ) ) 
              phi_d = MIN( 1.0_wp, ( 2.0_wp * ABS( rd ) ) /                   &
                                         ( rd**2 + 1.0_wp ) ) 
              phi_t = MIN( 1.0_wp, ( 2.0_wp * ABS( rt ) ) /                   &
                                         ( rt**2 + 1.0_wp ) ) 
!
!--           Calculate the resulting monotone flux. 
              swap_flux_x_local(k,j,tn) = fl_1 - phi_l *                      &
                                        ( fl_1 - swap_flux_x_local(k,j,tn)  )
              flux_r(k)                 = fr_1 - phi_r *                      &
                                        ( fr_1 - flux_r(k)                  )
              swap_flux_y_local(k,tn)   = fs_1 - phi_s *                      &
                                        ( fs_1 - swap_flux_y_local(k,tn)    )
              flux_n(k)                 = fn_1 - phi_n *                      &
                                        ( fn_1 - flux_n(k)                  )
              flux_d                    = fd_1 - phi_d *                      &
                                        ( fd_1 - flux_d                     )
              flux_t(k)                 = ft_1 - phi_t *                      &
                                        ( ft_1 - flux_t(k)                  )
!
!--          Moreover, modify dissipation flux according to the limiter. 
             swap_diss_x_local(k,j,tn) = swap_diss_x_local(k,j,tn) * phi_l
             diss_r(k)                 = diss_r(k)                 * phi_r
             swap_diss_y_local(k,tn)   = swap_diss_y_local(k,tn)   * phi_s
             diss_n(k)                 = diss_n(k)                 * phi_n
             diss_d                    = diss_d                    * phi_d
             diss_t(k)                 = diss_t(k)                 * phi_t

          ENDIF
!
!--       Calculate the divergence of the velocity field. A respective
!--       correction is needed to overcome numerical instabilities introduced
!--       by a not sufficient reduction of divergences near topography.
          div         =   ( u(k,j,i+1) - u(k,j,i)   ) * rho_air(k) * ddx      &
                        + ( v(k,j+1,i) - v(k,j,i)   ) * rho_air(k) * ddy      &
                        + ( w(k,j,i)   * rho_air_zw(k) -                      &
                            w(k-1,j,i) * rho_air_zw(k-1) ) * ddzw(k)

          tend(k,j,i) = tend(k,j,i) - (                                       &
                        ( flux_r(k) + diss_r(k) - swap_flux_x_local(k,j,tn) - &
                          swap_diss_x_local(k,j,tn)            ) * ddx        &
                      + ( flux_n(k) + diss_n(k) - swap_flux_y_local(k,tn)   - &
                          swap_diss_y_local(k,tn)              ) * ddy        &
                      + ( ( flux_t(k) + diss_t(k) ) -                         &
                          ( flux_d    + diss_d    )                           &
                                                    ) * drho_air(k) * ddzw(k) &
                                      ) + sk(k,j,i) * div


          swap_flux_y_local(k,tn)   = flux_n(k)
          swap_diss_y_local(k,tn)   = diss_n(k)
          swap_flux_x_local(k,j,tn) = flux_r(k)
          swap_diss_x_local(k,j,tn) = diss_r(k)
          flux_d                    = flux_t(k)
          diss_d                    = diss_t(k)

       ENDDO

!
!--    Evaluation of statistics.
       SELECT CASE ( sk_char )

          CASE ( 'pt' )

             DO  k = nzb, nzt
                sums_wspts_ws_l(k,tn) = sums_wspts_ws_l(k,tn) +                &
                    ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                    + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS( w(k,j,i) - hom(k,1,3,0)            )  &
                    ) * weight_substep(intermediate_timestep_count)
             ENDDO
            
          CASE ( 'sa' )

             DO  k = nzb, nzt
                sums_wssas_ws_l(k,tn) = sums_wssas_ws_l(k,tn) +                &
                    ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                    + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS( w(k,j,i) - hom(k,1,3,0)            )  &
                    ) * weight_substep(intermediate_timestep_count)
             ENDDO
            
          CASE ( 'q' )

             DO  k = nzb, nzt
                sums_wsqs_ws_l(k,tn)  = sums_wsqs_ws_l(k,tn) +                 &
                    ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                    + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS( w(k,j,i) - hom(k,1,3,0)            )  &
                    ) * weight_substep(intermediate_timestep_count)
             ENDDO

          CASE ( 'qr' )

             DO  k = nzb, nzt
                sums_wsqrs_ws_l(k,tn)  = sums_wsqrs_ws_l(k,tn) +               &
                    ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                    + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS( w(k,j,i) - hom(k,1,3,0)            )  &
                    ) * weight_substep(intermediate_timestep_count)
             ENDDO

          CASE ( 'nr' )

             DO  k = nzb, nzt
                sums_wsnrs_ws_l(k,tn)  = sums_wsnrs_ws_l(k,tn) +               &
                    ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                    + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS( w(k,j,i) - hom(k,1,3,0)            )  &
                    ) * weight_substep(intermediate_timestep_count)
             ENDDO
             
          CASE ( 's' )
          
             DO  k = nzb, nzt
                sums_wsss_ws_l(k,tn)  = sums_wsss_ws_l(k,tn) +                 &
                    ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                    + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS( w(k,j,i) - hom(k,1,3,0)            )  &
                    ) * weight_substep(intermediate_timestep_count)
             ENDDO
#ifdef KPP_CHEM
          CASE ( 'kc' )
          !kk Has to be implemented for kpp chemistry
#endif

         END SELECT
         
    END SUBROUTINE advec_s_ws_ij




!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of u-component - Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_u_ws_ij( i, j, i_omp, tn )

       USE arrays_3d,                                                         &
           ONLY:  ddzw, diss_l_u, diss_s_u, flux_l_u, flux_s_u, tend, u, v, w,&
                  drho_air, rho_air, rho_air_zw

       USE constants,                                                         &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                &
           ONLY:  intermediate_timestep_count, u_gtrans, v_gtrans

       USE grid_variables,                                                    &
           ONLY:  ddx, ddy

       USE indices,                                                           &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzb_max, nzt, wall_flags_0

       USE kinds

       USE statistics,                                                        &
           ONLY:  hom, sums_us2_ws_l, sums_wsus_ws_l, weight_substep

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit9  !<
       INTEGER(iwp) ::  ibit10 !<
       INTEGER(iwp) ::  ibit11 !<
       INTEGER(iwp) ::  ibit12 !<
       INTEGER(iwp) ::  ibit13 !<
       INTEGER(iwp) ::  ibit14 !<
       INTEGER(iwp) ::  ibit15 !<
       INTEGER(iwp) ::  ibit16 !<
       INTEGER(iwp) ::  ibit17 !<
       INTEGER(iwp) ::  i_omp  !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn     !<
       
       REAL(wp)    ::  diss_d   !<
       REAL(wp)    ::  div      !<
       REAL(wp)    ::  flux_d   !<
       REAL(wp)    ::  gu       !<
       REAL(wp)    ::  gv       !<
       REAL(wp)    ::  u_comp_l !<
       REAL(wp)    ::  v_comp   !<
       REAL(wp)    ::  w_comp   !<
       
       REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_n !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_r !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_t !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_n !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_r !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_t !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  u_comp !<

       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans
!
!--    Compute southside fluxes for the respective boundary of PE
       IF ( j == nys  )  THEN
       
          DO  k = nzb+1, nzb_max

             ibit14 = IBITS(wall_flags_0(k,j-1,i),14,1)
             ibit13 = IBITS(wall_flags_0(k,j-1,i),13,1)
             ibit12 = IBITS(wall_flags_0(k,j-1,i),12,1)

             v_comp      = v(k,j,i) + v(k,j,i-1) - gv
             flux_s_u(k,tn) = v_comp * (                                      &
                            ( 37.0_wp * ibit14 * adv_mom_5                    &
                         +     7.0_wp * ibit13 * adv_mom_3                    &
                         +              ibit12 * adv_mom_1                    &
                            ) *                                               &
                                        ( u(k,j,i)   + u(k,j-1,i) )           &
                     -      (  8.0_wp * ibit14 * adv_mom_5                    &
                         +              ibit13 * adv_mom_3                    &
                            ) *                                               &
                                        ( u(k,j+1,i) + u(k,j-2,i) )           &
                     +      (           ibit14 * adv_mom_5                    &
                            ) *                                               &
                                        ( u(k,j+2,i) + u(k,j-3,i) )           &
                                       )

             diss_s_u(k,tn) = - ABS ( v_comp ) * (                            &
                            ( 10.0_wp * ibit14 * adv_mom_5                    &
                         +     3.0_wp * ibit13 * adv_mom_3                    &
                         +              ibit12 * adv_mom_1                    &
                            ) *                                               &
                                        ( u(k,j,i)   - u(k,j-1,i) )           &
                     -      (  5.0_wp * ibit14 * adv_mom_5                    &
                         +              ibit13 * adv_mom_3                    &
                            ) *                                               &
                                        ( u(k,j+1,i) - u(k,j-2,i) )           &
                     +      (           ibit14 * adv_mom_5                    &
                            ) *                                               &
                                        ( u(k,j+2,i) - u(k,j-3,i) )           &
                                                 )

          ENDDO

          DO  k = nzb_max+1, nzt

             v_comp         = v(k,j,i) + v(k,j,i-1) - gv
             flux_s_u(k,tn) = v_comp * (                                      &
                           37.0_wp * ( u(k,j,i) + u(k,j-1,i)   )              &
                         -  8.0_wp * ( u(k,j+1,i) + u(k,j-2,i) )              &
                         +           ( u(k,j+2,i) + u(k,j-3,i) ) ) * adv_mom_5
             diss_s_u(k,tn) = - ABS(v_comp) * (                               &
                           10.0_wp * ( u(k,j,i) - u(k,j-1,i)   )              &
                         -  5.0_wp * ( u(k,j+1,i) - u(k,j-2,i) )              &
                         +           ( u(k,j+2,i) - u(k,j-3,i) ) ) * adv_mom_5

          ENDDO
          
       ENDIF
!
!--    Compute leftside fluxes for the respective boundary of PE
       IF ( i == i_omp )  THEN
       
          DO  k = nzb+1, nzb_max

             ibit11 = IBITS(wall_flags_0(k,j,i-1),11,1)
             ibit10 = IBITS(wall_flags_0(k,j,i-1),10,1)
             ibit9  = IBITS(wall_flags_0(k,j,i-1),9,1)

             u_comp_l         = u(k,j,i) + u(k,j,i-1) - gu
             flux_l_u(k,j,tn) = u_comp_l * (                                  &
                              ( 37.0_wp * ibit11 * adv_mom_5                  &
                           +     7.0_wp * ibit10 * adv_mom_3                  &
                           +              ibit9  * adv_mom_1                  &
                              ) *                                             &
                                          ( u(k,j,i)   + u(k,j,i-1) )         &
                       -      (  8.0_wp * ibit11 * adv_mom_5                  &
                           +              ibit10 * adv_mom_3                  &
                              ) *                                             &
                                          ( u(k,j,i+1) + u(k,j,i-2) )         &
                       +      (           ibit11 * adv_mom_5                  &
                              ) *                                             &
                                          ( u(k,j,i+2) + u(k,j,i-3) )         &
                                           )

              diss_l_u(k,j,tn) = - ABS( u_comp_l ) * (                        &
                              ( 10.0_wp * ibit11 * adv_mom_5                  &
                           +     3.0_wp * ibit10 * adv_mom_3                  &
                           +              ibit9  * adv_mom_1                  &
                              ) *                                             &
                                        ( u(k,j,i)   - u(k,j,i-1) )           &
                       -      (  5.0_wp * ibit11 * adv_mom_5                  &
                           +              ibit10 * adv_mom_3                  &
                              ) *                                             &
                                        ( u(k,j,i+1) - u(k,j,i-2) )           &
                       +      (           ibit11 * adv_mom_5                  &
                              ) *                                             &
                                        ( u(k,j,i+2) - u(k,j,i-3) )           &
                                                     )

          ENDDO

          DO  k = nzb_max+1, nzt

             u_comp_l         = u(k,j,i) + u(k,j,i-1) - gu
             flux_l_u(k,j,tn) = u_comp_l * (                                   &
                             37.0_wp * ( u(k,j,i) + u(k,j,i-1)   )             &
                           -  8.0_wp * ( u(k,j,i+1) + u(k,j,i-2) )             &
                           +           ( u(k,j,i+2) + u(k,j,i-3) ) ) * adv_mom_5
             diss_l_u(k,j,tn) = - ABS(u_comp_l) * (                            &
                             10.0_wp * ( u(k,j,i) - u(k,j,i-1)   )             &
                           -  5.0_wp * ( u(k,j,i+1) - u(k,j,i-2) )             &
                           +           ( u(k,j,i+2) - u(k,j,i-3) ) ) * adv_mom_5

          ENDDO
          
       ENDIF

       flux_t(0) = 0.0_wp
       diss_t(0) = 0.0_wp
       flux_d    = 0.0_wp
       diss_d    = 0.0_wp
!
!--    Now compute the fluxes tendency terms for the horizontal and
!--    vertical parts.
       DO  k = nzb+1, nzb_max

          ibit11 = IBITS(wall_flags_0(k,j,i),11,1)
          ibit10 = IBITS(wall_flags_0(k,j,i),10,1)
          ibit9  = IBITS(wall_flags_0(k,j,i),9,1)

          u_comp(k) = u(k,j,i+1) + u(k,j,i)
          flux_r(k) = ( u_comp(k) - gu ) * (                                  &
                     ( 37.0_wp * ibit11 * adv_mom_5                           &
                  +     7.0_wp * ibit10 * adv_mom_3                           &
                  +              ibit9  * adv_mom_1                           &
                     ) *                                                      &
                                    ( u(k,j,i+1) + u(k,j,i)   )               &
              -      (  8.0_wp * ibit11 * adv_mom_5                           &
                  +              ibit10 * adv_mom_3                           & 
                     ) *                                                      &
                                    ( u(k,j,i+2) + u(k,j,i-1) )               &
              +      (           ibit11 * adv_mom_5                           &
                     ) *                                                      &
                                    ( u(k,j,i+3) + u(k,j,i-2) )               &
                                           )

          diss_r(k) = - ABS( u_comp(k) - gu ) * (                             &
                     ( 10.0_wp * ibit11 * adv_mom_5                           &
                  +     3.0_wp * ibit10 * adv_mom_3                           &
                  +              ibit9  * adv_mom_1                           &
                     ) *                                                      &
                                    ( u(k,j,i+1) - u(k,j,i)   )               &
              -      (  5.0_wp * ibit11 * adv_mom_5                           &
                  +              ibit10 * adv_mom_3                           &
                     ) *                                                      &
                                    ( u(k,j,i+2) - u(k,j,i-1) )               &
              +      (           ibit11 * adv_mom_5                           &
                     ) *                                                      &
                                    ( u(k,j,i+3) - u(k,j,i-2) )               &
                                                )

          ibit14 = IBITS(wall_flags_0(k,j,i),14,1)
          ibit13 = IBITS(wall_flags_0(k,j,i),13,1)
          ibit12 = IBITS(wall_flags_0(k,j,i),12,1)

          v_comp    = v(k,j+1,i) + v(k,j+1,i-1) - gv
          flux_n(k) = v_comp * (                                              &
                     ( 37.0_wp * ibit14 * adv_mom_5                           &
                  +     7.0_wp * ibit13 * adv_mom_3                           &
                  +              ibit12 * adv_mom_1                           &
                     ) *                                                      &
                                    ( u(k,j+1,i) + u(k,j,i)   )               &
              -      (  8.0_wp * ibit14 * adv_mom_5                           &
                  +              ibit13 * adv_mom_3                           &
                     ) *                                                      &
                                    ( u(k,j+2,i) + u(k,j-1,i) )               &
              +      (           ibit14 * adv_mom_5                           &
                     ) *                                                      &
                                    ( u(k,j+3,i) + u(k,j-2,i) )               &
                               )

          diss_n(k) = - ABS ( v_comp ) * (                                    &
                     ( 10.0_wp * ibit14 * adv_mom_5                           &
                  +     3.0_wp * ibit13 * adv_mom_3                           &
                  +              ibit12 * adv_mom_1                           &
                     ) *                                                      &
                                    ( u(k,j+1,i) - u(k,j,i)   )               &
              -      (  5.0_wp * ibit14 * adv_mom_5                           &
                  +              ibit13 * adv_mom_3                           &
                     ) *                                                      &
                                    ( u(k,j+2,i) - u(k,j-1,i) )               &
              +      (           ibit14 * adv_mom_5                           &
                     ) *                                                      &
                                    ( u(k,j+3,i) - u(k,j-2,i) )               &
                                         )
!
!--       k index has to be modified near bottom and top, else array
!--       subscripts will be exceeded.
          ibit17 = IBITS(wall_flags_0(k,j,i),17,1)
          ibit16 = IBITS(wall_flags_0(k,j,i),16,1)
          ibit15 = IBITS(wall_flags_0(k,j,i),15,1)

          k_ppp = k + 3 * ibit17
          k_pp  = k + 2 * ( 1 - ibit15 )
          k_mm  = k - 2 * ibit17

          w_comp    = w(k,j,i) + w(k,j,i-1)
          flux_t(k) = w_comp * rho_air_zw(k) * (                              &
                     ( 37.0_wp * ibit17 * adv_mom_5                           &
                  +     7.0_wp * ibit16 * adv_mom_3                           &
                  +              ibit15 * adv_mom_1                           &
                     ) *                                                      &
                                ( u(k+1,j,i)  + u(k,j,i)     )                &
              -      (  8.0_wp * ibit17 * adv_mom_5                           &
                  +              ibit16 * adv_mom_3                           &
                     ) *                                                      &
                                ( u(k_pp,j,i) + u(k-1,j,i)   )                &
              +      (           ibit17 * adv_mom_5                           &
                     ) *                                                      &
                                ( u(k_ppp,j,i) + u(k_mm,j,i) )                &
                                 )

          diss_t(k) = - ABS( w_comp ) * rho_air_zw(k) * (                     &
                     ( 10.0_wp * ibit17 * adv_mom_5                           &
                  +     3.0_wp * ibit16 * adv_mom_3                           &
                  +              ibit15 * adv_mom_1                           &
                     ) *                                                      &
                                ( u(k+1,j,i)   - u(k,j,i)    )                &
              -      (  5.0_wp * ibit17 * adv_mom_5                           &
                  +              ibit16 * adv_mom_3                           &
                     ) *                                                      &
                                ( u(k_pp,j,i)  - u(k-1,j,i)  )                &
              +      (           ibit17 * adv_mom_5                           &
                     ) *                                                      &
                                ( u(k_ppp,j,i) - u(k_mm,j,i) )                &
                                         )
!
!--       Calculate the divergence of the velocity field. A respective
!--       correction is needed to overcome numerical instabilities introduced
!--       by a not sufficient reduction of divergences near topography. 
          div = ( ( u_comp(k)       * ( ibit9 + ibit10 + ibit11 )             &
                - ( u(k,j,i)   + u(k,j,i-1)   )                               &
                                    * ( IBITS(wall_flags_0(k,j,i-1),9,1)      &
                                      + IBITS(wall_flags_0(k,j,i-1),10,1)     &
                                      + IBITS(wall_flags_0(k,j,i-1),11,1)     &
                                      )                                       &
                  ) * rho_air(k) * ddx                                        &
               +  ( ( v_comp + gv ) * ( ibit12 + ibit13 + ibit14 )            &
                  - ( v(k,j,i)   + v(k,j,i-1 )  )                             &
                                    * ( IBITS(wall_flags_0(k,j-1,i),12,1)     &
                                      + IBITS(wall_flags_0(k,j-1,i),13,1)     &
                                      + IBITS(wall_flags_0(k,j-1,i),14,1)     &
                                      )                                       &
                  ) * rho_air(k) * ddy                                        &
               +  ( w_comp * rho_air_zw(k) * ( ibit15 + ibit16 + ibit17 )     &
                - ( w(k-1,j,i) + w(k-1,j,i-1) ) * rho_air_zw(k-1)             &
                                    * ( IBITS(wall_flags_0(k-1,j,i),15,1)     &
                                      + IBITS(wall_flags_0(k-1,j,i),16,1)     &
                                      + IBITS(wall_flags_0(k-1,j,i),17,1)     &
                                      )                                       &  
                  ) * ddzw(k)   &
                ) * 0.5_wp


          tend(k,j,i) = tend(k,j,i) - (                                       &
                            ( flux_r(k) + diss_r(k)                           &
                          -   flux_l_u(k,j,tn) - diss_l_u(k,j,tn) ) * ddx     &
                          + ( flux_n(k) + diss_n(k)                           &
                          -   flux_s_u(k,tn) - diss_s_u(k,tn)     ) * ddy     &
                          + ( ( flux_t(k) + diss_t(k) )                       &
                          -   ( flux_d    + diss_d )                          &
                                                    ) * drho_air(k) * ddzw(k) &
                                       ) + div * u(k,j,i)

           flux_l_u(k,j,tn) = flux_r(k)
           diss_l_u(k,j,tn) = diss_r(k)
           flux_s_u(k,tn)   = flux_n(k)
           diss_s_u(k,tn)   = diss_n(k)
           flux_d           = flux_t(k)
           diss_d           = diss_t(k)
!
!--        Statistical Evaluation of u'u'. The factor has to be applied for
!--        right evaluation when gallilei_trans = .T. .
           sums_us2_ws_l(k,tn) = sums_us2_ws_l(k,tn)                           &
                + ( flux_r(k)                                                  &
                    * ( u_comp(k) - 2.0_wp * hom(k,1,1,0)                   )  &
                    / ( u_comp(k) - gu + SIGN( 1.0E-20_wp, u_comp(k) - gu ) )  &
                  + diss_r(k)                                                  &
                    *   ABS( u_comp(k) - 2.0_wp * hom(k,1,1,0)              )  &
                    / ( ABS( u_comp(k) - gu ) + 1.0E-20_wp                  )  &
                  ) *   weight_substep(intermediate_timestep_count)
!
!--        Statistical Evaluation of w'u'.
           sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn)                         &
                + ( flux_t(k)                                                  &
                    * ( w_comp - 2.0_wp * hom(k,1,3,0)                   )     &
                    / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              )     &
                  + diss_t(k)                                                  &
                    *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              )     &
                    / ( ABS( w_comp ) + 1.0E-20_wp                       )     &
                  ) *   weight_substep(intermediate_timestep_count)
       ENDDO

       DO  k = nzb_max+1, nzt

          u_comp(k) = u(k,j,i+1) + u(k,j,i)
          flux_r(k) = ( u_comp(k) - gu ) * (                                  &
                         37.0_wp * ( u(k,j,i+1) + u(k,j,i)   )                &
                       -  8.0_wp * ( u(k,j,i+2) + u(k,j,i-1) )                &
                       +           ( u(k,j,i+3) + u(k,j,i-2) ) ) * adv_mom_5
          diss_r(k) = - ABS( u_comp(k) - gu ) * (                             &
                         10.0_wp * ( u(k,j,i+1) - u(k,j,i)   )                &
                       -  5.0_wp * ( u(k,j,i+2) - u(k,j,i-1) )                &
                       +           ( u(k,j,i+3) - u(k,j,i-2) ) ) * adv_mom_5

          v_comp    = v(k,j+1,i) + v(k,j+1,i-1) - gv
          flux_n(k) = v_comp * (                                              &
                         37.0_wp * ( u(k,j+1,i) + u(k,j,i)   )                &
                       -  8.0_wp * ( u(k,j+2,i) + u(k,j-1,i) )                &
                       +           ( u(k,j+3,i) + u(k,j-2,i) ) ) * adv_mom_5
          diss_n(k) = - ABS( v_comp ) * (                                     &
                         10.0_wp * ( u(k,j+1,i) - u(k,j,i)   )                &
                       -  5.0_wp * ( u(k,j+2,i) - u(k,j-1,i) )                &
                       +           ( u(k,j+3,i) - u(k,j-2,i) ) ) * adv_mom_5
!
!--       k index has to be modified near bottom and top, else array
!--       subscripts will be exceeded.
          ibit17 = IBITS(wall_flags_0(k,j,i),17,1)
          ibit16 = IBITS(wall_flags_0(k,j,i),16,1)
          ibit15 = IBITS(wall_flags_0(k,j,i),15,1)

          k_ppp = k + 3 * ibit17
          k_pp  = k + 2 * ( 1 - ibit15 )
          k_mm  = k - 2 * ibit17

          w_comp    = w(k,j,i) + w(k,j,i-1)
          flux_t(k) = w_comp * rho_air_zw(k) * (                              &
                     ( 37.0_wp * ibit17 * adv_mom_5                           &
                  +     7.0_wp * ibit16 * adv_mom_3                           &
                  +              ibit15 * adv_mom_1                           &
                     ) *                                                      &
                                ( u(k+1,j,i)  + u(k,j,i)     )                &
              -      (  8.0_wp * ibit17 * adv_mom_5                           &
                  +              ibit16 * adv_mom_3                           &
                     ) *                                                      &
                                ( u(k_pp,j,i) + u(k-1,j,i)   )                &
              +      (           ibit17 * adv_mom_5                           &
                     ) *                                                      &
                                ( u(k_ppp,j,i) + u(k_mm,j,i) )                &
                                 )

          diss_t(k) = - ABS( w_comp ) * rho_air_zw(k) * (                     &
                     ( 10.0_wp * ibit17 * adv_mom_5                           &
                  +     3.0_wp * ibit16 * adv_mom_3                           &
                  +              ibit15 * adv_mom_1                           &
                     ) *                                                      &
                                ( u(k+1,j,i)   - u(k,j,i)    )                &
              -      (  5.0_wp * ibit17 * adv_mom_5                           &
                  +              ibit16 * adv_mom_3                           &
                     ) *                                                      &
                                ( u(k_pp,j,i)  - u(k-1,j,i)  )                &
              +      (           ibit17 * adv_mom_5                           &
                     ) *                                                      &
                                ( u(k_ppp,j,i) - u(k_mm,j,i) )                &
                                         )
!
!--       Calculate the divergence of the velocity field. A respective
!--       correction is needed to overcome numerical instabilities introduced
!--       by a not sufficient reduction of divergences near topography.
          div = ( ( u_comp(k)   - ( u(k,j,i)   + u(k,j,i-1)   ) ) * ddx       &
                                                                  * rho_air(k)&
               +  ( v_comp + gv - ( v(k,j,i)   + v(k,j,i-1 )  ) ) * ddy       &
                                                                  * rho_air(k)&
               +  (   w_comp                      * rho_air_zw(k)   -         &
                    ( w(k-1,j,i) + w(k-1,j,i-1) ) * rho_air_zw(k-1)           &
                  ) * ddzw(k)                                                 &
                ) * 0.5_wp

          tend(k,j,i) = tend(k,j,i) - (                                       &
                            ( flux_r(k) + diss_r(k)                           &
                          -   flux_l_u(k,j,tn) - diss_l_u(k,j,tn) ) * ddx     &
                          + ( flux_n(k) + diss_n(k)                           &
                          -   flux_s_u(k,tn) - diss_s_u(k,tn)     ) * ddy     &
                          + ( ( flux_t(k) + diss_t(k) )                       &
                          -   ( flux_d    + diss_d    )                       &
                                                    ) * drho_air(k) * ddzw(k) &
                                       ) + div * u(k,j,i)

           flux_l_u(k,j,tn) = flux_r(k)
           diss_l_u(k,j,tn) = diss_r(k)
           flux_s_u(k,tn)   = flux_n(k)
           diss_s_u(k,tn)   = diss_n(k)
           flux_d           = flux_t(k)
           diss_d           = diss_t(k)
!
!--        Statistical Evaluation of u'u'. The factor has to be applied for
!--        right evaluation when gallilei_trans = .T. .
           sums_us2_ws_l(k,tn) = sums_us2_ws_l(k,tn)                           &
                + ( flux_r(k)                                                  &
                    * ( u_comp(k) - 2.0_wp * hom(k,1,1,0)                   )  &
                    / ( u_comp(k) - gu + SIGN( 1.0E-20_wp, u_comp(k) - gu ) )  &
                  + diss_r(k)                                                  &
                    *   ABS( u_comp(k) - 2.0_wp * hom(k,1,1,0)              )  &
                    / ( ABS( u_comp(k) - gu ) + 1.0E-20_wp                  )  &
                  ) *   weight_substep(intermediate_timestep_count)
!
!--        Statistical Evaluation of w'u'.
           sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn)                         &
                + ( flux_t(k)                                                  &
                    * ( w_comp - 2.0_wp * hom(k,1,3,0)                   )     &
                    / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              )     &
                  + diss_t(k)                                                  &
                    *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              )     &
                    / ( ABS( w_comp ) + 1.0E-20_wp                       )     &
                  ) *   weight_substep(intermediate_timestep_count)
       ENDDO

       sums_us2_ws_l(nzb,tn) = sums_us2_ws_l(nzb+1,tn)



    END SUBROUTINE advec_u_ws_ij



!-----------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of v-component - Call for grid point i,j
!-----------------------------------------------------------------------------!
   SUBROUTINE advec_v_ws_ij( i, j, i_omp, tn )

       USE arrays_3d,                                                          &
           ONLY:  ddzw, diss_l_v, diss_s_v, flux_l_v, flux_s_v, tend, u, v, w, &
                  drho_air, rho_air, rho_air_zw

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count, u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nysv, nzb, nzb_max, nzt, wall_flags_0

       USE kinds

       USE statistics,                                                         &
           ONLY:  hom, sums_vs2_ws_l, sums_wsvs_ws_l, weight_substep

       IMPLICIT NONE

       INTEGER(iwp)  ::  i      !<
       INTEGER(iwp)  ::  ibit18 !<
       INTEGER(iwp)  ::  ibit19 !<
       INTEGER(iwp)  ::  ibit20 !<
       INTEGER(iwp)  ::  ibit21 !<
       INTEGER(iwp)  ::  ibit22 !<
       INTEGER(iwp)  ::  ibit23 !<
       INTEGER(iwp)  ::  ibit24 !<
       INTEGER(iwp)  ::  ibit25 !<
       INTEGER(iwp)  ::  ibit26 !<
       INTEGER(iwp)  ::  i_omp  !<
       INTEGER(iwp)  ::  j      !<
       INTEGER(iwp)  ::  k      !<
       INTEGER(iwp)  ::  k_mm   !<
       INTEGER(iwp)  ::  k_pp   !<
       INTEGER(iwp)  ::  k_ppp  !<
       INTEGER(iwp)  ::  tn     !<
       
       REAL(wp)     ::  diss_d   !<
       REAL(wp)     ::  div      !<
       REAL(wp)     ::  flux_d   !<
       REAL(wp)     ::  gu       !<
       REAL(wp)     ::  gv       !<
       REAL(wp)     ::  u_comp   !<
       REAL(wp)     ::  v_comp_l !<
       REAL(wp)     ::  w_comp   !<
       
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_n !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_r !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_t !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_n !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_r !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_t !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  v_comp !<

       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans

!       
!--    Compute leftside fluxes for the respective boundary.
       IF ( i == i_omp )  THEN

          DO  k = nzb+1, nzb_max

             ibit20 = IBITS(wall_flags_0(k,j,i-1),20,1)
             ibit19 = IBITS(wall_flags_0(k,j,i-1),19,1)
             ibit18 = IBITS(wall_flags_0(k,j,i-1),18,1)

             u_comp           = u(k,j-1,i) + u(k,j,i) - gu
             flux_l_v(k,j,tn) = u_comp * (                                    &
                              ( 37.0_wp * ibit20 * adv_mom_5                  &
                           +     7.0_wp * ibit19 * adv_mom_3                  &
                           +              ibit18 * adv_mom_1                  &
                              ) *                                             &
                                        ( v(k,j,i)   + v(k,j,i-1) )           &
                       -      (  8.0_wp * ibit20 * adv_mom_5                  &
                           +              ibit19 * adv_mom_3                  &
                              ) *                                             &
                                        ( v(k,j,i+1) + v(k,j,i-2) )           &
                       +      (           ibit20 * adv_mom_5                  &
                              ) *                                             &
                                        ( v(k,j,i+2) + v(k,j,i-3) )           &
                                         )

              diss_l_v(k,j,tn) = - ABS( u_comp ) * (                          &
                              ( 10.0_wp * ibit20 * adv_mom_5                  &
                           +     3.0_wp * ibit19 * adv_mom_3                  &
                           +              ibit18 * adv_mom_1                  &
                              ) *                                             &
                                        ( v(k,j,i)   - v(k,j,i-1) )           &
                       -      (  5.0_wp * ibit20 * adv_mom_5                  &
                           +              ibit19 * adv_mom_3                  &
                              ) *                                             &
                                        ( v(k,j,i+1) - v(k,j,i-2) )           &
                       +      (           ibit20 * adv_mom_5                  &
                              ) *                                             &
                                        ( v(k,j,i+2) - v(k,j,i-3) )           &
                                                   )

          ENDDO

          DO  k = nzb_max+1, nzt

             u_comp           = u(k,j-1,i) + u(k,j,i) - gu
             flux_l_v(k,j,tn) = u_comp * (                                    &
                             37.0_wp * ( v(k,j,i) + v(k,j,i-1)   )            &
                           -  8.0_wp * ( v(k,j,i+1) + v(k,j,i-2) )            &
                           +           ( v(k,j,i+2) + v(k,j,i-3) ) ) * adv_mom_5
             diss_l_v(k,j,tn) = - ABS( u_comp ) * (                           &
                             10.0_wp * ( v(k,j,i) - v(k,j,i-1)   )            &
                           -  5.0_wp * ( v(k,j,i+1) - v(k,j,i-2) )            &
                           +           ( v(k,j,i+2) - v(k,j,i-3) ) ) * adv_mom_5

          ENDDO
          
       ENDIF
!
!--    Compute southside fluxes for the respective boundary.
       IF ( j == nysv )  THEN
       
          DO  k = nzb+1, nzb_max

             ibit23 = IBITS(wall_flags_0(k,j-1,i),23,1)
             ibit22 = IBITS(wall_flags_0(k,j-1,i),22,1)
             ibit21 = IBITS(wall_flags_0(k,j-1,i),21,1)

             v_comp_l       = v(k,j,i) + v(k,j-1,i) - gv
             flux_s_v(k,tn) = v_comp_l * (                                    &
                            ( 37.0_wp * ibit23 * adv_mom_5                    &
                         +     7.0_wp * ibit22 * adv_mom_3                    &
                         +              ibit21 * adv_mom_1                    &
                            ) *                                               &
                                        ( v(k,j,i)   + v(k,j-1,i) )           &
                     -      (  8.0_wp * ibit23 * adv_mom_5                    &
                         +              ibit22 * adv_mom_3                    &
                            ) *                                               &
                                        ( v(k,j+1,i) + v(k,j-2,i) )           &
                     +      (           ibit23 * adv_mom_5                    &
                            ) *                                               &
                                        ( v(k,j+2,i) + v(k,j-3,i) )           &
                                         )

             diss_s_v(k,tn) = - ABS( v_comp_l ) * (                           &
                            ( 10.0_wp * ibit23 * adv_mom_5                    &
                         +     3.0_wp * ibit22 * adv_mom_3                    &
                         +              ibit21 * adv_mom_1                    &
                            ) *                                               &
                                        ( v(k,j,i)   - v(k,j-1,i) )           &
                     -      (  5.0_wp * ibit23 * adv_mom_5                    &
                         +              ibit22 * adv_mom_3                    &
                            ) *                                               &
                                        ( v(k,j+1,i) - v(k,j-2,i) )           &
                     +      (           ibit23 * adv_mom_5                    &
                            ) *                                               &
                                        ( v(k,j+2,i) - v(k,j-3,i) )           &
                                                  )

          ENDDO

          DO  k = nzb_max+1, nzt

             v_comp_l       = v(k,j,i) + v(k,j-1,i) - gv
             flux_s_v(k,tn) = v_comp_l * (                                    &
                           37.0_wp * ( v(k,j,i) + v(k,j-1,i)   )              &
                         -  8.0_wp * ( v(k,j+1,i) + v(k,j-2,i) )              &
                         +           ( v(k,j+2,i) + v(k,j-3,i) ) ) * adv_mom_5
             diss_s_v(k,tn) = - ABS( v_comp_l ) * (                           &
                           10.0_wp * ( v(k,j,i) - v(k,j-1,i)   )              &
                         -  5.0_wp * ( v(k,j+1,i) - v(k,j-2,i) )              &
                         +           ( v(k,j+2,i) - v(k,j-3,i) ) ) * adv_mom_5

          ENDDO
          
       ENDIF

       flux_t(0) = 0.0_wp
       diss_t(0) = 0.0_wp
       flux_d    = 0.0_wp
       diss_d    = 0.0_wp
!
!--    Now compute the fluxes and tendency terms for the horizontal and
!--    verical parts.
       DO  k = nzb+1, nzb_max

          ibit20 = IBITS(wall_flags_0(k,j,i),20,1)
          ibit19 = IBITS(wall_flags_0(k,j,i),19,1)
          ibit18 = IBITS(wall_flags_0(k,j,i),18,1)
 
          u_comp    = u(k,j-1,i+1) + u(k,j,i+1) - gu
          flux_r(k) = u_comp * (                                              &
                     ( 37.0_wp * ibit20 * adv_mom_5                           &
                  +     7.0_wp * ibit19 * adv_mom_3                           &
                  +              ibit18 * adv_mom_1                           &
                     ) *                                                      &
                                    ( v(k,j,i+1) + v(k,j,i)   )               &
              -      (  8.0_wp * ibit20 * adv_mom_5                           &
                  +              ibit19 * adv_mom_3                           &
                     ) *                                                      &
                                    ( v(k,j,i+2) + v(k,j,i-1) )               &
              +      (           ibit20 * adv_mom_5                           &
                     ) *                                                      &
                                    ( v(k,j,i+3) + v(k,j,i-2) )               &
                               )

          diss_r(k) = - ABS( u_comp ) * (                                     &
                     ( 10.0_wp * ibit20 * adv_mom_5                           &
                  +     3.0_wp * ibit19 * adv_mom_3                           &
                  +              ibit18 * adv_mom_1                           &
                     ) *                                                      &
                                    ( v(k,j,i+1) - v(k,j,i)  )                &
              -      (  5.0_wp * ibit20 * adv_mom_5                           &
                  +              ibit19 * adv_mom_3                           &
                     ) *                                                      &
                                    ( v(k,j,i+2) - v(k,j,i-1) )               &
              +      (           ibit20 * adv_mom_5                           &
                     ) *                                                      &
                                    ( v(k,j,i+3) - v(k,j,i-2) )               &
                                        )

          ibit23 = IBITS(wall_flags_0(k,j,i),23,1)
          ibit22 = IBITS(wall_flags_0(k,j,i),22,1)
          ibit21 = IBITS(wall_flags_0(k,j,i),21,1)


          v_comp(k) = v(k,j+1,i) + v(k,j,i)
          flux_n(k) = ( v_comp(k) - gv ) * (                                  &
                     ( 37.0_wp * ibit23 * adv_mom_5                           &
                  +     7.0_wp * ibit22 * adv_mom_3                           &
                  +              ibit21 * adv_mom_1                           &
                     ) *                                                      &
                                    ( v(k,j+1,i) + v(k,j,i)   )               &
              -      (  8.0_wp * ibit23 * adv_mom_5                           &
                  +              ibit22 * adv_mom_3                           &
                     ) *                                                      &
                                    ( v(k,j+2,i) + v(k,j-1,i) )               &
              +      (           ibit23 * adv_mom_5                           &
                     ) *                                                      &
                                    ( v(k,j+3,i) + v(k,j-2,i) )               &
                                           )

          diss_n(k) = - ABS( v_comp(k) - gv ) * (                             &
                     ( 10.0_wp * ibit23 * adv_mom_5                           &
                  +     3.0_wp * ibit22 * adv_mom_3                           &
                  +              ibit21 * adv_mom_1                           &
                     ) *                                                      &
                                    ( v(k,j+1,i) - v(k,j,i)   )               &
              -      (  5.0_wp * ibit23 * adv_mom_5                           &
                  +              ibit22 * adv_mom_3                           &
                     ) *                                                      &
                                    ( v(k,j+2,i) - v(k,j-1,i) )               &
              +      (           ibit23 * adv_mom_5                           &
                     ) *                                                      &
                                    ( v(k,j+3,i) - v(k,j-2,i) )               &
                                                )
!
!--       k index has to be modified near bottom and top, else array
!--       subscripts will be exceeded.
          ibit26 = IBITS(wall_flags_0(k,j,i),26,1)
          ibit25 = IBITS(wall_flags_0(k,j,i),25,1)
          ibit24 = IBITS(wall_flags_0(k,j,i),24,1)

          k_ppp = k + 3 * ibit26
          k_pp  = k + 2 * ( 1 - ibit24  )
          k_mm  = k - 2 * ibit26

          w_comp    = w(k,j-1,i) + w(k,j,i)
          flux_t(k) = w_comp * rho_air_zw(k) * (                              &
                     ( 37.0_wp * ibit26 * adv_mom_5                           &
                  +     7.0_wp * ibit25 * adv_mom_3                           &
                  +              ibit24 * adv_mom_1                           &
                     ) *                                                      &
                                ( v(k+1,j,i)   + v(k,j,i)    )                &
              -      (  8.0_wp * ibit26 * adv_mom_5                           &
                  +              ibit25 * adv_mom_3                           &
                     ) *                                                      &
                                ( v(k_pp,j,i)  + v(k-1,j,i)  )                &
              +      (           ibit26 * adv_mom_5                           &
                     ) *                                                      &
                                ( v(k_ppp,j,i) + v(k_mm,j,i) )                &
                                 )

          diss_t(k) = - ABS( w_comp ) * rho_air_zw(k) * (                     &
                     ( 10.0_wp * ibit26 * adv_mom_5                           &
                  +     3.0_wp * ibit25 * adv_mom_3                           &
                  +              ibit24 * adv_mom_1                           &
                     ) *                                                      &
                                ( v(k+1,j,i)   - v(k,j,i)    )                &
              -      (  5.0_wp * ibit26 * adv_mom_5                           &
                  +              ibit25 * adv_mom_3                           &
                     ) *                                                      &
                                ( v(k_pp,j,i)  - v(k-1,j,i)  )                &
              +      (           ibit26 * adv_mom_5                           &
                     ) *                                                      &
                                ( v(k_ppp,j,i) - v(k_mm,j,i) )                &
                                         )
!
!--       Calculate the divergence of the velocity field. A respective
!--       correction is needed to overcome numerical instabilities introduced
!--       by a not sufficient reduction of divergences near topography. 
          div = ( ( ( u_comp     + gu )                                       &
                                       * ( ibit18 + ibit19 + ibit20 )         &
                  - ( u(k,j-1,i) + u(k,j,i) )                                 &
                                       * ( IBITS(wall_flags_0(k,j,i-1),18,1)  &
                                         + IBITS(wall_flags_0(k,j,i-1),19,1)  &
                                         + IBITS(wall_flags_0(k,j,i-1),20,1)  &
                                         )                                    &
                  ) * rho_air(k) * ddx                                        &
               +  ( v_comp(k)                                                 &
                                       * ( ibit21 + ibit22 + ibit23 )         &
                - ( v(k,j,i)     + v(k,j-1,i) )                               &
                                       * ( IBITS(wall_flags_0(k,j-1,i),21,1)  &
                                         + IBITS(wall_flags_0(k,j-1,i),22,1)  &
                                         + IBITS(wall_flags_0(k,j-1,i),23,1)  &
                                         )                                    &
                  ) * rho_air(k) * ddy                                        &
               +  ( w_comp * rho_air_zw(k) * ( ibit24 + ibit25 + ibit26 )     &
                - ( w(k-1,j-1,i) + w(k-1,j,i) ) * rho_air_zw(k-1)             &
                                       * ( IBITS(wall_flags_0(k-1,j,i),24,1)  &
                                         + IBITS(wall_flags_0(k-1,j,i),25,1)  &
                                         + IBITS(wall_flags_0(k-1,j,i),26,1)  &
                                         )                                    &
                   ) * ddzw(k)   &
                ) * 0.5_wp


          tend(k,j,i) = tend(k,j,i) - (                                       &
                         ( flux_r(k) + diss_r(k)                              &
                       -   flux_l_v(k,j,tn) - diss_l_v(k,j,tn)   ) * ddx      &
                       + ( flux_n(k) + diss_n(k)                              &
                       -   flux_s_v(k,tn) - diss_s_v(k,tn)       ) * ddy      &
                       + ( ( flux_t(k) + diss_t(k) )                          &
                       -   ( flux_d    + diss_d    )                          &
                                                   ) * drho_air(k) * ddzw(k)  &
                                      ) + v(k,j,i) * div

           flux_l_v(k,j,tn) = flux_r(k)
           diss_l_v(k,j,tn) = diss_r(k)
           flux_s_v(k,tn)   = flux_n(k)
           diss_s_v(k,tn)   = diss_n(k)
           flux_d           = flux_t(k)
           diss_d           = diss_t(k)

!
!--        Statistical Evaluation of v'v'. The factor has to be applied for
!--        right evaluation when gallilei_trans = .T. .
           sums_vs2_ws_l(k,tn) = sums_vs2_ws_l(k,tn)                           &
                + ( flux_n(k)                                                  &
                    * ( v_comp(k) - 2.0_wp * hom(k,1,2,0)                   )  &
                    / ( v_comp(k) - gv + SIGN( 1.0E-20_wp, v_comp(k) - gv ) )  &
                  + diss_n(k)                                                  &
                    *   ABS( v_comp(k) - 2.0_wp * hom(k,1,2,0)              )  &
                    / ( ABS( v_comp(k) - gv ) + 1.0E-20_wp                  )  &
                  ) *   weight_substep(intermediate_timestep_count)
!
!--        Statistical Evaluation of w'u'.
           sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn)                         &
                + ( flux_t(k)                                                  &
                    * ( w_comp - 2.0_wp * hom(k,1,3,0)                   )     &
                    / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              )     &
                  + diss_t(k)                                                  &
                    *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              )     &
                    / ( ABS( w_comp ) + 1.0E-20_wp                       )     &
                  ) *   weight_substep(intermediate_timestep_count)

       ENDDO

       DO  k = nzb_max+1, nzt

          u_comp    = u(k,j-1,i+1) + u(k,j,i+1) - gu
          flux_r(k) = u_comp * (                                              &
                      37.0_wp * ( v(k,j,i+1) + v(k,j,i)   )                   &
                    -  8.0_wp * ( v(k,j,i+2) + v(k,j,i-1) )                   &
                    +           ( v(k,j,i+3) + v(k,j,i-2) ) ) * adv_mom_5

          diss_r(k) = - ABS( u_comp ) * (                                     &
                      10.0_wp * ( v(k,j,i+1) - v(k,j,i) )                     &
                    -  5.0_wp * ( v(k,j,i+2) - v(k,j,i-1) )                   &
                    +           ( v(k,j,i+3) - v(k,j,i-2) ) ) * adv_mom_5


          v_comp(k) = v(k,j+1,i) + v(k,j,i)
          flux_n(k) = ( v_comp(k) - gv ) * (                                  &
                      37.0_wp * ( v(k,j+1,i) + v(k,j,i)   )                   &
                    -  8.0_wp * ( v(k,j+2,i) + v(k,j-1,i) )                   &
                      +         ( v(k,j+3,i) + v(k,j-2,i) ) ) * adv_mom_5

          diss_n(k) = - ABS( v_comp(k) - gv ) * (                             &
                      10.0_wp * ( v(k,j+1,i) - v(k,j,i)   )                   &
                    -  5.0_wp * ( v(k,j+2,i) - v(k,j-1,i) )                   &
                    +           ( v(k,j+3,i) - v(k,j-2,i) ) ) * adv_mom_5
!
!--       k index has to be modified near bottom and top, else array
!--       subscripts will be exceeded.
          ibit26 = IBITS(wall_flags_0(k,j,i),26,1)
          ibit25 = IBITS(wall_flags_0(k,j,i),25,1)
          ibit24 = IBITS(wall_flags_0(k,j,i),24,1)

          k_ppp = k + 3 * ibit26
          k_pp  = k + 2 * ( 1 - ibit24  )
          k_mm  = k - 2 * ibit26

          w_comp    = w(k,j-1,i) + w(k,j,i)
          flux_t(k) = w_comp * rho_air_zw(k) * (                              &
                     ( 37.0_wp * ibit26 * adv_mom_5                           &
                  +     7.0_wp * ibit25 * adv_mom_3                           &
                  +              ibit24 * adv_mom_1                           &
                     ) *                                                      &
                                ( v(k+1,j,i)   + v(k,j,i)    )                &
              -      (  8.0_wp * ibit26 * adv_mom_5                           &
                  +              ibit25 * adv_mom_3                           &
                     ) *                                                      &
                                ( v(k_pp,j,i)  + v(k-1,j,i)  )                &
              +      (           ibit26 * adv_mom_5                           &
                     ) *                                                      &
                                ( v(k_ppp,j,i) + v(k_mm,j,i) )                &
                                 )

          diss_t(k) = - ABS( w_comp ) * rho_air_zw(k) * (                     &
                     ( 10.0_wp * ibit26 * adv_mom_5                           &
                  +     3.0_wp * ibit25 * adv_mom_3                           &
                  +              ibit24 * adv_mom_1                           &
                     ) *                                                      &
                                ( v(k+1,j,i)   - v(k,j,i)    )                &
              -      (  5.0_wp * ibit26 * adv_mom_5                           &
                  +              ibit25 * adv_mom_3                           &
                     ) *                                                      &
                                ( v(k_pp,j,i)  - v(k-1,j,i)  )                &
              +      (           ibit26 * adv_mom_5                           &
                     ) *                                                      &
                                ( v(k_ppp,j,i) - v(k_mm,j,i) )                &
                                         )
!
!--       Calculate the divergence of the velocity field. A respective
!--       correction is needed to overcome numerical instabilities introduced
!--       by a not sufficient reduction of divergences near topography.
          div = ( ( u_comp + gu - ( u(k,j-1,i)   + u(k,j,i)   ) ) * ddx       &
                                                                  * rho_air(k)&
               +  ( v_comp(k)   - ( v(k,j,i)     + v(k,j-1,i) ) ) * ddy       &
                                                                  * rho_air(k)&
               +  (   w_comp                      * rho_air_zw(k)   -         &
                    ( w(k-1,j-1,i) + w(k-1,j,i) ) * rho_air_zw(k-1)           &
                  ) * ddzw(k)                                                 &
                ) * 0.5_wp

          tend(k,j,i) = tend(k,j,i) - (                                       &
                         ( flux_r(k) + diss_r(k)                              &
                       -   flux_l_v(k,j,tn) - diss_l_v(k,j,tn)   ) * ddx      &
                       + ( flux_n(k) + diss_n(k)                              &
                       -   flux_s_v(k,tn) - diss_s_v(k,tn)       ) * ddy      &
                       + ( ( flux_t(k) + diss_t(k) )                          &
                       -   ( flux_d    + diss_d    )                          &
                                                   ) * drho_air(k) * ddzw(k)  &
                                      ) + v(k,j,i) * div

           flux_l_v(k,j,tn) = flux_r(k)
           diss_l_v(k,j,tn) = diss_r(k)
           flux_s_v(k,tn)   = flux_n(k)
           diss_s_v(k,tn)   = diss_n(k)
           flux_d           = flux_t(k)
           diss_d           = diss_t(k)

!
!--        Statistical Evaluation of v'v'. The factor has to be applied for
!--        right evaluation when gallilei_trans = .T. .
           sums_vs2_ws_l(k,tn) = sums_vs2_ws_l(k,tn)                           &
                + ( flux_n(k)                                                  &
                    * ( v_comp(k) - 2.0_wp * hom(k,1,2,0)                   )  &
                    / ( v_comp(k) - gv + SIGN( 1.0E-20_wp, v_comp(k) - gv ) )  &
                  + diss_n(k)                                                  &
                    *   ABS( v_comp(k) - 2.0_wp * hom(k,1,2,0)              )  &
                    / ( ABS( v_comp(k) - gv ) + 1.0E-20_wp                  )  &
                  ) *   weight_substep(intermediate_timestep_count)
!
!--        Statistical Evaluation of w'u'.
           sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn)                         &
                + ( flux_t(k)                                                  &
                    * ( w_comp - 2.0_wp * hom(k,1,3,0)                   )     &
                    / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              )     &
                  + diss_t(k)                                                  &
                    *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              )     &
                    / ( ABS( w_comp ) + 1.0E-20_wp                       )     &
                  ) *   weight_substep(intermediate_timestep_count)

       ENDDO
       sums_vs2_ws_l(nzb,tn) = sums_vs2_ws_l(nzb+1,tn)


    END SUBROUTINE advec_v_ws_ij



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of w-component - Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_w_ws_ij( i, j, i_omp, tn )

       USE arrays_3d,                                                         &
           ONLY:  ddzu, diss_l_w, diss_s_w, flux_l_w, flux_s_w, tend, u, v, w,&
                  drho_air_zw, rho_air, rho_air_zw

       USE constants,                                                         &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                &
           ONLY:  intermediate_timestep_count, u_gtrans, v_gtrans

       USE grid_variables,                                                    &
           ONLY:  ddx, ddy

       USE indices,                                                           &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzb_max, nzt, wall_flags_0,        &
                  wall_flags_00

       USE kinds
       
       USE statistics,                                                        &
           ONLY:  hom, sums_ws2_ws_l, weight_substep

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit27 !<
       INTEGER(iwp) ::  ibit28 !<
       INTEGER(iwp) ::  ibit29 !<
       INTEGER(iwp) ::  ibit30 !<
       INTEGER(iwp) ::  ibit31 !<
       INTEGER(iwp) ::  ibit32 !<
       INTEGER(iwp) ::  ibit33 !<
       INTEGER(iwp) ::  ibit34 !<
       INTEGER(iwp) ::  ibit35 !<
       INTEGER(iwp) ::  i_omp  !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn     !<
       
       REAL(wp)    ::  diss_d  !<
       REAL(wp)    ::  div     !<
       REAL(wp)    ::  flux_d  !<
       REAL(wp)    ::  gu      !<
       REAL(wp)    ::  gv      !<
       REAL(wp)    ::  u_comp  !<
       REAL(wp)    ::  v_comp  !<
       REAL(wp)    ::  w_comp  !<
       
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_n !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_r !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_t !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_n !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_r !<
       REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_t !<

       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans

!
!--    Compute southside fluxes for the respective boundary.
       IF ( j == nys )  THEN

          DO  k = nzb+1, nzb_max
             ibit32 = IBITS(wall_flags_00(k,j-1,i),0,1)
             ibit31 = IBITS(wall_flags_0(k,j-1,i),31,1)
             ibit30 = IBITS(wall_flags_0(k,j-1,i),30,1)

             v_comp         = v(k+1,j,i) + v(k,j,i) - gv
             flux_s_w(k,tn) = v_comp * (                                      &
                            ( 37.0_wp * ibit32 * adv_mom_5                    &
                         +     7.0_wp * ibit31 * adv_mom_3                    &
                         +              ibit30 * adv_mom_1                    &
                            ) *                                               &
                                        ( w(k,j,i)   + w(k,j-1,i) )           &
                     -      (  8.0_wp * ibit32 * adv_mom_5                    &
                         +              ibit31 * adv_mom_3                    &
                            ) *                                               &
                                        ( w(k,j+1,i) + w(k,j-2,i) )           &
                     +      (           ibit32 * adv_mom_5                    &
                            ) *                                               &
                                        ( w(k,j+2,i) + w(k,j-3,i) )           &
                                       )

             diss_s_w(k,tn) = - ABS( v_comp ) * (                             &
                            ( 10.0_wp * ibit32 * adv_mom_5                    &
                         +     3.0_wp * ibit31 * adv_mom_3                    &
                         +              ibit30 * adv_mom_1                    &
                            ) *                                               &
                                        ( w(k,j,i)   - w(k,j-1,i) )           &
                     -      (  5.0_wp * ibit32 * adv_mom_5                    &
                         +              ibit31 * adv_mom_3                    &
                            ) *                                               &
                                        ( w(k,j+1,i) - w(k,j-2,i) )           &
                     +      (           ibit32 * adv_mom_5                    &
                            ) *                                               &
                                        ( w(k,j+2,i) - w(k,j-3,i) )           &
                                                )

          ENDDO

          DO  k = nzb_max+1, nzt

             v_comp         = v(k+1,j,i) + v(k,j,i) - gv
             flux_s_w(k,tn) = v_comp * (                                      &
                           37.0_wp * ( w(k,j,i) + w(k,j-1,i)   )              &
                         -  8.0_wp * ( w(k,j+1,i) +w(k,j-2,i)  )              &
                         +           ( w(k,j+2,i) + w(k,j-3,i) ) ) * adv_mom_5
             diss_s_w(k,tn) = - ABS( v_comp ) * (                             &
                           10.0_wp * ( w(k,j,i) - w(k,j-1,i)   )              &
                         -  5.0_wp * ( w(k,j+1,i) - w(k,j-2,i) )              &
                         +           ( w(k,j+2,i) - w(k,j-3,i) ) ) * adv_mom_5

          ENDDO

       ENDIF
!
!--    Compute leftside fluxes for the respective boundary.
       IF ( i == i_omp ) THEN

          DO  k = nzb+1, nzb_max

             ibit29 = IBITS(wall_flags_0(k,j,i-1),29,1)
             ibit28 = IBITS(wall_flags_0(k,j,i-1),28,1)
             ibit27 = IBITS(wall_flags_0(k,j,i-1),27,1)

             u_comp           = u(k+1,j,i) + u(k,j,i) - gu
             flux_l_w(k,j,tn) = u_comp * (                                    &
                             ( 37.0_wp * ibit29 * adv_mom_5                   &
                          +     7.0_wp * ibit28 * adv_mom_3                   &
                          +              ibit27 * adv_mom_1                   &
                             ) *                                              &
                                        ( w(k,j,i)   + w(k,j,i-1) )           &
                      -      (  8.0_wp * ibit29 * adv_mom_5                   &
                          +              ibit28 * adv_mom_3                   &
                             ) *                                              &
                                        ( w(k,j,i+1) + w(k,j,i-2) )           &
                      +      (           ibit29 * adv_mom_5                   &
                             ) *                                              &
                                        ( w(k,j,i+2) + w(k,j,i-3) )           &
                                         )

               diss_l_w(k,j,tn) = - ABS( u_comp ) * (                         &
                             ( 10.0_wp * ibit29 * adv_mom_5                   &
                          +     3.0_wp * ibit28 * adv_mom_3                   &
                          +              ibit27 * adv_mom_1                   &
                             ) *                                              &
                                        ( w(k,j,i)   - w(k,j,i-1) )           &
                      -      (  5.0_wp * ibit29 * adv_mom_5                   &
                          +              ibit28 * adv_mom_3                   &
                             ) *                                              &
                                        ( w(k,j,i+1) - w(k,j,i-2) )           &
                      +      (           ibit29 * adv_mom_5                   &
                             ) *                                              &
                                        ( w(k,j,i+2) - w(k,j,i-3) )           &
                                                    )

          ENDDO

          DO  k = nzb_max+1, nzt

             u_comp           = u(k+1,j,i) + u(k,j,i) - gu
             flux_l_w(k,j,tn) = u_comp * (                                    &
                            37.0_wp * ( w(k,j,i) + w(k,j,i-1)   )             &
                          -  8.0_wp * ( w(k,j,i+1) + w(k,j,i-2) )             &
                          +           ( w(k,j,i+2) + w(k,j,i-3) ) ) * adv_mom_5
             diss_l_w(k,j,tn) = - ABS( u_comp ) * (                           &
                            10.0_wp * ( w(k,j,i) - w(k,j,i-1)   )             &
                          -  5.0_wp * ( w(k,j,i+1) - w(k,j,i-2) )             &
                          +           ( w(k,j,i+2) - w(k,j,i-3) ) ) * adv_mom_5 

          ENDDO

       ENDIF
!
!--    The lower flux has to be calculated explicetely for the tendency at
!--    the first w-level. For topography wall this is done implicitely by
!--    wall_flags_0.
       k         = nzb + 1
       w_comp    = w(k,j,i) + w(k-1,j,i)
       flux_t(0) = w_comp       * ( w(k,j,i) + w(k-1,j,i) ) * adv_mom_1
       diss_t(0) = -ABS(w_comp) * ( w(k,j,i) - w(k-1,j,i) ) * adv_mom_1
       flux_d    = flux_t(0)
       diss_d    = diss_t(0)
!
!--    Now compute the fluxes and tendency terms for the horizontal
!--    and vertical parts.
       DO  k = nzb+1, nzb_max

          ibit29 = IBITS(wall_flags_0(k,j,i),29,1)
          ibit28 = IBITS(wall_flags_0(k,j,i),28,1)
          ibit27 = IBITS(wall_flags_0(k,j,i),27,1)

          u_comp    = u(k+1,j,i+1) + u(k,j,i+1) - gu
          flux_r(k) = u_comp * (                                              &
                     ( 37.0_wp * ibit29 * adv_mom_5                           &
                  +     7.0_wp * ibit28 * adv_mom_3                           &
                  +              ibit27 * adv_mom_1                           &
                     ) *                                                      &
                                    ( w(k,j,i+1) + w(k,j,i)   )               &
              -      (  8.0_wp * ibit29 * adv_mom_5                           &
                  +              ibit28 * adv_mom_3                           &
                     ) *                                                      &
                                    ( w(k,j,i+2) + w(k,j,i-1) )               &
              +      (           ibit29 * adv_mom_5                           &
                     ) *                                                      &
                                    ( w(k,j,i+3) + w(k,j,i-2) )               &
                               )

          diss_r(k) = - ABS( u_comp ) * (                                     &
                     ( 10.0_wp * ibit29 * adv_mom_5                           &
                  +     3.0_wp * ibit28 * adv_mom_3                           &
                  +              ibit27 * adv_mom_1                           &
                     ) *                                                      &
                                    ( w(k,j,i+1) - w(k,j,i)   )               &
              -      (  5.0_wp * ibit29 * adv_mom_5                           &
                  +              ibit28 * adv_mom_3                           &
                     ) *                                                      &
                                    ( w(k,j,i+2) - w(k,j,i-1) )               &
              +      (           ibit29 * adv_mom_5                           &
                     ) *                                                      &
                                    ( w(k,j,i+3) - w(k,j,i-2) )               &
                                        )

          ibit32 = IBITS(wall_flags_00(k,j,i),0,1)
          ibit31 = IBITS(wall_flags_0(k,j,i),31,1)
          ibit30 = IBITS(wall_flags_0(k,j,i),30,1)

          v_comp    = v(k+1,j+1,i) + v(k,j+1,i) - gv
          flux_n(k) = v_comp * (                                              &
                     ( 37.0_wp * ibit32 * adv_mom_5                           &
                  +     7.0_wp * ibit31 * adv_mom_3                           &
                  +              ibit30 * adv_mom_1                           &
                     ) *                                                      &
                                    ( w(k,j+1,i) + w(k,j,i)   )               &
              -      (  8.0_wp * ibit32 * adv_mom_5                           &
                  +              ibit31 * adv_mom_3                           &
                     ) *                                                      &
                                    ( w(k,j+2,i) + w(k,j-1,i) )               &
              +      (           ibit32 * adv_mom_5                           &
                     ) *                                                      &
                                    ( w(k,j+3,i) + w(k,j-2,i) )               &
                               )

          diss_n(k) = - ABS( v_comp ) * (                                     &
                     ( 10.0_wp * ibit32 * adv_mom_5                           &
                  +     3.0_wp * ibit31 * adv_mom_3                           &
                  +              ibit30 * adv_mom_1                           &
                     ) *                                                      &
                                    ( w(k,j+1,i) - w(k,j,i)  )                &
              -      (  5.0_wp * ibit32 * adv_mom_5                           &
                  +              ibit31 * adv_mom_3                           &
                     ) *                                                      &
                                   ( w(k,j+2,i) - w(k,j-1,i) )                &
              +      (           ibit32 * adv_mom_5                           &
                     ) *                                                      &
                                   ( w(k,j+3,i) - w(k,j-2,i) )                &
                                        )
!
!--       k index has to be modified near bottom and top, else array
!--       subscripts will be exceeded.
          ibit35 = IBITS(wall_flags_00(k,j,i),3,1)
          ibit34 = IBITS(wall_flags_00(k,j,i),2,1)
          ibit33 = IBITS(wall_flags_00(k,j,i),1,1)

          k_ppp = k + 3 * ibit35
          k_pp  = k + 2 * ( 1 - ibit33  )
          k_mm  = k - 2 * ibit35

          w_comp    = w(k+1,j,i) + w(k,j,i)
          flux_t(k) = w_comp * rho_air(k+1) * (                               &
                     ( 37.0_wp * ibit35 * adv_mom_5                           &
                  +     7.0_wp * ibit34 * adv_mom_3                           &
                  +              ibit33 * adv_mom_1                           &
                     ) *                                                      &
                                ( w(k+1,j,i)  + w(k,j,i)     )                &
              -      (  8.0_wp * ibit35 * adv_mom_5                           &
                  +              ibit34 * adv_mom_3                           &
                     ) *                                                      &
                                ( w(k_pp,j,i)  + w(k-1,j,i)  )                &
              +      (           ibit35 * adv_mom_5                           &
                     ) *                                                      &
                                ( w(k_ppp,j,i) + w(k_mm,j,i) )                &
                                )

          diss_t(k) = - ABS( w_comp ) * rho_air(k+1) * (                      &
                     ( 10.0_wp * ibit35 * adv_mom_5                           &
                  +     3.0_wp * ibit34 * adv_mom_3                           &
                  +              ibit33 * adv_mom_1                           &
                     ) *                                                      &
                                ( w(k+1,j,i)   - w(k,j,i)    )                &
              -      (  5.0_wp * ibit35 * adv_mom_5                           &
                  +              ibit34 * adv_mom_3                           &
                     ) *                                                      &
                                ( w(k_pp,j,i)  - w(k-1,j,i)  )                &
              +      (           ibit35 * adv_mom_5                           &
                     ) *                                                      &
                                ( w(k_ppp,j,i) - w(k_mm,j,i) )                &
                                        )

!
!--       Calculate the divergence of the velocity field. A respective
!--       correction is needed to overcome numerical instabilities introduced
!--       by a not sufficient reduction of divergences near topography. 
          div = ( ( ( u_comp + gu ) * ( ibit27 + ibit28 + ibit29 )            &
                  - ( u(k+1,j,i) + u(k,j,i)   )                               & 
                                    * ( IBITS(wall_flags_0(k,j,i-1),27,1)     &
                                      + IBITS(wall_flags_0(k,j,i-1),28,1)     &
                                      + IBITS(wall_flags_0(k,j,i-1),29,1)     &
                                      )                                       &
                  ) * rho_air_zw(k) * ddx                                     &
              +   ( ( v_comp + gv ) * ( ibit30 + ibit31 + ibit32 )            & 
                  - ( v(k+1,j,i) + v(k,j,i)   )                               &
                                    * ( IBITS(wall_flags_0(k,j-1,i),30,1)     &
                                      + IBITS(wall_flags_0(k,j-1,i),31,1)     &
                                      + IBITS(wall_flags_00(k,j-1,i),0,1)     &
                                      )                                       &
                  ) * rho_air_zw(k) * ddy                                     &
              +   ( w_comp * rho_air(k+1) * ( ibit33 + ibit34 + ibit35 )      &
                - ( w(k,j,i)   + w(k-1,j,i)   ) * rho_air(k)                  &
                                    * ( IBITS(wall_flags_00(k-1,j,i),1,1)     &
                                      + IBITS(wall_flags_00(k-1,j,i),2,1)     &
                                      + IBITS(wall_flags_00(k-1,j,i),3,1)     &
                                      )                                       & 
                  ) * ddzu(k+1)   &
                ) * 0.5_wp


          tend(k,j,i) = tend(k,j,i) - (                                       &
                      ( flux_r(k) + diss_r(k)                                 &
                    -   flux_l_w(k,j,tn) - diss_l_w(k,j,tn)   ) * ddx         &
                    + ( flux_n(k) + diss_n(k)                                 &
                    -   flux_s_w(k,tn) - diss_s_w(k,tn)       ) * ddy         &
                    + ( ( flux_t(k) + diss_t(k) )                             &
                    -   ( flux_d    + diss_d    )                             &
                                              ) * drho_air_zw(k) * ddzu(k+1)  &
                                      ) + div * w(k,j,i)

          flux_l_w(k,j,tn) = flux_r(k)
          diss_l_w(k,j,tn) = diss_r(k)
          flux_s_w(k,tn)   = flux_n(k)
          diss_s_w(k,tn)   = diss_n(k)
          flux_d           = flux_t(k)
          diss_d           = diss_t(k)
!
!--       Statistical Evaluation of w'w'.
          sums_ws2_ws_l(k,tn)  = sums_ws2_ws_l(k,tn)                          &
                      + ( flux_t(k)                                           &
                       * ( w_comp - 2.0_wp * hom(k,1,3,0)                   ) &
                       / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              ) &
                        + diss_t(k)                                           &
                       *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              ) &
                       / ( ABS( w_comp ) + 1.0E-20_wp                       ) &
                        ) *   weight_substep(intermediate_timestep_count)

       ENDDO

       DO  k = nzb_max+1, nzt

          u_comp    = u(k+1,j,i+1) + u(k,j,i+1) - gu
          flux_r(k) = u_comp * (                                              &
                      37.0_wp * ( w(k,j,i+1) + w(k,j,i)   )                   &
                    -  8.0_wp * ( w(k,j,i+2) + w(k,j,i-1) )                   &
                    +           ( w(k,j,i+3) + w(k,j,i-2) ) ) * adv_mom_5

          diss_r(k) = - ABS( u_comp ) * (                                     &
                      10.0_wp * ( w(k,j,i+1) - w(k,j,i)   )                   &
                    -  5.0_wp * ( w(k,j,i+2) - w(k,j,i-1) )                   &
                    +           ( w(k,j,i+3) - w(k,j,i-2) ) ) * adv_mom_5

          v_comp    = v(k+1,j+1,i) + v(k,j+1,i) - gv
          flux_n(k) = v_comp * (                                              &
                      37.0_wp * ( w(k,j+1,i) + w(k,j,i)   )                   &
                    -  8.0_wp * ( w(k,j+2,i) + w(k,j-1,i) )                   &
                    +           ( w(k,j+3,i) + w(k,j-2,i) ) ) * adv_mom_5

          diss_n(k) = - ABS( v_comp ) * (                                     &
                      10.0_wp * ( w(k,j+1,i) - w(k,j,i)   )                   &
                    -  5.0_wp * ( w(k,j+2,i) - w(k,j-1,i) )                   &
                    +           ( w(k,j+3,i) - w(k,j-2,i) ) ) * adv_mom_5
!
!--       k index has to be modified near bottom and top, else array
!--       subscripts will be exceeded.
          ibit35 = IBITS(wall_flags_00(k,j,i),3,1)
          ibit34 = IBITS(wall_flags_00(k,j,i),2,1)
          ibit33 = IBITS(wall_flags_00(k,j,i),1,1)

          k_ppp = k + 3 * ibit35
          k_pp  = k + 2 * ( 1 - ibit33  )
          k_mm  = k - 2 * ibit35

          w_comp    = w(k+1,j,i) + w(k,j,i)
          flux_t(k) = w_comp * rho_air(k+1) * (                               &
                     ( 37.0_wp * ibit35 * adv_mom_5                           &
                  +     7.0_wp * ibit34 * adv_mom_3                           &
                  +              ibit33 * adv_mom_1                           &
                     ) *                                                      &
                                ( w(k+1,j,i)  + w(k,j,i)     )                &
              -      (  8.0_wp * ibit35 * adv_mom_5                           &
                  +              ibit34 * adv_mom_3                           &
                     ) *                                                      &
                                ( w(k_pp,j,i)  + w(k-1,j,i)  )                &
              +      (           ibit35 * adv_mom_5                           &
                     ) *                                                      &
                                ( w(k_ppp,j,i) + w(k_mm,j,i) )                &
                                )

          diss_t(k) = - ABS( w_comp ) * rho_air(k+1) * (                      &
                     ( 10.0_wp * ibit35 * adv_mom_5                           &
                  +     3.0_wp * ibit34 * adv_mom_3                           &
                  +              ibit33 * adv_mom_1                           &
                     ) *                                                      &
                                ( w(k+1,j,i)   - w(k,j,i)    )                &
              -      (  5.0_wp * ibit35 * adv_mom_5                           &
                  +              ibit34 * adv_mom_3                           &
                     ) *                                                      &
                                ( w(k_pp,j,i)  - w(k-1,j,i)  )                &
              +      (           ibit35 * adv_mom_5                           &
                     ) *                                                      &
                                ( w(k_ppp,j,i) - w(k_mm,j,i) )                &
                                        )
!
!--       Calculate the divergence of the velocity field. A respective
!--       correction is needed to overcome numerical instabilities introduced
!--       by a not sufficient reduction of divergences near topography.
          div = ( ( u_comp + gu - ( u(k+1,j,i) + u(k,j,i)   ) ) * ddx         &
                                                              * rho_air_zw(k) &
              +   ( v_comp + gv - ( v(k+1,j,i) + v(k,j,i)   ) ) * ddy         &
                                                              * rho_air_zw(k) &
              +   (   w_comp                    * rho_air(k+1) -              &
                    ( w(k,j,i)   + w(k-1,j,i) ) * rho_air(k)                  &
                  ) * ddzu(k+1)                                               &
                ) * 0.5_wp

          tend(k,j,i) = tend(k,j,i) - (                                       &
                      ( flux_r(k) + diss_r(k)                                 &
                    -   flux_l_w(k,j,tn) - diss_l_w(k,j,tn)   ) * ddx         &
                    + ( flux_n(k) + diss_n(k)                                 &
                    -   flux_s_w(k,tn) - diss_s_w(k,tn)       ) * ddy         &
                    + ( ( flux_t(k) + diss_t(k) )                             &
                    -   ( flux_d    + diss_d    )                             &
                                              ) * drho_air_zw(k) * ddzu(k+1)  &
                                      ) + div * w(k,j,i)

          flux_l_w(k,j,tn) = flux_r(k)
          diss_l_w(k,j,tn) = diss_r(k)
          flux_s_w(k,tn)   = flux_n(k)
          diss_s_w(k,tn)   = diss_n(k)
          flux_d           = flux_t(k)
          diss_d           = diss_t(k)
!
!--       Statistical Evaluation of w'w'.
          sums_ws2_ws_l(k,tn)  = sums_ws2_ws_l(k,tn)                          &
                      + ( flux_t(k)                                           &
                       * ( w_comp - 2.0_wp * hom(k,1,3,0)                   ) &
                       / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              ) &
                        + diss_t(k)                                           &
                       *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              ) &
                       / ( ABS( w_comp ) + 1.0E-20_wp                       ) &
                        ) *   weight_substep(intermediate_timestep_count)

       ENDDO


    END SUBROUTINE advec_w_ws_ij
    

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Scalar advection - Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_ws( sk, sk_char )

       USE arrays_3d,                                                         &
           ONLY:  ddzw, drho_air, tend, u, v, w, rho_air, rho_air_zw

       USE constants,                                                         &
           ONLY:  adv_sca_1, adv_sca_3, adv_sca_5

       USE control_parameters,                                                &
           ONLY:  intermediate_timestep_count, monotonic_adjustment, u_gtrans,&
                  v_gtrans 

       USE grid_variables,                                                    &
           ONLY:  ddx, ddy

       USE indices,                                                           &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzb_max,   &
                  nzt, wall_flags_0
           
       USE kinds
       
       USE statistics,                                                        &
           ONLY:  hom, sums_wspts_ws_l, sums_wsqs_ws_l, sums_wssas_ws_l,      &
                  sums_wsqrs_ws_l, sums_wsnrs_ws_l, sums_wsss_ws_l,           &
                  weight_substep

       IMPLICIT NONE

       CHARACTER (LEN = *), INTENT(IN)    ::  sk_char !<
       
       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit0  !<
       INTEGER(iwp) ::  ibit1  !<
       INTEGER(iwp) ::  ibit2  !<
       INTEGER(iwp) ::  ibit3  !<
       INTEGER(iwp) ::  ibit4  !<
       INTEGER(iwp) ::  ibit5  !<
       INTEGER(iwp) ::  ibit6  !<
       INTEGER(iwp) ::  ibit7  !<
       INTEGER(iwp) ::  ibit8  !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_mmm  !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<
       
#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  sk !<
#endif

       REAL(wp) ::  diss_d !<
       REAL(wp) ::  div    !<
       REAL(wp) ::  flux_d !<
       REAL(wp) ::  fd_1   !<
       REAL(wp) ::  fl_1   !<
       REAL(wp) ::  fn_1   !<
       REAL(wp) ::  fr_1   !<
       REAL(wp) ::  fs_1   !<
       REAL(wp) ::  ft_1   !<
       REAL(wp) ::  phi_d  !<
       REAL(wp) ::  phi_l  !<
       REAL(wp) ::  phi_n  !<
       REAL(wp) ::  phi_r  !<
       REAL(wp) ::  phi_s  !<
       REAL(wp) ::  phi_t  !<
       REAL(wp) ::  rd     !<
       REAL(wp) ::  rl     !<
       REAL(wp) ::  rn     !<
       REAL(wp) ::  rr     !<
       REAL(wp) ::  rs     !<
       REAL(wp) ::  rt     !<
       REAL(wp) ::  u_comp !<
       REAL(wp) ::  v_comp !<
       
       REAL(wp), DIMENSION(nzb:nzt)   ::  diss_n !<
       REAL(wp), DIMENSION(nzb:nzt)   ::  diss_r !<
       REAL(wp), DIMENSION(nzb:nzt)   ::  diss_t !<
       REAL(wp), DIMENSION(nzb:nzt)   ::  flux_n !<
       REAL(wp), DIMENSION(nzb:nzt)   ::  flux_r !<
       REAL(wp), DIMENSION(nzb:nzt)   ::  flux_t !<
       
       REAL(wp), DIMENSION(nzb+1:nzt) ::  swap_diss_y_local !<
       REAL(wp), DIMENSION(nzb+1:nzt) ::  swap_flux_y_local !<
       
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_diss_x_local !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_flux_x_local !<
        

!
!--    Compute the fluxes for the whole left boundary of the processor domain.
       i = nxl
       DO  j = nys, nyn

          DO  k = nzb+1, nzb_max

             ibit2 = IBITS(wall_flags_0(k,j,i-1),2,1)
             ibit1 = IBITS(wall_flags_0(k,j,i-1),1,1)
             ibit0 = IBITS(wall_flags_0(k,j,i-1),0,1)

             u_comp                 = u(k,j,i) - u_gtrans
             swap_flux_x_local(k,j) = u_comp * (                              &
                                             ( 37.0_wp * ibit2 * adv_sca_5    &
                                          +     7.0_wp * ibit1 * adv_sca_3    &
                                          +              ibit0 * adv_sca_1    &
                                             ) *                              &
                                          ( sk(k,j,i)   + sk(k,j,i-1)    )    &
                                      -      (  8.0_wp * ibit2 * adv_sca_5    &
                                          +              ibit1 * adv_sca_3    &
                                             ) *                              &
                                          ( sk(k,j,i+1) + sk(k,j,i-2)    )    &
                                      +      (           ibit2 * adv_sca_5    & 
                                             ) *                              &
                                          ( sk(k,j,i+2) + sk(k,j,i-3)    )    &
                                               )

              swap_diss_x_local(k,j) = -ABS( u_comp ) * (                     &
                                             ( 10.0_wp * ibit2 * adv_sca_5    &
                                          +     3.0_wp * ibit1 * adv_sca_3    &
                                          +              ibit0 * adv_sca_1    &
                                             ) *                              &
                                          ( sk(k,j,i)   - sk(k,j,i-1) )       &
                                      -      (  5.0_wp * ibit2 * adv_sca_5    &
                                          +              ibit1 * adv_sca_3    &
                                             ) *                              &
                                         ( sk(k,j,i+1) - sk(k,j,i-2)  )       &
                                      +      (           ibit2 * adv_sca_5    &
                                             ) *                              &
                                          ( sk(k,j,i+2) - sk(k,j,i-3) )       &
                                                        )

          ENDDO

          DO  k = nzb_max+1, nzt

             u_comp                 = u(k,j,i) - u_gtrans
             swap_flux_x_local(k,j) = u_comp * (                              &
                                      37.0_wp * ( sk(k,j,i)   + sk(k,j,i-1) ) &
                                    -  8.0_wp * ( sk(k,j,i+1) + sk(k,j,i-2) ) &
                                    +           ( sk(k,j,i+2) + sk(k,j,i-3) ) &
                                               ) * adv_sca_5

             swap_diss_x_local(k,j) = -ABS( u_comp ) * (                      &
                                      10.0_wp * ( sk(k,j,i)   - sk(k,j,i-1) ) &
                                    -  5.0_wp * ( sk(k,j,i+1) - sk(k,j,i-2) ) &
                                    +           ( sk(k,j,i+2) - sk(k,j,i-3) ) &
                                                       ) * adv_sca_5

          ENDDO

       ENDDO

       DO  i = nxl, nxr

          j = nys
          DO  k = nzb+1, nzb_max

             ibit5 = IBITS(wall_flags_0(k,j-1,i),5,1)
             ibit4 = IBITS(wall_flags_0(k,j-1,i),4,1)
             ibit3 = IBITS(wall_flags_0(k,j-1,i),3,1)

             v_comp               = v(k,j,i) - v_gtrans
             swap_flux_y_local(k) = v_comp * (                                &
                                             ( 37.0_wp * ibit5 * adv_sca_5    &
                                          +     7.0_wp * ibit4 * adv_sca_3    &
                                          +              ibit3 * adv_sca_1    &
                                             ) *                              &
                                         ( sk(k,j,i)  + sk(k,j-1,i)     )     &
                                       -     (  8.0_wp * ibit5 * adv_sca_5    &
                                          +              ibit4 * adv_sca_3    &
                                              ) *                             &
                                         ( sk(k,j+1,i) + sk(k,j-2,i)    )     &
                                      +      (           ibit5 * adv_sca_5    &
                                             ) *                              &
                                        ( sk(k,j+2,i) + sk(k,j-3,i)     )     &
                                             )

             swap_diss_y_local(k) = -ABS( v_comp ) * (                        &
                                             ( 10.0_wp * ibit5 * adv_sca_5    &
                                          +     3.0_wp * ibit4 * adv_sca_3    &
                                          +              ibit3 * adv_sca_1    &
                                             ) *                              &
                                          ( sk(k,j,i)   - sk(k,j-1,i)    )    &
                                      -      (  5.0_wp * ibit5 * adv_sca_5    &
                                          +              ibit4 * adv_sca_3    &
                                             ) *                              &
                                          ( sk(k,j+1,i) - sk(k,j-2,i)    )    &
                                      +      (           ibit5 * adv_sca_5    &
                                             ) *                              &
                                          ( sk(k,j+2,i) - sk(k,j-3,i)    )    &
                                                     )

          ENDDO
!
!--       Above to the top of the highest topography. No degradation necessary.
          DO  k = nzb_max+1, nzt

             v_comp               = v(k,j,i) - v_gtrans
             swap_flux_y_local(k) = v_comp * (                               &
                                    37.0_wp * ( sk(k,j,i)   + sk(k,j-1,i) )  &
                                  -  8.0_wp * ( sk(k,j+1,i) + sk(k,j-2,i) )  &
                                  +           ( sk(k,j+2,i) + sk(k,j-3,i) )  &
                                             ) * adv_sca_5
              swap_diss_y_local(k) = -ABS( v_comp ) * (                      &
                                    10.0_wp * ( sk(k,j,i)   - sk(k,j-1,i) )  &
                                  -  5.0_wp * ( sk(k,j+1,i) - sk(k,j-2,i) )  &
                                  +             sk(k,j+2,i) - sk(k,j-3,i)    &
                                                      ) * adv_sca_5

          ENDDO

          DO  j = nys, nyn

             flux_t(0) = 0.0_wp
             diss_t(0) = 0.0_wp
             flux_d    = 0.0_wp
             diss_d    = 0.0_wp

             DO  k = nzb+1, nzb_max

                ibit2 = IBITS(wall_flags_0(k,j,i),2,1)
                ibit1 = IBITS(wall_flags_0(k,j,i),1,1)
                ibit0 = IBITS(wall_flags_0(k,j,i),0,1)

                u_comp    = u(k,j,i+1) - u_gtrans
                flux_r(k) = u_comp * (                                        &
                          ( 37.0_wp * ibit2 * adv_sca_5                       &
                      +      7.0_wp * ibit1 * adv_sca_3                       &
                      +               ibit0 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j,i+1) + sk(k,j,i)   )                    &
                   -      (  8.0_wp * ibit2 * adv_sca_5                       &
                       +              ibit1 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j,i+2) + sk(k,j,i-1) )                    &
                   +      (           ibit2 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j,i+3) + sk(k,j,i-2) )                    &
                                     )

                diss_r(k) = -ABS( u_comp ) * (                                &
                          ( 10.0_wp * ibit2 * adv_sca_5                       &
                       +     3.0_wp * ibit1 * adv_sca_3                       &
                       +              ibit0 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j,i+1) - sk(k,j,i)   )                    &
                   -      (  5.0_wp * ibit2 * adv_sca_5                       &
                       +              ibit1 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j,i+2) - sk(k,j,i-1) )                    &
                   +      (           ibit2 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j,i+3) - sk(k,j,i-2) )                    &
                                             )

                ibit5 = IBITS(wall_flags_0(k,j,i),5,1)
                ibit4 = IBITS(wall_flags_0(k,j,i),4,1)
                ibit3 = IBITS(wall_flags_0(k,j,i),3,1)

                v_comp    = v(k,j+1,i) - v_gtrans
                flux_n(k) = v_comp * (                                        &
                          ( 37.0_wp * ibit5 * adv_sca_5                       &
                       +     7.0_wp * ibit4 * adv_sca_3                       &
                       +              ibit3 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j+1,i) + sk(k,j,i)   )                    &
                   -      (  8.0_wp * ibit5 * adv_sca_5                       &
                       +              ibit4 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j+2,i) + sk(k,j-1,i) )                    &
                   +      (           ibit5 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j+3,i) + sk(k,j-2,i) )                    &
                                     )

                diss_n(k) = -ABS( v_comp ) * (                                &
                          ( 10.0_wp * ibit5 * adv_sca_5                       &
                       +     3.0_wp * ibit4 * adv_sca_3                       &
                       +              ibit3 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j+1,i) - sk(k,j,i)    )                   &
                   -      (  5.0_wp * ibit5 * adv_sca_5                       &
                       +              ibit4 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j+2,i) - sk(k,j-1,i)  )                   &
                   +      (           ibit5 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j+3,i) - sk(k,j-2,i) )                    &
                                             )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit8 = IBITS(wall_flags_0(k,j,i),8,1)
                ibit7 = IBITS(wall_flags_0(k,j,i),7,1)
                ibit6 = IBITS(wall_flags_0(k,j,i),6,1)

                k_ppp = k + 3 * ibit8
                k_pp  = k + 2 * ( 1 - ibit6  )
                k_mm  = k - 2 * ibit8


                flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                      &
                           ( 37.0_wp * ibit8 * adv_sca_5                      &
                        +     7.0_wp * ibit7 * adv_sca_3                      &
                        +           ibit6 * adv_sca_1                         &
                           ) *                                                &
                                   ( sk(k+1,j,i)  + sk(k,j,i)    )            &
                    -      (  8.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i) + sk(k-1,j,i)  )            &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *     ( sk(k_ppp,j,i)+ sk(k_mm,j,i) )            &
                                       )

                diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (              &
                           ( 10.0_wp * ibit8 * adv_sca_5                      &
                        +     3.0_wp * ibit7 * adv_sca_3                      &
                        +              ibit6 * adv_sca_1                      &
                           ) *                                                &
                                   ( sk(k+1,j,i)   - sk(k,j,i)    )           &
                    -      (  5.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i)  - sk(k-1,j,i)  )           &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *                                                &
                                   ( sk(k_ppp,j,i) - sk(k_mm,j,i) )           &
                                               )
!
!--             Apply monotonic adjustment.
                IF ( monotonic_adjustment )  THEN
!
!--                At first, calculate first order fluxes.
                   u_comp = u(k,j,i) - u_gtrans
                   fl_1   =  ( u_comp   * ( sk(k,j,i) + sk(k,j,i-1) )         &
                         -ABS( u_comp ) * ( sk(k,j,i) - sk(k,j,i-1) )         &
                             ) * adv_sca_1 

                   u_comp = u(k,j,i+1) - u_gtrans
                   fr_1   =  ( u_comp   * ( sk(k,j,i+1) + sk(k,j,i) )         &
                         -ABS( u_comp ) * ( sk(k,j,i+1) - sk(k,j,i) )         &
                             ) * adv_sca_1 

                   v_comp = v(k,j,i) - v_gtrans
                   fs_1   =  ( v_comp   * ( sk(k,j,i) + sk(k,j-1,i) )         &
                         -ABS( v_comp ) * ( sk(k,j,i) - sk(k,j-1,i) )         &
                             ) * adv_sca_1 

                   v_comp = v(k,j+1,i) - v_gtrans
                   fn_1   =  ( v_comp   * ( sk(k,j+1,i) + sk(k,j,i) )         &
                         -ABS( v_comp ) * ( sk(k,j+1,i) - sk(k,j,i) )         &
                             ) * adv_sca_1 

                   fd_1   = (  w(k-1,j,i)   * ( sk(k,j,i) + sk(k-1,j,i) )     &
                         -ABS( w(k-1,j,i) ) * ( sk(k,j,i) - sk(k-1,j,i) )     &
                            ) * adv_sca_1 * rho_air_zw(k)

                   ft_1   = (  w(k,j,i)   * ( sk(k+1,j,i) + sk(k,j,i) )       &
                         -ABS( w(k,j,i) ) * ( sk(k+1,j,i) - sk(k,j,i) )       &
                            ) * adv_sca_1 * rho_air_zw(k)
!
!--                Calculate ratio of upwind gradients. Note, Min/Max is just
!--                to avoid if statements.
                   rl     = ( MAX( 0.0_wp, u(k,j,i) - u_gtrans ) *            & 
                               ABS( ( sk(k,j,i-1) - sk(k,j,i-2)            ) /&
                                    ( sk(k,j,i)   - sk(k,j,i-1) + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, u(k,j,i) - u_gtrans ) *            &
                               ABS( ( sk(k,j,i)   - sk(k,j,i+1)            ) /&
                                    ( sk(k,j,i-1) - sk(k,j,i)   + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( u(k,j,i) - u_gtrans + 1E-20_wp )

                   rr     = ( MAX( 0.0_wp, u(k,j,i+1) - u_gtrans ) *          & 
                               ABS( ( sk(k,j,i)   - sk(k,j,i-1)            ) /&
                                    ( sk(k,j,i+1) - sk(k,j,i)   + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, u(k,j,i+1) - u_gtrans ) *          &
                               ABS( ( sk(k,j,i+1) - sk(k,j,i+2)            ) /&
                                    ( sk(k,j,i)   - sk(k,j,i+1) + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( u(k,j,i+1) - u_gtrans + 1E-20_wp )

                   rs     = ( MAX( 0.0_wp, v(k,j,i) - v_gtrans ) *            & 
                               ABS( ( sk(k,j-1,i) - sk(k,j-2,i)            ) /&
                                    ( sk(k,j,i)   - sk(k,j-1,i) + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, v(k,j,i) - v_gtrans ) *            &
                               ABS( ( sk(k,j,i)   - sk(k,j+1,i)            ) /&
                                    ( sk(k,j-1,i) - sk(k,j,i)   + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( v(k,j,i) - v_gtrans + 1E-20_wp )

                   rn     = ( MAX( 0.0_wp, v(k,j+1,i) - v_gtrans ) *          & 
                               ABS( ( sk(k,j,i)   - sk(k,j-1,i)            ) /&
                                    ( sk(k,j+1,i) - sk(k,j,i)   + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, v(k,j+1,i) - v_gtrans ) *          &
                               ABS( ( sk(k,j+1,i) - sk(k,j+2,i)            ) /&
                                    ( sk(k,j,i)   - sk(k,j+1,i) + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( v(k,j+1,i) - v_gtrans + 1E-20_wp )     
!   
!--                Reuse k_mm and compute k_mmm for the vertical gradient ratios. 
!--                Note, for vertical advection below the third grid point above 
!--                surface ( or below the model top) rd and rt are set to 0, i.e. 
!--                use of first order scheme is enforced.
                   k_mmm  = k - 3 * ibit8

                   rd     = ( MAX( 0.0_wp, w(k-1,j,i) ) *                     & 
                            ABS( ( sk(k_mm,j,i) - sk(k_mmm,j,i)           ) / &
                                 ( sk(k-1,j,i)  - sk(k_mm,j,i) + 1E-20_wp )   &
                               ) +                                            & 
                              MIN( 0.0_wp, w(k-1,j,i) ) *                     &
                            ABS( ( sk(k-1,j,i) - sk(k,j,i)            ) /     &
                                 ( sk(k_mm,j,i) - sk(k-1,j,i)   + 1E-20_wp )  &
                               )                                              &
                            ) * ibit8 / ABS( w(k-1,j,i) + 1E-20_wp ) 
 
                   rt     = ( MAX( 0.0_wp, w(k,j,i) ) *                       & 
                            ABS( ( sk(k,j,i)   - sk(k-1,j,i)            ) /   &
                                 ( sk(k+1,j,i) - sk(k,j,i)   + 1E-20_wp )     &
                               ) +                                            & 
                              MIN( 0.0_wp, w(k,j,i) ) *                       &
                            ABS( ( sk(k+1,j,i) - sk(k_pp,j,i)           ) /   &
                                 ( sk(k,j,i)   - sk(k+1,j,i) + 1E-20_wp )     &
                               )                                              &
                            ) * ibit8 / ABS( w(k,j,i) + 1E-20_wp )
!
!--                Calculate empirical limiter function (van Albada2 limiter).
                   phi_l = MIN( 1.0_wp, ( 2.0_wp * ABS( rl ) ) /              &
                                        ( rl**2 + 1.0_wp ) ) 
                   phi_r = MIN( 1.0_wp, ( 2.0_wp * ABS( rr ) ) /              &
                                        ( rr**2 + 1.0_wp ) ) 
                   phi_s = MIN( 1.0_wp, ( 2.0_wp * ABS( rs ) ) /              &
                                        ( rs**2 + 1.0_wp ) ) 
                   phi_n = MIN( 1.0_wp, ( 2.0_wp * ABS( rn ) ) /              &
                                        ( rn**2 + 1.0_wp ) ) 
                   phi_d = MIN( 1.0_wp, ( 2.0_wp * ABS( rd ) ) /              &
                                        ( rd**2 + 1.0_wp ) ) 
                   phi_t = MIN( 1.0_wp, ( 2.0_wp * ABS( rt ) ) /              &
                                        ( rt**2 + 1.0_wp ) ) 
!
!--                Calculate the resulting monotone flux. 
                   swap_flux_x_local(k,j) = fl_1 - phi_l *                    &
                                          ( fl_1 - swap_flux_x_local(k,j) )
                   flux_r(k)              = fr_1 - phi_r *                    &
                                          ( fr_1 - flux_r(k)              )
                   swap_flux_y_local(k)   = fs_1 - phi_s *                    &
                                          ( fs_1 - swap_flux_y_local(k)   )
                   flux_n(k)              = fn_1 - phi_n *                    &
                                          ( fn_1 - flux_n(k)              )
                   flux_d                 = fd_1 - phi_d *                    &
                                          ( fd_1 - flux_d                 )
                   flux_t(k)              = ft_1 - phi_t *                    &
                                          ( ft_1 - flux_t(k)              )
!
!--                Moreover, modify dissipation flux according to the limiter. 
                   swap_diss_x_local(k,j) = swap_diss_x_local(k,j) * phi_l
                   diss_r(k)              = diss_r(k)              * phi_r
                   swap_diss_y_local(k)   = swap_diss_y_local(k)   * phi_s
                   diss_n(k)              = diss_n(k)              * phi_n
                   diss_d                 = diss_d                 * phi_d
                   diss_t(k)              = diss_t(k)              * phi_t

                ENDIF
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div   =   ( u(k,j,i+1) * ( ibit0 + ibit1 + ibit2 )             &
                          - u(k,j,i)   * ( IBITS(wall_flags_0(k,j,i-1),0,1)    &
                                         + IBITS(wall_flags_0(k,j,i-1),1,1)    &
                                         + IBITS(wall_flags_0(k,j,i-1),2,1)    &
                                         )                                     &
                          ) * rho_air(k) * ddx                                 &
                        + ( v(k,j+1,i) * ( ibit3 + ibit4 + ibit5 )             &
                          - v(k,j,i)   * ( IBITS(wall_flags_0(k,j-1,i),3,1)    &
                                         + IBITS(wall_flags_0(k,j-1,i),4,1)    &
                                         + IBITS(wall_flags_0(k,j-1,i),5,1)    &
                                         )                                     &
                          ) * rho_air(k) * ddy                                 &
                        + ( w(k,j,i) * rho_air_zw(k) *                         &
                                         ( ibit6 + ibit7 + ibit8 )             &
                          - w(k-1,j,i) * rho_air_zw(k-1) *                     &
                                         ( IBITS(wall_flags_0(k-1,j,i),6,1)    &
                                         + IBITS(wall_flags_0(k-1,j,i),7,1)    &
                                         + IBITS(wall_flags_0(k-1,j,i),8,1)    &
                                         )                                     &      
                          ) * ddzw(k)


                tend(k,j,i) = tend(k,j,i) - (                                 &
                        ( flux_r(k) + diss_r(k) - swap_flux_x_local(k,j) -    &
                          swap_diss_x_local(k,j)            ) * ddx           &
                      + ( flux_n(k) + diss_n(k) - swap_flux_y_local(k)   -    &
                          swap_diss_y_local(k)              ) * ddy           &
                      + ( ( flux_t(k) + diss_t(k) ) -                         &
                          ( flux_d    + diss_d    )                           &
                                                    ) * drho_air(k) * ddzw(k) &
                                            ) + sk(k,j,i) * div

                swap_flux_y_local(k)   = flux_n(k)
                swap_diss_y_local(k)   = diss_n(k)
                swap_flux_x_local(k,j) = flux_r(k)
                swap_diss_x_local(k,j) = diss_r(k)
                flux_d                 = flux_t(k)
                diss_d                 = diss_t(k)

             ENDDO

             DO  k = nzb_max+1, nzt

                u_comp    = u(k,j,i+1) - u_gtrans
                flux_r(k) = u_comp * (                                        &
                      37.0_wp * ( sk(k,j,i+1) + sk(k,j,i)   )                 &
                    -  8.0_wp * ( sk(k,j,i+2) + sk(k,j,i-1) )                 &
                    +           ( sk(k,j,i+3) + sk(k,j,i-2) ) ) * adv_sca_5
                diss_r(k) = -ABS( u_comp ) * (                                &
                      10.0_wp * ( sk(k,j,i+1) - sk(k,j,i)   )                 &
                    -  5.0_wp * ( sk(k,j,i+2) - sk(k,j,i-1) )                 &
                    +           ( sk(k,j,i+3) - sk(k,j,i-2) ) ) * adv_sca_5

                v_comp    = v(k,j+1,i) - v_gtrans
                flux_n(k) = v_comp * (                                        &
                      37.0_wp * ( sk(k,j+1,i) + sk(k,j,i)   )                 &
                    -  8.0_wp * ( sk(k,j+2,i) + sk(k,j-1,i) )                 &
                    +           ( sk(k,j+3,i) + sk(k,j-2,i) ) ) * adv_sca_5
                diss_n(k) = -ABS( v_comp ) * (                                &
                      10.0_wp * ( sk(k,j+1,i) - sk(k,j,i)   )                 &
                    -  5.0_wp * ( sk(k,j+2,i) - sk(k,j-1,i) )                 &
                    +           ( sk(k,j+3,i) - sk(k,j-2,i) ) ) * adv_sca_5
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit8 = IBITS(wall_flags_0(k,j,i),8,1)
                ibit7 = IBITS(wall_flags_0(k,j,i),7,1)
                ibit6 = IBITS(wall_flags_0(k,j,i),6,1)

                k_ppp = k + 3 * ibit8
                k_pp  = k + 2 * ( 1 - ibit6  )
                k_mm  = k - 2 * ibit8


                flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                      &
                           ( 37.0_wp * ibit8 * adv_sca_5                      &
                        +     7.0_wp * ibit7 * adv_sca_3                      &
                        +              ibit6 * adv_sca_1                      &
                           ) *                                                &
                                   ( sk(k+1,j,i)  + sk(k,j,i)     )           &
                    -      (  8.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i) + sk(k-1,j,i)   )           &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *     ( sk(k_ppp,j,i)+ sk(k_mm,j,i)  )           &
                                       )

                diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (              &
                           ( 10.0_wp * ibit8 * adv_sca_5                      &
                        +     3.0_wp * ibit7 * adv_sca_3                      &
                        +              ibit6 * adv_sca_1                      &
                           ) *                                                &
                                   ( sk(k+1,j,i)   - sk(k,j,i)    )           &
                    -      (  5.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i)  - sk(k-1,j,i)  )           &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *                                                &
                                   ( sk(k_ppp,j,i) - sk(k_mm,j,i) )           &
                                               )
!
!--             Apply monotonic adjustment.
                IF ( monotonic_adjustment )  THEN
!
!--                At first, calculate first order fluxes.
                   u_comp = u(k,j,i) - u_gtrans
                   fl_1   =  ( u_comp   * ( sk(k,j,i) + sk(k,j,i-1) )         &
                         -ABS( u_comp ) * ( sk(k,j,i) - sk(k,j,i-1) )         &
                             ) * adv_sca_1 

                   u_comp = u(k,j,i+1) - u_gtrans
                   fr_1   =  ( u_comp   * ( sk(k,j,i+1) + sk(k,j,i) )         &
                         -ABS( u_comp ) * ( sk(k,j,i+1) - sk(k,j,i) )         &
                             ) * adv_sca_1 

                   v_comp = v(k,j,i) - v_gtrans
                   fs_1   =  ( v_comp   * ( sk(k,j,i) + sk(k,j-1,i) )         &
                         -ABS( v_comp ) * ( sk(k,j,i) - sk(k,j-1,i) )         &
                             ) * adv_sca_1 

                   v_comp = v(k,j+1,i) - v_gtrans
                   fn_1   =  ( v_comp   * ( sk(k,j+1,i) + sk(k,j,i) )         &
                         -ABS( v_comp ) * ( sk(k,j+1,i) - sk(k,j,i) )         &
                             ) * adv_sca_1 

                   fd_1   = (  w(k-1,j,i)   * ( sk(k,j,i) + sk(k-1,j,i) )     &
                         -ABS( w(k-1,j,i) ) * ( sk(k,j,i) - sk(k-1,j,i) )     &
                            ) * adv_sca_1 * rho_air_zw(k)

                   ft_1   = (  w(k,j,i)   * ( sk(k+1,j,i) + sk(k,j,i) )       &
                         -ABS( w(k,j,i) ) * ( sk(k+1,j,i) - sk(k,j,i) )       &
                            ) * adv_sca_1 * rho_air_zw(k)
!
!--                Calculate ratio of upwind gradients. Note, Min/Max is just
!--                to avoid if statements.
                   rl     = ( MAX( 0.0_wp, u(k,j,i) - u_gtrans ) *            & 
                               ABS( ( sk(k,j,i-1) - sk(k,j,i-2)            ) /&
                                    ( sk(k,j,i)   - sk(k,j,i-1) + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, u(k,j,i) - u_gtrans ) *            &
                               ABS( ( sk(k,j,i)   - sk(k,j,i+1)            ) /&
                                    ( sk(k,j,i-1) - sk(k,j,i)   + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( u(k,j,i) - u_gtrans + 1E-20_wp )

                   rr     = ( MAX( 0.0_wp, u(k,j,i+1) - u_gtrans ) *          & 
                               ABS( ( sk(k,j,i)   - sk(k,j,i-1)            ) /&
                                    ( sk(k,j,i+1) - sk(k,j,i)   + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, u(k,j,i+1) - u_gtrans ) *          &
                               ABS( ( sk(k,j,i+1) - sk(k,j,i+2)            ) /&
                                    ( sk(k,j,i)   - sk(k,j,i+1) + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( u(k,j,i+1) - u_gtrans + 1E-20_wp )

                   rs     = ( MAX( 0.0_wp, v(k,j,i) - v_gtrans ) *            & 
                               ABS( ( sk(k,j-1,i) - sk(k,j-2,i)            ) /&
                                    ( sk(k,j,i)   - sk(k,j-1,i) + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, v(k,j,i) - v_gtrans ) *            &
                               ABS( ( sk(k,j,i)   - sk(k,j+1,i)            ) /&
                                    ( sk(k,j-1,i) - sk(k,j,i)   + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( v(k,j,i) - v_gtrans + 1E-20_wp )

                   rn     = ( MAX( 0.0_wp, v(k,j+1,i) - v_gtrans ) *          & 
                               ABS( ( sk(k,j,i)   - sk(k,j-1,i)            ) /&
                                    ( sk(k,j+1,i) - sk(k,j,i)   + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, v(k,j+1,i) - v_gtrans ) *          &
                               ABS( ( sk(k,j+1,i) - sk(k,j+2,i)            ) /&
                                    ( sk(k,j,i)   - sk(k,j+1,i) + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( v(k,j+1,i) - v_gtrans + 1E-20_wp )     
!   
!--                Reuse k_mm and compute k_mmm for the vertical gradient ratios. 
!--                Note, for vertical advection below the third grid point above 
!--                surface ( or below the model top) rd and rt are set to 0, i.e. 
!--                use of first order scheme is enforced.
                   k_mmm  = k - 3 * ibit8

                   rd     = ( MAX( 0.0_wp, w(k-1,j,i) ) *                     & 
                            ABS( ( sk(k_mm,j,i) - sk(k_mmm,j,i)           ) / &
                                 ( sk(k-1,j,i)  - sk(k_mm,j,i) + 1E-20_wp )   &
                               ) +                                            & 
                              MIN( 0.0_wp, w(k-1,j,i) ) *                     &
                            ABS( ( sk(k-1,j,i) - sk(k,j,i)            ) /     &
                                 ( sk(k_mm,j,i) - sk(k-1,j,i)   + 1E-20_wp )  &
                               )                                              &
                            ) * ibit8 / ABS( w(k-1,j,i) + 1E-20_wp ) 
 
                   rt     = ( MAX( 0.0_wp, w(k,j,i) ) *                       & 
                            ABS( ( sk(k,j,i)   - sk(k-1,j,i)            ) /   &
                                 ( sk(k+1,j,i) - sk(k,j,i)   + 1E-20_wp )     &
                               ) +                                            & 
                              MIN( 0.0_wp, w(k,j,i) ) *                       &
                            ABS( ( sk(k+1,j,i) - sk(k_pp,j,i)           ) /   &
                                 ( sk(k,j,i)   - sk(k+1,j,i) + 1E-20_wp )     &
                               )                                              &
                            ) * ibit8 / ABS( w(k,j,i) + 1E-20_wp ) 
!
!--                Calculate empirical limiter function (van Albada2 limiter).
                   phi_l = MIN( 1.0_wp, ( 2.0_wp * ABS( rl ) ) /              &
                                        ( rl**2 + 1.0_wp ) ) 
                   phi_r = MIN( 1.0_wp, ( 2.0_wp * ABS( rr ) ) /              &
                                        ( rr**2 + 1.0_wp ) ) 
                   phi_s = MIN( 1.0_wp, ( 2.0_wp * ABS( rs ) ) /              &
                                        ( rs**2 + 1.0_wp ) ) 
                   phi_n = MIN( 1.0_wp, ( 2.0_wp * ABS( rn ) ) /              &
                                        ( rn**2 + 1.0_wp ) ) 
                   phi_d = MIN( 1.0_wp, ( 2.0_wp * ABS( rd ) ) /              &
                                        ( rd**2 + 1.0_wp ) ) 
                   phi_t = MIN( 1.0_wp, ( 2.0_wp * ABS( rt ) ) /              &
                                        ( rt**2 + 1.0_wp ) ) 
!
!--                Calculate the resulting monotone flux. 
                   swap_flux_x_local(k,j) = fl_1 - phi_l *                    &
                                          ( fl_1 - swap_flux_x_local(k,j) )
                   flux_r(k)              = fr_1 - phi_r *                    &
                                          ( fr_1 - flux_r(k)              )
                   swap_flux_y_local(k)   = fs_1 - phi_s *                    &
                                          ( fs_1 - swap_flux_y_local(k)   )
                   flux_n(k)              = fn_1 - phi_n *                    &
                                          ( fn_1 - flux_n(k)              )
                   flux_d                 = fd_1 - phi_d *                    &
                                          ( fd_1 - flux_d                 )
                   flux_t(k)              = ft_1 - phi_t *                    &
                                          ( ft_1 - flux_t(k)              )
!
!--                Moreover, modify dissipation flux according to the limiter. 
                   swap_diss_x_local(k,j) = swap_diss_x_local(k,j) * phi_l
                   diss_r(k)              = diss_r(k)              * phi_r
                   swap_diss_y_local(k)   = swap_diss_y_local(k)   * phi_s
                   diss_n(k)              = diss_n(k)              * phi_n
                   diss_d                 = diss_d                 * phi_d
                   diss_t(k)              = diss_t(k)              * phi_t

                ENDIF
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities introduced
!--             by a not sufficient reduction of divergences near topography.
                div         =   ( u(k,j,i+1) - u(k,j,i)   ) * rho_air(k) * ddx &
                              + ( v(k,j+1,i) - v(k,j,i)   ) * rho_air(k) * ddy &
                              + ( w(k,j,i)   * rho_air_zw(k) -                 &
                                  w(k-1,j,i) * rho_air_zw(k-1) ) * ddzw(k)

                tend(k,j,i) = tend(k,j,i) - (                                 &
                        ( flux_r(k) + diss_r(k) - swap_flux_x_local(k,j) -    &
                          swap_diss_x_local(k,j)            ) * ddx           &
                      + ( flux_n(k) + diss_n(k) - swap_flux_y_local(k)   -    &
                          swap_diss_y_local(k)              ) * ddy           &
                      + ( ( flux_t(k) + diss_t(k) ) -                         &
                          ( flux_d    + diss_d    )                           &
                                                    ) * drho_air(k) * ddzw(k) &
                                            ) + sk(k,j,i) * div

                swap_flux_y_local(k)   = flux_n(k)
                swap_diss_y_local(k)   = diss_n(k)
                swap_flux_x_local(k,j) = flux_r(k)
                swap_diss_x_local(k,j) = diss_r(k)
                flux_d                 = flux_t(k)
                diss_d                 = diss_t(k)

             ENDDO
!
!--          Evaluation of statistics. 
             SELECT CASE ( sk_char )

                 CASE ( 'pt' )
                    DO  k = nzb, nzt
                       sums_wspts_ws_l(k,tn) = sums_wspts_ws_l(k,tn)           &
                          + ( flux_t(k)                                        &
                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                            + diss_t(k)                                        &
                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )  &
                            ) * weight_substep(intermediate_timestep_count)
                    ENDDO
                 CASE ( 'sa' )
                    DO  k = nzb, nzt
                       sums_wssas_ws_l(k,tn) = sums_wssas_ws_l(k,tn)           &
                          + ( flux_t(k)                                        &
                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                            + diss_t(k)                                        &
                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )  &
                            ) * weight_substep(intermediate_timestep_count)
                    ENDDO
                 CASE ( 'q' )
                    DO  k = nzb, nzt
                       sums_wsqs_ws_l(k,tn)  = sums_wsqs_ws_l(k,tn)            &
                          + ( flux_t(k)                                        &
                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                            + diss_t(k)                                        &
                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )  &
                            ) * weight_substep(intermediate_timestep_count)
                    ENDDO
                 CASE ( 'qr' )
                    DO  k = nzb, nzt
                       sums_wsqrs_ws_l(k,tn)  = sums_wsqrs_ws_l(k,tn)          &
                          + ( flux_t(k)                                        &
                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                            + diss_t(k)                                        &
                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )  &
                            ) * weight_substep(intermediate_timestep_count)
                    ENDDO
                 CASE ( 'nr' )
                    DO  k = nzb, nzt
                       sums_wsnrs_ws_l(k,tn)  = sums_wsnrs_ws_l(k,tn)          &
                          + ( flux_t(k)                                        &
                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                            + diss_t(k)                                        &
                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )  &
                            ) * weight_substep(intermediate_timestep_count)
                    ENDDO
                 CASE ( 's' )
                    DO  k = nzb, nzt
                       sums_wsss_ws_l(k,tn)  = sums_wsss_ws_l(k,tn)            &
                          + ( flux_t(k)                                        &
                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3,0)                 )  &
                            + diss_t(k)                                        &
                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )  &
                            ) * weight_substep(intermediate_timestep_count)
                    ENDDO   
                                   

              END SELECT

         ENDDO
      ENDDO

    END SUBROUTINE advec_s_ws


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Scalar advection - Call for all grid points - accelerator version
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_ws_acc ( sk, sk_char )

       USE arrays_3d,                                                         &
           ONLY:  ddzw, drho_air, tend, u, v, w, rho_air, rho_air_zw

       USE constants,                                                         &
           ONLY:  adv_sca_1, adv_sca_3, adv_sca_5

       USE control_parameters,                                                &
           ONLY:  intermediate_timestep_count, monotonic_adjustment, u_gtrans,&
                  v_gtrans

       USE grid_variables,                                                    &
           ONLY:  ddx, ddy

       USE indices,                                                           &
           ONLY:  i_left, i_right, j_north, j_south, nxlg, nxrg, nyng, nysg,  &
                  nzb, nzb_max, nzt, wall_flags_0

       USE kinds
       
!        USE statistics,                                                       &
!            ONLY:  sums_wspts_ws_l, sums_wsqs_ws_l, sums_wssas_ws_l,          &
!                   sums_wsqrs_ws_l, sums_wsnrs_ws_l, weight_substep

       IMPLICIT NONE

       CHARACTER (LEN = *), INTENT(IN)    :: sk_char !<

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit0  !<
       INTEGER(iwp) ::  ibit1  !<
       INTEGER(iwp) ::  ibit2  !<
       INTEGER(iwp) ::  ibit3  !<
       INTEGER(iwp) ::  ibit4  !<
       INTEGER(iwp) ::  ibit5  !<
       INTEGER(iwp) ::  ibit6  !<
       INTEGER(iwp) ::  ibit7  !<
       INTEGER(iwp) ::  ibit8  !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_mmm  !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<

       REAL(wp)    ::  diss_d !<
       REAL(wp)    ::  diss_l !<
       REAL(wp)    ::  diss_n !<
       REAL(wp)    ::  diss_r !<
       REAL(wp)    ::  diss_s !<
       REAL(wp)    ::  diss_t !<
       REAL(wp)    ::  div    !<
       REAL(wp)    ::  flux_d !<
       REAL(wp)    ::  flux_l !<
       REAL(wp)    ::  flux_n !<
       REAL(wp)    ::  flux_r !<
       REAL(wp)    ::  flux_s !<
       REAL(wp)    ::  flux_t !<
       REAL(wp)    ::  fd_1   !<
       REAL(wp)    ::  fl_1   !<
       REAL(wp)    ::  fn_1   !<
       REAL(wp)    ::  fr_1   !<
       REAL(wp)    ::  fs_1   !<
       REAL(wp)    ::  ft_1   !<
       REAL(wp)    ::  phi_d  !<
       REAL(wp)    ::  phi_l  !<
       REAL(wp)    ::  phi_n  !<
       REAL(wp)    ::  phi_r  !<
       REAL(wp)    ::  phi_s  !<
       REAL(wp)    ::  phi_t  !<
       REAL(wp)    ::  rd     !<
       REAL(wp)    ::  rl     !<
       REAL(wp)    ::  rn     !<
       REAL(wp)    ::  rr     !<
       REAL(wp)    ::  rs     !<
       REAL(wp)    ::  rt     !<
       REAL(wp)    ::  u_comp !<
       REAL(wp)    ::  v_comp !<

       REAL(wp), INTENT(IN), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  ::  sk !<

!
!--    Computation of fluxes and tendency terms
       !$acc kernels present( ddzw, sk, tend, u, v, w, wall_flags_0 )
       DO  i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = nzb+1, nzt

                ibit2 = IBITS(wall_flags_0(k,j,i-1),2,1)
                ibit1 = IBITS(wall_flags_0(k,j,i-1),1,1)
                ibit0 = IBITS(wall_flags_0(k,j,i-1),0,1)

                u_comp              = u(k,j,i) - u_gtrans
                flux_l              = u_comp * (                              &
                                               ( 37.0_wp * ibit2 * adv_sca_5  &
                                            +     7.0_wp * ibit1 * adv_sca_3  &
                                            +              ibit0 * adv_sca_1  &
                                               ) *                            &
                                         ( sk(k,j,i)   + sk(k,j,i-1)    )     &
                                        -      (  8.0_wp * ibit2 * adv_sca_5  &
                                            +              ibit1 * adv_sca_3  &
                                               ) *                            &
                                         ( sk(k,j,i+1) + sk(k,j,i-2)    )     &
                                        +      (           ibit2 * adv_sca_5  &
                                               ) *                            &
                                         ( sk(k,j,i+2) + sk(k,j,i-3)    )     &
                                               )

                diss_l              = -ABS( u_comp ) * (                      &
                                               ( 10.0_wp * ibit2 * adv_sca_5  &
                                            +     3.0_wp * ibit1 * adv_sca_3  &
                                            +              ibit0 * adv_sca_1  &
                                               ) *                            &
                                         ( sk(k,j,i)   - sk(k,j,i-1)    )     &
                                        -      (  5.0_wp * ibit2 * adv_sca_5  &
                                            +              ibit1 * adv_sca_3  &
                                               ) *                            &
                                         ( sk(k,j,i+1) - sk(k,j,i-2)    )     &
                                        +      (           ibit2 * adv_sca_5  &
                                               ) *                            &
                                         ( sk(k,j,i+2) - sk(k,j,i-3)    )     &
                                                        )

                ibit2 = IBITS(wall_flags_0(k,j,i),2,1)
                ibit1 = IBITS(wall_flags_0(k,j,i),1,1)
                ibit0 = IBITS(wall_flags_0(k,j,i),0,1)

                u_comp    = u(k,j,i+1) - u_gtrans
                flux_r    = u_comp * (                                        &
                          ( 37.0_wp * ibit2 * adv_sca_5                       &
                      +      7.0_wp * ibit1 * adv_sca_3                       &
                      +               ibit0 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j,i+1) + sk(k,j,i)   )                    &
                   -      (  8.0_wp * ibit2 * adv_sca_5                       &
                       +              ibit1 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j,i+2) + sk(k,j,i-1) )                    &
                   +      (           ibit2 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j,i+3) + sk(k,j,i-2) )                    &
                                     )

                diss_r    = -ABS( u_comp ) * (                                &
                          ( 10.0_wp * ibit2 * adv_sca_5                       &
                       +     3.0_wp * ibit1 * adv_sca_3                       &
                       +              ibit0 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j,i+1) - sk(k,j,i)   )                    &
                   -      (  5.0_wp * ibit2 * adv_sca_5                       &
                       +              ibit1 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j,i+2) - sk(k,j,i-1) )                    &
                   +      (           ibit2 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j,i+3) - sk(k,j,i-2) )                    &
                                             )

                ibit5 = IBITS(wall_flags_0(k,j-1,i),5,1)
                ibit4 = IBITS(wall_flags_0(k,j-1,i),4,1)
                ibit3 = IBITS(wall_flags_0(k,j-1,i),3,1)

                v_comp    = v(k,j,i) - v_gtrans
                flux_s    = v_comp * (                                        &
                          ( 37.0_wp * ibit5 * adv_sca_5                       &
                       +     7.0_wp * ibit4 * adv_sca_3                       &
                       +              ibit3 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j,i)  + sk(k,j-1,i)     )                 &
                    -     (  8.0_wp * ibit5 * adv_sca_5                       &
                       +              ibit4 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j+1,i) + sk(k,j-2,i)    )                 &
                   +      (           ibit5 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j+2,i) + sk(k,j-3,i)    )                 &
                                     )

                diss_s    = -ABS( v_comp ) * (                                &
                          ( 10.0_wp * ibit5 * adv_sca_5                       &
                       +     3.0_wp * ibit4 * adv_sca_3                       &
                       +              ibit3 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j,i)   - sk(k,j-1,i)  )                   &
                   -      (  5.0_wp * ibit5 * adv_sca_5                       &
                       +              ibit4 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j+1,i) - sk(k,j-2,i)  )                   &
                   +      (           ibit5 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j+2,i) - sk(k,j-3,i)  )                   &
                                             )

                ibit5 = IBITS(wall_flags_0(k,j,i),5,1)
                ibit4 = IBITS(wall_flags_0(k,j,i),4,1)
                ibit3 = IBITS(wall_flags_0(k,j,i),3,1)

                v_comp    = v(k,j+1,i) - v_gtrans
                flux_n    = v_comp * (                                        &
                          ( 37.0_wp * ibit5 * adv_sca_5                       &
                       +     7.0_wp * ibit4 * adv_sca_3                       &
                       +              ibit3 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j+1,i) + sk(k,j,i)   )                    &
                   -      (  8.0_wp * ibit5 * adv_sca_5                       &
                       +              ibit4 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j+2,i) + sk(k,j-1,i) )                    &
                   +      (           ibit5 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j+3,i) + sk(k,j-2,i) )                    &
                                     )

                diss_n    = -ABS( v_comp ) * (                                &
                          ( 10.0_wp * ibit5 * adv_sca_5                       &
                       +     3.0_wp * ibit4 * adv_sca_3                       &
                       +              ibit3 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j+1,i) - sk(k,j,i)    )                   &
                   -      (  5.0_wp * ibit5 * adv_sca_5                       &
                       +              ibit4 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j+2,i) - sk(k,j-1,i)  )                   &
                   +      (           ibit5 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j+3,i) - sk(k,j-2,i)  )                   &
                                             )

!
!--             indizes k_m, k_mm, ... should be known at these point
                ibit8 = IBITS(wall_flags_0(k-1,j,i),8,1)
                ibit7 = IBITS(wall_flags_0(k-1,j,i),7,1)
                ibit6 = IBITS(wall_flags_0(k-1,j,i),6,1)

                k_pp  = k + 2 * ibit8
                k_mm  = k - 2 * ( ibit7 + ibit8 )
                k_mmm = k - 3 * ibit8

                flux_d    = w(k-1,j,i) * (                                    &
                           ( 37.0_wp * ibit8 * adv_sca_5                      &
                        +     7.0_wp * ibit7 * adv_sca_3                      &
                        +              ibit6 * adv_sca_1                      &
                           ) *                                                &
                                   ( sk(k,j,i)    + sk(k-1,j,i)  )            &
                    -      (  8.0_wp * ibit8 * adv_sca_5                      &
                          +            ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k+1,j,i) + sk(k_mm,j,i)  )            &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *     ( sk(k_pp,j,i)+ sk(k_mmm,j,i) )            &
                                         )

                diss_d    = -ABS( w(k-1,j,i) ) * (                            &
                           ( 10.0_wp * ibit8 * adv_sca_5                      &
                        +     3.0_wp * ibit7 * adv_sca_3                      &
                        +              ibit6 * adv_sca_1                      &
                           ) *                                                &
                                   ( sk(k,j,i)    - sk(k-1,j,i)   )           &
                    -      (  5.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k+1,j,i)  - sk(k_mm,j,i)  )           &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i) - sk(k_mmm,j,i) )           &
                                                 )

                ibit8 = IBITS(wall_flags_0(k,j,i),8,1)
                ibit7 = IBITS(wall_flags_0(k,j,i),7,1)
                ibit6 = IBITS(wall_flags_0(k,j,i),6,1)

                k_ppp = k + 3 * ibit8
                k_pp  = k + 2 * ( 1 - ibit6  )
                k_mm  = k - 2 * ibit8

                flux_t    = w(k,j,i) * rho_air_zw(k) * (                      &
                           ( 37.0_wp * ibit8 * adv_sca_5                      &
                        +     7.0_wp * ibit7 * adv_sca_3                      &
                        +              ibit6 * adv_sca_1                      &
                           ) *                                                &
                                   ( sk(k+1,j,i)  + sk(k,j,i)    )            &
                    -      (  8.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i) + sk(k-1,j,i)  )            &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *     ( sk(k_ppp,j,i)+ sk(k_mm,j,i) )            &
                                       )

                diss_t    = -ABS( w(k,j,i) ) * rho_air_zw(k) * (              &
                           ( 10.0_wp * ibit8 * adv_sca_5                      &
                        +     3.0_wp * ibit7 * adv_sca_3                      &
                        +              ibit6 * adv_sca_1                      &
                           ) *                                                &
                                   ( sk(k+1,j,i)   - sk(k,j,i)    )           &
                    -      (  5.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i)  - sk(k-1,j,i)  )           &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *                                                &
                                   ( sk(k_ppp,j,i) - sk(k_mm,j,i) )           &
                                         )
!
!--             Apply monotonic adjustment.
                IF ( monotonic_adjustment )  THEN
!
!--                At first, calculate first order fluxes.
                   u_comp = u(k,j,i) - u_gtrans
                   fl_1   =  ( u_comp   * ( sk(k,j,i) + sk(k,j,i-1) )         &
                         -ABS( u_comp ) * ( sk(k,j,i) - sk(k,j,i-1) )         &
                             ) * adv_sca_1 

                   u_comp = u(k,j,i+1) - u_gtrans
                   fr_1   =  ( u_comp   * ( sk(k,j,i+1) + sk(k,j,i) )         &
                         -ABS( u_comp ) * ( sk(k,j,i+1) - sk(k,j,i) )         &
                             ) * adv_sca_1 

                   v_comp = v(k,j,i) - v_gtrans
                   fs_1   =  ( v_comp   * ( sk(k,j,i) + sk(k,j-1,i) )         &
                         -ABS( v_comp ) * ( sk(k,j,i) - sk(k,j-1,i) )         &
                             ) * adv_sca_1 

                   v_comp = v(k,j+1,i) - v_gtrans
                   fn_1   =  ( v_comp   * ( sk(k,j+1,i) + sk(k,j,i) )         &
                         -ABS( v_comp ) * ( sk(k,j+1,i) - sk(k,j,i) )         &
                             ) * adv_sca_1 

                   fd_1   = (  w(k-1,j,i)   * ( sk(k,j,i) + sk(k-1,j,i) )     &
                        -ABS( w(k-1,j,i) ) * ( sk(k,j,i) - sk(k-1,j,i) )      &
                            ) * adv_sca_1 * rho_air_zw(k)

                   ft_1   = (  w(k,j,i)   * ( sk(k+1,j,i) + sk(k,j,i) )       &
                        -ABS( w(k,j,i) ) * ( sk(k+1,j,i) - sk(k,j,i) )        &
                            ) * adv_sca_1 * rho_air_zw(k)
!
!--                Calculate ratio of upwind gradients. Note, Min/Max is just
!--                to avoid if statements.
                   rl     = ( MAX( 0.0_wp, u(k,j,i) - u_gtrans ) *            & 
                               ABS( ( sk(k,j,i-1) - sk(k,j,i-2)            ) /&
                                    ( sk(k,j,i)   - sk(k,j,i-1) + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, u(k,j,i) - u_gtrans ) *            &
                               ABS( ( sk(k,j,i)   - sk(k,j,i+1)            ) /&
                                    ( sk(k,j,i-1) - sk(k,j,i)   + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( u(k,j,i) - u_gtrans + 1E-20_wp )

                   rr     = ( MAX( 0.0_wp, u(k,j,i+1) - u_gtrans ) *          & 
                               ABS( ( sk(k,j,i)   - sk(k,j,i-1)            ) /&
                                    ( sk(k,j,i+1) - sk(k,j,i)   + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, u(k,j,i+1) - u_gtrans ) *          &
                               ABS( ( sk(k,j,i+1) - sk(k,j,i+2)            ) /&
                                    ( sk(k,j,i)   - sk(k,j,i+1) + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( u(k,j,i+1) - u_gtrans + 1E-20_wp )

                   rs     = ( MAX( 0.0_wp, v(k,j,i) - v_gtrans ) *            & 
                               ABS( ( sk(k,j-1,i) - sk(k,j-2,i)            ) /&
                                    ( sk(k,j,i)   - sk(k,j-1,i) + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, v(k,j,i) - v_gtrans ) *            &
                               ABS( ( sk(k,j,i)   - sk(k,j+1,i)            ) /&
                                    ( sk(k,j-1,i) - sk(k,j,i)   + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( v(k,j,i) - v_gtrans + 1E-20_wp )

                   rn     = ( MAX( 0.0_wp, v(k,j+1,i) - v_gtrans ) *          & 
                               ABS( ( sk(k,j,i)   - sk(k,j-1,i)            ) /&
                                    ( sk(k,j+1,i) - sk(k,j,i)   + 1E-20_wp )  &
                                  ) +                                         & 
                              MIN( 0.0_wp, v(k,j+1,i) - v_gtrans ) *          &
                               ABS( ( sk(k,j+1,i) - sk(k,j+2,i)            ) /&
                                    ( sk(k,j,i)   - sk(k,j+1,i) + 1E-20_wp )  &
                                  )                                           &
                            ) / ABS( v(k,j+1,i) - v_gtrans + 1E-20_wp )     
!   
!--                Reuse k_mm and compute k_mmm for the vertical gradient ratios. 
!--                Note, for vertical advection below the third grid point above 
!--                surface ( or below the model top) rd and rt are set to 0, i.e. 
!--                use of first order scheme is enforced.
                   k_mmm  = k - 3 * ibit8

                   rd     = ( MAX( 0.0_wp, w(k-1,j,i) ) *                     & 
                            ABS( ( sk(k_mm,j,i) - sk(k_mmm,j,i)           ) / &
                                 ( sk(k-1,j,i)  - sk(k_mm,j,i) + 1E-20_wp )   &
                               ) +                                            & 
                              MIN( 0.0_wp, w(k-1,j,i) ) *                     &
                            ABS( ( sk(k-1,j,i) - sk(k,j,i)            ) /     &
                                 ( sk(k_mm,j,i) - sk(k-1,j,i)   + 1E-20_wp )  &
                               )                                              &
                            ) * ibit8 / ABS( w(k-1,j,i) + 1E-20_wp ) 
 
                   rt     = ( MAX( 0.0_wp, w(k,j,i) ) *                       & 
                            ABS( ( sk(k,j,i)   - sk(k-1,j,i)            ) /   &
                                 ( sk(k+1,j,i) - sk(k,j,i)   + 1E-20_wp )     &
                               ) +                                            & 
                              MIN( 0.0_wp, w(k,j,i) ) *                       &
                            ABS( ( sk(k+1,j,i) - sk(k_pp,j,i)           ) /   &
                                 ( sk(k,j,i)   - sk(k+1,j,i) + 1E-20_wp )     &
                               )                                              &
                            ) * ibit8 / ABS( w(k,j,i) + 1E-20_wp )
!
!--                Calculate empirical limiter function (van Albada2 limiter).
                   phi_l = MIN( 1.0_wp, ( 2.0_wp * ABS( rl ) ) /              &
                                        ( rl**2 + 1.0_wp ) ) 
                   phi_r = MIN( 1.0_wp, ( 2.0_wp * ABS( rr ) ) /              &
                                        ( rr**2 + 1.0_wp ) ) 
                   phi_s = MIN( 1.0_wp, ( 2.0_wp * ABS( rs ) ) /              &
                                        ( rs**2 + 1.0_wp ) ) 
                   phi_n = MIN( 1.0_wp, ( 2.0_wp * ABS( rn ) ) /              &
                                        ( rn**2 + 1.0_wp ) ) 
                   phi_d = MIN( 1.0_wp, ( 2.0_wp * ABS( rd ) ) /              &
                                        ( rd**2 + 1.0_wp ) ) 
                   phi_t = MIN( 1.0_wp, ( 2.0_wp * ABS( rt ) ) /              &
                                        ( rt**2 + 1.0_wp ) ) 
!
!--                Calculate the resulting monotone flux. 
                   flux_l = fl_1 - phi_l * ( fl_1 - flux_l )
                   flux_r = fr_1 - phi_r * ( fr_1 - flux_r )
                   flux_s = fs_1 - phi_s * ( fs_1 - flux_s )
                   flux_n = fn_1 - phi_n * ( fn_1 - flux_n )
                   flux_d = fd_1 - phi_d * ( fd_1 - flux_d )
                   flux_t = ft_1 - phi_t * ( ft_1 - flux_t )
!
!--                Moreover, modify dissipation flux according to the limiter. 
                   diss_l = diss_l * phi_l
                   diss_r = diss_r * phi_r
                   diss_s = diss_s * phi_s
                   diss_n = diss_n * phi_n
                   diss_d = diss_d * phi_d
                   diss_t = diss_t * phi_t

                ENDIF
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div   =   ( u(k,j,i+1) * ( ibit0 + ibit1 + ibit2 )             &
                          - u(k,j,i)   * ( IBITS(wall_flags_0(k,j,i-1),0,1)    &
                                         + IBITS(wall_flags_0(k,j,i-1),1,1)    &
                                         + IBITS(wall_flags_0(k,j,i-1),2,1)    &
                                         )                                     &
                          ) * rho_air(k) * ddx                                 &
                        + ( v(k,j+1,i) * ( ibit3 + ibit4 + ibit5 )             &
                          - v(k,j,i)   * ( IBITS(wall_flags_0(k,j-1,i),3,1)    &
                                         + IBITS(wall_flags_0(k,j-1,i),4,1)    &
                                         + IBITS(wall_flags_0(k,j-1,i),5,1)    &
                                         )                                     &
                          ) * rho_air(k) * ddy                                 &
                        + ( w(k,j,i) * rho_air_zw(k) *                         &
                                         ( ibit6 + ibit7 + ibit8 )             &
                          - w(k-1,j,i) * rho_air_zw(k-1) *                     &
                                         ( IBITS(wall_flags_0(k-1,j,i),6,1)    &
                                         + IBITS(wall_flags_0(k-1,j,i),7,1)    &
                                         + IBITS(wall_flags_0(k-1,j,i),8,1)    &
                                         )                                     &      
                          ) * ddzw(k)


                tend(k,j,i) = - (                                             &
                               ( flux_r + diss_r - flux_l - diss_l ) * ddx    &
                             + ( flux_n + diss_n - flux_s - diss_s ) * ddy    &
                             + ( ( flux_t + diss_t ) -                        &
                                 ( flux_d + diss_d )                          &
                                                    ) * drho_air(k) * ddzw(k) &
                                ) + div * sk(k,j,i)

!++
!--             Evaluation of statistics
!                SELECT CASE ( sk_char )
!
!                   CASE ( 'pt' )
!                      sums_wspts_ws_l(k,tn) = sums_wspts_ws_l(k,tn)         &
!                       + ( flux_t + diss_t )                                &
!                       *   weight_substep(intermediate_timestep_count)
!                   CASE ( 'sa' )
!                      sums_wssas_ws_l(k,tn) = sums_wssas_ws_l(k,tn)         &
!                       + ( flux_t + diss_t )                                &
!                       *   weight_substep(intermediate_timestep_count)
!                   CASE ( 'q' )
!                      sums_wsqs_ws_l(k,tn) = sums_wsqs_ws_l(k,tn)           &
!                      + ( flux_t + diss_t )                                 &
!                      *   weight_substep(intermediate_timestep_count)
!                   CASE ( 'qr' )
!                      sums_wsqrs_ws_l(k,tn) = sums_wsqrs_ws_l(k,tn)         &
!                      + ( flux_t + diss_t )                                 &
!                      *   weight_substep(intermediate_timestep_count)
!                   CASE ( 'nr' )
!                      sums_wsnrs_ws_l(k,tn) = sums_wsnrs_ws_l(k,tn)         &
!                      + ( flux_t + diss_t )                                 &
!                      *   weight_substep(intermediate_timestep_count)
!
!                END SELECT

             ENDDO
         ENDDO
      ENDDO
      !$acc end kernels

    END SUBROUTINE advec_s_ws_acc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of u - Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_u_ws

       USE arrays_3d,                                                          &
           ONLY:  ddzw, drho_air, tend, u, v, w, rho_air, rho_air_zw

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count, u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nzb, nzb_max, nzt, wall_flags_0
           
       USE kinds
       
       USE statistics,                                                         &
           ONLY:  hom, sums_us2_ws_l, sums_wsus_ws_l, weight_substep

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit9  !<
       INTEGER(iwp) ::  ibit10 !<
       INTEGER(iwp) ::  ibit11 !<
       INTEGER(iwp) ::  ibit12 !<
       INTEGER(iwp) ::  ibit13 !<
       INTEGER(iwp) ::  ibit14 !<
       INTEGER(iwp) ::  ibit15 !<
       INTEGER(iwp) ::  ibit16 !<
       INTEGER(iwp) ::  ibit17 !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<
       
       REAL(wp)    ::  diss_d !<
       REAL(wp)    ::  div    !<
       REAL(wp)    ::  flux_d !<
       REAL(wp)    ::  gu     !<
       REAL(wp)    ::  gv     !<
       REAL(wp)    ::  v_comp !<
       REAL(wp)    ::  w_comp !<
       
       REAL(wp), DIMENSION(nzb+1:nzt) ::  swap_diss_y_local_u !<
       REAL(wp), DIMENSION(nzb+1:nzt) ::  swap_flux_y_local_u !<
       
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_diss_x_local_u !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_flux_x_local_u !<
       
       REAL(wp), DIMENSION(nzb:nzt) ::  diss_n !<
       REAL(wp), DIMENSION(nzb:nzt) ::  diss_r !<
       REAL(wp), DIMENSION(nzb:nzt) ::  diss_t !<
       REAL(wp), DIMENSION(nzb:nzt) ::  flux_n !<
       REAL(wp), DIMENSION(nzb:nzt) ::  flux_r !<
       REAL(wp), DIMENSION(nzb:nzt) ::  flux_t !<
       REAL(wp), DIMENSION(nzb:nzt) ::  u_comp !<
 
       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans

!
!--    Compute the fluxes for the whole left boundary of the processor domain.
       i = nxlu
       DO  j = nys, nyn
          DO  k = nzb+1, nzb_max

             ibit11 = IBITS(wall_flags_0(k,j,i-1),11,1)
             ibit10 = IBITS(wall_flags_0(k,j,i-1),10,1)
             ibit9  = IBITS(wall_flags_0(k,j,i-1),9,1)

             u_comp(k)                = u(k,j,i) + u(k,j,i-1) - gu
             swap_flux_x_local_u(k,j) = u_comp(k) * (                          &
                                       ( 37.0_wp * ibit11 * adv_mom_5             &
                                    +     7.0_wp * ibit10 * adv_mom_3             &
                                    +              ibit9  * adv_mom_1             &
                                       ) *                                     &
                                     ( u(k,j,i)   + u(k,j,i-1) )               &
                                -      (  8.0_wp * ibit11 * adv_mom_5             &
                                    +              ibit10 * adv_mom_3             &
                                       ) *                                     &
                                     ( u(k,j,i+1) + u(k,j,i-2) )               &
                                +      (           ibit11 * adv_mom_5             &
                                       ) *                                     &
                                     ( u(k,j,i+2) + u(k,j,i-3) )               &
                                                   )

              swap_diss_x_local_u(k,j) = - ABS( u_comp(k) ) * (                &
                                       ( 10.0_wp * ibit11 * adv_mom_5             &
                                    +     3.0_wp * ibit10 * adv_mom_3             &
                                    +              ibit9  * adv_mom_1             &
                                       ) *                                     &
                                     ( u(k,j,i)   - u(k,j,i-1) )               &
                                -      (  5.0_wp * ibit11 * adv_mom_5             &
                                    +              ibit10 * adv_mom_3             &
                                       ) *                                     &
                                     ( u(k,j,i+1) - u(k,j,i-2) )               &
                                +      (           ibit11 * adv_mom_5             &
                                       ) *                                     &
                                     ( u(k,j,i+2) - u(k,j,i-3) )               &
                                                             )

          ENDDO

          DO  k = nzb_max+1, nzt

             u_comp(k)         = u(k,j,i) + u(k,j,i-1) - gu
             swap_flux_x_local_u(k,j) = u_comp(k) * (                          &
                             37.0_wp * ( u(k,j,i) + u(k,j,i-1)   )                &
                           -  8.0_wp * ( u(k,j,i+1) + u(k,j,i-2) )                &
                           +           ( u(k,j,i+2) + u(k,j,i-3) ) ) * adv_mom_5
             swap_diss_x_local_u(k,j) = - ABS(u_comp(k)) * (                   &
                             10.0_wp * ( u(k,j,i) - u(k,j,i-1)   )                &
                           -  5.0_wp * ( u(k,j,i+1) - u(k,j,i-2) )                &
                           +           ( u(k,j,i+2) - u(k,j,i-3) ) ) * adv_mom_5

          ENDDO
       ENDDO

       DO i = nxlu, nxr
!       
!--       The following loop computes the fluxes for the south boundary points
          j = nys
          DO  k = nzb+1, nzb_max

             ibit14 = IBITS(wall_flags_0(k,j-1,i),14,1)
             ibit13 = IBITS(wall_flags_0(k,j-1,i),13,1)
             ibit12 = IBITS(wall_flags_0(k,j-1,i),12,1)

             v_comp                 = v(k,j,i) + v(k,j,i-1) - gv
             swap_flux_y_local_u(k) = v_comp * (                              &
                                   ( 37.0_wp * ibit14 * adv_mom_5                &
                                +     7.0_wp * ibit13 * adv_mom_3                &
                                +              ibit12 * adv_mom_1                &
                                   ) *                                        &
                                     ( u(k,j,i)   + u(k,j-1,i) )              &
                            -      (  8.0_wp * ibit14 * adv_mom_5                &
                            +                  ibit13 * adv_mom_3                    &
                                   ) *                                        &
                                     ( u(k,j+1,i) + u(k,j-2,i) )              &
                        +      (               ibit14 * adv_mom_5                    &
                               ) *                                            &
                                     ( u(k,j+2,i) + u(k,j-3,i) )              &
                                               )

             swap_diss_y_local_u(k) = - ABS ( v_comp ) * (                    &
                                   ( 10.0_wp * ibit14 * adv_mom_5                &
                                +     3.0_wp * ibit13 * adv_mom_3                &
                                +              ibit12 * adv_mom_1                &
                                   ) *                                        &
                                     ( u(k,j,i)   - u(k,j-1,i) )              &
                            -      (  5.0_wp * ibit14 * adv_mom_5                &
                                +              ibit13 * adv_mom_3                &
                                   ) *                                        &
                                     ( u(k,j+1,i) - u(k,j-2,i) )              &
                            +      (           ibit14 * adv_mom_5                &
                                   ) *                                        &
                                     ( u(k,j+2,i) - u(k,j-3,i) )              &
                                                         )

          ENDDO

          DO  k = nzb_max+1, nzt

             v_comp                 = v(k,j,i) + v(k,j,i-1) - gv
             swap_flux_y_local_u(k) = v_comp * (                              &
                           37.0_wp * ( u(k,j,i) + u(k,j-1,i)   )                 &
                         -  8.0_wp * ( u(k,j+1,i) + u(k,j-2,i) )                 &
                         +           ( u(k,j+2,i) + u(k,j-3,i) ) ) * adv_mom_5
             swap_diss_y_local_u(k) = - ABS(v_comp) * (                       &
                           10.0_wp * ( u(k,j,i) - u(k,j-1,i)   )                 &
                         -  5.0_wp * ( u(k,j+1,i) - u(k,j-2,i) )                 &
                         +           ( u(k,j+2,i) - u(k,j-3,i) ) ) * adv_mom_5

          ENDDO
!
!--       Computation of interior fluxes and tendency terms
          DO  j = nys, nyn

             flux_t(0) = 0.0_wp
             diss_t(0) = 0.0_wp
             flux_d    = 0.0_wp
             diss_d    = 0.0_wp

             DO  k = nzb+1, nzb_max

                ibit11 = IBITS(wall_flags_0(k,j,i),11,1)
                ibit10 = IBITS(wall_flags_0(k,j,i),10,1)
                ibit9  = IBITS(wall_flags_0(k,j,i),9,1)

                u_comp(k) = u(k,j,i+1) + u(k,j,i)
                flux_r(k) = ( u_comp(k) - gu ) * (                           &
                          ( 37.0_wp * ibit11 * adv_mom_5                        &
                       +     7.0_wp * ibit10 * adv_mom_3                        &
                       +              ibit9  * adv_mom_1                        &
                          ) *                                                &
                                 ( u(k,j,i+1) + u(k,j,i)   )                 &
                   -      (  8.0_wp * ibit11 * adv_mom_5                        &
                       +              ibit10 * adv_mom_3                        &
                          ) *                                                &
                                 ( u(k,j,i+2) + u(k,j,i-1) )                 &
                   +      (           ibit11 * adv_mom_5                        &
                          ) *                                                &
                                 ( u(k,j,i+3) + u(k,j,i-2) )                 &
                                                 )

                diss_r(k) = - ABS( u_comp(k) - gu ) * (                      &
                          ( 10.0_wp * ibit11 * adv_mom_5                        &
                       +     3.0_wp * ibit10 * adv_mom_3                        & 
                       +              ibit9  * adv_mom_1                        &
                          ) *                                                &
                                 ( u(k,j,i+1) - u(k,j,i)  )                  &
                   -      (  5.0_wp * ibit11 * adv_mom_5                        &
                       +              ibit10 * adv_mom_3                        &
                          ) *                                                &
                                 ( u(k,j,i+2) - u(k,j,i-1) )                 &
                   +      (           ibit11 * adv_mom_5                        &
                          ) *                                                &
                                 ( u(k,j,i+3) - u(k,j,i-2) )                 &
                                                     )

                ibit14 = IBITS(wall_flags_0(k,j,i),14,1)
                ibit13 = IBITS(wall_flags_0(k,j,i),13,1)
                ibit12 = IBITS(wall_flags_0(k,j,i),12,1)

                v_comp    = v(k,j+1,i) + v(k,j+1,i-1) - gv
                flux_n(k) = v_comp * (                                       &
                          ( 37.0_wp * ibit14 * adv_mom_5                        &
                       +     7.0_wp * ibit13 * adv_mom_3                        &
                       +              ibit12 * adv_mom_1                        &
                          ) *                                                &
                                 ( u(k,j+1,i) + u(k,j,i)   )                 &
                   -      (  8.0_wp * ibit14 * adv_mom_5                        &
                       +              ibit13 * adv_mom_3                        &
                          ) *                                                &
                                 ( u(k,j+2,i) + u(k,j-1,i) )                 &
                   +      (           ibit14 * adv_mom_5                        & 
                          ) *                                                &
                                 ( u(k,j+3,i) + u(k,j-2,i) )                 &
                                                 )

                diss_n(k) = - ABS ( v_comp ) * (                             &
                          ( 10.0_wp * ibit14 * adv_mom_5                        &
                       +     3.0_wp * ibit13 * adv_mom_3                        &
                       +              ibit12 * adv_mom_1                        &
                          ) *                                                &
                                 ( u(k,j+1,i) - u(k,j,i)  )                  &
                   -      (  5.0_wp * ibit14 * adv_mom_5                        &
                       +              ibit13 * adv_mom_3                        &
                          ) *                                                &
                                 ( u(k,j+2,i) - u(k,j-1,i) )                 &
                   +      (           ibit14 * adv_mom_5                        &
                          ) *                                                &
                                 ( u(k,j+3,i) - u(k,j-2,i) )                 &
                                                      )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit17 = IBITS(wall_flags_0(k,j,i),17,1)
                ibit16 = IBITS(wall_flags_0(k,j,i),16,1)
                ibit15 = IBITS(wall_flags_0(k,j,i),15,1)

                k_ppp = k + 3 * ibit17
                k_pp  = k + 2 * ( 1 - ibit15  )
                k_mm  = k - 2 * ibit17

                w_comp    = w(k,j,i) + w(k,j,i-1)
                flux_t(k) = w_comp * rho_air_zw(k) * (                       &
                          ( 37.0_wp * ibit17 * adv_mom_5                        &
                       +     7.0_wp * ibit16 * adv_mom_3                        &
                       +              ibit15 * adv_mom_1                        & 
                          ) *                                                &
                             ( u(k+1,j,i)  + u(k,j,i)     )                  &
                   -      (  8.0_wp * ibit17 * adv_mom_5                        &
                       +              ibit16 * adv_mom_3                        &
                          ) *                                                &
                             ( u(k_pp,j,i) + u(k-1,j,i)   )                  &
                   +      (           ibit17 * adv_mom_5                        &
                          ) *                                                &
                             ( u(k_ppp,j,i) + u(k_mm,j,i) )                  &
                                      )

                diss_t(k) = - ABS( w_comp ) * rho_air_zw(k) * (              &
                          ( 10.0_wp * ibit17 * adv_mom_5                        &
                       +     3.0_wp * ibit16 * adv_mom_3                        &
                       +              ibit15 * adv_mom_1                        &
                          ) *                                                &
                             ( u(k+1,j,i)   - u(k,j,i)    )                  &
                   -      (  5.0_wp * ibit17 * adv_mom_5                        &
                       +              ibit16 * adv_mom_3                        &
                          ) *                                                &
                             ( u(k_pp,j,i)  - u(k-1,j,i)  )                  &
                   +      (           ibit17 * adv_mom_5                        &
                           ) *                                               &
                             ( u(k_ppp,j,i) - u(k_mm,j,i) )                  &
                                              )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( u_comp(k) * ( ibit9 + ibit10 + ibit11 )             &
                - ( u(k,j,i)   + u(k,j,i-1)   )                               &
                                    * ( IBITS(wall_flags_0(k,j,i-1),9,1)      &
                                      + IBITS(wall_flags_0(k,j,i-1),10,1)     &
                                      + IBITS(wall_flags_0(k,j,i-1),11,1)     &
                                      )                                       &
                  ) * rho_air(k) * ddx                                        &
               +  ( ( v_comp + gv ) * ( ibit12 + ibit13 + ibit14 )            &
                  - ( v(k,j,i)   + v(k,j,i-1 )  )                             &
                                    * ( IBITS(wall_flags_0(k,j-1,i),12,1)     &
                                      + IBITS(wall_flags_0(k,j-1,i),13,1)     &
                                      + IBITS(wall_flags_0(k,j-1,i),14,1)     &
                                      )                                       &
                  ) * rho_air(k) * ddy                                        &
               +  ( w_comp * rho_air_zw(k) * ( ibit15 + ibit16 + ibit17 )     &
                - ( w(k-1,j,i) + w(k-1,j,i-1) ) * rho_air_zw(k-1)             &
                                    * ( IBITS(wall_flags_0(k-1,j,i),15,1)     &
                                      + IBITS(wall_flags_0(k-1,j,i),16,1)     &
                                      + IBITS(wall_flags_0(k-1,j,i),17,1)     &
                                      )                                       &  
                  ) * ddzw(k)   &
                ) * 0.5_wp



                tend(k,j,i) = tend(k,j,i) - (                                  &
                 ( flux_r(k) + diss_r(k)                                       &
               -   swap_flux_x_local_u(k,j) - swap_diss_x_local_u(k,j) ) * ddx &
               + ( flux_n(k) + diss_n(k)                                       &
               -   swap_flux_y_local_u(k)   - swap_diss_y_local_u(k)   ) * ddy &
               + ( ( flux_t(k) + diss_t(k) )                                   &
               -   ( flux_d    + diss_d    )                                   &
                                                    ) * drho_air(k) * ddzw(k)  &
                                           ) + div * u(k,j,i)

                swap_flux_x_local_u(k,j) = flux_r(k)
                swap_diss_x_local_u(k,j) = diss_r(k)
                swap_flux_y_local_u(k)   = flux_n(k)
                swap_diss_y_local_u(k)   = diss_n(k)
                flux_d                   = flux_t(k)
                diss_d                   = diss_t(k)
!
!--             Statistical Evaluation of u'u'. The factor has to be applied
!--             for right evaluation when gallilei_trans = .T. .
                sums_us2_ws_l(k,tn) = sums_us2_ws_l(k,tn)                      &
                + ( flux_r(k)                                                  &
                    * ( u_comp(k) - 2.0_wp * hom(k,1,1,0)                   )  &
                    / ( u_comp(k) - gu + SIGN( 1.0E-20_wp, u_comp(k) - gu ) )  &
                  + diss_r(k)                                                  &
                    *   ABS( u_comp(k) - 2.0_wp * hom(k,1,1,0)              )  &
                    / ( ABS( u_comp(k) - gu ) + 1.0E-20_wp                  )  &
                  ) *   weight_substep(intermediate_timestep_count)
!
!--             Statistical Evaluation of w'u'.
                sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn)                    &
                + ( flux_t(k)                                                  &
                    * ( w_comp - 2.0_wp * hom(k,1,3,0)                   )     &
                    / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              )     &
                  + diss_t(k)                                                  &
                    *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              )     &
                    / ( ABS( w_comp ) + 1.0E-20_wp                       )     &
                  ) *   weight_substep(intermediate_timestep_count)

             ENDDO

             DO  k = nzb_max+1, nzt

                u_comp(k) = u(k,j,i+1) + u(k,j,i)
                flux_r(k) = ( u_comp(k) - gu ) * (                            &
                         37.0_wp * ( u(k,j,i+1) + u(k,j,i)   )                   &
                       -  8.0_wp * ( u(k,j,i+2) + u(k,j,i-1) )                   &
                       +           ( u(k,j,i+3) + u(k,j,i-2) ) ) * adv_mom_5
                diss_r(k) = - ABS( u_comp(k) - gu ) * (                       &
                         10.0_wp * ( u(k,j,i+1) - u(k,j,i)   )                   &
                       -  5.0_wp * ( u(k,j,i+2) - u(k,j,i-1) )                   &
                       +           ( u(k,j,i+3) - u(k,j,i-2) ) ) * adv_mom_5

                v_comp    = v(k,j+1,i) + v(k,j+1,i-1) - gv
                flux_n(k) = v_comp * (                                        &
                         37.0_wp * ( u(k,j+1,i) + u(k,j,i)   )                   &
                       -  8.0_wp * ( u(k,j+2,i) + u(k,j-1,i) )                   &
                       +           ( u(k,j+3,i) + u(k,j-2,i) ) ) * adv_mom_5
                diss_n(k) = - ABS( v_comp ) * (                               &
                         10.0_wp * ( u(k,j+1,i) - u(k,j,i)   )                   &
                       -  5.0_wp * ( u(k,j+2,i) - u(k,j-1,i) )                   &
                       +           ( u(k,j+3,i) - u(k,j-2,i) ) ) * adv_mom_5
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit17 = IBITS(wall_flags_0(k,j,i),17,1)
                ibit16 = IBITS(wall_flags_0(k,j,i),16,1)
                ibit15 = IBITS(wall_flags_0(k,j,i),15,1)

                k_ppp = k + 3 * ibit17
                k_pp  = k + 2 * ( 1 - ibit15  )
                k_mm  = k - 2 * ibit17

                w_comp    = w(k,j,i) + w(k,j,i-1)
                flux_t(k) = w_comp * rho_air_zw(k) * (                       &
                          ( 37.0_wp * ibit17 * adv_mom_5                        &
                       +     7.0_wp * ibit16 * adv_mom_3                        &
                       +              ibit15 * adv_mom_1                        &
                          ) *                                                &
                             ( u(k+1,j,i)  + u(k,j,i)     )                  &
                   -      (  8.0_wp * ibit17 * adv_mom_5                        &
                       +              ibit16 * adv_mom_3                        &
                          ) *                                                &
                             ( u(k_pp,j,i) + u(k-1,j,i)   )                  &
                   +      (           ibit17 * adv_mom_5                        &
                          ) *                                                &
                             ( u(k_ppp,j,i) + u(k_mm,j,i) )                  &
                                      )

                diss_t(k) = - ABS( w_comp ) * rho_air_zw(k) * (              &
                          ( 10.0_wp * ibit17 * adv_mom_5                        &
                       +     3.0_wp * ibit16 * adv_mom_3                        &
                       +              ibit15 * adv_mom_1                        &
                          ) *                                                &
                             ( u(k+1,j,i)   - u(k,j,i)    )                  &
                   -      (  5.0_wp * ibit17 * adv_mom_5                        &
                       +              ibit16 * adv_mom_3                        &
                          ) *                                                &
                             ( u(k_pp,j,i)  - u(k-1,j,i)  )                  &
                   +      (           ibit17 * adv_mom_5                        &
                           ) *                                               &
                             ( u(k_ppp,j,i) - u(k_mm,j,i) )                  &
                                              )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( u_comp(k)   - ( u(k,j,i)   + u(k,j,i-1)   ) ) * ddx &
                                                                  * rho_air(k)&
                     +  ( v_comp + gv - ( v(k,j,i)   + v(k,j,i-1 )  ) ) * ddy &
                                                                  * rho_air(k)&
                     +  (   w_comp                      * rho_air_zw(k) -     &
                          ( w(k-1,j,i) + w(k-1,j,i-1) ) * rho_air_zw(k-1)     &
                        ) * ddzw(k)                                           &
                      ) * 0.5_wp

                tend(k,j,i) = tend(k,j,i) - (                                  &
                 ( flux_r(k) + diss_r(k)                                       &
               -   swap_flux_x_local_u(k,j) - swap_diss_x_local_u(k,j) ) * ddx &
               + ( flux_n(k) + diss_n(k)                                       &
               -   swap_flux_y_local_u(k)   - swap_diss_y_local_u(k)   ) * ddy &
               + ( ( flux_t(k) + diss_t(k) )                                   &
               -   ( flux_d    + diss_d    )                                   &
                                                    ) * drho_air(k) * ddzw(k)  &
                                           ) + div * u(k,j,i)

                swap_flux_x_local_u(k,j) = flux_r(k)
                swap_diss_x_local_u(k,j) = diss_r(k)
                swap_flux_y_local_u(k)   = flux_n(k)
                swap_diss_y_local_u(k)   = diss_n(k)
                flux_d                   = flux_t(k)
                diss_d                   = diss_t(k)
!
!--             Statistical Evaluation of u'u'. The factor has to be applied
!--             for right evaluation when gallilei_trans = .T. .
                sums_us2_ws_l(k,tn) = sums_us2_ws_l(k,tn)                      &
                + ( flux_r(k)                                                  &
                    * ( u_comp(k) - 2.0_wp * hom(k,1,1,0)                   )  &
                    / ( u_comp(k) - gu + SIGN( 1.0E-20_wp, u_comp(k) - gu ) )  &
                  + diss_r(k)                                                  &
                    *   ABS( u_comp(k) - 2.0_wp * hom(k,1,1,0)              )  &
                    / ( ABS( u_comp(k) - gu ) + 1.0E-20_wp                  )  &
                  ) *   weight_substep(intermediate_timestep_count)
!
!--             Statistical Evaluation of w'u'.
                sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn)                    &
                + ( flux_t(k)                                                  &
                    * ( w_comp - 2.0_wp * hom(k,1,3,0)                   )     &
                    / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              )     &
                  + diss_t(k)                                                  &
                    *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              )     &
                    / ( ABS( w_comp ) + 1.0E-20_wp                       )     &
                  ) *   weight_substep(intermediate_timestep_count)
             ENDDO
          ENDDO
       ENDDO
       sums_us2_ws_l(nzb,tn) = sums_us2_ws_l(nzb+1,tn)


    END SUBROUTINE advec_u_ws
    
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of u - Call for all grid points - accelerator version
!------------------------------------------------------------------------------!
    SUBROUTINE advec_u_ws_acc

       USE arrays_3d,                                                          &
           ONLY:  ddzw, drho_air, tend, u, v, w, rho_air, rho_air_zw

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count, u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  i_left, i_right, j_north, j_south, nxl, nxr, nyn, nys, nzb,  &
                  nzb_max, nzt, wall_flags_0
           
       USE kinds
       
!        USE statistics,                                                       &
!            ONLY:  hom, sums_us2_ws_l, sums_wsus_ws_l, weight_substep

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit9  !<
       INTEGER(iwp) ::  ibit10 !<
       INTEGER(iwp) ::  ibit11 !<
       INTEGER(iwp) ::  ibit12 !<
       INTEGER(iwp) ::  ibit13 !<
       INTEGER(iwp) ::  ibit14 !<
       INTEGER(iwp) ::  ibit15 !<
       INTEGER(iwp) ::  ibit16 !<
       INTEGER(iwp) ::  ibit17 !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mmm  !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<

       REAL(wp)    ::  diss_d   !<
       REAL(wp)    ::  diss_l   !<
       REAL(wp)    ::  diss_n   !<
       REAL(wp)    ::  diss_r   !<
       REAL(wp)    ::  diss_s   !<
       REAL(wp)    ::  diss_t   !<
       REAL(wp)    ::  div      !<
       REAL(wp)    ::  flux_d   !<
       REAL(wp)    ::  flux_l   !<
       REAL(wp)    ::  flux_n   !<
       REAL(wp)    ::  flux_r   !<
       REAL(wp)    ::  flux_s   !<
       REAL(wp)    ::  flux_t   !<
       REAL(wp)    ::  gu       !<
       REAL(wp)    ::  gv       !<
       REAL(wp)    ::  u_comp   !<
       REAL(wp)    ::  u_comp_l !<
       REAL(wp)    ::  v_comp   !<
       REAL(wp)    ::  v_comp_s !<
       REAL(wp)    ::  w_comp   !<


       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans

!
!--    Computation of fluxes and tendency terms
       !$acc  kernels present( ddzw, tend, u, v, w, wall_flags_0 )
       DO i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = nzb+1, nzt

                ibit11 = IBITS(wall_flags_0(k,j,i-1),11,1)
                ibit10 = IBITS(wall_flags_0(k,j,i-1),10,1)
                ibit9  = IBITS(wall_flags_0(k,j,i-1),9,1)

                u_comp_l           = u(k,j,i) + u(k,j,i-1) - gu
                flux_l             = u_comp_l * (                          &
                                    ( 37.0_wp * ibit11 * adv_mom_5             &
                                 +     7.0_wp * ibit10 * adv_mom_3             &
                                 +              ibit9  * adv_mom_1             &
                                    ) *                                     &
                                  ( u(k,j,i)   + u(k,j,i-1) )               &
                             -      (  8.0_wp * ibit11 * adv_mom_5             &
                                 +              ibit10 * adv_mom_3             &
                                    ) *                                     &
                                  ( u(k,j,i+1) + u(k,j,i-2) )               &
                             +      (           ibit11 * adv_mom_5             &
                                    ) *                                     &
                                  ( u(k,j,i+2) + u(k,j,i-3) )               &
                                                )

                diss_l             = - ABS( u_comp_l ) * (                &
                                   ( 10.0_wp * ibit11 * adv_mom_5             &
                                +     3.0_wp * ibit10 * adv_mom_3             &
                                +              ibit9  * adv_mom_1             &
                                   ) *                                     &
                                 ( u(k,j,i)   - u(k,j,i-1) )               &
                            -      (  5.0_wp * ibit11 * adv_mom_5             &
                                +              ibit10 * adv_mom_3             &
                                   ) *                                     &
                                 ( u(k,j,i+1) - u(k,j,i-2) )               &
                            +      (           ibit11 * adv_mom_5             &
                                   ) *                                     &
                                 ( u(k,j,i+2) - u(k,j,i-3) )               &
                                                         )

                ibit11 = IBITS(wall_flags_0(k,j,i),11,1)
                ibit10 = IBITS(wall_flags_0(k,j,i),10,1)
                ibit9  = IBITS(wall_flags_0(k,j,i),9,1)

                u_comp    = u(k,j,i+1) + u(k,j,i)
                flux_r    = ( u_comp   - gu ) * (                           &
                          ( 37.0_wp * ibit11 * adv_mom_5                        &
                       +     7.0_wp * ibit10 * adv_mom_3                        &
                       +              ibit9  * adv_mom_1                        &
                          ) *                                                &
                                 ( u(k,j,i+1) + u(k,j,i)   )                 &
                   -      (  8.0_wp * ibit11 * adv_mom_5                        &
                       +              ibit10 * adv_mom_3                        &
                          ) *                                                &
                                 ( u(k,j,i+2) + u(k,j,i-1) )                 &
                   +      (           ibit11 * adv_mom_5                        &
                          ) *                                                &
                                 ( u(k,j,i+3) + u(k,j,i-2) )                 &
                                                 )

                diss_r    = - ABS( u_comp    - gu ) * (                      &
                          ( 10.0_wp * ibit11 * adv_mom_5                        &
                       +     3.0_wp * ibit10 * adv_mom_3                        &
                       +              ibit9  * adv_mom_1                        &
                          ) *                                                &
                                 ( u(k,j,i+1) - u(k,j,i)  )                  &
                   -      (  5.0_wp * ibit11 * adv_mom_5                        &
                       +              ibit10 * adv_mom_3                        &
                          ) *                                                &
                                 ( u(k,j,i+2) - u(k,j,i-1) )                 &
                   +      (           ibit11 * adv_mom_5                        &
                          ) *                                                &
                                 ( u(k,j,i+3) - u(k,j,i-2) )                 &
                                                     )

                ibit14 = IBITS(wall_flags_0(k,j-1,i),14,1)
                ibit13 = IBITS(wall_flags_0(k,j-1,i),13,1)
                ibit12 = IBITS(wall_flags_0(k,j-1,i),12,1)

                v_comp_s                 = v(k,j,i) + v(k,j,i-1) - gv
                flux_s                   = v_comp_s * (                       &
                                   ( 37.0_wp * ibit14 * adv_mom_5                &
                                +     7.0_wp * ibit13 * adv_mom_3                &
                                +              ibit12 * adv_mom_1                &
                                   ) *                                         &
                                     ( u(k,j,i)   + u(k,j-1,i) )              &
                            -      (  8.0_wp * ibit14 * adv_mom_5                &
                            +                  ibit13 * adv_mom_3                    &
                                   ) *                                        &
                                     ( u(k,j+1,i) + u(k,j-2,i) )              &
                        +      (               ibit14 * adv_mom_5                    &
                               ) *                                            &
                                     ( u(k,j+2,i) + u(k,j-3,i) )              &
                                               )

                diss_s                  = - ABS ( v_comp_s ) * (              &
                                   ( 10.0_wp * ibit14 * adv_mom_5                &
                                +     3.0_wp * ibit13 * adv_mom_3                &
                                +              ibit12 * adv_mom_1                &
                                   ) *                                        &
                                     ( u(k,j,i)   - u(k,j-1,i) )              &
                            -      (  5.0_wp * ibit14 * adv_mom_5                &
                                +              ibit13 * adv_mom_3                &
                                   ) *                                        &
                                     ( u(k,j+1,i) - u(k,j-2,i) )              &
                            +      (           ibit14 * adv_mom_5                &
                                   ) *                                        &
                                     ( u(k,j+2,i) - u(k,j-3,i) )              &
                                                         )


                ibit14 = IBITS(wall_flags_0(k,j,i),14,1)
                ibit13 = IBITS(wall_flags_0(k,j,i),13,1)
                ibit12 = IBITS(wall_flags_0(k,j,i),12,1)

                v_comp    = v(k,j+1,i) + v(k,j+1,i-1) - gv
                flux_n    = v_comp * (                                       &
                          ( 37.0_wp * ibit14 * adv_mom_5                        &
                       +     7.0_wp * ibit13 * adv_mom_3                        &
                       +              ibit12 * adv_mom_1                        &
                          ) *                                                &
                                 ( u(k,j+1,i) + u(k,j,i)   )                 &
                   -      (  8.0_wp * ibit14 * adv_mom_5                        &
                       +              ibit13 * adv_mom_3                        &
                          ) *                                                &
                                 ( u(k,j+2,i) + u(k,j-1,i) )                 &
                   +      (           ibit14 * adv_mom_5                        &
                          ) *                                                &
                                 ( u(k,j+3,i) + u(k,j-2,i) )                 &
                                                 )

                diss_n    = - ABS ( v_comp ) * (                             &
                          ( 10.0_wp * ibit14 * adv_mom_5                        &
                       +     3.0_wp * ibit13 * adv_mom_3                        &
                       +              ibit12 * adv_mom_1                        &
                          ) *                                                &
                                 ( u(k,j+1,i) - u(k,j,i)  )                  &
                   -      (  5.0_wp * ibit14 * adv_mom_5                        &
                       +              ibit13 * adv_mom_3                        &
                          ) *                                                &
                                 ( u(k,j+2,i) - u(k,j-1,i) )                 &
                   +      (           ibit14 * adv_mom_5                        &
                          ) *                                                &
                                 ( u(k,j+3,i) - u(k,j-2,i) )                 &
                                                      )

                ibit17 = IBITS(wall_flags_0(k-1,j,i),17,1)
                ibit16 = IBITS(wall_flags_0(k-1,j,i),16,1)
                ibit15 = IBITS(wall_flags_0(k-1,j,i),15,1)

                k_pp  = k + 2 * ibit17
                k_mm  = k - 2 * ( ibit16 + ibit17 )
                k_mmm = k - 3 * ibit17

                w_comp    = w(k-1,j,i) + w(k-1,j,i-1)
                flux_d    = w_comp  * (                                      &
                          ( 37.0_wp * ibit17 * adv_mom_5                        &
                       +     7.0_wp * ibit16 * adv_mom_3                        &
                       +              ibit15 * adv_mom_1                        &
                          ) *                                                &
                             ( u(k,j,i)    + u(k-1,j,i)   )                  &
                   -      (  8.0_wp * ibit17 * adv_mom_5                        &
                       +              ibit16 * adv_mom_3                        &
                          ) *                                                &
                             ( u(k+1,j,i) + u(k_mm,j,i)   )                  &
                   +      (           ibit17 * adv_mom_5                        &
                          ) *                                                 &
                             ( u(k_pp,j,i) + u(k_mmm,j,i) )                  &
                                      )

                diss_d    = - ABS( w_comp ) * (                              &
                          ( 10.0_wp * ibit17 * adv_mom_5                        &
                       +     3.0_wp * ibit16 * adv_mom_3                        &
                       +              ibit15 * adv_mom_1                        &
                          ) *                                                &
                             ( u(k,j,i)     - u(k-1,j,i)  )                  &
                   -      (  5.0_wp * ibit17 * adv_mom_5                        &
                       +              ibit16 * adv_mom_3                        &
                          ) *                                                &
                             ( u(k+1,j,i)  - u(k_mm,j,i)  )                  &
                   +      (           ibit17 * adv_mom_5                        &
                           ) *                                               &
                             ( u(k_pp,j,i) - u(k_mmm,j,i) )                  &
                                              )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit17 = IBITS(wall_flags_0(k,j,i),17,1)
                ibit16 = IBITS(wall_flags_0(k,j,i),16,1)
                ibit15 = IBITS(wall_flags_0(k,j,i),15,1)

                k_ppp = k + 3 * ibit17
                k_pp  = k + 2 * ( 1 - ibit15  )
                k_mm  = k - 2 * ibit17

                w_comp    = w(k,j,i) + w(k,j,i-1)
                flux_t    = w_comp * rho_air_zw(k) * (                       &
                          ( 37.0_wp * ibit17 * adv_mom_5                        &
                       +     7.0_wp * ibit16 * adv_mom_3                        &
                       +              ibit15 * adv_mom_1                        &
                          ) *                                                &
                             ( u(k+1,j,i)  + u(k,j,i)     )                  &
                   -      (  8.0_wp * ibit17 * adv_mom_5                        &
                       +              ibit16 * adv_mom_3                        &
                          ) *                                                &
                             ( u(k_pp,j,i) + u(k-1,j,i)   )                  &
                   +      (           ibit17 * adv_mom_5                        &
                          ) *                                                &
                             ( u(k_ppp,j,i) + u(k_mm,j,i) )                  &
                                      )

                diss_t    = - ABS( w_comp ) * rho_air_zw(k) * (              &
                          ( 10.0_wp * ibit17 * adv_mom_5                        &
                       +     3.0_wp * ibit16 * adv_mom_3                        &
                       +              ibit15 * adv_mom_1                        &
                          ) *                                                &
                             ( u(k+1,j,i)   - u(k,j,i)    )                  &
                   -      (  5.0_wp * ibit17 * adv_mom_5                        &
                       +              ibit16 * adv_mom_3                        &
                          ) *                                                &
                             ( u(k_pp,j,i)  - u(k-1,j,i)  )                  &
                   +      (           ibit17 * adv_mom_5                        &
                           ) *                                               &
                             ( u(k_ppp,j,i) - u(k_mm,j,i) )                  &
                                              )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( u_comp    * ( ibit9 + ibit10 + ibit11 )             &
                - ( u(k,j,i)   + u(k,j,i-1)   )                               &
                                    * ( IBITS(wall_flags_0(k,j,i-1),9,1)      &
                                      + IBITS(wall_flags_0(k,j,i-1),10,1)     &
                                      + IBITS(wall_flags_0(k,j,i-1),11,1)     &
                                      )                                       &
                  ) * rho_air(k) * ddx                                        &
               +  ( ( v_comp + gv ) * ( ibit12 + ibit13 + ibit14 )            &
                  - ( v(k,j,i)   + v(k,j,i-1 )  )                             &
                                    * ( IBITS(wall_flags_0(k,j-1,i),12,1)     &
                                      + IBITS(wall_flags_0(k,j-1,i),13,1)     &
                                      + IBITS(wall_flags_0(k,j-1,i),14,1)     &
                                      )                                       &
                  ) * rho_air(k) * ddy                                        &
               +  ( w_comp * rho_air_zw(k) * ( ibit15 + ibit16 + ibit17 )     &
                - ( w(k-1,j,i) + w(k-1,j,i-1) ) * rho_air_zw(k-1)             &
                                    * ( IBITS(wall_flags_0(k-1,j,i),15,1)     &
                                      + IBITS(wall_flags_0(k-1,j,i),16,1)     &
                                      + IBITS(wall_flags_0(k-1,j,i),17,1)     &
                                      )                                       &
                  ) * ddzw(k)   &
                ) * 0.5_wp


                tend(k,j,i) = - (                                              &
                               ( flux_r + diss_r - flux_l - diss_l ) * ddx     &
                             + ( flux_n + diss_n - flux_s - diss_s ) * ddy     &
                             + ( ( flux_t + diss_t ) -                         &
                                 ( flux_d + diss_d )                           &
                                                     ) * drho_air(k) * ddzw(k) &
                                ) + div * u(k,j,i)

!++
!--             Statistical Evaluation of u'u'. The factor has to be applied
!--             for right evaluation when gallilei_trans = .T. .
!                sums_us2_ws_l(k,tn) = sums_us2_ws_l(k,tn)                     &
!                              + ( flux_r    *                                 &
!                                ( u_comp    - 2.0_wp * hom(k,1,1,0) )            &
!                              / ( u_comp    - gu + 1.0E-20_wp   )             &
!                              +   diss_r    *                                 &
!                                  ABS( u_comp    - 2.0_wp * hom(k,1,1,0) )       &
!                              / ( ABS( u_comp    - gu ) + 1.0E-20_wp ) )      &
!                              *   weight_substep(intermediate_timestep_count)
!
!--             Statistical Evaluation of w'u'.
!                sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn)                   &
!                              + ( flux_t    + diss_t    )                     &
!                              *   weight_substep(intermediate_timestep_count)
             ENDDO
          ENDDO
       ENDDO
       !$acc end kernels

!++
!       sums_us2_ws_l(nzb,tn) = sums_us2_ws_l(nzb+1,tn)

    END SUBROUTINE advec_u_ws_acc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of v - Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_v_ws

       USE arrays_3d,                                                          &
           ONLY:  ddzw, drho_air, tend, u, v, w, rho_air, rho_air_zw

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count, u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nysv, nzb, nzb_max, nzt, wall_flags_0

       USE kinds

       USE statistics,                                                         &
           ONLY:  hom, sums_vs2_ws_l, sums_wsvs_ws_l, weight_substep

       IMPLICIT NONE


       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit18 !<
       INTEGER(iwp) ::  ibit19 !<
       INTEGER(iwp) ::  ibit20 !<
       INTEGER(iwp) ::  ibit21 !<
       INTEGER(iwp) ::  ibit22 !<
       INTEGER(iwp) ::  ibit23 !<
       INTEGER(iwp) ::  ibit24 !<
       INTEGER(iwp) ::  ibit25 !<
       INTEGER(iwp) ::  ibit26 !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<
       
       REAL(wp)    ::  diss_d !<
       REAL(wp)    ::  div    !<
       REAL(wp)    ::  flux_d !<
       REAL(wp)    ::  gu     !<
       REAL(wp)    ::  gv     !<
       REAL(wp)    ::  u_comp !<
       REAL(wp)    ::  w_comp !<
       
       REAL(wp), DIMENSION(nzb+1:nzt) ::  swap_diss_y_local_v !<
       REAL(wp), DIMENSION(nzb+1:nzt) ::  swap_flux_y_local_v !<
       
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_diss_x_local_v !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_flux_x_local_v !<
       
       REAL(wp), DIMENSION(nzb:nzt) ::  diss_n !<
       REAL(wp), DIMENSION(nzb:nzt) ::  diss_r !<
       REAL(wp), DIMENSION(nzb:nzt) ::  diss_t !<
       REAL(wp), DIMENSION(nzb:nzt) ::  flux_n !<
       REAL(wp), DIMENSION(nzb:nzt) ::  flux_r !<
       REAL(wp), DIMENSION(nzb:nzt) ::  flux_t !<
       REAL(wp), DIMENSION(nzb:nzt) ::  v_comp !<

       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans
!
!--    First compute the whole left boundary of the processor domain
       i = nxl
       DO  j = nysv, nyn
          DO  k = nzb+1, nzb_max

             ibit20 = IBITS(wall_flags_0(k,j,i-1),20,1)
             ibit19 = IBITS(wall_flags_0(k,j,i-1),19,1)
             ibit18 = IBITS(wall_flags_0(k,j,i-1),18,1)

             u_comp                   = u(k,j-1,i) + u(k,j,i) - gu
             swap_flux_x_local_v(k,j) = u_comp * (                             &
                                      ( 37.0_wp * ibit20 * adv_mom_5              &
                                   +     7.0_wp * ibit19 * adv_mom_3              &
                                   +              ibit18 * adv_mom_1              &
                                      ) *                                      &
                                     ( v(k,j,i)   + v(k,j,i-1) )               &
                               -      (  8.0_wp * ibit20 * adv_mom_5              &
                                   +              ibit19 * adv_mom_3              &
                                      ) *                                      &
                                     ( v(k,j,i+1) + v(k,j,i-2) )               &
                               +      (           ibit20 * adv_mom_5              &
                                      ) *                                      &
                                     ( v(k,j,i+2) + v(k,j,i-3) )               &
                                                 )

              swap_diss_x_local_v(k,j) = - ABS( u_comp ) * (                   &
                                      ( 10.0_wp * ibit20 * adv_mom_5              &
                                   +     3.0_wp * ibit19 * adv_mom_3              &
                                   +              ibit18 * adv_mom_1              &
                                      ) *                                      &
                                     ( v(k,j,i)   - v(k,j,i-1) )               &
                               -      (  5.0_wp * ibit20 * adv_mom_5              &
                                   +              ibit19 * adv_mom_3              &
                                      ) *                                      &
                                     ( v(k,j,i+1) - v(k,j,i-2) )               &
                               +      (           ibit20 * adv_mom_5              &
                                      ) *                                      &
                                     ( v(k,j,i+2) - v(k,j,i-3) )               &
                                                           )

          ENDDO

          DO  k = nzb_max+1, nzt

             u_comp                   = u(k,j-1,i) + u(k,j,i) - gu
             swap_flux_x_local_v(k,j) = u_comp * (                            &
                             37.0_wp * ( v(k,j,i) + v(k,j,i-1)   )               &
                           -  8.0_wp * ( v(k,j,i+1) + v(k,j,i-2) )               &
                           +           ( v(k,j,i+2) + v(k,j,i-3) ) ) * adv_mom_5
             swap_diss_x_local_v(k,j) = - ABS( u_comp ) * (                   &
                             10.0_wp * ( v(k,j,i) - v(k,j,i-1)   )               &
                           -  5.0_wp * ( v(k,j,i+1) - v(k,j,i-2) )               &
                           +           ( v(k,j,i+2) - v(k,j,i-3) ) ) * adv_mom_5

          ENDDO

       ENDDO

       DO i = nxl, nxr

          j = nysv
          DO  k = nzb+1, nzb_max

             ibit23 = IBITS(wall_flags_0(k,j-1,i),23,1)
             ibit22 = IBITS(wall_flags_0(k,j-1,i),22,1)
             ibit21 = IBITS(wall_flags_0(k,j-1,i),21,1)

             v_comp(k)              = v(k,j,i) + v(k,j-1,i) - gv
             swap_flux_y_local_v(k) = v_comp(k) * (                           &
                                   ( 37.0_wp * ibit23 * adv_mom_5                &
                                +     7.0_wp * ibit22 * adv_mom_3                &
                                +              ibit21 * adv_mom_1                &
                                   ) *                                        &
                                     ( v(k,j,i)   + v(k,j-1,i) )              &
                            -      (  8.0_wp * ibit23 * adv_mom_5                &
                                +              ibit22 * adv_mom_3                &
                                   ) *                                        &
                                     ( v(k,j+1,i) + v(k,j-2,i) )              &
                            +      (           ibit23 * adv_mom_5                &
                                   ) *                                        &
                                     ( v(k,j+2,i) + v(k,j-3,i) )              &
                                                 )

             swap_diss_y_local_v(k) = - ABS( v_comp(k) ) * (                  &
                                   ( 10.0_wp * ibit23 * adv_mom_5                &
                                +     3.0_wp * ibit22 * adv_mom_3                &
                                +              ibit21 * adv_mom_1                &
                                   ) *                                        &
                                     ( v(k,j,i)   - v(k,j-1,i) )              &
                            -      (  5.0_wp * ibit23 * adv_mom_5                &
                                +              ibit22 * adv_mom_3                &
                                   ) *                                        &
                                     ( v(k,j+1,i) - v(k,j-2,i) )              &
                            +      (           ibit23 * adv_mom_5                &
                                   ) *                                        &
                                     ( v(k,j+2,i) - v(k,j-3,i) )              &
                                                          )

          ENDDO

          DO  k = nzb_max+1, nzt

             v_comp(k)              = v(k,j,i) + v(k,j-1,i) - gv
             swap_flux_y_local_v(k) = v_comp(k) * (                           &
                           37.0_wp * ( v(k,j,i) + v(k,j-1,i)   )                 &
                         -  8.0_wp * ( v(k,j+1,i) + v(k,j-2,i) )                 &
                         +           ( v(k,j+2,i) + v(k,j-3,i) ) ) * adv_mom_5
             swap_diss_y_local_v(k) = - ABS( v_comp(k) ) * (                  &
                           10.0_wp * ( v(k,j,i) - v(k,j-1,i)   )                 &
                         -  5.0_wp * ( v(k,j+1,i) - v(k,j-2,i) )                 &
                         +           ( v(k,j+2,i) - v(k,j-3,i) ) ) * adv_mom_5

          ENDDO

          DO  j = nysv, nyn

             flux_t(0) = 0.0_wp
             diss_t(0) = 0.0_wp
             flux_d    = 0.0_wp
             diss_d    = 0.0_wp

             DO  k = nzb+1, nzb_max

                ibit20 = IBITS(wall_flags_0(k,j,i),20,1)
                ibit19 = IBITS(wall_flags_0(k,j,i),19,1)
                ibit18 = IBITS(wall_flags_0(k,j,i),18,1)

                u_comp    = u(k,j-1,i+1) + u(k,j,i+1) - gu
                flux_r(k) = u_comp * (                                       &
                          ( 37.0_wp * ibit20 * adv_mom_5                        &
                       +     7.0_wp * ibit19 * adv_mom_3                        &
                       +              ibit18 * adv_mom_1                        &
                          ) *                                                &
                                 ( v(k,j,i+1) + v(k,j,i)   )                 &
                   -      (  8.0_wp * ibit20 * adv_mom_5                        &
                       +              ibit19 * adv_mom_3                        &
                          ) *                                                &
                                 ( v(k,j,i+2) + v(k,j,i-1) )                 &
                   +      (           ibit20 * adv_mom_5                        &
                          ) *                                                &
                                 ( v(k,j,i+3) + v(k,j,i-2) )                 &
                                     )

                diss_r(k) = - ABS( u_comp ) * (                              &
                          ( 10.0_wp * ibit20 * adv_mom_5                        &
                       +     3.0_wp * ibit19 * adv_mom_3                        &
                       +              ibit18 * adv_mom_1                        &
                          ) *                                                &
                                 ( v(k,j,i+1) - v(k,j,i)  )                  &
                   -      (  5.0_wp * ibit20 * adv_mom_5                        &
                       +              ibit19 * adv_mom_3                        &
                          ) *                                                &
                                 ( v(k,j,i+2) - v(k,j,i-1) )                 &
                   +      (           ibit20 * adv_mom_5                        &
                          ) *                                                &
                                 ( v(k,j,i+3) - v(k,j,i-2) )                 &
                                              )

                ibit23 = IBITS(wall_flags_0(k,j,i),23,1)
                ibit22 = IBITS(wall_flags_0(k,j,i),22,1)
                ibit21 = IBITS(wall_flags_0(k,j,i),21,1)

                v_comp(k) = v(k,j+1,i) + v(k,j,i)
                flux_n(k) = ( v_comp(k) - gv ) * (                           &
                          ( 37.0_wp * ibit23 * adv_mom_5                        &
                       +     7.0_wp * ibit22 * adv_mom_3                        &
                       +              ibit21 * adv_mom_1                        &
                          ) *                                                &
                                 ( v(k,j+1,i) + v(k,j,i)   )                 &
                   -      (  8.0_wp * ibit23 * adv_mom_5                        &
                       +              ibit22 * adv_mom_3                        &
                          ) *                                                &
                                 ( v(k,j+2,i) + v(k,j-1,i) )                 &
                   +      (           ibit23 * adv_mom_5                        &
                          ) *                                                &
                                 ( v(k,j+3,i) + v(k,j-2,i) )                 &
                                     )

                diss_n(k) = - ABS( v_comp(k) - gv ) * (                      &
                          ( 10.0_wp * ibit23 * adv_mom_5                        &
                       +     3.0_wp * ibit22 * adv_mom_3                        &
                       +              ibit21 * adv_mom_1                        &
                          ) *                                                &
                                 ( v(k,j+1,i) - v(k,j,i)  )                  &
                   -      (  5.0_wp * ibit23 * adv_mom_5                        &
                       +              ibit22 * adv_mom_3                        &
                          ) *                                                &
                                 ( v(k,j+2,i) - v(k,j-1,i) )                 &
                   +      (           ibit23 * adv_mom_5                        &
                          ) *                                                &
                                 ( v(k,j+3,i) - v(k,j-2,i) )                 &
                                                      )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit26 = IBITS(wall_flags_0(k,j,i),26,1)
                ibit25 = IBITS(wall_flags_0(k,j,i),25,1)
                ibit24 = IBITS(wall_flags_0(k,j,i),24,1)

                k_ppp = k + 3 * ibit26
                k_pp  = k + 2 * ( 1 - ibit24  )
                k_mm  = k - 2 * ibit26

                w_comp    = w(k,j-1,i) + w(k,j,i)
                flux_t(k) = w_comp * rho_air_zw(k) * (                       &
                          ( 37.0_wp * ibit26 * adv_mom_5                        &
                       +     7.0_wp * ibit25 * adv_mom_3                        &
                       +              ibit24 * adv_mom_1                        &
                          ) *                                                &
                             ( v(k+1,j,i)   + v(k,j,i)    )                  &
                   -      (  8.0_wp * ibit26 * adv_mom_5                        &
                       +              ibit25 * adv_mom_3                        &
                          ) *                                                &
                             ( v(k_pp,j,i)  + v(k-1,j,i)  )                  &
                   +      (           ibit26 * adv_mom_5                        &
                          ) *                                                &
                             ( v(k_ppp,j,i) + v(k_mm,j,i) )                  &
                                      )

                diss_t(k) = - ABS( w_comp ) * rho_air_zw(k) * (              &
                          ( 10.0_wp * ibit26 * adv_mom_5                        &
                       +     3.0_wp * ibit25 * adv_mom_3                        &
                       +              ibit24 * adv_mom_1                        &
                          ) *                                                &
                             ( v(k+1,j,i)   - v(k,j,i)    )                  &
                   -      (  5.0_wp * ibit26 * adv_mom_5                        &
                       +              ibit25 * adv_mom_3                        &
                          ) *                                                &
                             ( v(k_pp,j,i)  - v(k-1,j,i)  )                  &
                   +      (           ibit26 * adv_mom_5                        &
                          ) *                                                &
                             ( v(k_ppp,j,i) - v(k_mm,j,i) )                  &
                                               )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( ( u_comp     + gu )                                 &
                                       * ( ibit18 + ibit19 + ibit20 )         &
                - ( u(k,j-1,i)   + u(k,j,i) )                                 &
                                       * ( IBITS(wall_flags_0(k,j,i-1),18,1)  &
                                         + IBITS(wall_flags_0(k,j,i-1),19,1)  &
                                         + IBITS(wall_flags_0(k,j,i-1),20,1)  &
                                         )                                    &
                  ) * rho_air(k) * ddx                                        &
               +  ( v_comp(k)                                                 &
                                       * ( ibit21 + ibit22 + ibit23 )         &
                - ( v(k,j,i)     + v(k,j-1,i) )                               &
                                       * ( IBITS(wall_flags_0(k,j-1,i),21,1)  &
                                         + IBITS(wall_flags_0(k,j-1,i),22,1)  &
                                         + IBITS(wall_flags_0(k,j-1,i),23,1)  &
                                         )                                    &
                  ) * rho_air(k) * ddy                                        &
               +  ( w_comp * rho_air_zw(k)                                    &
                                       * ( ibit24 + ibit25 + ibit26 )         &
                - ( w(k-1,j-1,i) + w(k-1,j,i) ) * rho_air_zw(k-1)             &
                                       * ( IBITS(wall_flags_0(k-1,j,i),24,1)  &
                                         + IBITS(wall_flags_0(k-1,j,i),25,1)  &
                                         + IBITS(wall_flags_0(k-1,j,i),26,1)  &
                                         )                                    &
                   ) * ddzw(k)   &
                ) * 0.5_wp


                tend(k,j,i) = tend(k,j,i) - (                                 &
                       ( flux_r(k) + diss_r(k)                                &
                     -   swap_flux_x_local_v(k,j) - swap_diss_x_local_v(k,j)  &
                       ) * ddx                                                &
                     + ( flux_n(k) + diss_n(k)                                &
                     -   swap_flux_y_local_v(k) - swap_diss_y_local_v(k)      &
                       ) * ddy                                                &
                     + ( ( flux_t(k) + diss_t(k) )                            &
                     -   ( flux_d    + diss_d    )                            &
                       ) * drho_air(k) * ddzw(k)                              &
                                            )  + v(k,j,i) * div

                swap_flux_x_local_v(k,j) = flux_r(k)
                swap_diss_x_local_v(k,j) = diss_r(k)
                swap_flux_y_local_v(k)   = flux_n(k)
                swap_diss_y_local_v(k)   = diss_n(k)
                flux_d                   = flux_t(k)
                diss_d                   = diss_t(k)

!
!--             Statistical Evaluation of v'v'. The factor has to be applied
!--             for right evaluation when gallilei_trans = .T. .
                sums_vs2_ws_l(k,tn) = sums_vs2_ws_l(k,tn)                      &
                + ( flux_n(k)                                                  &
                    * ( v_comp(k) - 2.0_wp * hom(k,1,2,0)                   )  &
                    / ( v_comp(k) - gv + SIGN( 1.0E-20_wp, v_comp(k) - gv ) )  &
               +   diss_n(k)                                                   &
                    *   ABS( v_comp(k) - 2.0_wp * hom(k,1,2,0)              )  &
                    / ( ABS( v_comp(k) - gv ) + 1.0E-20_wp                  )  &
                  ) *   weight_substep(intermediate_timestep_count)
!
!--             Statistical Evaluation of w'u'.
                sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn)                    &
                + ( flux_t(k)                                                  &
                    * ( w_comp - 2.0_wp * hom(k,1,3,0)                   )     &
                    / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              )     &
               +   diss_t(k)                                                   &
                    *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              )     &
                    / ( ABS( w_comp ) + 1.0E-20_wp                       )     &
                  ) *   weight_substep(intermediate_timestep_count)

             ENDDO

             DO  k = nzb_max+1, nzt

                u_comp    = u(k,j-1,i+1) + u(k,j,i+1) - gu
                flux_r(k) = u_comp * (                                        &
                      37.0_wp * ( v(k,j,i+1) + v(k,j,i)   )                      &
                    -  8.0_wp * ( v(k,j,i+2) + v(k,j,i-1) )                      &
                    +           ( v(k,j,i+3) + v(k,j,i-2) ) ) * adv_mom_5

                diss_r(k) = - ABS( u_comp ) * (                               &
                      10.0_wp * ( v(k,j,i+1) - v(k,j,i) )                        &
                    -  5.0_wp * ( v(k,j,i+2) - v(k,j,i-1) )                      &
                    +           ( v(k,j,i+3) - v(k,j,i-2) ) ) * adv_mom_5


                v_comp(k) = v(k,j+1,i) + v(k,j,i)
                flux_n(k) = ( v_comp(k) - gv ) * (                            &
                      37.0_wp * ( v(k,j+1,i) + v(k,j,i)   )                      &
                    -  8.0_wp * ( v(k,j+2,i) + v(k,j-1,i) )                      &
                      +         ( v(k,j+3,i) + v(k,j-2,i) ) ) * adv_mom_5

                diss_n(k) = - ABS( v_comp(k) - gv ) * (                       &
                      10.0_wp * ( v(k,j+1,i) - v(k,j,i)   )                      &
                    -  5.0_wp * ( v(k,j+2,i) - v(k,j-1,i) )                      &
                    +           ( v(k,j+3,i) - v(k,j-2,i) ) ) * adv_mom_5
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit26 = IBITS(wall_flags_0(k,j,i),26,1)
                ibit25 = IBITS(wall_flags_0(k,j,i),25,1)
                ibit24 = IBITS(wall_flags_0(k,j,i),24,1)

                k_ppp = k + 3 * ibit26
                k_pp  = k + 2 * ( 1 - ibit24  )
                k_mm  = k - 2 * ibit26

                w_comp    = w(k,j-1,i) + w(k,j,i)
                flux_t(k) = w_comp * rho_air_zw(k) * (                       &
                          ( 37.0_wp * ibit26 * adv_mom_5                        &
                       +     7.0_wp * ibit25 * adv_mom_3                        &
                       +              ibit24 * adv_mom_1                        &
                          ) *                                                &
                             ( v(k+1,j,i)   + v(k,j,i)    )                  &
                   -      (  8.0_wp * ibit26 * adv_mom_5                        &
                       +              ibit25 * adv_mom_3                        &
                          ) *                                                &
                             ( v(k_pp,j,i)  + v(k-1,j,i)  )                  &
                   +      (           ibit26 * adv_mom_5                        &
                          ) *                                                &
                             ( v(k_ppp,j,i) + v(k_mm,j,i) )                  &
                                      )

                diss_t(k) = - ABS( w_comp ) * rho_air_zw(k) * (              &
                          ( 10.0_wp * ibit26 * adv_mom_5                        &
                       +     3.0_wp * ibit25 * adv_mom_3                        &
                       +              ibit24 * adv_mom_1                        &
                          ) *                                                &
                             ( v(k+1,j,i)   - v(k,j,i)    )                  &
                   -      (  5.0_wp * ibit26 * adv_mom_5                        &
                       +              ibit25 * adv_mom_3                        &
                          ) *                                                &
                             ( v(k_pp,j,i)  - v(k-1,j,i)  )                  &
                   +      (           ibit26 * adv_mom_5                        &
                          ) *                                                &
                             ( v(k_ppp,j,i) - v(k_mm,j,i) )                  &
                                               )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( u_comp + gu - ( u(k,j-1,i)   + u(k,j,i)   ) ) * ddx &
                                                                  * rho_air(k)&
                     +  ( v_comp(k)   - ( v(k,j,i)     + v(k,j-1,i) ) ) * ddy &
                                                                  * rho_air(k)&
                     +  (   w_comp                      * rho_air_zw(k) -     &
                          ( w(k-1,j-1,i) + w(k-1,j,i) ) * rho_air_zw(k-1)     &
                        ) * ddzw(k)                                           &
                      ) * 0.5_wp
 
                tend(k,j,i) = tend(k,j,i) - (                                 &
                       ( flux_r(k) + diss_r(k)                                &
                     -   swap_flux_x_local_v(k,j) - swap_diss_x_local_v(k,j)  &
                       ) * ddx                                                &
                     + ( flux_n(k) + diss_n(k)                                &
                     -   swap_flux_y_local_v(k) - swap_diss_y_local_v(k)      &
                       ) * ddy                                                &
                     + ( ( flux_t(k) + diss_t(k) )                            &
                     -   ( flux_d    + diss_d    )                            &
                       ) * drho_air(k) * ddzw(k)                              &
                                            )  + v(k,j,i) * div

                swap_flux_x_local_v(k,j) = flux_r(k)
                swap_diss_x_local_v(k,j) = diss_r(k)
                swap_flux_y_local_v(k)   = flux_n(k)
                swap_diss_y_local_v(k)   = diss_n(k)
                flux_d                   = flux_t(k)
                diss_d                   = diss_t(k)

!
!--             Statistical Evaluation of v'v'. The factor has to be applied
!--             for right evaluation when gallilei_trans = .T. .
                sums_vs2_ws_l(k,tn) = sums_vs2_ws_l(k,tn)                      &
                + ( flux_n(k)                                                  &
                    * ( v_comp(k) - 2.0_wp * hom(k,1,2,0)                   )  &
                    / ( v_comp(k) - gv + SIGN( 1.0E-20_wp, v_comp(k) - gv ) )  &
               +   diss_n(k)                                                   &
                    *   ABS( v_comp(k) - 2.0_wp * hom(k,1,2,0)              )  &
                    / ( ABS( v_comp(k) - gv ) + 1.0E-20_wp                  )  &
                  ) *   weight_substep(intermediate_timestep_count)
!
!--             Statistical Evaluation of w'u'.
                sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn)                    &
                + ( flux_t(k)                                                  &
                    * ( w_comp - 2.0_wp * hom(k,1,3,0)                   )     &
                    / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              )     &
               +   diss_t(k)                                                   &
                    *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              )     &
                    / ( ABS( w_comp ) + 1.0E-20_wp                       )     &
                  ) *   weight_substep(intermediate_timestep_count)

             ENDDO
          ENDDO
       ENDDO
       sums_vs2_ws_l(nzb,tn) = sums_vs2_ws_l(nzb+1,tn)


    END SUBROUTINE advec_v_ws
    
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of v - Call for all grid points - accelerator version
!------------------------------------------------------------------------------!
    SUBROUTINE advec_v_ws_acc

       USE arrays_3d,                                                          &
           ONLY:  ddzw, drho_air, tend, u, v, w, rho_air, rho_air_zw

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count, u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  i_left, i_right, j_north, j_south, nxl, nxr, nyn, nys, nzb,  &
                  nzb_max, nzt, wall_flags_0
           
       USE kinds
       
!        USE statistics,                                                       &
!            ONLY:  hom, sums_vs2_ws_l, sums_wsvs_ws_l, weight_substep

       IMPLICIT NONE


       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit18 !<
       INTEGER(iwp) ::  ibit19 !<
       INTEGER(iwp) ::  ibit20 !<
       INTEGER(iwp) ::  ibit21 !<
       INTEGER(iwp) ::  ibit22 !<
       INTEGER(iwp) ::  ibit23 !<
       INTEGER(iwp) ::  ibit24 !<
       INTEGER(iwp) ::  ibit25 !<
       INTEGER(iwp) ::  ibit26 !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_mmm  !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<

       REAL(wp)    ::  diss_d   !<
       REAL(wp)    ::  diss_l   !<
       REAL(wp)    ::  diss_n   !<
       REAL(wp)    ::  diss_r   !<
       REAL(wp)    ::  diss_s   !<
       REAL(wp)    ::  diss_t   !<
       REAL(wp)    ::  div      !<
       REAL(wp)    ::  flux_d   !<
       REAL(wp)    ::  flux_l   !<
       REAL(wp)    ::  flux_n   !<
       REAL(wp)    ::  flux_r   !<
       REAL(wp)    ::  flux_s   !<
       REAL(wp)    ::  flux_t   !<
       REAL(wp)    ::  gu       !<
       REAL(wp)    ::  gv       !<
       REAL(wp)    ::  u_comp   !<
       REAL(wp)    ::  u_comp_l !<
       REAL(wp)    ::  v_comp   !<
       REAL(wp)    ::  v_comp_s !<
       REAL(wp)    ::  w_comp   !<

       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans

!
!--    Computation of fluxes and tendency terms
       !$acc kernels present( ddzw, tend, u, v, w, wall_flags_0 )
       DO  i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = nzb+1, nzt

                ibit20 = IBITS(wall_flags_0(k,j,i-1),20,1)
                ibit19 = IBITS(wall_flags_0(k,j,i-1),19,1)
                ibit18 = IBITS(wall_flags_0(k,j,i-1),18,1)

                u_comp_l                 = u(k,j-1,i) + u(k,j,i) - gu
                flux_l                   = u_comp_l * (                          &
                                      ( 37.0_wp * ibit20 * adv_mom_5              &
                                   +     7.0_wp * ibit19 * adv_mom_3              &
                                   +              ibit18 * adv_mom_1              &
                                      ) *                                      &
                                     ( v(k,j,i)   + v(k,j,i-1) )               &
                               -      (  8.0_wp * ibit20 * adv_mom_5              &
                                   +              ibit19 * adv_mom_3              &
                                      ) *                                      &
                                     ( v(k,j,i+1) + v(k,j,i-2) )               &
                               +      (           ibit20 * adv_mom_5              &
                                      ) *                                      &
                                     ( v(k,j,i+2) + v(k,j,i-3) )               &
                                                 )

                diss_l                   = - ABS( u_comp_l ) * (                 &
                                      ( 10.0_wp * ibit20 * adv_mom_5              &
                                   +     3.0_wp * ibit19 * adv_mom_3              &
                                   +              ibit18 * adv_mom_1              &
                                      ) *                                      &
                                     ( v(k,j,i)   - v(k,j,i-1) )               &
                               -      (  5.0_wp * ibit20 * adv_mom_5              &
                                   +              ibit19 * adv_mom_3              &
                                      ) *                                      &
                                     ( v(k,j,i+1) - v(k,j,i-2) )               &
                               +      (           ibit20 * adv_mom_5              &
                                      ) *                                      &
                                     ( v(k,j,i+2) - v(k,j,i-3) )               &
                                                           )

                ibit20 = IBITS(wall_flags_0(k,j,i),20,1)
                ibit19 = IBITS(wall_flags_0(k,j,i),19,1)
                ibit18 = IBITS(wall_flags_0(k,j,i),18,1)

                u_comp    = u(k,j-1,i+1) + u(k,j,i+1) - gu
                flux_r    = u_comp * (                                       &
                          ( 37.0_wp * ibit20 * adv_mom_5                        &
                       +     7.0_wp * ibit19 * adv_mom_3                        &
                       +              ibit18 * adv_mom_1                        &
                          ) *                                                &
                                 ( v(k,j,i+1) + v(k,j,i)   )                 &
                   -      (  8.0_wp * ibit20 * adv_mom_5                        &
                       +              ibit19 * adv_mom_3                        &
                          ) *                                                &
                                 ( v(k,j,i+2) + v(k,j,i-1) )                 &
                   +      (           ibit20 * adv_mom_5                        &
                          ) *                                                &
                                 ( v(k,j,i+3) + v(k,j,i-2) )                 &
                                     )

                diss_r    = - ABS( u_comp ) * (                              &
                          ( 10.0_wp * ibit20 * adv_mom_5                        &
                       +     3.0_wp * ibit19 * adv_mom_3                        &
                       +              ibit18 * adv_mom_1                        &
                          ) *                                                &
                                 ( v(k,j,i+1) - v(k,j,i)  )                  &
                   -      (  5.0_wp * ibit20 * adv_mom_5                        &
                       +              ibit19 * adv_mom_3                        &
                          ) *                                                &
                                 ( v(k,j,i+2) - v(k,j,i-1) )                 &
                   +      (           ibit20 * adv_mom_5                        &
                          ) *                                                &
                                 ( v(k,j,i+3) - v(k,j,i-2) )                 &
                                              )

                ibit23 = IBITS(wall_flags_0(k,j-1,i),23,1)
                ibit22 = IBITS(wall_flags_0(k,j-1,i),22,1)
                ibit21 = IBITS(wall_flags_0(k,j-1,i),21,1)


                v_comp_s              = v(k,j,i) + v(k,j-1,i) - gv
                flux_s                = v_comp_s    * (                       &
                                   ( 37.0_wp * ibit23 * adv_mom_5                &
                                +     7.0_wp * ibit22 * adv_mom_3                &
                                +              ibit21 * adv_mom_1                &
                                   ) *                                        &
                                     ( v(k,j,i)   + v(k,j-1,i) )              &
                            -      (  8.0_wp * ibit23 * adv_mom_5                &
                                +              ibit22 * adv_mom_3                &
                                   ) *                                        &
                                     ( v(k,j+1,i) + v(k,j-2,i) )              &
                            +      (           ibit23 * adv_mom_5                &
                                   ) *                                        &
                                     ( v(k,j+2,i) + v(k,j-3,i) )              &
                                                 )

                diss_s                = - ABS( v_comp_s ) * (                 &
                                   ( 10.0_wp * ibit23 * adv_mom_5                &
                                +     3.0_wp * ibit22 * adv_mom_3                &
                                +              ibit21 * adv_mom_1                &
                                   ) *                                        &
                                     ( v(k,j,i)   - v(k,j-1,i) )              &
                            -      (  5.0_wp * ibit23 * adv_mom_5                &
                                +              ibit22 * adv_mom_3                &
                                   ) *                                        &
                                     ( v(k,j+1,i) - v(k,j-2,i) )              &
                            +      (           ibit23 * adv_mom_5                &
                                   ) *                                        &
                                     ( v(k,j+2,i) - v(k,j-3,i) )              &
                                                          )

                ibit23 = IBITS(wall_flags_0(k,j,i),23,1)
                ibit22 = IBITS(wall_flags_0(k,j,i),22,1)
                ibit21 = IBITS(wall_flags_0(k,j,i),21,1)

                v_comp = v(k,j+1,i) + v(k,j,i)
                flux_n = ( v_comp - gv ) * (                                 &
                          ( 37.0_wp * ibit23 * adv_mom_5                        &
                       +     7.0_wp * ibit22 * adv_mom_3                        &
                       +              ibit21 * adv_mom_1                        &
                          ) *                                                &
                                 ( v(k,j+1,i) + v(k,j,i)   )                 &
                   -      (  8.0_wp * ibit23 * adv_mom_5                        &
                       +              ibit22 * adv_mom_3                        &
                          ) *                                                &
                                 ( v(k,j+2,i) + v(k,j-1,i) )                 &
                   +      (           ibit23 * adv_mom_5                        &
                          ) *                                                &
                                 ( v(k,j+3,i) + v(k,j-2,i) )                 &
                                     )

                diss_n = - ABS( v_comp - gv ) * (                         &
                          ( 10.0_wp * ibit23 * adv_mom_5                        &
                       +     3.0_wp * ibit22 * adv_mom_3                        &
                       +              ibit21 * adv_mom_1                        &
                          ) *                                                &
                                 ( v(k,j+1,i) - v(k,j,i)  )                  &
                   -      (  5.0_wp * ibit23 * adv_mom_5                        &
                       +              ibit22 * adv_mom_3                        &
                          ) *                                                &
                                 ( v(k,j+2,i) - v(k,j-1,i) )                 &
                   +      (           ibit23 * adv_mom_5                        &
                          ) *                                                &
                                 ( v(k,j+3,i) - v(k,j-2,i) )                 &
                                                     )

                ibit26 = IBITS(wall_flags_0(k-1,j,i),26,1)
                ibit25 = IBITS(wall_flags_0(k-1,j,i),25,1)
                ibit24 = IBITS(wall_flags_0(k-1,j,i),24,1)

                k_pp  = k + 2 * ibit26
                k_mm  = k - 2 * ( ibit25 + ibit26 )
                k_mmm = k - 3 * ibit26

                w_comp    = w(k-1,j-1,i) + w(k-1,j,i)
                flux_d    = w_comp  * (                                      &
                          ( 37.0_wp * ibit26 * adv_mom_5                        &
                       +     7.0_wp * ibit25 * adv_mom_3                        &
                       +              ibit24 * adv_mom_1                        &
                          ) *                                                &
                             ( v(k,j,i)     + v(k-1,j,i)  )                  &
                   -      (  8.0_wp * ibit26 * adv_mom_5                        &
                       +              ibit25 * adv_mom_3                        &
                          ) *                                                &
                             ( v(k+1,j,i)  + v(k_mm,j,i)  )                  &
                   +      (           ibit26 * adv_mom_5                        &
                          ) *                                                &
                             ( v(k_pp,j,i) + v(k_mmm,j,i) )                  &
                                      )

                diss_d    = - ABS( w_comp ) * (                              &
                          ( 10.0_wp * ibit26 * adv_mom_5                        &
                       +     3.0_wp * ibit25 * adv_mom_3                        &
                       +              ibit24 * adv_mom_1                        &
                          ) *                                                &
                             ( v(k,j,i)     - v(k-1,j,i)  )                  &
                   -      (  5.0_wp * ibit26 * adv_mom_5                        &
                       +              ibit25 * adv_mom_3                        &
                          ) *                                                &
                             ( v(k+1,j,i)  - v(k_mm,j,i)  )                  &
                   +      (           ibit26 * adv_mom_5                        &
                          ) *                                                &
                             ( v(k_pp,j,i) - v(k_mmm,j,i) )                  &
                                               )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit26 = IBITS(wall_flags_0(k,j,i),26,1)
                ibit25 = IBITS(wall_flags_0(k,j,i),25,1)
                ibit24 = IBITS(wall_flags_0(k,j,i),24,1)

                k_ppp = k + 3 * ibit26
                k_pp  = k + 2 * ( 1 - ibit24  )
                k_mm  = k - 2 * ibit26

                w_comp    = w(k,j-1,i) + w(k,j,i)
                flux_t    = w_comp * rho_air_zw(k) * (                       &
                          ( 37.0_wp * ibit26 * adv_mom_5                        &
                       +     7.0_wp * ibit25 * adv_mom_3                        &
                       +              ibit24 * adv_mom_1                        &
                          ) *                                                &
                             ( v(k+1,j,i)   + v(k,j,i)    )                  &
                   -      (  8.0_wp * ibit26 * adv_mom_5                        &
                       +              ibit25 * adv_mom_3                        &
                          ) *                                                &
                             ( v(k_pp,j,i)  + v(k-1,j,i)  )                  &
                   +      (           ibit26 * adv_mom_5                        &
                          ) *                                                &
                             ( v(k_ppp,j,i) + v(k_mm,j,i) )                  &
                                      )

                diss_t    = - ABS( w_comp ) * rho_air_zw(k) * (              &
                          ( 10.0_wp * ibit26 * adv_mom_5                        &
                       +     3.0_wp * ibit25 * adv_mom_3                        &
                       +              ibit24 * adv_mom_1                        &
                          ) *                                                &
                             ( v(k+1,j,i)   - v(k,j,i)    )                  &
                   -      (  5.0_wp * ibit26 * adv_mom_5                        &
                       +              ibit25 * adv_mom_3                        &
                          ) *                                                &
                             ( v(k_pp,j,i)  - v(k-1,j,i)  )                  &
                   +      (           ibit26 * adv_mom_5                        &
                          ) *                                                &
                             ( v(k_ppp,j,i) - v(k_mm,j,i) )                  &
                                               )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( ( u_comp     + gu )                                 &
                                       * ( ibit18 + ibit19 + ibit20 )         &
                - ( u(k,j-1,i)   + u(k,j,i) )                                 &
                                       * ( IBITS(wall_flags_0(k,j,i-1),18,1)  &
                                         + IBITS(wall_flags_0(k,j,i-1),19,1)  &
                                         + IBITS(wall_flags_0(k,j,i-1),20,1)  &
                                         )                                    &
                  ) * rho_air(k) * ddx                                        &
               +  ( v_comp                                                    &
                                       * ( ibit21 + ibit22 + ibit23 )         &
                - ( v(k,j,i)     + v(k,j-1,i) )                               &
                                       * ( IBITS(wall_flags_0(k,j-1,i),21,1)  &
                                         + IBITS(wall_flags_0(k,j-1,i),22,1)  &
                                         + IBITS(wall_flags_0(k,j-1,i),23,1)  &
                                         )                                    &
                  ) * rho_air(k) * ddy                                        &
               +  ( w_comp * rho_air_zw(k)                                    &
                                       * ( ibit24 + ibit25 + ibit26 )         &
                - ( w(k-1,j-1,i) + w(k-1,j,i) ) * rho_air_zw(k-1)             &
                                       * ( IBITS(wall_flags_0(k-1,j,i),24,1)  &
                                         + IBITS(wall_flags_0(k-1,j,i),25,1)  &
                                         + IBITS(wall_flags_0(k-1,j,i),26,1)  &
                                         )                                    &
                   ) * ddzw(k)   &
                ) * 0.5_wp


                tend(k,j,i) = - (                                              &
                               ( flux_r + diss_r - flux_l - diss_l ) * ddx     &
                             + ( flux_n + diss_n - flux_s - diss_s ) * ddy     &
                             + ( ( flux_t + diss_t ) -                         &
                                 ( flux_d + diss_d )                           &
                               ) * drho_air(k) * ddzw(k)                       &
                                ) + div * v(k,j,i)


!++
!--             Statistical Evaluation of v'v'. The factor has to be applied
!--             for right evaluation when gallilei_trans = .T. .
!                sums_vs2_ws_l(k,tn) = sums_vs2_ws_l(k,tn)                  &
!                      + ( flux_n                                           &
!                      * ( v_comp - 2.0_wp * hom(k,1,2,0) )                    &
!                      / ( v_comp - gv + 1.0E-20_wp )                       &
!                      +   diss_n                                           &
!                      *   ABS( v_comp - 2.0_wp * hom(k,1,2,0) )               &
!                      / ( ABS( v_comp - gv ) +1.0E-20_wp ) )               &
!                      *   weight_substep(intermediate_timestep_count)
!
!--              Statistical Evaluation of w'v'.
!                 sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn)                &
!                              + ( flux_t + diss_t )                         &
!                              *   weight_substep(intermediate_timestep_count)

             ENDDO
          ENDDO
       ENDDO
       !$acc end kernels

!++
!       sums_vs2_ws_l(nzb,tn) = sums_vs2_ws_l(nzb+1,tn)

    END SUBROUTINE advec_v_ws_acc
    
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of w - Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_w_ws

       USE arrays_3d,                                                          &
           ONLY:  ddzu, drho_air_zw, tend, u, v, w, rho_air, rho_air_zw

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count, u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzb_max, nzt, wall_flags_0,         &
                  wall_flags_00

       USE kinds
       
       USE statistics,                                                         &
           ONLY:  hom, sums_ws2_ws_l, weight_substep

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit27 !<
       INTEGER(iwp) ::  ibit28 !<
       INTEGER(iwp) ::  ibit29 !<
       INTEGER(iwp) ::  ibit30 !<
       INTEGER(iwp) ::  ibit31 !<
       INTEGER(iwp) ::  ibit32 !<
       INTEGER(iwp) ::  ibit33 !<
       INTEGER(iwp) ::  ibit34 !<
       INTEGER(iwp) ::  ibit35 !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<
       
       REAL(wp)    ::  diss_d !<
       REAL(wp)    ::  div    !<
       REAL(wp)    ::  flux_d !<
       REAL(wp)    ::  gu     !<
       REAL(wp)    ::  gv     !<
       REAL(wp)    ::  u_comp !<
       REAL(wp)    ::  v_comp !<
       REAL(wp)    ::  w_comp !<
       
       REAL(wp), DIMENSION(nzb:nzt)    ::  diss_t !<
       REAL(wp), DIMENSION(nzb:nzt)    ::  flux_t !<
       
       REAL(wp), DIMENSION(nzb+1:nzt)  ::  diss_n !<
       REAL(wp), DIMENSION(nzb+1:nzt)  ::  diss_r !<
       REAL(wp), DIMENSION(nzb+1:nzt)  ::  flux_n !<
       REAL(wp), DIMENSION(nzb+1:nzt)  ::  flux_r !<
       REAL(wp), DIMENSION(nzb+1:nzt)  ::  swap_diss_y_local_w !<
       REAL(wp), DIMENSION(nzb+1:nzt)  ::  swap_flux_y_local_w !<
       
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_diss_x_local_w !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_flux_x_local_w !<
 
       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans
!
!--   compute the whole left boundary of the processor domain
       i = nxl
       DO  j = nys, nyn
          DO  k = nzb+1, nzb_max

             ibit29 = IBITS(wall_flags_0(k,j,i-1),29,1)
             ibit28 = IBITS(wall_flags_0(k,j,i-1),28,1)
             ibit27 = IBITS(wall_flags_0(k,j,i-1),27,1)

             u_comp                   = u(k+1,j,i) + u(k,j,i) - gu
             swap_flux_x_local_w(k,j) = u_comp * (                             &
                                      ( 37.0_wp * ibit29 * adv_mom_5              &
                                   +     7.0_wp * ibit28 * adv_mom_3              &
                                   +              ibit27 * adv_mom_1              &
                                      ) *                                      &
                                     ( w(k,j,i)   + w(k,j,i-1) )               &
                               -      (  8.0_wp * ibit29 * adv_mom_5              &
                                   +              ibit28 * adv_mom_3              &
                                      ) *                                      &
                                     ( w(k,j,i+1) + w(k,j,i-2) )               &
                               +      (           ibit29 * adv_mom_5              &
                                      ) *                                      &
                                     ( w(k,j,i+2) + w(k,j,i-3) )               &
                                                 )

               swap_diss_x_local_w(k,j) = - ABS( u_comp ) * (                  &
                                        ( 10.0_wp * ibit29 * adv_mom_5            &
                                     +     3.0_wp * ibit28 * adv_mom_3            &
                                     +              ibit27 * adv_mom_1            &
                                        ) *                                    &
                                     ( w(k,j,i)   - w(k,j,i-1) )               &
                                 -      (  5.0_wp * ibit29 * adv_mom_5            &
                                     +              ibit28 * adv_mom_3            &
                                        ) *                                    &
                                     ( w(k,j,i+1) - w(k,j,i-2) )               &
                                 +      (           ibit29 * adv_mom_5            &
                                        ) *                                    &
                                     ( w(k,j,i+2) - w(k,j,i-3) )               &
                                                            )

          ENDDO

          DO  k = nzb_max+1, nzt

             u_comp                   = u(k+1,j,i) + u(k,j,i) - gu
             swap_flux_x_local_w(k,j) = u_comp * (                             &
                            37.0_wp * ( w(k,j,i) + w(k,j,i-1)   )                 &
                          -  8.0_wp * ( w(k,j,i+1) + w(k,j,i-2) )                 &
                          +           ( w(k,j,i+2) + w(k,j,i-3) ) ) * adv_mom_5
             swap_diss_x_local_w(k,j) = - ABS( u_comp ) * (                    &
                            10.0_wp * ( w(k,j,i) - w(k,j,i-1)   )                 &
                          -  5.0_wp * ( w(k,j,i+1) - w(k,j,i-2) )                 &
                          +           ( w(k,j,i+2) - w(k,j,i-3) ) ) * adv_mom_5

          ENDDO

       ENDDO

       DO i = nxl, nxr

          j = nys
          DO  k = nzb+1, nzb_max

             ibit32 = IBITS(wall_flags_00(k,j-1,i),0,1)
             ibit31 = IBITS(wall_flags_0(k,j-1,i),31,1)
             ibit30 = IBITS(wall_flags_0(k,j-1,i),30,1)

             v_comp                 = v(k+1,j,i) + v(k,j,i) - gv
             swap_flux_y_local_w(k) = v_comp * (                              &
                                    ( 37.0_wp * ibit32 * adv_mom_5               &
                                 +     7.0_wp * ibit31 * adv_mom_3               &
                                 +              ibit30 * adv_mom_1               &
                                    ) *                                        &
                                     ( w(k,j,i)   + w(k,j-1,i) )              &
                             -      (  8.0_wp * ibit32 * adv_mom_5               &
                                 +              ibit31 * adv_mom_3               &
                                    ) *                                       &
                                     ( w(k,j+1,i) + w(k,j-2,i) )              &
                             +      (           ibit32 * adv_mom_5               &
                                    ) *                                       &
                                     ( w(k,j+2,i) + w(k,j-3,i) )              &
                                               )

             swap_diss_y_local_w(k) = - ABS( v_comp ) * (                     &
                                    ( 10.0_wp * ibit32 * adv_mom_5               &
                                 +     3.0_wp * ibit31 * adv_mom_3               &
                                 +              ibit30 * adv_mom_1               &
                                    ) *                                       &
                                     ( w(k,j,i)   - w(k,j-1,i) )              &
                             -      (  5.0_wp * ibit32 * adv_mom_5               &
                                 +              ibit31 * adv_mom_3               &
                                    ) *                                       &
                                     ( w(k,j+1,i) - w(k,j-2,i) )              &
                             +      (           ibit32 * adv_mom_5               &
                                    ) *                                       &
                                     ( w(k,j+2,i) - w(k,j-3,i) )              &
                                                        )

          ENDDO

          DO  k = nzb_max+1, nzt

             v_comp                 = v(k+1,j,i) + v(k,j,i) - gv
             swap_flux_y_local_w(k) = v_comp * (                              &
                           37.0_wp * ( w(k,j,i) + w(k,j-1,i)   )                 &
                         -  8.0_wp * ( w(k,j+1,i) +w(k,j-2,i)  )                 &
                         +           ( w(k,j+2,i) + w(k,j-3,i) ) ) * adv_mom_5
             swap_diss_y_local_w(k) = - ABS( v_comp ) * (                     &
                           10.0_wp * ( w(k,j,i) - w(k,j-1,i)   )                 &
                         -  5.0_wp * ( w(k,j+1,i) - w(k,j-2,i) )                 &
                         +           ( w(k,j+2,i) - w(k,j-3,i) ) ) * adv_mom_5

          ENDDO

          DO  j = nys, nyn

!
!--          The lower flux has to be calculated explicetely for the tendency
!--          at the first w-level. For topography wall this is done implicitely
!--          by wall_flags_0.
             k         = nzb + 1
             w_comp    = w(k,j,i) + w(k-1,j,i)
             flux_t(0) = w_comp       * ( w(k,j,i) + w(k-1,j,i) ) * adv_mom_1
             diss_t(0) = -ABS(w_comp) * ( w(k,j,i) - w(k-1,j,i) ) * adv_mom_1
             flux_d    = flux_t(0)
             diss_d    = diss_t(0)

             DO  k = nzb+1, nzb_max

                ibit29 = IBITS(wall_flags_0(k,j,i),29,1)
                ibit28 = IBITS(wall_flags_0(k,j,i),28,1)
                ibit27 = IBITS(wall_flags_0(k,j,i),27,1)

                u_comp    = u(k+1,j,i+1) + u(k,j,i+1) - gu
                flux_r(k) = u_comp * (                                       &
                          ( 37.0_wp * ibit29 * adv_mom_5                        &
                       +     7.0_wp * ibit28 * adv_mom_3                        &
                       +              ibit27 * adv_mom_1                        &
                          ) *                                                &
                                 ( w(k,j,i+1) + w(k,j,i)   )                 &
                   -      (  8.0_wp * ibit29 * adv_mom_5                        &
                       +              ibit28 * adv_mom_3                        &
                          ) *                                                &
                                 ( w(k,j,i+2) + w(k,j,i-1) )                 &
                   +      (           ibit29 * adv_mom_5                        &
                          ) *                                                &
                                 ( w(k,j,i+3) + w(k,j,i-2) )                 &
                                     )

                diss_r(k) = - ABS( u_comp ) * (                              &
                          ( 10.0_wp * ibit29 * adv_mom_5                        &
                       +     3.0_wp * ibit28 * adv_mom_3                        &
                       +              ibit27 * adv_mom_1                        &
                          ) *                                                &
                                 ( w(k,j,i+1) - w(k,j,i)  )                  &
                   -      (  5.0_wp * ibit29 * adv_mom_5                        &
                       +              ibit28 * adv_mom_3                        &
                          ) *                                                &
                                 ( w(k,j,i+2) - w(k,j,i-1) )                 &
                   +      (           ibit29 * adv_mom_5                        &
                          ) *                                                &
                                 ( w(k,j,i+3) - w(k,j,i-2) )                 &
                                              )

                ibit32 = IBITS(wall_flags_00(k,j,i),0,1)
                ibit31 = IBITS(wall_flags_0(k,j,i),31,1)
                ibit30 = IBITS(wall_flags_0(k,j,i),30,1)

                v_comp    = v(k+1,j+1,i) + v(k,j+1,i) - gv
                flux_n(k) = v_comp * (                                       &
                          ( 37.0_wp * ibit32 * adv_mom_5                        &
                       +     7.0_wp * ibit31 * adv_mom_3                        &
                       +              ibit30 * adv_mom_1                        &
                          ) *                                                &
                                 ( w(k,j+1,i) + w(k,j,i)   )                 &
                   -      (  8.0_wp * ibit32 * adv_mom_5                        &
                       +              ibit31 * adv_mom_3                        &
                          ) *                                                 &
                                 ( w(k,j+2,i) + w(k,j-1,i) )                 &
                   +      (           ibit32 * adv_mom_5                        &
                          ) *                                                &
                                 ( w(k,j+3,i) + w(k,j-2,i) )                 &
                                     )

                diss_n(k) = - ABS( v_comp ) * (                              &
                          ( 10.0_wp * ibit32 * adv_mom_5                        &
                       +     3.0_wp * ibit31 * adv_mom_3                        &
                       +              ibit30 * adv_mom_1                        &
                          ) *                                                &
                                 ( w(k,j+1,i) - w(k,j,i)  )                  &
                   -      (  5.0_wp * ibit32 * adv_mom_5                        &
                       +              ibit31 * adv_mom_3                        &
                          ) *                                                &
                                 ( w(k,j+2,i) - w(k,j-1,i) )                 &
                   +      (           ibit32 * adv_mom_5                        &
                          ) *                                                &
                                 ( w(k,j+3,i) - w(k,j-2,i) )                 &
                                              )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit35 = IBITS(wall_flags_00(k,j,i),3,1)
                ibit34 = IBITS(wall_flags_00(k,j,i),2,1)
                ibit33 = IBITS(wall_flags_00(k,j,i),1,1)

                k_ppp = k + 3 * ibit35
                k_pp  = k + 2 * ( 1 - ibit33  )
                k_mm  = k - 2 * ibit35

                w_comp    = w(k+1,j,i) + w(k,j,i)
                flux_t(k) = w_comp * rho_air(k+1) * (                        &
                          ( 37.0_wp * ibit35 * adv_mom_5                        &
                       +     7.0_wp * ibit34 * adv_mom_3                        &
                       +              ibit33 * adv_mom_1                        &
                          ) *                                                &
                             ( w(k+1,j,i)  + w(k,j,i)     )                  &
                   -      (  8.0_wp * ibit35 * adv_mom_5                        &
                       +              ibit34 * adv_mom_3                        &
                          ) *                                                &
                             ( w(k_pp,j,i)  + w(k-1,j,i)  )                  &
                   +      (           ibit35 * adv_mom_5                        &
                          ) *                                                &
                             ( w(k_ppp,j,i) + w(k_mm,j,i) )                  &
                                       )

                diss_t(k) = - ABS( w_comp ) * rho_air(k+1) * (               &
                          ( 10.0_wp * ibit35 * adv_mom_5                        &
                       +     3.0_wp * ibit34 * adv_mom_3                        &
                       +              ibit33 * adv_mom_1                        &
                          ) *                                                &
                             ( w(k+1,j,i)   - w(k,j,i)    )                  &
                   -      (  5.0_wp * ibit35 * adv_mom_5                        &
                       +              ibit34 * adv_mom_3                        &
                          ) *                                                &
                             ( w(k_pp,j,i)  - w(k-1,j,i)  )                  &
                   +      (           ibit35 * adv_mom_5                        &
                          ) *                                                &
                             ( w(k_ppp,j,i) - w(k_mm,j,i) )                  &
                                               )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( ( u_comp + gu ) * ( ibit27 + ibit28 + ibit29 )      &
                  - ( u(k+1,j,i) + u(k,j,i)   )                               & 
                                    * ( IBITS(wall_flags_0(k,j,i-1),27,1)     &
                                      + IBITS(wall_flags_0(k,j,i-1),28,1)     &
                                      + IBITS(wall_flags_0(k,j,i-1),29,1)     &
                                      )                                       &
                  ) * rho_air_zw(k) * ddx                                     &
              +   ( ( v_comp + gv ) * ( ibit30 + ibit31 + ibit32 )            &
                  - ( v(k+1,j,i) + v(k,j,i)   )                               &
                                    * ( IBITS(wall_flags_0(k,j-1,i),30,1)     &
                                      + IBITS(wall_flags_0(k,j-1,i),31,1)     &
                                      + IBITS(wall_flags_00(k,j-1,i),0,1)     &
                                      )                                       &
                  ) * rho_air_zw(k) * ddy                                     &
              +   ( w_comp * rho_air(k+1) * ( ibit33 + ibit34 + ibit35 )      &
                - ( w(k,j,i)   + w(k-1,j,i)   ) * rho_air(k)                  &
                                    * ( IBITS(wall_flags_00(k-1,j,i),1,1)     &
                                      + IBITS(wall_flags_00(k-1,j,i),2,1)     &
                                      + IBITS(wall_flags_00(k-1,j,i),3,1)     &
                                      )                                       & 
                  ) * ddzu(k+1)   &
                ) * 0.5_wp



                tend(k,j,i) = tend(k,j,i) - (                                 &
                      ( flux_r(k) + diss_r(k)                                 &
                    -   swap_flux_x_local_w(k,j) - swap_diss_x_local_w(k,j)   &
                      ) * ddx                                                 &
                    + ( flux_n(k) + diss_n(k)                                 &
                    -   swap_flux_y_local_w(k)   - swap_diss_y_local_w(k)     &
                      ) * ddy                                                 &
                    + ( ( flux_t(k) + diss_t(k) )                             &
                    -   ( flux_d    + diss_d    )                             &
                      ) * drho_air_zw(k) * ddzu(k+1)                          &
                                            )  + div * w(k,j,i)

                swap_flux_x_local_w(k,j) = flux_r(k)
                swap_diss_x_local_w(k,j) = diss_r(k)
                swap_flux_y_local_w(k)   = flux_n(k)
                swap_diss_y_local_w(k)   = diss_n(k)
                flux_d                   = flux_t(k)
                diss_d                   = diss_t(k)

                sums_ws2_ws_l(k,tn)  = sums_ws2_ws_l(k,tn)                    &
                      + ( flux_t(k)                                           &
                       * ( w_comp - 2.0_wp * hom(k,1,3,0)                   ) &
                       / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              ) &
                        + diss_t(k)                                           &
                       *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              ) &
                       / ( ABS( w_comp ) + 1.0E-20_wp                       ) &
                        ) *   weight_substep(intermediate_timestep_count)

             ENDDO

             DO  k = nzb_max+1, nzt

                u_comp    = u(k+1,j,i+1) + u(k,j,i+1) - gu
                flux_r(k) = u_comp * (                                      &
                      37.0_wp * ( w(k,j,i+1) + w(k,j,i)   )                    &
                    -  8.0_wp * ( w(k,j,i+2) + w(k,j,i-1) )                    &
                    +           ( w(k,j,i+3) + w(k,j,i-2) ) ) * adv_mom_5

                diss_r(k) = - ABS( u_comp ) * (                             &
                      10.0_wp * ( w(k,j,i+1) - w(k,j,i)   )                    &
                    -  5.0_wp * ( w(k,j,i+2) - w(k,j,i-1) )                    &
                    +           ( w(k,j,i+3) - w(k,j,i-2) ) ) * adv_mom_5

                v_comp    = v(k+1,j+1,i) + v(k,j+1,i) - gv
                flux_n(k) = v_comp * (                                      &
                      37.0_wp * ( w(k,j+1,i) + w(k,j,i)   )                    &
                    -  8.0_wp * ( w(k,j+2,i) + w(k,j-1,i) )                    &
                    +           ( w(k,j+3,i) + w(k,j-2,i) ) ) * adv_mom_5

                diss_n(k) = - ABS( v_comp ) * (                             &
                      10.0_wp * ( w(k,j+1,i) - w(k,j,i)   )                    &
                    -  5.0_wp * ( w(k,j+2,i) - w(k,j-1,i) )                    &
                    +           ( w(k,j+3,i) - w(k,j-2,i) ) ) * adv_mom_5
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit35 = IBITS(wall_flags_00(k,j,i),3,1)
                ibit34 = IBITS(wall_flags_00(k,j,i),2,1)
                ibit33 = IBITS(wall_flags_00(k,j,i),1,1)

                k_ppp = k + 3 * ibit35
                k_pp  = k + 2 * ( 1 - ibit33  )
                k_mm  = k - 2 * ibit35

                w_comp    = w(k+1,j,i) + w(k,j,i)
                flux_t(k) = w_comp * rho_air(k+1) * (                        &
                          ( 37.0_wp * ibit35 * adv_mom_5                        &
                       +     7.0_wp * ibit34 * adv_mom_3                        &
                       +              ibit33 * adv_mom_1                        &
                          ) *                                                &
                             ( w(k+1,j,i)  + w(k,j,i)     )                  &
                   -      (  8.0_wp * ibit35 * adv_mom_5                        &
                       +              ibit34 * adv_mom_3                        &
                          ) *                                                &
                             ( w(k_pp,j,i)  + w(k-1,j,i)  )                  &
                   +      (           ibit35 * adv_mom_5                        &
                          ) *                                                &
                             ( w(k_ppp,j,i) + w(k_mm,j,i) )                  &
                                       )

                diss_t(k) = - ABS( w_comp ) * rho_air(k+1) * (               &
                          ( 10.0_wp * ibit35 * adv_mom_5                        &
                       +     3.0_wp * ibit34 * adv_mom_3                        &
                       +              ibit33 * adv_mom_1                        &
                          ) *                                                &
                             ( w(k+1,j,i)   - w(k,j,i)    )                  &
                   -      (  5.0_wp * ibit35 * adv_mom_5                        &
                       +              ibit34 * adv_mom_3                        &
                          ) *                                                &
                             ( w(k_pp,j,i)  - w(k-1,j,i)  )                  &
                   +      (           ibit35 * adv_mom_5                        &
                          ) *                                                &
                             ( w(k_ppp,j,i) - w(k_mm,j,i) )                  &
                                               )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( u_comp + gu - ( u(k+1,j,i) + u(k,j,i)   ) ) * ddx  &
                                                             * rho_air_zw(k) &
                    +   ( v_comp + gv - ( v(k+1,j,i) + v(k,j,i)   ) ) * ddy  &
                                                             * rho_air_zw(k) &
                    +   (   w_comp                    * rho_air(k+1) -       &
                          ( w(k,j,i)   + w(k-1,j,i) ) * rho_air(k)           &
                        ) * ddzu(k+1)                                        &
                      ) * 0.5_wp

                tend(k,j,i) = tend(k,j,i) - (                                 &
                      ( flux_r(k) + diss_r(k)                                 &
                    -   swap_flux_x_local_w(k,j) - swap_diss_x_local_w(k,j)   &
                      ) * ddx                                                 &
                    + ( flux_n(k) + diss_n(k)                                 &
                    -   swap_flux_y_local_w(k)   - swap_diss_y_local_w(k)     &
                      ) * ddy                                                 &
                    + ( ( flux_t(k) + diss_t(k) )                             &
                    -   ( flux_d    + diss_d    )                             &
                      ) * drho_air_zw(k) * ddzu(k+1)                          &
                                            )  + div * w(k,j,i)

                swap_flux_x_local_w(k,j) = flux_r(k)
                swap_diss_x_local_w(k,j) = diss_r(k)
                swap_flux_y_local_w(k)   = flux_n(k)
                swap_diss_y_local_w(k)   = diss_n(k)
                flux_d                   = flux_t(k)
                diss_d                   = diss_t(k)

                sums_ws2_ws_l(k,tn)  = sums_ws2_ws_l(k,tn)                    &
                      + ( flux_t(k)                                           &
                       * ( w_comp - 2.0_wp * hom(k,1,3,0)                   ) &
                       / ( w_comp + SIGN( 1.0E-20_wp, w_comp )              ) &
                        + diss_t(k)                                           &
                       *   ABS( w_comp - 2.0_wp * hom(k,1,3,0)              ) &
                       / ( ABS( w_comp ) + 1.0E-20_wp                       ) &
                        ) *   weight_substep(intermediate_timestep_count)

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE advec_w_ws


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of w - Call for all grid points - accelerator version
!------------------------------------------------------------------------------!
    SUBROUTINE advec_w_ws_acc

       USE arrays_3d,                                                          &
           ONLY:  ddzu, drho_air_zw, tend, u, v, w, rho_air, rho_air_zw

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count, u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  i_left, i_right, j_north, j_south, nxl, nxr, nyn, nys, nzb,  &
                  nzb_max, nzt, wall_flags_0, wall_flags_00
           
       USE kinds
       
!        USE statistics,                                                       &
!            ONLY:  hom, sums_ws2_ws_l, weight_substep

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit27 !<
       INTEGER(iwp) ::  ibit28 !<
       INTEGER(iwp) ::  ibit29 !<
       INTEGER(iwp) ::  ibit30 !<
       INTEGER(iwp) ::  ibit31 !<
       INTEGER(iwp) ::  ibit32 !<
       INTEGER(iwp) ::  ibit33 !<
       INTEGER(iwp) ::  ibit34 !<
       INTEGER(iwp) ::  ibit35 !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mmm  !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<

       REAL(wp)    ::  diss_d   !<
       REAL(wp)    ::  diss_l   !<
       REAL(wp)    ::  diss_n   !<
       REAL(wp)    ::  diss_r   !<
       REAL(wp)    ::  diss_s   !<
       REAL(wp)    ::  diss_t   !<
       REAL(wp)    ::  div      !<
       REAL(wp)    ::  flux_d   !<
       REAL(wp)    ::  flux_l   !<
       REAL(wp)    ::  flux_n   !<
       REAL(wp)    ::  flux_r   !<
       REAL(wp)    ::  flux_s   !<
       REAL(wp)    ::  flux_t   !<
       REAL(wp)    ::  gu       !<
       REAL(wp)    ::  gv       !<
       REAL(wp)    ::  u_comp   !<
       REAL(wp)    ::  u_comp_l !<
       REAL(wp)    ::  v_comp   !<
       REAL(wp)    ::  v_comp_s !<
       REAL(wp)    ::  w_comp   !<

       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans


!
!--    Computation of fluxes and tendency terms
       !$acc kernels present( ddzu, tend, u, v, w, wall_flags_0, wall_flags_00 )
       DO i = i_left, i_right
          DO  j = j_south, j_north
             DO  k = nzb+1, nzt

                ibit27 = IBITS(wall_flags_0(k,j,i-1),27,1)
                ibit28 = IBITS(wall_flags_0(k,j,i-1),28,1)
                ibit29 = IBITS(wall_flags_0(k,j,i-1),29,1)

                u_comp_l                 = u(k+1,j,i) + u(k,j,i) - gu
                flux_l                   = u_comp_l * (                        &
                                      ( 37.0_wp * ibit29 * adv_mom_5              &
                                   +     7.0_wp * ibit28 * adv_mom_3              &
                                   +              ibit27 * adv_mom_1              &
                                      ) *                                      &
                                     ( w(k,j,i)   + w(k,j,i-1) )               &
                               -      (  8.0_wp * ibit29 * adv_mom_5              &
                                   +              ibit28 * adv_mom_3              &
                                      ) *                                      &
                                     ( w(k,j,i+1) + w(k,j,i-2) )               &
                               +      (           ibit29 * adv_mom_5              &
                                      ) *                                      &
                                     ( w(k,j,i+2) + w(k,j,i-3) )               &
                                                 )

                diss_l                    = - ABS( u_comp_l ) * (              &
                                        ( 10.0_wp * ibit29 * adv_mom_5            &
                                     +     3.0_wp * ibit28 * adv_mom_3            &
                                     +              ibit27 * adv_mom_1            &
                                        ) *                                    &
                                     ( w(k,j,i)   - w(k,j,i-1) )               &
                                 -      (  5.0_wp * ibit29 * adv_mom_5            &
                                     +              ibit28 * adv_mom_3            &
                                        ) *                                    &
                                     ( w(k,j,i+1) - w(k,j,i-2) )               &
                                 +      (           ibit29 * adv_mom_5            &
                                        ) *                                    &
                                     ( w(k,j,i+2) - w(k,j,i-3) )               &
                                                            )

                ibit27 = IBITS(wall_flags_0(k,j,i),27,1)
                ibit28 = IBITS(wall_flags_0(k,j,i),28,1)
                ibit29 = IBITS(wall_flags_0(k,j,i),29,1)

                u_comp    = u(k+1,j,i+1) + u(k,j,i+1) - gu
                flux_r    = u_comp * (                                       &
                          ( 37.0_wp * ibit29 * adv_mom_5                        &
                       +     7.0_wp * ibit28 * adv_mom_3                        &
                       +              ibit27 * adv_mom_1                        &
                          ) *                                                &
                                 ( w(k,j,i+1) + w(k,j,i)   )                 &
                   -      (  8.0_wp * ibit29 * adv_mom_5                        &
                       +              ibit28 * adv_mom_3                        &
                          ) *                                                &
                                 ( w(k,j,i+2) + w(k,j,i-1) )                 &
                   +      (           ibit29 * adv_mom_5                        &
                          ) *                                                &
                                 ( w(k,j,i+3) + w(k,j,i-2) )                 &
                                     )

                diss_r    = - ABS( u_comp ) * (                              &
                          ( 10.0_wp * ibit29 * adv_mom_5                        &
                       +     3.0_wp * ibit28 * adv_mom_3                        &
                       +              ibit27 * adv_mom_1                        &
                          ) *                                                &
                                 ( w(k,j,i+1) - w(k,j,i)  )                  &
                   -      (  5.0_wp * ibit29 * adv_mom_5                        &
                       +              ibit28 * adv_mom_3                        &
                          ) *                                                &
                                 ( w(k,j,i+2) - w(k,j,i-1) )                 &
                   +      (           ibit29 * adv_mom_5                        &
                          ) *                                                &
                                 ( w(k,j,i+3) - w(k,j,i-2) )                 &
                                              )
                ibit32 = IBITS(wall_flags_00(k,j-1,i),0,1)
                ibit31 = IBITS(wall_flags_0(k,j-1,i),31,1)
                ibit30 = IBITS(wall_flags_0(k,j-1,i),30,1)

                v_comp_s               = v(k+1,j,i) + v(k,j,i) - gv
                flux_s                 = v_comp_s * (                         &
                                    ( 37.0_wp * ibit32 * adv_mom_5               &
                                 +     7.0_wp * ibit31 * adv_mom_3               &
                                 +              ibit30 * adv_mom_1               &
                                    ) *                                       &
                                     ( w(k,j,i)   + w(k,j-1,i) )              &
                             -      (  8.0_wp * ibit32 * adv_mom_5               &
                                 +              ibit31 * adv_mom_3               &
                                    ) *                                       &
                                     ( w(k,j+1,i) + w(k,j-2,i) )              &
                             +      (           ibit32 * adv_mom_5               &
                                    ) *                                       &
                                     ( w(k,j+2,i) + w(k,j-3,i) )              &
                                               )

                diss_s                 = - ABS( v_comp_s ) * (                &
                                    ( 10.0_wp * ibit32 * adv_mom_5               &
                                 +     3.0_wp * ibit31 * adv_mom_3               &
                                 +              ibit30 * adv_mom_1               &
                                    ) *                                       &
                                     ( w(k,j,i)   - w(k,j-1,i) )              &
                             -      (  5.0_wp * ibit32 * adv_mom_5               &
                                 +              ibit31 * adv_mom_3               &
                                    ) *                                       &
                                     ( w(k,j+1,i) - w(k,j-2,i) )              &
                             +      (           ibit32 * adv_mom_5               &
                                    ) *                                       &
                                     ( w(k,j+2,i) - w(k,j-3,i) )              &
                                                        )

                ibit32 = IBITS(wall_flags_00(k,j,i),0,1)
                ibit31 = IBITS(wall_flags_0(k,j,i),31,1)
                ibit30 = IBITS(wall_flags_0(k,j,i),30,1)

                v_comp    = v(k+1,j+1,i) + v(k,j+1,i) - gv
                flux_n    = v_comp * (                                       &
                          ( 37.0_wp * ibit32 * adv_mom_5                        &
                       +     7.0_wp * ibit31 * adv_mom_3                        &
                       +              ibit30 * adv_mom_1                        &
                          ) *                                                 &
                                 ( w(k,j+1,i) + w(k,j,i)   )                 &
                   -      (  8.0_wp * ibit32 * adv_mom_5                        &
                       +              ibit31 * adv_mom_3                        &
                          ) *                                                &
                                 ( w(k,j+2,i) + w(k,j-1,i) )                 &
                   +      (           ibit32 * adv_mom_5                        &
                          ) *                                                &
                                 ( w(k,j+3,i) + w(k,j-2,i) )                 &
                                     )

                diss_n    = - ABS( v_comp ) * (                              &
                          ( 10.0_wp * ibit32 * adv_mom_5                        &
                       +     3.0_wp * ibit31 * adv_mom_3                        &
                       +              ibit30 * adv_mom_1                        &
                          ) *                                                &
                                 ( w(k,j+1,i) - w(k,j,i)  )                  &
                   -      (  5.0_wp * ibit32 * adv_mom_5                        &
                       +              ibit31 * adv_mom_3                        &
                          ) *                                                &
                                 ( w(k,j+2,i) - w(k,j-1,i) )                 &
                   +      (           ibit32 * adv_mom_5                        &
                          ) *                                                &
                                 ( w(k,j+3,i) - w(k,j-2,i) )                 &
                                              )

                ibit35 = IBITS(wall_flags_00(k-1,j,i),3,1)
                ibit34 = IBITS(wall_flags_00(k-1,j,i),2,1)
                ibit33 = IBITS(wall_flags_00(k-1,j,i),1,1)

                k_pp  = k + 2 * ibit35
                k_mm  = k - 2 * ( ibit34 + ibit35 )
                k_mmm = k - 3 * ibit35

                w_comp    = w(k,j,i) + w(k-1,j,i)
                flux_d    = w_comp  * (                                      &
                          ( 37.0_wp * ibit35 * adv_mom_5                        &
                       +     7.0_wp * ibit34 * adv_mom_3                        &
                       +              ibit33 * adv_mom_1                        &
                          ) *                                                &
                             ( w(k,j,i)    + w(k-1,j,i)   )                  &
                   -      (  8.0_wp * ibit35 * adv_mom_5                        &
                       +              ibit34 * adv_mom_3                        &
                          ) *                                                &
                             ( w(k+1,j,i)  + w(k_mm,j,i)  )                  &
                   +      (           ibit35 * adv_mom_5                        &
                          ) *                                                &
                             ( w(k_pp,j,i) + w(k_mmm,j,i) )                  &
                                       )

                diss_d    = - ABS( w_comp ) * (                              &
                          ( 10.0_wp * ibit35 * adv_mom_5                        &
                       +     3.0_wp * ibit34 * adv_mom_3                        &
                       +              ibit33 * adv_mom_1                        &
                          ) *                                                &
                             ( w(k,j,i)    - w(k-1,j,i)   )                  &
                   -      (  5.0_wp * ibit35 * adv_mom_5                        &
                       +              ibit34 * adv_mom_3                        &
                          ) *                                                &
                             ( w(k+1,j,i)  - w(k_mm,j,i)  )                  &
                   +      (           ibit35 * adv_mom_5                        &
                          ) *                                                &
                             ( w(k_pp,j,i) - w(k_mmm,j,i) )                  &
                                               )

!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ibit35 = IBITS(wall_flags_00(k,j,i),3,1)
                ibit34 = IBITS(wall_flags_00(k,j,i),2,1)
                ibit33 = IBITS(wall_flags_00(k,j,i),1,1)

                k_ppp = k + 3 * ibit35
                k_pp  = k + 2 * ( 1 - ibit33  )
                k_mm  = k - 2 * ibit35

                w_comp    = w(k+1,j,i) + w(k,j,i)
                flux_t    = w_comp * rho_air(k+1) * (                        &
                          ( 37.0_wp * ibit35 * adv_mom_5                        &
                       +     7.0_wp * ibit34 * adv_mom_3                        &
                       +              ibit33 * adv_mom_1                        &
                          ) *                                                &
                             ( w(k+1,j,i)  + w(k,j,i)     )                  &
                   -      (  8.0_wp * ibit35 * adv_mom_5                        &
                       +              ibit34 * adv_mom_3                        &
                          ) *                                                &
                             ( w(k_pp,j,i)  + w(k-1,j,i)  )                  &
                   +      (           ibit35 * adv_mom_5                        &
                          ) *                                                &
                             ( w(k_ppp,j,i) + w(k_mm,j,i) )                  &
                                       )

                diss_t    = - ABS( w_comp ) * rho_air(k+1) * (               &
                          ( 10.0_wp * ibit35 * adv_mom_5                        &
                       +     3.0_wp * ibit34 * adv_mom_3                        &
                       +              ibit33 * adv_mom_1                        &
                          ) *                                                &
                             ( w(k+1,j,i)   - w(k,j,i)    )                  &
                   -      (  5.0_wp * ibit35 * adv_mom_5                        &
                       +              ibit34 * adv_mom_3                        &
                          ) *                                                &
                             ( w(k_pp,j,i)  - w(k-1,j,i)  )                  &
                   +      (           ibit35 * adv_mom_5                        &
                          ) *                                                &
                             ( w(k_ppp,j,i) - w(k_mm,j,i) )                  &
                                               )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( ( u_comp + gu ) * ( ibit27 + ibit28 + ibit29 )      &
                  - ( u(k+1,j,i) + u(k,j,i)   )                               & 
                                    * ( IBITS(wall_flags_0(k,j,i-1),27,1)     &
                                      + IBITS(wall_flags_0(k,j,i-1),28,1)     &
                                      + IBITS(wall_flags_0(k,j,i-1),29,1)     &
                                      )                                       &
                  ) * rho_air_zw(k) * ddx                                     &
              +   ( ( v_comp + gv ) * ( ibit30 + ibit31 + ibit32 )            &
                  - ( v(k+1,j,i) + v(k,j,i)   )                               &
                                    * ( IBITS(wall_flags_0(k,j-1,i),30,1)     &
                                      + IBITS(wall_flags_0(k,j-1,i),31,1)     &
                                      + IBITS(wall_flags_00(k,j-1,i),0,1)     &
                                      )                                       &
                  ) * rho_air_zw(k) * ddy                                     &
              +   ( w_comp * rho_air(k+1) * ( ibit33 + ibit34 + ibit35 )      &
                - ( w(k,j,i)   + w(k-1,j,i)   ) * rho_air(k)                  &
                                    * ( IBITS(wall_flags_00(k-1,j,i),1,1)     &
                                      + IBITS(wall_flags_00(k-1,j,i),2,1)     &
                                      + IBITS(wall_flags_00(k-1,j,i),3,1)     &
                                      )                                       & 
                  ) * ddzu(k+1)   &
                ) * 0.5_wp


                tend(k,j,i) = - (                                                &
                               ( flux_r + diss_r - flux_l - diss_l ) * ddx       &
                             + ( flux_n + diss_n - flux_s - diss_s ) * ddy       &
                             + ( ( flux_t + diss_t ) -                           &
                                 ( flux_d + diss_d ) * rho_air(k)                &
                               ) * drho_air_zw(k) * ddzu(k+1)                    &
                                 ) + div * w(k,j,i)


!++
!--             Statistical Evaluation of w'w'.
!                sums_ws2_ws_l(k,tn)  = sums_ws2_ws_l(k,tn)                    &
!                               + ( flux_t + diss_t )                    &
!                               *   weight_substep(intermediate_timestep_count)

             ENDDO
          ENDDO
       ENDDO
       !$acc end kernels

    END SUBROUTINE advec_w_ws_acc

 END MODULE advec_ws
