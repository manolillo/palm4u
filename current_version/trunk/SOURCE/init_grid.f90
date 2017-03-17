!> @file init_grid.f90
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
! $Id: init_grid.f90 2038 2016-10-26 11:16:56Z knoop $
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
! 
! 2021 2016-10-07 14:08:57Z suehring
! Bugfix: setting Neumann boundary conditions for topography required for 
! topography flags in multigrid_noopt solver
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1994 2016-08-15 09:52:21Z suehring
! Bugfix in definition of generic topography
! 
! 1982 2016-08-01 11:04:48Z suehring
! Bugfix concering consistency check for topography
! 
! 1968 2016-07-18 12:01:49Z suehring
! Changed: PE-wise reading of topography file in order to avoid global definition
! of arrays nzb_local and nzb_tmp. Thereby, topography definition for single
! buildings and street canyons has changed, as well as flag setting for 
! multigrid scheme.
!
! Bugfix in checking l_grid anisotropy. 
! Simplify initial computation of lwall and vertical_influence, i.e. remove 
! nzb_s_inner as it is still zero at this point.
! 
! 1942 2016-06-14 12:18:18Z suehring
! Topography filter implemented to fill holes resolved by only one grid point.
! Initialization of flags for ws-scheme moved to advec_ws.  
! 
! 1931 2016-06-10 12:06:59Z suehring
! Rename multigrid into multigrid_noopt and multigrid_fast into multigrid
!
! 1910 2016-05-26 06:49:46Z raasch
! Bugfix: if topography is read from file, Neumann conditions are used for the
! nzb_local array (instead of cyclic conditions) in case that non-cyclic
! boundary conditions are switched on for the run
!
! 1902 2016-05-09 11:18:56Z suehring
! Set topography flags for multigrid solver only (not for multigrid_fast)
!
! 1886 2016-04-21 11:20:47Z suehring
! Bugfix: setting advection flags near walls
! reformulated index values for nzb_v_inner
! variable discriptions added in declaration block
!
! 1845 2016-04-08 08:29:13Z raasch
! nzb_2d removed
! 
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1779 2016-03-03 08:01:28Z raasch
! coupling_char is trimmed at every place it occurs, because it can have
! different length now
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1743 2016-01-13 10:23:51Z raasch
! Bugfix for calculation of nzb_s_outer and nzb_u_outer at north boundary of
! total domain
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer.
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1677 2015-10-02 13:25:23Z boeske
! Bugfix: Ghost points are included in wall_flags_0 and wall_flags_00
! 
! 1675 2015-10-02 08:28:59Z gronemeier
! Bugfix: Definition of topography grid levels
! 
! 1660 2015-09-21 08:15:16Z gronemeier
! Bugfix: Definition of topography grid levels if vertical grid stretching
!         starts below the maximum topography height.
!
! 1580 2015-04-10 13:43:49Z suehring
! Bugfix: setting flags for 5th order scheme near buildings
!
! 1575 2015-03-27 09:56:27Z raasch
! adjustments for psolver-queries
!
! 1557 2015-03-05 16:43:04Z suehring
! Adjustment for monotoinic limiter
!
! 1418 2014-06-06 13:05:08Z fricke
! Bugfix: Change if-condition for stretched grid in the ocean, with the old
!          condition and a negative value for dz_stretch_level the condition
!          was always true for the whole model domain
!
! 1409 2014-05-23 12:11:32Z suehring
! Bugfix: set wall_flags_0 at inflow and outflow boundary also for i <= nxlu 
! j <= nysv
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
!
! 1322 2014-03-20 16:38:49Z raasch
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
! 1221 2013-09-10 08:59:13Z raasch
! wall_flags_00 introduced to hold bits 32-63,
! additional 3D-flag arrays for replacing the 2D-index array nzb_s_inner in
! loops optimized for openACC (pres + flow_statistics)
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1069 2012-11-28 16:18:43Z maronga
! bugfix: added coupling_char to TOPOGRAPHY_DATA to allow topography in the
!         ocean model in case of coupled runs
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! lower index for calculating wall_flags_0 set to nzb_w_inner instead of
! nzb_w_inner+1
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting
!
! 978 2012-08-09 08:28:32Z fricke
! Bugfix: nzb_max is set to nzt at non-cyclic lateral boundaries
! Bugfix: Set wall_flags_0 for inflow boundary
!
! 927 2012-06-06 19:15:04Z raasch
! Wall flags are not set for multigrid method in case of masking method
!
! 864 2012-03-27 15:10:33Z gryschka 
! In case of ocean and Dirichlet bottom bc for u and v dzu_mg and ddzu_pres
! were not correctly defined for k=1.
!
! 861 2012-03-26 14:18:34Z suehring
! Set wall_flags_0. The array is needed for degradation in ws-scheme near walls,
! inflow and outflow boundaries as well as near the bottom and the top of the
! model domain.!
! Initialization of nzb_s_inner and nzb_w_inner.
! gls has to be at least nbgp to do not exceed the array bounds of nzb_local
! while setting wall_flags_0
!
! 843 2012-02-29 15:16:21Z gryschka
! In case of ocean and dirichlet bc for u and v at the bottom
! the first u-level ist defined at same height as the first w-level
!
! 818 2012-02-08 16:11:23Z maronga
! Bugfix: topo_height is only required if topography is used. It is thus now
! allocated in the topography branch
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! Revision 1.1  1997/08/11 06:17:45  raasch
! Initial revision (Testversion)
!
!
! Description:
! ------------
!> Creating grid depending constants
!------------------------------------------------------------------------------!
 SUBROUTINE init_grid
 
    USE advec_ws,                                                              &
        ONLY:  ws_init_flags

    USE arrays_3d,                                                             &
        ONLY:  dd2zu, ddzu, ddzu_pres, ddzw, dzu, dzu_mg, dzw, dzw_mg, f1_mg,  &
               f2_mg, f3_mg, l_grid, l_wall, zu, zw
        
    USE control_parameters,                                                    &
        ONLY:  bc_lr_cyc, bc_ns_cyc, building_height, building_length_x,       &
               building_length_y, building_wall_left, building_wall_south,     &
               canyon_height, canyon_wall_left, canyon_wall_south,             &
               canyon_width_x, canyon_width_y, constant_flux_layer,            &
               coupling_char, dp_level_ind_b, dz, dz_max, dz_stretch_factor,   &
               dz_stretch_level, dz_stretch_level_index, grid_level, ibc_uv_b, &
               io_blocks, io_group, inflow_l, inflow_n, inflow_r, inflow_s,    &
               masking_method, maximum_grid_level, message_string,             &
               momentum_advec, nest_domain, nest_bound_l, nest_bound_n,        &
               nest_bound_r, nest_bound_s, ocean, outflow_l, outflow_n,        &
               outflow_r, outflow_s, psolver, scalar_advec, topography,        &
               topography_grid_convention, use_surface_fluxes, use_top_fluxes, &
               wall_adjustment_factor
         
    USE grid_variables,                                                        &
        ONLY:  ddx, ddx2, ddy, ddy2, dx, dx2, dy, dy2, fwxm,                   &
               fwxp, fwym, fwyp, fxm, fxp, fym, fyp, wall_e_x, wall_e_y,       &
               wall_u, wall_v, wall_w_x, wall_w_y, zu_s_inner, zw_w_inner
        
    USE indices,                                                               &
        ONLY:  flags, nbgp, nx, nxl, nxlg, nxl_mg, nxr, nxrg, nxr_mg,          &
               ny, nyn, nyng, nyn_mg, nys, nys_mg, nysg, nz, nzb,              &
               nzb_diff, nzb_diff_s_inner, nzb_diff_s_outer, nzb_diff_u,       &
               nzb_diff_v, nzb_max, nzb_s_inner, nzb_s_outer, nzb_u_inner,     &
               nzb_u_outer, nzb_v_inner, nzb_v_outer, nzb_w_inner,             &
               nzb_w_outer, nzt, nzt_diff, nzt_mg, rflags_invers,              &
               rflags_s_inner, wall_flags_0, wall_flags_00, wall_flags_1,      &
               wall_flags_10, wall_flags_2, wall_flags_3,  wall_flags_4,       &
               wall_flags_5, wall_flags_6, wall_flags_7, wall_flags_8,         &
               wall_flags_9
    
    USE kinds
    
    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  bh            !< temporary vertical index of building height
    INTEGER(iwp) ::  blx           !< grid point number of building size along x
    INTEGER(iwp) ::  bly           !< grid point number of building size along y
    INTEGER(iwp) ::  bxl           !< index for left building wall
    INTEGER(iwp) ::  bxr           !< index for right building wall
    INTEGER(iwp) ::  byn           !< index for north building wall
    INTEGER(iwp) ::  bys           !< index for south building wall
    INTEGER(iwp) ::  ch            !< temporary vertical index for canyon height
    INTEGER(iwp) ::  cwx           !< grid point number of canyon size along x
    INTEGER(iwp) ::  cwy           !< grid point number of canyon size along y
    INTEGER(iwp) ::  cxl           !< index for left canyon wall
    INTEGER(iwp) ::  cxr           !< index for right canyon wall
    INTEGER(iwp) ::  cyn           !< index for north canyon wall
    INTEGER(iwp) ::  cys           !< index for south canyon wall
    INTEGER(iwp) ::  i             !< index variable along x
    INTEGER(iwp) ::  ii            !< loop variable for reading topography file
    INTEGER(iwp) ::  inc           !< incremental parameter for coarsening grid level
    INTEGER(iwp) ::  j             !< index variable along y
    INTEGER(iwp) ::  k             !< index variable along z
    INTEGER(iwp) ::  l             !< loop variable 
    INTEGER(iwp) ::  nxl_l         !< index of left PE boundary for multigrid level
    INTEGER(iwp) ::  nxr_l         !< index of right PE boundary for multigrid level
    INTEGER(iwp) ::  nyn_l         !< index of north PE boundary for multigrid level
    INTEGER(iwp) ::  nys_l         !< index of south PE boundary for multigrid level
    INTEGER(iwp) ::  nzb_local_max !< vertical grid index of maximum topography height
    INTEGER(iwp) ::  nzb_local_min !< vertical grid index of minimum topography height
    INTEGER(iwp) ::  nzb_si        !< dummy index for local nzb_s_inner
    INTEGER(iwp) ::  nzt_l         !< index of top PE boundary for multigrid level
    INTEGER(iwp) ::  num_hole      !< number of holes (in topography) resolved by only one grid point 
    INTEGER(iwp) ::  num_hole_l    !< number of holes (in topography) resolved by only one grid point on local PE     
    INTEGER(iwp) ::  num_wall      !< number of surrounding vertical walls for a single grid point
    INTEGER(iwp) ::  skip_n_rows   !< counting variable to skip rows while reading topography file    
    INTEGER(iwp) ::  vi            !< dummy for vertical influence

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE   ::                               &
                     vertical_influence  !< number of vertical grid points above obstacle where adjustment of near-wall mixing length is required
                                         
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  corner_nl      !< index of north-left corner location to limit near-wall mixing length
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  corner_nr      !< north-right
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  corner_sl      !< south-left
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  corner_sr      !< south-right
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nzb_local      !< index for topography top at cell-center
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nzb_tmp        !< dummy to calculate topography indices on u- and v-grid
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  wall_l         !< distance to adjacent left-facing wall
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  wall_n         !< north-facing
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  wall_r         !< right-facing
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  wall_s         !< right-facing

    REAL(wp) ::  dum           !< dummy variable to skip columns while reading topography file    
    REAL(wp) ::  dz_stretched  !< stretched vertical grid spacing

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  topo_height   !< input variable for topography height
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zu_s_inner_l  !< dummy array on global scale to write topography output array
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  zw_w_inner_l  !< dummy array on global scale to write topography output array

    
!
!-- Calculation of horizontal array bounds including ghost layers
    nxlg = nxl - nbgp
    nxrg = nxr + nbgp
    nysg = nys - nbgp
    nyng = nyn + nbgp

!
!-- Allocate grid arrays
    ALLOCATE( ddzu(1:nzt+1), ddzw(1:nzt+1), dd2zu(1:nzt), dzu(1:nzt+1),        &
              dzw(1:nzt+1), l_grid(1:nzt), zu(nzb:nzt+1), zw(nzb:nzt+1) )

!
!-- Compute height of u-levels from constant grid length and dz stretch factors
    IF ( dz == -1.0_wp )  THEN
       message_string = 'missing dz'
       CALL message( 'init_grid', 'PA0200', 1, 2, 0, 6, 0 ) 
    ELSEIF ( dz <= 0.0_wp )  THEN
       WRITE( message_string, * ) 'dz=',dz,' <= 0.0'
       CALL message( 'init_grid', 'PA0201', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Define the vertical grid levels
    IF ( .NOT. ocean )  THEN
!
!--    Grid for atmosphere with surface at z=0 (k=0, w-grid).
!--    The second u-level (k=1) corresponds to the top of the 
!--    Prandtl-layer.

       IF ( ibc_uv_b == 0 .OR. ibc_uv_b == 2 ) THEN
          zu(0) = 0.0_wp
      !    zu(0) = - dz * 0.5_wp
       ELSE
          zu(0) = - dz * 0.5_wp
       ENDIF
       zu(1) =   dz * 0.5_wp

       dz_stretch_level_index = nzt+1
       dz_stretched = dz
       DO  k = 2, nzt+1
          IF ( dz_stretch_level <= zu(k-1)  .AND.  dz_stretched < dz_max )  THEN
             dz_stretched = dz_stretched * dz_stretch_factor
             dz_stretched = MIN( dz_stretched, dz_max )
             IF ( dz_stretch_level_index == nzt+1 ) dz_stretch_level_index = k-1
          ENDIF
          zu(k) = zu(k-1) + dz_stretched
       ENDDO

!
!--    Compute the w-levels. They are always staggered half-way between the 
!--    corresponding u-levels. In case of dirichlet bc for u and v at the 
!--    ground the first u- and w-level (k=0) are defined at same height (z=0). 
!--    The top w-level is extrapolated linearly.
       zw(0) = 0.0_wp
       DO  k = 1, nzt
          zw(k) = ( zu(k) + zu(k+1) ) * 0.5_wp
       ENDDO
       zw(nzt+1) = zw(nzt) + 2.0_wp * ( zu(nzt+1) - zw(nzt) )

    ELSE
!
!--    Grid for ocean with free water surface is at k=nzt (w-grid).
!--    In case of neumann bc at the ground the first first u-level (k=0) lies
!--    below the first w-level (k=0). In case of dirichlet bc the first u- and
!--    w-level are defined at same height, but staggered from the second level.
!--    The second u-level (k=1) corresponds to the top of the Prandtl-layer.
       zu(nzt+1) =   dz * 0.5_wp
       zu(nzt)   = - dz * 0.5_wp

       dz_stretch_level_index = 0
       dz_stretched = dz
       DO  k = nzt-1, 0, -1
!
!--       The default value of dz_stretch_level is positive, thus the first 
!--       condition is always true. Hence, the second condition is necessary.
          IF ( dz_stretch_level >= zu(k+1)  .AND.  dz_stretch_level <= 0.0  &
               .AND.  dz_stretched < dz_max )  THEN
             dz_stretched = dz_stretched * dz_stretch_factor
             dz_stretched = MIN( dz_stretched, dz_max )
             IF ( dz_stretch_level_index == 0 ) dz_stretch_level_index = k+1
          ENDIF
          zu(k) = zu(k+1) - dz_stretched
       ENDDO

!
!--    Compute the w-levels. They are always staggered half-way between the 
!--    corresponding u-levels, except in case of dirichlet bc for u and v
!--    at the ground. In this case the first u- and w-level are defined at 
!--    same height. The top w-level (nzt+1) is not used but set for 
!--    consistency, since w and all scalar variables are defined up tp nzt+1.
       zw(nzt+1) = dz
       zw(nzt)   = 0.0_wp
       DO  k = 0, nzt
          zw(k) = ( zu(k) + zu(k+1) ) * 0.5_wp
       ENDDO

!
!--    In case of dirichlet bc for u and v the first u- and w-level are defined
!--    at same height.
       IF ( ibc_uv_b == 0 ) THEN
          zu(0) = zw(0)
       ENDIF

    ENDIF

!
!-- Compute grid lengths.
    DO  k = 1, nzt+1
       dzu(k)  = zu(k) - zu(k-1)
       ddzu(k) = 1.0_wp / dzu(k)
       dzw(k)  = zw(k) - zw(k-1)
       ddzw(k) = 1.0_wp / dzw(k)
    ENDDO

    DO  k = 1, nzt
       dd2zu(k) = 1.0_wp / ( dzu(k) + dzu(k+1) )
    ENDDO
    
!    
!-- The FFT- SOR-pressure solvers assume grid spacings of a staggered grid
!-- everywhere. For the actual grid, the grid spacing at the lowest level
!-- is only dz/2, but should be dz. Therefore, an additional array
!-- containing with appropriate grid information is created for these
!-- solvers.
    IF ( psolver(1:9) /= 'multigrid' )  THEN
       ALLOCATE( ddzu_pres(1:nzt+1) )
       ddzu_pres = ddzu
       ddzu_pres(1) = ddzu_pres(2)  ! change for lowest level
    ENDIF

!
!-- Compute the reciprocal values of the horizontal grid lengths.
    ddx = 1.0_wp / dx
    ddy = 1.0_wp / dy
    dx2 = dx * dx
    dy2 = dy * dy
    ddx2 = 1.0_wp / dx2
    ddy2 = 1.0_wp / dy2

!
!-- Compute the grid-dependent mixing length.
    DO  k = 1, nzt
       l_grid(k)  = ( dx * dy * dzw(k) )**0.33333333333333_wp
    ENDDO

!
!-- Allocate outer and inner index arrays for topography and set
!-- defaults.

    ALLOCATE( corner_nl(nys:nyn,nxl:nxr), corner_nr(nys:nyn,nxl:nxr),       &
              corner_sl(nys:nyn,nxl:nxr), corner_sr(nys:nyn,nxl:nxr),       &
              wall_l(nys:nyn,nxl:nxr), wall_n(nys:nyn,nxl:nxr),             &
              wall_r(nys:nyn,nxl:nxr), wall_s(nys:nyn,nxl:nxr) )                     
     
    ALLOCATE( fwxm(nysg:nyng,nxlg:nxrg), fwxp(nysg:nyng,nxlg:nxrg),         &
              fwym(nysg:nyng,nxlg:nxrg), fwyp(nysg:nyng,nxlg:nxrg),         &
              fxm(nysg:nyng,nxlg:nxrg), fxp(nysg:nyng,nxlg:nxrg),           &
              fym(nysg:nyng,nxlg:nxrg), fyp(nysg:nyng,nxlg:nxrg),           &
              nzb_s_inner(nysg:nyng,nxlg:nxrg),                             &
              nzb_s_outer(nysg:nyng,nxlg:nxrg),                             &
              nzb_u_inner(nysg:nyng,nxlg:nxrg),                             &
              nzb_u_outer(nysg:nyng,nxlg:nxrg),                             &
              nzb_v_inner(nysg:nyng,nxlg:nxrg),                             &
              nzb_v_outer(nysg:nyng,nxlg:nxrg),                             &
              nzb_w_inner(nysg:nyng,nxlg:nxrg),                             &
              nzb_w_outer(nysg:nyng,nxlg:nxrg),                             &
              nzb_diff_s_inner(nysg:nyng,nxlg:nxrg),                        &
              nzb_diff_s_outer(nysg:nyng,nxlg:nxrg),                        &
              nzb_diff_u(nysg:nyng,nxlg:nxrg),                              &
              nzb_diff_v(nysg:nyng,nxlg:nxrg),                              &
              nzb_local(nysg:nyng,nxlg:nxrg),                               &
              nzb_tmp(nysg:nyng,nxlg:nxrg),                                 &
              rflags_s_inner(nzb:nzt+2,nysg:nyng,nxlg:nxrg),                &
              rflags_invers(nysg:nyng,nxlg:nxrg,nzb:nzt+2),                 &
              wall_e_x(nysg:nyng,nxlg:nxrg),                                &
              wall_e_y(nysg:nyng,nxlg:nxrg),                                &
              wall_u(nysg:nyng,nxlg:nxrg),                                  &
              wall_v(nysg:nyng,nxlg:nxrg),                                  &
              wall_w_x(nysg:nyng,nxlg:nxrg),                                &
              wall_w_y(nysg:nyng,nxlg:nxrg) )



    ALLOCATE( l_wall(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )


    nzb_s_inner = nzb;  nzb_s_outer = nzb
    nzb_u_inner = nzb;  nzb_u_outer = nzb
    nzb_v_inner = nzb;  nzb_v_outer = nzb
    nzb_w_inner = nzb;  nzb_w_outer = nzb

    rflags_s_inner = 1.0_wp
    rflags_invers  = 1.0_wp

!
!-- Define vertical gridpoint from (or to) which on the usual finite difference
!-- form (which does not use surface fluxes) is applied 
    IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
       nzb_diff = nzb + 2
    ELSE
       nzb_diff = nzb + 1
    ENDIF
    IF ( use_top_fluxes )  THEN
       nzt_diff = nzt - 1
    ELSE
       nzt_diff = nzt
    ENDIF

    nzb_diff_s_inner = nzb_diff;  nzb_diff_s_outer = nzb_diff
    nzb_diff_u = nzb_diff;  nzb_diff_v = nzb_diff

    wall_e_x = 0.0_wp;  wall_e_y = 0.0_wp;  wall_u = 0.0_wp;  wall_v = 0.0_wp
    wall_w_x = 0.0_wp;  wall_w_y = 0.0_wp
    fwxp = 1.0_wp;  fwxm = 1.0_wp;  fwyp = 1.0_wp;  fwym = 1.0_wp
    fxp  = 1.0_wp;  fxm  = 1.0_wp;  fyp  = 1.0_wp;  fym  = 1.0_wp

!
!-- Initialize near-wall mixing length l_wall only in the vertical direction
!-- for the moment,
!-- multiplication with wall_adjustment_factor near the end of this routine
    l_wall(nzb,:,:)   = l_grid(1)
    DO  k = nzb+1, nzt
       l_wall(k,:,:)  = l_grid(k)
    ENDDO
    l_wall(nzt+1,:,:) = l_grid(nzt)

    ALLOCATE ( vertical_influence(nzb:nzt) )
    DO  k = 1, nzt
       vertical_influence(k) = MIN ( INT( l_grid(k) / &
                     ( wall_adjustment_factor * dzw(k) ) + 0.5_wp ), nzt - k )
    ENDDO

    DO  k = 1, nzt
       IF ( l_grid(k) > 1.5_wp * dx * wall_adjustment_factor .OR.  &
            l_grid(k) > 1.5_wp * dy * wall_adjustment_factor )  THEN
          WRITE( message_string, * ) 'grid anisotropy exceeds ', &
                                     'threshold given by only local', &
                                     ' &horizontal reduction of near_wall ', &
                                     'mixing length l_wall', &
                                     ' &starting from height level k = ', k, '.'
          CALL message( 'init_grid', 'PA0202', 0, 1, 0, 6, 0 )
          EXIT
       ENDIF
    ENDDO
    vertical_influence(0) = vertical_influence(1)

    DO  k = nzb + 1, nzb + vertical_influence(nzb)
       l_wall(k,:,:) = zu(k) - zw(nzb)
    ENDDO

!
!-- Set outer and inner index arrays for non-flat topography.
!-- Here consistency checks concerning domain size and periodicity are
!-- necessary.
!-- Within this SELECT CASE structure only nzb_local is initialized
!-- individually depending on the chosen topography type, all other index 
!-- arrays are initialized further below.
    SELECT CASE ( TRIM( topography ) )

       CASE ( 'flat' )
!
!--       nzb_local is required for the multigrid solver
          nzb_local = 0

       CASE ( 'single_building' )
!
!--       Single rectangular building, by default centered in the middle of the
!--       total domain
          blx = NINT( building_length_x / dx )
          bly = NINT( building_length_y / dy )
          bh  = MINLOC( ABS( zw - building_height ), 1 ) - 1
          IF ( ABS( zw(bh  ) - building_height ) == &
               ABS( zw(bh+1) - building_height )    )  bh = bh + 1

          IF ( building_wall_left == 9999999.9_wp )  THEN
             building_wall_left = ( nx + 1 - blx ) / 2 * dx
          ENDIF
          bxl = NINT( building_wall_left / dx )
          bxr = bxl + blx

          IF ( building_wall_south == 9999999.9_wp )  THEN
             building_wall_south = ( ny + 1 - bly ) / 2 * dy
          ENDIF
          bys = NINT( building_wall_south / dy )
          byn = bys + bly

!
!--       Building size has to meet some requirements
          IF ( ( bxl < 1 ) .OR. ( bxr > nx-1 ) .OR. ( bxr < bxl+3 ) .OR.  &
               ( bys < 1 ) .OR. ( byn > ny-1 ) .OR. ( byn < bys+3 ) )  THEN
             WRITE( message_string, * ) 'inconsistent building parameters:',   &
                                      '& bxl=', bxl, 'bxr=', bxr, 'bys=', bys, &
                                      'byn=', byn, 'nx=', nx, 'ny=', ny
             CALL message( 'init_grid', 'PA0203', 1, 2, 0, 6, 0 )
          ENDIF

!
!--       Define the building. 
          nzb_local = 0
          IF ( bxl <= nxr  .AND.  bxr >= nxl  .AND.                            &
               bys <= nyn  .AND.  byn >= nys )                                 &        
             nzb_local(MAX(nys,bys):MIN(nyn,byn),MAX(nxl,bxl):MIN(nxr,bxr)) = bh

          CALL exchange_horiz_2d_int( nzb_local, nys, nyn, nxl, nxr, nbgp )

       CASE ( 'single_street_canyon' )
!
!--       Single quasi-2D street canyon of infinite length in x or y direction.
!--       The canyon is centered in the other direction by default.
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
!
!--          Street canyon in y direction
             cwx = NINT( canyon_width_x / dx )
             IF ( canyon_wall_left == 9999999.9_wp )  THEN
                canyon_wall_left = ( nx + 1 - cwx ) / 2 * dx
             ENDIF
             cxl = NINT( canyon_wall_left / dx )
             cxr = cxl + cwx

          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
!
!--          Street canyon in x direction
             cwy = NINT( canyon_width_y / dy )
             IF ( canyon_wall_south == 9999999.9_wp )  THEN
                canyon_wall_south = ( ny + 1 - cwy ) / 2 * dy
             ENDIF
             cys = NINT( canyon_wall_south / dy )
             cyn = cys + cwy

          ELSE
             
             message_string = 'no street canyon width given'
             CALL message( 'init_grid', 'PA0204', 1, 2, 0, 6, 0 )
  
          ENDIF

          ch  = MINLOC( ABS( zw - canyon_height ), 1 ) - 1
          IF ( ABS( zw(ch  ) - canyon_height ) == &
               ABS( zw(ch+1) - canyon_height )    )  ch = ch + 1

          dp_level_ind_b = ch
!
!--       Street canyon size has to meet some requirements
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
             IF ( ( cxl < 1 ) .OR. ( cxr > nx-1 ) .OR. ( cwx < 3 ) .OR.        &
               ( ch < 3 ) )  THEN
                WRITE( message_string, * ) 'inconsistent canyon parameters:',  &
                                           '&cxl=', cxl, 'cxr=', cxr,          &
                                           'cwx=', cwx,                        &
                                           'ch=', ch, 'nx=', nx, 'ny=', ny
                CALL message( 'init_grid', 'PA0205', 1, 2, 0, 6, 0 ) 
             ENDIF
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
             IF ( ( cys < 1 ) .OR. ( cyn > ny-1 ) .OR. ( cwy < 3 ) .OR.        &
               ( ch < 3 ) )  THEN
                WRITE( message_string, * ) 'inconsistent canyon parameters:',  &
                                           '&cys=', cys, 'cyn=', cyn,          &
                                           'cwy=', cwy,                        &
                                           'ch=', ch, 'nx=', nx, 'ny=', ny
                CALL message( 'init_grid', 'PA0206', 1, 2, 0, 6, 0 ) 
             ENDIF
          ENDIF
          IF ( canyon_width_x /= 9999999.9_wp .AND.                            &                 
               canyon_width_y /= 9999999.9_wp )  THEN
             message_string = 'inconsistent canyon parameters:' //             &   
                              '&street canyon can only be oriented' //         &
                              '&either in x- or in y-direction'
             CALL message( 'init_grid', 'PA0207', 1, 2, 0, 6, 0 )
          ENDIF

          nzb_local = ch
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
             IF ( cxl <= nxr  .AND.  cxr >= nxl )                              &
                nzb_local(:,MAX(nxl,cxl+1):MIN(nxr,cxr-1)) = 0
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
             IF ( cys <= nyn  .AND.  cyn >= nys )                              &          
                nzb_local(MAX(nys,cys+1):MIN(nyn,cyn-1),:) = 0
          ENDIF

          CALL exchange_horiz_2d_int( nzb_local, nys, nyn, nxl, nxr, nbgp )

       CASE ( 'read_from_file' )

          ALLOCATE ( topo_height(nys:nyn,nxl:nxr) )

          DO  ii = 0, io_blocks-1
             IF ( ii == io_group )  THEN

!
!--             Arbitrary irregular topography data in PALM format (exactly
!--             matching the grid size and total domain size)
                OPEN( 90, FILE='TOPOGRAPHY_DATA'//TRIM( coupling_char ),       &
                          STATUS='OLD', FORM='FORMATTED', ERR=10 )
!
!--             Read topography PE-wise. Rows are read from nyn to nys, columns
!--             are read from nxl to nxr. At first, ny-nyn rows need to be skipped.
                skip_n_rows = 0
                DO WHILE ( skip_n_rows < ny - nyn )
                   READ( 90, * )  
                   skip_n_rows = skip_n_rows + 1
                ENDDO
!
!--             Read data from nyn to nys and nxl to nxr. Therefore, skip 
!--             column until nxl-1 is reached
                DO  j = nyn, nys, -1
                   READ( 90, *, ERR=11, END=11 )                               &
                                              ( dum, i = 0, nxl-1 ),           &
                                              ( topo_height(j,i), i = nxl, nxr )
                ENDDO

                GOTO 12
          
 10             message_string = 'file TOPOGRAPHY'//TRIM( coupling_char )//    &
                                 ' does not exist'
                CALL message( 'init_grid', 'PA0208', 1, 2, 0, 6, 0 )

 11             message_string = 'errors in file TOPOGRAPHY_DATA'//            &
                                 TRIM( coupling_char )
                CALL message( 'init_grid', 'PA0209', 1, 2, 0, 6, 0 )

 12             CLOSE( 90 )

             ENDIF
#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
          ENDDO

!
!--       Calculate the index height of the topography
          nzb_local = 0
          DO  i = nxl, nxr
             DO  j = nys, nyn
                nzb_local(j,i) = MINLOC( ABS( zw - topo_height(j,i) ), 1 ) - 1
                IF ( ABS( zw(nzb_local(j,i)  ) - topo_height(j,i) ) == &
                     ABS( zw(nzb_local(j,i)+1) - topo_height(j,i) )    )  &
                   nzb_local(j,i) = nzb_local(j,i) + 1
             ENDDO
          ENDDO

          DEALLOCATE ( topo_height )
!
!--       Filter topography, i.e. fill holes resolved by only one grid point.  
!--       Such holes are suspected to lead to velocity blow-ups as continuity 
!--       equation on discrete grid cannot be fulfilled in such case.
!--       For now, check only for holes and fill them to the lowest height level
!--       of the directly adjoining grid points along x- and y- direction. 
!--       Before checking for holes, set lateral boundary conditions for 
!--       topography. After hole-filling, boundary conditions must be set again!
          CALL exchange_horiz_2d_int( nzb_local, nys, nyn, nxl, nxr, nbgp )
          
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  nzb_local(-1,:)   = nzb_local(0,:)
             IF ( nyn == ny )  nzb_local(ny+1,:) = nzb_local(ny,:)
          ENDIF

          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  nzb_local(:,-1)   = nzb_local(:,0)
             IF ( nxr == nx )  nzb_local(:,nx+1) = nzb_local(:,nx)          
          ENDIF

          num_hole_l = 0
          DO i = nxl, nxr
             DO j = nys, nyn

                num_wall = 0

                IF ( nzb_local(j-1,i) > nzb_local(j,i) )                       &
                   num_wall = num_wall + 1
                IF ( nzb_local(j+1,i) > nzb_local(j,i) )                       &
                   num_wall = num_wall + 1
                IF ( nzb_local(j,i-1) > nzb_local(j,i) )                       &
                   num_wall = num_wall + 1
                IF ( nzb_local(j,i+1) > nzb_local(j,i) )                       &
                   num_wall = num_wall + 1

                IF ( num_wall == 4 )  THEN
                   nzb_local(j,i) = MIN( nzb_local(j-1,i), nzb_local(j+1,i),   &
                                         nzb_local(j,i-1), nzb_local(j,i+1) )
                   num_hole_l     = num_hole_l + 1
                ENDIF
             ENDDO
          ENDDO
!
!--       Count the total number of holes, required for informative message.
#if defined( __parallel )
          CALL MPI_ALLREDUCE( num_hole_l, num_hole, 1, MPI_INTEGER, MPI_SUM,   &
                              comm2d, ierr )
#else
          num_hole = num_hole_l
#endif    
!
!--       Create an informative message if any hole was removed.
          IF ( num_hole > 0 )  THEN
             WRITE( message_string, * ) num_hole, 'hole(s) resolved by only '//&
                                                  'one grid point were filled'
             CALL message( 'init_grid', 'PA0430', 0, 0, 0, 6, 0 )
          ENDIF
!
!--       Exchange ghost-points, as well as add cyclic or Neumann boundary 
!--       conditions.
          CALL exchange_horiz_2d_int( nzb_local, nys, nyn, nxl, nxr, nbgp )
          
          IF ( .NOT. bc_ns_cyc )  THEN
             IF ( nys == 0  )  nzb_local(-1,:)   = nzb_local(0,:)
             IF ( nyn == ny )  nzb_local(ny+1,:) = nzb_local(ny,:)
          ENDIF

          IF ( .NOT. bc_lr_cyc )  THEN
             IF ( nxl == 0  )  nzb_local(:,-1)   = nzb_local(:,0)
             IF ( nxr == nx )  nzb_local(:,nx+1) = nzb_local(:,nx)          
          ENDIF

       CASE DEFAULT
!
!--       The DEFAULT case is reached either if the parameter topography
!--       contains a wrong character string or if the user has defined a special
!--       case in the user interface. There, the subroutine user_init_grid 
!--       checks which of these two conditions applies.
          CALL user_init_grid( nzb_local )

    END SELECT
!
!-- Determine the maximum level of topography. Furthermore it is used for
!-- steering the degradation of order of the applied advection scheme.
!-- In case of non-cyclic lateral boundaries, the order of the advection
!-- scheme has to be reduced up to nzt (required at the lateral boundaries).
#if defined( __parallel )
    CALL MPI_ALLREDUCE( MAXVAL( nzb_local ) + 1, nzb_max, 1, MPI_INTEGER,      &
                        MPI_MAX, comm2d, ierr )
#else
    nzb_max = MAXVAL( nzb_local ) + 1
#endif
    IF ( inflow_l .OR. outflow_l .OR. inflow_r .OR. outflow_r .OR.             &
         inflow_n .OR. outflow_n .OR. inflow_s .OR. outflow_s .OR.             &
         nest_domain )                                                         &
    THEN
       nzb_max = nzt
    ENDIF

!
!-- Consistency checks and index array initialization are only required for
!-- non-flat topography, also the initialization of topography height arrays
!-- zu_s_inner and zw_w_inner
    IF ( TRIM( topography ) /= 'flat' )  THEN
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MAXVAL( nzb_local ), nzb_local_max, 1, MPI_INTEGER, &
                           MPI_MAX, comm2d, ierr )
       CALL MPI_ALLREDUCE( MINVAL( nzb_local ), nzb_local_min, 1, MPI_INTEGER, &
                           MPI_MIN, comm2d, ierr )                           
#else
       nzb_local_max = MAXVAL( nzb_local )
       nzb_local_min = MINVAL( nzb_local )
#endif

!
!--    Consistency checks
       IF ( nzb_local_min < 0  .OR.  nzb_local_max  > nz + 1 )  THEN
          WRITE( message_string, * ) 'nzb_local values are outside the',       &
                                'model domain',                                &
                                '&MINVAL( nzb_local ) = ', nzb_local_min,      &
                                '&MAXVAL( nzb_local ) = ', nzb_local_max
          CALL message( 'init_grid', 'PA0210', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( topography_grid_convention == 'cell_edge' )  THEN
! 
!--       The array nzb_local as defined using the 'cell_edge' convention 
!--       describes the actual total size of topography which is defined at the 
!--       cell edges where u=0 on the topography walls in x-direction and v=0 
!--       on the topography walls in y-direction. However, PALM uses individual
!--       arrays nzb_u|v|w|s_inner|outer that are based on nzb_s_inner.
!--       Therefore, the extent of topography in nzb_local is now reduced by 
!--       1dx at the E topography walls and by 1dy at the N topography walls 
!--       to form the basis for nzb_s_inner. 
!--       Note, the reverse memory access (i-j instead of j-i) is absolutely
!--       required at this point.
          DO  j = nys+1, nyn+1
             DO  i = nxl-1, nxr
                nzb_local(j,i) = MIN( nzb_local(j,i), nzb_local(j,i+1) )
             ENDDO
          ENDDO
!
!--       Exchange ghost points
          CALL exchange_horiz_2d_int( nzb_local, nys, nyn, nxl, nxr, nbgp )

          DO  i = nxl, nxr+1
             DO  j = nys-1, nyn
                nzb_local(j,i) = MIN( nzb_local(j,i), nzb_local(j+1,i) )
             ENDDO
          ENDDO
!
!--       Exchange ghost points          
          CALL exchange_horiz_2d_int( nzb_local, nys, nyn, nxl, nxr, nbgp )
       ENDIF
!
!--    Initialize index arrays nzb_s_inner and nzb_w_inner
       nzb_s_inner = nzb_local
       nzb_w_inner = nzb_local

!
!--    Initialize remaining index arrays:
!--    first pre-initialize them with nzb_s_inner...
       nzb_u_inner = nzb_s_inner
       nzb_u_outer = nzb_s_inner
       nzb_v_inner = nzb_s_inner
       nzb_v_outer = nzb_s_inner
       nzb_w_outer = nzb_s_inner
       nzb_s_outer = nzb_s_inner

!
!--    ...then extend pre-initialized arrays in their according directions
!--    based on nzb_local using nzb_tmp as a temporary global index array

!
!--    nzb_s_outer: 
!--    extend nzb_local east-/westwards first, then north-/southwards
       nzb_tmp = nzb_local
       DO  j = nys, nyn
          DO  i = nxl, nxr
             nzb_tmp(j,i) = MAX( nzb_local(j,i-1), nzb_local(j,i),             &
                                 nzb_local(j,i+1) )
          ENDDO
       ENDDO
       
       CALL exchange_horiz_2d_int( nzb_tmp, nys, nyn, nxl, nxr, nbgp )
       
       DO  i = nxl, nxr
          DO  j = nys, nyn
             nzb_s_outer(j,i) = MAX( nzb_tmp(j-1,i), nzb_tmp(j,i),             &
                                     nzb_tmp(j+1,i) )
          ENDDO
!
!--       non-cyclic boundary conditions (overwritten by call of
!--       exchange_horiz_2d_int below in case of cyclic boundary conditions)
          IF ( nys == 0 )  THEN
             j = -1
             nzb_s_outer(j,i) = MAX( nzb_tmp(j+1,i), nzb_tmp(j,i) )
          ENDIF
          IF ( nyn == ny )  THEN
             j = ny + 1
             nzb_s_outer(j,i) = MAX( nzb_tmp(j-1,i), nzb_tmp(j,i) )
          ENDIF
       ENDDO
!
!--    nzb_w_outer: 
!--    identical to nzb_s_outer
       nzb_w_outer = nzb_s_outer

!
!--    nzb_u_inner: 
!--    extend nzb_local rightwards only
       nzb_tmp = nzb_local
       DO  j = nys, nyn
          DO  i = nxl, nxr
             nzb_tmp(j,i) = MAX( nzb_local(j,i-1), nzb_local(j,i) )
          ENDDO
       ENDDO
       
       CALL exchange_horiz_2d_int( nzb_tmp, nys, nyn, nxl, nxr, nbgp )
       
       nzb_u_inner = nzb_tmp
!
!--    nzb_u_outer: 
!--    extend current nzb_tmp (nzb_u_inner) north-/southwards
       DO  i = nxl, nxr
          DO  j = nys, nyn
             nzb_u_outer(j,i) = MAX( nzb_tmp(j-1,i), nzb_tmp(j,i),             &
                                     nzb_tmp(j+1,i) )
          ENDDO
!
!--       non-cyclic boundary conditions (overwritten by call of
!--       exchange_horiz_2d_int below in case of cyclic boundary conditions)
          IF ( nys == 0 )  THEN
             j = -1
             nzb_u_outer(j,i) = MAX( nzb_tmp(j+1,i), nzb_tmp(j,i) )
          ENDIF
          IF ( nyn == ny )  THEN
             j = ny + 1
             nzb_u_outer(j,i) = MAX( nzb_tmp(j-1,i), nzb_tmp(j,i) )
          ENDIF
       ENDDO

!
!--    nzb_v_inner:
!--    extend nzb_local northwards only
       nzb_tmp = nzb_local
       DO  i = nxl, nxr
          DO  j = nys, nyn
             nzb_tmp(j,i) = MAX( nzb_local(j-1,i), nzb_local(j,i) )
          ENDDO
       ENDDO
       
       CALL exchange_horiz_2d_int( nzb_tmp, nys, nyn, nxl, nxr, nbgp )      
       nzb_v_inner = nzb_tmp

!
!--    nzb_v_outer: 
!--    extend current nzb_tmp (nzb_v_inner) right-/leftwards
       DO  j = nys, nyn
          DO  i = nxl, nxr
             nzb_v_outer(j,i) = MAX( nzb_tmp(j,i-1), nzb_tmp(j,i),             &
                                     nzb_tmp(j,i+1) )
          ENDDO
!
!--       non-cyclic boundary conditions (overwritten by call of
!--       exchange_horiz_2d_int below in case of cyclic boundary conditions)
          IF ( nxl == 0 )  THEN
             i = -1
             nzb_v_outer(j,i) = MAX( nzb_tmp(j,i+1), nzb_tmp(j,i) )
          ENDIF
          IF ( nxr == nx )  THEN
             i = nx + 1
             nzb_v_outer(j,i) = MAX( nzb_tmp(j,i-1), nzb_tmp(j,i) )
          ENDIF
       ENDDO

!
!--    Exchange of lateral boundary values (parallel computers) and cyclic
!--    boundary conditions, if applicable.
!--    Since nzb_s_inner and nzb_w_inner are derived directly from nzb_local
!--    they do not require exchange and are not included here.
       CALL exchange_horiz_2d_int( nzb_u_inner, nys, nyn, nxl, nxr, nbgp )
       CALL exchange_horiz_2d_int( nzb_u_outer, nys, nyn, nxl, nxr, nbgp )
       CALL exchange_horiz_2d_int( nzb_v_inner, nys, nyn, nxl, nxr, nbgp )
       CALL exchange_horiz_2d_int( nzb_v_outer, nys, nyn, nxl, nxr, nbgp )
       CALL exchange_horiz_2d_int( nzb_w_outer, nys, nyn, nxl, nxr, nbgp )
       CALL exchange_horiz_2d_int( nzb_s_outer, nys, nyn, nxl, nxr, nbgp )

!
!--    Allocate and set the arrays containing the topography height
       ALLOCATE( zu_s_inner(0:nx+1,0:ny+1), zw_w_inner(0:nx+1,0:ny+1),         &
                 zu_s_inner_l(0:nx+1,0:ny+1), zw_w_inner_l(0:nx+1,0:ny+1) )
                 
       zu_s_inner   = 0.0_wp
       zw_w_inner   = 0.0_wp
       zu_s_inner_l = 0.0_wp
       zw_w_inner_l = 0.0_wp
       
       DO  i = nxl, nxr
          DO  j = nys, nyn
             zu_s_inner_l(i,j) = zu(nzb_local(j,i))
             zw_w_inner_l(i,j) = zw(nzb_local(j,i))
          ENDDO
       ENDDO
       
#if defined( __parallel )
       CALL MPI_REDUCE( zu_s_inner_l, zu_s_inner, (nx+2)*(ny+2),         &
                           MPI_REAL, MPI_SUM, 0, comm2d, ierr )       
       CALL MPI_REDUCE( zw_w_inner_l, zw_w_inner, (nx+2)*(ny+2),         &
                           MPI_REAL, MPI_SUM, 0, comm2d, ierr )  
#else
       zu_s_inner = zu_s_inner_l
       zw_w_inner = zw_w_inner_l
#endif

      DEALLOCATE( zu_s_inner_l, zw_w_inner_l )
      IF ( myid /= 0 )  DEALLOCATE( zu_s_inner, zw_w_inner )
!
!--   Set south and left ghost points, required for netcdf output
      IF ( myid == 0 )  THEN
         IF( bc_lr_cyc )  THEN
            zu_s_inner(nx+1,:) = zu_s_inner(0,:)
            zw_w_inner(nx+1,:) = zw_w_inner(0,:)
         ELSE
            zu_s_inner(nx+1,:) = zu_s_inner(nx,:)
            zw_w_inner(nx+1,:) = zw_w_inner(nx,:)
         ENDIF
         IF( bc_ns_cyc )  THEN
            zu_s_inner(:,ny+1) = zu_s_inner(:,0)
            zw_w_inner(:,ny+1) = zw_w_inner(:,0)
         ELSE
            zu_s_inner(:,ny+1) = zu_s_inner(:,ny)
            zw_w_inner(:,ny+1) = zw_w_inner(:,ny)
         ENDIF
      ENDIF
!
!--    Set flag arrays to be used for masking of grid points
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb, nzt+1
                IF ( k <= nzb_s_inner(j,i) )  rflags_s_inner(k,j,i) = 0.0_wp
                IF ( k <= nzb_s_inner(j,i) )  rflags_invers(j,i,k)  = 0.0_wp
             ENDDO
          ENDDO
       ENDDO

    ENDIF
!
!-- Deallocate temporary array, as it might be reused for different
!-- grid-levels further below.
    DEALLOCATE( nzb_tmp )

!
!-- Set the individual index arrays which define the k index from which on
!-- the usual finite difference form (which does not use surface fluxes) is
!-- applied 
    IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
       nzb_diff_u         = nzb_u_inner + 2
       nzb_diff_v         = nzb_v_inner + 2
       nzb_diff_s_inner   = nzb_s_inner + 2
       nzb_diff_s_outer   = nzb_s_outer + 2
    ELSE
       nzb_diff_u         = nzb_u_inner + 1
       nzb_diff_v         = nzb_v_inner + 1
       nzb_diff_s_inner   = nzb_s_inner + 1
       nzb_diff_s_outer   = nzb_s_outer + 1
    ENDIF

!
!-- Calculation of wall switches and factors required by diffusion_u/v.f90 and 
!-- for limitation of near-wall mixing length l_wall further below
    corner_nl = 0
    corner_nr = 0
    corner_sl = 0
    corner_sr = 0
    wall_l    = 0
    wall_n    = 0
    wall_r    = 0
    wall_s    = 0

    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       u-component
          IF ( nzb_u_outer(j,i) > nzb_u_outer(j+1,i) )  THEN
             wall_u(j,i) = 1.0_wp   ! north wall (location of adjacent fluid)
             fym(j,i)    = 0.0_wp
             fyp(j,i)    = 1.0_wp
          ELSEIF ( nzb_u_outer(j,i) > nzb_u_outer(j-1,i) )  THEN
             wall_u(j,i) = 1.0_wp   ! south wall (location of adjacent fluid)
             fym(j,i)    = 1.0_wp
             fyp(j,i)    = 0.0_wp
          ENDIF
!
!--       v-component
          IF ( nzb_v_outer(j,i) > nzb_v_outer(j,i+1) )  THEN
             wall_v(j,i) = 1.0_wp   ! rigth wall (location of adjacent fluid)
             fxm(j,i)    = 0.0_wp
             fxp(j,i)    = 1.0_wp
          ELSEIF ( nzb_v_outer(j,i) > nzb_v_outer(j,i-1) )  THEN
             wall_v(j,i) = 1.0_wp   ! left wall (location of adjacent fluid)
             fxm(j,i)    = 1.0_wp
             fxp(j,i)    = 0.0_wp
          ENDIF
!
!--       w-component, also used for scalars, separate arrays for shear 
!--       production of tke
          IF ( nzb_w_outer(j,i) > nzb_w_outer(j+1,i) )  THEN
             wall_e_y(j,i) =  1.0_wp   ! north wall (location of adjacent fluid)
             wall_w_y(j,i) =  1.0_wp
             fwym(j,i)     =  0.0_wp
             fwyp(j,i)     =  1.0_wp
          ELSEIF ( nzb_w_outer(j,i) > nzb_w_outer(j-1,i) )  THEN
             wall_e_y(j,i) = -1.0_wp   ! south wall (location of adjacent fluid)
             wall_w_y(j,i) =  1.0_wp
             fwym(j,i)     =  1.0_wp
             fwyp(j,i)     =  0.0_wp
          ENDIF
          IF ( nzb_w_outer(j,i) > nzb_w_outer(j,i+1) )  THEN
             wall_e_x(j,i) =  1.0_wp   ! right wall (location of adjacent fluid)
             wall_w_x(j,i) =  1.0_wp
             fwxm(j,i)     =  0.0_wp
             fwxp(j,i)     =  1.0_wp
          ELSEIF ( nzb_w_outer(j,i) > nzb_w_outer(j,i-1) )  THEN
             wall_e_x(j,i) = -1.0_wp   ! left wall (location of adjacent fluid)
             wall_w_x(j,i) =  1.0_wp
             fwxm(j,i)     =  1.0_wp
             fwxp(j,i)     =  0.0_wp
          ENDIF
!
!--       Wall and corner locations inside buildings for limitation of
!--       near-wall mixing length l_wall
          IF ( nzb_s_inner(j,i) > nzb_s_inner(j+1,i) )  THEN

             wall_n(j,i) = nzb_s_inner(j+1,i) + 1            ! North wall

             IF ( nzb_s_inner(j,i) > nzb_s_inner(j,i-1) )  THEN
                corner_nl(j,i) = MAX( nzb_s_inner(j+1,i),  & ! Northleft corner
                                      nzb_s_inner(j,i-1) ) + 1
             ENDIF

             IF ( nzb_s_inner(j,i) > nzb_s_inner(j,i+1) )  THEN
                corner_nr(j,i) = MAX( nzb_s_inner(j+1,i),  & ! Northright corner
                                      nzb_s_inner(j,i+1) ) + 1
             ENDIF

          ENDIF

          IF ( nzb_s_inner(j,i) > nzb_s_inner(j-1,i) )  THEN

             wall_s(j,i) = nzb_s_inner(j-1,i) + 1            ! South wall
             IF ( nzb_s_inner(j,i) > nzb_s_inner(j,i-1) )  THEN
                corner_sl(j,i) = MAX( nzb_s_inner(j-1,i),  & ! Southleft corner
                                      nzb_s_inner(j,i-1) ) + 1
             ENDIF

             IF ( nzb_s_inner(j,i) > nzb_s_inner(j,i+1) )  THEN
                corner_sr(j,i) = MAX( nzb_s_inner(j-1,i),  & ! Southright corner
                                      nzb_s_inner(j,i+1) ) + 1
             ENDIF

          ENDIF

          IF ( nzb_s_inner(j,i) > nzb_s_inner(j,i-1) )  THEN
             wall_l(j,i) = nzb_s_inner(j,i-1) + 1            ! Left wall
          ENDIF

          IF ( nzb_s_inner(j,i) > nzb_s_inner(j,i+1) )  THEN
             wall_r(j,i) = nzb_s_inner(j,i+1) + 1            ! Right wall
          ENDIF

       ENDDO
    ENDDO
!
!-- Calculate wall flag arrays for the multigrid method.
!-- Please note, wall flags are only applied in the not cache-optimized 
!-- version.
    IF ( psolver == 'multigrid_noopt' )  THEN

!
!--    Gridpoint increment of the current level. 
       inc = 1
       DO  l = maximum_grid_level, 1 , -1
!
!--       Set grid_level as it is required for exchange_horiz_2d_int
          grid_level = l

          nxl_l = nxl_mg(l)
          nxr_l = nxr_mg(l)
          nys_l = nys_mg(l)
          nyn_l = nyn_mg(l)
          nzt_l = nzt_mg(l)
!
!--       Assign the flag level to be calculated
          SELECT CASE ( l )
             CASE ( 1 )
                flags => wall_flags_1
             CASE ( 2 )
                flags => wall_flags_2
             CASE ( 3 )
                flags => wall_flags_3
             CASE ( 4 )
                flags => wall_flags_4
             CASE ( 5 )
                flags => wall_flags_5
             CASE ( 6 )
                flags => wall_flags_6
             CASE ( 7 )
                flags => wall_flags_7
             CASE ( 8 )
                flags => wall_flags_8
             CASE ( 9 )
                flags => wall_flags_9
             CASE ( 10 )
                flags => wall_flags_10
          END SELECT

!
!--       Depending on the grid level, set the respective bits in case of
!--       neighbouring walls
!--       Bit 0:  wall to the bottom
!--       Bit 1:  wall to the top (not realized in remaining PALM code so far)
!--       Bit 2:  wall to the south
!--       Bit 3:  wall to the north
!--       Bit 4:  wall to the left
!--       Bit 5:  wall to the right
!--       Bit 6:  inside building

          flags = 0

!
!--       In case of masking method, flags are not set and multigrid method
!--       works like FFT-solver
          IF ( .NOT. masking_method )  THEN

!
!--          Allocate temporary array for topography heights on coarser grid
!--          level. Please note, 2 ghoist points are required, in order to
!--          calculate flags() on the interior ghost point. 
             ALLOCATE( nzb_tmp(nys_l-2:nyn_l+2,nxl_l-2:nxr_l+2) )
             nzb_tmp = 0
             
             DO  i = nxl_l, nxr_l
                DO  j = nys_l, nyn_l
                   nzb_tmp(j,i) = nzb_local(j*inc,i*inc)
                ENDDO
             ENDDO
!
!--          Exchange ghost points on respective multigrid level. 2 ghost points
!--          are required, in order to calculate flags on 
!--          nys_l-1 / nyn_l+1 / nxl_l-1 / nxr_l+1. The alternative would be to 
!--          exchange 3D-INTEGER array flags on the respective multigrid level.
             CALL exchange_horiz_2d_int( nzb_tmp, nys_l, nyn_l, nxl_l, nxr_l, 2 )
!
!--          Set non-cyclic boundary conditions on respective multigrid level
             IF ( .NOT. bc_ns_cyc )  THEN
                IF ( inflow_s .OR. outflow_s .OR. nest_bound_s  )  THEN
                   nzb_tmp(-2,:) = nzb_tmp(0,:)
                   nzb_tmp(-1,:) = nzb_tmp(0,:)
                ENDIF
                IF ( inflow_n .OR. outflow_n .OR. nest_bound_n )  THEN
                   nzb_tmp(nyn_l+2,:) = nzb_tmp(nyn_l,:)
                   nzb_tmp(nyn_l+1,:) = nzb_tmp(nyn_l,:)
                ENDIF
             ENDIF
             IF ( .NOT. bc_lr_cyc )  THEN
                IF ( inflow_l .OR. outflow_l .OR. nest_bound_l  )  THEN
                   nzb_tmp(:,-2) = nzb_tmp(:,0)
                   nzb_tmp(:,-1) = nzb_tmp(:,0)
                ENDIF
                IF ( inflow_r .OR. outflow_r .OR. nest_bound_r )  THEN
                   nzb_tmp(:,nxr_l+1) = nzb_tmp(:,nxr_l)    
                   nzb_tmp(:,nxr_l+2) = nzb_tmp(:,nxr_l)      
                ENDIF        
             ENDIF
                       
             DO  i = nxl_l-1, nxr_l+1
                DO  j = nys_l-1, nyn_l+1
                   DO  k = nzb, nzt_l+1     
!
!--                   Inside/outside building (inside building does not need
!--                   further tests for walls)
                      IF ( k*inc <= nzb_tmp(j,i) )  THEN

                         flags(k,j,i) = IBSET( flags(k,j,i), 6 )

                      ELSE
!
!--                      Bottom wall
                         IF ( (k-1)*inc <= nzb_tmp(j,i) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 0 )
                         ENDIF
!
!--                      South wall
                         IF ( k*inc <= nzb_tmp(j-1,i) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 2 )
                         ENDIF
!
!--                      North wall
                         IF ( k*inc <= nzb_tmp(j+1,i) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 3 )
                         ENDIF
!
!--                      Left wall
                         IF ( k*inc <= nzb_tmp(j,i-1) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 4 )
                         ENDIF
!
!--                      Right wall
                         IF ( k*inc <= nzb_tmp(j,i+1) )  THEN
                            flags(k,j,i) = IBSET( flags(k,j,i), 5 )
                         ENDIF

                      ENDIF
                            
                   ENDDO
                ENDDO
             ENDDO

             DEALLOCATE( nzb_tmp )

          ENDIF

          inc = inc * 2

       ENDDO
!
!--    Reset grid_level to "normal" grid
       grid_level = 0
       
    ENDIF
!
!-- Allocate flags needed for masking walls. Even though these flags are only
!-- required in the ws-scheme, the arrays need to be allocated here as they are 
!-- used in OpenACC directives. 
    ALLOCATE( wall_flags_0(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                     &
              wall_flags_00(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    wall_flags_0  = 0
    wall_flags_00 = 0
!
!-- Init flags for ws-scheme to degrade order of the numerics near walls, i.e. 
!-- to decrease the numerical stencil appropriately.
    IF ( momentum_advec == 'ws-scheme'  .OR.  scalar_advec == 'ws-scheme'  .OR.&
         scalar_advec   == 'ws-scheme-mono' )  THEN
       CALL ws_init_flags
    ENDIF

!
!-- In case of topography: limit near-wall mixing length l_wall further:
!-- Go through all points of the subdomain one by one and look for the closest
!-- surface
    IF ( TRIM(topography) /= 'flat' )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn

             nzb_si = nzb_s_inner(j,i)
             vi     = vertical_influence(nzb_si)

             IF ( wall_n(j,i) > 0 )  THEN
!
!--             North wall (y distance)
                DO  k = wall_n(j,i), nzb_si
                   l_wall(k,j+1,i) = MIN( l_wall(k,j+1,i), 0.5_wp * dy )
                ENDDO
!
!--             Above North wall (yz distance)
                DO  k = nzb_si + 1, nzb_si + vi
                   l_wall(k,j+1,i) = MIN( l_wall(k,j+1,i),                     &
                                          SQRT( 0.25_wp * dy**2 +              &
                                          ( zu(k) - zw(nzb_si) )**2 ) )
                ENDDO
!
!--             Northleft corner (xy distance)
                IF ( corner_nl(j,i) > 0 )  THEN
                   DO  k = corner_nl(j,i), nzb_si
                      l_wall(k,j+1,i-1) = MIN( l_wall(k,j+1,i-1), &
                                               0.5_wp * SQRT( dx**2 + dy**2 ) )
                   ENDDO
!
!--                Above Northleft corner (xyz distance)
                   DO  k = nzb_si + 1, nzb_si + vi
                      l_wall(k,j+1,i-1) = MIN( l_wall(k,j+1,i-1),              &
                                            SQRT( 0.25_wp * (dx**2 + dy**2) +  &
                                            ( zu(k) - zw(nzb_si) )**2 ) )
                   ENDDO
                ENDIF
!
!--             Northright corner (xy distance)
                IF ( corner_nr(j,i) > 0 )  THEN
                   DO  k = corner_nr(j,i), nzb_si
                       l_wall(k,j+1,i+1) = MIN( l_wall(k,j+1,i+1),             &
                                                0.5_wp * SQRT( dx**2 + dy**2 ) )
                   ENDDO
!
!--                Above northright corner (xyz distance)
                   DO  k = nzb_si + 1, nzb_si + vi
                      l_wall(k,j+1,i+1) = MIN( l_wall(k,j+1,i+1),              &
                                            SQRT( 0.25_wp * (dx**2 + dy**2) +  &
                                            ( zu(k) - zw(nzb_si) )**2 ) )
                   ENDDO
                ENDIF
             ENDIF

             IF ( wall_s(j,i) > 0 )  THEN
!
!--             South wall (y distance)
                DO  k = wall_s(j,i), nzb_si
                   l_wall(k,j-1,i) = MIN( l_wall(k,j-1,i), 0.5_wp * dy )
                ENDDO
!
!--             Above south wall (yz distance)
                DO  k = nzb_si + 1, nzb_si + vi
                   l_wall(k,j-1,i) = MIN( l_wall(k,j-1,i),                     &
                                          SQRT( 0.25_wp * dy**2 +              &
                                          ( zu(k) - zw(nzb_si) )**2 ) )
                ENDDO
!
!--             Southleft corner (xy distance)
                IF ( corner_sl(j,i) > 0 )  THEN
                   DO  k = corner_sl(j,i), nzb_si
                      l_wall(k,j-1,i-1) = MIN( l_wall(k,j-1,i-1),              &
                                               0.5_wp * SQRT( dx**2 + dy**2 ) )
                   ENDDO
!
!--                Above southleft corner (xyz distance)
                   DO  k = nzb_si + 1, nzb_si + vi
                      l_wall(k,j-1,i-1) = MIN( l_wall(k,j-1,i-1),              &
                                            SQRT( 0.25_wp * (dx**2 + dy**2) +  &
                                            ( zu(k) - zw(nzb_si) )**2 ) )
                   ENDDO
                ENDIF
!
!--             Southright corner (xy distance)
                IF ( corner_sr(j,i) > 0 )  THEN
                   DO  k = corner_sr(j,i), nzb_si
                      l_wall(k,j-1,i+1) = MIN( l_wall(k,j-1,i+1),              &
                                               0.5_wp * SQRT( dx**2 + dy**2 ) )
                   ENDDO
!
!--                Above southright corner (xyz distance)
                   DO  k = nzb_si + 1, nzb_si + vi
                      l_wall(k,j-1,i+1) = MIN( l_wall(k,j-1,i+1),              &
                                            SQRT( 0.25_wp * (dx**2 + dy**2) +  &
                                            ( zu(k) - zw(nzb_si) )**2 ) )
                   ENDDO
                ENDIF

             ENDIF

             IF ( wall_l(j,i) > 0 )  THEN
!
!--             Left wall (x distance)
                DO  k = wall_l(j,i), nzb_si
                   l_wall(k,j,i-1) = MIN( l_wall(k,j,i-1), 0.5_wp * dx )
                ENDDO
!
!--             Above left wall (xz distance)
                DO  k = nzb_si + 1, nzb_si + vi
                   l_wall(k,j,i-1) = MIN( l_wall(k,j,i-1),                     &
                                       SQRT( 0.25_wp * dx**2 +                 &
                                       ( zu(k) - zw(nzb_si) )**2 ) )
                ENDDO
             ENDIF

             IF ( wall_r(j,i) > 0 )  THEN
!
!--             Right wall (x distance)
                DO  k = wall_r(j,i), nzb_si
                   l_wall(k,j,i+1) = MIN( l_wall(k,j,i+1), 0.5_wp * dx )
                ENDDO
!
!--             Above right wall (xz distance)
                DO  k = nzb_si + 1, nzb_si + vi
                   l_wall(k,j,i+1) = MIN( l_wall(k,j,i+1),                     &
                                          SQRT( 0.25_wp * dx**2 +              &
                                          ( zu(k) - zw(nzb_si) )**2 ) )
                ENDDO

             ENDIF

          ENDDO
       ENDDO

    ENDIF

!
!-- Multiplication with wall_adjustment_factor
    l_wall = wall_adjustment_factor * l_wall

!
!-- Set lateral boundary conditions for l_wall
    CALL exchange_horiz( l_wall, nbgp )

    DEALLOCATE( corner_nl, corner_nr, corner_sl, corner_sr, nzb_local, &
                vertical_influence, wall_l, wall_n, wall_r, wall_s )


 END SUBROUTINE init_grid
