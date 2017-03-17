!> @file lpm_init.f90
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
! $Id: lpm_init.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 2016-06-09 16:25:25Z suehring
! Bugfix in determining initial particle height and grid index in case of 
! seed_follows_topography.
! Bugfix concerning random positions, ensure that particles do not move more 
! than one grid length.
! Bugfix logarithmic interpolation.
! Initial setting of sgs_wf_part.
!
! 1890 2016-04-22 08:52:11Z hoffmann
! Initialization of aerosol equilibrium radius not possible in supersaturated 
! environments. Therefore, a maximum supersaturation of -1 % is assumed during
! initialization.
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod

! 
! 1871 2016-04-15 11:46:09Z hoffmann
! Initialization of aerosols added.
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 1831 2016-04-07 13:15:51Z hoffmann
! curvature_solution_effects moved to particle_attributes
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module added
!
! 1725 2015-11-17 13:01:51Z hoffmann 
! Bugfix: Processor-dependent seed for random function is generated before it is 
! used.
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer.
!
! 1685 2015-10-08 07:32:13Z raasch
! bugfix concerning vertical index offset in case of ocean
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1575 2015-03-27 09:56:27Z raasch
! initial vertical particle position is allowed to follow the topography
!
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
! Kind definition added to all floating point numbers.
! lpm_init changed form a subroutine to a module.
! 
! 1327 2014-03-21 11:00:16Z raasch
! -netcdf_output
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
! bugfix: #if defined( __parallel ) added
!
! 1314 2014-03-14 18:25:17Z suehring
! Vertical logarithmic interpolation of horizontal particle speed for particles 
! between roughness height and first vertical grid level.
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! routine renamed: init_particles -> lpm_init
! de_dx, de_dy, de_dz are allocated here (instead of automatic arrays in
! advec_particles),
! sort_particles renamed lpm_sort_arrays, user_init_particles renamed lpm_init
!
! 828 2012-02-21 12:00:36Z raasch
! call of init_kernels, particle feature color renamed class
!
! 824 2012-02-17 09:09:57Z raasch
! particle attributes speed_x|y|z_sgs renamed rvar1|2|3,
! array particles implemented as pointer
!
! 667 2010-12-23 12:06:00Z suehring/gryschka
! nxl-1, nxr+1, nys-1, nyn+1 replaced by nxlg, nxrg, nysg, nyng for allocation 
! of arrays.
!
! Revision 1.1  1999/11/25 16:22:38  raasch
! Initial revision
!
!
! Description:
! ------------
!> This routine initializes a set of particles and their attributes (position,
!> radius, ..) which are used by the Lagrangian particle model (see lpm).
!------------------------------------------------------------------------------!
 MODULE lpm_init_mod
 

    USE arrays_3d,                                                             &
        ONLY:  de_dx, de_dy, de_dz, zu, zw, z0

    USE control_parameters,                                                    &
        ONLY:  cloud_droplets, constant_flux_layer, current_timestep_number,   &
               dz, initializing_actions, message_string, ocean, simulated_time

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy

    USE indices,                                                               &
        ONLY:  nx, nxl, nxlg, nxrg, nxr, ny, nyn, nys, nyng, nysg, nz, nzb,    &
               nzb_w_inner, nzt

    USE kinds

    USE lpm_collision_kernels_mod,                                             &
        ONLY:  init_kernels

    USE netcdf_interface,                                                      &
        ONLY:  netcdf_data_format

    USE particle_attributes,                                                   &
        ONLY:   alloc_factor, bc_par_b, bc_par_lr, bc_par_ns, bc_par_t,        &
                block_offset, block_offset_def, collision_kernel,              &
                curvature_solution_effects,                                    &
                density_ratio, grid_particles,                                 &
                initial_weighting_factor, ibc_par_b, ibc_par_lr, ibc_par_ns,   &
                ibc_par_t, iran_part, log_z_z0,                                &
                max_number_of_particle_groups, maximum_number_of_particles,    &
                min_nr_particle, mpi_particle_type,                            &
                number_of_particles,                                           &
                number_of_particle_groups, number_of_sublayers,                &
                offset_ocean_nzt, offset_ocean_nzt_m1,                         &
                particles, particle_advection_start, particle_groups,          &
                particle_groups_type, particles_per_point,                     &
                particle_type, pdx, pdy, pdz,                                  &
                prt_count, psb, psl, psn, psr, pss, pst,                       &
                radius, random_start_position, read_particles_from_restartfile,&
                seed_follows_topography, sgs_wf_part, sort_count,              &
                total_number_of_particles,                                     &
                use_sgs_for_particles,                                         &
                write_particle_statistics, uniform_particles, zero_particle,   &
                z0_av_global

    USE pegrid

    USE random_function_mod,                                                   &
        ONLY:  random_function

    IMPLICIT NONE

    PRIVATE

    INTEGER(iwp), PARAMETER         :: PHASE_INIT    = 1  !<
    INTEGER(iwp), PARAMETER, PUBLIC :: PHASE_RELEASE = 2  !<

    INTERFACE lpm_init
       MODULE PROCEDURE lpm_init
    END INTERFACE lpm_init

    INTERFACE lpm_create_particle
       MODULE PROCEDURE lpm_create_particle
    END INTERFACE lpm_create_particle

    PUBLIC lpm_init, lpm_create_particle

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_init

    USE lpm_collision_kernels_mod,                                             &
        ONLY:  init_kernels

    IMPLICIT NONE

    INTEGER(iwp) ::  i                           !<
    INTEGER(iwp) ::  j                           !<
    INTEGER(iwp) ::  k                           !<

#if defined( __parallel )
    INTEGER(iwp), DIMENSION(3) ::  blocklengths  !<
    INTEGER(iwp), DIMENSION(3) ::  displacements !<
    INTEGER(iwp), DIMENSION(3) ::  types         !<
#endif

    REAL(wp) ::  height_int                      !<
    REAL(wp) ::  height_p                        !<
    REAL(wp) ::  z_p                             !<
    REAL(wp) ::  z0_av_local                     !<

#if defined( __parallel )
!
!-- Define MPI derived datatype for FORTRAN datatype particle_type (see module
!-- particle_attributes). Integer length is 4 byte, Real is 8 byte
    blocklengths(1)  = 19;  blocklengths(2)  =   6;  blocklengths(3)  =   1
    displacements(1) =  0;  displacements(2) = 152;  displacements(3) = 176

    types(1) = MPI_REAL
    types(2) = MPI_INTEGER
    types(3) = MPI_UB
    CALL MPI_TYPE_STRUCT( 3, blocklengths, displacements, types, &
                          mpi_particle_type, ierr )
    CALL MPI_TYPE_COMMIT( mpi_particle_type, ierr )
#endif

!
!-- In case of oceans runs, the vertical index calculations need an offset,
!-- because otherwise the k indices will become negative
    IF ( ocean )  THEN
       offset_ocean_nzt    = nzt
       offset_ocean_nzt_m1 = nzt - 1
    ENDIF

!
!-- Define block offsets for dividing a gridcell in 8 sub cells

    block_offset(0) = block_offset_def (-1,-1,-1)
    block_offset(1) = block_offset_def (-1,-1, 0)
    block_offset(2) = block_offset_def (-1, 0,-1)
    block_offset(3) = block_offset_def (-1, 0, 0)
    block_offset(4) = block_offset_def ( 0,-1,-1)
    block_offset(5) = block_offset_def ( 0,-1, 0)
    block_offset(6) = block_offset_def ( 0, 0,-1)
    block_offset(7) = block_offset_def ( 0, 0, 0)
!
!-- Check the number of particle groups.
    IF ( number_of_particle_groups > max_number_of_particle_groups )  THEN
       WRITE( message_string, * ) 'max_number_of_particle_groups =',      &
                                  max_number_of_particle_groups ,         &
                                  '&number_of_particle_groups reset to ', &
                                  max_number_of_particle_groups
       CALL message( 'lpm_init', 'PA0213', 0, 1, 0, 6, 0 )
       number_of_particle_groups = max_number_of_particle_groups
    ENDIF

!
!-- Set default start positions, if necessary
    IF ( psl(1) == 9999999.9_wp )  psl(1) = -0.5_wp * dx
    IF ( psr(1) == 9999999.9_wp )  psr(1) = ( nx + 0.5_wp ) * dx
    IF ( pss(1) == 9999999.9_wp )  pss(1) = -0.5_wp * dy
    IF ( psn(1) == 9999999.9_wp )  psn(1) = ( ny + 0.5_wp ) * dy
    IF ( psb(1) == 9999999.9_wp )  psb(1) = zu(nz/2)
    IF ( pst(1) == 9999999.9_wp )  pst(1) = psb(1)

    IF ( pdx(1) == 9999999.9_wp  .OR.  pdx(1) == 0.0_wp )  pdx(1) = dx
    IF ( pdy(1) == 9999999.9_wp  .OR.  pdy(1) == 0.0_wp )  pdy(1) = dy
    IF ( pdz(1) == 9999999.9_wp  .OR.  pdz(1) == 0.0_wp )  pdz(1) = zu(2) - zu(1)

    DO  j = 2, number_of_particle_groups
       IF ( psl(j) == 9999999.9_wp )  psl(j) = psl(j-1)
       IF ( psr(j) == 9999999.9_wp )  psr(j) = psr(j-1)
       IF ( pss(j) == 9999999.9_wp )  pss(j) = pss(j-1)
       IF ( psn(j) == 9999999.9_wp )  psn(j) = psn(j-1)
       IF ( psb(j) == 9999999.9_wp )  psb(j) = psb(j-1)
       IF ( pst(j) == 9999999.9_wp )  pst(j) = pst(j-1)
       IF ( pdx(j) == 9999999.9_wp  .OR.  pdx(j) == 0.0_wp )  pdx(j) = pdx(j-1)
       IF ( pdy(j) == 9999999.9_wp  .OR.  pdy(j) == 0.0_wp )  pdy(j) = pdy(j-1)
       IF ( pdz(j) == 9999999.9_wp  .OR.  pdz(j) == 0.0_wp )  pdz(j) = pdz(j-1)
    ENDDO

!
!-- Allocate arrays required for calculating particle SGS velocities. 
!-- Initialize prefactor required for stoachastic Weil equation.
    IF ( use_sgs_for_particles  .AND.  .NOT. cloud_droplets )  THEN
       ALLOCATE( de_dx(nzb:nzt+1,nysg:nyng,nxlg:nxrg), &
                 de_dy(nzb:nzt+1,nysg:nyng,nxlg:nxrg), &
                 de_dz(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

       sgs_wf_part = 1.0_wp / 3.0_wp    
    ENDIF

!
!-- Allocate array required for logarithmic vertical interpolation of 
!-- horizontal particle velocities between the surface and the first vertical
!-- grid level. In order to avoid repeated CPU cost-intensive CALLS of 
!-- intrinsic FORTRAN procedure LOG(z/z0), LOG(z/z0) is precalculated for 
!-- several heights. Splitting into 20 sublayers turned out to be sufficient. 
!-- To obtain exact height levels of particles, linear interpolation is applied
!-- (see lpm_advec.f90).
    IF ( constant_flux_layer )  THEN
       
       ALLOCATE ( log_z_z0(0:number_of_sublayers) ) 
       z_p         = zu(nzb+1) - zw(nzb)

!
!--    Calculate horizontal mean value of z0 used for logartihmic
!--    interpolation. Note: this is not exact for heterogeneous z0. 
!--    However, sensitivity studies showed that the effect is 
!--    negligible. 
       z0_av_local  = SUM( z0(nys:nyn,nxl:nxr) )
       z0_av_global = 0.0_wp

#if defined( __parallel )
       CALL MPI_ALLREDUCE(z0_av_local, z0_av_global, 1, MPI_REAL, MPI_SUM, &
                          comm2d, ierr )
#else
       z0_av_global = z0_av_local
#endif

       z0_av_global = z0_av_global  / ( ( ny + 1 ) * ( nx + 1 ) )
!
!--    Horizontal wind speed is zero below and at z0
       log_z_z0(0) = 0.0_wp
!
!--    Calculate vertical depth of the sublayers
       height_int  = ( z_p - z0_av_global ) / REAL( number_of_sublayers, KIND=wp )
!
!--    Precalculate LOG(z/z0)
       height_p    = z0_av_global
       DO  k = 1, number_of_sublayers

          height_p    = height_p + height_int
          log_z_z0(k) = LOG( height_p / z0_av_global )

       ENDDO

    ENDIF

!
!-- Check boundary condition and set internal variables
    SELECT CASE ( bc_par_b )
    
       CASE ( 'absorb' )
          ibc_par_b = 1

       CASE ( 'reflect' )
          ibc_par_b = 2
          
       CASE DEFAULT
          WRITE( message_string, * )  'unknown boundary condition ',   &
                                       'bc_par_b = "', TRIM( bc_par_b ), '"'
          CALL message( 'lpm_init', 'PA0217', 1, 2, 0, 6, 0 )
          
    END SELECT
    SELECT CASE ( bc_par_t )
    
       CASE ( 'absorb' )
          ibc_par_t = 1

       CASE ( 'reflect' )
          ibc_par_t = 2
          
       CASE DEFAULT
          WRITE( message_string, * ) 'unknown boundary condition ',   &
                                     'bc_par_t = "', TRIM( bc_par_t ), '"'
          CALL message( 'lpm_init', 'PA0218', 1, 2, 0, 6, 0 )
          
    END SELECT
    SELECT CASE ( bc_par_lr )

       CASE ( 'cyclic' )
          ibc_par_lr = 0

       CASE ( 'absorb' )
          ibc_par_lr = 1

       CASE ( 'reflect' )
          ibc_par_lr = 2
          
       CASE DEFAULT
          WRITE( message_string, * ) 'unknown boundary condition ',   &
                                     'bc_par_lr = "', TRIM( bc_par_lr ), '"'
          CALL message( 'lpm_init', 'PA0219', 1, 2, 0, 6, 0 )
          
    END SELECT
    SELECT CASE ( bc_par_ns )

       CASE ( 'cyclic' )
          ibc_par_ns = 0

       CASE ( 'absorb' )
          ibc_par_ns = 1

       CASE ( 'reflect' )
          ibc_par_ns = 2
          
       CASE DEFAULT
          WRITE( message_string, * ) 'unknown boundary condition ',   &
                                     'bc_par_ns = "', TRIM( bc_par_ns ), '"'
          CALL message( 'lpm_init', 'PA0220', 1, 2, 0, 6, 0 )
          
    END SELECT

!
!-- Initialize collision kernels
    IF ( collision_kernel /= 'none' )  CALL init_kernels

!
!-- For the first model run of a possible job chain initialize the
!-- particles, otherwise read the particle data from restart file.
    IF ( TRIM( initializing_actions ) == 'read_restart_data'  &
         .AND.  read_particles_from_restartfile )  THEN

       CALL lpm_read_restart_file

    ELSE

!
!--    Allocate particle arrays and set attributes of the initial set of
!--    particles, which can be also periodically released at later times.
       ALLOCATE( prt_count(nzb:nzt+1,nysg:nyng,nxlg:nxrg), &
                 grid_particles(nzb+1:nzt,nys:nyn,nxl:nxr) )

       maximum_number_of_particles = 0
       number_of_particles         = 0

       sort_count = 0
       prt_count  = 0

!
!--    Initialize all particles with dummy values (otherwise errors may
!--    occur within restart runs). The reason for this is still not clear
!--    and may be presumably caused by errors in the respective user-interface.
       zero_particle = particle_type( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  &
                                      0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0, 0, 0, &
                                      0, .FALSE., -1 )

       particle_groups = particle_groups_type( 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp )

!
!--    Set values for the density ratio and radius for all particle
!--    groups, if necessary
       IF ( density_ratio(1) == 9999999.9_wp )  density_ratio(1) = 0.0_wp
       IF ( radius(1)        == 9999999.9_wp )  radius(1) = 0.0_wp
       DO  i = 2, number_of_particle_groups
          IF ( density_ratio(i) == 9999999.9_wp )  THEN
             density_ratio(i) = density_ratio(i-1)
          ENDIF
          IF ( radius(i) == 9999999.9_wp )  radius(i) = radius(i-1)
       ENDDO

       DO  i = 1, number_of_particle_groups
          IF ( density_ratio(i) /= 0.0_wp  .AND.  radius(i) == 0 )  THEN
             WRITE( message_string, * ) 'particle group #', i, 'has a', &
                                        'density ratio /= 0 but radius = 0'
             CALL message( 'lpm_init', 'PA0215', 1, 2, 0, 6, 0 )
          ENDIF
          particle_groups(i)%density_ratio = density_ratio(i)
          particle_groups(i)%radius        = radius(i)
       ENDDO

!
!--    Set a seed value for the random number generator to be exclusively
!--    used for the particle code. The generated random numbers should be
!--    different on the different PEs.
       iran_part = iran_part + myid

       CALL lpm_create_particle (PHASE_INIT)
!
!--    User modification of initial particles
       CALL user_lpm_init

!
!--    Open file for statistical informations about particle conditions
       IF ( write_particle_statistics )  THEN
          CALL check_open( 80 )
          WRITE ( 80, 8000 )  current_timestep_number, simulated_time,         &
                              number_of_particles,                             &
                              maximum_number_of_particles
          CALL close_file( 80 )
       ENDIF

    ENDIF

!
!-- To avoid programm abort, assign particles array to the local version of 
!-- first grid cell
    number_of_particles = prt_count(nzb+1,nys,nxl)
    particles => grid_particles(nzb+1,nys,nxl)%particles(1:number_of_particles)
!
!-- Formats
8000 FORMAT (I6,1X,F7.2,4X,I10,71X,I10)

 END SUBROUTINE lpm_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_create_particle (phase)

    USE lpm_exchange_horiz_mod,                                                &
        ONLY: lpm_exchange_horiz, lpm_move_particle, realloc_particles_array

    USE lpm_pack_arrays_mod,                                                   &
        ONLY: lpm_pack_all_arrays

    USE particle_attributes,                                                   &
        ONLY: deleted_particles, monodisperse_aerosols

    IMPLICIT  NONE

    INTEGER(iwp)               ::  alloc_size  !< relative increase of allocated memory for particles
    INTEGER(iwp)               ::  i           !< loop variable ( particle groups )
    INTEGER(iwp)               ::  ip          !< index variable along x
    INTEGER(iwp)               ::  j           !< loop variable ( particles per point )
    INTEGER(iwp)               ::  jp          !< index variable along y
    INTEGER(iwp)               ::  kp          !< index variable along z
    INTEGER(iwp)               ::  loop_stride !< loop variable for initialization
    INTEGER(iwp)               ::  n           !< loop variable ( number of particles )
    INTEGER(iwp)               ::  new_size    !< new size of allocated memory for particles

    INTEGER(iwp), INTENT(IN)   ::  phase       !< mode of inititialization

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  local_count !< start address of new particle
    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  local_start !< start address of new particle

    LOGICAL                    ::  first_stride !< flag for initialization

    REAL(wp)                   ::  pos_x      !< increment for particle position in x      
    REAL(wp)                   ::  pos_y      !< increment for particle position in y  
    REAL(wp)                   ::  pos_z      !< increment for particle position in z     
    REAL(wp)                   ::  rand_contr !< dummy argument for random position

    TYPE(particle_type),TARGET ::  tmp_particle !< temporary particle used for initialization

!
!-- Calculate particle positions and store particle attributes, if
!-- particle is situated on this PE
    DO  loop_stride = 1, 2
       first_stride = (loop_stride == 1)
       IF ( first_stride )   THEN
          local_count = 0           ! count number of particles
       ELSE
          local_count = prt_count   ! Start address of new particles
       ENDIF

       n = 0
       DO  i = 1, number_of_particle_groups

          pos_z = psb(i)

          DO WHILE ( pos_z <= pst(i) )

             pos_y = pss(i)

             DO WHILE ( pos_y <= psn(i) )

                IF ( pos_y >= ( nys - 0.5_wp ) * dy  .AND.  &
                     pos_y <  ( nyn + 0.5_wp ) * dy )  THEN

                   pos_x = psl(i)

            xloop: DO WHILE ( pos_x <= psr(i) )

                      IF ( pos_x >= ( nxl - 0.5_wp ) * dx  .AND.  &
                           pos_x <  ( nxr + 0.5_wp ) * dx )  THEN

                         DO  j = 1, particles_per_point

                            n = n + 1
                            tmp_particle%x             = pos_x
                            tmp_particle%y             = pos_y
                            tmp_particle%z             = pos_z
                            tmp_particle%age           = 0.0_wp
                            tmp_particle%age_m         = 0.0_wp
                            tmp_particle%dt_sum        = 0.0_wp
                            tmp_particle%dvrp_psize    = 0.0_wp !unused
                            tmp_particle%e_m           = 0.0_wp
                            IF ( curvature_solution_effects )  THEN
!
!--                            Initial values (internal timesteps, derivative)
!--                            for Rosenbrock method
                               tmp_particle%rvar1      = 1.0E-6_wp     !last Rosenbrock timestep
                               tmp_particle%rvar2      = 0.1E-6_wp     !dry aerosol radius
                               tmp_particle%rvar3      = -9999999.9_wp !unused
                            ELSE
!
!--                            Initial values for SGS velocities
                               tmp_particle%rvar1      = 0.0_wp
                               tmp_particle%rvar2      = 0.0_wp
                               tmp_particle%rvar3      = 0.0_wp
                            ENDIF
                            tmp_particle%speed_x       = 0.0_wp
                            tmp_particle%speed_y       = 0.0_wp
                            tmp_particle%speed_z       = 0.0_wp
                            tmp_particle%origin_x      = pos_x
                            tmp_particle%origin_y      = pos_y
                            tmp_particle%origin_z      = pos_z
                            tmp_particle%radius        = particle_groups(i)%radius
                            tmp_particle%weight_factor = initial_weighting_factor
                            tmp_particle%class         = 1
                            tmp_particle%group         = i
                            tmp_particle%tailpoints    = 0     !unused
                            tmp_particle%particle_mask = .TRUE.
                            tmp_particle%tail_id       = 0     !unused


!
!--                         Determine the grid indices of the particle position
                            ip = ( tmp_particle%x + 0.5_wp * dx ) * ddx
                            jp = ( tmp_particle%y + 0.5_wp * dy ) * ddy
                            kp = tmp_particle%z / dz + 1 + offset_ocean_nzt

                            IF ( seed_follows_topography )  THEN
!
!--                            Particle height is given relative to topography
                               kp = kp + nzb_w_inner(jp,ip)
                               tmp_particle%z = tmp_particle%z +               &
                                                         zw(nzb_w_inner(jp,ip))
                               IF ( kp > nzt )  THEN
                                  pos_x = pos_x + pdx(i)
                                  CYCLE xloop
                               ENDIF
                            ELSEIF ( .NOT. seed_follows_topography .AND.       &
                                     tmp_particle%z <= zw(nzb_w_inner(jp,ip)) )  THEN
                               pos_x = pos_x + pdx(i)
                               CYCLE xloop                               
                            ENDIF

                            local_count(kp,jp,ip) = local_count(kp,jp,ip) + 1
                            IF ( .NOT. first_stride )  THEN
                               IF ( ip < nxl  .OR.  jp < nys  .OR.  kp < nzb+1 )  THEN
                                  write(6,*) 'xl ',ip,jp,kp,nxl,nys,nzb+1
                               ENDIF
                               IF ( ip > nxr  .OR.  jp > nyn  .OR.  kp > nzt )  THEN
                                  write(6,*) 'xu ',ip,jp,kp,nxr,nyn,nzt
                               ENDIF
                               grid_particles(kp,jp,ip)%particles(local_count(kp,jp,ip)) = tmp_particle

                            ENDIF
                         ENDDO

                      ENDIF

                      pos_x = pos_x + pdx(i)

                   ENDDO xloop

                ENDIF

                pos_y = pos_y + pdy(i)

             ENDDO

             pos_z = pos_z + pdz(i)

          ENDDO

       ENDDO

       IF ( first_stride )  THEN
          DO  ip = nxl, nxr
             DO  jp = nys, nyn
                DO  kp = nzb+1, nzt
                   IF ( phase == PHASE_INIT )  THEN
                      IF ( local_count(kp,jp,ip) > 0 )  THEN
                         alloc_size = MAX( INT( local_count(kp,jp,ip) *        &
                            ( 1.0_wp + alloc_factor / 100.0_wp ) ),            &
                            min_nr_particle )
                      ELSE
                         alloc_size = min_nr_particle
                      ENDIF
                      ALLOCATE(grid_particles(kp,jp,ip)%particles(1:alloc_size))
                      DO  n = 1, alloc_size
                         grid_particles(kp,jp,ip)%particles(n) = zero_particle
                      ENDDO
                   ELSEIF ( phase == PHASE_RELEASE )  THEN
                      IF ( local_count(kp,jp,ip) > 0 )  THEN
                         new_size   = local_count(kp,jp,ip) + prt_count(kp,jp,ip)
                         alloc_size = MAX( INT( new_size * ( 1.0_wp +          &
                            alloc_factor / 100.0_wp ) ), min_nr_particle )
                         IF( alloc_size > SIZE( grid_particles(kp,jp,ip)%particles) )  THEN
                           CALL realloc_particles_array(ip,jp,kp,alloc_size)
                         ENDIF
                      ENDIF
                   ENDIF

                ENDDO
             ENDDO
          ENDDO
       ENDIF

    ENDDO

    local_start = prt_count+1
    prt_count   = local_count

!
!-- Initialize aerosol background spectrum
    IF ( curvature_solution_effects  .AND.  .NOT. monodisperse_aerosols )  THEN
       CALL lpm_init_aerosols(local_start)
    ENDIF

!
!-- Add random fluctuation to particle positions.
    IF ( random_start_position )  THEN
       DO  ip = nxl, nxr
          DO  jp = nys, nyn
             DO  kp = nzb+1, nzt
                number_of_particles = prt_count(kp,jp,ip)
                IF ( number_of_particles <= 0 )  CYCLE
                particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
!
!--             Move only new particles. Moreover, limit random fluctuation 
!--             in order to prevent that particles move more than one grid box, 
!--             which would lead to problems concerning particle exchange 
!--             between processors in case pdx/pdy are larger than dx/dy, 
!--             respectively.  
                DO  n = local_start(kp,jp,ip), number_of_particles
                   IF ( psl(particles(n)%group) /= psr(particles(n)%group) )  THEN
                      rand_contr = ( random_function( iran_part ) - 0.5_wp ) * &
                                     pdx(particles(n)%group)
                      particles(n)%x = particles(n)%x +                        &
                              MERGE( rand_contr, SIGN( dx, rand_contr ), &
                                     ABS( rand_contr ) < dx                    &
                                   ) 
                   ENDIF
                   IF ( pss(particles(n)%group) /= psn(particles(n)%group) )  THEN
                      rand_contr = ( random_function( iran_part ) - 0.5_wp ) * &
                                     pdy(particles(n)%group)
                      particles(n)%y = particles(n)%y +                        &
                              MERGE( rand_contr, SIGN( dy, rand_contr ), &
                                     ABS( rand_contr ) < dy                    &
                                   ) 
                   ENDIF
                   IF ( psb(particles(n)%group) /= pst(particles(n)%group) )  THEN
                      rand_contr = ( random_function( iran_part ) - 0.5_wp ) * &
                                     pdz(particles(n)%group)
                      particles(n)%z = particles(n)%z +                        &
                              MERGE( rand_contr, SIGN( dz, rand_contr ), &
                                     ABS( rand_contr ) < dz                    &
                                   ) 
                   ENDIF
                ENDDO
!
!--             Identify particles located outside the model domain and reflect
!--             or absorb them if necessary.
                CALL lpm_boundary_conds( 'bottom/top' )
!
!--             Furthermore, remove particles located in topography. Note, as
!--             the particle speed is still zero at this point, wall 
!--             reflection boundary conditions will not work in this case.
                particles =>                                                   &
                       grid_particles(kp,jp,ip)%particles(1:number_of_particles)
                DO  n = local_start(kp,jp,ip), number_of_particles
                   i = ( particles(n)%x + 0.5_wp * dx ) * ddx
                   j = ( particles(n)%y + 0.5_wp * dy ) * ddy
                   IF ( particles(n)%z <= zw(nzb_w_inner(j,i)) )  THEN
                      particles(n)%particle_mask = .FALSE.
                      deleted_particles = deleted_particles + 1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
!
!--    Exchange particles between grid cells and processors
       CALL lpm_move_particle
       CALL lpm_exchange_horiz

    ENDIF
!
!-- In case of random_start_position, delete particles identified by 
!-- lpm_exchange_horiz and lpm_boundary_conds. Then sort particles into blocks, 
!-- which is needed for a fast interpolation of the LES fields on the particle 
!-- position.
    CALL lpm_pack_all_arrays

!
!-- Determine maximum number of particles (i.e., all possible particles that 
!-- have been allocated) and the current number of particles
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             maximum_number_of_particles = maximum_number_of_particles         &
                                           + SIZE(grid_particles(kp,jp,ip)%particles)
             number_of_particles         = number_of_particles                 &
                                           + prt_count(kp,jp,ip)
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculate the number of particles of the total domain
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( number_of_particles, total_number_of_particles, 1, &
    MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
    total_number_of_particles = number_of_particles
#endif

    RETURN

 END SUBROUTINE lpm_create_particle

 SUBROUTINE lpm_init_aerosols(local_start)

    USE arrays_3d,                                                             &
        ONLY: hyp, pt, q  

    USE cloud_parameters,                                                      &
        ONLY: l_d_rv, rho_l

    USE constants,                                                             &
        ONLY: pi

    USE kinds

    USE particle_attributes,                                                   &
        ONLY: init_aerosol_probabilistic, molecular_weight_of_solute,          &
              molecular_weight_of_water, n1, n2, n3, rho_s, rm1, rm2, rm3,     &
              s1, s2, s3, vanthoff

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  cdf     !< CDF of aerosol spectrum 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  r_temp  !< dry aerosol radius spectrum

    REAL(wp)  :: bfactor            !< solute effects
    REAL(wp)  :: dr                 !< width of radius bin
    REAL(wp)  :: e_a                !< vapor pressure
    REAL(wp)  :: e_s                !< saturation vapor pressure
    REAL(wp)  :: n_init             !< sum of all aerosol concentrations
    REAL(wp)  :: pdf                !< PDF of aerosol spectrum
    REAL(wp)  :: rmin = 1.0e-8_wp   !< minimum aerosol radius
    REAL(wp)  :: rmax = 1.0e-6_wp   !< maximum aerosol radius
    REAL(wp)  :: rs_rand            !< random number
    REAL(wp)  :: r_mid              !< mean radius
    REAL(wp)  :: t_int              !< temperature
    REAL(wp)  :: weight_sum         !< sum of all weighting factors

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  local_start !<

    INTEGER(iwp)  :: n              !<
    INTEGER(iwp)  :: nn             !<
    INTEGER(iwp)  :: no_bins = 999  !< number of bins
    INTEGER(iwp)  :: ip             !<
    INTEGER(iwp)  :: jp             !<
    INTEGER(iwp)  :: kp             !<

    LOGICAL ::  new_pdf = .FALSE.   !< check if aerosol PDF has to be recalculated

!
!-- Compute aerosol background distribution
    IF ( init_aerosol_probabilistic )  THEN
       ALLOCATE( cdf(0:no_bins), r_temp(0:no_bins) )
       DO n = 0, no_bins
          r_temp(n) = EXP( LOG(rmin) + ( LOG(rmax) - LOG(rmin ) ) /            &
                           REAL(no_bins, KIND=wp) * REAL(n, KIND=wp) )

          cdf(n) = 0.0_wp
          n_init = n1 + n2 + n3
          IF ( n1 > 0.0_wp )  THEN 
             cdf(n) = cdf(n) + n1 / n_init * ( 0.5_wp + 0.5_wp *        &
                                  ERF( LOG( r_temp(n) / rm1 ) /         &
                                       ( SQRT(2.0_wp) * LOG(s1) )       &
                                     ) )
          ENDIF
          IF ( n2 > 0.0_wp )  THEN 
             cdf(n) = cdf(n) + n2 / n_init * ( 0.5_wp + 0.5_wp *        &
                                  ERF( LOG( r_temp(n) / rm2 ) /         &
                                       ( SQRT(2.0_wp) * LOG(s2) )       &
                                     ) )
          ENDIF
          IF ( n3 > 0.0_wp )  THEN 
             cdf(n) = cdf(n) + n3 / n_init * ( 0.5_wp + 0.5_wp *        &
                                  ERF( LOG( r_temp(n) / rm3 ) /         &
                                       ( SQRT(2.0_wp) * LOG(s3) )       &
                                     ) )
          ENDIF

       ENDDO
    ENDIF

    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt

             number_of_particles = prt_count(kp,jp,ip)
             IF ( number_of_particles <= 0 )  CYCLE
             particles => grid_particles(kp,jp,ip)%particles(1:number_of_particles)
!
!--          Initialize the aerosols with a predefined spectral distribution
!--          of the dry radius (logarithmically increasing bins) and a varying 
!--          weighting factor
             IF ( .NOT. init_aerosol_probabilistic )  THEN

                new_pdf = .FALSE.
                IF ( .NOT. ALLOCATED( r_temp ) )  THEN
                   new_pdf = .TRUE.
                ELSE
                   IF ( SIZE( r_temp ) .NE. &
                        number_of_particles - local_start(kp,jp,ip) + 2 )  THEN
                      new_pdf = .TRUE.
                      DEALLOCATE( r_temp )
                   ENDIF
                ENDIF

                IF ( new_pdf )  THEN

                   no_bins = number_of_particles + 1 - local_start(kp,jp,ip)
                   ALLOCATE( r_temp(0:no_bins) )

                   DO n = 0, no_bins
                      r_temp(n) = EXP( LOG(rmin) + ( LOG(rmax) - LOG(rmin ) ) / &
                                       REAL(no_bins, KIND=wp) *                 &
                                       REAL(n, KIND=wp) )
                   ENDDO

                ENDIF

!
!--             Calculate radius and concentration of each aerosol
                DO n = local_start(kp,jp,ip), number_of_particles

                   nn = n - local_start(kp,jp,ip)

                   r_mid = SQRT( r_temp(nn) * r_temp(nn+1) )
                   dr    = r_temp(nn+1) - r_temp(nn)

                   pdf    = 0.0_wp
                   n_init = n1 + n2 + n3
                   IF ( n1 > 0.0_wp )  THEN 
                      pdf = pdf + n1 / n_init * ( 1.0_wp / ( r_mid * LOG(s1) *      &
                                                             SQRT( 2.0_wp * pi )    &
                                                           ) *                      &
                                                  EXP( -( LOG( r_mid / rm1 ) )**2 / &
                                                       ( 2.0_wp * LOG(s1)**2 )      &
                                                     )                              &
                                                )
                   ENDIF
                   IF ( n2 > 0.0_wp )  THEN
                      pdf = pdf + n2 / n_init * ( 1.0_wp / ( r_mid * LOG(s2) *      &
                                                             SQRT( 2.0_wp * pi )    &
                                                           ) *                      &
                                                  EXP( -( LOG( r_mid / rm2 ) )**2 / &
                                                       ( 2.0_wp * LOG(s2)**2 )      &
                                                     )                              &
                                                )
                   ENDIF
                   IF ( n3 > 0.0_wp )  THEN
                      pdf = pdf + n3 / n_init * ( 1.0_wp / ( r_mid * LOG(s3) *      &
                                                             SQRT( 2.0_wp * pi )    &
                                                           ) *                      &
                                                  EXP( -( LOG( r_mid / rm3 ) )**2 / &
                                                       ( 2.0_wp * LOG(s3)**2 )      &
                                                     )                              &
                                                )
                   ENDIF

                   particles(n)%rvar2         = r_mid
                   particles(n)%weight_factor = pdf * dr

                END DO
!
!--             Adjust weighting factors to initialize the same number of aerosols
!--             in every grid box
                weight_sum = SUM(particles(local_start(kp,jp,ip):number_of_particles)%weight_factor)

                particles(local_start(kp,jp,ip):number_of_particles)%weight_factor =     &
                   particles(local_start(kp,jp,ip):number_of_particles)%weight_factor /  &
                   weight_sum * initial_weighting_factor * ( no_bins + 1 )

             ENDIF
!
!--          Initialize the aerosols with a predefined weighting factor but 
!--          a randomly choosen dry radius
             IF ( init_aerosol_probabilistic )  THEN

                DO  n = local_start(kp,jp,ip), number_of_particles  !only new particles 

                   rs_rand = -1.0_wp
                   DO WHILE ( rs_rand .LT. cdf(0)  .OR.  rs_rand .GE. cdf(no_bins)  )
                      rs_rand = random_function( iran_part )
                   ENDDO
!
!--                Determine aerosol dry radius by a random number generator
                   DO nn = 0, no_bins-1
                      IF ( cdf(nn) .LE. rs_rand  .AND.  cdf(nn+1) .GT. rs_rand )  THEN
                         particles(n)%rvar2 = r_temp(nn) + ( r_temp(nn+1) - r_temp(nn) ) / &
                                              ( cdf(nn+1) - cdf(nn) ) * ( rs_rand - cdf(nn) )
                         EXIT
                      ENDIF
                   ENDDO

                ENDDO

             ENDIF

!
!--          Set particle radius to equilibrium radius based on the environmental
!--          supersaturation (Khvorostyanov and Curry, 2007, JGR). This avoids
!--          the sometimes lengthy growth toward their equilibrium radius within 
!--          the simulation.
             t_int  = pt(kp,jp,ip) * ( hyp(kp) / 100000.0_wp )**0.286_wp

             e_s = 611.0_wp * EXP( l_d_rv * ( 3.6609E-3_wp - 1.0_wp / t_int ) )
             e_a = q(kp,jp,ip) * hyp(kp) / ( 0.378_wp * q(kp,jp,ip) + 0.622_wp )

!
!--          The formula is only valid for subsaturated environments. For 
!--          supersaturations higher than -1 %, the supersaturation is set to -1%.
             IF ( e_a / e_s < 0.99_wp )  THEN

                DO  n = local_start(kp,jp,ip), number_of_particles  !only new particles 

                   bfactor             = vanthoff * molecular_weight_of_water *    &
                                         rho_s * particles(n)%rvar2**3 /           &
                                         ( molecular_weight_of_solute * rho_l )
                   particles(n)%radius = particles(n)%rvar2 * ( bfactor /          &
                                         particles(n)%rvar2**3 )**(1.0_wp/3.0_wp) *&
                                         ( 1.0_wp - e_a / e_s )**(-1.0_wp/3.0_wp)

                ENDDO

             ELSE

                DO  n = local_start(kp,jp,ip), number_of_particles  !only new particles 

                   bfactor             = vanthoff * molecular_weight_of_water *    &
                                         rho_s * particles(n)%rvar2**3 /           &
                                         ( molecular_weight_of_solute * rho_l )
                   particles(n)%radius = particles(n)%rvar2 * ( bfactor /          &
                                         particles(n)%rvar2**3 )**(1.0_wp/3.0_wp) *&
                                         0.01_wp**(-1.0_wp/3.0_wp)

                ENDDO

             ENDIF

          ENDDO
       ENDDO
    ENDDO
!
!-- Deallocate used arrays
    IF ( ALLOCATED(r_temp) )  DEALLOCATE( r_temp )
    IF ( ALLOCATED(cdf) )     DEALLOCATE( cdf )

 END SUBROUTINE lpm_init_aerosols

END MODULE lpm_init_mod
