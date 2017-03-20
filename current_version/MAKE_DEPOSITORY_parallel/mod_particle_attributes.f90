!> @file mod_particle_attributes.f90
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
! $Id: mod_particle_attributes.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1936 2016-06-13 13:37:44Z suehring
! +deallocate_memory, step_dealloc
!
! 1929 2016-06-09 16:25:25Z suehring
! -sgs_wfu_par, sgs_wfv_par, sgs_wfw_par
! + sgs_wf_par
!
! 1871 2016-04-15 11:46:09Z hoffmann
! Initialization of aerosols added.
!
! 1849 2016-04-08 11:33:18Z hoffmann
! bfactor, mass_of_solute, molecular_weight_of_solute, molecular_weight_of_water,
! vanthoff added from modules
!
! 1831 2016-04-07 13:15:51Z hoffmann
! palm_kernel removed, curvature_solution_effects added
!
! 1822 2016-04-07 07:49:42Z hoffmann
! +collision_algorithm, all_or_nothing, average_impact
! Tails removed.
!
! 1727 2015-11-20 07:22:02Z knoop
! Bugfix: Cause of syntax warning gfortran preprocessor removed 
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1575 2015-03-27 09:56:27Z raasch
! +seed_follows_topography
!
! 1359 2014-04-11 17:15:14Z hoffmann
! new module containing all particle related variables
! -dt_sort_particles
!
! Description:
! ------------
!> Definition of variables used to compute particle transport
!------------------------------------------------------------------------------!
MODULE particle_attributes
 

    USE kinds

    CHARACTER(LEN=15) ::  bc_par_lr = 'cyclic'                    !< left/right boundary condition
    CHARACTER(LEN=15) ::  bc_par_ns = 'cyclic'                    !< north/south boundary condition
    CHARACTER(LEN=15) ::  bc_par_b  = 'reflect'                   !< bottom boundary condition
    CHARACTER(LEN=15) ::  bc_par_t  = 'absorb'                    !< top boundary condition
    CHARACTER(LEN=15) ::  collision_algorithm = 'all_or_nothing'  !< collision algorithm
    CHARACTER(LEN=15) ::  collision_kernel = 'none'               !< collision kernel

    INTEGER(iwp) ::  deleted_particles = 0,                                    &
                     dissipation_classes = 10, ibc_par_lr,                     &
                     ibc_par_ns, ibc_par_b, ibc_par_t, iran_part = -1234567,   &
                     maximum_number_of_particles = 0,                          &
                     min_nr_particle = 50,                                     &
                     mpi_particle_type,                                        &
                     number_of_particles = 0,                                  &
                     number_of_particle_groups = 1,                            &
                     number_of_sublayers = 20,                                 &
                     offset_ocean_nzt = 0,                                     &
                     offset_ocean_nzt_m1 = 0, particles_per_point = 1,         &
                     particle_file_count = 0, radius_classes = 20,             &
                     sort_count = 0, step_dealloc = 100,                       &
                     total_number_of_particles,                                &
                     trlp_count_sum, trlp_count_recv_sum, trrp_count_sum,      &
                     trrp_count_recv_sum, trsp_count_sum, trsp_count_recv_sum, &
                     trnp_count_sum, trnp_count_recv_sum

    INTEGER(iwp), PARAMETER ::  max_number_of_particle_groups = 10

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  prt_count

    LOGICAL ::  all_or_nothing = .FALSE., average_impact = .FALSE.,            &
                curvature_solution_effects = .FALSE.,                          &
                deallocate_memory = .TRUE.,                                    &
                hall_kernel = .FALSE., particle_advection = .FALSE.,           &
                random_start_position = .FALSE.,                               &
                read_particles_from_restartfile = .TRUE.,                      &
                seed_follows_topography = .FALSE.,                             &
                uniform_particles = .TRUE., use_kernel_tables = .FALSE.,       &
                use_sgs_for_particles = .FALSE., wang_kernel = .FALSE.,        &
                write_particle_statistics = .FALSE.

    LOGICAL, DIMENSION(max_number_of_particle_groups) ::                       &
                vertical_particle_advection = .TRUE.

    REAL(wp) ::  alloc_factor = 20.0_wp, c_0 = 3.0_wp,                         &
                 dt_min_part = 0.0002_wp, dt_prel = 9999999.9_wp,              &
                 dt_write_particle_data = 9999999.9_wp,                        &
                 end_time_prel = 9999999.9_wp,                                 &
                 initial_weighting_factor = 1.0_wp,                            &
                 particle_advection_start = 0.0_wp,                            &
                 sgs_wf_part, time_prel = 0.0_wp, time_sort_particles = 0.0_wp,&
                 time_write_particle_data = 0.0_wp, z0_av_global

    REAL(wp), DIMENSION(max_number_of_particle_groups) ::                      &
                 density_ratio = 9999999.9_wp, pdx = 9999999.9_wp,             &
                 pdy = 9999999.9_wp, pdz = 9999999.9_wp, psb = 9999999.9_wp,   &
                 psl = 9999999.9_wp, psn = 9999999.9_wp, psr = 9999999.9_wp,   &
                 pss = 9999999.9_wp, pst = 9999999.9_wp, radius = 9999999.9_wp

    REAL(wp), DIMENSION(:), ALLOCATABLE     ::  log_z_z0

    TYPE particle_type
        SEQUENCE
        REAL(wp)     ::  radius, age, age_m, dt_sum, dvrp_psize, e_m,          &
                         origin_x, origin_y, origin_z, rvar1, rvar2, rvar3,    &
                         speed_x, speed_y, speed_z, weight_factor, x, y, z
        INTEGER(iwp) ::  class, group, tailpoints, tail_id
        LOGICAL      ::  particle_mask
        INTEGER(iwp) ::  block_nr
    END TYPE particle_type

    TYPE(particle_type), DIMENSION(:), POINTER ::  particles
    TYPE(particle_type)                        ::  zero_particle

    TYPE particle_groups_type
        SEQUENCE
        REAL(wp) ::  density_ratio, radius, exp_arg, exp_term
    END TYPE particle_groups_type

    TYPE(particle_groups_type), DIMENSION(max_number_of_particle_groups) ::    &
       particle_groups

    TYPE  grid_particle_def
        INTEGER(iwp), DIMENSION(0:7)               ::  start_index
        INTEGER(iwp), DIMENSION(0:7)               ::  end_index
        LOGICAL                                    ::  time_loop_done
        TYPE(particle_type), POINTER, DIMENSION(:) ::  particles                !Particle array for this grid cell
    END TYPE grid_particle_def

    TYPE(grid_particle_def), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  grid_particles

    TYPE block_offset_def
        INTEGER(iwp) ::  i_off
        INTEGER(iwp) ::  j_off
        INTEGER(iwp) ::  k_off
    END TYPE block_offset_def

    TYPE(block_offset_def), DIMENSION(0:7)         ::  block_offset

!
!-- Lagrangian cloud model constants (especially for steering aerosols)
    REAL(wp) ::  molecular_weight_of_solute = 0.05844_wp    !< mol. m. NaCl (kg mol-1)
    REAL(wp) ::  molecular_weight_of_water = 0.01801528_wp  !< mol. m. H2O (kg mol-1)
    REAL(wp) ::  rho_s = 2165.0_wp                          !< density of NaCl (kg m-3) 
    REAL(wp) ::  vanthoff = 2.0_wp                          !< van't Hoff factor for NaCl

    REAL(wp) ::  n1 = 100.0_wp, s1 = 2.0_wp, rm1 = 0.05E-6_wp, &
                 n2 =   0.0_wp, s2 = 2.0_wp, rm2 = 0.05E-6_wp, &
                 n3 =   0.0_wp, s3 = 2.0_wp, rm3 = 0.05E-6_wp

    LOGICAL  ::  monodisperse_aerosols      = .FALSE.       !< initialize monodisperse aerosols
    LOGICAL  ::  init_aerosol_probabilistic = .FALSE.       !< how to initialize aerosol spectra

    SAVE


END MODULE particle_attributes

