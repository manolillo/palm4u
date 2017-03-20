!> @file urban_surface_mod.f90
!--------------------------------------------------------------------------------!
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
! Copyright 2015-2016 Czech Technical University in Prague
! Copyright 1997-2016 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
! 
! 
! Former revisions:
! -----------------
! $Id: urban_surface_mod.f90 2032 2016-10-21 15:13:51Z knoop $
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
! 
! 2024 2016-10-12 16:42:37Z kanani
! Bugfixes in deallocation of array plantt and reading of csf/csfsurf,
! optimization of MPI-RMA operations,
! declaration of pcbl as integer,
! renamed usm_radnet -> usm_rad_net, usm_canopy_khf -> usm_canopy_hr,
! splitted arrays svf -> svf & csf, svfsurf -> svfsurf & csfsurf,
! use of new control parameter varnamelength,
! added output variables usm_rad_ressw, usm_rad_reslw,
! minor formatting changes,
! minor optimizations.
! 
! 2011 2016-09-19 17:29:57Z kanani
! Major reformatting according to PALM coding standard (comments, blanks,
! alphabetical ordering, etc.),
! removed debug_prints, 
! removed auxiliary SUBROUTINE get_usm_info, instead, USM flag urban_surface is
! defined in MODULE control_parameters (modules.f90) to avoid circular 
! dependencies,
! renamed canopy_heat_flux to pc_heating_rate, as meaning of quantity changed.
! 
! 2007 2016-08-24 15:47:17Z kanani
! Initial revision
!
!
! Description:
! ------------
! 2016/6/9 - Initial version of the USM (Urban Surface Model)
!            authors: Jaroslav Resler, Pavel Krc
!                     (Czech Technical University in Prague and Institute of
!                      Computer Science of the Czech Academy of Sciences, Prague)
!            with contributions: Michal Belda, Nina Benesova, Ondrej Vlcek
!            partly inspired by PALM LSM (B. Maronga)
!            parameterizations of Ra checked with TUF3D (E. S. Krayenhoff)
!> Module for Urban Surface Model (USM)
!> The module includes:
!>    1. radiation model with direct/diffuse radiation, shading, reflections
!>       and integration with plant canopy
!>    2. wall and wall surface model
!>    3. surface layer energy balance
!>    4. anthropogenic heat (only from transportation so far)
!>    5. necessary auxiliary subroutines (reading inputs, writing outputs,
!>       restart simulations, ...)
!> It also make use of standard radiation and integrates it into
!> urban surface model.
!>
!> Further work:
!> -------------
!> 1. Reduce number of shape view factors by merging factors for distant surfaces
!>    under shallow angles. Idea: Iteratively select the smallest shape view
!>    factor by value (among all sources and targets) which has a similarly
!>    oriented source neighbor (or near enough) SVF and merge them by adding
!>    value of the smaller SVF to the larger one and deleting the smaller one.
!>    This will allow for better scaling at higher resolutions.
!>
!> 2. Remove global arrays surfouts, surfoutl and only keep track of radiosity
!>    from surfaces that are visible from local surfaces (i.e. there is a SVF
!>    where target is local). To do that, radiosity will be exchanged after each
!>    reflection step using MPI_Alltoall instead of current MPI_Allgather.
!>
!> @todo Check optimizations for RMA operations
!> @todo Alternatives for MPI_WIN_ALLOCATE? (causes problems with openmpi)
!> @todo Check for load imbalances in CPU measures, e.g. for exchange_horiz_prog
!>       factor 3 between min and max time
!------------------------------------------------------------------------------!
 MODULE urban_surface_mod

    USE arrays_3d,                                                             &
        ONLY:  zu, pt, pt_1, pt_2, p, ol, shf, ts, us, u, v, w, hyp, tend

    USE cloud_parameters,                                                      &
        ONLY:  cp, r_d

    USE constants,                                                             &
        ONLY:  pi
    
    USE control_parameters,                                                    &
        ONLY:  dz, topography, dt_3d, intermediate_timestep_count,             &
               initializing_actions, intermediate_timestep_count_max,          &
               simulated_time, end_time, timestep_scheme, tsc,                 &
               coupling_char, io_blocks, io_group, message_string,             &
               time_since_reference_point, surface_pressure,                   &
               g, pt_surface, large_scale_forcing, lsf_surf,                   &
               time_do3d, dt_do3d, average_count_3d, varnamelength,            &
               urban_surface

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s
     
    USE grid_variables,                                                        &
        ONLY:  dx, dy, ddx, ddy, ddx2, ddy2
    
    USE indices,                                                               &
        ONLY:  nx, ny, nnx, nny, nnz, nxl, nxlg, nxr, nxrg, nyn, nyng, nys,    &
               nysg, nzb_s_inner, nzb_s_outer, nzb, nzt, nbgp

    USE, INTRINSIC :: iso_c_binding 

    USE kinds
              
    USE pegrid
    
    USE plant_canopy_model_mod,                                                &
        ONLY:  plant_canopy, pch_index,                                        &
               pc_heating_rate, lad_s
    
    USE radiation_model_mod,                                                   &
        ONLY:  radiation, calc_zenith, zenith, day_init, time_utc_init,        &
               rad_net, rad_sw_in, rad_lw_in, rad_sw_out, rad_lw_out,          &
               sigma_sb, sun_direction, sun_dir_lat, sun_dir_lon,              &
               force_radiation_call

    USE statistics,                                                            &
        ONLY:  hom, statistic_regions

               

    IMPLICIT NONE

!-- configuration parameters (they can be setup in PALM config)
    LOGICAL                                        ::  split_diffusion_radiation = .TRUE. !< split direct and diffusion dw radiation
                                                                                          !< (.F. in case the radiation model already does it)    
    LOGICAL                                        ::  usm_energy_balance_land = .TRUE.   !< flag parameter indicating wheather the energy balance is calculated for land and roofs
    LOGICAL                                        ::  usm_energy_balance_wall = .TRUE.   !< flag parameter indicating wheather the energy balance is calculated for land and roofs
    LOGICAL                                        ::  usm_material_model = .TRUE.        !< flag parameter indicating wheather the  model of heat in materials is used
    LOGICAL                                        ::  usm_anthropogenic_heat = .FALSE.   !< flag parameter indicating wheather the anthropogenic heat sources (e.g.transportation) are used
    LOGICAL                                        ::  force_radiation_call_l = .FALSE.   !< flag parameter for unscheduled radiation model calls
    LOGICAL                                        ::  mrt_factors = .FALSE.              !< whether to generate MRT factor files during init
    LOGICAL                                        ::  write_svf_on_init = .FALSE.
    LOGICAL                                        ::  read_svf_on_init = .FALSE.
    LOGICAL                                        ::  usm_lad_rma = .TRUE.               !< use MPI RMA to access LAD for raytracing (instead of global array)
    
    INTEGER(iwp)                                   ::  nrefsteps = 0                      !< number of reflection steps to perform
    
    INTEGER(iwp)                                   ::  land_category = 2                  !< default category for land surface
    INTEGER(iwp)                                   ::  wall_category = 2                  !< default category for wall surface over pedestrian zone
    INTEGER(iwp)                                   ::  pedestrant_category = 2            !< default category for wall surface in pedestrian zone
    INTEGER(iwp)                                   ::  roof_category = 2                  !< default category for root surface
    REAL(wp)                                       ::  roof_height_limit = 4._wp          !< height for distinguish between land surfaces and roofs

    REAL(wp), PARAMETER                            ::  ext_coef = 0.6_wp                  !< extinction coefficient (a.k.a. alpha)
    REAL(wp)                                       ::  ra_horiz_coef = 5.0_wp             !< mysterious coefficient for correction of overestimation
                                                                                          !< of r_a for horizontal surfaces -> TODO
    
!-- parameters of urban surface model
    INTEGER(iwp), PARAMETER                        ::  usm_version_len = 10               !< length of identification string of usm version
    CHARACTER(usm_version_len), PARAMETER          ::  usm_version = 'USM v. 1.0'         !< identification of version of binary svf and restart files
    INTEGER(iwp), PARAMETER                        ::  svf_code_len = 15                  !< length of code for verification of the end of svf file
    CHARACTER(svf_code_len), PARAMETER             ::  svf_code = '*** end svf ***'       !< code for verification of the end of svf file
    INTEGER(iwp)                                   ::  nzu                                !< number of layers of urban surface (will be calculated)
    INTEGER(iwp)                                   ::  nzub,nzut                          !< bottom and top layer of urban surface (will be calculated)
    INTEGER(iwp), PARAMETER                        ::  nzut_free = 3                      !< number of free layers in urban surface layer above top of buildings 
    INTEGER(iwp), PARAMETER                        ::  ndsvf = 2                          !< number of dimensions of real values in SVF
    INTEGER(iwp), PARAMETER                        ::  idsvf = 2                          !< number of dimensions of integer values in SVF
    INTEGER(iwp), PARAMETER                        ::  ndcsf = 2                          !< number of dimensions of real values in CSF
    INTEGER(iwp), PARAMETER                        ::  idcsf = 2                          !< number of dimensions of integer values in CSF
    INTEGER(iwp), PARAMETER                        ::  kdcsf = 4                          !< number of dimensions of integer values in CSF calculation array
    INTEGER(iwp), PARAMETER                        ::  id = 1                             !< position of d-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  iz = 2                             !< position of k-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  iy = 3                             !< position of j-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  ix = 4                             !< position of i-index in surfl and surf
    INTEGER(iwp), PARAMETER                        ::  iroof = 0                          !< 0 - index of ground or roof
    INTEGER(iwp), PARAMETER                        ::  isouth = 1                         !< 1 - index of south facing wall
    INTEGER(iwp), PARAMETER                        ::  inorth = 2                         !< 2 - index of north facing wall
    INTEGER(iwp), PARAMETER                        ::  iwest  = 3                         !< 3 - index of west facing wall
    INTEGER(iwp), PARAMETER                        ::  ieast  = 4                         !< 4 - index of east facing wall
    INTEGER(iwp), PARAMETER                        ::  isky = 5                           !< 5 - index of top border of the urban surface layer ("urban sky")
    INTEGER(iwp), PARAMETER                        ::  inorthb = 6                        !< 6 - index of free north border of the domain (south facing)
    INTEGER(iwp), PARAMETER                        ::  isouthb = 7                        !< 7 - index of north south border of the domain (north facing)
    INTEGER(iwp), PARAMETER                        ::  ieastb  = 8                        !< 8 - index of east border of the domain (west facing)
    INTEGER(iwp), PARAMETER                        ::  iwestb  = 9                        !< 9 - index of wast border of the domain (east facing)
    INTEGER(iwp), DIMENSION(0:9), PARAMETER        ::  idir = (/0,0,0,-1,1,0,0,0,-1,1/)   !< surface normal direction x indices
    INTEGER(iwp), DIMENSION(0:9), PARAMETER        ::  jdir = (/0,-1,1,0,0,0,-1,1,0,0/)   !< surface normal direction y indices
    INTEGER(iwp), DIMENSION(0:9), PARAMETER        ::  kdir = (/1,0,0,0,0,-1,0,0,0,0/)    !< surface normal direction z indices
    REAL(wp), DIMENSION(1:4)                       ::  ddxy2                              !< 1/dx^2 or 1/dy^2 (in surface normal direction)
    INTEGER(iwp), DIMENSION(1:4,6:9)               ::  ijdb                               !< start and end of the local domain border coordinates (set in code)
    LOGICAL, DIMENSION(6:9)                        ::  isborder                           !< is PE on the border of the domain in four corresponding directions
                                                                                          !< parameter but set in the code

!-- indices and sizes of urban surface model
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  surfl            !< coordinates of i-th local surface in local grid - surfl[:,k] = [d, z, y, x]
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  surf             !< coordinates of i-th surface in grid - surf[:,k] = [d, z, y, x]
    INTEGER(iwp)                                   ::  nsurfl           !< number of all surfaces in local processor
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  nsurfs           !< array of number of all surfaces in individual processors
    INTEGER(iwp)                                   ::  startsky         !< start index of block of sky
    INTEGER(iwp)                                   ::  endsky           !< end index of block of sky
    INTEGER(iwp)                                   ::  nskys            !< number of sky surfaces in local processor
    INTEGER(iwp)                                   ::  startland        !< start index of block of land and roof surfaces
    INTEGER(iwp)                                   ::  endland          !< end index of block of land and roof surfaces
    INTEGER(iwp)                                   ::  nlands           !< number of land and roof surfaces in local processor
    INTEGER(iwp)                                   ::  startwall        !< start index of block of wall surfaces
    INTEGER(iwp)                                   ::  endwall          !< end index of block of wall surfaces
    INTEGER(iwp)                                   ::  nwalls           !< number of wall surfaces in local processor
    INTEGER(iwp)                                   ::  startenergy      !< start index of block of real surfaces (land, walls and roofs)
    INTEGER(iwp)                                   ::  endenergy        !< end index of block of real surfaces (land, walls and roofs)
    INTEGER(iwp)                                   ::  nenergy          !< number of real surfaces in local processor
    INTEGER(iwp)                                   ::  nsurf            !< global number of surfaces in index array of surfaces (nsurf = Σproc nsurfs)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  surfstart        !< starts of blocks of surfaces for individual processors in array surf
                                                                        !< respective block for particular processor is surfstart[iproc]+1 : surfstart[iproc+1]
    INTEGER(iwp)                                   ::  nsvfl            !< number of svf for local processor
    INTEGER(iwp)                                   ::  ncsfl            !< no. of csf in local processor
                                                                        !< needed only during calc_svf but must be here because it is
                                                                        !< shared between subroutines usm_calc_svf and usm_raytrace

!-- type for calculation of svf
    TYPE t_svf
        INTEGER(iwp)                               :: isurflt           !< 
        INTEGER(iwp)                               :: isurfs            !< 
        REAL(wp)                                   :: rsvf              !< 
        REAL(wp)                                   :: rtransp           !< 
    END TYPE

!-- type for calculation of csf
    TYPE t_csf
        INTEGER(iwp)                               :: ip                !< 
        INTEGER(iwp)                               :: itx               !< 
        INTEGER(iwp)                               :: ity               !< 
        INTEGER(iwp)                               :: itz               !< 
        INTEGER(iwp)                               :: isurfs            !< 
        REAL(wp)                                   :: rsvf              !< 
        REAL(wp)                                   :: rtransp           !< 
    END TYPE

!-- arrays for calculation of svf and csf
    TYPE(t_svf), DIMENSION(:), POINTER             ::  asvf             !< pointer to growing svc array
    TYPE(t_csf), DIMENSION(:), POINTER             ::  acsf             !< pointer to growing csf array
    TYPE(t_svf), DIMENSION(:), ALLOCATABLE, TARGET ::  asvf1, asvf2     !< realizations of svf array
    TYPE(t_csf), DIMENSION(:), ALLOCATABLE, TARGET ::  acsf1, acsf2     !< realizations of csf array
    INTEGER(iwp)                                   ::  nsvfla           !< dimmension of array allocated for storage of svf in local processor
    INTEGER(iwp)                                   ::  ncsfla           !< dimmension of array allocated for storage of csf in local processor
    INTEGER(iwp)                                   ::  msvf, mcsf       !< mod for swapping the growing array
    INTEGER(iwp), PARAMETER                        ::  gasize = 10000   !< initial size of growing arrays
!-- temporary arrays for calculation of csf in raytracing
    INTEGER(iwp)                                   ::  maxboxesg        !< max number of boxes ray can cross in the domain
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  boxes            !< coordinates of gridboxes being crossed by ray
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  crlens           !< array of crossing lengths of ray for particular grid boxes
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  lad_ip           !< array of numbers of process where lad is stored 
    INTEGER(kind=MPI_ADDRESS_KIND), &
                  DIMENSION(:), ALLOCATABLE        ::  lad_disp         !< array of displaycements of lad in local array of proc lad_ip
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  lad_s_ray        !< array of received lad_s for appropriate gridboxes crossed by ray 

!-- arrays storing the values of USM
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  svfsurf          !< svfsurf[:,isvf] = index of source and target surface for svf[isvf]
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  svf              !< array of shape view factors+direct irradiation factors for local surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfins          !< array of sw radiation falling to local surface after i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinl          !< array of lw radiation for local surface after i-th reflection
    
                                                                        !< Inward radiation is also valid for virtual surfaces (radiation leaving domain)
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinsw         !< array of sw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlw         !< array of lw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdir      !< array of direct sw radiation falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdif      !< array of diffuse sw radiation from sky and model boundary falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlwdif      !< array of diffuse lw radiation from sky and model boundary falling to local surface
    
                                                                        !< Outward radiation is only valid for nonvirtual surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutsl        !< array of reflected sw radiation for local surface in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutll        !< array of reflected + emitted lw radiation for local surface in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfouts         !< array of reflected sw radiation for all surfaces in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutl         !< array of reflected + emitted lw radiation for all surfaces in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutsw        !< array of total sw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutlw        !< array of total lw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfhf           !< array of total radiation flux incoming to minus outgoing from local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  rad_net_l        !< local copy of rad_net (net radiation at surface)

!-- arrays for time averages
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  rad_net_av       !< average of rad_net_l
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinsw_av      !< average of sw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlw_av      !< average of lw radiation falling to local surface including radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdir_av   !< average of direct sw radiation falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswdif_av   !< average of diffuse sw radiation from sky and model boundary falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlwdif_av   !< average of diffuse lw radiation from sky and model boundary falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinswref_av   !< average of sw radiation falling to surface from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinlwref_av   !< average of lw radiation falling to surface from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutsw_av     !< average of total sw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfoutlw_av     !< average of total lw radiation outgoing from nonvirtual surfaces surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfins_av       !< average of array of residua of sw radiation absorbed in surface after last reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfinl_av       !< average of array of residua of lw radiation absorbed in surface after last reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  surfhf_av        !< average of total radiation flux incoming to minus outgoing from local surface    
    
!-- block variables needed for calculation of the plant canopy model inside the urban surface model
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  csfsurf          !< csfsurf[:,icsf] = index of target surface and csf grid index for csf[icsf]
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  csf              !< array of plant canopy sink fators + direct irradiation factors (transparency)
                                                                        !< for local surfaces
    INTEGER(wp), DIMENSION(:,:), ALLOCATABLE       ::  pcbl             !< k,j,i coordinates of l-th local plant canopy box pcbl[:,l] = [k, j, i]
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE    ::  gridpcbl         !< index of local pcb[k,j,i]
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinsw          !< array of absorbed sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  pcbinlw          !< array of absorbed lw radiation for local plant canopy box
    INTEGER(iwp)                                   ::  npcbl            !< number of the plant canopy gridboxes in local processor
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  pch              !< heights of the plant canopy
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE      ::  pct              !< top layer of the plant canopy
    REAL(wp), DIMENSION(:,:,:), POINTER            ::  usm_lad          !< subset of lad_s within urban surface, transformed to plain Z coordinate
    REAL(wp), DIMENSION(:), POINTER                ::  usm_lad_g        !< usm_lad globalized (used to avoid MPI RMA calls in raytracing)
    REAL(wp)                                       ::  prototype_lad    !< prototype leaf area density for computing effective optical depth
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        ::  nzterr, plantt   !< temporary global arrays for raytracing
    
!-- radiation related arrays (it should be better in interface of radiation module of PALM
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rad_sw_in_dir    !< direct sw radiation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rad_sw_in_diff   !< diffusion sw radiation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  rad_lw_in_diff   !< diffusion lw radiation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- anthropogenic heat sources
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          ::  aheat             !< daily average of anthropogenic heat (W/m2)
    REAL(wp), DIMENSION(:), ALLOCATABLE            ::  aheatprof         !< diurnal profile of anthropogenic heat 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- wall surface model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- wall surface model constants
    INTEGER(iwp), PARAMETER                        :: nzb_wall = 0       !< inner side of the wall model (to be switched)
    INTEGER(iwp), PARAMETER                        :: nzt_wall = 3       !< outer side of the wall model (to be switched)
    INTEGER(iwp), PARAMETER                        :: nzw = 4            !< number of wall layers (fixed for now)

    REAL(wp), DIMENSION(nzb_wall:nzt_wall)         :: zwn_default = (/0.0242_wp, 0.0969_wp, 0.346_wp, 1.0_wp /)
                                                                         !< normalized soil, wall and roof layer depths (m/m)
                                                                        
    REAL(wp)                                       ::   wall_inner_temperature = 296.0_wp    !< temperature of the inner wall surface (~23 degrees C) (K)
    REAL(wp)                                       ::   roof_inner_temperature = 296.0_wp    !< temperature of the inner roof surface (~23 degrees C) (K)
    REAL(wp)                                       ::   soil_inner_temperature = 283.0_wp    !< temperature of the deep soil (~10 degrees C) (K)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- surface and material model variables for walls, ground, roofs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        :: surface_types      !< array of types of wall parameters
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: zwn                !< normalized wall layer depths (m)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: ddz_wall           !< 1/dz_wall
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: ddz_wall_stag      !< 1/dz_wall_stag
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: dz_wall            !< wall grid spacing (center-center)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: dz_wall_stag       !< wall grid spacing (edge-edge) 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: zw                 !< wall layer depths (m)

#if defined( __nopointer )
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf             !< wall surface temperature (K)
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_p           !< progn. wall surface temperature (K)
#else
    REAL(wp), DIMENSION(:), POINTER                :: t_surf
    REAL(wp), DIMENSION(:), POINTER                :: t_surf_p 

    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_1
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_2
#endif
    REAL(wp), DIMENSION(:), ALLOCATABLE, TARGET    :: t_surf_av          !< average of wall surface temperature (K)

!-- Temporal tendencies for time stepping            
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: tt_surface_m       !< surface temperature tendency (K)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- Energy balance variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- parameters of the land, roof and wall surfaces
    LOGICAL,  DIMENSION(:), ALLOCATABLE            :: isroof_surf        !< is the surface the part of a roof
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: albedo_surf        !< albedo of the surface
!-- parameters of the wall surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: c_surface          !< heat capacity of the wall surface skin ( J m−2 K−1 )
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: emiss_surf         !< emissivity of the wall surface
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: lambda_surf        !< heat conductivity λS between air and surface ( W m−2 K−1 )
    
!-- parameters of the walls material
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: thickness_wall     !< thickness of the wall, roof and soil layers
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: rho_c_wall         !< volumetric heat capacity of the material ( J m-3 K-1 ) (= 2.19E6)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: lambda_h           !< heat conductivity λT of the material ( W m-1 K-1 )
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: roughness_wall     !< roughness relative to concrete
    
!-- output wall heat flux arrays
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: wshf               !< kinematic wall heat flux of sensible heat (needed for diffusion_s!<)
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: wshf_eb            !< wall heat flux of sensible heat in wall normal direction
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: wshf_eb_av         !< average of wshf_eb
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: wghf_eb            !< wall ground heat flux
    REAL(wp), DIMENSION(:), ALLOCATABLE            :: wghf_eb_av         !< average of wghf_eb

#if defined( __nopointer )
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET  :: t_wall             !< Wall temperature (K)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET  :: t_wall_av          !< Average of t_wall
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET  :: t_wall_p           !< Prog. wall temperature (K)
#else
    REAL(wp), DIMENSION(:,:), POINTER              :: t_wall, t_wall_p
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, TARGET  :: t_wall_av, t_wall_1, t_wall_2
#endif

!-- Wall temporal tendencies for time stepping
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: tt_wall_m          !< t_wall prognostic array 

!-- Surface and material parameters classes (surface_type)
!-- albedo, emissivity, lambda_surf, roughness, thickness, volumetric heat capacity, thermal conductivity
    INTEGER(iwp)                                   :: n_surface_types      !< number of the wall type categories
    INTEGER(iwp), PARAMETER                        :: n_surface_params = 8 !< number of parameters for each type of the wall
    INTEGER(iwp), PARAMETER                        :: ialbedo  = 1         !< albedo of the surface
    INTEGER(iwp), PARAMETER                        :: iemiss   = 2         !< emissivity of the surface
    INTEGER(iwp), PARAMETER                        :: ilambdas = 3         !< heat conductivity λS between air and surface ( W m−2 K−1 )
    INTEGER(iwp), PARAMETER                        :: irough   = 4         !< roughness relative to concrete
    INTEGER(iwp), PARAMETER                        :: icsurf   = 5         !< Surface skin layer heat capacity (J m−2 K−1 )
    INTEGER(iwp), PARAMETER                        :: ithick   = 6         !< thickness of the surface (wall, roof, land)  ( m )
    INTEGER(iwp), PARAMETER                        :: irhoC    = 7         !< volumetric heat capacity rho_ocean*C of the material ( J m−3 K−1 )
    INTEGER(iwp), PARAMETER                        :: ilambdah = 8         !< thermal conductivity λH of the wall (W m−1 K−1 )
    CHARACTER(12), DIMENSION(:), ALLOCATABLE       :: surface_type_names   !< names of wall types (used only for reports)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE        :: surface_type_codes   !< codes of wall types
    REAL(wp), DIMENSION(:,:), ALLOCATABLE          :: surface_params       !< parameters of wall types
    
    CHARACTER(len=*), PARAMETER                    :: svf_file_name='usm_svf'
    
!-- interfaces of subroutines accessed from outside of this module
    INTERFACE usm_check_data_output
       MODULE PROCEDURE usm_check_data_output
    END INTERFACE usm_check_data_output
    
    INTERFACE usm_check_parameters
       MODULE PROCEDURE usm_check_parameters
    END INTERFACE usm_check_parameters
    
    INTERFACE usm_data_output_3d
       MODULE PROCEDURE usm_data_output_3d
    END INTERFACE usm_data_output_3d
    
    INTERFACE usm_define_netcdf_grid
       MODULE PROCEDURE usm_define_netcdf_grid
    END INTERFACE usm_define_netcdf_grid

    INTERFACE usm_init_urban_surface
       MODULE PROCEDURE usm_init_urban_surface
    END INTERFACE usm_init_urban_surface

    INTERFACE usm_material_heat_model
       MODULE PROCEDURE usm_material_heat_model
    END INTERFACE usm_material_heat_model
    
    INTERFACE usm_parin
       MODULE PROCEDURE usm_parin
    END INTERFACE usm_parin

    INTERFACE usm_radiation
       MODULE PROCEDURE usm_radiation
    END INTERFACE usm_radiation
    
    INTERFACE usm_read_restart_data 
       MODULE PROCEDURE usm_read_restart_data
    END INTERFACE usm_read_restart_data

    INTERFACE usm_surface_energy_balance
       MODULE PROCEDURE usm_surface_energy_balance
    END INTERFACE usm_surface_energy_balance
    
    INTERFACE usm_swap_timelevel
       MODULE PROCEDURE usm_swap_timelevel
    END INTERFACE usm_swap_timelevel
    
    INTERFACE usm_wall_heat_flux
       MODULE PROCEDURE usm_wall_heat_flux
       MODULE PROCEDURE usm_wall_heat_flux_ij
    END INTERFACE usm_wall_heat_flux
    
    INTERFACE usm_write_restart_data
       MODULE PROCEDURE usm_write_restart_data
    END INTERFACE usm_write_restart_data
    
    SAVE

    PRIVATE
    
!-- Public parameters, constants and initial values
    PUBLIC split_diffusion_radiation,                                          &
           usm_anthropogenic_heat, usm_material_model, mrt_factors,            &
           usm_check_parameters,                                               &
           usm_energy_balance_land, usm_energy_balance_wall, nrefsteps,        &
           usm_init_urban_surface, usm_radiation, usm_read_restart_data,       &
           usm_wall_heat_flux,                                                 &
           usm_surface_energy_balance, usm_material_heat_model,                &
           usm_swap_timelevel, usm_check_data_output, usm_average_3d_data,     &
           usm_data_output_3d, usm_define_netcdf_grid, usm_parin,              &
           usm_write_restart_data,                                             &
           nzub, nzut, ra_horiz_coef, usm_lad_rma,                             &
           land_category, pedestrant_category, wall_category, roof_category,   &
           write_svf_on_init, read_svf_on_init


 CONTAINS

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine creates the necessary indices of the urban surfaces
!> and plant canopy and it allocates the needed arrays for USM
!------------------------------------------------------------------------------!
    SUBROUTINE usm_allocate_urban_surface
    
        IMPLICIT NONE
       
        INTEGER(iwp)                            :: i, j, k, d, l, ir, jr, ids
        INTEGER(iwp)                            :: nzubl, nzutl, isurf, ipcgb
        INTEGER(iwp)                            :: procid

        

        
!--     auxiliary vars
        ddxy2 = (/ddy2,ddy2,ddx2,ddx2/)      !< 1/dx^2 or 1/dy^2 (in surface normal direction)
        
        CALL location_message( '', .TRUE. )
        CALL location_message( '    allocation of needed arrays', .TRUE. )
!--     find nzub, nzut, nzu
        nzubl = minval(nzb_s_inner(nys:nyn,nxl:nxr))
        nzutl = maxval(nzb_s_inner(nys:nyn,nxl:nxr))
        nzubl = max(nzubl,nzb)
        
        IF ( plant_canopy )  THEN
!--         allocate needed arrays
            ALLOCATE( pct(nys:nyn,nxl:nxr) )
            ALLOCATE( pch(nys:nyn,nxl:nxr) )

!--         calculate plant canopy height
            npcbl = 0
            pct = 0.0_wp
            pch = 0.0_wp
            DO i = nxl, nxr
                DO j = nys, nyn
                    DO k = nzt+1, 0, -1
                        IF ( lad_s(k,j,i) /= 0.0_wp )  THEN
!--                         we are at the top of the pcs
                            pct(j,i) = k + nzb_s_inner(j,i)
                            pch(j,i) = k
                            npcbl = npcbl + pch(j,i)
                            EXIT
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
            
            nzutl = max(nzutl, maxval(pct))
!--         code of plant canopy model uses parameter pch_index
!--         we need to setup it here to right value
!--         (pch_index, lad_s and other arrays in PCM are defined flat)
            pch_index = maxval(pch)

            prototype_lad = maxval(lad_s) * .9_wp  !< better be *1.0 if lad is either 0 or maxval(lad) everywhere
            IF ( prototype_lad <= 0._wp ) prototype_lad = .3_wp
            !WRITE(message_string, '(a,f6.3)') 'Precomputing effective box optical ' &
            !    // 'depth using prototype leaf area density = ', prototype_lad
            !CALL message('usm_init_urban_surface', 'PA0520', 0, 0, -1, 6, 0)
        ENDIF
        
        nzutl = min(nzutl+nzut_free, nzt)
                  
#if defined( __parallel )
        CALL MPI_AllReduce(nzubl,nzub,1,MPI_INTEGER,MPI_MIN,comm2d,ierr);
        CALL MPI_AllReduce(nzutl,nzut,1,MPI_INTEGER,MPI_MAX,comm2d,ierr);
#else
        nzub = nzubl
        nzut = nzutl
#endif

!--     global number of urban layers
        nzu = nzut - nzub + 1
        
!--     allocate urban surfaces grid
!--     calc number of surfaces in local proc
        CALL location_message( '    calculation of indices for surfaces', .TRUE. )
        nsurfl = 0
!--     calculate land surface and roof
        startland = nsurfl+1
        nsurfl = nsurfl+(nxr-nxl+1)*(nyn-nys+1)
        endland = nsurfl
        nlands = endland-startland+1

!--     calculation of the walls
        startwall = nsurfl+1
        DO i = nxl, nxr
            DO j = nys, nyn
!--             test for walls
!--             (we don't use array flags because it isn't calculated in case of masking_method=.T.)
                DO ids = 1, 4  !-- four wall directions
                    jr = min(max(j-jdir(ids),0),ny)
                    ir = min(max(i-idir(ids),0),nx)
                    nsurfl = nsurfl + max(0, nzb_s_inner(jr,ir)-nzb_s_inner(j,i))
                ENDDO
            ENDDO
        ENDDO
        endwall = nsurfl
        nwalls = endwall-startwall+1
        
!--     range of energy balance surfaces
        nenergy = 0
        IF ( usm_energy_balance_land )  THEN
            startenergy = startland
            nenergy = nenergy + nlands
        ELSE
            startenergy = startwall
        ENDIF
        IF ( usm_energy_balance_wall )  THEN
            endenergy = endwall
            nenergy = nenergy + nwalls
        ELSE
            endenergy = endland
        ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     block of virtual surfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     calculate sky surfaces
        startsky = nsurfl+1
        nsurfl = nsurfl+(nxr-nxl+1)*(nyn-nys+1)
        endsky = nsurfl
        nskys = endsky-startsky+1
        
!--     border flags
#if defined( __parallel )
        isborder = (/ north_border_pe, south_border_pe, right_border_pe, left_border_pe /)
#else
        isborder = (/.TRUE.,.TRUE.,.TRUE.,.TRUE./)
#endif
!--     fill array of the limits of the local domain borders
        ijdb = RESHAPE( (/ nxl,nxr,nyn,nyn,nxl,nxr,nys,nys,nxr,nxr,nys,nyn,nxl,nxl,nys,nyn /), (/4, 4/) )
!--     calulation of the free borders of the domain
        DO ids = 6,9
            IF ( isborder(ids) )  THEN
!--             free border of the domain in direction ids
                DO i = ijdb(1,ids), ijdb(2,ids)
                    DO j = ijdb(3,ids), ijdb(4,ids)
                        k = nzut - max(nzb_s_inner(j,i), nzb_s_inner(j-jdir(ids),i-idir(ids)))
                        nsurfl = nsurfl + k
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
        
!--     fill gridpcbl and pcbl
        IF ( plant_canopy )  THEN
            ALLOCATE( pcbl(iz:ix, 1:npcbl) )
            ALLOCATE( gridpcbl(nzub:nzut,nys:nyn,nxl:nxr) )
            gridpcbl(:,:,:) = 0
            ipcgb = 0
            DO i = nxl, nxr
                DO j = nys, nyn
                    DO k = nzb_s_inner(j,i)+1, pct(j,i)
                        ipcgb = ipcgb + 1
                        gridpcbl(k,j,i) = ipcgb
                        pcbl(:,ipcgb) = (/ k, j, i /)
                    ENDDO
                ENDDO
            ENDDO

            ALLOCATE( pcbinsw( 1:npcbl ) )
            ALLOCATE( pcbinlw( 1:npcbl ) )
        ENDIF

!--     fill surfl
        ALLOCATE(surfl(4,nsurfl))
        isurf = 0
        
!--     add land surfaces or roofs
        DO i = nxl, nxr
            DO j = nys, nyn
                isurf = isurf + 1
                k = nzb_s_inner(j,i)+1
                surfl(:,isurf) = (/iroof,k,j,i/)
            ENDDO
        ENDDO

!--     add walls
        DO i = nxl, nxr
            DO j = nys, nyn
                DO ids = 1, 4  !> four wall directions
                    jr = min(max(j-jdir(ids),0),ny)
                    ir = min(max(i-idir(ids),0),nx)
                    DO k = nzb_s_inner(j,i)+1, nzb_s_inner(jr,ir)
                        isurf = isurf + 1
                        surfl(:,isurf) = (/ids,k,j,i/)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

!--     add sky
        DO i = nxl, nxr
            DO j = nys, nyn
                isurf = isurf + 1
                k = nzut
                surfl(:,isurf) = (/isky,k,j,i/)
            ENDDO
        ENDDO
        
!--     calulation of the free borders of the domain
        DO ids = 6,9
            IF ( isborder(ids) )  THEN
!--             free border of the domain in direction ids
                DO i = ijdb(1,ids), ijdb(2,ids)
                    DO j = ijdb(3,ids), ijdb(4,ids)
                        DO k = max(nzb_s_inner(j,i),nzb_s_inner(j-jdir(ids),i-idir(ids)))+1, nzut
                            isurf = isurf + 1
                            surfl(:,isurf) = (/ids,k,j,i/)
                        ENDDO
                    ENDDO
                ENDDO
            ENDIF
        ENDDO
        
!--     global array surf of indices of surfaces and displacement index array surfstart
        ALLOCATE(nsurfs(0:numprocs-1))
        
#if defined( __parallel )
        CALL MPI_Allgather(nsurfl,1,MPI_INTEGER,nsurfs,1,MPI_INTEGER,comm2d,ierr)
#else
        nsurfs(0) = nsurfl
#endif
        ALLOCATE(surfstart(0:numprocs))
        k = 0
        DO i=0,numprocs-1
            surfstart(i) = k
            k = k+nsurfs(i)
        ENDDO
        surfstart(numprocs) = k
        nsurf = k
        ALLOCATE(surf(4,nsurf))
        
#if defined( __parallel )
        CALL MPI_AllGatherv(surfl, nsurfl*4, MPI_INTEGER, surf, nsurfs*4, surfstart*4, MPI_INTEGER, comm2d, ierr)
#else
        surf = surfl
#endif
        
!--
!--     allocation of the arrays for direct and diffusion radiation
        CALL location_message( '    allocation of radiation arrays', .TRUE. )
!--     rad_sw_in, rad_lw_in are computed in radiation model,
!--     splitting of direct and diffusion part is done
!--     in usm_calc_diffusion_radiation for now
        ALLOCATE( rad_sw_in_dir(nysg:nyng,nxlg:nxrg) )
        ALLOCATE( rad_sw_in_diff(nysg:nyng,nxlg:nxrg) )
        ALLOCATE( rad_lw_in_diff(nysg:nyng,nxlg:nxrg) )
        
!--     allocate radiation arrays
        ALLOCATE( surfins(nsurfl) )
        ALLOCATE( surfinl(nsurfl) )
        ALLOCATE( surfinsw(nsurfl) )
        ALLOCATE( surfinlw(nsurfl) )
        ALLOCATE( surfinswdir(nsurfl) )
        ALLOCATE( surfinswdif(nsurfl) )
        ALLOCATE( surfinlwdif(nsurfl) )
        ALLOCATE( surfoutsl(startenergy:endenergy) )
        ALLOCATE( surfoutll(startenergy:endenergy) )
        ALLOCATE( surfoutsw(startenergy:endenergy) )
        ALLOCATE( surfoutlw(startenergy:endenergy) )
        ALLOCATE( surfouts(nsurf) ) !TODO: global surfaces without virtual
        ALLOCATE( surfoutl(nsurf) ) !TODO: global surfaces without virtual
        ALLOCATE( surfhf(startenergy:endenergy) )
        ALLOCATE( rad_net_l(startenergy:endenergy) )

!--     Wall surface model
!--     allocate arrays for wall surface model and define pointers
       
!--     allocate array of wall types and wall parameters
        ALLOCATE ( surface_types(startenergy:endenergy) )
       
!--     broadband albedo of the land, roof and wall surface
!--     for domain border and sky set artifically to 1.0
!--     what allows us to calculate heat flux leaving over
!--     side and top borders of the domain
        ALLOCATE ( albedo_surf(nsurfl) )
        albedo_surf = 1.0_wp
       
!--     wall and roof surface parameters
        ALLOCATE ( isroof_surf(startenergy:endenergy) )
        ALLOCATE ( emiss_surf(startenergy:endenergy) )
        ALLOCATE ( lambda_surf(startenergy:endenergy) )
        ALLOCATE ( c_surface(startenergy:endenergy) )
        ALLOCATE ( roughness_wall(startenergy:endenergy) )
       
!--     allocate wall and roof material parameters
        ALLOCATE ( thickness_wall(startenergy:endenergy) )
        ALLOCATE ( lambda_h(nzb_wall:nzt_wall,startenergy:endenergy) )
        ALLOCATE ( rho_c_wall(nzb_wall:nzt_wall,startenergy:endenergy) )

!--     allocate wall and roof layers sizes
        ALLOCATE ( zwn(nzb_wall:nzt_wall) )
        ALLOCATE ( dz_wall(nzb_wall:nzt_wall+1, startenergy:endenergy) )
        ALLOCATE ( ddz_wall(nzb_wall:nzt_wall+1, startenergy:endenergy) )
        ALLOCATE ( dz_wall_stag(nzb_wall:nzt_wall, startenergy:endenergy) )
        ALLOCATE ( ddz_wall_stag(nzb_wall:nzt_wall, startenergy:endenergy) )
        ALLOCATE ( zw(nzb_wall:nzt_wall, startenergy:endenergy) )

!--     allocate wall and roof temperature arrays
#if defined( __nopointer )
        ALLOCATE ( t_surf(startenergy:endenergy) )
        ALLOCATE ( t_surf_p(startenergy:endenergy) )
        ALLOCATE ( t_wall(nzb_wall:nzt_wall+1,startenergy:endenergy) )
        ALLOCATE ( t_wall_p(nzb_wall:nzt_wall+1,startenergy:endenergy) )
#else
        ALLOCATE ( t_surf_1(startenergy:endenergy) )
        ALLOCATE ( t_surf_2(startenergy:endenergy) )
        ALLOCATE ( t_wall_1(nzb_wall:nzt_wall+1,startenergy:endenergy) )
        ALLOCATE ( t_wall_2(nzb_wall:nzt_wall+1,startenergy:endenergy) )

!--     initial assignment of the pointers
        t_wall    => t_wall_1;    t_wall_p    => t_wall_2
        t_surf => t_surf_1; t_surf_p => t_surf_2
#endif

!--     allocate intermediate timestep arrays
        ALLOCATE ( tt_surface_m(startenergy:endenergy) )
        ALLOCATE ( tt_wall_m(nzb_wall:nzt_wall+1,startenergy:endenergy) )

!--     allocate wall heat flux output array
        ALLOCATE ( wshf(startwall:endwall) )
        ALLOCATE ( wshf_eb(startenergy:endenergy) )
        ALLOCATE ( wghf_eb(startenergy:endenergy) )

!--     set inital values for prognostic quantities
        tt_surface_m = 0.0_wp
        tt_wall_m    = 0.0_wp

        wshf = 0.0_wp
        wshf_eb = 0.0_wp
        wghf_eb = 0.0_wp
        
    END SUBROUTINE usm_allocate_urban_surface



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average urban surface output quantities as well as allocate
!> the array necessary for storing the average.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_average_3d_data( mode, variable )

        IMPLICIT NONE

        CHARACTER (len=*), INTENT(IN) ::  mode
        CHARACTER (len=*), INTENT(IN) :: variable
  
        INTEGER(iwp)                                       :: i, j, k, l, ids, iwl,istat
        CHARACTER (len=varnamelength)                      :: var, surfid
        INTEGER(iwp), PARAMETER                            :: nd = 5
        CHARACTER(len=6), DIMENSION(0:nd-1), PARAMETER     :: dirname = (/ '_roof ', '_south', '_north', '_west ', '_east ' /)

!--     find the real name of the variable
        var = TRIM(variable)
        DO i = 0, nd-1
            k = len(TRIM(var))
            j = len(TRIM(dirname(i)))
            IF ( var(k-j+1:k) == dirname(i) )  THEN
                ids = i
                var = var(:k-j)
                EXIT
            ENDIF
        ENDDO
        IF ( ids == -1 )  THEN
            var = TRIM(variable)
        ENDIF
        IF ( var(1:11) == 'usm_t_wall_'  .AND.  len(TRIM(var)) >= 12 )  THEN
!--          wall layers
            READ(var(12:12), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:10)
            ELSE
!--             wrong wall layer index
                RETURN
            ENDIF
        ENDIF

        IF ( mode == 'allocate' )  THEN
           
           SELECT CASE ( TRIM( var ) )
                
                CASE ( 'usm_rad_net' )
!--                 array of complete radiation balance
                    IF ( .NOT.  ALLOCATED(rad_net_av) )  THEN
                        ALLOCATE( rad_net_av(startenergy:endenergy) )
                        rad_net_av = 0.0_wp
                    ENDIF
                    
                CASE ( 'usm_rad_insw' )
!--                 array of sw radiation falling to surface after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfinsw_av) )  THEN
                        ALLOCATE( surfinsw_av(startenergy:endenergy) )
                        surfinsw_av = 0.0_wp
                    ENDIF
                                    
                CASE ( 'usm_rad_inlw' )
!--                 array of lw radiation falling to surface after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfinlw_av) )  THEN
                        ALLOCATE( surfinlw_av(startenergy:endenergy) )
                        surfinlw_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inswdir' )
!--                 array of direct sw radiation falling to surface from sun
                    IF ( .NOT.  ALLOCATED(surfinswdir_av) )  THEN
                        ALLOCATE( surfinswdir_av(startenergy:endenergy) )
                        surfinswdir_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inswdif' )
!--                 array of difusion sw radiation falling to surface from sky and borders of the domain
                    IF ( .NOT.  ALLOCATED(surfinswdif_av) )  THEN
                        ALLOCATE( surfinswdif_av(startenergy:endenergy) )
                        surfinswdif_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inswref' )
!--                 array of sw radiation falling to surface from reflections
                    IF ( .NOT.  ALLOCATED(surfinswref_av) )  THEN
                        ALLOCATE( surfinswref_av(startenergy:endenergy) )
                        surfinswref_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inlwdif' )
!--                 array of sw radiation falling to surface after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfinlwdif_av) )  THEN
                        ALLOCATE( surfinlwdif_av(startenergy:endenergy) )
                        surfinlwdif_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_inlwref' )
!--                 array of lw radiation falling to surface from reflections
                    IF ( .NOT.  ALLOCATED(surfinlwref_av) )  THEN
                        ALLOCATE( surfinlwref_av(startenergy:endenergy) )
                        surfinlwref_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_outsw' )
!--                 array of sw radiation emitted from surface after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfoutsw_av) )  THEN
                        ALLOCATE( surfoutsw_av(startenergy:endenergy) )
                        surfoutsw_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_outlw' )
!--                 array of lw radiation emitted from surface after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfoutlw_av) )  THEN
                        ALLOCATE( surfoutlw_av(startenergy:endenergy) )
                        surfoutlw_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_rad_ressw' )
!--                 array of residua of sw radiation absorbed in surface after last reflection
                    IF ( .NOT.  ALLOCATED(surfins_av) )  THEN
                        ALLOCATE( surfins_av(startenergy:endenergy) )
                        surfins_av = 0.0_wp
                    ENDIF
                                    
                CASE ( 'usm_rad_reslw' )
!--                 array of residua of lw radiation absorbed in surface after last reflection
                    IF ( .NOT.  ALLOCATED(surfinl_av) )  THEN
                        ALLOCATE( surfinl_av(startenergy:endenergy) )
                        surfinl_av = 0.0_wp
                    ENDIF
                                    
                CASE ( 'usm_rad_hf' )
!--                 array of heat flux from radiation for surfaces after i-th reflection
                    IF ( .NOT.  ALLOCATED(surfhf_av) )  THEN
                        ALLOCATE( surfhf_av(startenergy:endenergy) )
                        surfhf_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_wshf' )
!--                 array of sensible heat flux from surfaces
!--                 land surfaces
                    IF ( .NOT.  ALLOCATED(wshf_eb_av) )  THEN
                        ALLOCATE( wshf_eb_av(startenergy:endenergy) )
                        wshf_eb_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_wghf' )
!--                 array of heat flux from ground (wall, roof, land)
                    IF ( .NOT.  ALLOCATED(wghf_eb_av) )  THEN
                        ALLOCATE( wghf_eb_av(startenergy:endenergy) )
                        wghf_eb_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_t_surf' )
!--                 surface temperature for surfaces
                    IF ( .NOT.  ALLOCATED(t_surf_av) )  THEN
                        ALLOCATE( t_surf_av(startenergy:endenergy) )
                        t_surf_av = 0.0_wp
                    ENDIF

                CASE ( 'usm_t_wall' )
!--                 wall temperature for iwl layer of walls and land
                    IF ( .NOT.  ALLOCATED(t_wall_av) )  THEN
                        ALLOCATE( t_wall_av(nzb_wall:nzt_wall,startenergy:endenergy) )
                        t_wall_av = 0.0_wp
                    ENDIF

               CASE DEFAULT
                   CONTINUE

           END SELECT

        ELSEIF ( mode == 'sum' )  THEN
           
           SELECT CASE ( TRIM( var ) )
                
                CASE ( 'usm_rad_net' )
!--                 array of complete radiation balance
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            rad_net_av(l) = rad_net_av(l) + rad_net_l(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_insw' )
!--                 array of sw radiation falling to surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinsw_av(l) = surfinsw_av(l) + surfinsw(l)
                        ENDIF
                    ENDDO
                             
                CASE ( 'usm_rad_inlw' )
!--                 array of lw radiation falling to surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinlw_av(l) = surfinlw_av(l) + surfinlw(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inswdir' )
!--                 array of direct sw radiation falling to surface from sun
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinswdir_av(l) = surfinswdir_av(l) + surfinswdir(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inswdif' )
!--                 array of difusion sw radiation falling to surface from sky and borders of the domain
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinswdif_av(l) = surfinswdif_av(l) + surfinswdif(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inswref' )
!--                 array of sw radiation falling to surface from reflections
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinswref_av(l) = surfinswref_av(l) + surfinsw(l) - &
                                                surfinswdir(l) - surfinswdif(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inlwdif' )
!--                 array of sw radiation falling to surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinlwdif_av(l) = surfinlwdif_av(l) + surfinlwdif(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_inlwref' )
!--                 array of lw radiation falling to surface from reflections
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinlwref_av(l) = surfinlwref_av(l) + &
                                                surfinlw(l) - surfinlwdif(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_outsw' )
!--                 array of sw radiation emitted from surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfoutsw_av(l) = surfoutsw_av(l) + surfoutsw(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_outlw' )
!--                 array of lw radiation emitted from surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfoutlw_av(l) = surfoutlw_av(l) + surfoutlw(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_ressw' )
!--                 array of residua of sw radiation absorbed in surface after last reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfins_av(l) = surfins_av(l) + surfins(l)
                        ENDIF
                    ENDDO
                                    
                CASE ( 'usm_rad_reslw' )
!--                 array of residua of lw radiation absorbed in surface after last reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinl_av(l) = surfinl_av(l) + surfinl(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_hf' )
!--                 array of heat flux from radiation for surfaces after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfhf_av(l) = surfhf_av(l) + surfhf(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_wshf' )
!--                 array of sensible heat flux from surfaces (land, roof, wall)
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            wshf_eb_av(l) = wshf_eb_av(l) + wshf_eb(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_wghf' )
!--                 array of heat flux from ground (wall, roof, land)
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            wghf_eb_av(l) = wghf_eb_av(l) + wghf_eb(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_t_surf' )
!--                 surface temperature for surfaces
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            t_surf_av(l) = t_surf_av(l) + t_surf(l)
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_t_wall' )
!--                 wall temperature for  iwl layer of walls and land
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            t_wall_av(iwl, l) = t_wall_av(iwl,l) + t_wall(iwl, l)
                        ENDIF
                    ENDDO
                    
                CASE DEFAULT
                    CONTINUE

           END SELECT

        ELSEIF ( mode == 'average' )  THEN
           
           SELECT CASE ( TRIM( var ) )
                
                CASE ( 'usm_rad_net' )
!--                 array of complete radiation balance
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            rad_net_av(l) = rad_net_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_insw' )
!--                 array of sw radiation falling to surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinsw_av(l) = surfinsw_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                                    
                CASE ( 'usm_rad_inlw' )
!--                 array of lw radiation falling to surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinlw_av(l) = surfinlw_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_rad_inswdir' )
!--                 array of direct sw radiation falling to surface from sun
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinswdir_av(l) = surfinswdir_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_rad_inswdif' )
!--                 array of difusion sw radiation falling to surface from sky and borders of the domain
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinswdif_av(l) = surfinswdif_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_rad_inswref' )
!--                 array of sw radiation falling to surface from reflections
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinswref_av(l) = surfinswref_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_rad_inlwdif' )
!--                 array of sw radiation falling to surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinlwdif_av(l) = surfinlwdif_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_rad_inlwref' )
!--                 array of lw radiation falling to surface from reflections
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinlwref_av(l) = surfinlwref_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_rad_outsw' )
!--                 array of sw radiation emitted from surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfoutsw_av(l) = surfoutsw_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_rad_outlw' )
!--                 array of lw radiation emitted from surface after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfoutlw_av(l) = surfoutlw_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_rad_ressw' )
!--                 array of residua of sw radiation absorbed in surface after last reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfins_av(l) = surfins_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                                    
                CASE ( 'usm_rad_reslw' )
!--                 array of residua of lw radiation absorbed in surface after last reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfinl_av(l) = surfinl_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO
                    
                CASE ( 'usm_rad_hf' )
!--                 array of heat flux from radiation for surfaces after i-th reflection
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            surfhf_av(l) = surfhf_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_wshf' )
!--                 array of sensible heat flux from surfaces (land, roof, wall)
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            wshf_eb_av(l) = wshf_eb_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_wghf' )
!--                 array of heat flux from ground (wall, roof, land)
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            wghf_eb_av(l) = wghf_eb_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_t_surf' )
!--                 surface temperature for surfaces
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            t_surf_av(l) = t_surf_av(l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

                CASE ( 'usm_t_wall' )
!--                 wall temperature for  iwl layer of walls and land
                    DO l = startenergy, endenergy
                        IF ( surfl(id,l) == ids )  THEN
                            t_wall_av(iwl, l) = t_wall_av(iwl,l) / REAL( average_count_3d, kind=wp )
                        ENDIF
                    ENDDO

           END SELECT

        ENDIF

    END SUBROUTINE usm_average_3d_data


!------------------------------------------------------------------------------!
!> Calculates radiation absorbed by box with given size and LAD.
!>
!> Simulates resol**2 rays (by equally spacing a bounding horizontal square
!> conatining all possible rays that would cross the box) and calculates
!> average transparency per ray. Returns fraction of absorbed radiation flux
!> and area for which this fraction is effective.
!------------------------------------------------------------------------------!
    PURE SUBROUTINE usm_box_absorb(boxsize, resol, dens, uvec, area, absorb)
        IMPLICIT NONE

        REAL(wp), DIMENSION(3), INTENT(in) :: &
            boxsize, &      !< z, y, x size of box in m
            uvec            !< z, y, x unit vector of incoming flux
        INTEGER(iwp), INTENT(in) :: &
            resol           !< No. of rays in x and y dimensions
        REAL(wp), INTENT(in) :: &
            dens            !< box density (e.g. Leaf Area Density)
        REAL(wp), INTENT(out) :: &
            area, &         !< horizontal area for flux absorbtion
            absorb          !< fraction of absorbed flux
        REAL(wp) :: &
            xshift, yshift, &
            xmin, xmax, ymin, ymax, &
            xorig, yorig, &
            dx1, dy1, dz1, dx2, dy2, dz2, &
            crdist, &
            transp
        INTEGER(iwp) :: &
            i, j

        xshift = uvec(3) / uvec(1) * boxsize(1)
        xmin = min(0._wp, -xshift)
        xmax = boxsize(3) + max(0._wp, -xshift)
        yshift = uvec(2) / uvec(1) * boxsize(1)
        ymin = min(0._wp, -yshift)
        ymax = boxsize(2) + max(0._wp, -yshift)

        transp = 0._wp
        DO i = 1, resol
            xorig = xmin + (xmax-xmin) * (i-.5_wp) / resol
            DO j = 1, resol
                yorig = ymin + (ymax-ymin) * (j-.5_wp) / resol

                dz1 = 0._wp
                dz2 = boxsize(1)/uvec(1)

                IF ( uvec(2) > 0._wp )  THEN
                    dy1 = -yorig             / uvec(2) !< crossing with y=0
                    dy2 = (boxsize(2)-yorig) / uvec(2) !< crossing with y=boxsize(2)
                ELSE IF ( uvec(2) < 0._wp )  THEN
                    dy1 = (boxsize(2)-yorig) / uvec(2) !< crossing with y=boxsize(2)
                    dy2 = -yorig             / uvec(2) !< crossing with y=0
                ELSE !uvec(2)==0
                    dy1 = -huge(1._wp)
                    dy2 = huge(1._wp)
                ENDIF

                IF ( uvec(3) > 0._wp )  THEN
                    dx1 = -xorig             / uvec(3) !< crossing with x=0
                    dx2 = (boxsize(3)-xorig) / uvec(3) !< crossing with x=boxsize(3)
                ELSE IF ( uvec(3) < 0._wp )  THEN
                    dx1 = (boxsize(3)-xorig) / uvec(3) !< crossing with x=boxsize(3)
                    dx2 = -xorig             / uvec(3) !< crossing with x=0
                ELSE !uvec(1)==0
                    dx1 = -huge(1._wp)
                    dx2 = huge(1._wp)
                ENDIF

                crdist = max(0._wp, (min(dz2, dy2, dx2) - max(dz1, dy1, dx1)))
                transp = transp + exp(-ext_coef * dens * crdist)
            ENDDO
        ENDDO
        transp = transp / resol**2
        area = (boxsize(3)+xshift)*(boxsize(2)+yshift)
        absorb = 1._wp - transp
        
    END SUBROUTINE usm_box_absorb
    
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine splits direct and diffusion dw radiation
!> It sould not be called in case the radiation model already does it
!> It follows <CITATION>
!------------------------------------------------------------------------------!
    SUBROUTINE usm_calc_diffusion_radiation 
    
        REAL(wp), PARAMETER                          ::  sol_const = 1367.0_wp   !< solar conbstant
        REAL(wp), PARAMETER                          :: lowest_solarUp = 0.1_wp  !< limit the sun elevation to protect stability of the calculation
        INTEGER(iwp)                                 :: i, j
        REAL(wp), PARAMETER                          ::  year_seconds = 86400._wp * 365._wp
        REAL(wp)                                     ::  year_angle              !< angle
        REAL(wp)                                     ::  etr                     !< extraterestrial radiation
        REAL(wp)                                     ::  corrected_solarUp       !< corrected solar up radiation
        REAL(wp)                                     ::  horizontalETR           !< horizontal extraterestrial radiation
        REAL(wp)                                     ::  clearnessIndex          !< clearness index
        REAL(wp)                                     ::  diff_frac               !< diffusion fraction of the radiation

        
!--     Calculate current day and time based on the initial values and simulation time
        year_angle = ((day_init*86400) + time_utc_init+time_since_reference_point) &
                       / year_seconds * 2.0_wp * pi
        
        etr = sol_const * (1.00011_wp +                                            &
                          0.034221_wp * cos(year_angle) +                          &
                          0.001280_wp * sin(year_angle) +                          &
                          0.000719_wp * cos(2.0_wp * year_angle) +                 &
                          0.000077_wp * sin(2.0_wp * year_angle))
        
!--    
!--     Under a very low angle, we keep extraterestrial radiation at
!--     the last small value, therefore the clearness index will be pushed
!--     towards 0 while keeping full continuity.
!--    
        IF ( zenith(0) <= lowest_solarUp )  THEN
            corrected_solarUp = lowest_solarUp
        ELSE
            corrected_solarUp = zenith(0)
        ENDIF
        
        horizontalETR = etr * corrected_solarUp
        
        DO i = nxlg, nxrg
            DO j = nysg, nyng
                clearnessIndex = rad_sw_in(0,j,i) / horizontalETR
                diff_frac = 1.0_wp / (1.0_wp + exp(-5.0033_wp + 8.6025_wp * clearnessIndex))
                rad_sw_in_diff(j,i) = rad_sw_in(0,j,i) * diff_frac
                rad_sw_in_dir(j,i)  = rad_sw_in(0,j,i) * (1.0_wp - diff_frac)
                rad_lw_in_diff(j,i) = rad_lw_in(0,j,i)
            ENDDO
        ENDDO
        
    END SUBROUTINE usm_calc_diffusion_radiation
    

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates shape view factors SVF and plant sink canopy factors PSCF
!> !!!!!DESCRIPTION!!!!!!!!!!
!------------------------------------------------------------------------------!
    SUBROUTINE usm_calc_svf
    
        IMPLICIT NONE
        
        INTEGER(iwp)                                :: i, j, k, l, d, ip, jp
        INTEGER(iwp)                                :: isvf, ksvf, icsf, kcsf, npcsfl, isvf_surflt, imrtt, imrtf
        INTEGER(iwp)                                :: sd, td, ioln, iproc
        REAL(wp),     DIMENSION(0:9)                :: facearea
        INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE   :: nzterrl, planthl
        REAL(wp),     DIMENSION(:,:), ALLOCATABLE   :: csflt, pcsflt
        INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE   :: kcsflt,kpcsflt
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE     :: icsflt,dcsflt,ipcsflt,dpcsflt
        REAL(wp), DIMENSION(3)                      :: uv
        LOGICAL                                     :: visible
        REAL(wp), DIMENSION(3)                      :: sa, ta          !< real coordinates z,y,x of source and target
        REAL(wp)                                    :: transparency, rirrf, sqdist, svfsum
        INTEGER(iwp)                                :: isurflt, isurfs, isurflt_prev
        INTEGER(iwp)                                :: itx, ity, itz
        CHARACTER(len=7)                            :: pid_char = ''
        INTEGER(iwp)                                :: win_lad, minfo
        REAL(wp), DIMENSION(:,:,:), POINTER         :: lad_s_rma       !< fortran pointer, but lower bounds are 1
        TYPE(c_ptr)                                 :: lad_s_rma_p     !< allocated c pointer
        INTEGER(kind=MPI_ADDRESS_KIND)              :: size_lad_rma
    
!--     calculation of the SVF
        CALL location_message( '    calculation of SVF and CSF', .TRUE. )

!--     precalculate face areas for different face directions using normal vector
        DO d = 0, 9
            facearea(d) = 1._wp
            IF ( idir(d) == 0 ) facearea(d) = facearea(d) * dx
            IF ( jdir(d) == 0 ) facearea(d) = facearea(d) * dy
            IF ( kdir(d) == 0 ) facearea(d) = facearea(d) * dz
        ENDDO

!--     initialize variables and temporary arrays for calculation of svf and csf
        nsvfl  = 0
        ncsfl  = 0
        nsvfla = gasize
        msvf   = 1
        ALLOCATE( asvf1(nsvfla) )
        asvf => asvf1
        IF ( plant_canopy )  THEN
            ncsfla = gasize
            mcsf   = 1
            ALLOCATE( acsf1(ncsfla) )
            acsf => acsf1
        ENDIF
        
!--     initialize temporary terrain and plant canopy height arrays (global 2D array!)
        ALLOCATE( nzterr(0:(nx+1)*(ny+1)-1) )
#if defined( __parallel )
        ALLOCATE( nzterrl(nys:nyn,nxl:nxr) )
        nzterrl = nzb_s_inner(nys:nyn,nxl:nxr)
        CALL MPI_AllGather( nzterrl, nnx*nny, MPI_INTEGER, &
                            nzterr, nnx*nny, MPI_INTEGER, comm2d, ierr )
        DEALLOCATE(nzterrl)
#else
        nzterr = RESHAPE( nzb_s_inner(nys:nyn,nxl:nxr), (/(nx+1)*(ny+1)/) )
#endif
        IF ( plant_canopy )  THEN
            ALLOCATE( plantt(0:(nx+1)*(ny+1)-1) )
            maxboxesg = nx + ny + nzu + 1
!--         temporary arrays storing values for csf calculation during raytracing
            ALLOCATE( boxes(3, maxboxesg) )
            ALLOCATE( crlens(maxboxesg) )

#if defined( __parallel )
            ALLOCATE( planthl(nys:nyn,nxl:nxr) )
            planthl = pch(nys:nyn,nxl:nxr)
        
            CALL MPI_AllGather( planthl, nnx*nny, MPI_INTEGER, &
                                plantt, nnx*nny, MPI_INTEGER, comm2d, ierr )
            DEALLOCATE( planthl )
            
!--         temporary arrays storing values for csf calculation during raytracing
            ALLOCATE( lad_ip(maxboxesg) )
            ALLOCATE( lad_disp(maxboxesg) )

            IF ( usm_lad_rma )  THEN
                ALLOCATE( lad_s_ray(maxboxesg) )
                
                ! set conditions for RMA communication
                CALL MPI_Info_create(minfo, ierr)
                CALL MPI_Info_set(minfo, 'accumulate_ordering', '', ierr)
                CALL MPI_Info_set(minfo, 'accumulate_ops', 'same_op', ierr)
                CALL MPI_Info_set(minfo, 'same_size', 'true', ierr)
                CALL MPI_Info_set(minfo, 'same_disp_unit', 'true', ierr)

!--             Allocate and initialize the MPI RMA window
!--             must be in accordance with allocation of lad_s in plant_canopy_model
!--             optimization of memory should be done
!--             Argument X of function c_sizeof(X) needs arbitrary REAL(wp) value, set to 1.0_wp for now
                size_lad_rma = c_sizeof(1.0_wp)*nnx*nny*nzu
                CALL MPI_Win_allocate(size_lad_rma, c_sizeof(1.0_wp), minfo, comm2d, &
                                        lad_s_rma_p, win_lad, ierr)
                CALL c_f_pointer(lad_s_rma_p, lad_s_rma, (/ nzu, nny, nnx /))
                usm_lad(nzub:, nys:, nxl:) => lad_s_rma(:,:,:)
            ELSE
                ALLOCATE(usm_lad(nzub:nzut, nys:nyn, nxl:nxr))
            ENDIF
#else
            plantt = RESHAPE( pct(nys:nyn,nxl:nxr), (/(nx+1)*(ny+1)/) )
            ALLOCATE(usm_lad(nzub:nzut, nys:nyn, nxl:nxr))
#endif
            usm_lad(:,:,:) = 0._wp
            DO i = nxl, nxr
                DO j = nys, nyn
                    k = nzb_s_inner(j, i)
                    usm_lad(k:nzut, j, i) = lad_s(0:nzut-k, j, i)
                ENDDO
            ENDDO

#if defined( __parallel )
            IF ( usm_lad_rma )  THEN
                CALL MPI_Info_free(minfo, ierr)
                CALL MPI_Win_lock_all(0, win_lad, ierr)
            ELSE
                ALLOCATE( usm_lad_g(0:(nx+1)*(ny+1)*nzu-1) )
                CALL MPI_AllGather( usm_lad, nnx*nny*nzu, MPI_REAL, &
                                    usm_lad_g, nnx*nny*nzu, MPI_REAL, comm2d, ierr )
            ENDIF
#endif
        ENDIF

        IF ( mrt_factors )  THEN
            OPEN(153, file='MRT_TARGETS', access='SEQUENTIAL', &
                    action='READ', status='OLD', form='FORMATTED', err=524)
            OPEN(154, file='MRT_FACTORS'//myid_char, access='DIRECT', recl=(5*4+2*8), &
                    action='WRITE', status='REPLACE', form='UNFORMATTED', err=525)
            imrtf = 1
            DO
                READ(153, *, end=526, err=524) imrtt, i, j, k
                IF ( i < nxl  .OR.  i > nxr &
                     .OR.  j < nys  .OR.  j > nyn ) CYCLE
                ta = (/ REAL(k), REAL(j), REAL(i) /)

                DO isurfs = 1, nsurf
                    IF ( .NOT.  usm_facing(i, j, k, -1, &
                        surf(ix, isurfs), surf(iy, isurfs), &
                        surf(iz, isurfs), surf(id, isurfs)) )  THEN
                        CYCLE
                    ENDIF
                      
                    sd = surf(id, isurfs)
                    sa = (/ REAL(surf(iz, isurfs), wp) - 0.5_wp * kdir(sd), &
                            REAL(surf(iy, isurfs), wp) - 0.5_wp * jdir(sd), &
                            REAL(surf(ix, isurfs), wp) - 0.5_wp * idir(sd) /)

!--                 unit vector source -> target
                    uv = (/ (ta(1)-sa(1))*dz, (ta(2)-sa(2))*dy, (ta(3)-sa(3))*dx /)
                    sqdist = SUM(uv(:)**2)
                    uv = uv / SQRT(sqdist)

!--                 irradiance factor - see svf. Here we consider that target face is always normal,
!--                 i.e. the second dot product equals 1
                    rirrf = dot_product((/ kdir(sd), jdir(sd), idir(sd) /), uv) &
                        / (pi * sqdist) * facearea(sd)

!--                 raytrace while not creating any canopy sink factors
                    CALL usm_raytrace(sa, ta, isurfs, rirrf, 1._wp, .FALSE., &
                            visible, transparency, win_lad)
                    IF ( .NOT.  visible ) CYCLE

                    !rsvf = rirrf * transparency
                    WRITE(154, rec=imrtf, err=525) INT(imrtt, kind=4), &
                        INT(surf(id, isurfs), kind=4), &
                        INT(surf(iz, isurfs), kind=4), &
                        INT(surf(iy, isurfs), kind=4), &
                        INT(surf(ix, isurfs), kind=4), &
                        REAL(rirrf, kind=8), REAL(transparency, kind=8)
                    imrtf = imrtf + 1

                ENDDO !< isurfs
            ENDDO !< MRT_TARGETS record

524         message_string = 'error reading file MRT_TARGETS'
            CALL message( 'usm_calc_svf', 'PA0524', 1, 2, 0, 6, 0 )

525         message_string = 'error writing file MRT_FACTORS'//myid_char
            CALL message( 'usm_calc_svf', 'PA0525', 1, 2, 0, 6, 0 )

526         CLOSE(153)
            CLOSE(154)
        ENDIF  !< mrt_factors

        
        DO isurflt = 1, nsurfl
!--         determine face centers
            td = surfl(id, isurflt)
            IF ( td >= isky  .AND.  .NOT.  plant_canopy ) CYCLE
            ta = (/ REAL(surfl(iz, isurflt), wp) - 0.5_wp * kdir(td),  &
                      REAL(surfl(iy, isurflt), wp) - 0.5_wp * jdir(td),  &
                      REAL(surfl(ix, isurflt), wp) - 0.5_wp * idir(td)  /)
            DO isurfs = 1, nsurf
                IF ( .NOT.  usm_facing(surfl(ix, isurflt), surfl(iy, isurflt), &
                    surfl(iz, isurflt), surfl(id, isurflt), &
                    surf(ix, isurfs), surf(iy, isurfs), &
                    surf(iz, isurfs), surf(id, isurfs)) )  THEN
                    CYCLE
                ENDIF
                  
                sd = surf(id, isurfs)
                sa = (/ REAL(surf(iz, isurfs), wp) - 0.5_wp * kdir(sd),  &
                        REAL(surf(iy, isurfs), wp) - 0.5_wp * jdir(sd),  &
                        REAL(surf(ix, isurfs), wp) - 0.5_wp * idir(sd)  /)

!--             unit vector source -> target
                uv = (/ (ta(1)-sa(1))*dz, (ta(2)-sa(2))*dy, (ta(3)-sa(3))*dx /)
                sqdist = SUM(uv(:)**2)
                uv = uv / SQRT(sqdist)
                
!--             irradiance factor (our unshaded shape view factor) = view factor per differential target area * source area
                rirrf = dot_product((/ kdir(sd), jdir(sd), idir(sd) /), uv) & ! cosine of source normal and direction
                    * dot_product((/ kdir(td), jdir(td), idir(td) /), -uv) &  ! cosine of target normal and reverse direction
                    / (pi * sqdist) & ! square of distance between centers
                    * facearea(sd)

!--             raytrace + process plant canopy sinks within
                CALL usm_raytrace(sa, ta, isurfs, rirrf, facearea(td), .TRUE., &
                        visible, transparency, win_lad)
                
                IF ( .NOT.  visible ) CYCLE
                IF ( td >= isky ) CYCLE !< we calculated these only for raytracing
                                        !< to find plant canopy sinks, we don't need svf for them
                ! rsvf = rirrf * transparency

!--             write to the svf array
                nsvfl = nsvfl + 1
!--             check dimmension of asvf array and enlarge it if needed
                IF ( nsvfla < nsvfl )  THEN
                    k = nsvfla * 2
                    IF ( msvf == 0 )  THEN
                        msvf = 1
                        ALLOCATE( asvf1(k) )
                        asvf => asvf1
                        asvf1(1:nsvfla) = asvf2
                        DEALLOCATE( asvf2 )
                    ELSE
                        msvf = 0
                        ALLOCATE( asvf2(k) )
                        asvf => asvf2
                        asvf2(1:nsvfla) = asvf1
                        DEALLOCATE( asvf1 )
                    ENDIF
                    nsvfla = k
                ENDIF
!--             write svf values into the array
                asvf(nsvfl)%isurflt = isurflt
                asvf(nsvfl)%isurfs = isurfs
                asvf(nsvfl)%rsvf = rirrf !we postopne multiplication by transparency
                asvf(nsvfl)%rtransp = transparency !a.k.a. Direct Irradiance Factor
            ENDDO
        ENDDO

        CALL location_message( '    waiting for completion of SVF and CSF calculation in all processes', .TRUE. )
!--     deallocate temporary global arrays
        DEALLOCATE(nzterr)
        
        IF ( plant_canopy )  THEN
!--         finalize mpi_rma communication and deallocate temporary arrays
#if defined( __parallel )
            IF ( usm_lad_rma )  THEN
                CALL MPI_Win_flush_all(win_lad, ierr)
!--             unlock MPI window
                CALL MPI_Win_unlock_all(win_lad, ierr)
!--             free MPI window
                CALL MPI_Win_free(win_lad, ierr)
                
!--             deallocate temporary arrays storing values for csf calculation during raytracing
                DEALLOCATE( lad_s_ray )
!--             usm_lad is the pointer to lad_s_rma in case of usm_lad_rma
!--             and must not be deallocated here
            ELSE
                DEALLOCATE(usm_lad)
                DEALLOCATE(usm_lad_g)
            ENDIF
#else
            DEALLOCATE(usm_lad)
#endif
            DEALLOCATE( boxes )
            DEALLOCATE( crlens )
            DEALLOCATE( plantt )
        ENDIF

        CALL location_message( '    calculation of the complete SVF array', .TRUE. )

!--     sort svf ( a version of quicksort )
        CALL quicksort_svf(asvf,1,nsvfl)

        ALLOCATE( svf(ndsvf,nsvfl) )
        ALLOCATE( svfsurf(idsvf,nsvfl) )

        !< load svf from the structure array to plain arrays
        isurflt_prev = -1
        ksvf = 1
        svfsum = 0._wp
        DO isvf = 1, nsvfl
!--         normalize svf per target face
            IF ( asvf(ksvf)%isurflt /= isurflt_prev )  THEN
                IF ( isurflt_prev /= -1  .AND.  svfsum /= 0._wp )  THEN
!--                 TODO detect and log when normalization differs too much from 1
                    svf(1, isvf_surflt:isvf-1) = svf(1, isvf_surflt:isvf-1) / svfsum
                ENDIF
                isurflt_prev = asvf(ksvf)%isurflt
                isvf_surflt = isvf
                svfsum = asvf(ksvf)%rsvf !?? / asvf(ksvf)%rtransp
            ELSE
                svfsum = svfsum + asvf(ksvf)%rsvf !?? / asvf(ksvf)%rtransp
            ENDIF

            svf(:, isvf) = (/ asvf(ksvf)%rsvf, asvf(ksvf)%rtransp /)
            svfsurf(:, isvf) = (/ asvf(ksvf)%isurflt, asvf(ksvf)%isurfs /)

!--         next element
            ksvf = ksvf + 1
        ENDDO

        IF ( isurflt_prev /= -1  .AND.  svfsum /= 0._wp )  THEN
!--         TODO detect and log when normalization differs too much from 1
            svf(1, isvf_surflt:nsvfl) = svf(1, isvf_surflt:nsvfl) / svfsum
        ENDIF

!--     deallocate temporary asvf array
!--     DEALLOCATE(asvf) - ifort has a problem with deallocation of allocatable target
!--     via pointing pointer - we need to test original targets
        IF ( ALLOCATED(asvf1) )  THEN
            DEALLOCATE(asvf1)
        ENDIF
        IF ( ALLOCATED(asvf2) )  THEN
            DEALLOCATE(asvf2)
        ENDIF

        npcsfl = 0
        IF ( plant_canopy )  THEN

            CALL location_message( '    calculation of the complete CSF array', .TRUE. )

!--         sort and merge csf for the last time, keeping the array size to minimum
            CALL usm_merge_and_grow_csf(-1)
            
!--         aggregate csb among processors
!--         allocate necessary arrays
            ALLOCATE( csflt(ndcsf,max(ncsfl,ndcsf)) )
            ALLOCATE( kcsflt(kdcsf,max(ncsfl,kdcsf)) )
            ALLOCATE( icsflt(0:numprocs-1) )
            ALLOCATE( dcsflt(0:numprocs-1) )
            ALLOCATE( ipcsflt(0:numprocs-1) )
            ALLOCATE( dpcsflt(0:numprocs-1) )
            
!--         fill out arrays of csf values and 
!--         arrays of number of elements and displacements
!--         for particular precessors
            icsflt = 0
            dcsflt = 0
            ip = -1
            j = -1
            d = 0
            DO kcsf = 1, ncsfl
                j = j+1
                IF ( acsf(kcsf)%ip /= ip )  THEN
!--                 new block of the processor
!--                 number of elements of previous block
                    IF ( ip>=0) icsflt(ip) = j
                    d = d+j
!--                 blank blocks
                    DO jp = ip+1, acsf(kcsf)%ip-1
!--                     number of elements is zero, displacement is equal to previous
                        icsflt(jp) = 0
                        dcsflt(jp) = d
                    ENDDO
!--                 the actual block
                    ip = acsf(kcsf)%ip
                    dcsflt(ip) = d
                    j = 0
                ENDIF
!--             fill out real values of rsvf, rtransp
                csflt(1,kcsf) = acsf(kcsf)%rsvf
                csflt(2,kcsf) = acsf(kcsf)%rtransp
!--             fill out integer values of itz,ity,itx,isurfs
                kcsflt(1,kcsf) = acsf(kcsf)%itz
                kcsflt(2,kcsf) = acsf(kcsf)%ity
                kcsflt(3,kcsf) = acsf(kcsf)%itx
                kcsflt(4,kcsf) = acsf(kcsf)%isurfs
            ENDDO
!--         last blank blocks at the end of array
            j = j+1
            IF ( ip>=0 ) icsflt(ip) = j
            d = d+j
            DO jp = ip+1, numprocs-1
!--             number of elements is zero, displacement is equal to previous
                icsflt(jp) = 0
                dcsflt(jp) = d
            ENDDO
            
!--         deallocate temporary acsf array
!--         DEALLOCATE(acsf) - ifort has a problem with deallocation of allocatable target
!--         via pointing pointer - we need to test original targets
            IF ( ALLOCATED(acsf1) )  THEN
                DEALLOCATE(acsf1)
            ENDIF
            IF ( ALLOCATED(acsf2) )  THEN
                DEALLOCATE(acsf2)
            ENDIF
                    
#if defined( __parallel )
!--         scatter and gather the number of elements to and from all processor
!--         and calculate displacements
            CALL MPI_AlltoAll(icsflt,1,MPI_INTEGER,ipcsflt,1,MPI_INTEGER,comm2d, ierr)
            
            npcsfl = SUM(ipcsflt)
            d = 0
            DO i = 0, numprocs-1
                dpcsflt(i) = d
                d = d + ipcsflt(i)
            ENDDO
        
!--         exchange csf fields between processors
            ALLOCATE( pcsflt(ndcsf,max(npcsfl,ndcsf)) )
            ALLOCATE( kpcsflt(kdcsf,max(npcsfl,kdcsf)) )
            CALL MPI_AlltoAllv(csflt, ndcsf*icsflt, ndcsf*dcsflt, MPI_REAL, &
                pcsflt, ndcsf*ipcsflt, ndcsf*dpcsflt, MPI_REAL, comm2d, ierr)
            CALL MPI_AlltoAllv(kcsflt, kdcsf*icsflt, kdcsf*dcsflt, MPI_INTEGER, &
                kpcsflt, kdcsf*ipcsflt, kdcsf*dpcsflt, MPI_INTEGER, comm2d, ierr)
            
#else
            npcsfl = ncsfl
            ALLOCATE( pcsflt(ndcsf,max(npcsfl,ndcsf)) )
            ALLOCATE( kpcsflt(kdcsf,max(npcsfl,kdcsf)) )
            pcsflt = csflt
            kpcsflt = kcsflt
#endif

!--         deallocate temporary arrays
            DEALLOCATE( csflt )
            DEALLOCATE( kcsflt )
            DEALLOCATE( icsflt )
            DEALLOCATE( dcsflt )
            DEALLOCATE( ipcsflt )
            DEALLOCATE( dpcsflt )

!--         sort csf ( a version of quicksort )
            CALL quicksort_csf2(kpcsflt, pcsflt, 1, npcsfl)

!--         aggregate canopy sink factor records with identical box & source
!--         againg across all values from all processors
            IF ( npcsfl > 0 )  THEN
                icsf = 1 !< reading index
                kcsf = 1 !< writing index
                DO while (icsf < npcsfl)
!--                 here kpcsf(kcsf) already has values from kpcsf(icsf)
                    IF ( kpcsflt(3,icsf) == kpcsflt(3,icsf+1)  .AND.  &
                         kpcsflt(2,icsf) == kpcsflt(2,icsf+1)  .AND.  &
                         kpcsflt(1,icsf) == kpcsflt(1,icsf+1)  .AND.  &
                         kpcsflt(4,icsf) == kpcsflt(4,icsf+1) )  THEN
!--                     We could simply take either first or second rtransp, both are valid. As a very simple heuristic about which ray
!--                     probably passes nearer the center of the target box, we choose DIF from the entry with greater CSF, since that
!--                     might mean that the traced beam passes longer through the canopy box.
                        IF ( pcsflt(1,kcsf) < pcsflt(1,icsf+1) )  THEN
                            pcsflt(2,kcsf) = pcsflt(2,icsf+1)
                        ENDIF
                        pcsflt(1,kcsf) = pcsflt(1,kcsf) + pcsflt(1,icsf+1)

!--                     advance reading index, keep writing index
                        icsf = icsf + 1
                    ELSE
!--                     not identical, just advance and copy
                        icsf = icsf + 1
                        kcsf = kcsf + 1
                        kpcsflt(:,kcsf) = kpcsflt(:,icsf)
                        pcsflt(:,kcsf) = pcsflt(:,icsf)
                    ENDIF
                ENDDO
!--             last written item is now also the last item in valid part of array
                npcsfl = kcsf
            ENDIF

            ncsfl = npcsfl
            IF ( ncsfl > 0 )  THEN
                ALLOCATE( csf(ndcsf,ncsfl) )
                ALLOCATE( csfsurf(idcsf,ncsfl) )
                DO icsf = 1, ncsfl
                    csf(:,icsf) = pcsflt(:,icsf)
                    csfsurf(1,icsf) =  gridpcbl(kpcsflt(1,icsf),kpcsflt(2,icsf),kpcsflt(3,icsf))
                    csfsurf(2,icsf) =  kpcsflt(4,icsf)
                ENDDO
            ENDIF
            
!--         deallocation of temporary arrays
            DEALLOCATE( pcsflt )
            DEALLOCATE( kpcsflt )
            
        ENDIF
        
        RETURN
        
301     WRITE( message_string, * )  &
            'I/O error when processing shape view factors / ',  &
            'plant canopy sink factors / direct irradiance factors.'
        CALL message( 'init_urban_surface', 'PA0502', 2, 2, 0, 6, 0 )
       
    END SUBROUTINE usm_calc_svf


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine checks variables and assigns units.
!> It is caaled out from subroutine check_parameters.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_check_data_output( variable, unit )
        
        IMPLICIT NONE
 
        CHARACTER (len=*),INTENT(IN)    ::  variable !:
        CHARACTER (len=*),INTENT(OUT)   ::  unit     !:
        
        CHARACTER (len=varnamelength)   :: var

        var = TRIM(variable)
        IF ( var(1:12) == 'usm_rad_net_'  .OR.  var(1:13) == 'usm_rad_insw_'  .OR.        &
             var(1:13) == 'usm_rad_inlw_'  .OR.  var(1:16) == 'usm_rad_inswdir_'  .OR.    &
             var(1:16) == 'usm_rad_inswdif_'  .OR.  var(1:16) == 'usm_rad_inswref_'  .OR. &
             var(1:16) == 'usm_rad_inlwdif_'  .OR.  var(1:16) == 'usm_rad_inlwref_'  .OR. &
             var(1:14) == 'usm_rad_outsw_'  .OR.  var(1:14) == 'usm_rad_outlw_'  .OR.     &
             var(1:14) == 'usm_rad_ressw_'  .OR.  var(1:14) == 'usm_rad_reslw_'  .OR.     &
             var(1:11) == 'usm_rad_hf_'  .OR.                                             &
             var(1:9)  == 'usm_wshf_'  .OR.  var(1:9) == 'usm_wghf_' )  THEN
            unit = 'W/m2'
        ELSE IF ( var(1:10) == 'usm_t_surf'  .OR.  var(1:10) == 'usm_t_wall' )  THEN
            unit = 'K'
        ELSE IF ( var(1:9) == 'usm_surfz'  .OR.  var(1:7) == 'usm_svf'  .OR.              & 
                  var(1:7) == 'usm_dif'  .OR.  var(1:11) == 'usm_surfcat'  .OR.           &
                  var(1:11) == 'usm_surfalb'  .OR.  var(1:12) == 'usm_surfemis')  THEN
            unit = '1'
        ELSE IF ( plant_canopy  .AND.  var(1:7) == 'usm_lad' )  THEN
            unit = 'm2/m3'
        ELSE IF ( plant_canopy  .AND.  var(1:13) == 'usm_canopy_hr' )  THEN
            unit = 'K/s'
        ELSE
            unit = 'illegal'
        ENDIF

    END SUBROUTINE usm_check_data_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for urban surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_check_parameters
    
       USE control_parameters,                                                 &
           ONLY:  bc_pt_b, bc_q_b, constant_flux_layer, large_scale_forcing,   &
                  lsf_surf, topography

!
!--    Dirichlet boundary conditions are required as the surface fluxes are
!--    calculated from the temperature/humidity gradients in the urban surface
!--    model
       IF ( bc_pt_b == 'neumann'   .OR.   bc_q_b == 'neumann' )  THEN
          message_string = 'urban surface model requires setting of '//        &
                           'bc_pt_b = "dirichlet" and '//                      &
                           'bc_q_b  = "dirichlet"'
          CALL message( 'check_parameters', 'PA0590', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( .NOT.  constant_flux_layer )  THEN
          message_string = 'urban surface model requires '//                   &
                           'constant_flux_layer = .T.'
          CALL message( 'check_parameters', 'PA0591', 1, 2, 0, 6, 0 )
       ENDIF
!        
!--    Surface forcing has to be disabled for LSF in case of enabled 
!--    urban surface module
       IF ( large_scale_forcing )  THEN
          lsf_surf = .FALSE.
       ENDIF
!
!--    Topography
       IF ( topography == 'flat' )  THEN
          message_string = 'topography /= "flat" is required '//               &
                           'when using the urban surface model'
          CALL message( 'check_parameters', 'PA0592', 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE usm_check_parameters


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Output of the 3D-arrays in netCDF and/or AVS format 
!> for variables of urban_surface model.
!> It resorts the urban surface module output quantities from surf style
!> indexing into temporary 3D array with indices (i,j,k).
!> It is called from subroutine data_output_3d.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
        
        IMPLICIT NONE

        INTEGER(iwp), INTENT(IN)       ::  av        !< 
        CHARACTER (len=*), INTENT(IN)  ::  variable  !< 
        INTEGER(iwp), INTENT(IN)       ::  nzb_do    !< lower limit of the data output (usually 0)
        INTEGER(iwp), INTENT(IN)       ::  nzt_do    !< vertical upper limit of the data output (usually nz_do3d)
        LOGICAL, INTENT(OUT)           ::  found     !< 
        REAL(sp), DIMENSION(nxlg:nxrg,nysg:nyng,nzb_do:nzt_do) ::  local_pf   !< sp - it has to correspond to module data_output_3d
        REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg)     ::  temp_pf    !< temp array for urban surface output procedure
        
        CHARACTER (len=varnamelength)                          :: var, surfid
        INTEGER(iwp), PARAMETER                                :: nd = 5
        CHARACTER(len=6), DIMENSION(0:nd-1), PARAMETER         :: dirname = (/ '_roof ', '_south', '_north', '_west ', '_east ' /)
        INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER             :: dirint = (/ iroof, isouth, inorth, iwest, ieast /)
        INTEGER(iwp), DIMENSION(0:nd-1)                        :: dirstart
        INTEGER(iwp), DIMENSION(0:nd-1)                        :: dirend
        INTEGER(iwp)                                           :: ids,isurf,isvf,isurfs,isurflt
        INTEGER(iwp)                                           :: is,js,ks,i,j,k,iwl,istat

        dirstart = (/ startland, startwall, startwall, startwall, startwall /)
        dirend = (/ endland, endwall, endwall, endwall, endwall /)

        found = .TRUE.
        temp_pf = -1._wp
        
        ids = -1
        var = TRIM(variable)
        DO i = 0, nd-1
            k = len(TRIM(var))
            j = len(TRIM(dirname(i)))
            IF ( var(k-j+1:k) == dirname(i) )  THEN
                ids = i
                var = var(:k-j)
                EXIT
            ENDIF
        ENDDO
        IF ( ids == -1 )  THEN
            var = TRIM(variable)
        ENDIF
        IF ( var(1:11) == 'usm_t_wall_'  .AND.  len(TRIM(var)) >= 12 )  THEN
!--         wall layers
            READ(var(12:12), '(I1)', iostat=istat ) iwl
            IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
                var = var(1:10)
            ENDIF
        ENDIF
        IF ( (var(1:8) == 'usm_svf_'  .OR.  var(1:8) == 'usm_dif_')  .AND.  len(TRIM(var)) >= 13 )  THEN
!--         svf values to particular surface
            surfid = var(9:)
            i = index(surfid,'_')
            j = index(surfid(i+1:),'_')
            READ(surfid(1:i-1),*, iostat=istat ) is
            IF ( istat == 0 )  THEN
                READ(surfid(i+1:i+j-1),*, iostat=istat ) js
            ENDIF
            IF ( istat == 0 )  THEN
                READ(surfid(i+j+1:),*, iostat=istat ) ks
            ENDIF
            IF ( istat == 0 )  THEN
                var = var(1:7)
            ENDIF
        ENDIF
        
        SELECT CASE ( TRIM(var) )

          CASE ( 'usm_surfz' )
!--           array of lw radiation falling to local surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                  IF ( surfl(id,isurf) == ids )  THEN
                      IF ( surfl(id,isurf) == iroof )  THEN
                          temp_pf(0,surfl(iy,isurf),surfl(ix,isurf)) =             &
                                  max(temp_pf(0,surfl(iy,isurf),surfl(ix,isurf)),  &
                                      REAL(surfl(iz,isurf),wp))
                      ELSE
                          temp_pf(0,surfl(iy,isurf),surfl(ix,isurf)) =             &
                                  max(temp_pf(0,surfl(iy,isurf),surfl(ix,isurf)),  &
                                      REAL(surfl(iz,isurf),wp)+1.0_wp)
                      ENDIF
                  ENDIF
              ENDDO

          CASE ( 'usm_surfcat' )
!--           surface category
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surface_types(isurf)
                 ENDIF
              ENDDO
              
          CASE ( 'usm_surfalb' )
!--           surface albedo
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = albedo_surf(isurf)
                 ENDIF
              ENDDO
              
          CASE ( 'usm_surfemis' )
!--           surface albedo
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = emiss_surf(isurf)
                 ENDIF
              ENDDO
              
          CASE ( 'usm_svf', 'usm_dif' )
!--           shape view factors or iradiance factors to selected surface
              IF ( TRIM(var)=='usm_svf' )  THEN
                  k = 1
              ELSE
                  k = 2
              ENDIF
              DO isvf = 1, nsvfl
                  isurflt = svfsurf(1, isvf)
                  isurfs = svfsurf(2, isvf)
                             
                  IF ( surf(ix,isurfs) == is  .AND.  surf(iy,isurfs) == js  .AND.       &
                       surf(iz,isurfs) == ks  .AND.  surf(id,isurfs) == ids )  THEN
  !--                 correct source surface
                      temp_pf(surfl(iz,isurflt),surfl(iy,isurflt),surfl(ix,isurflt)) = svf(k,isvf)
                  ENDIF
              ENDDO

          CASE ( 'usm_rad_net' )
!--           array of complete radiation balance
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = rad_net_l(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = rad_net_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_insw' )
!--           array of sw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinsw(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinsw_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inlw' )
!--           array of lw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlw(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlw_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inswdir' )
!--           array of direct sw radiation falling to surface from sun
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswdir(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswdir_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inswdif' )
!--           array of difusion sw radiation falling to surface from sky and borders of the domain
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswdif(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswdif_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inswref' )
!--           array of sw radiation falling to surface from reflections
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = &
                       surfinsw(isurf) - surfinswdir(isurf) - surfinswdif(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinswref_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inlwdif' )
!--           array of sw radiation falling to surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlwdif(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlwdif_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_inlwref' )
!--           array of lw radiation falling to surface from reflections
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlw(isurf) - surfinlwdif(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinlwref_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_outsw' )
!--           array of sw radiation emitted from surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfoutsw(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfoutsw_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_outlw' )
!--           array of lw radiation emitted from surface after i-th reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfoutlw(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfoutlw_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_ressw' )
!--           average of array of residua of sw radiation absorbed in surface after last reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfins(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfins_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_reslw' )
!--           average of array of residua of lw radiation absorbed in surface after last reflection
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinl(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfinl_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_rad_hf' )
!--           array of heat flux from radiation for surfaces after all reflections
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfhf(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = surfhf_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_wshf' )
!--           array of sensible heat flux from surfaces
!--           horizontal surfaces
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = wshf_eb(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = wshf_eb_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_wghf' )
!--           array of heat flux from ground (land, wall, roof)
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = wghf_eb(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = wghf_eb_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_t_surf' )
!--           surface temperature for surfaces
              DO isurf = max(startenergy,dirstart(ids)), min(endenergy,dirend(ids))
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = t_surf(isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = t_surf_av(isurf)
                   ENDIF
                 ENDIF
              ENDDO
              
          CASE ( 'usm_t_wall' )
!--           wall temperature for  iwl layer of walls and land
              DO isurf = dirstart(ids), dirend(ids)
                 IF ( surfl(id,isurf) == ids )  THEN
                   IF ( av == 0 )  THEN
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = t_wall(iwl,isurf)
                   ELSE
                     temp_pf(surfl(iz,isurf),surfl(iy,isurf),surfl(ix,isurf)) = t_wall_av(iwl,isurf)
                   ENDIF
                 ENDIF
              ENDDO

          CASE ( 'usm_lad' )
!--           leaf area density
              DO i = nxl, nxr
                 DO j = nys, nyn
                     DO k = nzb_s_inner(j,i), nzut
                         temp_pf(k,j,i) = lad_s(k-nzb_s_inner(j,i),j,i)
                     ENDDO
                 ENDDO
              ENDDO
              
          CASE ( 'usm_canopy_hr' )
!--           canopy heating rate
              DO i = nxl, nxr
                 DO j = nys, nyn
                     DO k = nzb_s_inner(j,i), nzut
                         temp_pf(k,j,i) = pc_heating_rate(k-nzb_s_inner(j,i),j,i)
                     ENDDO
                 ENDDO
              ENDDO
             
          CASE DEFAULT
              found = .FALSE.
              
        END SELECT
        
!--     fill out array local_pf which is subsequently treated by data_output_3d
        CALL exchange_horiz( temp_pf, nbgp )
        DO j = nysg,nyng
            DO i = nxlg,nxrg
                DO k = nzb_do, nzt_do
                    local_pf(i,j,k) = temp_pf(k,j,i)
                ENDDO
            ENDDO
        ENDDO
        
    END SUBROUTINE usm_data_output_3d
    

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Soubroutine defines appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )
    
        IMPLICIT NONE

        CHARACTER (len=*), INTENT(IN)  ::  variable    !< 
        LOGICAL, INTENT(OUT)           ::  found       !< 
        CHARACTER (len=*), INTENT(OUT) ::  grid_x      !< 
        CHARACTER (len=*), INTENT(OUT) ::  grid_y      !< 
        CHARACTER (len=*), INTENT(OUT) ::  grid_z      !< 

        CHARACTER (len=varnamelength)  :: var

        var = TRIM(variable)
        IF ( var(1:12) == 'usm_rad_net_'  .OR.  var(1:13) == 'usm_rad_insw_'  .OR.          &
             var(1:13) == 'usm_rad_inlw_'  .OR.  var(1:16) == 'usm_rad_inswdir_'  .OR.      &
             var(1:16) == 'usm_rad_inswdif_'  .OR.  var(1:16) == 'usm_rad_inswref_'  .OR.   &
             var(1:16) == 'usm_rad_inlwdif_'  .OR.  var(1:16) == 'usm_rad_inlwref_'  .OR.   &
             var(1:14) == 'usm_rad_outsw_'  .OR.  var(1:14) == 'usm_rad_outlw_'  .OR.       &
             var(1:14) == 'usm_rad_ressw_'  .OR.  var(1:14) == 'usm_rad_reslw_'  .OR.       &
             var(1:11) == 'usm_rad_hf_'  .OR.                                               &
             var(1:9) == 'usm_wshf_'  .OR.  var(1:9) == 'usm_wghf_'  .OR.                   &
             var(1:10) == 'usm_t_surf'  .OR.  var(1:10) == 'usm_t_wall'  .OR.               &
             var(1:9) == 'usm_surfz'  .OR.  var(1:7) == 'usm_svf'  .OR.                     & 
             var(1:7) == 'usm_dif'  .OR.  var(1:11) == 'usm_surfcat'  .OR.                  &
             var(1:11) == 'usm_surfalb'  .OR.  var(1:12) == 'usm_surfemis'  .OR.            &
             var(1:7) == 'usm_lad'  .OR.  var(1:13) == 'usm_canopy_hr' )  THEN

            found = .TRUE.
            grid_x = 'x'
            grid_y = 'y'
            grid_z = 'zu'
        ELSE
            found  = .FALSE.
            grid_x = 'none'
            grid_y = 'none'
            grid_z = 'none'
        ENDIF

    END SUBROUTINE usm_define_netcdf_grid
    
    
!------------------------------------------------------------------------------!
!> Finds first model boundary crossed by a ray
!------------------------------------------------------------------------------!
    PURE SUBROUTINE usm_find_boundary_face(origin, uvect, bdycross)
        IMPLICIT NONE
        REAL(wp), DIMENSION(3), INTENT(in)      :: origin    !< ray origin
        REAL(wp), DIMENSION(3), INTENT(in)      :: uvect     !< ray unit vector
        INTEGER(iwp), DIMENSION(4), INTENT(out) :: bdycross  !< found boundary crossing (d, z, y, x)
        REAL(wp), DIMENSION(3)                  :: crossdist !< crossing distance
        INTEGER(iwp), DIMENSION(3)              :: bdyd      !< boundary direction
        REAL(wp)                                :: bdydim    !< 
        REAL(wp)                                :: dist      !< 
        INTEGER(iwp)                            :: seldim    !< found fist crossing index
        INTEGER(iwp)                            :: d         !< 

        bdydim = nzut + .5_wp !< top boundary
        bdyd(1) = isky
        crossdist(1) = (bdydim - origin(1)) / uvect(1)

        IF ( uvect(2) >= 0._wp )  THEN
            bdydim = ny + .5_wp !< north global boundary
            bdyd(2) = inorthb
        ELSE
            bdydim = -.5_wp !< south global boundary
            bdyd(2) = isouthb
        ENDIF
        crossdist(2) = (bdydim - origin(2)) / uvect(2)

        IF ( uvect(3) >= 0._wp )  THEN
            bdydim = nx + .5_wp !< east global boundary
            bdyd(3) = ieastb
        ELSE
            bdydim = -.5_wp !< west global boundary
            bdyd(3) = iwestb
        ENDIF
        crossdist(3) = (bdydim - origin(3)) / uvect(3)

        seldim = minloc(crossdist, 1)
        dist = crossdist(seldim)
        d = bdyd(seldim)

        bdycross(1) = d
        bdycross(2:4) = NINT( origin(:) + uvect(:)*dist &
                        + .5_wp * (/ kdir(d), jdir(d), idir(d) /) )
    END SUBROUTINE


!------------------------------------------------------------------------------!
!> Determines whether two faces are oriented towards each other
!------------------------------------------------------------------------------!
    PURE LOGICAL FUNCTION usm_facing(x, y, z, d, x2, y2, z2, d2)
        IMPLICIT NONE
        INTEGER(iwp),   INTENT(in)  :: x, y, z, d, x2, y2, z2, d2
      
        usm_facing = .FALSE.
        IF ( d==iroof  .AND.  d2==iroof ) RETURN
        IF ( d==isky  .AND.  d2==isky ) RETURN
        IF ( (d==isouth  .OR.  d==inorthb)  .AND.  (d2==isouth.OR.d2==inorthb) ) RETURN
        IF ( (d==inorth  .OR.  d==isouthb)  .AND.  (d2==inorth.OR.d2==isouthb) ) RETURN
        IF ( (d==iwest  .OR.  d==ieastb)  .AND.  (d2==iwest.OR.d2==ieastb) ) RETURN
        IF ( (d==ieast  .OR.  d==iwestb)  .AND.  (d2==ieast.OR.d2==iwestb) ) RETURN

        SELECT CASE (d)
            CASE (iroof)                   !< ground, roof
                IF ( z2 < z ) RETURN
            CASE (isky)                    !< sky
                IF ( z2 > z ) RETURN
            CASE (isouth, inorthb)         !< south facing
                IF ( y2 > y ) RETURN
            CASE (inorth, isouthb)         !< north facing
                IF ( y2 < y ) RETURN
            CASE (iwest, ieastb)           !< west facing
                IF ( x2 > x ) RETURN
            CASE (ieast, iwestb)           !< east facing
                IF ( x2 < x ) RETURN
        END SELECT

        SELECT CASE (d2)
            CASE (iroof)                   !< ground, roof
                IF ( z < z2 ) RETURN
            CASE (isky)                    !< sky
                IF ( z > z2 ) RETURN
            CASE (isouth, inorthb)         !< south facing
                IF ( y > y2 ) RETURN
            CASE (inorth, isouthb)         !< north facing
                IF ( y < y2 ) RETURN
            CASE (iwest, ieastb)           !< west facing
                IF ( x > x2 ) RETURN
            CASE (ieast, iwestb)           !< east facing
                IF ( x < x2 ) RETURN
            CASE (-1)
                CONTINUE
        END SELECT

        usm_facing = .TRUE.
        
    END FUNCTION usm_facing
    

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the wall surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_init_material_model

        IMPLICIT NONE

        INTEGER(iwp) ::  k, l            !< running indices
        
        CALL location_message( '    initialization of wall surface model', .TRUE. )
       
!--     Calculate wall grid spacings. 
!--     Temperature is defined at the center of the wall layers,
!--     whereas gradients/fluxes are defined at the edges (_stag)
        DO l = nzb_wall, nzt_wall
           zwn(l) = zwn_default(l)
        ENDDO
        
!--     apply for all particular wall grids
        DO l = startenergy, endenergy
           zw(:,l) = zwn(:) * thickness_wall(l)
           dz_wall(nzb_wall,l) = zw(nzb_wall,l)
           DO k = nzb_wall+1, nzt_wall
               dz_wall(k,l) = zw(k,l) - zw(k-1,l)
           ENDDO
           
           dz_wall(nzt_wall+1,l) = dz_wall(nzt_wall,l)

           DO k = nzb_wall, nzt_wall-1
               dz_wall_stag(k,l) = 0.5 * (dz_wall(k+1,l) + dz_wall(k,l))
           ENDDO
           dz_wall_stag(nzt_wall,l) = dz_wall(nzt_wall,l)
        ENDDO
        
        ddz_wall      = 1.0_wp / dz_wall
        ddz_wall_stag = 1.0_wp / dz_wall_stag
        
        CALL location_message( '    wall structures filed out', .TRUE. )

        CALL location_message( '    initialization of wall surface model finished', .TRUE. )

    END SUBROUTINE usm_init_material_model

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the urban surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_init_urban_surface
    
        IMPLICIT NONE

        INTEGER(iwp) ::  i, j, k, l            !< running indices
        REAL(wp)     ::  c, d, tin, exn
        
        CALL cpu_log( log_point_s(78), 'usm_init', 'start' )
!--     surface forcing have to be disabled for LSF 
!--     in case of enabled urban surface module
        IF ( large_scale_forcing )  THEN
            lsf_surf = .FALSE.
        ENDIF
       
!--     init anthropogenic sources of heat
        CALL usm_allocate_urban_surface()
       
!--     read the surface_types array somewhere
        CALL usm_read_urban_surface_types()
        
!--     init material heat model
        CALL usm_init_material_model()
        
        IF ( usm_anthropogenic_heat )  THEN
!--         init anthropogenic sources of heat (from transportation for now)
            CALL usm_read_anthropogenic_heat()
        ENDIF
        
        IF ( read_svf_on_init )  THEN
!--         read svf, csf, svfsurf and csfsurf data from file
            CALL location_message( '    Start reading SVF from file', .TRUE. )
            CALL usm_read_svf_from_file()
            CALL location_message( '    Reading SVF from file has finished', .TRUE. )
        ELSE
!--         calculate SFV and CSF
            CALL location_message( '    Start calculation of SVF', .TRUE. )
            CALL cpu_log( log_point_s(79), 'usm_calc_svf', 'start' )
            CALL usm_calc_svf()
            CALL cpu_log( log_point_s(79), 'usm_calc_svf', 'stop' )
            CALL location_message( '    Calculation of SVF has finished', .TRUE. )
        ENDIF

        IF ( write_svf_on_init )  THEN
!--         write svf, csf svfsurf and csfsurf data to file
            CALL location_message( '    Store SVF and CSF to file', .TRUE. )
            CALL usm_write_svf_to_file()
        ENDIF
        
        IF ( plant_canopy )  THEN
!--         gridpcbl was only necessary for initialization
            DEALLOCATE( gridpcbl )
            IF ( .NOT.  ALLOCATED(pc_heating_rate) )  THEN
!--             then pc_heating_rate is allocated in init_plant_canopy
!--             in case of cthf /= 0 => we need to allocate it for our use here
                ALLOCATE( pc_heating_rate(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
            ENDIF
        ENDIF

!--     Intitialization of the surface and wall/ground/roof temperature

!--     Initialization for restart runs
        IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN

!--         restore data from restart file
            CALL usm_read_restart_data()
        ELSE
       
!--         Calculate initial surface temperature
            exn = ( surface_pressure / 1000.0_wp )**0.286_wp

            DO l = startenergy, endenergy
                k = surfl(iz,l)
                j = surfl(iy,l)
                i = surfl(ix,l)

!--              Initial surface temperature set from pt of adjacent gridbox
                t_surf(l) = pt(k,j,i) * exn
            ENDDO
      
!--         initial values for t_wall
!--         outer value is set to surface temperature
!--         inner value is set to wall_inner_temperature
!--         and profile is logaritmic (linear in nz)
            DO l = startenergy, endenergy
                IF ( isroof_surf(l) )  THEN
                    tin = roof_inner_temperature
                ELSE IF ( surf(id,l)==iroof )  THEN
                    tin = soil_inner_temperature
                ELSE
                    tin = wall_inner_temperature
                ENDIF
                DO k = nzb_wall, nzt_wall+1
                    c = REAL(k-nzb_wall,wp)/REAL(nzt_wall+1-nzb_wall,wp)
                    t_wall(k,:) = (1.0_wp-c)*t_surf(:) + c*tin
                ENDDO
            ENDDO
        ENDIF
        
!--    
!--        Possibly DO user-defined actions (e.g. define heterogeneous wall surface)
        CALL user_init_urban_surface

!--     initialize prognostic values for the first timestep
        t_surf_p = t_surf
        t_wall_p = t_wall
        
!--     Adjust radiative fluxes for urban surface at model start
        CALL usm_radiation
        
        CALL cpu_log( log_point_s(78), 'usm_init', 'stop' )
        
        
    END SUBROUTINE usm_init_urban_surface


!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Wall model as part of the urban surface model. The model predicts wall
!> temperature.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_material_heat_model


        IMPLICIT NONE

        INTEGER(iwp) ::  i,j,k,l,kw                      !< running indices

        REAL(wp), DIMENSION(nzb_wall:nzt_wall) :: wtend  !< tendency

                                               
        DO l = startenergy, endenergy
!--         calculate frequently used parameters
            k = surfl(iz,l)
            j = surfl(iy,l)
            i = surfl(ix,l)

            !
!--         prognostic equation for ground/wall/roof temperature t_wall
            wtend(:) = 0.0_wp
            wtend(nzb_wall) = (1.0_wp/rho_c_wall(nzb_wall,l)) *                     &
                       ( lambda_h(nzb_wall,l) * ( t_wall(nzb_wall+1,l)              &
                         - t_wall(nzb_wall,l) ) * ddz_wall(nzb_wall+1,l)            &
                         + wghf_eb(l) ) * ddz_wall_stag(nzb_wall,l)
            
            DO  kw = nzb_wall+1, nzt_wall
                wtend(kw) = (1.0_wp/rho_c_wall(kw,l))                               &
                              * (   lambda_h(kw,l)                                  &
                                 * ( t_wall(kw+1,l) - t_wall(kw,l) )                &
                                 * ddz_wall(kw+1,l)                                 &
                              - lambda_h(kw-1,l)                                    &
                                 * ( t_wall(kw,l) - t_wall(kw-1,l) )                &
                                 * ddz_wall(kw,l)                                   &
                              ) * ddz_wall_stag(kw,l)
            ENDDO

            t_wall_p(nzb_wall:nzt_wall,l) = t_wall(nzb_wall:nzt_wall,l)             &
                                             + dt_3d * ( tsc(2)                     &
                                             * wtend(nzb_wall:nzt_wall) + tsc(3)    &
                                             * tt_wall_m(nzb_wall:nzt_wall,l) )   
            
            !
!--         calculate t_wall tendencies for the next Runge-Kutta step
            IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   DO  kw = nzb_wall, nzt_wall
                      tt_wall_m(kw,l) = wtend(kw)
                   ENDDO
                ELSEIF ( intermediate_timestep_count <                              &
                         intermediate_timestep_count_max )  THEN
                    DO  kw = nzb_wall, nzt_wall
                        tt_wall_m(kw,l) = -9.5625_wp * wtend(kw) + 5.3125_wp        &
                                         * tt_wall_m(kw,l)
                    ENDDO
                ENDIF
            ENDIF
        ENDDO

    END SUBROUTINE usm_material_heat_model


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &usm_par for urban surface model
!------------------------------------------------------------------------------!
    SUBROUTINE usm_parin

       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< string containing current line of file PARIN

       NAMELIST /urban_surface_par/                                            &
                           land_category,                                      &
                           mrt_factors,                                        &
                           nrefsteps,                                          &
                           pedestrant_category,                                &
                           ra_horiz_coef,                                      &
                           read_svf_on_init,                                   &
                           roof_category,                                      &
                           split_diffusion_radiation,                          &
                           urban_surface,                                      &
                           usm_anthropogenic_heat,                             &
                           usm_energy_balance_land,                            &
                           usm_energy_balance_wall,                            &
                           usm_material_model,                                 &
                           usm_lad_rma,                                        &
                           wall_category,                                      &
                           write_svf_on_init

       line = ' '

!
!--    Try to find urban surface model package
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&urban_surface_par' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, urban_surface_par )

!
!--    Set flag that indicates that the land surface model is switched on
       urban_surface = .TRUE.


 10    CONTINUE

    END SUBROUTINE usm_parin


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine calculates interaction of the solar radiation
!> with urban surface and updates surface, roofs and walls heatfluxes.
!> It also updates rad_sw_out and rad_lw_out.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_radiation
    
        IMPLICIT NONE
       
        INTEGER(iwp)               :: i, j, k, kk, is, js, d, ku, refstep
        INTEGER(iwp)               :: nzubl, nzutl, isurf, isurfsrc, isurf1, isvf, icsf, ipcgb
        INTEGER(iwp), DIMENSION(4) :: bdycross
        REAL(wp), DIMENSION(3,3)   :: mrot            !< grid rotation matrix (xyz)
        REAL(wp), DIMENSION(3,0:9) :: vnorm           !< face direction normal vectors (xyz)
        REAL(wp), DIMENSION(3)     :: sunorig         !< grid rotated solar direction unit vector (xyz)
        REAL(wp), DIMENSION(3)     :: sunorig_grid    !< grid squashed solar direction unit vector (zyx)
        REAL(wp), DIMENSION(0:9)   :: costheta        !< direct irradiance factor of solar angle
        REAL(wp), DIMENSION(nzub:nzut) :: pchf_prep   !< precalculated factor for canopy temp tendency
        REAL(wp), PARAMETER        :: alpha = 0._wp   !< grid rotation (TODO: add to namelist or remove)
        REAL(wp)                   :: rx, ry, rz
        REAL(wp)                   :: pc_box_area, pc_abs_frac, pc_abs_eff
        INTEGER(iwp)               :: pc_box_dimshift !< transform for best accuracy
        
        
        IF ( plant_canopy )  THEN
            pchf_prep(:) = r_d * (hyp(nzub:nzut) / 100000.0_wp)**0.286_wp &
                        / (cp * hyp(nzub:nzut) * dx*dy*dz) !< equals to 1 / (rho_ocean * c_p * Vbox * T)
        ENDIF

        sun_direction = .TRUE.
        CALL calc_zenith  !< required also for diffusion radiation

!--     prepare rotated normal vectors and irradiance factor
        vnorm(1,:) = idir(:)
        vnorm(2,:) = jdir(:)
        vnorm(3,:) = kdir(:)
        mrot(1, :) = (/ cos(alpha), -sin(alpha), 0._wp /)
        mrot(2, :) = (/ sin(alpha),  cos(alpha), 0._wp /)
        mrot(3, :) = (/ 0._wp,       0._wp,      1._wp /)
        sunorig = (/ sun_dir_lon, sun_dir_lat, zenith(0) /)
        sunorig = matmul(mrot, sunorig)
        DO d = 0, 9
            costheta(d) = dot_product(sunorig, vnorm(:,d))
        ENDDO
        
        IF ( zenith(0) > 0 )  THEN
!--         now we will "squash" the sunorig vector by grid box size in
!--         each dimension, so that this new direction vector will allow us
!--         to traverse the ray path within grid coordinates directly
            sunorig_grid = (/ sunorig(3)/dz, sunorig(2)/dy, sunorig(1)/dx /)
!--         sunorig_grid = sunorig_grid / norm2(sunorig_grid)
            sunorig_grid = sunorig_grid / SQRT(SUM(sunorig_grid**2))

            IF ( plant_canopy )  THEN
!--            precompute effective box depth with prototype Leaf Area Density
               pc_box_dimshift = maxloc(sunorig, 1) - 1
               CALL usm_box_absorb(cshift((/dx,dy,dz/), pc_box_dimshift),      &
                                   60, prototype_lad,                          &
                                   cshift(sunorig, pc_box_dimshift),           &
                                   pc_box_area, pc_abs_frac)
               pc_box_area = pc_box_area * sunorig(pc_box_dimshift+1) / sunorig(3)
               pc_abs_eff = log(1._wp - pc_abs_frac) / prototype_lad
            ENDIF
        ENDIF
        
!--     split diffusion and direct part of the solar downward radiation
!--     comming from radiation model and store it in 2D arrays
!--     rad_sw_in_diff, rad_sw_in_dir and rad_lw_in_diff
        IF ( split_diffusion_radiation )  THEN
            CALL usm_calc_diffusion_radiation
        ELSE
            rad_sw_in_diff = 0.0_wp
            rad_sw_in_dir(:,:)  = rad_sw_in(0,:,:)
            rad_lw_in_diff(:,:) = rad_lw_in(0,:,:)
        ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     First pass: direct + diffuse irradiance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        surfinswdir   = 0._wp
        surfinswdif   = 0._wp
        surfinlwdif   = 0._wp
        surfins   = 0._wp
        surfinl   = 0._wp
        surfoutsl    = 0._wp
        surfoutll    = 0._wp
        
!--     Set up thermal radiation from surfaces
!--     emiss_surf is defined only for surfaces for which energy balance is calculated
        surfoutll(startenergy:endenergy) = emiss_surf(startenergy:endenergy) * sigma_sb   &
                                           * t_surf(startenergy:endenergy)**4
        
#if defined( __parallel )
!--     might be optimized and gather only values relevant for current processor
        CALL MPI_AllGatherv(surfoutll, nenergy, MPI_REAL, &
                            surfoutl, nsurfs, surfstart, MPI_REAL, comm2d, ierr)
#else
        surfoutl(:) = surfoutll(:)
#endif
        
        isurf1 = -1   !< previous processed surface
        DO isvf = 1, nsvfl
            isurf = svfsurf(1, isvf)
            k = surfl(iz, isurf)
            j = surfl(iy, isurf)
            i = surfl(ix, isurf)
            isurfsrc = svfsurf(2, isvf)
            IF ( zenith(0) > 0  .AND.  isurf /= isurf1 )  THEN
!--             locate the virtual surface where the direct solar ray crosses domain boundary
!--             (once per target surface)
                d = surfl(id, isurf)
                rz = REAL(k, wp) - 0.5_wp * kdir(d)
                ry = REAL(j, wp) - 0.5_wp * jdir(d)
                rx = REAL(i, wp) - 0.5_wp * idir(d)
                
                CALL usm_find_boundary_face( (/ rz, ry, rx /), sunorig_grid, bdycross)
                
                isurf1 = isurf
            ENDIF

            IF ( surf(id, isurfsrc) >= isky )  THEN
!--             diffuse rad from boundary surfaces. Since it is a simply
!--             calculated value, it is not assigned to surfref(s/l),
!--             instead it is used directly here
!--             we consider the radiation from the radiation model falling on surface
!--             as the radiation falling on the top of urban layer into the place of the source surface
!--             we consider it as a very reasonable simplification which allow as avoid
!--             necessity of other global range arrays and some all to all mpi communication
                surfinswdif(isurf) = surfinswdif(isurf) + rad_sw_in_diff(j,i) * svf(1,isvf) * svf(2,isvf)
                                                                !< canopy shading is applied only to shortwave
                surfinlwdif(isurf) = surfinlwdif(isurf) + rad_lw_in_diff(j,i) * svf(1,isvf)
            ELSE
!--             for surface-to-surface factors we calculate thermal radiation in 1st pass
                surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutl(isurfsrc)
            ENDIF
            
            IF ( zenith(0) > 0  .AND.  all( surf(:, isurfsrc) == bdycross ) )  THEN
!--             found svf between model boundary and the face => face isn't shaded
                surfinswdir(isurf) = rad_sw_in_dir(j, i) &
                    * costheta(surfl(id, isurf)) * svf(2,isvf) / zenith(0)

            ENDIF
        ENDDO

        IF ( plant_canopy )  THEN
        
            pcbinsw(:) = 0._wp
            pcbinlw(:) = 0._wp  !< will stay always 0 since we don't absorb lw anymore
            !
!--         pcsf first pass
            isurf1 = -1  !< previous processed pcgb
            DO icsf = 1, ncsfl
                ipcgb = csfsurf(1, icsf)
                i = pcbl(ix,ipcgb)
                j = pcbl(iy,ipcgb)
                k = pcbl(iz,ipcgb)
                isurfsrc = csfsurf(2, icsf)

                IF ( zenith(0) > 0  .AND.  ipcgb /= isurf1 )  THEN
!--                 locate the virtual surface where the direct solar ray crosses domain boundary
!--                 (once per target PC gridbox)
                    rz = REAL(k, wp)
                    ry = REAL(j, wp)
                    rx = REAL(i, wp)
                    CALL usm_find_boundary_face( (/ rz, ry, rx /), &
                        sunorig_grid, bdycross)

                    isurf1 = ipcgb
                ENDIF

                IF ( surf(id, isurfsrc) >= isky )  THEN
!--                 Diffuse rad from boundary surfaces. See comments for svf above.
                    pcbinsw(ipcgb) = pcbinsw(ipcgb) + svf(1,isvf) * svf(2,isvf) * rad_sw_in_diff(j,i)
!--                 canopy shading is applied only to shortwave, therefore no absorbtion for lw
!--                 pcbinlw(ipcgb) = pcbinlw(ipcgb) + svf(1,isvf) * rad_lw_in_diff(j,i)
                !ELSE
!--                 Thermal radiation in 1st pass
!--                 pcbinlw(ipcgb) = pcbinlw(ipcgb) + svf(1,isvf) * surfoutl(isurfsrc)
                ENDIF

                IF ( zenith(0) > 0  .AND.  all( surf(:, isurfsrc) == bdycross ) )  THEN
!--                 found svf between model boundary and the pcgb => pcgb isn't shaded
                    pc_abs_frac = 1._wp - exp(pc_abs_eff * lad_s(k,j,i))
                    pcbinsw(ipcgb) = pcbinsw(ipcgb) &
                        + rad_sw_in_dir(j, i) * pc_box_area * svf(2,isvf) * pc_abs_frac
                ENDIF
            ENDDO
        ENDIF
        surfins(startenergy:endenergy) = surfinswdir(startenergy:endenergy) + surfinswdif(startenergy:endenergy)
        surfinl(startenergy:endenergy) = surfinl(startenergy:endenergy) + surfinlwdif(startenergy:endenergy)
        surfinsw(:) = surfins(:)
        surfinlw(:) = surfinl(:)
        surfoutsw(:) = 0.0_wp
        surfoutlw(:) = surfoutll(:)
        surfhf(startenergy:endenergy) = surfinsw(startenergy:endenergy) + surfinlw(startenergy:endenergy) &
                                      - surfoutsw(startenergy:endenergy) - surfoutlw(startenergy:endenergy)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     Next passes - reflections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO refstep = 1, nrefsteps
        
            surfoutsl(startenergy:endenergy) = albedo_surf(startenergy:endenergy) * surfins(startenergy:endenergy)
!--         for non-transparent surfaces, longwave albedo is 1 - emissivity
            surfoutll(startenergy:endenergy) = (1._wp - emiss_surf(startenergy:endenergy)) * surfinl(startenergy:endenergy)

#if defined( __parallel )
            CALL MPI_AllGatherv(surfoutsl, nsurfl, MPI_REAL, &
                surfouts, nsurfs, surfstart, MPI_REAL, comm2d, ierr)
            CALL MPI_AllGatherv(surfoutll, nsurfl, MPI_REAL, &
                surfoutl, nsurfs, surfstart, MPI_REAL, comm2d, ierr)
#else
            surfouts(:) = surfoutsl(:)
            surfoutl(:) = surfoutll(:)
#endif

!--         reset for next pass input
            surfins(:) = 0._wp
            surfinl(:) = 0._wp
            
!--         reflected radiation
            DO isvf = 1, nsvfl
                isurf = svfsurf(1, isvf)
                isurfsrc = svfsurf(2, isvf)

!--             TODO: to remove if, use start+end for isvf
                IF ( surf(id, isurfsrc) < isky )  THEN
                    surfins(isurf) = surfins(isurf) + svf(1,isvf) * svf(2,isvf) * surfouts(isurfsrc)
                    surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutl(isurfsrc)
                ENDIF
            ENDDO

!--         radiation absorbed by plant canopy
            DO icsf = 1, ncsfl
                ipcgb = csfsurf(1, icsf)
                isurfsrc = csfsurf(2, icsf)

                IF ( surf(id, isurfsrc) < isky )  THEN
                    pcbinsw(ipcgb) = pcbinsw(ipcgb) + csf(1,icsf) * csf(2,icsf) * surfouts(isurfsrc)
!--                 pcbinlw(ipcgb) = pcbinlw(ipcgb) + csf(1,icsf) * surfoutl(isurfsrc)
                ENDIF
            ENDDO
            
            surfinsw(:) = surfinsw(:)  + surfins(:)
            surfinlw(:) = surfinlw(:)  + surfinl(:)
            surfoutsw(startenergy:endenergy) = surfoutsw(startenergy:endenergy) + surfoutsl(startenergy:endenergy)
            surfoutlw(startenergy:endenergy) = surfoutlw(startenergy:endenergy) + surfoutll(startenergy:endenergy)
            surfhf(startenergy:endenergy) = surfinsw(startenergy:endenergy) + surfinlw(startenergy:endenergy) &
                                          - surfoutsw(startenergy:endenergy) - surfoutlw(startenergy:endenergy)
        
        ENDDO

!--     push heat flux absorbed by plant canopy to respective 3D arrays
        IF ( plant_canopy )  THEN
            pc_heating_rate(:,:,:) = 0._wp
            DO ipcgb = 1, npcbl
                j = pcbl(iy, ipcgb)
                i = pcbl(ix, ipcgb)
                k = pcbl(iz, ipcgb)
                kk = k - nzb_s_inner(j,i)  !- lad arrays are defined flat
                pc_heating_rate(kk, j, i) = (pcbinsw(ipcgb) + pcbinlw(ipcgb)) &
                    * pchf_prep(k) * pt(k, j, i) !-- = dT/dt
            ENDDO
        ENDIF

!--     return surface radiation to horizontal surfaces 
!--     to rad_sw_in, rad_lw_in and rad_net for outputs
        !!!!!!!!!!
!--     we need the original radiation on urban top layer
!--     for calculation of MRT so we can't do adjustment here for now
        !!!!!!!!!!
        !!!DO isurf = 1, nsurfl
        !!!    i = surfl(ix,isurf)
        !!!    j = surfl(iy,isurf)
        !!!    k = surfl(iz,isurf)
        !!!    d = surfl(id,isurf)
        !!!    IF ( d==iroof )  THEN
        !!!        rad_sw_in(:,j,i) = surfinsw(isurf)
        !!!        rad_lw_in(:,j,i) = surfinlw(isurf)
        !!!        rad_net(j,i) = rad_sw_in(k,j,i) - rad_sw_out(k,j,i) + rad_lw_in(k,j,i) - rad_lw_out(k,j,i)
        !!!    ENDIF
        !!!ENDDO

    END SUBROUTINE usm_radiation

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Raytracing for detecting obstacles and calculating compound canopy sink
!> factors. (A simple obstacle detection would only need to process faces in
!> 3 dimensions without any ordering.)
!> Assumtions:
!> -----------
!> 1. The ray always originates from a face midpoint (only one coordinate equals
!>    *.5, i.e. wall) and doesn't travel parallel to the surface (that would mean
!>    shape factor=0). Therefore, the ray may never travel exactly along a face
!>    or an edge.
!> 2. From grid bottom to urban surface top the grid has to be *equidistant*
!>    within each of the dimensions, including vertical (but the resolution
!>    doesn't need to be the same in all three dimensions).
!------------------------------------------------------------------------------!
    SUBROUTINE usm_raytrace(src, targ, isrc, rirrf, atarg, create_csf, visible, transparency, win_lad)
        IMPLICIT NONE

        REAL(wp), DIMENSION(3), INTENT(in)     :: src, targ    !< real coordinates z,y,x
        INTEGER(iwp), INTENT(in)               :: isrc         !< index of source face for csf
        REAL(wp), INTENT(in)                   :: rirrf        !< irradiance factor for csf
        REAL(wp), INTENT(in)                   :: atarg        !< target surface area for csf
        LOGICAL, INTENT(in)                    :: create_csf   !< whether to generate new CSFs during raytracing
        LOGICAL, INTENT(out)                   :: visible
        REAL(wp), INTENT(out)                  :: transparency !< along whole path
        INTEGER(iwp), INTENT(in)               :: win_lad
        INTEGER(iwp)                           :: i, j, k, d
        INTEGER(iwp)                           :: seldim       !< dimension to be incremented
        INTEGER(iwp)                           :: ncsb         !< no of written plant canopy sinkboxes
        INTEGER(iwp)                           :: maxboxes     !< max no of gridboxes visited
        REAL(wp)                               :: distance     !< euclidean along path
        REAL(wp)                               :: crlen        !< length of gridbox crossing
        REAL(wp)                               :: lastdist     !< beginning of current crossing
        REAL(wp)                               :: nextdist     !< end of current crossing
        REAL(wp)                               :: realdist     !< distance in meters per unit distance
        REAL(wp)                               :: crmid        !< midpoint of crossing
        REAL(wp)                               :: cursink      !< sink factor for current canopy box
        REAL(wp), DIMENSION(3)                 :: delta        !< path vector
        REAL(wp), DIMENSION(3)                 :: uvect        !< unit vector
        REAL(wp), DIMENSION(3)                 :: dimnextdist  !< distance for each dimension increments
        INTEGER(iwp), DIMENSION(3)             :: box          !< gridbox being crossed
        INTEGER(iwp), DIMENSION(3)             :: dimnext      !< next dimension increments along path
        INTEGER(iwp), DIMENSION(3)             :: dimdelta     !< dimension direction = +- 1
        INTEGER(iwp)                           :: px, py       !< number of processors in x and y dir before 
                                                               !< the processor in the question
        INTEGER(iwp)                           :: ip           !< number of processor where gridbox reside
        INTEGER(iwp)                           :: ig           !< 1D index of gridbox in global 2D array
        REAL(wp)                               :: lad_s_target !< recieved lad_s of particular grid box
        REAL(wp), PARAMETER                    :: grow_factor = 1.5_wp !< factor of expansion of grow arrays

!--     Maximum number of gridboxes visited equals to maximum number of boundaries crossed in each dimension plus one. That's also
!--     the maximum number of plant canopy boxes written. We grow the acsf array accordingly using exponential factor.
        maxboxes = SUM(ABS(NINT(targ) - NINT(src))) + 1
        IF ( plant_canopy  .AND.  ncsfl + maxboxes > ncsfla )  THEN
!--         use this code for growing by fixed exponential increments (equivalent to case where ncsfl always increases by 1)
!--         k = CEILING(grow_factor ** real(CEILING(log(real(ncsfl + maxboxes, kind=wp)) &
!--                                                / log(grow_factor)), kind=wp))
!--         or use this code to simply always keep some extra space after growing
            k = CEILING(REAL(ncsfl + maxboxes, kind=wp) * grow_factor)

            CALL usm_merge_and_grow_csf(k)
        ENDIF
        
        transparency = 1._wp
        ncsb = 0

        delta(:) = targ(:) - src(:)
        distance = SQRT(SUM(delta(:)**2))
        IF ( distance == 0._wp )  THEN
            visible = .TRUE.
            RETURN
        ENDIF
        uvect(:) = delta(:) / distance
        realdist = SQRT(SUM( (uvect(:)*(/dz,dy,dx/))**2 ))

        lastdist = 0._wp

!--     Since all face coordinates have values *.5 and we'd like to use
!--     integers, all these have .5 added
        DO d = 1, 3
            IF ( uvect(d) == 0._wp )  THEN
                dimnext(d) = 999999999
                dimdelta(d) = 999999999
                dimnextdist(d) = 1.0E20_wp
            ELSE IF ( uvect(d) > 0._wp )  THEN
                dimnext(d) = CEILING(src(d) + .5_wp)
                dimdelta(d) = 1
                dimnextdist(d) = (dimnext(d) - .5_wp - src(d)) / uvect(d)
            ELSE
                dimnext(d) = FLOOR(src(d) + .5_wp)
                dimdelta(d) = -1
                dimnextdist(d) = (dimnext(d) - .5_wp - src(d)) / uvect(d)
            ENDIF
        ENDDO

        DO
!--         along what dimension will the next wall crossing be?
            seldim = minloc(dimnextdist, 1)
            nextdist = dimnextdist(seldim)
            IF ( nextdist > distance ) nextdist = distance

            crlen = nextdist - lastdist
            IF ( crlen > .001_wp )  THEN
                crmid = (lastdist + nextdist) * .5_wp
                box = NINT(src(:) + uvect(:) * crmid)

!--             calculate index of the grid with global indices (box(2),box(3))
!--             in the array nzterr and plantt and id of the coresponding processor
                px = box(3)/nnx
                py = box(2)/nny
                ip = px*pdims(2)+py
                ig = ip*nnx*nny + (box(3)-px*nnx)*nny + box(2)-py*nny
                IF ( box(1) <= nzterr(ig) )  THEN
                    visible = .FALSE.
                    RETURN
                ENDIF

                IF ( plant_canopy )  THEN
                    IF ( box(1) <= plantt(ig) )  THEN
                        ncsb = ncsb + 1
                        boxes(:,ncsb) = box
                        crlens(ncsb) = crlen
#if defined( __parallel )
                        lad_ip(ncsb) = ip
                        lad_disp(ncsb) = (box(3)-px*nnx)*(nny*nzu) + (box(2)-py*nny)*nzu + box(1)-nzub
#endif
                    ENDIF
                ENDIF
            ENDIF

            IF ( nextdist >= distance ) EXIT
            lastdist = nextdist
            dimnext(seldim) = dimnext(seldim) + dimdelta(seldim)
            dimnextdist(seldim) = (dimnext(seldim) - .5_wp - src(seldim)) / uvect(seldim)
        ENDDO
        
        IF ( plant_canopy )  THEN
#if defined( __parallel )
            IF ( usm_lad_rma )  THEN
!--             send requests for lad_s to appropriate processor
                CALL cpu_log( log_point_s(77), 'usm_init_rma', 'start' )
                DO i = 1, ncsb
                    CALL MPI_Get(lad_s_ray(i), 1, MPI_REAL, lad_ip(i), lad_disp(i), &
                                 1, MPI_REAL, win_lad, ierr)
                    IF ( ierr /= 0 )  THEN
                        WRITE(message_string, *) 'MPI error ', ierr, ' at MPI_Get'
                        CALL message( 'usm_raytrace', 'PA0519', 1, 2, 0, 6, 0 )
                    ENDIF
                ENDDO
                
!--             wait for all pending local requests complete
                CALL MPI_Win_flush_local_all(win_lad, ierr)
                IF ( ierr /= 0 )  THEN
                    WRITE(message_string, *) 'MPI error ', ierr, ' at MPI_Win_flush_local_all'
                    CALL message( 'usm_raytrace', 'PA0519', 1, 2, 0, 6, 0 )
                ENDIF
                CALL cpu_log( log_point_s(77), 'usm_init_rma', 'stop' )
                
            ENDIF
#endif

!--         calculate csf and transparency
            DO i = 1, ncsb
#if defined( __parallel )
                IF ( usm_lad_rma )  THEN
                    lad_s_target = lad_s_ray(i)
                ELSE
                    lad_s_target = usm_lad_g(lad_ip(i)*nnx*nny*nzu + lad_disp(i))
                ENDIF
#else
                lad_s_target = usm_lad(boxes(1,i),boxes(2,i),boxes(3,i))
#endif
                cursink = 1._wp - exp(-ext_coef * lad_s_target * crlens(i)*realdist)

                IF ( create_csf )  THEN
!--                 write svf values into the array
                    ncsfl = ncsfl + 1
                    acsf(ncsfl)%ip = lad_ip(i)
                    acsf(ncsfl)%itx = boxes(3,i)
                    acsf(ncsfl)%ity = boxes(2,i)
                    acsf(ncsfl)%itz = boxes(1,i)
                    acsf(ncsfl)%isurfs = isrc
                    acsf(ncsfl)%rsvf = REAL(cursink*rirrf*atarg, wp) !-- we postpone multiplication by transparency
                    acsf(ncsfl)%rtransp = REAL(transparency, wp)
                ENDIF  !< create_csf

                transparency = transparency * (1._wp - cursink)
                
            ENDDO
        ENDIF
        
        visible = .TRUE.
        
    END SUBROUTINE usm_raytrace
    
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This subroutine is part of the urban surface model.
!> It reads daily heat produced by anthropogenic sources
!> and the diurnal cycle of the heat.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_anthropogenic_heat
    
        INTEGER(iwp)                  :: i,j,ii
        REAL(wp)                      :: heat
        
!--     allocation of array of sources of anthropogenic heat and their diural profile
        ALLOCATE( aheat(nys:nyn,nxl:nxr) )
        ALLOCATE( aheatprof(0:24) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     read daily amount of heat and its daily cycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        aheat = 0.0_wp
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

!--             open anthropogenic heat file 
                OPEN( 151, file='ANTHROPOGENIC_HEAT'//TRIM(coupling_char), action='read', &
                           status='old', form='formatted', err=11 )
                i = 0
                j = 0
                DO
                    READ( 151, *, err=12, end=13 )  i, j, heat
                    IF ( i >= nxl  .AND.  i <= nxr  .AND.  j >= nys  .AND.  j <= nyn )  THEN
!--                     write heat into the array
                        aheat(j,i) = heat
                    ENDIF
                    CYCLE
 12                 WRITE(message_string,'(a,2i4)') 'error in file ANTHROPOGENIC_HEAT'//TRIM(coupling_char)//' after line ',i,j
                    CALL message( 'usm_read_anthropogenic_heat', 'PA0515', 0, 1, 0, 6, 0 )
                ENDDO
 13             CLOSE(151)
                CYCLE
 11             message_string = 'file ANTHROPOGENIC_HEAT'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_anthropogenic_heat', 'PA0516', 1, 2, 0, 6, 0 )
            ENDIF
            
#if defined( __parallel ) && ! defined ( __check )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     read diurnal profiles of heat sources
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        aheatprof = 0.0_wp
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

!--             open anthropogenic heat profile file 
                OPEN( 151, file='ANTHROPOGENIC_HEAT_PROFILE'//TRIM(coupling_char), action='read', &
                           status='old', form='formatted', err=21 )
                i = 0
                DO
                    READ( 151, *, err=22, end=23 )  i, heat
                    IF ( i >= 0  .AND.  i <= 24 )  THEN
!--                     write heat into the array
                        aheatprof(i) = heat
                    ENDIF
                    CYCLE
 22                 WRITE(message_string,'(a,i4)') 'error in file ANTHROPOGENIC_HEAT_PROFILE'// &
                                                     TRIM(coupling_char)//' after line ',i
                    CALL message( 'usm_read_anthropogenic_heat', 'PA0517', 0, 1, 0, 6, 0 )
                ENDDO
                aheatprof(24) = aheatprof(0)
 23             CLOSE(151)
                CYCLE
 21             message_string = 'file ANTHROPOGENIC_HEAT_PROFILE'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_anthropogenic_heat', 'PA0518', 1, 2, 0, 6, 0 )
            ENDIF
            
#if defined( __parallel ) && ! defined ( __check )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO
        
    END SUBROUTINE usm_read_anthropogenic_heat
   

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Soubroutine reads t_surf and t_wall data from restart files
!kanani: Renamed this routine according to corresponging routines in PALM
!kanani: Modified the routine to match read_var_list, from where usm_read_restart_data
!        shall be called in the future. This part has not been tested yet. (see virtual_flight_mod)
!        Also, I had some trouble with the allocation of t_surf, since this is a pointer.
!        So, I added some directives here.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_restart_data


       IMPLICIT NONE
       
       CHARACTER (LEN=30) ::  variable_chr  !< dummy variable to read string
       
       INTEGER(iwp)       ::  i             !< running index


       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN
             READ ( 13 )  variable_chr
             DO   WHILE ( TRIM( variable_chr ) /= '*** end usm ***' )

                SELECT CASE ( TRIM( variable_chr ) )
                
                   CASE ( 't_surf' )
#if defined( __nopointer )                   
                      IF ( .NOT.  ALLOCATED( t_surf ) )                         &
                         ALLOCATE( t_surf(startenergy:endenergy) )
                      READ ( 13 )  t_surf
#else                      
                      IF ( .NOT.  ALLOCATED( t_surf_1 ) )                         &
                         ALLOCATE( t_surf_1(startenergy:endenergy) )
                      READ ( 13 )  t_surf_1
#endif

                   CASE ( 't_wall' )
#if defined( __nopointer )
                      IF ( .NOT.  ALLOCATED( t_wall ) )                         &
                         ALLOCATE( t_wall(nzb_wall:nzt_wall+1,startenergy:endenergy) )
                      READ ( 13 )  t_wall
#else
                      IF ( .NOT.  ALLOCATED( t_wall_1 ) )                         &
                         ALLOCATE( t_wall_1(nzb_wall:nzt_wall+1,startenergy:endenergy) )
                      READ ( 13 )  t_wall_1
#endif

                   CASE DEFAULT
                      WRITE ( message_string, * )  'unknown variable named "', &
                                        TRIM( variable_chr ), '" found in',    &
                                        '&data from prior run on PE ', myid
                      CALL message( 'user_read_restart_data', 'UI0012', 1, 2, 0, 6, 0 )

                END SELECT

                READ ( 13 )  variable_chr

             ENDDO
          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

    END SUBROUTINE usm_read_restart_data


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Soubroutine reads svf and svfsurf data from saved file
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_svf_from_file

        IMPLICIT NONE
        INTEGER(iwp)                 :: fsvf = 89
        INTEGER(iwp)                 :: i
        CHARACTER(usm_version_len)   :: usm_version_field
        CHARACTER(svf_code_len)      :: svf_code_field

        DO  i = 0, io_blocks-1
            IF ( i == io_group )  THEN
                OPEN ( fsvf, file=TRIM(svf_file_name)//TRIM(coupling_char)//myid_char,               &
                    form='unformatted', status='old' )

!--             read and check version
                READ ( fsvf ) usm_version_field
                IF ( TRIM(usm_version_field) /= TRIM(usm_version) )  THEN
                    WRITE( message_string, * ) 'Version of binary SVF file "',           &
                                            TRIM(usm_version_field), '" does not match ',            &
                                            'the version of model "', TRIM(usm_version), '"'
                    CALL message( 'usm_read_svf_from_file', 'UI0012', 1, 2, 0, 6, 0 )
                ENDIF
                
!--             read nsvfl, ncsfl
                READ ( fsvf ) nsvfl, ncsfl
                IF ( nsvfl <= 0  .OR.  ncsfl < 0 )  THEN
                    WRITE( message_string, * ) 'Wrong number of SVF or CSF'
                    CALL message( 'usm_read_svf_from_file', 'UI0012', 1, 2, 0, 6, 0 )
                ELSE
                    WRITE(message_string,*) '    Number of SVF and CSF to read', nsvfl, ncsfl
                    CALL location_message( message_string, .TRUE. )
                ENDIF
                
                ALLOCATE(svf(ndsvf,nsvfl))
                ALLOCATE(svfsurf(idsvf,nsvfl))
                READ(fsvf) svf
                READ(fsvf) svfsurf
                IF ( plant_canopy )  THEN
                    ALLOCATE(csf(ndcsf,ncsfl))
                    ALLOCATE(csfsurf(idcsf,ncsfl))
                    READ(fsvf) csf
                    READ(fsvf) csfsurf
                ENDIF
                READ ( fsvf ) svf_code_field
                
                IF ( TRIM(svf_code_field) /= TRIM(svf_code) )  THEN
                    WRITE( message_string, * ) 'Wrong structure of binary svf file'
                    CALL message( 'usm_read_svf_from_file', 'UI0012', 1, 2, 0, 6, 0 )
                ENDIF
                
                CLOSE (fsvf)
                
            ENDIF
#if defined( __parallel )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO

    END SUBROUTINE usm_read_svf_from_file

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This subroutine reads walls, roofs and land categories and it parameters
!> from input files.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_read_urban_surface_types
    
        CHARACTER(12)                                         :: wtn
        INTEGER(iwp)                                          :: wtc
        REAL(wp), DIMENSION(n_surface_params)                 :: wtp
    
        INTEGER(iwp), DIMENSION(0:17, nysg:nyng, nxlg:nxrg)   :: usm_par
        REAL(wp), DIMENSION(1:14, nysg:nyng, nxlg:nxrg)       :: usm_val
        INTEGER(iwp)                                          :: k, l, d, iw, jw, kw, it, ip, ii, ij
        INTEGER(iwp)                                          :: i, j
        INTEGER(iwp)                                          :: nz, roof, dirwe, dirsn
        INTEGER(iwp)                                          :: category
        INTEGER(iwp)                                          :: weheight1, wecat1, snheight1, sncat1
        INTEGER(iwp)                                          :: weheight2, wecat2, snheight2, sncat2
        INTEGER(iwp)                                          :: weheight3, wecat3, snheight3, sncat3
        REAL(wp)                                              :: height, albedo, thick
        REAL(wp)                                              :: wealbedo1, wethick1, snalbedo1, snthick1
        REAL(wp)                                              :: wealbedo2, wethick2, snalbedo2, snthick2
        REAL(wp)                                              :: wealbedo3, wethick3, snalbedo3, snthick3
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     read categories of walls and their parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

!--             open urban surface file 
                OPEN( 151, file='SURFACE_PARAMETERS'//coupling_char, action='read', &
                           status='old', form='formatted', err=15 ) 
!--             first test and get n_surface_types
                k = 0
                l = 0
                DO
                    l = l+1
                    READ( 151, *, err=11, end=12 )  wtc, wtp, wtn
                    k = k+1
                    CYCLE
 11                 CONTINUE
                ENDDO
 12             n_surface_types = k
                ALLOCATE( surface_type_names(n_surface_types) )
                ALLOCATE( surface_type_codes(n_surface_types) )
                ALLOCATE( surface_params(n_surface_params, n_surface_types) )
!--             real reading
                rewind( 151 )
                k = 0
                DO
                    READ( 151, *, err=13, end=14 )  wtc, wtp, wtn
                    k = k+1
                    surface_type_codes(k) = wtc
                    surface_params(:,k) = wtp
                    surface_type_names(k) = wtn
                    CYCLE
13                  WRITE(6,'(i3,a,2i5)') myid, 'readparams2 error k=', k
                    FLUSH(6)
                    CONTINUE
                ENDDO
 14             CLOSE(151)
                CYCLE
 15             message_string = 'file SURFACE_PARAMETERS'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_urban_surface_types', 'PA0513', 1, 2, 0, 6, 0 )
            ENDIF
        ENDDO
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--     read types of surfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        usm_par = 0
        DO  ii = 0, io_blocks-1
            IF ( ii == io_group )  THEN

                !
!--             open csv urban surface file 
                OPEN( 151, file='URBAN_SURFACE'//TRIM(coupling_char), action='read', &
                      status='old', form='formatted', err=23 )
                
                l = 0
                DO
                    l = l+1
!--                 i, j, height, nz, roof, dirwe, dirsn, category, soilcat,
!--                 weheight1, wecat1, snheight1, sncat1, weheight2, wecat2, snheight2, sncat2,
!--                 weheight3, wecat3, snheight3, sncat3
                    READ( 151, *, err=21, end=25 )  i, j, height, nz, roof, dirwe, dirsn,            &
                                            category, albedo, thick,                                 &
                                            weheight1, wecat1, wealbedo1, wethick1,                  &
                                            weheight2, wecat2, wealbedo2, wethick2,                  &
                                            weheight3, wecat3, wealbedo3, wethick3,                  &
                                            snheight1, sncat1, snalbedo1, snthick1,                  &
                                            snheight2, sncat2, snalbedo2, snthick2,                  &
                                            snheight3, sncat3, snalbedo3, snthick3

                    IF ( i >= nxlg  .AND.  i <= nxrg  .AND.  j >= nysg  .AND.  j <= nyng )  THEN
!--                     write integer variables into array
                        usm_par(:,j,i) = (/1, nz, roof, dirwe, dirsn, category,                      &
                                          weheight1, wecat1, weheight2, wecat2, weheight3, wecat3,   &
                                          snheight1, sncat1, snheight2, sncat2, snheight3, sncat3 /)
!--                     write real values into array
                        usm_val(:,j,i) = (/ albedo, thick,                                           &
                                           wealbedo1, wethick1, wealbedo2, wethick2,                 &
                                           wealbedo3, wethick3, snalbedo1, snthick1,                 &
                                           snalbedo2, snthick2, snalbedo3, snthick3 /)
                    ENDIF
                    CYCLE
 21                 WRITE (message_string, "(A,I5)") 'errors in file URBAN_SURFACE'//TRIM(coupling_char)//' on line ', l
                    CALL message( 'usm_read_urban_surface_types', 'PA0512', 0, 1, 0, 6, 0 )
                ENDDO
          
 23             message_string = 'file URBAN_SURFACE'//TRIM(coupling_char)//' does not exist'
                CALL message( 'usm_read_urban_surface_types', 'PA0514', 1, 2, 0, 6, 0 )

 25             CLOSE( 90 )

            ENDIF
#if defined( __parallel ) && ! defined ( __check )
            CALL MPI_BARRIER( comm2d, ierr )
#endif
        ENDDO
        
        !
!--     check completeness and formal correctness of the data
        DO i = nxlg, nxrg
            DO j = nysg, nyng
                IF ( usm_par(0,j,i) /= 0  .AND.  (        &  !< incomplete data,supply default values later
                     usm_par(1,j,i) < nzb  .OR.           &
                     usm_par(1,j,i) > nzt  .OR.           &  !< incorrect height (nz < nzb  .OR.  nz > nzt)
                     usm_par(2,j,i) < 0  .OR.             &
                     usm_par(2,j,i) > 1  .OR.             &  !< incorrect roof sign
                     usm_par(3,j,i) < nzb-nzt  .OR.       & 
                     usm_par(3,j,i) > nzt-nzb  .OR.       &  !< incorrect west-east wall direction sign
                     usm_par(4,j,i) < nzb-nzt  .OR.       &
                     usm_par(4,j,i) > nzt-nzb  .OR.       &  !< incorrect south-north wall direction sign
                     usm_par(6,j,i) < nzb  .OR.           & 
                     usm_par(6,j,i) > nzt  .OR.           &  !< incorrect pedestrian level height for west-east wall
                     usm_par(8,j,i) > nzt  .OR.           &
                     usm_par(10,j,i) > nzt  .OR.          &  !< incorrect wall or roof level height for west-east wall
                     usm_par(12,j,i) < nzb  .OR.          & 
                     usm_par(12,j,i) > nzt  .OR.          &  !< incorrect pedestrian level height for south-north wall
                     usm_par(14,j,i) > nzt  .OR.          &
                     usm_par(16,j,i) > nzt                &  !< incorrect wall or roof level height for south-north wall
                    ) )  THEN
!--                 incorrect input data
                    WRITE (message_string, "(A,2I5)") 'missing or incorrect data in file URBAN_SURFACE'// &
                                                       TRIM(coupling_char)//' for i,j=', i,j
                    CALL message( 'usm_read_urban_surface', 'PA0504', 1, 2, 0, 6, 0 )
                ENDIF
                
            ENDDO
        ENDDO
        
!--     assign the surface types to local surface array
        DO  l = startenergy, endenergy
            
            d = surfl(id,l)
            kw = surfl(iz,l)
            j = surfl(iy,l)
            i = surfl(ix,l)
            IF ( d == iroof )  THEN
!--             horizontal surface - land or roof
                iw = i
                jw = j
                IF ( usm_par(5,jw,iw) == 0 )  THEN
                    IF ( zu(kw) >= roof_height_limit )  THEN
                        isroof_surf(l) = .TRUE.
                        surface_types(l) = roof_category         !< default category for root surface
                    ELSE
                        isroof_surf(l) = .FALSE.
                        surface_types(l) = land_category         !< default category for land surface
                    ENDIF
                    albedo_surf(l) = -1.0_wp
                    thickness_wall(l) = -1.0_wp
                ELSE
                    IF ( usm_par(2,jw,iw)==0 )  THEN
                        isroof_surf(l) = .FALSE.
                        thickness_wall(l) = -1.0_wp
                    ELSE
                        isroof_surf(l) = .TRUE.
                        thickness_wall(l) = usm_val(2,jw,iw)
                    ENDIF
                    surface_types(l) = usm_par(5,jw,iw)
                    albedo_surf(l) = usm_val(1,jw,iw)
                ENDIF
            ELSE
                SELECT CASE (d)
                    CASE (iwest)
                        iw = i
                        jw = j
                        ii = 6
                        ij = 3
                    CASE (ieast)
                        iw = i-1
                        jw = j
                        ii = 6
                        ij = 3
                    CASE (isouth)
                        iw = i
                        jw = j
                        ii = 12
                        ij = 9
                    CASE (inorth)
                        iw = i
                        jw = j-1
                        ii = 12
                        ij = 9
                END SELECT
                
                IF ( kw <= usm_par(ii,jw,iw) )  THEN
!--                 pedestrant zone
                    isroof_surf(l) = .FALSE.
                    IF ( usm_par(ii+1,jw,iw) == 0 )  THEN
                        surface_types(l) = pedestrant_category   !< default category for wall surface in pedestrant zone
                        albedo_surf(l) = -1.0_wp
                        thickness_wall(l) = -1.0_wp
                    ELSE
                        surface_types(l) = usm_par(ii+1,jw,iw)
                        albedo_surf(l) = usm_val(ij,jw,iw)
                        thickness_wall(l) = usm_val(ij+1,jw,iw)
                    ENDIF
                ELSE IF ( kw <= usm_par(ii+2,jw,iw) )  THEN
!--                 wall zone
                    isroof_surf(l) = .FALSE.
                    IF ( usm_par(ii+3,jw,iw) == 0 )  THEN
                        surface_types(l) = wall_category         !< default category for wall surface
                        albedo_surf(l) = -1.0_wp
                        thickness_wall(l) = -1.0_wp
                    ELSE
                        surface_types(l) = usm_par(ii+3,jw,iw)
                        albedo_surf(l) = usm_val(ij+2,jw,iw)
                        thickness_wall(l) = usm_val(ij+3,jw,iw)
                    ENDIF
                ELSE IF ( kw <= usm_par(ii+4,jw,iw) )  THEN
!--                 roof zone
                    isroof_surf(l) = .TRUE.
                    IF ( usm_par(ii+5,jw,iw) == 0 )  THEN
                        surface_types(l) = roof_category         !< default category for roof surface
                        albedo_surf(l) = -1.0_wp
                        thickness_wall(l) = -1.0_wp
                    ELSE
                        surface_types(l) = usm_par(ii+5,jw,iw)
                        albedo_surf(l) = usm_val(ij+4,jw,iw)
                        thickness_wall(l) = usm_val(ij+5,jw,iw)
                    ENDIF
                ELSE
!--                 something wrong
                    CALL message( 'usm_read_urban_surface', 'PA0505', 1, 2, 0, 6, 0 )
                ENDIF
            ENDIF
            
!--         find the type position
            it = surface_types(l)
            ip = -99999
            DO k = 1, n_surface_types
                IF ( surface_type_codes(k) == it )  THEN
                    ip = k
                    EXIT
                ENDIF
            ENDDO
            IF ( ip == -99999 )  THEN
!--             wall category not found
                WRITE (message_string, "(A,I5,A,3I5)") 'wall category ', it, ' not found  for i,j,k=', iw,jw,kw
                CALL message( 'usm_read_urban_surface', 'PA0506', 1, 2, 0, 6, 0 )
            ENDIF
            
!--         Fill out the parameters of the wall
!--         wall surface:
            
!--         albedo
            IF ( albedo_surf(l) < 0.0_wp )  THEN
                albedo_surf(l) = surface_params(ialbedo, ip)
            ENDIF
            
!--         emissivity of the wall
            emiss_surf(l) = surface_params(iemiss, ip)
            
!--         heat conductivity λS between air and wall ( W m−2 K−1 )
            lambda_surf(l) = surface_params(ilambdas, ip)
            
!--         roughness relative to concrete
            roughness_wall(l) = surface_params(irough, ip)
            
!--         Surface skin layer heat capacity (J m−2 K−1 )
            c_surface(l) = surface_params(icsurf, ip)
            
!--         wall material parameters:
            
!--         thickness of the wall (m)
!--         missing values are replaced by default value for category
            IF ( thickness_wall(l) <= 0.001_wp )  THEN
                thickness_wall(l) = surface_params(ithick, ip)
            ENDIF
            
!--         volumetric heat capacity rho_ocean*C of the wall ( J m−3 K−1 )
            rho_c_wall(:,l) = surface_params(irhoC, ip)
            
!--         thermal conductivity λH of the wall (W m−1 K−1 )
            lambda_h(:,l) = surface_params(ilambdah, ip)
            
        ENDDO

        CALL location_message( '    types and parameters of urban surfaces read', .TRUE. )
   
    END SUBROUTINE usm_read_urban_surface_types


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Solver for the energy balance at the ground/roof/wall surface.
!> It follows basic ideas and structure of lsm_energy_balance
!> with many simplifications and adjustments.
!> TODO better description
!------------------------------------------------------------------------------!
    SUBROUTINE usm_surface_energy_balance

        IMPLICIT NONE

        INTEGER(iwp)                          :: i, j, k, l, d      !< running indices
        
        REAL(wp)                              :: pt1                !< temperature at first grid box adjacent to surface
        REAL(wp)                              :: u1,v1,w1           !< near wall u,v,w
        REAL(wp)                              :: stend              !< surface tendency
        REAL(wp)                              :: coef_1             !< first coeficient for prognostic equation
        REAL(wp)                              :: coef_2             !< second  coeficient for prognostic equation
        REAL(wp)                              :: rho_cp             !< rho_wall_surface * cp
        REAL(wp)                              :: r_a                !< aerodynamic resistance for horizontal and vertical surfaces
        REAL(wp)                              :: f_shf              !< factor for shf_eb
        REAL(wp)                              :: lambda_surface     !< current value of lambda_surface (heat conductivity between air and wall)
        REAL(wp)                              :: Ueff               !< effective wind speed for calculation of heat transfer coefficients
        REAL(wp)                              :: httc               !< heat transfer coefficient
        REAL(wp), DIMENSION(nzub:nzut)        :: exn                !< value of the Exner function in layers
        
        REAL(wp), DIMENSION(0:4)              :: dxdir              !< surface normal direction gridbox length
        REAL(wp)                              :: dtime              !< simulated time of day (in UTC)
        INTEGER(iwp)                          :: dhour              !< simulated hour of day (in UTC)
        REAL(wp)                              :: acoef              !< actual coefficient of diurnal profile of anthropogenic heat

        dxdir = (/dz,dy,dy,dx,dx/)
        
        exn(:) = (hyp(nzub:nzut) / 100000.0_wp )**0.286_wp          !< Exner function
            
!--    
        DO l = startenergy, endenergy
!--         Calculate frequently used parameters
            d = surfl(id,l)
            k = surfl(iz,l)
            j = surfl(iy,l)
            i = surfl(ix,l)

!--         TODO - how to calculate lambda_surface for horizontal surfaces
!--         (lambda_surface is set according to stratification in land surface model)
            IF ( ol(j,i) >= 0.0_wp )  THEN
                lambda_surface = lambda_surf(l)
            ELSE
                lambda_surface = lambda_surf(l)
            ENDIF
            
            pt1  = pt(k,j,i)

!--         calculate rho_ocean * cp coefficient at surface layer
            rho_cp  = cp * hyp(k) / ( r_d * pt1 * exn(k) )

!--         calculate aerodyamic resistance. 
            IF ( d == iroof )  THEN
!--             calculation for horizontal surfaces follows LSM formulation
!--             pt, us, ts are not available for the prognostic time step,
!--             data from the last time step is used here.
                
                r_a = (pt1 - t_surf(l)/exn(k)) / (ts(j,i) * us(j,i) + 1.0E-10_wp)
                
!--             make sure that the resistance does not drop to zero
                IF ( ABS(r_a) < 1.0E-10_wp )  r_a = 1.0E-10_wp
                
!--             the parameterization is developed originally for larger scales
!--             (compare with remark in TUF-3D)
!--             our first experiences show that the parameterization underestimates
!--             r_a in meter resolution.
!--             temporary solution - multiplication by magic constant :-(.
                r_a = r_a * ra_horiz_coef
                
!--             factor for shf_eb
                f_shf  = rho_cp / r_a
            ELSE
!--             calculation of r_a for vertical surfaces
!--
!--             heat transfer coefficient for forced convection along vertical walls
!--             follows formulation in TUF3d model (Krayenhoff & Voogt, 2006)
!--            
!--             H = httc (Tsfc - Tair)
!--             httc = rw * (11.8 + 4.2 * Ueff) - 4.0
!--            
!--                   rw: wall patch roughness relative to 1.0 for concrete
!--                   Ueff: effective wind speed
!--                   - 4.0 is a reduction of Rowley et al (1930) formulation based on
!--                   Cole and Sturrock (1977)
!--            
!--                   Ucan: Canyon wind speed
!--                   wstar: convective velocity
!--                   Qs: surface heat flux
!--                   zH: height of the convective layer
!--                   wstar = (g/Tcan*Qs*zH)**(1./3.)
                
!--             staggered grid needs to be taken into consideration
                IF ( d == inorth )  THEN
                    u1 = (u(k,j,i)+u(k,j,i+1))*0.5_wp
                    v1 = v(k,j+1,i)
                ELSE IF ( d == isouth )  THEN
                    u1 = (u(k,j,i)+u(k,j,i+1))*0.5_wp
                    v1 = v(k,j,i)
                ELSE IF ( d == ieast )  THEN
                    u1 = u(k,j,i+1)
                    v1 = (v(k,j,i)+v(k,j+1,i))*0.5_wp
                ELSE IF ( d == iwest )  THEN
                    u1 = u(k,j,i)
                    v1 = (v(k,j,i)+v(k,j+1,i))*0.5_wp
                ELSE
                    STOP
                ENDIF
                w1 = (w(k,j,i)+w(k-1,j,i))*0.5_wp
                
                Ueff = SQRT(u1**2 + v1**2 + w1**2)
                httc = roughness_wall(l) * (11.8 + 4.2 * Ueff) - 4.0
                f_shf  = httc
            ENDIF
        
!--         add LW up so that it can be removed in prognostic equation
            rad_net_l(l) = surfinsw(l) - surfoutsw(l) + surfinlw(l) - surfoutlw(l)

!--         numerator of the prognostic equation
            coef_1 = rad_net_l(l) +    &    ! coef +1 corresponds to -lwout included in calculation of radnet_l
                     (3.0_wp+1.0_wp) * emiss_surf(l) * sigma_sb * t_surf(l) ** 4 +      &  
                     f_shf  * pt1 +                                                     &
                     lambda_surface * t_wall(nzb_wall,l)

!--         denominator of the prognostic equation
            coef_2 = 4.0_wp * emiss_surf(l) * sigma_sb * t_surf(l) ** 3                 &
                         + lambda_surface + f_shf / exn(k)

!--         implicit solution when the surface layer has no heat capacity,
!--         otherwise use RK3 scheme.
            t_surf_p(l) = ( coef_1 * dt_3d * tsc(2) + c_surface(l) * t_surf(l) ) /      & 
                              ( c_surface(l) + coef_2 * dt_3d * tsc(2) ) 

!--         add RK3 term
            t_surf_p(l) = t_surf_p(l) + dt_3d * tsc(3) * tt_surface_m(l)
            
!--         calculate true tendency
            stend = (t_surf_p(l) - t_surf(l) - dt_3d * tsc(3) * tt_surface_m(l)) / (dt_3d  * tsc(2))

!--         calculate t_surf tendencies for the next Runge-Kutta step
            IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                    tt_surface_m(l) = stend
                ELSEIF ( intermediate_timestep_count <                                  &
                         intermediate_timestep_count_max )  THEN
                    tt_surface_m(l) = -9.5625_wp * stend + 5.3125_wp                    &
                                       * tt_surface_m(l)
                ENDIF
            ENDIF

!--         in case of fast changes in the skin temperature, it is required to
!--         update the radiative fluxes in order to keep the solution stable
            IF ( ABS( t_surf_p(l) - t_surf(l) ) > 1.0_wp )  THEN
               force_radiation_call_l = .TRUE.
            ENDIF
            
!--         for horizontal surfaces is pt(nzb_s_inner(j,i),j,i) = pt_surf.
!--         there is no equivalent surface gridpoint for vertical surfaces.
!--         pt(k,j,i) is calculated for all directions in diffusion_s
!--         using surface and wall heat fluxes
            IF ( d == iroof )  THEN
               pt(nzb_s_inner(j,i),j,i) = t_surf_p(l) / exn(k)
            ENDIF

!--         calculate fluxes
!--         rad_net_l is never used!           
            rad_net_l(l)     = rad_net_l(l) + 3.0_wp * sigma_sb                         &
                                * t_surf(l)**4 - 4.0_wp * sigma_sb                      &
                                * t_surf(l)**3 * t_surf_p(l)
            wghf_eb(l)       = lambda_surface * (t_surf_p(l) - t_wall(nzb_wall,l))

!--         ground/wall/roof surface heat flux
            wshf_eb(l)  = - f_shf  * ( pt1 - t_surf_p(l) )
            
!--         store kinematic surface heat fluxes for utilization in other processes
!--         diffusion_s, surface_layer_fluxes,...
            IF ( d == iroof )  THEN
!--             shf is used in diffusion_s and also
!--             for calculation of surface layer fluxes
!--             update for horizontal surfaces
                shf(j,i) = wshf_eb(l) / rho_cp
            ELSE
!--             surface heat flux for vertical surfaces
!--             used in diffusion_s
                wshf(l) = wshf_eb(l) / rho_cp
            ENDIF

        ENDDO
        
        
        IF ( usm_anthropogenic_heat  .AND.  &
             intermediate_timestep_count == intermediate_timestep_count_max )  THEN
!--         application of the additional anthropogenic heat sources
!--         we considere the traffic for now so all heat is absorbed
!--         to the first layer, generalization would be worth
            
!--         calculation of actual profile coefficient
!--         ??? check time_since_reference_point ???
            dtime = mod(simulated_time + time_utc_init, 24.0_wp*3600.0_wp)
            dhour = INT(dtime/3600.0_wp)
!--         linear interpolation of coeficient
            acoef = (REAL(dhour+1,wp)-dtime/3600.0_wp)*aheatprof(dhour) + (dtime/3600.0_wp-REAL(dhour,wp))*aheatprof(dhour+1)
            DO i = nxl, nxr
                DO j = nys, nyn
                    IF ( aheat(j,i) > 0.0_wp )  THEN
!--                     TODO the increase of pt in box i,j,nzb_s_inner(j,i)+1 in time dt_3d 
!--                     given to anthropogenic heat aheat*acoef (W*m-2)
!--                     k = nzb_s_inner(j,i)+1
!--                     pt(k,j,i) = pt(k,j,i) + aheat(j,i)*acoef*dt_3d/(exn(k)*rho_cp*dz)
!--                     Instead of this, we can adjust shf in case AH only at surface
                        shf(j,i) = shf(j,i) + aheat(j,i)*acoef * ddx * ddy / rho_cp
                    ENDIF
                ENDDO
            ENDDO
        ENDIF
        
!--     pt and shf are defined on nxlg:nxrg,nysg:nyng
!--     get the borders from neighbours
        CALL exchange_horiz( pt, nbgp )
        CALL exchange_horiz_2d( shf )


!--    calculation of force_radiation_call:
!--    Make logical OR for all processes.
!--    Force radiation call if at least one processor forces it.
       IF ( intermediate_timestep_count == intermediate_timestep_count_max-1 )          &
       THEN
#if defined( __parallel )
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
              CALL mpi_allreduce( force_radiation_call_l, force_radiation_call,         &
                                  1, MPI_LOGICAL, MPI_LOR, comm2d, ierr )
#else
          force_radiation_call = force_radiation_call_l
#endif
          force_radiation_call_l = .FALSE.
       ENDIF

    END SUBROUTINE usm_surface_energy_balance


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels for t_surf and t_wall
!> called out from subroutine swap_timelevel
!------------------------------------------------------------------------------!
    SUBROUTINE usm_swap_timelevel ( mod_count )

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) :: mod_count
       INTEGER(iwp)             :: i
      
#if defined( __nopointer )
       t_surf    = t_surf_p
       t_wall    = t_wall_p
#else
       SELECT CASE ( mod_count )
          CASE ( 0 )
             t_surf  => t_surf_1; t_surf_p  => t_surf_2
             t_wall     => t_wall_1;    t_wall_p     => t_wall_2
          CASE ( 1 )
             t_surf  => t_surf_2; t_surf_p  => t_surf_1
             t_wall     => t_wall_2;    t_wall_p     => t_wall_1
       END SELECT
#endif
        
    END SUBROUTINE usm_swap_timelevel


!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This function applies the kinematic wall heat fluxes
!> for walls in four directions for all gridboxes in urban layer.
!> It is called out from subroutine prognostic_equations.
!> TODO Compare performance with cycle runnig l=startwall,endwall...
!------------------------------------------------------------------------------!
    SUBROUTINE usm_wall_heat_flux
    
        IMPLICIT NONE

        INTEGER(iwp)              ::  i,j,k,d,l             !< running indices
        
        DO l = startenergy, endenergy
            j = surfl(iy,l)
            i = surfl(ix,l)
            k = surfl(iz,l)
            d = surfl(id,l)
            tend(k,j,i) = tend(k,j,i) + wshf(l) * ddxy2(d)
        ENDDO

    END SUBROUTINE usm_wall_heat_flux
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> This function applies the kinematic wall heat fluxes
!> for walls in four directions around the gridbox i,j.
!> It is called out from subroutine prognostic_equations.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_wall_heat_flux_ij(i,j) 
    
        IMPLICIT NONE

        INTEGER(iwp), INTENT(in)  ::  i,j                   !< indices of grid box
        INTEGER(iwp)              ::  ii,jj,k,d,l
        
        DO l = startenergy, endenergy
            jj = surfl(iy,l)
            ii = surfl(ix,l)
            IF ( ii == i  .AND.  jj == j ) THEN
               k = surfl(iz,l)
               IF ( k >=  nzb_s_inner(j,i)+1  .AND.  k <=  nzb_s_outer(j,i) ) THEN
                  d = surfl(id,l)
                  IF ( d >= 1 .and. d <= 4 )   THEN
                     tend(k,j,i) = tend(k,j,i) + wshf(l) * ddxy2(d)
                  ENDIF
               ENDIF
            ENDIF
        ENDDO

    END SUBROUTINE usm_wall_heat_flux_ij
 

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine writes t_surf and t_wall data into restart files
!kanani: Renamed this routine according to corresponging routines in PALM
!kanani: Modified the routine to match write_var_list, from where usm_write_restart_data
!        shall be called in the future. This part has not been tested yet. (see virtual_flight_mod)
!        Also, I had some trouble with the allocation of t_surf, since this is a pointer.
!        So, I added some directives here.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_write_restart_data
    
       IMPLICIT NONE
       
       INTEGER(iwp)  ::  i

       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN
             WRITE ( 14 )  't_surf                        '
#if defined( __nopointer )             
             WRITE ( 14 )  t_surf
#else
             WRITE ( 14 )  t_surf_1
#endif
             WRITE ( 14 )  't_wall                        '
#if defined( __nopointer )             
             WRITE ( 14 )  t_wall
#else
             WRITE ( 14 )  t_wall_1
#endif
             WRITE ( 14 )  '*** end usm ***               '
          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

       
    END SUBROUTINE usm_write_restart_data


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine stores svf, svfsurf, csf and csfsurf data to a file.
!------------------------------------------------------------------------------!
    SUBROUTINE usm_write_svf_to_file

        IMPLICIT NONE
        INTEGER(iwp)        :: fsvf = 89
        INTEGER(iwp)        :: i

        DO  i = 0, io_blocks-1
            IF ( i == io_group )  THEN
                OPEN ( fsvf, file=TRIM(svf_file_name)//TRIM(coupling_char)//myid_char,               &
                    form='unformatted', status='new' )

                WRITE ( fsvf )  usm_version
                WRITE ( fsvf )  nsvfl, ncsfl
                WRITE ( fsvf )  svf
                WRITE ( fsvf )  svfsurf
                IF ( plant_canopy )  THEN
                    WRITE ( fsvf )  csf
                    WRITE ( fsvf )  csfsurf
                ENDIF
                WRITE ( fsvf )  TRIM(svf_code)

                CLOSE (fsvf)
#if defined( __parallel )
                CALL MPI_BARRIER( comm2d, ierr )
#endif
            ENDIF
        ENDDO
    END SUBROUTINE usm_write_svf_to_file


!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Block of auxiliary subroutines:
!> 1. quicksort and corresponding comparison
!> 2. usm_merge_and_grow_csf for implementation of "dynamical growing"
!>    array for csf
!------------------------------------------------------------------------------!    
    PURE FUNCTION svf_lt(svf1,svf2) result (res)
      TYPE (t_svf), INTENT(in) :: svf1,svf2
      LOGICAL                  :: res
      IF ( svf1%isurflt < svf2%isurflt  .OR.    &
          (svf1%isurflt == svf2%isurflt  .AND.  svf1%isurfs < svf2%isurfs) )  THEN
          res = .TRUE.
      ELSE
          res = .FALSE.
      ENDIF
    END FUNCTION svf_lt
    
 
!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_svf(svfl, first, last)
        IMPLICIT NONE
        TYPE(t_svf), DIMENSION(:), INTENT(INOUT)  :: svfl
        INTEGER(iwp), INTENT(IN)                  :: first, last
        TYPE(t_svf)                               :: x, t
        INTEGER(iwp)                              :: i, j

        IF ( first>=last ) RETURN
        x = svfl( (first+last) / 2 )
        i = first
        j = last
        DO
            DO while ( svf_lt(svfl(i),x) )
	        i=i+1
            ENDDO
            DO while ( svf_lt(x,svfl(j)) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t = svfl(i);  svfl(i) = svfl(j);  svfl(j) = t
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_svf(svfl, first, i-1)
        IF ( j+1 < last )  CALL quicksort_svf(svfl, j+1, last)
    END SUBROUTINE quicksort_svf

    
    PURE FUNCTION csf_lt(csf1,csf2) result (res)
      TYPE (t_csf), INTENT(in) :: csf1,csf2
      LOGICAL                  :: res
      IF ( csf1%ip < csf2%ip  .OR.    &
           (csf1%ip == csf2%ip  .AND.  csf1%itx < csf2%itx)  .OR.  &
           (csf1%ip == csf2%ip  .AND.  csf1%itx == csf2%itx  .AND.  csf1%ity < csf2%ity)  .OR.  &
           (csf1%ip == csf2%ip  .AND.  csf1%itx == csf2%itx  .AND.  csf1%ity == csf2%ity  .AND.   &
            csf1%itz < csf2%itz)  .OR.  &
           (csf1%ip == csf2%ip  .AND.  csf1%itx == csf2%itx  .AND.  csf1%ity == csf2%ity  .AND.   &
            csf1%itz == csf2%itz  .AND.  csf1%isurfs < csf2%isurfs) )  THEN
          res = .TRUE.
      ELSE
          res = .FALSE.
      ENDIF
    END FUNCTION csf_lt


!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_csf(csfl, first, last)
        IMPLICIT NONE
        TYPE(t_csf), DIMENSION(:), INTENT(INOUT)  :: csfl
        INTEGER(iwp), INTENT(IN)                  :: first, last
        TYPE(t_csf)                               :: x, t
        INTEGER(iwp)                              :: i, j

        IF ( first>=last ) RETURN
        x = csfl( (first+last)/2 )
        i = first
        j = last
        DO
            DO while ( csf_lt(csfl(i),x) )
                i=i+1
            ENDDO
            DO while ( csf_lt(x,csfl(j)) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t = csfl(i);  csfl(i) = csfl(j);  csfl(j) = t
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_csf(csfl, first, i-1)
        IF ( j+1 < last )  CALL quicksort_csf(csfl, j+1, last)
    END SUBROUTINE quicksort_csf

    
    SUBROUTINE usm_merge_and_grow_csf(newsize)
        INTEGER(iwp), INTENT(in)                :: newsize  !< new array size after grow, must be >= ncsfl
                                                            !< or -1 to shrink to minimum
        INTEGER(iwp)                            :: iread, iwrite
        TYPE(t_csf), DIMENSION(:), POINTER      :: acsfnew

        IF ( newsize == -1 )  THEN
!--         merge in-place
            acsfnew => acsf
        ELSE
!--         allocate new array
            IF ( mcsf == 0 )  THEN
                ALLOCATE( acsf1(newsize) )
                acsfnew => acsf1
            ELSE
                ALLOCATE( acsf2(newsize) )
                acsfnew => acsf2
            ENDIF
        ENDIF

        IF ( ncsfl >= 1 )  THEN
!--         sort csf in place (quicksort)
            CALL quicksort_csf(acsf,1,ncsfl)

!--         while moving to a new array, aggregate canopy sink factor records with identical box & source
            acsfnew(1) = acsf(1)
            iwrite = 1
            DO iread = 2, ncsfl
!--             here acsf(kcsf) already has values from acsf(icsf)
                IF ( acsfnew(iwrite)%itx == acsf(iread)%itx &
                         .AND.  acsfnew(iwrite)%ity == acsf(iread)%ity &
                         .AND.  acsfnew(iwrite)%itz == acsf(iread)%itz &
                         .AND.  acsfnew(iwrite)%isurfs == acsf(iread)%isurfs )  THEN
!--                 We could simply take either first or second rtransp, both are valid. As a very simple heuristic about which ray
!--                 probably passes nearer the center of the target box, we choose DIF from the entry with greater CSF, since that
!--                 might mean that the traced beam passes longer through the canopy box.
                    IF ( acsfnew(iwrite)%rsvf < acsf(iread)%rsvf )  THEN
                        acsfnew(iwrite)%rtransp = acsf(iread)%rtransp
                    ENDIF
                    acsfnew(iwrite)%rsvf = acsfnew(iwrite)%rsvf + acsf(iread)%rsvf
!--                 advance reading index, keep writing index
                ELSE
!--                 not identical, just advance and copy
                    iwrite = iwrite + 1
                    acsfnew(iwrite) = acsf(iread)
                ENDIF
            ENDDO
            ncsfl = iwrite
        ENDIF

        IF ( newsize == -1 )  THEN
!--         allocate new array and copy shrinked data
            IF ( mcsf == 0 )  THEN
                ALLOCATE( acsf1(ncsfl) )
                acsf1(1:ncsfl) = acsf2(1:ncsfl)
            ELSE
                ALLOCATE( acsf2(ncsfl) )
                acsf2(1:ncsfl) = acsf1(1:ncsfl)
            ENDIF
        ENDIF

!--     deallocate old array
        IF ( mcsf == 0 )  THEN
            mcsf = 1
            acsf => acsf1
            DEALLOCATE( acsf2 )
        ELSE
            mcsf = 0
            acsf => acsf2
            DEALLOCATE( acsf1 )
        ENDIF
        ncsfla = newsize
    END SUBROUTINE usm_merge_and_grow_csf

    
!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    RECURSIVE SUBROUTINE quicksort_csf2(kpcsflt, pcsflt, first, last)
        IMPLICIT NONE
        INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT)  :: kpcsflt
        REAL(wp), DIMENSION(:,:), INTENT(INOUT)      :: pcsflt
        INTEGER(iwp), INTENT(IN)                     :: first, last
        REAL(wp), DIMENSION(ndcsf)                   :: t2
        INTEGER(iwp), DIMENSION(kdcsf)               :: x, t1
        INTEGER(iwp)                                 :: i, j

        IF ( first>=last ) RETURN
        x = kpcsflt(:, (first+last)/2 )
        i = first
        j = last
        DO
            DO while ( csf_lt2(kpcsflt(:,i),x) )
                i=i+1
            ENDDO
            DO while ( csf_lt2(x,kpcsflt(:,j)) )
                j=j-1
            ENDDO
            IF ( i >= j ) EXIT
            t1 = kpcsflt(:,i);  kpcsflt(:,i) = kpcsflt(:,j);  kpcsflt(:,j) = t1
            t2 = pcsflt(:,i);  pcsflt(:,i) = pcsflt(:,j);  pcsflt(:,j) = t2
            i=i+1
            j=j-1
        ENDDO
        IF ( first < i-1 ) CALL quicksort_csf2(kpcsflt, pcsflt, first, i-1)
        IF ( j+1 < last )  CALL quicksort_csf2(kpcsflt, pcsflt, j+1, last)
    END SUBROUTINE quicksort_csf2
    

    PURE FUNCTION csf_lt2(item1, item2) result(res)
        INTEGER(iwp), DIMENSION(kdcsf), INTENT(in)  :: item1, item2
        LOGICAL                                     :: res
        res = ( (item1(3) < item2(3))                                                        &
             .OR.  (item1(3) == item2(3)  .AND.  item1(2) < item2(2))                            &
             .OR.  (item1(3) == item2(3)  .AND.  item1(2) == item2(2)  .AND.  item1(1) < item2(1)) &
             .OR.  (item1(3) == item2(3)  .AND.  item1(2) == item2(2)  .AND.  item1(1) == item2(1) &
                 .AND.  item1(4) < item2(4)) )
    END FUNCTION csf_lt2

    
 END MODULE urban_surface_mod
