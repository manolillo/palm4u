!> @file plant_canopy_model_mod.f90
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
! $Id: plant_canopy_model_mod.f90 2025 2016-10-12 16:44:03Z kanani $
!
! 2024 2016-10-12 16:42:37Z kanani
! Added missing lad_s initialization
! 
! 2011 2016-09-19 17:29:57Z kanani
! Renamed canopy_heat_flux to pc_heating_rate, since the original meaning/
! calculation of the quantity has changed, related to the urban surface model
! and similar future applications. 
! 
! 2007 2016-08-24 15:47:17Z kanani
! Added SUBROUTINE pcm_read_plant_canopy_3d for reading 3d plant canopy data
! from file (new case canopy_mode=read_from_file_3d) in the course of
! introduction of urban surface model,
! introduced variable ext_coef,
! resorted SUBROUTINEs to alphabetical order
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! 
! 1953 2016-06-21 09:28:42Z suehring
! Bugfix, lad_s and lad must be public
! 
! 1826 2016-04-07 12:01:39Z maronga
! Further modularization
!
! 1721 2015-11-16 12:56:48Z raasch
! bugfixes: shf is reduced in areas covered with canopy only,
!           canopy is set on top of topography
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   module plant_canopy_model_mod now contains a subroutine for the 
!   initialization of the canopy model (pcm_init),
!   limitation of the canopy drag (previously accounted for by calculation of 
!   a limiting canopy timestep for the determination of the maximum LES timestep
!   in subroutine timestep) is now realized by the calculation of pre-tendencies
!   and preliminary velocities in subroutine pcm_tendency,
!   some redundant MPI communication removed in subroutine pcm_init
!   (was previously in init_3d_model),
!   unnecessary 3d-arrays lad_u, lad_v, lad_w removed - lad information on the
!   respective grid is now provided only by lad_s (e.g. in the calculation of 
!   the tendency terms or of cum_lai_hf),
!   drag_coefficient, lai, leaf_surface_concentration, 
!   scalar_exchange_coefficient, sec and sls renamed to canopy_drag_coeff, 
!   cum_lai_hf, leaf_surface_conc, leaf_scalar_exch_coeff, lsec and lsc,
!   respectively,
!   unnecessary 3d-arrays cdc, lsc and lsec now defined as single-value constants,
!   USE-statements and ONLY-lists modified accordingly
! 
! 1340 2014-03-25 19:45:13Z kanani
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
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 138 2007-11-28 10:03:58Z letzel
! Initial revision
!
! Description:
! ------------
!> 1) Initialization of the canopy model, e.g. construction of leaf area density 
!> profile (subroutine pcm_init).
!> 2) Calculation of sinks and sources of momentum, heat and scalar concentration 
!> due to canopy elements (subroutine pcm_tendency).
!------------------------------------------------------------------------------!
 MODULE plant_canopy_model_mod
 
    USE arrays_3d,                                                             &
        ONLY:  dzu, dzw, e, q, s, shf, tend, u, v, w, zu, zw 

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxlu, nxr, nxrg, nyn, nyng, nys, nysg, nysv,   &
               nz, nzb, nzb_s_inner, nzb_u_inner, nzb_v_inner, nzb_w_inner, nzt

    USE kinds


    IMPLICIT NONE


    CHARACTER (LEN=20)   ::  canopy_mode = 'block' !< canopy coverage

    INTEGER(iwp) ::  pch_index = 0                 !< plant canopy height/top index
    INTEGER(iwp) ::                                                            &
       lad_vertical_gradient_level_ind(10) = -9999 !< lad-profile levels (index)

    LOGICAL ::  calc_beta_lad_profile = .FALSE. !< switch for calc. of lad from beta func.
    LOGICAL ::  plant_canopy = .FALSE.          !< switch for use of canopy model

    REAL(wp) ::  alpha_lad = 9999999.9_wp   !< coefficient for lad calculation
    REAL(wp) ::  beta_lad = 9999999.9_wp    !< coefficient for lad calculation
    REAL(wp) ::  canopy_drag_coeff = 0.0_wp !< canopy drag coefficient (parameter)
    REAL(wp) ::  cdc = 0.0_wp               !< canopy drag coeff. (abbreviation used in equations)
    REAL(wp) ::  cthf = 0.0_wp              !< canopy top heat flux
    REAL(wp) ::  dt_plant_canopy = 0.0_wp   !< timestep account. for canopy drag
    REAL(wp) ::  ext_coef = 0.6_wp          !< extinction coefficient
    REAL(wp) ::  lad_surface = 0.0_wp       !< lad surface value
    REAL(wp) ::  lai_beta = 0.0_wp          !< leaf area index (lai) for lad calc.
    REAL(wp) ::                                                                &
       leaf_scalar_exch_coeff = 0.0_wp      !< canopy scalar exchange coeff.
    REAL(wp) ::                                                                &
       leaf_surface_conc = 0.0_wp           !< leaf surface concentration
    REAL(wp) ::  lsec = 0.0_wp              !< leaf scalar exchange coeff.
    REAL(wp) ::  lsc = 0.0_wp               !< leaf surface concentration

    REAL(wp) ::                                                                &
       lad_vertical_gradient(10) = 0.0_wp              !< lad gradient
    REAL(wp) ::                                                                &
       lad_vertical_gradient_level(10) = -9999999.9_wp !< lad-prof. levels (in m)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  lad            !< leaf area density
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pre_lad        !< preliminary lad
    
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::                                 &
       pc_heating_rate                                    !< plant canopy heating rate
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  cum_lai_hf !< cumulative lai for heatflux calc.
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  lad_s      !< lad on scalar-grid


    SAVE


    PRIVATE
  
!
!-- Public functions
    PUBLIC pcm_check_parameters, pcm_header, pcm_init, pcm_parin, pcm_tendency

!
!-- Public variables and constants
    PUBLIC pc_heating_rate, canopy_mode, cthf, dt_plant_canopy, lad, lad_s,   &
           pch_index, plant_canopy
           


    INTERFACE pcm_check_parameters
       MODULE PROCEDURE pcm_check_parameters
    END INTERFACE pcm_check_parameters     
    
     INTERFACE pcm_header
       MODULE PROCEDURE pcm_header
    END INTERFACE pcm_header        
    
    INTERFACE pcm_init
       MODULE PROCEDURE pcm_init
    END INTERFACE pcm_init

    INTERFACE pcm_parin
       MODULE PROCEDURE pcm_parin
    END INTERFACE pcm_parin

    INTERFACE pcm_read_plant_canopy_3d
       MODULE PROCEDURE pcm_read_plant_canopy_3d
    END INTERFACE pcm_read_plant_canopy_3d
    
    INTERFACE pcm_tendency
       MODULE PROCEDURE pcm_tendency
       MODULE PROCEDURE pcm_tendency_ij
    END INTERFACE pcm_tendency


 CONTAINS

 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for plant canopy model
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_check_parameters

       USE control_parameters,                                                 &
           ONLY: cloud_physics, message_string, microphysics_seifert
                 
    
       IMPLICIT NONE

    
       IF ( canopy_drag_coeff == 0.0_wp )  THEN
          message_string = 'plant_canopy = .TRUE. requires a non-zero drag '// &
                           'coefficient & given value is canopy_drag_coeff = 0.0'
          CALL message( 'check_parameters', 'PA0041', 1, 2, 0, 6, 0 )
       ENDIF
    
       IF ( ( alpha_lad /= 9999999.9_wp  .AND.  beta_lad == 9999999.9_wp )  .OR.&
              beta_lad /= 9999999.9_wp   .AND.  alpha_lad == 9999999.9_wp )  THEN
          message_string = 'using the beta function for the construction ' //  &
                           'of the leaf area density profile requires '    //  &
                           'both alpha_lad and beta_lad to be /= 9999999.9'
          CALL message( 'check_parameters', 'PA0118', 1, 2, 0, 6, 0 )
       ENDIF
    
       IF ( calc_beta_lad_profile  .AND.  lai_beta == 0.0_wp )  THEN
          message_string = 'using the beta function for the construction ' //  &
                           'of the leaf area density profile requires '    //  &
                           'a non-zero lai_beta, but given value is '      //  &
                           'lai_beta = 0.0'
          CALL message( 'check_parameters', 'PA0119', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( calc_beta_lad_profile  .AND.  lad_surface /= 0.0_wp )  THEN
          message_string = 'simultaneous setting of alpha_lad /= 9999999.9' // &
                           'and lad_surface /= 0.0 is not possible, '       // &
                           'use either vertical gradients or the beta '     // &
                           'function for the construction of the leaf area '// &
                           'density profile'
          CALL message( 'check_parameters', 'PA0120', 1, 2, 0, 6, 0 )
       ENDIF 

       IF ( cloud_physics  .AND.  microphysics_seifert )  THEN
          message_string = 'plant_canopy = .TRUE. requires cloud_scheme /=' // &
                          ' seifert_beheng'
          CALL message( 'check_parameters', 'PA0360', 1, 2, 0, 6, 0 )
       ENDIF

 
    END SUBROUTINE pcm_check_parameters 
 

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for plant canopy model
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_header ( io )

       USE control_parameters,                                                 &
           ONLY: dz, passive_scalar


       IMPLICIT NONE
 
       CHARACTER (LEN=10) ::  coor_chr            !<

       CHARACTER (LEN=86) ::  coordinates         !<
       CHARACTER (LEN=86) ::  gradients           !<
       CHARACTER (LEN=86) ::  leaf_area_density   !<
       CHARACTER (LEN=86) ::  slices              !<
 
       INTEGER(iwp) :: i                !< 
       INTEGER(iwp),  INTENT(IN) ::  io !< Unit of the output file
       INTEGER(iwp) :: k                !<       
   
       REAL(wp) ::  canopy_height       !< canopy height (in m)
       
       canopy_height = pch_index * dz

       WRITE ( io, 1 )  canopy_mode, canopy_height, pch_index,                 &
                          canopy_drag_coeff
       IF ( passive_scalar )  THEN
          WRITE ( io, 2 )  leaf_scalar_exch_coeff,                             &
                             leaf_surface_conc
       ENDIF

!
!--    Heat flux at the top of vegetation
       WRITE ( io, 3 )  cthf

!
!--    Leaf area density profile, calculated either from given vertical 
!--    gradients or from beta probability density function.
       IF (  .NOT.  calc_beta_lad_profile )  THEN

!--       Building output strings, starting with surface value
          WRITE ( leaf_area_density, '(F7.4)' )  lad_surface
          gradients = '------'
          slices = '     0'
          coordinates = '   0.0'
          i = 1
          DO  WHILE ( i < 11  .AND.  lad_vertical_gradient_level_ind(i)        &
                      /= -9999 )

             WRITE (coor_chr,'(F7.2)')  lad(lad_vertical_gradient_level_ind(i))
             leaf_area_density = TRIM( leaf_area_density ) // ' ' //           &
                                 TRIM( coor_chr )
 
             WRITE (coor_chr,'(F7.2)')  lad_vertical_gradient(i)
             gradients = TRIM( gradients ) // ' ' // TRIM( coor_chr )

             WRITE (coor_chr,'(I7)')  lad_vertical_gradient_level_ind(i)
             slices = TRIM( slices ) // ' ' // TRIM( coor_chr )

             WRITE (coor_chr,'(F7.1)')  lad_vertical_gradient_level(i)
             coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

             i = i + 1
          ENDDO

          WRITE ( io, 4 )  TRIM( coordinates ), TRIM( leaf_area_density ),     &
                             TRIM( gradients ), TRIM( slices )

       ELSE
       
          WRITE ( leaf_area_density, '(F7.4)' )  lad_surface
          coordinates = '   0.0'
          
          DO  k = 1, pch_index

             WRITE (coor_chr,'(F7.2)')  lad(k)
             leaf_area_density = TRIM( leaf_area_density ) // ' ' //           &
                                 TRIM( coor_chr )
 
             WRITE (coor_chr,'(F7.1)')  zu(k)
             coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

          ENDDO       

          WRITE ( io, 5 ) TRIM( coordinates ), TRIM( leaf_area_density ),      &
                          alpha_lad, beta_lad, lai_beta

       ENDIF  

1 FORMAT (//' Vegetation canopy (drag) model:'/                                &
              ' ------------------------------'//                              &
              ' Canopy mode: ', A /                                            &
              ' Canopy height: ',F6.2,'m (',I4,' grid points)' /               &
              ' Leaf drag coefficient: ',F6.2 /)
2 FORMAT (/ ' Scalar exchange coefficient: ',F6.2 /                            &
              ' Scalar concentration at leaf surfaces in kg/m**3: ',F6.2 /)
3 FORMAT (' Predefined constant heatflux at the top of the vegetation: ',F6.2, &
          ' K m/s')
4 FORMAT (/ ' Characteristic levels of the leaf area density:'//               &
              ' Height:              ',A,'  m'/                                &
              ' Leaf area density:   ',A,'  m**2/m**3'/                        &
              ' Gradient:            ',A,'  m**2/m**4'/                        &
              ' Gridpoint:           ',A)
5 FORMAT (//' Characteristic levels of the leaf area density and coefficients:'&
          //  ' Height:              ',A,'  m'/                                &
              ' Leaf area density:   ',A,'  m**2/m**3'/                        &
              ' Coefficient alpha: ',F6.2 /                                    &
              ' Coefficient beta: ',F6.2 /                                     &
              ' Leaf area index: ',F6.2,'  m**2/m**2' /)   
       
    END SUBROUTINE pcm_header
 
 
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the plant canopy model
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_init
    

       USE control_parameters,                                                 &
           ONLY: coupling_char, dz, humidity, io_blocks, io_group,             &
                 message_string, ocean, passive_scalar 


       IMPLICIT NONE

       CHARACTER(10) :: pct
       
       INTEGER(iwp) ::  i   !< running index
       INTEGER(iwp) ::  ii  !< index       
       INTEGER(iwp) ::  j   !< running index
       INTEGER(iwp) ::  k   !< running index

       REAL(wp) ::  int_bpdf        !< vertical integral for lad-profile construction
       REAL(wp) ::  dzh             !< vertical grid spacing in units of canopy height
       REAL(wp) ::  gradient        !< gradient for lad-profile construction
       REAL(wp) ::  canopy_height   !< canopy height for lad-profile construction
       REAL(wp) ::  pcv(nzb:nzt+1)  !<
       
!
!--    Allocate one-dimensional arrays for the computation of the 
!--    leaf area density (lad) profile
       ALLOCATE( lad(0:nz+1), pre_lad(0:nz+1) )
       lad = 0.0_wp
       pre_lad = 0.0_wp

!
!--    Set flag that indicates that the lad-profile shall be calculated by using
!--    a beta probability density function
       IF ( alpha_lad /= 9999999.9_wp  .AND.  beta_lad /= 9999999.9_wp )  THEN
          calc_beta_lad_profile = .TRUE.
       ENDIF
       
       
!
!--    Compute the profile of leaf area density used in the plant
!--    canopy model. The profile can either be constructed from
!--    prescribed vertical gradients of the leaf area density or by
!--    using a beta probability density function (see e.g. Markkanen et al., 
!--    2003: Boundary-Layer Meteorology, 106, 437-459) 
       IF (  .NOT.  calc_beta_lad_profile )  THEN   

!
!--       Use vertical gradients for lad-profile construction    
          i = 1
          gradient = 0.0_wp

          IF (  .NOT.  ocean )  THEN

             lad(0) = lad_surface
             lad_vertical_gradient_level_ind(1) = 0
 
             DO k = 1, pch_index
                IF ( i < 11 )  THEN
                   IF ( lad_vertical_gradient_level(i) < zu(k)  .AND.          &
                        lad_vertical_gradient_level(i) >= 0.0_wp )  THEN
                      gradient = lad_vertical_gradient(i)
                      lad_vertical_gradient_level_ind(i) = k - 1
                      i = i + 1
                   ENDIF
                ENDIF
                IF ( gradient /= 0.0_wp )  THEN
                   IF ( k /= 1 )  THEN
                      lad(k) = lad(k-1) + dzu(k) * gradient
                   ELSE
                      lad(k) = lad_surface + dzu(k) * gradient
                   ENDIF
                ELSE
                   lad(k) = lad(k-1)
                ENDIF
             ENDDO

          ENDIF

!
!--       In case of no given leaf area density gradients, choose a vanishing
!--       gradient. This information is used for the HEADER and the RUN_CONTROL
!--       file.
          IF ( lad_vertical_gradient_level(1) == -9999999.9_wp )  THEN
             lad_vertical_gradient_level(1) = 0.0_wp
          ENDIF

       ELSE

! 
!--       Use beta function for lad-profile construction
          int_bpdf = 0.0_wp
          canopy_height = pch_index * dz

          DO k = nzb, pch_index
             int_bpdf = int_bpdf +                                             &
                      ( ( ( zw(k) / canopy_height )**( alpha_lad-1.0_wp ) ) *  &
                      ( ( 1.0_wp - ( zw(k) / canopy_height ) )**(              &
                          beta_lad-1.0_wp ) )                                  &
                      * ( ( zw(k+1)-zw(k) ) / canopy_height ) )
          ENDDO

!
!--       Preliminary lad profile (defined on w-grid)
          DO k = nzb, pch_index
             pre_lad(k) =  lai_beta *                                          &
                        ( ( ( zw(k) / canopy_height )**( alpha_lad-1.0_wp ) )  &
                        * ( ( 1.0_wp - ( zw(k) / canopy_height ) )**(          &
                              beta_lad-1.0_wp ) ) / int_bpdf                   &
                        ) / canopy_height
          ENDDO

!
!--       Final lad profile (defined on scalar-grid level, since most prognostic
!--       quantities are defined there, hence, less interpolation is required 
!--       when calculating the canopy tendencies)
          lad(0) = pre_lad(0)
          DO k = nzb+1, pch_index
             lad(k) = 0.5 * ( pre_lad(k-1) + pre_lad(k) )
          ENDDO          

       ENDIF

!
!--    Allocate 3D-array for the leaf area density (lad_s). In case of a 
!--    prescribed canopy-top heat flux (cthf), allocate 3D-arrays for
!--    the cumulative leaf area index (cum_lai_hf) and the canopy heat flux.
       ALLOCATE( lad_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

       IF ( cthf /= 0.0_wp )  THEN
          ALLOCATE( cum_lai_hf(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                 &
                    pc_heating_rate(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

!
!--    Initialize canopy parameters cdc (canopy drag coefficient), 
!--    lsec (leaf scalar exchange coefficient), lsc (leaf surface concentration)
!--    with the prescribed values
       cdc = canopy_drag_coeff
       lsec = leaf_scalar_exch_coeff
       lsc = leaf_surface_conc

!
!--    Initialization of the canopy coverage in the model domain:
!--    Setting the parameter canopy_mode = 'block' initializes a canopy, which
!--    fully covers the domain surface
       SELECT CASE ( TRIM( canopy_mode ) )

          CASE( 'block' )

             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   lad_s(:,j,i) = lad(:)
                ENDDO
             ENDDO

          CASE ( 'read_from_file_3d' )
!
!--          Initialize canopy parameters cdc (canopy drag coefficient),
!--          lsec (leaf scalar exchange coefficient), lsc (leaf surface concentration)
!--          from file which contains complete 3D data (separate vertical profiles for
!--          each location).
             CALL pcm_read_plant_canopy_3d

          CASE DEFAULT
!
!--          The DEFAULT case is reached either if the parameter 
!--          canopy mode contains a wrong character string or if the 
!--          user has coded a special case in the user interface. 
!--          There, the subroutine user_init_plant_canopy checks 
!--          which of these two conditions applies.
             CALL user_init_plant_canopy
 
       END SELECT

!
!--    Initialization of the canopy heat source distribution due to heating
!--    of the canopy layers by incoming solar radiation, in case that a non-zero
!--    value is set for the canopy top heat flux (cthf), which equals the
!--    available net radiation at canopy top.
!--    The heat source distribution is calculated by a decaying exponential 
!--    function of the downward cumulative leaf area index (cum_lai_hf), 
!--    assuming that the foliage inside the plant canopy is heated by solar
!--    radiation penetrating the canopy layers according to the distribution
!--    of net radiation as suggested by Brown & Covey (1966; Agric. Meteorol. 3, 
!--    73–96). This approach has been applied e.g. by Shaw & Schumann (1992;
!--    Bound.-Layer Meteorol. 61, 47–64)
       IF ( cthf /= 0.0_wp )  THEN
!
!--       Piecewise calculation of the cumulative leaf area index by vertical
!--       integration of the leaf area density
          cum_lai_hf(:,:,:) = 0.0_wp
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = pch_index-1, 0, -1
                   IF ( k == pch_index-1 )  THEN
                      cum_lai_hf(k,j,i) = cum_lai_hf(k+1,j,i) +                &
                         ( 0.5_wp * lad_s(k+1,j,i) *                           &
                           ( zw(k+1) - zu(k+1) ) )  +                          &
                         ( 0.5_wp * ( 0.5_wp * ( lad_s(k+1,j,i) +              &
                                                 lad_s(k,j,i) ) +              &
                                      lad_s(k+1,j,i) ) *                       &
                           ( zu(k+1) - zw(k) ) ) 
                   ELSE
                      cum_lai_hf(k,j,i) = cum_lai_hf(k+1,j,i) +                &
                         ( 0.5_wp * ( 0.5_wp * ( lad_s(k+2,j,i) +              &
                                                 lad_s(k+1,j,i) ) +            &
                                      lad_s(k+1,j,i) ) *                       &
                           ( zw(k+1) - zu(k+1) ) )  +                          &
                         ( 0.5_wp * ( 0.5_wp * ( lad_s(k+1,j,i) +              &
                                                 lad_s(k,j,i) ) +              &
                                      lad_s(k+1,j,i) ) *                       &
                           ( zu(k+1) - zw(k) ) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculation of the heating rate (K/s) within the different layers of 
!--       the plant canopy
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
!
!--             Calculation only necessary in areas covered with canopy
                IF ( cum_lai_hf(0,j,i) /= 0.0_wp )  THEN
!--             
!--                In areas with canopy the surface value of the canopy heat 
!--                flux distribution overrides the surface heat flux (shf) 
                   shf(j,i) = cthf * exp( -ext_coef * cum_lai_hf(0,j,i) )
!
!--                Within the different canopy layers the plant-canopy heating
!--                rate (pc_heating_rate) is calculated as the vertical 
!--                divergence of the canopy heat fluxes at the top and bottom
!--                of the respective layer
                   DO  k = 1, pch_index
                      pc_heating_rate(k,j,i) = cthf *                         &
                                                ( exp(-ext_coef*cum_lai_hf(k,j,i)) -  &
                                                  exp(-ext_coef*cum_lai_hf(k-1,j,i)) ) / dzw(k)
                   ENDDO
                ENDIF
             ENDDO
          ENDDO

       ENDIF



    END SUBROUTINE pcm_init


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &canopy_par for plant canopy model
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_parin


       IMPLICIT NONE

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file 
       
       NAMELIST /canopy_par/      alpha_lad, beta_lad, canopy_drag_coeff,      &
                                  canopy_mode, cthf,                           &
                                  lad_surface,                                 & 
                                  lad_vertical_gradient,                       &
                                  lad_vertical_gradient_level,                 &
                                  lai_beta,                                    &
                                  leaf_scalar_exch_coeff,                      &
                                  leaf_surface_conc, pch_index
       
       line = ' '
       
!
!--    Try to find radiation model package
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&canopy_par' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, canopy_par )

!
!--    Set flag that indicates that the radiation model is switched on
       plant_canopy = .TRUE.

 10    CONTINUE
       

    END SUBROUTINE pcm_parin



!------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Loads 3D plant canopy data from file. File format is as follows:
!>
!> num_levels
!> dtype,x,y,value(nzb),value(nzb+1), ... ,value(nzb+num_levels-1)
!> dtype,x,y,value(nzb),value(nzb+1), ... ,value(nzb+num_levels-1)
!> dtype,x,y,value(nzb),value(nzb+1), ... ,value(nzb+num_levels-1)
!> ...
!>
!> i.e. first line determines number of levels and further lines represent plant
!> canopy data, one line per column and variable. In each data line,
!> dtype represents variable to be set:
!>
!> dtype=1: leaf area density (lad_s)
!> dtype=2: canopy drag coefficient (cdc)
!> dtype=3: leaf scalar exchange coefficient (lsec)
!> dtype=4: leaf surface concentration (lsc)
!>
!> Zeros are added automatically above num_levels until top of domain.  Any
!> non-specified (x,y) columns have zero values as default.
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_read_plant_canopy_3d
        USE control_parameters, &
            ONLY: passive_scalar, message_string
        IMPLICIT NONE

        INTEGER(iwp)                            :: i, j, dtype, nzp, nzpltop, nzpl, kk
        REAL(wp), DIMENSION(:), ALLOCATABLE     :: col

        lad_s = 0.0_wp
        OPEN(152, file='PLANT_CANOPY_DATA_3D', access='SEQUENTIAL', &
                action='READ', status='OLD', form='FORMATTED', err=515)
        READ(152, *, err=516, end=517) nzp   !< read first line = number of vertical layers
        ALLOCATE(col(1:nzp))
        nzpltop = MIN(nzt+1, nzb+nzp-1)
        nzpl = nzpltop - nzb + 1    !< no. of layers to assign

        DO
            READ(152, *, err=516, end=517) dtype, i, j, col(:)
            IF ( i < nxlg .or. i > nxrg .or. j < nysg .or. j > nyng ) CYCLE

            SELECT CASE (dtype)
              CASE( 1 ) !< leaf area density
                !-- only lad_s has flat z-coordinate, others have regular
                kk = nzb_s_inner(j, i)
                lad_s(nzb:nzpltop-kk, j, i) = col(1+kk:nzpl)
!               CASE( 2 ) !< canopy drag coefficient
!                 cdc(nzb:nzpltop, j, i) = col(1:nzpl)
!               CASE( 3 ) !< leaf scalar exchange coefficient
!                 lsec(nzb:nzpltop, j, i) = col(1:nzpl)
!               CASE( 4 ) !< leaf surface concentration
!                 lsc(nzb:nzpltop, j, i) = col(1:nzpl)
              CASE DEFAULT
                write(message_string, '(a,i2,a)')   &
                    'Unknown record type in file PLANT_CANOPY_DATA_3D: "', dtype, '"'
                CALL message( 'pcm_read_plant_canopy_3d', 'PA0530', 1, 2, 0, 6, 0 )
            END SELECT
        ENDDO

515     message_string = 'error opening file PLANT_CANOPY_DATA_3D'
        CALL message( 'pcm_read_plant_canopy_3d', 'PA0531', 1, 2, 0, 6, 0 )

516     message_string = 'error reading file PLANT_CANOPY_DATA_3D'
        CALL message( 'pcm_read_plant_canopy_3d', 'PA0532', 1, 2, 0, 6, 0 )

517     CLOSE(152)
        DEALLOCATE(col)
        
    END SUBROUTINE pcm_read_plant_canopy_3d
    
    

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the tendency terms, accounting for the effect of the plant
!> canopy on momentum and scalar quantities.
!>
!> The canopy is located where the leaf area density lad_s(k,j,i) > 0.0 
!> (defined on scalar grid), as initialized in subroutine pcm_init. 
!> The lad on the w-grid is vertically interpolated from the surrounding
!> lad_s. The upper boundary of the canopy is defined on the w-grid at 
!> k = pch_index. Here, the lad is zero.
!>
!> The canopy drag must be limited (previously accounted for by calculation of 
!> a limiting canopy timestep for the determination of the maximum LES timestep
!> in subroutine timestep), since it is physically impossible that the canopy
!> drag alone can locally change the sign of a velocity component. This 
!> limitation is realized by calculating preliminary tendencies and velocities.
!> It is subsequently checked if the preliminary new velocity has a different
!> sign than the current velocity. If so, the tendency is limited in a way that
!> the velocity can at maximum be reduced to zero by the canopy drag. 
!>
!>
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_tendency( component )


       USE control_parameters,                                                 &
           ONLY:  dt_3d, message_string

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component !< prognostic variable (u,v,w,pt,q,e)
       INTEGER(iwp) ::  i         !< running index
       INTEGER(iwp) ::  j         !< running index
       INTEGER(iwp) ::  k         !< running index
       INTEGER(iwp) ::  kk        !< running index for flat lad arrays

       REAL(wp) ::  ddt_3d    !< inverse of the LES timestep (dt_3d)
       REAL(wp) ::  lad_local !< local lad value
       REAL(wp) ::  pre_tend  !< preliminary tendency
       REAL(wp) ::  pre_u     !< preliminary u-value 
       REAL(wp) ::  pre_v     !< preliminary v-value 
       REAL(wp) ::  pre_w     !< preliminary w-value 


       ddt_3d = 1.0_wp / dt_3d
 
!
!--    Compute drag for the three velocity components and the SGS-TKE:
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             DO  i = nxlu, nxr
                DO  j = nys, nyn
                   DO  k = nzb_u_inner(j,i)+1, nzb_u_inner(j,i)+pch_index

                      kk = k - nzb_u_inner(j,i)  !- lad arrays are defined flat
!
!--                   In order to create sharp boundaries of the plant canopy, 
!--                   the lad on the u-grid at index (k,j,i) is equal to 
!--                   lad_s(k,j,i), rather than being interpolated from the 
!--                   surrounding lad_s, because this would yield smaller lad 
!--                   at the canopy boundaries than inside of the canopy.
!--                   For the same reason, the lad at the rightmost(i+1)canopy 
!--                   boundary on the u-grid equals lad_s(k,j,i).
                      lad_local = lad_s(kk,j,i)
                      IF ( lad_local == 0.0_wp .AND. lad_s(kk,j,i-1) > 0.0_wp )&
                      THEN
                         lad_local = lad_s(kk,j,i-1)
                      ENDIF

                      pre_tend = 0.0_wp
                      pre_u = 0.0_wp
!
!--                   Calculate preliminary value (pre_tend) of the tendency
                      pre_tend = - cdc *                                       &
                                   lad_local *                                 &
                                   SQRT( u(k,j,i)**2 +                         &
                                         ( 0.25_wp * ( v(k,j,i-1) +            &
                                                       v(k,j,i)   +            &
                                                       v(k,j+1,i) +            &
                                                       v(k,j+1,i-1) )          &
                                         )**2 +                                &
                                         ( 0.25_wp * ( w(k-1,j,i-1) +          &
                                                       w(k-1,j,i)   +          &
                                                       w(k,j,i-1)   +          &
                                                       w(k,j,i) )              &
                                         )**2                                  &
                                       ) *                                     &
                                   u(k,j,i)

!
!--                   Calculate preliminary new velocity, based on pre_tend
                      pre_u = u(k,j,i) + dt_3d * pre_tend
!
!--                   Compare sign of old velocity and new preliminary velocity,
!--                   and in case the signs are different, limit the tendency
                      IF ( SIGN(pre_u,u(k,j,i)) /= pre_u )  THEN
                         pre_tend = - u(k,j,i) * ddt_3d
                      ELSE
                         pre_tend = pre_tend
                      ENDIF
!
!--                   Calculate final tendency
                      tend(k,j,i) = tend(k,j,i) + pre_tend

                   ENDDO
                ENDDO
             ENDDO

!
!--       v-component
          CASE ( 2 )
             DO  i = nxl, nxr
                DO  j = nysv, nyn
                   DO  k = nzb_v_inner(j,i)+1, nzb_v_inner(j,i)+pch_index

                      kk = k - nzb_v_inner(j,i)  !- lad arrays are defined flat
!
!--                   In order to create sharp boundaries of the plant canopy, 
!--                   the lad on the v-grid at index (k,j,i) is equal to 
!--                   lad_s(k,j,i), rather than being interpolated from the 
!--                   surrounding lad_s, because this would yield smaller lad 
!--                   at the canopy boundaries than inside of the canopy.
!--                   For the same reason, the lad at the northmost(j+1) canopy 
!--                   boundary on the v-grid equals lad_s(k,j,i).
                      lad_local = lad_s(kk,j,i)
                      IF ( lad_local == 0.0_wp .AND. lad_s(kk,j-1,i) > 0.0_wp )&
                      THEN
                         lad_local = lad_s(kk,j-1,i)
                      ENDIF

                      pre_tend = 0.0_wp
                      pre_v = 0.0_wp
!
!--                   Calculate preliminary value (pre_tend) of the tendency
                      pre_tend = - cdc *                                       &
                                   lad_local *                                 &
                                   SQRT( ( 0.25_wp * ( u(k,j-1,i)   +          &
                                                       u(k,j-1,i+1) +          &
                                                       u(k,j,i)     +          &
                                                       u(k,j,i+1) )            &
                                         )**2 +                                &
                                         v(k,j,i)**2 +                         &
                                         ( 0.25_wp * ( w(k-1,j-1,i) +          &
                                                       w(k-1,j,i)   +          &
                                                       w(k,j-1,i)   +          &
                                                       w(k,j,i) )              &
                                         )**2                                  &
                                       ) *                                     &
                                   v(k,j,i)

!
!--                   Calculate preliminary new velocity, based on pre_tend
                      pre_v = v(k,j,i) + dt_3d * pre_tend
!
!--                   Compare sign of old velocity and new preliminary velocity,
!--                   and in case the signs are different, limit the tendency
                      IF ( SIGN(pre_v,v(k,j,i)) /= pre_v )  THEN
                         pre_tend = - v(k,j,i) * ddt_3d
                      ELSE
                         pre_tend = pre_tend
                      ENDIF
!
!--                   Calculate final tendency
                      tend(k,j,i) = tend(k,j,i) + pre_tend

                   ENDDO
                ENDDO
             ENDDO

!
!--       w-component
          CASE ( 3 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_w_inner(j,i)+1, nzb_w_inner(j,i)+pch_index-1

                      kk = k - nzb_w_inner(j,i)  !- lad arrays are defined flat

                      pre_tend = 0.0_wp
                      pre_w = 0.0_wp
!
!--                   Calculate preliminary value (pre_tend) of the tendency
                      pre_tend = - cdc *                                       &
                                   (0.5_wp *                                   &
                                      ( lad_s(kk+1,j,i) + lad_s(kk,j,i) )) *   &
                                   SQRT( ( 0.25_wp * ( u(k,j,i)   +            &
                                                       u(k,j,i+1) +            &
                                                       u(k+1,j,i) +            &
                                                       u(k+1,j,i+1) )          &
                                         )**2 +                                &
                                         ( 0.25_wp * ( v(k,j,i)   +            &
                                                       v(k,j+1,i) +            &
                                                       v(k+1,j,i) +            &
                                                       v(k+1,j+1,i) )          &
                                         )**2 +                                &
                                         w(k,j,i)**2                           &
                                       ) *                                     &
                                   w(k,j,i)
!
!--                   Calculate preliminary new velocity, based on pre_tend
                      pre_w = w(k,j,i) + dt_3d * pre_tend
!
!--                   Compare sign of old velocity and new preliminary velocity,
!--                   and in case the signs are different, limit the tendency
                      IF ( SIGN(pre_w,w(k,j,i)) /= pre_w )  THEN
                         pre_tend = - w(k,j,i) * ddt_3d
                      ELSE
                         pre_tend = pre_tend
                      ENDIF
!
!--                   Calculate final tendency
                      tend(k,j,i) = tend(k,j,i) + pre_tend

                   ENDDO
                ENDDO
             ENDDO

!
!--       potential temperature
          CASE ( 4 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzb_s_inner(j,i)+pch_index
                      kk = k - nzb_s_inner(j,i)  !- lad arrays are defined flat
                      tend(k,j,i) = tend(k,j,i) + pc_heating_rate(kk,j,i)
                   ENDDO
                ENDDO
             ENDDO

!
!--       humidity
          CASE ( 5 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzb_s_inner(j,i)+pch_index
                      kk = k - nzb_s_inner(j,i)  !- lad arrays are defined flat
                      tend(k,j,i) = tend(k,j,i) -                              &
                                       lsec *                                  &
                                       lad_s(kk,j,i) *                         &
                                       SQRT( ( 0.5_wp * ( u(k,j,i) +           &
                                                          u(k,j,i+1) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( v(k,j,i) +           &
                                                          v(k,j+1,i) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( w(k-1,j,i) +         & 
                                                          w(k,j,i) )           &
                                             )**2                              &
                                           ) *                                 &
                                       ( q(k,j,i) - lsc )
                   ENDDO
                ENDDO
             ENDDO

!
!--       sgs-tke
          CASE ( 6 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzb_s_inner(j,i)+pch_index
                      kk = k - nzb_s_inner(j,i)  !- lad arrays are defined flat
                      tend(k,j,i) = tend(k,j,i) -                              &
                                       2.0_wp * cdc *                          &
                                       lad_s(kk,j,i) *                         &
                                       SQRT( ( 0.5_wp * ( u(k,j,i) +           &
                                                          u(k,j,i+1) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( v(k,j,i) +           &
                                                          v(k,j+1,i) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( w(k,j,i) +           &
                                                          w(k+1,j,i) )         &
                                             )**2                              &
                                           ) *                                 &
                                       e(k,j,i)
                   ENDDO
                ENDDO
             ENDDO 
!
!--       scalar concentration
          CASE ( 7 )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_s_inner(j,i)+1, nzb_s_inner(j,i)+pch_index
                      kk = k - nzb_s_inner(j,i)  !- lad arrays are defined flat
                      tend(k,j,i) = tend(k,j,i) -                              &
                                       lsec *                                  &
                                       lad_s(kk,j,i) *                         &
                                       SQRT( ( 0.5_wp * ( u(k,j,i) +           &
                                                          u(k,j,i+1) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( v(k,j,i) +           &
                                                          v(k,j+1,i) )         &
                                             )**2 +                            &
                                             ( 0.5_wp * ( w(k-1,j,i) +         & 
                                                          w(k,j,i) )           &
                                             )**2                              &
                                           ) *                                 &
                                       ( s(k,j,i) - lsc )
                   ENDDO
                ENDDO
             ENDDO    



          CASE DEFAULT

             WRITE( message_string, * ) 'wrong component: ', component
             CALL message( 'pcm_tendency', 'PA0279', 1, 2, 0, 6, 0 ) 

       END SELECT

    END SUBROUTINE pcm_tendency


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the tendency terms, accounting for the effect of the plant
!> canopy on momentum and scalar quantities.
!>
!> The canopy is located where the leaf area density lad_s(k,j,i) > 0.0 
!> (defined on scalar grid), as initialized in subroutine pcm_init. 
!> The lad on the w-grid is vertically interpolated from the surrounding
!> lad_s. The upper boundary of the canopy is defined on the w-grid at 
!> k = pch_index. Here, the lad is zero.
!>
!> The canopy drag must be limited (previously accounted for by calculation of 
!> a limiting canopy timestep for the determination of the maximum LES timestep
!> in subroutine timestep), since it is physically impossible that the canopy
!> drag alone can locally change the sign of a velocity component. This 
!> limitation is realized by calculating preliminary tendencies and velocities.
!> It is subsequently checked if the preliminary new velocity has a different
!> sign than the current velocity. If so, the tendency is limited in a way that
!> the velocity can at maximum be reduced to zero by the canopy drag. 
!>
!>
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE pcm_tendency_ij( i, j, component )


       USE control_parameters,                                                 &
           ONLY:  dt_3d, message_string

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  component !< prognostic variable (u,v,w,pt,q,e)
       INTEGER(iwp) ::  i         !< running index
       INTEGER(iwp) ::  j         !< running index
       INTEGER(iwp) ::  k         !< running index
       INTEGER(iwp) ::  kk        !< running index for flat lad arrays

       REAL(wp) ::  ddt_3d    !< inverse of the LES timestep (dt_3d)
       REAL(wp) ::  lad_local !< local lad value
       REAL(wp) ::  pre_tend  !< preliminary tendency
       REAL(wp) ::  pre_u     !< preliminary u-value 
       REAL(wp) ::  pre_v     !< preliminary v-value 
       REAL(wp) ::  pre_w     !< preliminary w-value 


       ddt_3d = 1.0_wp / dt_3d

!
!--    Compute drag for the three velocity components and the SGS-TKE
       SELECT CASE ( component )

!
!--       u-component
          CASE ( 1 )
             DO  k = nzb_u_inner(j,i)+1, nzb_u_inner(j,i)+pch_index

                kk = k - nzb_u_inner(j,i)  !- lad arrays are defined flat
!
!--             In order to create sharp boundaries of the plant canopy, 
!--             the lad on the u-grid at index (k,j,i) is equal to lad_s(k,j,i),
!--             rather than being interpolated from the surrounding lad_s, 
!--             because this would yield smaller lad at the canopy boundaries 
!--             than inside of the canopy.
!--             For the same reason, the lad at the rightmost(i+1)canopy 
!--             boundary on the u-grid equals lad_s(k,j,i).
                lad_local = lad_s(kk,j,i)
                IF ( lad_local == 0.0_wp .AND. lad_s(kk,j,i-1) > 0.0_wp )  THEN
                   lad_local = lad_s(kk,j,i-1)
                ENDIF

                pre_tend = 0.0_wp
                pre_u = 0.0_wp
!
!--             Calculate preliminary value (pre_tend) of the tendency
                pre_tend = - cdc *                                             &
                             lad_local *                                       &   
                             SQRT( u(k,j,i)**2 +                               &
                                   ( 0.25_wp * ( v(k,j,i-1)  +                 &
                                                 v(k,j,i)    +                 &
                                                 v(k,j+1,i)  +                 &
                                                 v(k,j+1,i-1) )                &
                                   )**2 +                                      &
                                   ( 0.25_wp * ( w(k-1,j,i-1) +                &
                                                 w(k-1,j,i)   +                &
                                                 w(k,j,i-1)   +                &
                                                 w(k,j,i) )                    &
                                   )**2                                        &
                                 ) *                                           &
                             u(k,j,i)

!
!--             Calculate preliminary new velocity, based on pre_tend
                pre_u = u(k,j,i) + dt_3d * pre_tend
!
!--             Compare sign of old velocity and new preliminary velocity,
!--             and in case the signs are different, limit the tendency
                IF ( SIGN(pre_u,u(k,j,i)) /= pre_u )  THEN
                   pre_tend = - u(k,j,i) * ddt_3d
                ELSE
                   pre_tend = pre_tend
                ENDIF
!
!--             Calculate final tendency
                tend(k,j,i) = tend(k,j,i) + pre_tend
             ENDDO


!
!--       v-component
          CASE ( 2 )
             DO  k = nzb_v_inner(j,i)+1, nzb_v_inner(j,i)+pch_index

                kk = k - nzb_v_inner(j,i)  !- lad arrays are defined flat
!
!--             In order to create sharp boundaries of the plant canopy, 
!--             the lad on the v-grid at index (k,j,i) is equal to lad_s(k,j,i),
!--             rather than being interpolated from the surrounding lad_s, 
!--             because this would yield smaller lad at the canopy boundaries 
!--             than inside of the canopy.
!--             For the same reason, the lad at the northmost(j+1)canopy 
!--             boundary on the v-grid equals lad_s(k,j,i).
                lad_local = lad_s(kk,j,i)
                IF ( lad_local == 0.0_wp .AND. lad_s(kk,j-1,i) > 0.0_wp )  THEN
                   lad_local = lad_s(kk,j-1,i)
                ENDIF

                pre_tend = 0.0_wp
                pre_v = 0.0_wp
!
!--             Calculate preliminary value (pre_tend) of the tendency
                pre_tend = - cdc *                                             &
                             lad_local *                                       &
                             SQRT( ( 0.25_wp * ( u(k,j-1,i)   +                &
                                                 u(k,j-1,i+1) +                &
                                                 u(k,j,i)     +                &
                                                 u(k,j,i+1) )                  &
                                   )**2 +                                      &
                                   v(k,j,i)**2 +                               &
                                   ( 0.25_wp * ( w(k-1,j-1,i) +                &
                                                 w(k-1,j,i)   +                &
                                                 w(k,j-1,i)   +                &
                                                 w(k,j,i) )                    &
                                   )**2                                        &
                                 ) *                                           &
                             v(k,j,i)

!
!--             Calculate preliminary new velocity, based on pre_tend
                pre_v = v(k,j,i) + dt_3d * pre_tend
!
!--             Compare sign of old velocity and new preliminary velocity,
!--             and in case the signs are different, limit the tendency
                IF ( SIGN(pre_v,v(k,j,i)) /= pre_v )  THEN
                   pre_tend = - v(k,j,i) * ddt_3d
                ELSE
                   pre_tend = pre_tend
                ENDIF
!
!--             Calculate final tendency
                tend(k,j,i) = tend(k,j,i) + pre_tend
             ENDDO


!
!--       w-component
          CASE ( 3 )
             DO  k = nzb_w_inner(j,i)+1, nzb_w_inner(j,i)+pch_index-1

                kk = k - nzb_w_inner(j,i)  !- lad arrays are defined flat

                pre_tend = 0.0_wp
                pre_w = 0.0_wp
!
!--             Calculate preliminary value (pre_tend) of the tendency
                pre_tend = - cdc *                                             &
                             (0.5_wp *                                         &
                                ( lad_s(kk+1,j,i) + lad_s(kk,j,i) )) *         &
                             SQRT( ( 0.25_wp * ( u(k,j,i)    +                 &  
                                                 u(k,j,i+1)  +                 &
                                                 u(k+1,j,i)  +                 &
                                                 u(k+1,j,i+1) )                &
                                   )**2 +                                      &
                                   ( 0.25_wp * ( v(k,j,i)    +                 &
                                                 v(k,j+1,i)  +                 &
                                                 v(k+1,j,i)  +                 &
                                                 v(k+1,j+1,i) )                &
                                   )**2 +                                      &
                                   w(k,j,i)**2                                 &
                                 ) *                                           &
                             w(k,j,i)
!
!--             Calculate preliminary new velocity, based on pre_tend
                pre_w = w(k,j,i) + dt_3d * pre_tend
!
!--             Compare sign of old velocity and new preliminary velocity,
!--             and in case the signs are different, limit the tendency
                IF ( SIGN(pre_w,w(k,j,i)) /= pre_w )  THEN
                   pre_tend = - w(k,j,i) * ddt_3d
                ELSE
                   pre_tend = pre_tend
                ENDIF
!
!--             Calculate final tendency
                tend(k,j,i) = tend(k,j,i) + pre_tend
             ENDDO

!
!--       potential temperature
          CASE ( 4 )
             DO  k = nzb_s_inner(j,i)+1, nzb_s_inner(j,i)+pch_index
                kk = k - nzb_s_inner(j,i)  !- lad arrays are defined flat
                tend(k,j,i) = tend(k,j,i) + pc_heating_rate(kk,j,i)
             ENDDO


!
!--       humidity
          CASE ( 5 )
             DO  k = nzb_s_inner(j,i)+1, nzb_s_inner(j,i)+pch_index
                kk = k - nzb_s_inner(j,i)  !- lad arrays are defined flat
                tend(k,j,i) = tend(k,j,i) -                                    &
                                 lsec *                                        &
                                 lad_s(kk,j,i) *                               &
                                 SQRT( ( 0.5_wp * ( u(k,j,i) +                 &
                                                    u(k,j,i+1) )               &
                                       )**2  +                                 &
                                       ( 0.5_wp * ( v(k,j,i) +                 &
                                                    v(k,j+1,i) )               &
                                       )**2 +                                  &
                                       ( 0.5_wp * ( w(k-1,j,i) +               &
                                                    w(k,j,i) )                 &
                                       )**2                                    &
                                     ) *                                       &
                                 ( q(k,j,i) - lsc )
             ENDDO   

!
!--       sgs-tke
          CASE ( 6 )
             DO  k = nzb_s_inner(j,i)+1, nzb_s_inner(j,i)+pch_index
                kk = k - nzb_s_inner(j,i)  !- lad arrays are defined flat
                tend(k,j,i) = tend(k,j,i) -                                    &
                                 2.0_wp * cdc *                                &
                                 lad_s(kk,j,i) *                               &
                                 SQRT( ( 0.5_wp * ( u(k,j,i) +                 &
                                                    u(k,j,i+1) )               &
                                       )**2 +                                  &  
                                       ( 0.5_wp * ( v(k,j,i) +                 &
                                                    v(k,j+1,i) )               &
                                       )**2 +                                  &
                                       ( 0.5_wp * ( w(k,j,i) +                 &
                                                    w(k+1,j,i) )               &
                                       )**2                                    &
                                     ) *                                       &
                                 e(k,j,i)
             ENDDO
             
!
!--       scalar concentration 
          CASE ( 7 )
             DO  k = nzb_s_inner(j,i)+1, nzb_s_inner(j,i)+pch_index
                kk = k - nzb_s_inner(j,i)  !- lad arrays are defined flat
                tend(k,j,i) = tend(k,j,i) -                                    &
                                 lsec *                                        &
                                 lad_s(kk,j,i) *                               &
                                 SQRT( ( 0.5_wp * ( u(k,j,i) +                 &
                                                    u(k,j,i+1) )               &
                                       )**2  +                                 &
                                       ( 0.5_wp * ( v(k,j,i) +                 &
                                                    v(k,j+1,i) )               &
                                       )**2 +                                  &
                                       ( 0.5_wp * ( w(k-1,j,i) +               &
                                                    w(k,j,i) )                 &
                                       )**2                                    &
                                     ) *                                       &
                                 ( s(k,j,i) - lsc )
             ENDDO                

       CASE DEFAULT

          WRITE( message_string, * ) 'wrong component: ', component
          CALL message( 'pcm_tendency', 'PA0279', 1, 2, 0, 6, 0 ) 

       END SELECT

    END SUBROUTINE pcm_tendency_ij



 END MODULE plant_canopy_model_mod
