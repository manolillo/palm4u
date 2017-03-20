!> @file kchem_driver.f90
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
! kk: Intial version (Klaus Ketelsen)
! bK: Changed initial values of chem_species
! RFo: Added tmp_fact
! FKa: Some formatting, todos added below, "_wp" added to REAL quantities' values,
!      the "if(var(1:3) == 'kc_')" is now included in check_parameters (as done
!      for "usm_" quantities
!
! Former revisions:
! -----------------
!
! Description:
! ------------
!> Interface PALM <-> kpp chemistry
!>
!> @todo Format this module according to PALM coding standard (see e.g. module
!>       template under http://palm.muk.uni-hannover.de/mosaik/downloads/8 or
!>       D3_coding_standard.pdf under https://palm.muk.uni-hannover.de/trac/downloads/16)
!> @todo Rename ssws, sswst to csws, cswst to avoid confusion with ssws/sswst 
!>       used together with passive scalar s
!> @todo Use same run index for all loops over NSPEC, e.g. nsp = 1, NSPEC
!>       (and: NSPEC-->nspec, i.e. no capital letters in variable names 
!>       [see formatting rules])
!> @todo "#ifdef KPP_CHEM" --> "#if defined( __chem)" (see other examples 
!>       in PALM)               #if defined( __parallel )
!> @todo „write(6...“ printing should eventually be handled in a separate 
!>       subroutine „kchem_header“ (see e.g. how it is done in routine 
!>       lsm_header in land_surface_model_mod.f90), which prints this 
!>       information to a MONITORING file named <jobname>_header (in this file, 
!>       all PALM setup information are included)
!> @todo Not all routines have meaningful names yet, e.g. "kchem_integrate" only
!>       refers to chemical reactions part, but the name itself could also refer
!>       to the other part of the chemistry integration (advection/diffusion)
!>       --> Please clarify
!> @todo Rename flag „use_kpp_chemistry“ to „chemistry“ (see e.g. flag 
!>       „plant_canopy“ used for the plant canopy model
!> @todo Please enter all of your TODOs here in this list
!> @todo Every subroutine needs a "Description" header 
!>       (as e.g. in urban_surface_mod.f90)
!> @todo Clear circular dependencies between check_parameters, kchem_initialize,
!>       and kchem_check_data_output

!------------------------------------------------------------------------------!


MODULE kchem_driver
   USE kinds,              ONLY: wp, iwp
   USE indices,            ONLY: nzb,nzt,nysg,nyng,nxlg,nxrg,nzb_s_inner,nys,nyn
   USE pegrid,             ONLY: myid, threads_per_task
   USE control_parameters, ONLY: dt_3d, ws_scheme_sca  !ws_sch... added by bK
   USE arrays_3d,          ONLY: pt
   USE kchem_kpp,          ONLY: NSPEC, SPC_NAMES, NKPPCTRL, NMAXFIXSTEPS, t_steps, FILL_TEMP, kpp_integrate,     &
                                 NVAR, ATOL, RTOL
   USE cpulog,             ONLY: cpu_log, log_point_s

   IMPLICIT   none
   PRIVATE
   SAVE

!- Define chemical variables

   TYPE   species_def
      CHARACTER(LEN=8)                                   :: name
      CHARACTER(LEN=16)                                  :: unit
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: conc
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: conc_p
      REAL(kind=wp),POINTER,DIMENSION(:,:,:)             :: tconc_m
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:)           :: ssws, sswst
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:)           :: flux_s, diss_s
      REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:)         :: flux_l, diss_l
   END TYPE species_def

   PUBLIC  species_def

   logical, PUBLIC                                               :: use_kpp_chemistry = .FALSE.

   TYPE(species_def),ALLOCATABLE,DIMENSION(:),TARGET, PUBLIC     :: chem_species

   REAL(kind=wp),ALLOCATABLE,DIMENSION(:,:,:,:),TARGET   :: spec_conc_1,spec_conc_2,spec_conc_3

   INTEGER,DIMENSION(NKPPCTRL)                           :: icntrl                            ! Fine tuning kpp
   REAL(kind=wp),DIMENSION(NKPPCTRL)                     :: rcntrl                            ! Fine tuning kpp

   PUBLIC NSPEC
   PUBLIC NVAR       ! added NVAR for pe  bK, kd3
   PUBLIC SPC_NAMES  ! added for pe bk, kd4
!- Interface section

   INTERFACE kchem_initialize
      MODULE PROCEDURE kchem_initialize
   END INTERFACE kchem_initialize

   INTERFACE kchem_parin
      MODULE PROCEDURE kchem_parin
   END INTERFACE kchem_parin

   INTERFACE kchem_integrate
      MODULE PROCEDURE kchem_integrate_ij
   END INTERFACE kchem_integrate

   INTERFACE kchem_swap_timelevel
      MODULE PROCEDURE kchem_swap_timelevel
   END INTERFACE kchem_swap_timelevel

   INTERFACE kchem_define_netcdf_grid
      MODULE PROCEDURE kchem_define_netcdf_grid
   END INTERFACE kchem_define_netcdf_grid

   INTERFACE kchem_data_output_3d
      MODULE PROCEDURE kchem_data_output_3d
   END INTERFACE kchem_data_output_3d

   INTERFACE kchem_check_data_output
      MODULE PROCEDURE kchem_check_data_output
   END INTERFACE kchem_check_data_output


   PUBLIC kchem_check_data_output, kchem_data_output_3d,                       &
          kchem_define_netcdf_grid, kchem_initialize, kchem_integrate,         &
          kchem_parin, kchem_swap_timelevel
            


 CONTAINS

   SUBROUTINE kchem_initialize
      IMPLICIT   none
!--   local variables
      INTEGER                  :: i

!--   Allocate Memory for chemical species

      ALLOCATE(chem_species(NSPEC))
      ALLOCATE(spec_conc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC))
      ALLOCATE(spec_conc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC))
      ALLOCATE(spec_conc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg,NSPEC))

      if(myid == 0)  then
         write(6,'(/,a,/)')  'kpp chemics >>>> List of species: '
      end if

      do i=1,NSPEC
         chem_species(i)%name = SPC_NAMES(i)
         if(myid == 0)  then
            write(6,'(a,i4,3x,a)')  '   Species: ',i,trim(SPC_NAMES(i))
         end if
         chem_species(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_1(:,:,:,i)
         chem_species(i)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_2(:,:,:,i)
         chem_species(i)%tconc_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) => spec_conc_3(:,:,:,i)

         ALLOCATE (chem_species(i)%ssws(nysg:nyng,nxlg:nxrg))
         ALLOCATE (chem_species(i)%sswst(nysg:nyng,nxlg:nxrg))

         chem_species(i)%ssws  = 0.0_wp
         chem_species(i)%sswst = 0.0_wp


!          IF ( ws_scheme_sca )  THEN                                           !fk (needs to be revisited. due to
            ALLOCATE (chem_species(i)%flux_s(nzb+1:nzt,0:threads_per_task-1))   !circular dependencies between check_parameters,
            ALLOCATE (chem_species(i)%diss_s(nzb+1:nzt,0:threads_per_task-1))   !kchem_initialize, kchem_check_data_output,
            ALLOCATE (chem_species(i)%flux_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1)) !this IF clause creates problems,
            ALLOCATE (chem_species(i)%diss_l(nzb+1:nzt,nys:nyn,0:threads_per_task-1)) !since ws_scheme_sca is first set "true"  
            chem_species(i)%flux_s = 0.0_wp                          !bK        !in check_parameters, which is so far called   
            chem_species(i)%flux_l = 0.0_wp                                     !after kchem_initialize
            chem_species(i)%diss_s = 0.0_wp
            chem_species(i)%diss_l = 0.0_wp
!          ENDIF         
      end do

!--   Set initial values

      CALL set_const_initial_values

      if(myid == 0)  write(6,*) ' '

      RETURN

     CONTAINS
      SUBROUTINE set_const_initial_values
         IMPLICIT   none

!--      local variables
         INTEGER                  :: i

         if(myid == 0)  then
            write(6,'(/,a,/)')  'kpp chemics >>>> Set constant Initial Values: '
         end if

!        Default values are taken from smog.def from supplied kpp example
         do i=1,NSPEC
            if(trim(chem_species(i)%name) == 'NO')   then
!              chem_species(i)%conc = 8.725*1.0E+08
               chem_species(i)%conc = 0.5_wp                           !added by bK
            else if(trim(chem_species(i)%name) == 'NO2') then
!              chem_species(i)%conc = 2.240*1.0E+08             
               chem_species(i)%conc = 0.1_wp                           !added by bK
            else if(trim(chem_species(i)%name) == 'O3') then
               chem_species(i)%conc = 0.05_wp                          !added by bK
            else if(trim(chem_species(i)%name) == 'H2O') then
!              chem_species(i)%conc = 5.326*1.0E+11
               chem_species(i)%conc = 1.30*1.0E+4_wp                   !added by bK
            else if(trim(chem_species(i)%name) == 'O2') then
!              chem_species(i)%conc = 1.697*1.0E+16
               chem_species(i)%conc = 2.0*1.0E+5_wp                    !added by bK
            else if(trim(chem_species(i)%name) == 'RH') then
!               chem_species(i)%conc = 9.906*1.0E+01
                chem_species(i)%conc = 2.0_wp                          !added by bK
            else if(trim(chem_species(i)%name) == 'RCHO') then
!              chem_species(i)%conc = 6.624*1.0E+08
               chem_species(i)%conc = 2.0_wp                           !added by bK
            else if(trim(chem_species(i)%name) == 'OH') then
               chem_species(i)%conc = 1.0*1.0E-07_wp                   !added by bK
            else if(trim(chem_species(i)%name) == 'HO2') then
               chem_species(i)%conc = 1*1.0E-7_wp                      !added by bK
            else if(trim(chem_species(i)%name) == 'RCCO2') then
               chem_species(i)%conc = 1.0*1.0E-7_wp                    !added by bK
            else if(trim(chem_species(i)%name) == 'RCOO2NO2') then
               chem_species(i)%conc = 1.0*1.0E-7_wp                   !added by bK
           else
!   H2O = 2.0e+04;
               chem_species(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg) = 0.0_wp
            end if
!             chem_species(i)%conc_p  = 0.0_wp
            chem_species(i)%conc_p  = chem_species(i)%conc       ! added bK
            chem_species(i)%tconc_m = 0.0_wp

            if(myid == 0)  then
               write(6,'(a,3x,a,3x,a,e12.4)')  '   Species:     ',chem_species(i)%name(1:7),' Initial Value = ',chem_species(i)%conc(nzb,nysg,nxlg)
            end if
         end do

#if defined( __nopointer )
!kk      Hier mit message abbrechen
         if(myid == 0)  then
            write(6,*)  '   KPP does only run with POINTER Version'
         end if
         stop 'error'
#endif

         RETURN
      END SUBROUTINE set_const_initial_values
   END SUBROUTINE kchem_initialize

   SUBROUTINE kchem_parin
      IMPLICIT   none

      CHARACTER (LEN=80)                             ::  line                   ! < dummy string that contains the current line of the parameter file
      REAL(kind=wp), DIMENSION(NMAXFIXSTEPS)         ::  my_steps               ! List of fixed timesteps   my_step(1) = 0.0 automatic stepping

      NAMELIST /kpp_chem/ icntrl, rcntrl, my_steps

!--   Read kpp_chem namelist
      icntrl    = 0
      rcntrl    = 0.0_wp
      my_steps  = 0.0_wp
      line      = ' '
      icntrl(2) = 1                                   !Atol and Rtol are Scalar

      ATOL = 1.0_wp
      RTOL = 0.01_wp
!
!--   Try to find kpp chemistry package
      REWIND ( 11 )
      line = ' '
      DO   WHILE ( INDEX( line, '&kpp_chem' ) == 0 )
        READ ( 11, '(A)', END=10 )  line
      ENDDO
      BACKSPACE ( 11 )

!
!--   Read user-defined namelist
      READ ( 11, kpp_chem )

      use_kpp_chemistry = .TRUE.

 10   CONTINUE

      t_steps = my_steps

      if(myid <= 1)  then
         write(6,*) 'KPP Parin ',icntrl(2),icntrl(3),ATOL(1)
      end if


      RETURN
   END SUBROUTINE kchem_parin

   SUBROUTINE kchem_integrate_ij (i, j)
      IMPLICIT   none
      INTEGER,INTENT(IN)       :: i,j

!--   local variables
      INTEGER                  :: k,m,istat
      INTEGER,dimension(20)    :: istatus
      REAL(kind=wp),dimension(nzb_s_inner(j,i)+1:nzt,NSPEC)                :: tmp_conc           !bK NSPEC repl with NVAR kd6
      REAL(kind=wp),dimension(nzb_s_inner(j,i)+1:nzt)                      :: tmp_temp
      REAL(kind=wp),dimension(nzb_s_inner(j,i)+1:nzt)                      :: tmp_fact   !! RFo

! ppm to molecules/cm**3
       tmp_fact = 10e-6_iwp*6.022e23_wp/22.414_wp/1000._wp     !*T/273.15*1013/p
!      tmp_temp = 292
!      CALL fill_temp (istat, tmp_temp)                             ! Load constant temperature into kpp context
      CALL fill_temp (istat, pt(nzb_s_inner(j,i)+1:nzt,j,i))                             ! Load temperature into kpp context


      do m=1,NSPEC
         tmp_conc(:,m) = chem_species(m)%conc (nzb_s_inner(j,i)+1:nzt,j,i) * tmp_fact(:) ! RFo
      end do

      if(myid == 0 .and. i == 10 .and. j == 10)  then
         write(0,*) 'begin KPP step ',dt_3d
      end if

      CALL cpu_log( log_point_s(80), 'kpp_integrate', 'start' )

      CALL kpp_integrate (dt_3d, tmp_conc, istatus=istatus)

      CALL cpu_log( log_point_s(80), 'kpp_integrate', 'stop' )

      do m=1,NSPEC
         chem_species(m)%conc (nzb_s_inner(j,i)+1:nzt,j,i) = tmp_conc(:,m) / tmp_fact(:)  ! RFo
      end do

      if(myid == 0 .and. i == 10 .and. j == 10)  then
         write(6,'(a,8i7)') ' KPP Status ',istatus(1:8)
      end if

      RETURN
   END SUBROUTINE kchem_integrate_ij

   SUBROUTINE kchem_swap_timelevel (level)
      IMPLICIT   none

      INTEGER,INTENT(IN)                  :: level

!--   local variables

      INTEGER               :: i

      if(level == 0)  then
         do i=1,NVAR                                        ! NSPEC replaced with NVAR bK kd1  
            chem_species(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_1(:,:,:,i)
            chem_species(i)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_2(:,:,:,i)
         end do
      else
         do i=1,NVAR                                        ! NSPEC replaced with NVAR bk  kd2
            chem_species(i)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_2(:,:,:,i)
            chem_species(i)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_1(:,:,:,i)
         end do
      end if

      RETURN
   END SUBROUTINE kchem_swap_timelevel

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!------------------------------------------------------------------------------!
   SUBROUTINE kchem_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

      IMPLICIT NONE

      CHARACTER (LEN=*), INTENT(IN)  ::  var         !<
      LOGICAL, INTENT(OUT)           ::  found       !<
      CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !<
      CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !<
      CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !<

      found  = .TRUE.

      if(var(1:3) == 'kc_')   then                    ! always the same grid for chemistry variables
            grid_x = 'x'
            grid_y = 'y'
            grid_z = 'zu'                             !kk Use same z axis as u variables. Has to be checked if OK
      else
            found  = .FALSE.
            grid_x = 'none'
            grid_y = 'none'
            grid_z = 'none'
      end if

      write(6,*) 'kchem_define_netcdf_grid ',TRIM(var),' ',trim(grid_x),' ',found

   END SUBROUTINE kchem_define_netcdf_grid

   SUBROUTINE kchem_check_data_output( var, unit, i, ilen, k )


      USE control_parameters,                                                 &
         ONLY: data_output, message_string

      IMPLICIT NONE

      CHARACTER (LEN=*) ::  unit     !<
      CHARACTER (LEN=*) ::  var      !<

      INTEGER(iwp) :: i
      INTEGER(iwp) :: ilen
      INTEGER(iwp) :: k

      INTEGER              :: n
      CHARACTER(len=16)    ::  spec_name

      unit = 'illegal'

      spec_name = TRIM(var(4:))

      do n=1,NSPEC
         if(TRIM(spec_name) == TRIM(chem_species(n)%name))   Then
            unit = 'ppm'
         end if
      end do

      RETURN
   END SUBROUTINE kchem_check_data_output

   SUBROUTINE kchem_data_output_3d( av, variable, found, local_pf )


      USE indices

      USE kinds


      IMPLICIT NONE

      CHARACTER (LEN=*) ::  variable !<
      LOGICAL      ::  found !<
      INTEGER(iwp) ::  av    !<
      REAL(sp), DIMENSION(nxlg:nxrg,nysg:nyng,nzb:nzt+1) ::  local_pf !<

      !-- local variables

      INTEGER              ::  i, j, k,n
      CHARACTER(len=16)    ::  spec_name


      found = .FALSE.

      spec_name = TRIM(variable(4:))

      do n=1,NSPEC
         if(TRIM(spec_name) == TRIM(chem_species(n)%name))   Then
            write(6,*) 'Output of species ',TRIM(variable),' ',TRIM(chem_species(n)%name)
            DO  i = nxlg, nxrg
               DO  j = nysg, nyng
                  DO  k = nzb, nzt+1
                     local_pf(i,j,k) = chem_species(n)%conc(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
            found = .TRUE.
         end if
      end do


      RETURN
   END SUBROUTINE kchem_data_output_3d


END MODULE kchem_driver
