!> @file lpm_write_exchange_statistics.f90
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
! $Id: lpm_write_exchange_statistics.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! Unused variables removed.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated. 
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 849 2012-03-15 10:35:09Z raasch
! initial revision (former part of advec_particles)
!
!
! Description:
! ------------
!> Write particle statistics (total particle numbers and number of particles
!> exchanged between subdomains) on ASCII file.
!>
!> @attention Output format of this file could be further improved! At current
!>            stage it is only a test output.
!------------------------------------------------------------------------------!
 SUBROUTINE lpm_write_exchange_statistics
 

    USE control_parameters,                                                    &
        ONLY:  current_timestep_number, dt_3d, simulated_time

    USE indices,                                                               &
        ONLY:  nxl, nxr, nys, nyn, nzb, nzt

    USE particle_attributes,                                                   &
        ONLY:  grid_particles, maximum_number_of_particles,                    &
               number_of_particles, prt_count,                                 &
               trlp_count_sum, trlp_count_recv_sum, trnp_count_sum,            &
               trnp_count_recv_sum, trrp_count_sum, trrp_count_recv_sum,       &
               trsp_count_sum, trsp_count_recv_sum

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) :: ip         !<
    INTEGER(iwp) :: jp         !<
    INTEGER(iwp) :: kp         !<

!
!-- Determine maximum number of particles (i.e., all possible particles that 
!-- have been allocated) and the current number of particles
    number_of_particles         = 0
    maximum_number_of_particles = 0
    DO  ip = nxl, nxr
       DO  jp = nys, nyn
          DO  kp = nzb+1, nzt
             number_of_particles = number_of_particles                         &
                                     + prt_count(kp,jp,ip)
             maximum_number_of_particles = maximum_number_of_particles         &
                                     + SIZE(grid_particles(kp,jp,ip)%particles)
          ENDDO
       ENDDO
    ENDDO

    CALL check_open( 80 )
#if defined( __parallel )
    WRITE ( 80, 8000 )  current_timestep_number+1, simulated_time+dt_3d, &
                        number_of_particles, pleft, trlp_count_sum,      &
                        trlp_count_recv_sum, pright, trrp_count_sum,     &
                        trrp_count_recv_sum, psouth, trsp_count_sum,     &
                        trsp_count_recv_sum, pnorth, trnp_count_sum,     &
                        trnp_count_recv_sum, maximum_number_of_particles
    CALL close_file( 80 )
#else
    WRITE ( 80, 8000 )  current_timestep_number+1, simulated_time+dt_3d, &
                        number_of_particles, maximum_number_of_particles
#endif

!
!-- Formats
8000 FORMAT (I6,1X,F7.2,4X,I10,5X,4(I3,1X,I4,'/',I4,2X),6X,I10)


 END SUBROUTINE lpm_write_exchange_statistics
