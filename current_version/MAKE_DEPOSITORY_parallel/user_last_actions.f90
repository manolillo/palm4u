!> @file user_last_actions.f90
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
! $Id: user_last_actions.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
!
! 1320 2014-03-20 08:40:49Z raasch
! revision history before 2012 removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 211 2008-11-11 04:46:24Z raasch
! Former file user_interface.f90 split into one file per subroutine
!
! Description:
! ------------
!> Execution of user-defined actions at the end of a job.
!------------------------------------------------------------------------------!
 SUBROUTINE user_last_actions
 

    USE control_parameters
        
    USE kinds
    
    USE user

    IMPLICIT NONE

!
!-- Here the user-defined actions at the end of a job follow.
!-- Sample for user-defined output:
    IF ( write_binary(1:4) == 'true' )  THEN
!       IF ( ALLOCATED( u2_av ) )  THEN
!          WRITE ( 14 )  'u2_av               ';  WRITE ( 14 )  u2_av
!       ENDIF

       WRITE ( 14 )  '*** end user ***    '

    ENDIF

 END SUBROUTINE user_last_actions

