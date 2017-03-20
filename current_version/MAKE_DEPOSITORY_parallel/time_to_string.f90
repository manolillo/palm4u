!> @file time_to_string.f90
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
! $Id: time_to_string.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1342 2014-03-26 17:04:47Z kanani
! REAL constants defined as wp-kind
! 
! 1320 2014-03-20 08:40:49Z raasch
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
! Revision 1.1  1997/08/11 06:26:08  raasch
! Initial revision
!
!
! Description:
! ------------
!> Transforming the time from real to character-string hh:mm:ss
!------------------------------------------------------------------------------!
 FUNCTION time_to_string( time )
 

    USE kinds

    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string !< 

    INTEGER(iwp)      ::  hours   !< 
    INTEGER(iwp)      ::  minutes !< 
    INTEGER(iwp)      ::  seconds !< 

    REAL(wp)          ::  rest_time !< 
    REAL(wp)          ::  time      !< 

!
!-- Calculate the number of hours, minutes, and seconds
    hours     = INT( time / 3600.0_wp )
    rest_time = time - hours * 3600_wp
    minutes   = INT( rest_time / 60.0_wp )
    seconds   = rest_time - minutes * 60

!
!-- Build the string
    IF ( hours < 100 )  THEN
       WRITE (time_to_string,'(I2.2,'':'',I2.2,'':'',I2.2)')  hours, minutes, &
                                                              seconds
    ELSE
       WRITE (time_to_string,'(I3.3,'':'',I2.2,'':'',I2.2)')  hours, minutes, &
                                                              seconds
    ENDIF

 END FUNCTION time_to_string
