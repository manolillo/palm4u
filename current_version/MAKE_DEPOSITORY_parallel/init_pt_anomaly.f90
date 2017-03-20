!> @file init_pt_anomaly.f90
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
! $Id: init_pt_anomaly.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp_kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 861 2012-03-26 14:18:34Z suehring
! Modification of the amplitude to obtain a visible temperature perturbation.
!
! Revision 1.1  1997/08/29 08:58:56  raasch
! Initial revision
!
!
! Description:
! ------------
!> Impose a temperature perturbation for an advection test.
!------------------------------------------------------------------------------!
 SUBROUTINE init_pt_anomaly
 

    USE arrays_3d,                                                             &
        ONLY:  pt, zu

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxr, nyn, nys, nzb, nzt
        
    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  ic !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  jc !<
    INTEGER(iwp) ::  k  !<
    INTEGER(iwp) ::  kc !<
    
    REAL(wp)     ::  betrag !<
    REAL(wp)     ::  radius !<
    REAL(wp)     ::  rc     !<
    REAL(wp)     ::  x      !<
    REAL(wp)     ::  y      !<
    REAL(wp)     ::  z      !<
    
!
!-- Defaults: radius rc, strength z,
!--           position of centre: ic, jc, kc
    rc =  10.0_wp * dx
    ic =  ( nx+1 ) / 2
    jc =  ic
    kc =  nzt / 2

!
!-- Compute the perturbation.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             x = ( i - ic ) * dx
             y = ( j - jc ) * dy
             z = ABS( zu(k) - zu(kc) )
             radius = SQRT( x**2 + y**2 + z**2 )
             IF ( radius <= rc )  THEN
                betrag = 5.0_wp * EXP( -( radius * 0.001_wp / 2.0_wp )**2 )
             ELSE
                betrag = 0.0_wp
             ENDIF

             pt(k,j,i) = pt(k,j,i) + betrag

          ENDDO
       ENDDO
    ENDDO

!
!-- Exchange of boundary values for temperature
    CALL exchange_horiz( pt, nbgp )


 END SUBROUTINE init_pt_anomaly
