!> @file advec_v_pw.f90
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
! $Id: advec_v_pw.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
! 
! 
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
! 
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
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
! Revision 1.1  1997/08/11 06:09:57  raasch
! Initial revision
!
!
! Description:
! ------------
!> Advection term for v velocity-component using Piacsek and Williams.
!> Vertical advection at the first grid point above the surface is done with
!> normal centred differences, because otherwise no information from the surface
!> would be communicated upwards due to w=0 at K=nzb.
!------------------------------------------------------------------------------!
 MODULE advec_v_pw_mod
 

    PRIVATE
    PUBLIC advec_v_pw

    INTERFACE advec_v_pw
       MODULE PROCEDURE advec_v_pw
       MODULE PROCEDURE advec_v_pw_ij
    END INTERFACE advec_v_pw
 
 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_v_pw

       USE arrays_3d,                                                          &
           ONLY:  ddzw, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nysv, nzb_v_inner, nzt

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<
       
       REAL(wp)    ::  gu !<
       REAL(wp)    ::  gv !<
 

       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans
       DO  i = nxl, nxr
          DO  j = nysv, nyn
             DO  k = nzb_v_inner(j,i)+1, nzt
                tend(k,j,i) = tend(k,j,i) - 0.25_wp * (                        &
                         ( v(k,j,i+1) * ( u(k,j-1,i+1) + u(k,j,i+1) - gu )     &
                         - v(k,j,i-1) * ( u(k,j-1,i) + u(k,j,i) - gu ) ) * ddx &
                       + ( v(k,j+1,i) * ( v(k,j+1,i) + v(k,j,i) - gv )         &
                         - v(k,j-1,i) * ( v(k,j,i) + v(k,j-1,i) - gv ) ) * ddy &
                       + ( v(k+1,j,i) * ( w(k,j-1,i) + w(k,j,i) )              &
                         - v(k-1,j,i) * ( w(k-1,j-1,i) + w(k-1,j,i) ) )        &
                                                                  * ddzw(k)    &
                                                      )
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE advec_v_pw


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_v_pw_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  ddzw, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nzb_v_inner, nzt

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<
       
       REAL(wp)    ::  gu !<
       REAL(wp)    ::  gv !<


       gu = 2.0_wp * u_gtrans
       gv = 2.0_wp * v_gtrans
       DO  k = nzb_v_inner(j,i)+1, nzt
          tend(k,j,i) = tend(k,j,i) - 0.25_wp * (                              &
                         ( v(k,j,i+1) * ( u(k,j-1,i+1) + u(k,j,i+1) - gu )     &
                         - v(k,j,i-1) * ( u(k,j-1,i) + u(k,j,i) - gu ) ) * ddx &
                       + ( v(k,j+1,i) * ( v(k,j+1,i) + v(k,j,i) - gv )         &
                         - v(k,j-1,i) * ( v(k,j,i) + v(k,j-1,i) - gv ) ) * ddy &
                       + ( v(k+1,j,i) * ( w(k,j-1,i) + w(k,j,i) )              &
                         - v(k-1,j,i) * ( w(k-1,j-1,i) + w(k-1,j,i) ) )        &
                                                                  * ddzw(k)    &
                                                )
       ENDDO
 
    END SUBROUTINE advec_v_pw_ij

 END MODULE advec_v_pw_mod
 
