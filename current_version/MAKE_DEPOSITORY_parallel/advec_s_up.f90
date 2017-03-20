!> @file advec_s_up.f90
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
! $Id: advec_s_up.f90 2001 2016-08-20 18:41:22Z knoop $
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
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY list
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
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! 981 2012-08-09 14:57:44Z maronga
! Typo removed
!
! Revision 1.1  1997/08/29 08:54:33  raasch
! Initial revision
!
!
! Description:
! ------------
!> Advection term for scalar quantities using the Upstream scheme.
!> NOTE: vertical advection at k=1 still has wrong grid spacing for w>0!
!>       The same problem occurs for all topography boundaries!
!------------------------------------------------------------------------------!
 MODULE advec_s_up_mod
 

    PRIVATE
    PUBLIC advec_s_up

    INTERFACE advec_s_up
       MODULE PROCEDURE advec_s_up
       MODULE PROCEDURE advec_s_up_ij
    END INTERFACE advec_s_up

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_up( sk )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzb_s_inner,&
                  nzt

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<

       REAL(wp) ::  ukomp !<
       REAL(wp) ::  vkomp !<
       REAL(wp) ::  wkomp !<
#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  sk
#endif


       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_s_inner(j,i)+1, nzt
!
!--             x-direction
                ukomp = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - u_gtrans
                IF ( ukomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - ukomp *                         &
                                         ( sk(k,j,i) - sk(k,j,i-1) ) * ddx
                ELSE
                   tend(k,j,i) = tend(k,j,i) - ukomp *                         &
                                         ( sk(k,j,i+1) - sk(k,j,i) ) * ddx
                ENDIF
!
!--             y-direction
                vkomp = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - v_gtrans
                IF ( vkomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - vkomp *                         &
                                         ( sk(k,j,i) - sk(k,j-1,i) ) * ddy
                ELSE
                   tend(k,j,i) = tend(k,j,i) - vkomp *                         &
                                         ( sk(k,j+1,i) - sk(k,j,i) ) * ddy
                ENDIF
!
!--             z-direction
                wkomp = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )
                IF ( wkomp > 0.0_wp )  THEN
                   tend(k,j,i) = tend(k,j,i) - wkomp *                         &
                                         ( sk(k,j,i) - sk(k-1,j,i) ) * ddzu(k)
                ELSE
                   tend(k,j,i) = tend(k,j,i) - wkomp *                         &
                                         ( sk(k+1,j,i)-sk(k,j,i) ) * ddzu(k+1)
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE advec_s_up


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_up_ij( i, j, sk )

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

       USE control_parameters,                                                 &
           ONLY:  u_gtrans, v_gtrans

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzb_s_inner, nzt

       USE kinds


       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<

       REAL(wp) ::  ukomp !<
       REAL(wp) ::  vkomp !<
       REAL(wp) ::  wkomp !<
       
#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  sk
#endif


       DO  k = nzb_s_inner(j,i)+1, nzt
!
!--       x-direction
          ukomp = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - u_gtrans
          IF ( ukomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - ukomp *                               &
                                         ( sk(k,j,i) - sk(k,j,i-1) ) * ddx
          ELSE
             tend(k,j,i) = tend(k,j,i) - ukomp *                               &
                                         ( sk(k,j,i+1) - sk(k,j,i) ) * ddx
          ENDIF
!
!--       y-direction
          vkomp = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - v_gtrans
          IF ( vkomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - vkomp *                               &
                                         ( sk(k,j,i) - sk(k,j-1,i) ) * ddy
          ELSE
             tend(k,j,i) = tend(k,j,i) - vkomp *                               &
                                         ( sk(k,j+1,i) - sk(k,j,i) ) * ddy
          ENDIF
!
!--       z-direction
          wkomp = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )
          IF ( wkomp > 0.0_wp )  THEN
             tend(k,j,i) = tend(k,j,i) - wkomp *                               &
                                         ( sk(k,j,i) - sk(k-1,j,i) ) * ddzu(k)
          ELSE
             tend(k,j,i) = tend(k,j,i) - wkomp *                               &
                                         ( sk(k+1,j,i)-sk(k,j,i) ) * ddzu(k+1)
          ENDIF

       ENDDO

    END SUBROUTINE advec_s_up_ij

 END MODULE advec_s_up_mod
