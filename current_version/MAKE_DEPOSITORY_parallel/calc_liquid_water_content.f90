!> @file calc_liquid_water_content.f90
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
! $Id: calc_liquid_water_content.f90 2001 2016-08-20 18:41:22Z knoop $
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme removed. microphysics_seifert added.
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
! 1253 2013-11-07 10:48:12Z fricke
! Bugfix: q is set to qr in case that q is smaller than qr
!
! 1115 2013-03-26 18:16:16Z hoffmann
! drizzle can be used independently from precipitation
!
! 1053 2012-11-13 17:11:03Z hoffmann
! description expanded to the two-moment cloud scheme
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  2000/04/13 14:50:45  schroeter
! Initial revision
!
!
!
! Description:
! ------------
!> Calculation of the liquid water content (0%-or-100%-scheme). This scheme is 
!> used by the one and the two moment cloud physics scheme. Using the two moment 
!> scheme, this calculation results in the cloud water content.
!------------------------------------------------------------------------------!
 SUBROUTINE calc_liquid_water_content
 


    USE arrays_3d,                                                             &
        ONLY:  hyp, pt, q, qc, ql, qr

    USE cloud_parameters,                                                      &
        ONLY:  l_d_cp, l_d_r, t_d_pt

    USE control_parameters,                                                    &
        ONLY:  microphysics_seifert

    USE indices,                                                               &
        ONLY:  nxlg, nxrg, nyng, nysg, nzb_s_inner, nzt

    USE kinds

    USE pegrid


    IMPLICIT NONE

    INTEGER(iwp) ::  i !<
    INTEGER(iwp) ::  j !<
    INTEGER(iwp) ::  k !<

    REAL(wp) ::  alpha !<
    REAL(wp) ::  e_s   !<
    REAL(wp) ::  q_s   !<
    REAL(wp) ::  t_l   !<

    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          DO  k = nzb_s_inner(j,i)+1, nzt

!
!--          Compute the liquid water temperature
             t_l = t_d_pt(k) * pt(k,j,i)

!
!--          Compute saturation water vapor pressure at t_l
             e_s = 610.78_wp * EXP( 17.269_wp * ( t_l - 273.16_wp ) /          &
                                                ( t_l - 35.86_wp ) )

!
!--          Compute approximation of saturation humidity
             q_s = 0.622_wp * e_s / ( hyp(k) - 0.378_wp * e_s )

!
!--          Correction factor
             alpha = 0.622_wp * l_d_r * l_d_cp / ( t_l * t_l )

!
!--          Correction of the approximated value
!--          (see: Cuijpers + Duynkerke, 1993, JAS, 23)
             q_s = q_s * ( 1.0_wp + alpha * q(k,j,i) ) / ( 1.0_wp + alpha * q_s )

!
!--          Compute the liquid water content
             IF ( microphysics_seifert )  THEN
                IF ( ( q(k,j,i) - q_s - qr(k,j,i) ) > 0.0_wp ) THEN
                   qc(k,j,i) = q(k,j,i) - q_s - qr(k,j,i)
                   ql(k,j,i) = qc(k,j,i) + qr(k,j,i)
                ELSE
                   IF ( q(k,j,i) < qr(k,j,i) )  q(k,j,i) = qr(k,j,i)
                   qc(k,j,i) = 0.0_wp 
                   ql(k,j,i) = qr(k,j,i)
                ENDIF
             ELSE
                IF ( ( q(k,j,i) - q_s ) > 0.0_wp ) THEN
                   qc(k,j,i) = q(k,j,i) - q_s
                   ql(k,j,i) = qc(k,j,i)
                ELSE
                   qc(k,j,i) = 0.0_wp
                   ql(k,j,i) = 0.0_wp
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
 END SUBROUTINE calc_liquid_water_content
