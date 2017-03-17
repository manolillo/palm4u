!> @file eqn_state_seawater.f90
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
! $Id: eqn_state_seawater.f90 2032 2016-10-21 15:13:51Z knoop $
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
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
! 97 2007-06-21 08:23:15Z raasch
! Initial revision
!
!
! Description:
! ------------
!> Equation of state for seawater as a function of potential temperature,
!> salinity, and pressure.
!> For coefficients see Jackett et al., 2006: J. Atm. Ocean Tech.
!> eqn_state_seawater calculates the potential density referred at hyp(0).
!> eqn_state_seawater_func calculates density.
!------------------------------------------------------------------------------!
 MODULE eqn_state_seawater_mod
 
    
    USE kinds

    IMPLICIT NONE

    PRIVATE
    PUBLIC eqn_state_seawater, eqn_state_seawater_func

    REAL(wp), DIMENSION(12), PARAMETER ::  nom =                               &
                          (/ 9.9984085444849347D2,   7.3471625860981584D0,     &
                            -5.3211231792841769D-2,  3.6492439109814549D-4,    &
                             2.5880571023991390D0,  -6.7168282786692354D-3,    &
                             1.9203202055760151D-3,  1.1798263740430364D-2,    &
                             9.8920219266399117D-8,  4.6996642771754730D-6,    &
                            -2.5862187075154352D-8, -3.2921414007960662D-12 /)
                          !<

    REAL(wp), DIMENSION(13), PARAMETER ::  den =                               &
                          (/ 1.0D0,                  7.2815210113327091D-3,    &
                            -4.4787265461983921D-5,  3.3851002965802430D-7,    &
                             1.3651202389758572D-10, 1.7632126669040377D-3,    &
                            -8.8066583251206474D-6, -1.8832689434804897D-10,   &
                             5.7463776745432097D-6,  1.4716275472242334D-9,    &
                             6.7103246285651894D-6, -2.4461698007024582D-17,   &
                            -9.1534417604289062D-18 /)
                          !<

    INTERFACE eqn_state_seawater
       MODULE PROCEDURE eqn_state_seawater
       MODULE PROCEDURE eqn_state_seawater_ij
    END INTERFACE eqn_state_seawater
 
    INTERFACE eqn_state_seawater_func
       MODULE PROCEDURE eqn_state_seawater_func
    END INTERFACE eqn_state_seawater_func
 
 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE eqn_state_seawater

       USE arrays_3d,                                                          &
           ONLY:  hyp, prho, pt_p, rho_ocean, sa_p
       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb_s_inner, nzt

       IMPLICIT NONE

       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
       INTEGER(iwp) ::  k  !<

       REAL(wp) ::  pden  !<
       REAL(wp) ::  pnom  !<
       REAL(wp) ::  p1    !<
       REAL(wp) ::  p2    !<
       REAL(wp) ::  p3    !<
       REAL(wp) ::  pt1   !<
       REAL(wp) ::  pt2   !<
       REAL(wp) ::  pt3   !<
       REAL(wp) ::  pt4   !<
       REAL(wp) ::  sa1   !<
       REAL(wp) ::  sa15  !<
       REAL(wp) ::  sa2   !<
       
                       

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_s_inner(j,i)+1, nzt
!
!--             Pressure is needed in dbar
                p1 = hyp(k) * 1E-4_wp
                p2 = p1 * p1
                p3 = p2 * p1

!
!--             Temperature needed in degree Celsius
                pt1 = pt_p(k,j,i) - 273.15_wp
                pt2 = pt1 * pt1
                pt3 = pt1 * pt2
                pt4 = pt2 * pt2

                sa1  = sa_p(k,j,i)
                sa15 = sa1 * SQRT( sa1 )
                sa2  = sa1 * sa1

                pnom = nom(1)           + nom(2)*pt1     + nom(3)*pt2     +    &
                       nom(4)*pt3       + nom(5)*sa1     + nom(6)*sa1*pt1 +    &
                       nom(7)*sa2

                pden = den(1)           + den(2)*pt1     + den(3)*pt2     +    &
                       den(4)*pt3       + den(5)*pt4     + den(6)*sa1     +    &
                       den(7)*sa1*pt1   + den(8)*sa1*pt3 + den(9)*sa15    +    &
                       den(10)*sa15*pt2

!
!--             Potential density (without pressure terms)
                prho(k,j,i) = pnom / pden

                pnom = pnom +             nom(8)*p1      + nom(9)*p1*pt2  +    &
                       nom(10)*p1*sa1   + nom(11)*p2     + nom(12)*p2*pt2

                pden = pden +             den(11)*p1     + den(12)*p2*pt3 +    &
                       den(13)*p3*pt1

!
!--             In-situ density
                rho_ocean(k,j,i) = pnom / pden

             ENDDO
!
!--          Neumann conditions are assumed at bottom and top boundary
             prho(nzt+1,j,i)            = prho(nzt,j,i)
             prho(nzb_s_inner(j,i),j,i) = prho(nzb_s_inner(j,i)+1,j,i)
             rho_ocean(nzt+1,j,i)             = rho_ocean(nzt,j,i)
             rho_ocean(nzb_s_inner(j,i),j,i)  = rho_ocean(nzb_s_inner(j,i)+1,j,i)

          ENDDO
       ENDDO

    END SUBROUTINE eqn_state_seawater


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE eqn_state_seawater_ij( i, j )

       USE arrays_3d,                                                          &
           ONLY:  hyp, prho, pt_p, rho_ocean, sa_p
           
       USE indices,                                                            &
           ONLY:  nzb_s_inner, nzt

       IMPLICIT NONE

       INTEGER(iwp) ::  i, j, k

       REAL(wp)     ::  pden, pnom, p1, p2, p3, pt1, pt2, pt3, pt4, sa1, sa15, &
                        sa2

       DO  k = nzb_s_inner(j,i)+1, nzt
!
!--       Pressure is needed in dbar
          p1 = hyp(k) * 1E-4_wp
          p2 = p1 * p1
          p3 = p2 * p1

!
!--       Temperature needed in degree Celsius
          pt1 = pt_p(k,j,i) - 273.15_wp
          pt2 = pt1 * pt1
          pt3 = pt1 * pt2
          pt4 = pt2 * pt2

          sa1  = sa_p(k,j,i)
          sa15 = sa1 * SQRT( sa1 )
          sa2  = sa1 * sa1

          pnom = nom(1)           + nom(2)*pt1     + nom(3)*pt2     +          &
                 nom(4)*pt3       + nom(5)*sa1     + nom(6)*sa1*pt1 +          &
                 nom(7)*sa2

          pden = den(1)           + den(2)*pt1     + den(3)*pt2     +          &
                 den(4)*pt3       + den(5)*pt4     + den(6)*sa1     +          &
                 den(7)*sa1*pt1   + den(8)*sa1*pt3 + den(9)*sa15    +          &
                 den(10)*sa15*pt2

!
!--       Potential density (without pressure terms)
          prho(k,j,i) = pnom / pden

          pnom = pnom +             nom(8)*p1      + nom(9)*p1*pt2  +          &
                 nom(10)*p1*sa1   + nom(11)*p2     + nom(12)*p2*pt2
          pden = pden +             den(11)*p1     + den(12)*p2*pt3 +          &
                 den(13)*p3*pt1

!
!--       In-situ density
          rho_ocean(k,j,i) = pnom / pden


       ENDDO

!
!--    Neumann conditions are assumed at bottom and top boundary
       prho(nzt+1,j,i)            = prho(nzt,j,i)
       prho(nzb_s_inner(j,i),j,i) = prho(nzb_s_inner(j,i)+1,j,i)
       rho_ocean(nzt+1,j,i)             = rho_ocean(nzt,j,i)
       rho_ocean(nzb_s_inner(j,i),j,i)  = rho_ocean(nzb_s_inner(j,i)+1,j,i)

    END SUBROUTINE eqn_state_seawater_ij


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Equation of state as a function
!------------------------------------------------------------------------------!
    REAL(wp) FUNCTION eqn_state_seawater_func( p, pt, sa )

       IMPLICIT NONE

       REAL(wp) ::  p      !<
       REAL(wp) ::  p1     !<
       REAL(wp) ::  p2     !<
       REAL(wp) ::  p3     !<
       REAL(wp) ::  pt     !<
       REAL(wp) ::  pt1    !<
       REAL(wp) ::  pt2    !<
       REAL(wp) ::  pt3    !<
       REAL(wp) ::  pt4    !<
       REAL(wp) ::  sa     !<
       REAL(wp) ::  sa15   !<
       REAL(wp) ::  sa2    !<

!
!--    Pressure is needed in dbar
       p1 = p  * 1E-4_wp
       p2 = p1 * p1
       p3 = p2 * p1

!
!--    Temperature needed in degree Celsius
       pt1 = pt - 273.15_wp
       pt2 = pt1 * pt1
       pt3 = pt1 * pt2
       pt4 = pt2 * pt2

       sa15 = sa * SQRT( sa )
       sa2  = sa * sa


       eqn_state_seawater_func =                                               &
         ( nom(1)        + nom(2)*pt1       + nom(3)*pt2    + nom(4)*pt3     + &
           nom(5)*sa     + nom(6)*sa*pt1    + nom(7)*sa2    + nom(8)*p1      + &
           nom(9)*p1*pt2 + nom(10)*p1*sa    + nom(11)*p2    + nom(12)*p2*pt2   &
         ) /                                                                   &
         ( den(1)        + den(2)*pt1       + den(3)*pt2    + den(4)*pt3     + &
           den(5)*pt4    + den(6)*sa        + den(7)*sa*pt1 + den(8)*sa*pt3  + &
           den(9)*sa15   + den(10)*sa15*pt2 + den(11)*p1    + den(12)*p2*pt3 + &
           den(13)*p3*pt1                                                      &
         )


    END FUNCTION eqn_state_seawater_func

 END MODULE eqn_state_seawater_mod
