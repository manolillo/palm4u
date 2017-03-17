 PROGRAM interpret_config

!--------------------------------------------------------------------------------!
! This file is part of PALM.
!
! PALM is free software: you can redistribute it and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation,
! either version 3 of the License, or (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2014  Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
!
!
! Former revisions:
! -----------------
! $Id: interpret_config.f90 2018 2016-09-29 06:23:32Z raasch $
!
! 2017 2016-09-29 06:22:33Z raasch
! Bugfix in case of # of output files > 100
!
! 1289 2014-03-04 07:12:34Z raasch
! routine local_getenv removed
!
! 1094 2013-02-03 01:52:12Z raasch
! unused variables removed
!
! 1046 2012-11-09 14:38:45Z maronga
! code put under GPL (PALM 3.9)
!
! 785 2011-11-28 09:47:19Z raasch
! character length extended to 1000 in order to account for long lines in the
! configuration file
!
! 07/09/10 - Siggi - bugfix in 2201 statement: closing " was missing
! unknown  - Siggi - mrun environment variables are read from NAMELIST instead
!                    of using GETENV. Variables are always assigned a value,
!                    also if they already got one. These values are re-assigned
!                    later in mrun.
! 28/02/07 - Siggi - empty lines in configuration file are accepted
! 01/11/05 - Siggi - s2b-Feld erlaubt den Wert locopt
! 29/06/05 - Siggi - Fehlermeldung ins englische uebertragen und ergaenzt
! 29/04/05 - Siggi - extin wird auch fuer Input-Dateien ausgegeben
! 18/11/97 - Siggi - Komma in 2010-FORMAT hinzugefuegt
! 21/07/97 - Siggi - Erste Fassung
!
! Description:
! -------------
! This program reads the mrun-configuration file .mrun.config and outputs
! its content in form of ksh-commands, which can then be executed by mrun.
! mrun is also able to directly read from the configuration file by using
! the option "-S" (but with much slower speed).
!------------------------------------------------------------------------------!

    IMPLICIT NONE

    CHARACTER (LEN=1)    ::  bs = ACHAR( 92 )   ! backslash (auf vpp sonst n.
                                               ! druckbar)
    CHARACTER (LEN=20)   ::  do_remote, do_trace, host, localhost
    CHARACTER (LEN=100)  ::  config_file, icf
    CHARACTER (LEN=1000) ::  cond1, cond2, empty = REPEAT( ' ', 240 ),       &
                             for_cond1, for_cond2, for_host, input_list,     &
                             iolist, output_list, s1, s2, s2a, s2b, s2c, s3, &
                             s3cond, s4, s5, s6, value, var, zeile

    INTEGER ::  i, icomment = 0, icond1, icond2, idatver = 0, iec = 0,         &
                ienvvar = 0, ifor_cond1, ifor_cond2, ifor_host, ihost,         &
                iic = 0, iicf, iin = 0, iinput_list, il, ilocalhost,  ioc = 0, &
                ios, iout = 0, ioutput_list, is1, is2, is2a, is2b, is2c,       &
                is3, is3cond, is4, is5, is6, ivalue, ivar, izeile, linenr

    LOGICAL ::  found

    NAMELIST /mrun_environment/  cond1, cond2, config_file, do_remote,       &
                                 do_trace, host, input_list, icf, localhost, &
                                 output_list


    OPEN ( 1, FILE='.mrun_environment', FORM='FORMATTED' )
    READ ( 1, mrun_environment )
    CLOSE ( 1 )

    icond1       = LEN_TRIM( cond1 )
    icond2       = LEN_TRIM( cond2 )
    il           = LEN_TRIM( config_file )
    ihost        = LEN_TRIM( host )
    iinput_list  = LEN_TRIM( input_list )
    iicf         = LEN_TRIM( icf )
    ilocalhost   = LEN_TRIM( localhost )
    ioutput_list = LEN_TRIM( output_list )

    iolist = input_list(1:iinput_list) // output_list(1:ioutput_list)

    IF ( do_trace(1:4) == 'true' )  THEN
       PRINT*,'*** cond1="',cond1(1:icond1),'"'
       PRINT*,'*** cond2="',cond2(1:icond2),'"'
       PRINT*,'*** config_file="',config_file(1:il),'"'
       PRINT*,'*** do_remote="',do_remote,'"'
       PRINT*,'*** do_trace="',do_trace,'"'
       PRINT*,'*** host="',host(1:ihost),'"'
       PRINT*,'*** input_list="',input_list(1:iinput_list),'"'
       PRINT*,'*** interpreted_config_file="',icf(1:iicf),'"'
       PRINT*,'*** localhost="',localhost(1:ilocalhost),'"'
       PRINT*,'*** output_list="',output_list(1:ioutput_list),'"'
    ENDIF

    OPEN ( 1, FILE=config_file(1:il), FORM='formatted' )
    OPEN ( 2, FILE=icf(1:iicf), FORM='formatted' )

    READ ( 1, '(A)', IOSTAT=ios )  zeile
    linenr = 1


    DO WHILE ( ios == 0 )

       izeile = LEN_TRIM( zeile )

       IF ( LEN_TRIM( zeile ) == 0 )  THEN

          CONTINUE

       ELSEIF ( zeile(1:1) == '#' )  THEN

          icomment = icomment + 1

       ELSEIF ( zeile(1:1) == '%' )  THEN

          ienvvar = ienvvar + 1
          i = INDEX( zeile, ' ' )
          var = zeile(2:i-1)
          ivar = i - 2

          zeile(1:i) = empty(1:i)
          zeile = ADJUSTL( zeile )
          i = INDEX( zeile, ' ' )
          value = zeile(1:i-1)
          ivalue = i - 1

          zeile(1:i) = empty(1:i)
          zeile = ADJUSTL( zeile )
          i = INDEX( zeile, ' ' )

          IF ( i /= 1 )  THEN
             for_host = zeile(1:i-1)
             ifor_host = i - 1

             zeile(1:i) = empty(1:i)
             zeile = ADJUSTL( zeile )
             i = INDEX( zeile, ' ' )

             IF ( i /= 1 )  THEN
                for_cond1 = zeile(1:i-1)
                ifor_cond1 = i - 1

                zeile(1:i) = empty(1:i)
                zeile = ADJUSTL( zeile )
                i = INDEX( zeile, ' ' )

                IF ( i /= 1 )  THEN
                   for_cond2 = zeile(1:i-1)
                   ifor_cond2 = i - 1
                ELSE
                   for_cond2 = ''
                   ifor_cond2 = 0
                ENDIF
             ELSE
                for_cond1 = ''
                ifor_cond1 = 0
                for_cond2 = ''
                ifor_cond2 = 0
             ENDIF
          ELSE
             for_host = ' '
             ifor_host = 1
             for_cond1 = ''
             ifor_cond1 = 0
             for_cond2 = ''
             ifor_cond2 = 0
          ENDIF
          IF ( do_trace(1:4) == 'true' )  THEN
             PRINT*,'var="',var(1:ivar),'"'
             PRINT*,'value="',value(1:ivalue),'"'
             PRINT*,'for_host="',for_host(1:ifor_host),'"'
             PRINT*,'for_cond1="',for_cond1(1:ifor_cond1),'"'
             PRINT*,'for_cond2="',for_cond2(1:ifor_cond2),'"'
          ENDIF
!
!--       Geltungsbereich pruefen und evtl. Variable ausgeben
          IF ( for_host == ' '  .OR.  ( &
               for_host(1:ifor_host) == host(1:ihost)  .AND. &
               for_cond1(1:ifor_cond1) == cond1(1:icond1)  .AND. &
               for_cond2(1:ifor_cond2) == cond2(1:icond2) &
                                      )  .OR. ( &
               INDEX( iolist, for_host(1:ifor_host) ) /= 0 &
                                              ) )  THEN

!
!--          Zuerst Doppelpunkte durch Blanks ersetzen (aber doppelt
!--          auftretende Doppelpunkte durch einen Doppelpunkt)
             i = 0
             DO
                i = i + 1
                IF ( i > ivalue )  EXIT
                IF ( value(i:i) == ':' )  THEN
                   IF ( value(i+1:i+1) == ':' )  THEN
                      value = value(1:i) // value(i+2:ivalue)
                      ivalue = ivalue - 1
                   ELSE
                      value(i:i) = ' '
                   ENDIF
                ENDIF
             ENDDO

!
!--          Variable ausgeben
             WRITE (2,2200)  var(1:ivar), bs, value(1:ivalue), bs, &
                             var(1:ivar)
 2200        FORMAT ('eval ',A,'=',A,'"',A,A,'"'/'export ',A)

             IF ( do_trace(1:4) == 'true' )  THEN
                WRITE (2,2201)  bs, var(1:ivar), value(1:ivalue)
 2201           FORMAT ('printf "',A,'n*** ENVIRONMENT-VARIABLE ',A,' = ',A,'"')
             ENDIF

          ENDIF

!
!--       Variable "host" muss gleich ausgewertet werden, da mit ihr ein
!--       neuer Geltungsbereich festgelegt wird
          IF ( var(1:ivar) == 'host' )  THEN

             host  = value(1:ivalue)
             ihost = ivalue

          ENDIF

       ELSEIF ( zeile(1:3) == 'EC:' )  THEN
!
!--       Error-Kommandos
          iec = iec + 1
          IF ( iec < 10 )  THEN
             WRITE (2,'(''err_command['',I1,'']="'',A,''"'')')  iec, &
                                                                zeile(4:izeile)
          ELSEIF ( iec < 100 )  THEN
             WRITE (2,'(''err_command['',I2,'']="'',A,''"'')')  iec, &
                                                                zeile(4:izeile)
          ELSE
             WRITE (2,'(''err_command['',I3,'']="'',A,''"'')')  iec, &
                                                                zeile(4:izeile)
          ENDIF

       ELSEIF ( zeile(1:3) == 'IC:' )  THEN
!
!--       Input-Kommandos
          iic = iic + 1
          IF ( iic < 10 )  THEN
             WRITE (2,'(''in_command['',I1,'']="'',A,''"'')')  iic, &
                                                               zeile(4:izeile)
          ELSEIF ( iic < 100 )  THEN
             WRITE (2,'(''in_command['',I2,'']="'',A,''"'')')  iic, &
                                                               zeile(4:izeile)
          ELSE
             WRITE (2,'(''in_command['',I3,'']="'',A,''"'')')  iic, &
                                                               zeile(4:izeile)
          ENDIF

       ELSEIF ( zeile(1:3) == 'OC:' )  THEN
!
!--       Output-Kommandos
          ioc = ioc + 1
          IF ( ioc < 10 )  THEN
             WRITE (2,'(''out_command['',I1,'']="'',A,''"'')')  ioc, &
                                                                zeile(4:izeile)
          ELSEIF ( ioc < 100 )  THEN
             WRITE (2,'(''out_command['',I2,'']="'',A,''"'')')  ioc, &
                                                                zeile(4:izeile)
          ELSE
             WRITE (2,'(''out_command['',I3,'']="'',A,''"'')')  ioc, &
                                                                zeile(4:izeile)
          ENDIF

       ELSE
!
!--       Dateiverbindungen
          idatver = idatver + 1
!
!--       Lokaler Name
          i   = INDEX( zeile , ' ' )
          s1  = zeile(1:i-1)
          is1 = i-1
!
!--       Dateieigenschaften
          zeile = ADJUSTL( zeile(i:izeile) )
          i   = INDEX( zeile , ' ' )
          s2  = zeile(1:i-1)
          is2 = i-1
!
!--       Geltungsbereich
          zeile = ADJUSTL( zeile(i:izeile) )
          i   = INDEX( zeile , ' ' )
          s3  = zeile(1:i-1)
          is3 = i-1
!
!--       Pfadname
          zeile = ADJUSTL( zeile(i:izeile) )
          i   = INDEX( zeile , ' ' )
          s4  = zeile(1:i-1)
          is4 = i-1
!
!--       evtl. Extension
          zeile = ADJUSTL( zeile(i:izeile) )
          i = INDEX( zeile , ' ' )
          IF ( i == 1 )  THEN
             s5  = ' '
             is5 = 1
             s6  = ' '
             is6 = 1
          ELSE
             s5  = zeile(1:i-1)
             is5 = i-1
!
!--          evtl. 2. Extension
             zeile = ADJUSTL( zeile(i:izeile) )
             i = INDEX( zeile , ' ' )
             IF ( i == 1 )  THEN
                s6  = ' '
                is6 = 1
             ELSE
                s6  = zeile(1:i-1)
                is6 = i-1
             ENDIF
          ENDIF

!
!--       Dateieigenschaften aufspalten
          i = INDEX( s2 , ':' )
          IF ( i == 0 )  THEN
             s2a  = s2
             is2a = is2
             s2b  = ''
             is2b = 0
             s2c  = ''
             is2c = 0
          ELSE
             s2a  = s2(1:i-1)
             is2a = i-1
             s2   = s2(i+1:is2)

             i = INDEX( s2 , ':' )
             IF ( i == 0 )  THEN
                s2b  = s2
                is2b = LEN_TRIM( s2 )
                s2c  = ''
                is2c = 0
             ELSE
                s2b  = s2(1:i-1)
                is2b = i-1
                s2c  = s2(i+1:)
                is2c = LEN_TRIM( s2c )
             ENDIF
          ENDIF
!
!--       Pruefung, ob Eingabedateiverbindung abgespeichert werden soll
          IF ( s2a(1:is2a) == 'in'  .AND.  .NOT. (                     &
               do_remote(1:4) == 'true'  .AND.                         &
               ( s2b(1:is2b) == 'loc'  .OR.  s2b(1:is2b) == 'locopt' ) &
                                                 ) )  THEN
             found = .FALSE.
             i = INDEX( s3 , ':' )
             IF ( i == 0 )  THEN
                s3cond  = s3
                is3cond = LEN_TRIM( s3cond )
             ELSE
                s3cond  = s3(1:i-1)
                is3cond = i-1
                s3      = s3(i+1:)
             ENDIF

             DO WHILE ( s3cond(1:1) /= ' ' )

                IF ( INDEX( input_list(1:iinput_list) , s3cond(1:is3cond) ) /= 0 &
                     .OR.  s3cond(1:is3cond) == '-' )  THEN
                   found = .TRUE.
                ENDIF

                IF ( s3(1:1) == ' ' )  THEN
                   s3cond = ' '
                ELSE
                   i = INDEX( s3 , ':' )
                   IF ( i == 0 )  THEN
                      s3cond  = s3
                      is3cond = LEN_TRIM( s3cond )
                      s3      = ' '
                   ELSE
                      s3cond  = s3(1:i-1)
                      is3cond = i-1
                      s3      = s3(i+1:)
                   ENDIF
                ENDIF

             ENDDO

!
!--          Wenn Geltungsbereich erfuellt, dann Dateiverbindung abspeichern
             IF ( found )  THEN

                iin = iin + 1
                IF ( iin < 10 )  THEN
                   WRITE (2,2000)  iin, s1(1:is1), iin, s2b(1:is2b), &
                                   iin, s2c(1:is2c), &
                                   iin, s3(1:is3), iin, s4(1:is4), &
                                   iin, s5(1:is5), iin, s6(1:is6)
2000               FORMAT ('localin[',I1,']="',A,'"; transin[',I1,']="',A, &
                           '"; actionin[',I1,']="',A, &
                           '"; typein[',I1,']="',A,'"'/'pathin[',I1,']="',A, &
                           '"; endin[',I1,']="',A,'"; extin[',I1,']="',A,'"')
                ELSEIF ( iin < 100 )  THEN
                   WRITE (2,2001)  iin, s1(1:is1), iin, s2b(1:is2b), &
                                   iin, s2c(1:is2c), &
                                   iin, s3(1:is3), iin, s4(1:is4), &
                                   iin, s5(1:is5), iin, s6(1:is6)
2001               FORMAT ('localin[',I2,']="',A,'"; transin[',I2,']="',A, &
                           '"; actionin[',I2,']="',A, &
                           '"; typein[',I2,']="',A,'"'/'pathin[',I2,']="',A, &
                           '"; endin[',I2,']="',A,'"; extin[',I2,']="',A,'"')
                ELSE
                   WRITE (2,2002)  iin, s1(1:is1), iin, s2b(1:is2b), &
                                   iin, s2c(1:is2c), &
                                   iin, s3(1:is3), iin, s4(1:is4), &
                                   iin, s5(1:is5), iin, s6(1:is6)
2002               FORMAT ('localin[',I3,']="',A,'"; transin[',I3,']="',A, &
                           '"; actionin[',I3,']="',A, &
                           '"; typein[',I3,']="',A,'"'/'pathin[',I3,']="',A, &
                           '"; endin[',I3,']="',A,'"; extin[',I3,']="',A,'"')
                ENDIF
             ENDIF

          ELSEIF ( s2a(1:is2a) == 'out'  .AND.  .NOT. ( &
                   do_remote(1:4) == 'true'  .AND.  s2b(1:is2b) == 'loc' &
                                                      ) )  THEN
!
!--          Pruefung, ob Ausgabedateiverbindung abgespeichert werden soll
             found = .FALSE.
             i = INDEX( s3 , ':' )
             IF ( i == 0 )  THEN
                s3cond  = s3
                is3cond = LEN_TRIM( s3cond )
             ELSE
                s3cond  = s3(1:i-1)
                is3cond = i-1
                s3      = s3(i+1:)
             ENDIF

             DO WHILE ( s3cond(1:1) /= ' ' )

                IF ( INDEX( output_list(1:ioutput_list) , s3cond(1:is3cond) ) /= 0 &
                     .OR.  s3cond(1:is3cond) == '-' )  THEN
                   found = .TRUE.
                ENDIF

                IF ( s3(1:1) == ' ' )  THEN
                   s3cond = ' '
                ELSE
                   i = INDEX( s3 , ':' )
                   IF ( i == 0 )  THEN
                      s3cond  = s3
                      is3cond = LEN_TRIM( s3cond )
                      s3      = ' '
                   ELSE
                      s3cond  = s3(1:i-1)
                      is3cond = i-1
                      s3      = s3(i+1:)
                   ENDIF
                ENDIF

             ENDDO
!
!--          Wenn Geltungsbereich erfuellt, dann Dateiverbindung abspeichern
             IF ( found )  THEN

                iout = iout + 1
                IF ( iout < 10 )  THEN
                   WRITE (2,2003)  iout, s1(1:is1), iout, s2c(1:is2c), &
                                   iout, s3(1:is3), iout, s4(1:is4), &
                                   iout, s5(1:is5), iout, s6(1:is6)
 2003              FORMAT ('localout[',I1,']="',A,'"; actionout[',I1,']="',A, &
                           '"; typeout[',I1,']="',A,'"'/'pathout[',I1,']="',A, &
                           '"; endout[',I1,']="',A,'"; extout[',I1,']="',A,'"')
                ELSEIF ( iout < 100 )  THEN
                      WRITE (2,2004)  iout, s1(1:is1), iout, s2c(1:is2c), &
                                      iout, s3(1:is3), iout, s4(1:is4), &
                                      iout, s5(1:is5), iout, s6(1:is6)
 2004              FORMAT ('localout[',I2,']="',A,'"; actionout[',I2,']="',A, &
                           '"; typeout[',I2,']="',A,'"'/'pathout[',I2,']="',A, &
                           '"; endout[',I2,']="',A,'"; extout[',I2,']="',A,'"')
                ELSE
                      WRITE (2,2005)  iout, s1(1:is1), iout, s2c(1:is2c), &
                                      iout, s3(1:is3), iout, s4(1:is4), &
                                      iout, s5(1:is5), iout, s6(1:is6)
 2005              FORMAT ('localout[',I3,']="',A,'"; actionout[',I3,']="',A, &
                           '"; typeout[',I3,']="',A,'"'/'pathout[',I3,']="',A, &
                           '"; endout[',I3,']="',A,'"; extout[',I3,']="',A,'"')
                ENDIF
             ENDIF

          ELSEIF ( s2a(1:is2a) /= 'in'  .AND.  s2a(1:is2a) /= 'out' )  THEN
!
!--          Kein gueltiger Wert fuer I/O-Feld
             WRITE (2,2010)  bs, bs, config_file(1:il), linenr, bs, bs, &
                             s2a(1:is2a), bs, bs, bs, bs, bs
 2010        FORMAT ('printf "',A,'n',A,'n +++ I/O-field in configuration ', &
                     'file ',A, ', line ', I5, ' has the illegal"'/          &
                     'printf "',A,'n     value ',A,'"',A,A,'". Only ',       &
                     A,'"in',A,'" or ',A,'"out',A,'" are allowed!"'          &
                    )
             WRITE (2,'(''locat=connect; exit'')')
             STOP
          ENDIF
          
       ENDIF

       READ( 1, '(A)', IOSTAT=ios )  zeile
       linenr = linenr + 1

    ENDDO

!
!-- Ausgabe der Anzahl von gefundenen Zeilen
    IF ( iec > 0 )  WRITE (2,'(''(( iec = '',I3,'' ))'')')  iec
    IF ( iic > 0 )  WRITE (2,'(''(( iic = '',I3,'' ))'')')  iic
    IF ( ioc > 0 )  WRITE (2,'(''(( ioc = '',I3,'' ))'')')  ioc
    IF ( iin > 0 )  WRITE (2,'(''(( iin = '',I3,'' ))'')')  iin
    IF ( iout > 0 )  WRITE (2,'(''(( iout = '',I3,'' ))'')')  iout

    IF ( do_trace(1:4) == 'true' )  THEN
       PRINT*,' '
       PRINT*,'*** Inhalt von: ',config_file(1:il)
       PRINT*,icomment,' Kommentarzeilen'
       PRINT*,ienvvar,' Environment-Variablen-Vereinbarungen'
       PRINT*,iec,' Error-Kommandos'
       PRINT*,iic,' Input-Kommandos'
       PRINT*,ioc,' Output-Kommandos'
       PRINT*,idatver,' Dateiverbindungs-Anweisungen'
       PRINT*,'Davon interpretiert:'
       PRINT*,iin,' Eingabedateien'
       PRINT*,iout,' Ausgabedateien'
    ENDIF

 END PROGRAM interpret_config
