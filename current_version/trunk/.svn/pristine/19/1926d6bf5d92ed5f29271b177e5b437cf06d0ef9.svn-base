#!/bin/bash

# subjob - script for automatic generation and submission of batch-job files
#          for various batch queuing systems

#--------------------------------------------------------------------------------#
# This file is part of PALM.
#
# PALM is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# PALM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM. If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 1997-2014  Leibniz Universitaet Hannover
#--------------------------------------------------------------------------------#
#
# Current revisions:
# ------------------
# 
# 
# Former revisions:
# -----------------
# $Id$
#
# 1944 2016-06-15 06:29:00Z raasch
# adjustments for using HLRN ssh-keys
#
# 1940 2016-06-14 05:15:20Z raasch
# adjustments for lckiaps
#
# 1866 2016-04-15 06:50:59Z raasch
# adjusted for lcocean
#
# 1841 2016-04-07 19:14:06Z raasch
# script now running under bash
# 
# 1701 2015-11-02 07:43:04Z maronga
# Bugfix: added missing init_cmds for lccrayh/lccrayb
# 
# 1621 2015-07-17 11:39:33Z heinze
# adjustments for Mistral at DKRZ Hamburg (lcbullhh)
#
# 1575 2015-03-27 09:56:27Z raasch
# mpp2-queues added to lccrayh
#
# 1547 2015-01-29 15:09:12Z witha
# adjustments for ForWind computing cluster (lcflow)
# 
# 1545 2015-01-29 06:52:23Z heinze
# local host name for blizzard further specified
#
# 1480 2014-10-17 14:41:49Z raasch
# adjustments for 2nd stage of HLRNIII
#
# 1468 2014-09-24 14:06:57Z maronga
# Typo removed (addres->address)
# Adjustments for lcxe6
# 
# 1452 2014-08-22 09:41:06Z heinze
# local hosts for blizzard added
# 
# 1450 2014-08-21 07:31:51Z heinze
# HLRN-III (lccrayb): testq queue adjusted to mpp1testq
# 
# 1442 2014-07-28 07:09:10Z raasch
# HLRN-III (lccrayb/lccrayh) queues adjusted
#
# 1378 2014-04-28 06:04:58Z raasch
# -et option added for lctit
#
# 1350 2014-04-04 13:01:30Z maronga
# location of qsub updated for lcxe6
#
# 1289 2014-03-04 07:12:34Z raasch
# German comments translated to English
# fimm-, necriam-, scirocco-, ibmy-, and sgi-specific code removed
#
# 1279 2014-01-28 12:10:14Z raasch
# node calculation modified due to changes in mrun (tasks_per_node must not be
# an integral divisor of numprocs any more)
#
# 1274 2014-01-09 13:14:54Z heinze
# adjustments for lccrayh
#
# 1266 2013-12-11 12:07:34Z heinze
# further adjustments for lccrayb (use msub instead of qsub)
#
# 1264 2013-12-09 12:46:09Z fricke
# Bugfix: Using number of nodes instead of number of processors (lccrayb)
#
# 1262 2013-12-09 10:57:20Z fricke
# further adjustments for lccrayb
#
# 1260 2013-12-04 12:48:04Z raasch
# jaboticaba admitted
#
# 1255 2013-11-07 14:43:35Z raasch
# further adjustments for lccrayb
#
# 1224 2013-09-16 07:27:23Z raasch
# first adjustments for lccrayb
#
# 1202 2013-07-10 16:22:07Z witha
# adjustments for Forwind cluster (lcflow)
#
# 1199 2013-07-05 14:52:22Z raasch
# adjustments for CSC Helsinki (lccrayf)
#
# use of cluster/express queue enabled (ibmh)
# vinessa added (imuk)
#
# 1103 2013-02-20 02:15:53Z raasch
# bash compatibility adjustments (usage of OPTIND, output formatting with printf
# instead typeset -L/R),
# further adjustments for lckyuh
#
# 2013-02-10 01:47:43Z raasch
# adjustments for Kyushu-Univeristy computing center (lckyuh - hayaka)
# and for Forwind cluster (lcflow)
#
# 1094 2013-02-03 01:52:12Z raasch
# new option -P for explicit setting of ssh/scp port,
# decalpha parts (yonsei) removed
#
# 2013-02-02 07:06:13Z raasch
# adjustments for Kyushu-University computing center (lckyut - tatara)
# old changelog messages removed
#
# 1046 2012-11-09 14:38:45Z maronga
# code put under GPL (PALM 3.9)
#
# 08/07/94 - Siggi - first version finished
# 29/06/94 - Siggi - script development started
#--------------------------------------------------------------------------------#
# subjob - script for automatic generation and submission of batch-job files
#          for various batch queuing systems
#--------------------------------------------------------------------------------#


    # VARIABLE-DECLARATIONS AND DEFAULT VALUES
 delete_dayfile=false
 email_notification=none
 group_number=none
 locat=normal
 no_default_queue=none
 no_submit=false
 job_catalog="~/job_queue"
 job_name=none
 local_user=$LOGNAME
 node_usage=shared
 numprocs=0
 punkte="..........................................................."
 submcom=qsub
 queue=default
 remote_host=none
 remote_user=""
 verify=true

 typeset  -i   cputime=0  memory=0  Memory=0  minuten  resttime  sekunden  stunden
 typeset  -i   numprocs  mpi_tasks=0  nodes=0  processes_per_node=0 tasks_per_node=0  threads_per_task=1



    # ERROR HANDLING
    # IN CASE OF EXIT:
 trap 'if [[ $locat != normal ]]
       then
          case  $locat  in
             (option)  printf "\n  --> available optios can be displayed"
                       printf " by typing:"
                       printf "\n      \"subjob ?\" \n";;
             (ftpcopy|parameter|scp|verify)  printf "\n";;
             (*)       printf "\n  +++ unknown error"
                       printf "\n      please inform S. Raasch!\n"
          esac
          [[ -f $job_to_send ]]  &&  rm  $job_to_send
          printf "\n\n+++ SUBJOB killed \n\n"
       fi' exit


    # IN CASE OF TERMINAL-BREAK:
 trap '[[ -f $job_to_send ]]  &&  rm  $job_to_send
       printf "\n\n+++ SUBJOB killed \n\n"
       exit
      ' 2


    # DETERMINE NAME OF LOCAL HOST
 local_host=$(hostname)

    # SET HOST-SPECIFIC VARIABLES VEREINBAREN (CHECK, IF LOCAL HOST
    # IS ADMITTED AT ALL)
    # NOTE: ONE OF THE ENTRIES FOR "lck" OR "lckordi" ALWAYS HAS TO BE
    # COMMENT OUT, BECAUSE THE HOSTNAME (node*) IS SAME FOR BOTH MACHINES
 case  $local_host  in
     (ambiel-lx)             local_address=134.106.74.48;  local_host=lcfor;;
     (atmos)                 local_address=172.20.25.35;   local_host=lcide;;
     (austru)                local_address=130.75.105.128; local_host=lcmuk;;
     (autan)                 local_address=130.75.105.57;  local_host=lcmuk;;
     (bora)                  local_address=130.75.105.103; local_host=lcmuk;;
     (b04*)                  local_address=133.5.4.33;     local_host=lckyuh;;
     (blizzard1|p0*|p1*|p2*|p3*|p4*|p5*|p6*|p7*|p8*|p9*)   local_address=136.172.40.15;  local_host=ibmh;;
     (blizzard2|p0*|p1*|p2*|p3*|p4*|p5*|p6*|p7*|p8*|p9*)   local_address=136.172.40.16;  local_host=ibmh;;
     (blogin*|bxc*)          local_address=130.73.233.1;   local_host=lccrayb;;
     (hlogin*|hxc*)          local_address=130.75.4.1;     local_host=lccrayh;;
     (breva)                 local_address=130.75.105.98;  local_host=lcmuk;;
     (buran)                 local_address=130.75.105.58;  local_host=lcmuk;;
     (caurus)                local_address=130.75.105.19;  local_host=lcmuk;;
     (climate*)              local_address=165.132.26.68;  local_host=lcyon;;
     (clogin*)               local_address=86.50.166.21;   local_host=lccrayf;;
     (cs*)                   local_address=136.172.44.131; local_host=nech;;
     (elephanta)             local_address=130.75.105.6;   local_host=lcmuk;;
     (flow01)                local_address=10.141.255.71;  local_host=lcflow;;
     (flow02)                local_address=10.141.255.72;  local_host=lcflow;;
     (node*)                 local_address=165.132.26.61   local_host=lck;;
   #  (node*)                 local_address=210.219.61.8    local_host=lckordi;;
     (gaia*)                 local_address=150.183.146.24; local_host=ibmkisti;;
     (gallego)               local_address=130.75.105.10;  local_host=lcmuk;;
     (gregale)               local_address=130.75.105.109; local_host=lcmuk;;
     (hababai)               local_address=130.75.105.108; local_host=lcmuk;;
     (hayaka*)               local_address=133.5.4.33;     local_host=lckyuh;;
     (hexagon.bccs.uib.no)   local_address=129.177.20.113; local_host=lcxe6;;
     (hx*)                   local_address=133.3.51.11;    local_host=lckyoto;;
     (inferno)               local_address=130.75.105.5;   local_host=lcmuk;;
     (irifi)                 local_address=130.75.105.104; local_host=lcmuk;;
     (jaboticaba)            local_address=150.163.25.181; local_host=lcbr;;
     (sno)                   local_address=130.75.105.113; local_host=lcmuk;;
     (levanto)               local_address=130.75.105.45;  local_host=lcmuk;;
     (login*)                local_address=118.128.66.201; local_host=lckiaps;;
     (maestro)               local_address=130.75.105.2;   local_host=lcmuk;;
     (meller)                local_address=134.106.74.155; local_host=lcfor;;
     (meteo-login*)          local_address=193.166.211.144;local_host=lcxt5m;;
     (mlogin1*|m1*)          local_address=136.172.50.13;  local_host=lcbullhh;;
     (hexagon*)		     local_address=129.177.20.113; local_host=lcxe6;;
     (nobel*)                local_address=150.183.5.101;  local_host=ibms;;
     (ocean)                 local_address="ocean";        local_host=lcocean;; 
     (orkan)                 local_address=130.75.105.3;   local_host=lcmuk;;
     (ostria)                local_address=130.75.105.106; local_host=lcmuk;;
     (paesano)               local_address=130.75.105.46;  local_host=lcmuk;;
     (pcj*)                  local_address=172.31.120.1;   local_host=lckyut;;
     (pingui)                local_address=134.106.74.118; local_host=lcfor;;
     (quanero)               local_address=130.75.105.107; local_host=lcmuk;;
     (rte*)                  local_address=133.5.185.60;   local_host=lcrte;;
     (schultzl-Latitude-E6540)  local_address="schultzl-Latitude-E6540"; local_host=lcsch;;
     (shiokaze-lx)           local_address=134.106.74.123; local_host=lcfor;;
     (sisu-login*)           local_address=86.50.166.21;   local_host=lccrayf;;
     (solano)                local_address=130.75.105.110; local_host=lcmuk;;
     (sugoka*)               local_address=172.31.120.1;   local_host=lckyut;;
     (tc*)                   local_address="ocean";        local_host=lcocean;; 
     (t2a*)                  local_address=10.1.6.165;     local_host=lctit;;
     (urban*)                local_address=147.46.30.151   local_host=lcsb;;
     (vinessa)               local_address=130.75.105.112; local_host=lcmuk;;
     (vorias)                local_address=172.20.25.43;   local_host=lcmuk;;
     (*.cc.kyushu-u.ac.jp)   local_address=133.5.4.129;    local_host=ibmku;;
     (*)                     printf "\n  +++ \"$local_host\" unknown";
                             printf "\n      please contact the PALM group at IMUK";
                             locat=parameter; exit;;
 esac



    # BY DEFAULT, THE REMOTE HOST IS THE LOCAL HOST
 remote_host=$local_host




    # READ THE SHELLSCRIPT-OPTIONS
 while  getopts  :c:dDe:g:h:m:n:N:O:P:q:t:T:u:vX:  option
 do
   case  $option  in
       (c)   job_catalog=$OPTARG;;
       (d)   delete_dayfile=true;;
       (D)   no_submit=true;;
       (e)   email_notification=$OPTARG;;
       (g)   group_number=$OPTARG;;
       (h)   remote_host=$OPTARG;;
       (m)   memory=$OPTARG;;
       (n)   job_name=$OPTARG;;
       (N)   node_usage=$OPTARG;;
       (O)   threads_per_task=$OPTARG;;
       (P)   scp_port=$OPTARG;;
       (q)   no_default_queue=$OPTARG;;
       (t)   cputime=$OPTARG;;
       (T)   tasks_per_node=$OPTARG;;
       (u)   remote_user=$OPTARG;;
       (v)   verify=false;;
       (X)   numprocs=$OPTARG;;
       (\?)  printf "\n  +++ Option $OPTARG unknown \n";
             locat=option; exit;;
   esac
 done


    # GET THE NAME OF THE JOBFILE AS NEXT ARGUMENT
 (( to_shift = $OPTIND - 1 ))
 shift $to_shift; file_to_send=$1


    # OUTPUT OF SHORT DESCRIPTION OF SCRIPT-OPTIONS
 if [ "$1" = "?" ]
 then
   (printf "\n  *** subjob can be called as follows:\n"
    printf "\n      subjob -c.. -d -D -h.. -m.. -q.. -t.. -u.. -v  <jobfile>\n"
    printf "\n      Description of available options:\n"
    printf "\n      Option  Description                         Default-Value"
    printf "\n        -c    job-input- and output-catalog       ~/job_queue"
    printf "\n        -d    no job-protocol will be created     ---"
    printf "\n        -D    only the job-file will be created   ---"
    printf "\n        -h    execution host, available hosts:    $remote_host"
    printf "\n              ibm, ibmh, ibmkisti, ibmku, ibms, lc...,"
    printf "\n              lckiaps, lctit, nech"
    printf "\n        -m    memory demand per process in MByte  ---"
    printf "\n        -n    jobname                             <jobdatei>"
    printf "\n        -O    threads per task (for OpenMP usage) 1"
    printf "\n        -P    ssh/scp port                        default port"
    printf "\n        -q    job-queue to be used                default"
    printf "\n        -t    allowed cpu-time in seconds         ---"
    printf "\n        -T    tasks per node (on parallel hosts)  ---"
    printf "\n        -u    username on execution host          from .netrc"
    printf "\n        -v    no prompt for confirmation          ---"
    printf "\n        -X    # of processors (on parallel hosts) 1"
    printf "\n "
    printf "\n      The only possible positional parameter is <jobfile>:"
    printf "\n      The complete NQS-job must be provided here."
    printf "\n      <jobfile>=? creates this outline\n\n") | more
    exit
 fi



    # CHECK, IF JOB-FILE HAS BEEN GIVEN AS ARGUMENT AND IF THE FILE ITSELF EXISTS
 if [[ "$file_to_send" = "" ]]
 then
    printf "\n  +++ job-file missing"
    locat=parameter; exit
 else
    if [[ -f $file_to_send ]]
    then
       true
    else
       printf "\n  +++ job-file: "
       printf "\n           $file_to_send"
       printf "\n      does not exist"
       locat=parameter; exit
    fi
 fi



    # IF NO JOBNAME HAS BEEN GIVEN, JOBNAME IS SET TO THE NAME OF THE JOB-FILE,
    # PROVIDED THAT THE JOB-FILE NAME DOES NOT CONTAIN ANY PATH
 if [[ $job_name = none ]]
 then
    job_name=$file_to_send
 fi
 if [[ $(echo $job_name | grep -c "/") != 0 ]]
 then
    printf "\n  +++ job-file name: "
    printf "\n           $job_name"
    printf "\n      must not contain \"/\"-characters"
    locat=parameter; exit
 fi




    # SET HOST-SPECIFIC QUANTITIES, OR TERMINATE IN CASE OF UNKNOWN HOST,
    # OR IF NO HOST HAS BEEN GIVEN
 if [[ $remote_host = none ]]
 then
    printf "\n  +++ host missing"
    locat=option; exit
 else
    case  $remote_host  in
        (ibm)     queue=p690_standard; remote_address=134.76.99.81; submcom=/usr/lpp/LoadL/full/bin/llsubmit;;
        (ibmh)    queue=cluster; remote_address=136.172.40.15; submcom=/usr/lpp/LoadL/full/bin/llsubmit;;
        (ibmkisti) queue=class.32plus; remote_address=150.183.146.24; submcom=/usr/lpp/LoadL/full/bin/llsubmit;;
        (ibmku)   queue=s4; remote_address=133.5.4.129; submcom=/usr/local/bin/llsubmit;;
        (ibms)    queue=p_normal; remote_address=150.183.5.101; submcom=/usr/lpp/LoadL/full/bin/llsubmit;;
        (lcbullhh)    queue=compute; remote_address=136.172.50.13; submcom=/usr/bin/sbatch;;
        (lccrayb) queue=mpp1testq; remote_address=130.73.233.1; submcom=/opt/moab/default/bin/msub;;
        (lccrayh) queue=mpp1testq; remote_address=130.75.4.1; submcom=/opt/moab/default/bin/msub;;
        (lccrayf) queue=small; remote_address=86.50.166.21; submcom=/opt/slurm/default/bin/sbatch;;
        (lcflow)  remote_address=10.140.1.72; submcom=qsub;;
        (lckyoto) remote_address=133.3.51.11; submcom=/thin/local/bin/qsub;;
        (lck)     remote_address=165.132.26.61; submcom=/usr/torque/bin/qsub;;
        (lckiaps) remote_address=118.128.66.201; submcom=/opt/pbs/default/bin/qsub;;
        (lckordi) remote_address=210.219.61.8; submcom=/usr/torque/bin/qsub;;
        (lckyuh)  remote_address=133.5.4.33; submcom=/usr/bin/pjsub;;
        (lckyut)  remote_address=133.5.4.37; submcom=/usr/bin/pjsub;;
        (lcocean) remote_address="ocean"; submcom=qsub;;
        (lcsb)    remote_address=147.46.30.151; submcom=/usr/torque/bin/qsub;;
        (lctit)   queue=S; remote_address=10.1.6.165; submcom=/opt/pbs/tools/bin/t2sub;;
        (lcxe6)	  remote_address=129.177.20.113; submcom=/opt/torque/default/bin/qsub;;
        (lcxt5m)  remote_address=193.166.211.144; submcom=/opt/pbs/10.1.0.91350/bin/qsub;;
        (lcyon)   remote_address=165.132.26.68; submcom=/usr/torque/bin/qsub;;
        (nech)    qsubmem=memsz_job; qsubtime=cputim_job; remote_address=136.172.44.147; submcom="/usr/local/bin/qsub";;
        (*)       printf "\n  +++ hostname \"$remote_host\" not allowed";
                  locat=parameter; exit;;
    esac
 fi


    # CHECK, IF A VALID QUEUE HAS BEEN GIVEN
 if [[ $no_default_queue != none ]]
 then
    error=false
    ndq=$no_default_queue
    case  $remote_host  in
        (ibm)    case  $ndq  in
                     (p690_express|p690_standard|p690_long)  error=false;;
                     (*)                                     error=true;;
                 esac;;
        (ibmh)   case  $ndq  in
                     (cluster|express)  error=false;;
                     (*)                                     error=true;;
                 esac;;
        (ibmkisti)   case  $ndq  in
                     (class.32plus|class.1-2|class.2-32)  error=false;;
                     (*)                                     error=true;;
                 esac;;
        (ibmku)  case  $ndq  in
                     (sdbg1|sdbg2|sdbg4|s4|s16|s32|s32-s)    error=false;;
                     (*)                                     error=true;;
                 esac;;
        (ibms)   case  $ndq  in
                     (express|normal|p_express|p_normal|p_normal_1.3|p_normal_1.7|grand)     error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lcbullhh) case  $ndq  in
                     (compute|shared)  error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lccrayb) case  $ndq  in
                     (dataq|mpp1q|mpp1testq|mpp2q|mpp2testq|smp1q|smp1testq|specialm1q)   error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lccrayh) case  $ndq  in
                     (dataq|mpp1q|mpp1testq|mpp2q|mpp2testq|smp1q|smp1testq|specialm1q)   error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lccrayf) case  $ndq  in
                     (usup|test*|small|large)                error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lcflow) case  $ndq  in
                     (cfd_lom_long.q|cfd_him_long.q|cfd_lom_serl.q|cfd_lom_shrt.q|cfd_him_shrt.q|cfd_ivy_shrt.q)  error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lckiaps) case  $ndq  in
                     (express|normal|normal20|quickq)        error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lckyoto) case  $ndq  in
                     (eh|ph)                                 error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lckyuh) case  $ndq  in
                     (fx-dbg|fx-single|fx-small|fx-middle|fx-large)  error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lckyut) case  $ndq  in
                     (cx-dbg|cx-single|cx-small|cx-middle|cx-large)  error=false;;
                     (*)                                     error=true;;
                 esac;;
        (lctit)  case  $ndq  in
                     (G|L128|L256|L512H|S|S96|V)             error=false;;
                     (*)                                     error=true;;
                 esac;;
        (t3eb)   case  $ndq  in
                     (berte|p50|p100|p392|forfree|p25himem)  error=false;;
                     (*)    error=true;;
                 esac;;
        (t3eh)   case  $ndq  in
                     (para_t3e|em|k|l|lm|comp_t3e|c|p|ht)  error=false;;
                     (*)    error=true;;
                 esac;;
        (t3ej2|t3ej5)  case  $ndq  in
                     (low|normal|high)  error=false;;
                     (*)    error=true;;
                 esac;;
        (t3es)  case  $ndq  in
                     (batch|serial-4|pe4|p48|pe16|pe32|pe64|pe128)  error=false;;
                     (*)    error=true;;
                 esac;;
    esac
    if [[ $error = true ]]
    then
       printf "\n  +++ queue \"$no_default_queue\" on host \"$remote_host\" not allowed"
       locat=parameter; exit
    else
       queue=$no_default_queue
    fi
 fi



    # CHECK THE CPU-TIME
    # SPLIT TIME INTO HOURS, MINUTES, AND SECONDS
 done=false
 while [[ $done = false ]]
 do
    if (( $cputime <= 0 ))
    then
       printf "\n  +++ wrong cpu-time or cpu-time missing"
       printf "\n  >>> Please type cpu-time in seconds as INTEGER:"
       printf "\n  >>> "
       read  cputime  1>/dev/null  2>&1
    else
       done=true
    fi
 done
 if [[ $remote_host = nech ]]
 then
    if (( tasks_per_node != 0 ))
    then
       (( cputime = cputime * tasks_per_node ))
    elif [[ $numprocs != 0 ]]
    then
       (( cputime = cputime * numprocs ))
    fi
 fi
 (( stunden  = cputime / 3600 ))
 (( resttime = cputime - stunden * 3600 ))
 (( minuten  = resttime / 60 ))
 (( sekunden = resttime - minuten * 60 ))
 timestring=${stunden}:${minuten}:${sekunden}



    # CHECK THE MEMORY DEMAND
 done=false
 while [[ $done = false ]]
 do
    if (( memory <= 0 ))
    then
       printf "\n  +++ wrong memory demand or memory demand missing"
       printf "\n  >>> Please type memory in  MByte per process  as INTEGER:"
       printf "\n  >>> "
       read  memory  1>/dev/null  2>&1
    else
       done=true
    fi
 done

 if [[ $remote_host = nech ]]
 then
    if (( tasks_per_node != 0 ))
    then
       (( Memory = memory * tasks_per_node / 1000 ))
    elif [[ $numprocs != 0 ]]
    then
       (( Memory = memory * numprocs / 1000 ))
    else
       (( Memory = memory / 1000 ))
    fi
 elif [[ $remote_host = lctit ]]
 then
    (( Memory = memory * tasks_per_node / 1000 ))
 fi


    # MEMORY DEMAND IN CASE OF OPENMP-USAGE ON IBM-SYSTEMS
 if [[ $(echo $remote_host | cut -c1-3) = ibm ]]
 then
    (( memory = memory * threads_per_task ))
 fi


    # CALCULATE NUMBER OF REQUIRED NODES
 if (( tasks_per_node != 0 ))
 then
    (( nodes = ( numprocs - 1 ) / ( tasks_per_node * threads_per_task ) + 1 ))
 fi


    # CALCULATE NUMBER OF PROCESSES PER NODE
 (( processes_per_node = tasks_per_node * threads_per_task ))


    # CALCULATE NUMBER OF MPI TASKS
 (( mpi_tasks = numprocs / threads_per_task ))


    # SET PORT NUMBER OPTION FOR CALLS OF ssh/scp, subjob AND batch_scp SCRIPTS
 if [[ "$scp_port" != "" ]]
 then
    PORTOPT="-P $scp_port"
    SSH_PORTOPT="-p $scp_port"
 fi


    # HEADER-OUTPUT
 if [[ $verify = true ]]
 then
    printf "\n\n"
    printf "#--------------------------------------------------------------# \n"
    spalte1=SUBJOB;spalte2=$(date)
    printf "| %-20s%40s | \n" "$spalte1" "$spalte2"
    printf "|                                                              | \n"
    printf "| values of parameters/options:                                | \n"
    spalte1=$(echo local_host$punkte | cut -c-20)
    spalte2=$punkte$local_host
    printf "| %-20s%40s | \n" "$spalte1" "${spalte2: -40}"
    spalte1=$(echo remote_host$punkte | cut -c-20)
    spalte2=$punkte$remote_host
    printf "| %-20s%40s | \n" "$spalte1" "${spalte2: -40}"
    spalte1=$(echo queue$punkte | cut -c-20)
    spalte2=$punkte$queue
    printf "| %-20s%40s | \n" "$spalte1" "${spalte2: -40}"
    spalte1=$(echo memory$punkte | cut -c-20)
    spalte2="$punkte$memory mb"
    printf "| %-20s%40s | \n" "$spalte1" "${spalte2: -40}"
    spalte1=$(echo cputime$punkte | cut -c-20)
    spalte2="$punkte$cputime sec"
    printf "| %-20s%40s | \n" "$spalte1" "${spalte2: -40}"
    spalte1=$(echo job_name$punkte | cut -c-20)
    spalte2="$punkte$job_name"
    printf "| %-20s%40s | \n" "$spalte1" "${spalte2: -40}"
    printf "#--------------------------------------------------------------# \n\n"


       # QUERY CHECK
    antwort="dummy"
    while [[ $antwort != y  &&  $antwort != Y  &&  $antwort != n  &&  $antwort != N ]]
    do
       read antwort?" >>> continue (y/n) ? "
    done
    if [[ $antwort = n  ||  $antwort = N ]]
    then
       locat=verify; exit
    fi
    printf "\n"
 fi

    # GENERATE RANDOM IDENTIFIER, AND DETERMINE THE JOBNAME ON THE TARGET HOST
 identifier=$RANDOM
 job_on_remhost=${job_name}_${identifier}_$local_host
 job_to_send=job_to_send_$identifier
 if [[ $delete_dayfile = false ]]
 then
    remote_dayfile=${local_host}_${job_name}_result_$identifier
    local_dayfile=${remote_host}_${job_name}
 else
    remote_dayfile=/dev/null
 fi


    # GENERATE THE BATCH-JOB SCRIPTS (FOR QUEUEING-SYSTEMS qsub/msub/LoadLeveler)
 if [[ $(echo $remote_host | cut -c1-3) = ibm  &&  $numprocs != 0 ]]
 then

       # GENERAL LOADLEVELER SETTINGS
    execute_in_shell="#!/bin/ksh"
    use_shell="# @ shell = /bin/ksh"
    consumable_memory="ConsumableMemory($memory mb)"
    class="# @ class = $queue"
    environment="# @ environment = OMP_NUM_THREADS=$threads_per_task; MP_SHARED_MEMORY=yes"
    network_to_use="# @ network.mpi = sn_all,shared,us"
    data_limit="# @ data_limit = 1.76gb"
    image_size="# @ image_size = 50"
    wall_clock_limit="# @ wall_clock_limit = ${timestring},$timestring"

    if [[ $email_notification = none ]]
    then
       notify_user=""
    else
       notify_user="# @ notify_user = $email_notification"
       if [[ $delete_dayfile = true ]]
       then
          notification='# @ notification = never'
       fi
    fi

    if [[ $remote_host = ibmh ]]
    then
       data_limit=""
       network_to_use=""
       class="# @ class = $queue"
       environment=""
       rset="# @ rset = RSET_MCM_AFFINITY"
       task_affinity="# @ task_affinity = core(1)"
    elif [[ $remote_host = ibmkisti ]]
    then
       network_to_use="# @ network.MPI = sn_all,shared,US"
       wall_clock_limit="# @ wall_clock_limit = $timestring"
       if [[ $threads_per_task = 1 ]]
       then
          rset="# @ rset = RSET_MCM_AFFINITY"
          mcm_affinity_options="# @ mcm_affinity_options = mcm_mem_pref mcm_sni_none mcm_distribute"
       fi
       environment=""
       use_shell=""
       data_limit=""
       image_size=""
    elif [[ $remote_host = ibmku ]]
    then
       execute_in_shell="#!/usr/bin/ksh"
       use_shell="# @ shell = /usr/bin/ksh"
       consumable_memory=""
       environment=""
       network_to_use="# @ network.mpi = sn_all,shared,us"
       data_limit=""
       image_size=""
    elif [[ $remote_host = ibms ]]
    then
       network_to_use="# @ network.mpi = csss,shared,us"
    fi

    cat > $job_to_send << %%END%%
$execute_in_shell
$use_shell

# @ job_type = parallel
# @ job_name = $job_name
# @ resources = ConsumableCpus($threads_per_task) $consumable_memory
# @ output = $remote_dayfile
# @ error = $remote_dayfile
$wall_clock_limit
$image_size
$class
$environment
$network_to_use
$data_limit
$rset
$mcm_affinity_options
$task_affinity
$notification
$notify_user

%%END%%

    if (( nodes > 0 ))
    then

       if [[ $remote_host != ibmkisti ]]
       then

          cat >> $job_to_send << %%END%%
# @ node = $nodes
# @ tasks_per_node = $processes_per_node
# @ node_usage = $node_usage
# @ queue

%%END%%

       else

          cat >> $job_to_send << %%END%%
# @ total_tasks = $mpi_tasks
# @ blocking = unlimited
# @ queue

%%END%%

       fi

    else

       cat >> $job_to_send << %%END%%
# @ blocking = unlimited
# @ total_tasks = $numprocs
# @ node_usage = $node_usage
# @ queue

%%END%%

    fi

       # WORKAROUND BECAUSE OF SILLY JOB FILTER ON ibmkisti
    if [[ $remote_host = ibmkisti  &&  $threads_per_task != 1 ]]
    then
       echo  "export OMP_NUM_THREADS=$threads_per_task"  >>  $job_to_send
    fi

 elif [[ $(echo $remote_host | cut -c1-3) = ibm  &&  $numprocs = 0 ]]
 then

    cat > $job_to_send << %%END%%
#!/bin/ksh

# @ job_type = serial
# @ node_usage = $node_usage
# @ job_name = palm
# @ wall_clock_limit = ${timestring},$timestring
# @ resources = ConsumableCpus(1) ConsumableMemory(1 gb)
# @ output = $remote_dayfile
# @ error = $remote_dayfile
$class
$notification

# @ queue

%%END%%

 elif [[ $remote_host = lcbullhh ]]
 then
    if [[ $numprocs != 0 ]]
    then
       cat > $job_to_send << %%END%%
#!/bin/bash -l
#SBATCH -J $job_name
#SBATCH -t $timestring
#SBATCH -N $nodes
#SBATCH --ntasks-per-node=$processes_per_node
#SBATCH -p $queue
#SBATCH -o $remote_dayfile
#SBATCH -e $remote_dayfile
#SBATCH -A $project_account

$init_cmds
$module_calls

%%END%%

    else
       cat > $job_to_send << %%END%%
#!/bin/bash -l
#SBATCH -J $job_name
#SBATCH -t $timestring
#SBATCH -l ncpus=1
#SBATCH -l pmem=${memory}mb
#SBATCH -m abe
#SBATCH -o $remote_dayfile
#SBATCH -e $remote_dayfile
#SBATCH -A $project_account

$init_cmds
$module_calls

%%END%%

    fi

 elif [[ $remote_host = lccrayb || $remote_host = lccrayh ]]
 then

    if [[ "$feature" != "" ]]
    then
       featuredir="#PBS -l feature=$feature"
    fi

    if [[ $numprocs != 0 ]]
    then
       cat > $job_to_send << %%END%%
#!/bin/bash
#PBS -N $job_name
#PBS -l walltime=$timestring
#PBS -l nodes=$nodes:ppn=$processes_per_node
#PBS -o $remote_dayfile
#PBS -j oe
#PBS -q $queue
$featuredir

$init_cmds
$module_calls

%%END%%

    else

       continue

    fi

 elif [[ $remote_host = lccrayf ]]
 then

    if [[ $numprocs != 0 ]]
    then
       cat > $job_to_send << %%END%%
#!/bin/bash -l
#SBATCH -J $job_name
#SBATCH -t $timestring
#SBATCH -N $nodes
#SBATCH --ntasks-per-node=$processes_per_node
#SBATCH -p $queue
#SBATCH -o $remote_dayfile
#SBATCH -e $remote_dayfile

$init_cmds
$module_calls

%%END%%

    else
       cat > $job_to_send << %%END%%
#!/bin/bash -l
#SBATCH -J $job_name
#SBATCH -t $timestring
#SBATCH -l ncpus=1
#SBATCH -l pmem=${memory}mb
#SBATCH -m abe
#SBATCH -o $remote_dayfile
#SBATCH -e $remote_dayfile

$init_cmds
$module_calls

%%END%%

    fi

 elif [[ $remote_host = lcflow ]]
 then
    if [[ $numprocs != 0 ]]
    then
      pe_set="#$ -pe impi $numprocs"
    else
      pe_set="#$ -pe impi 1"
    fi
    if [[ $queue = default ]]
    then
      queue_set=""
    else
      queue_set="#$ -q $queue"
    fi
    [[ "$disc_space" = "" ]]  &&  disc_space=50

       cat > $job_to_send << %%END%%
#!/bin/bash
#$ -S /bin/bash
#$ -N $job_name
#$ -cwd
#$ -l h_rt=$timestring
#$ -l h_vmem=${memory}M
#$ -o $remote_dayfile
#$ -j y
$pe_set
#$ -R y
#$ -l h_fsize=${disc_space}G
$queue_set

%%END%%

 elif [[ $remote_host = lck || $remote_host = lckordi || $remote_host = lcsb ]]
 then

    if [[ $numprocs != 0 ]]
    then
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -N $job_name
#PBS -l walltime=$timestring
#PBS -l ncpus=$numprocs
#PBS -l pmem=${memory}mb
#PBS -o $remote_dayfile
#PBS -l nodes=$nodes:ppn=${processes_per_node}
#PBS -j oe

mpd &

%%END%%

    else
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -N $job_name
#PBS -l walltime=$timestring
#PBS -l ncpus=1
#PBS -l pmem=${memory}mb
#PBS -o $remote_dayfile
#PBS -j oe

%%END%%

    fi

 elif [[ $remote_host = lckiaps ]]
 then

    if [[ $numprocs != 0 ]]
    then
       cat > $job_to_send << %%END%%
#!/bin/bash
#PBS -N $job_name
#PBS -l walltime=$timestring
#PBS -l select=1:ncpus=$numprocs
#PBS -l pmem=${memory}mb
#PBS -q $queue
#PBS -o $remote_dayfile
#PBS -j oe
#PBS -V

%%END%%

    else
       cat > $job_to_send << %%END%%
#!/bin/bash
#PBS -N $job_name
#PBS -l walltime=$timestring
#PBS -l ncpus=1
#PBS -l pmem=${memory}mb
#PBS -o $remote_dayfile
#PBS -j oe

%%END%%

    fi

 elif [[ $remote_host = lcyon ]]
 then

    if [[ $numprocs != 0 ]]
    then
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -N $job_name
#PBS -l walltime=$timestring
#PBS -l ncpus=$numprocs
#PBS -l pmem=${memory}mb
#PBS -o $remote_dayfile
#PBS -j oe

%%END%%

    else
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -N $job_name
#PBS -l walltime=$timestring
#PBS -l ncpus=1
#PBS -l pmem=${memory}mb
#PBS -o $remote_dayfile
#PBS -j oe

%%END%%

    fi

 elif [[ $remote_host = lcxe6 ]]
 then

    if [[ $numprocs != 0 ]]
    then
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -S /bin/ksh
#PBS -N $job_name
#PBS -A $project_account
#PBS -j oe
#PBS -l walltime=$timestring
#PBS -l mppwidth=${numprocs}
#PBS -l mppnppn=${processes_per_node}
#PBS -m abe
#PBS -o $remote_dayfile
$email_directive

$init_cmds
$module_calls

%%END%%

    else
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -S /bin/ksh
#PBS -N $job_name
#PBS -A $project_account
#PBS -j oe
#PBS -l walltime=$timestring
#PBS -l ncpus=1
#PBS -l pmem=${memory}mb
#PBS -m abe
$email_directive
#PBS -o $remote_dayfile

$init_cmds
$module_calls

%%END%%

    fi

 elif [[ $remote_host = lckyoto ]]
 then

       cat > $job_to_send << %%END%%
#!/bin/ksh
# @\$-o $remote_dayfile
# @\$-eo -oi
# @\$-lP 16
# @\$-lp 1
# @\$-lm 28gb  -llm unlimited -ls unlimited
# @\$-q $queue
# @\$-Pvn abs_pack
##for intel? @\$-Pvn abs_unpack -Pvs unpack -Pvc unpack
#. /thin/local/etc/setprofile/intel-11.0.sh
#. /thin/local/etc/setprofile/mvapich2-1.4+intel-11.0.sh
. ~/.myprofile
#. /home2/t/t51254/palm/current_version/myprofile
#. /thin/apps/pgi/mpi.sh
#
env
#
set -x

%%END%%

 elif [[ $remote_host = lcxt5m ]]
 then

    if [[ $numprocs != 0 ]]
    then
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -S /bin/ksh
#PBS -N $job_name
#PBS -j oe
#PBS -l walltime=$timestring
#PBS -l mppwidth=${numprocs}
#PBS -l mppnppn=${processes_per_node}
#PBS -m abe
#PBS -o $remote_dayfile

$init_cmds
$module_calls

%%END%%

    else
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -S /bin/ksh
#PBS -N $job_name
#PBS -j oe
#PBS -l walltime=$timestring
#PBS -l ncpus=1
#PBS -l pmem=${memory}mb
#PBS -m abe
#PBS -o $remote_dayfile

$init_cmds
$module_calls

%%END%%

    fi

 elif [[ $remote_host = lckyuh ]]
 then
    cat > $job_to_send << %%END%%
#!/bin/bash
#PJM -L "rscgrp=$queue"
#PJM -L "node=$nodes"
#PJM --mpi "proc=$numprocs"
#PJM -L "elapse=$timestring"
#PJM -o $remote_dayfile
#PJM -j
#PJM -X
#PJM --no-stging

export LANG=en_US.UTF-8
%%END%%

 elif [[ $remote_host = lckyut ]]
 then
    cat > $job_to_send << %%END%%
#!/bin/bash
#PJM -L "rscgrp=$queue"
#PJM -L "vnode=$numprocs"
#PJM -L "vnode-core=1"
#PJM -L "elapse=$timestring"
#PJM --mpi proc=$numprocs
#PJM -o $remote_dayfile
#PJM -j
#PJM -X
#PJM --no-stging

export LANG=en_US.UTF-8
%%END%%

 elif [[ $remote_host = lcocean ]]
 then
   cat > $job_to_send << %%END%%
#!/bin/bash
#$ -cwd
#$ -V 
#$ -N $job_name
#$ -pe orte $numprocs
#$ -o $remote_dayfile
#$ -j y
#$ -R y
$init_cmds
$module_calls

%%END%%

 elif [[ $remote_host = nech ]]
 then

    if (( nodes > 1 ))
    then
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -l cpunum_prc=$processes_per_node,cputim_job=$cputime
#PBS -l ${qsubmem}=${Memory}gb
#PBS -b $nodes
#PBS -o $remote_dayfile
#PBS -N palm
#PBS -j o
#PBS -T mpisx

%%END%%

    elif [[ $numprocs != 0 ]]
    then
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -l cpunum_prc=$processes_per_node,cputim_job=$cputime
#PBS -l ${qsubmem}=${Memory}gb
#PBS -o $remote_dayfile
#PBS -N palm
#PBS -j o

%%END%%

    else
       cat > $job_to_send << %%END%%
#!/bin/ksh
#PBS -l ${qsubmem}=${Memory}gb,${qsubtime}=$cputime
#PBS -o $remote_dayfile
#PBS -j o

%%END%%

    fi

 elif [[ $remote_host = lctit ]]
 then
    cat > $job_to_send << %%END%%
#!/bin/ksh
$init_cmds
$module_calls

%%END%%

       # SET OPTIONS FOR SUBMIT-COMMAND
    if [[ $tasks_per_node != $processes_per_node ]]
    then
       submcom="$submcom -W group_list=$group_number -N $job_name -l walltime=$timestring -l select=$nodes:ncpus=$processes_per_node:mpiprocs=$tasks_per_node:mem=${Memory}gb -o $remote_dayfile -j oe -et 1 -q $queue "
    else
       submcom="$submcom -W group_list=$group_number -N $job_name -l walltime=$timestring -l select=$nodes:ncpus=$processes_per_node:mpiprocs=$tasks_per_node:mem=${Memory}gb -l place=scatter -o $remote_dayfile -j oe -et 1 -q $queue "
    fi

 else

    cat > $job_to_send << %%END%%
# @\$-q ${queue}
# @\$-l${qsubtime} $timestring
# @\$-l${qsubmem} ${memory}mb
# @\$-o $remote_dayfile
# @\$-eo

%%END%%

 fi


    # IN CASE OF JOBS EXECUTING ON REMOTE-HOSTS, THE TRANSFER OF THE DAYFILES
    # TO THE LOCAL HOSTS WILL BE INITIATED BY TRAP ON EXIT
    # NO TRANSFER POSSIBLE ON IBM IN SEOUL
 if [[ $delete_dayfile = false  &&  $remote_host != $local_host ]]
 then
    echo "set +vx"                              >>  $job_to_send
    echo "trap '"                               >>  $job_to_send
    echo "set +vx"                              >>  $job_to_send
    if [[ $(echo $remote_host | cut -c1-3) = ibm  ||  $remote_host = lcbullhh  ||  $remote_host = lccrayb  ||  $remote_host = lccrayh  ||  $(echo $remote_host | cut -c1-3) = nec  ||  $remote_host = lcflow  ||  $remote_host = lckiaps  ||  $remote_host = lckyu* || $remote_host = lcxe6  ||  $remote_host = lcocean ]]
    then
       if [[ $remote_host = ibmh ]]
       then
          return_queue=c1
       elif [[ $remote_host = ibmkisti ]]
       then
          return_queue=class.1-2
       elif [[ $remote_host = ibmku ]]
       then
          return_queue=sdbg2
       elif [[ $remote_host = ibms ]]
       then
          return_queue=p_normal
       elif [[ $remote_host = lcbullhh ]]
       then
          return_queue=shared
       elif [[ $remote_host = lccrayb || $remote_host = lccrayh ]]
       then
          return_queue=dataq
       elif [[ $remote_host = lcxe6 ]]
       then
          return_queue=debug
       elif [[ $remote_host = lckiaps ]]
       then
          return_queue=express
       elif [[ $remote_host = lckyuh ]]
       then
          return_queue=cx-single
       elif [[ $remote_host = lckyut ]]
       then
          return_queue=cx-single
       else
          return_queue=unknown
       fi

       if [[ $(echo $remote_host | cut -c1-3) = ibm ]]
       then

          if [[ $remote_host = ibmku ]]
          then
             echo "echo \"#!/usr/bin/ksh\" >> scpjob.$identifier"            >>  $job_to_send
             echo "echo \"# @ shell = /usr/bin/ksh\" >> scpjob.$identifier"  >>  $job_to_send
          else
             echo "echo \"#!/bin/ksh\" >> scpjob.$identifier"                >>  $job_to_send
          fi
          echo "echo \"# @ job_type = serial\" >> scpjob.$identifier"    >>  $job_to_send
          echo "echo \"# @ job_name = transfer\" >> scpjob.$identifier"  >>  $job_to_send
          echo "echo \"# @ resources = ConsumableCpus(1) ConsumableMemory(1 gb)\" >> scpjob.$identifier"  >>  $job_to_send
          echo "echo \"# @ wall_clock_limit = 00:10:00,00:10:00\" >> scpjob.$identifier "  >>  $job_to_send
          echo "echo \"# @ output = job_queue/last_job_transfer_protocol\" >> scpjob.$identifier"  >>  $job_to_send
          echo "echo \"# @ error = job_queue/last_job_transfer_protocol\" >> scpjob.$identifier"  >>  $job_to_send
          if [[ $host != "ibmh" ]]
          then
             echo "echo \"# @ class = $return_queue\" >> scpjob.$identifier"  >>  $job_to_send
          fi
          echo "echo \"# @ image_size = 10\" >> scpjob.$identifier"      >>  $job_to_send
          echo "echo \"# @ notification = never\" >> scpjob.$identifier" >>  $job_to_send

          echo "echo \"# @ queue\" >> scpjob.$identifier"                >>  $job_to_send
          echo "echo \" \" >> scpjob.$identifier"                        >>  $job_to_send

          echo "echo \"set -x\" >> scpjob.$identifier"                   >>  $job_to_send
          echo "echo \"batch_scp  $PORTOPT  -d  -w 10  -u $local_user  $local_address  ${job_catalog}/$remote_dayfile  \\\"$job_catalog\\\"  $local_dayfile\" >> scpjob.$identifier"  >>  $job_to_send
          if [[ $remote_host = ibmku ]]
          then
             echo "echo \"rm  scpjob.$identifier\" >> scpjob.$identifier"   >>  $job_to_send
          fi
          echo "echo \"exit\" >> scpjob.$identifier"                     >>  $job_to_send

       elif [[ $remote_host = nech ]]
       then
          echo "cd /pf/b/${remote_user}/job_queue" >>  $job_to_send
          echo "cat > scpjob.$identifier << %%END%%"  >>  $job_to_send
          echo "#PBS -l ${qsubmem}=1GB,${qsubtime}=100"  >>  $job_to_send
          echo "#PBS -o last_job_transfer_protocol"      >>  $job_to_send
          echo "#PBS -j o"                         >>  $job_to_send
          echo " "                                 >>  $job_to_send
          echo "set -x"                            >>  $job_to_send
          echo "cd /pf/b/${remote_user}/job_queue" >>  $job_to_send
          echo "batch_scp  $PORTOPT  -d  -w 10  -u $local_user $local_address  $remote_dayfile  \"$job_catalog\"  $local_dayfile"  >>  $job_to_send
          echo "[[ \"\$for_subjob_to_do\" != \"\" ]]  &&  eval \$for_subjob_to_do"  >>  $job_to_send
          echo "%%END%%"                           >>  $job_to_send

       elif [[ $remote_host = lcbullhh ]]
       then
          echo "cat > scpjob.$identifier << %%END%%"        >>  $job_to_send
          echo "#!/bin/bash"                             >>  $job_to_send
          echo "#SBATCH --job-name=job_protocol_transfer" >>  $job_to_send
          echo "#SBATCH -t 00:20:00"                     >>  $job_to_send
          echo "#SBATCH -N 1"                            >>  $job_to_send
          echo "#SBATCH -n 1"                            >>  $job_to_send
          echo "#SBATCH -o \$HOME/job_queue/last_job_transfer_protocol"      >>  $job_to_send
          echo "#SBATCH -o $remote_dayfile"              >>  $job_to_send
          echo "#SBATCH -e $remote_dayfile"              >>  $job_to_send
          echo "#SBATCH -A $project_account"             >>  $job_to_send
          echo "#SBATCH -p $return_queue"                >>  $job_to_send
          echo " "                                       >>  $job_to_send
          echo "set -x"                                  >>  $job_to_send
          echo "batch_scp  $PORTOPT  -d  -w 10  -u $local_user $local_address  ${job_catalog}/$remote_dayfile  \"$job_catalog\"  $local_dayfile"  >>  $job_to_send
          echo "[[ \"\$for_subjob_to_do\" != \"\" ]]  &&  eval \$for_subjob_to_do"  >>  $job_to_send
          echo "%%END%%"                                 >>  $job_to_send

       elif [[ $remote_host = lckyuh ]]
       then
          echo "cat > scpjob.$identifier << %%END%%"  >>  $job_to_send
          echo "#!/bin/bash"                       >>  $job_to_send
          echo "#PJM -L \"node=1\""                >>  $job_to_send
          echo "#PJM -L \"rscgrp=$return_queue\""  >>  $job_to_send
          echo "#PJM --no-stging"                  >>  $job_to_send
          echo "#PJM -L \"elapse=30:00\""          >>  $job_to_send
          echo "#PJM -o \$HOME/job_queue/last_job_transfer_protocol"  >>  $job_to_send
          echo "#PJM -j"                           >>  $job_to_send
          echo " "                                 >>  $job_to_send
          echo "export LANG=en_US.UTF-8"           >>  $job_to_send
          echo "set -x"                            >>  $job_to_send
          echo "batch_scp  $PORTOPT  -d  -w 10  -u $local_user $local_address  $remote_dayfile  \"$job_catalog\"  $local_dayfile"  >>  $job_to_send
          echo "[[ \"\$for_subjob_to_do\" != \"\" ]]  &&  eval \$for_subjob_to_do"  >>  $job_to_send
          echo "%%END%%"                           >>  $job_to_send

       elif [[ $remote_host = lckyut ]]
       then
          echo "cat > scpjob.$identifier << %%END%%"  >>  $job_to_send
          echo "#!/bin/bash"                       >>  $job_to_send
          echo "#PJM -L \"vnode=1\""               >>  $job_to_send
          echo "#PJM -L \"rscgrp=$return_queue\""  >>  $job_to_send
          echo "#PJM --no-stging"                  >>  $job_to_send
          echo "#PJM -L \"elapse=30:00\""          >>  $job_to_send
          echo "#PJM -o \$HOME/job_queue/last_job_transfer_protocol"  >>  $job_to_send
          echo "#PJM -j"                           >>  $job_to_send
          echo " "                                 >>  $job_to_send
          echo "export LANG=en_US.UTF-8"           >>  $job_to_send
          echo "set -x"                            >>  $job_to_send
          echo "batch_scp  $PORTOPT  -d  -w 10  -u $local_user $local_address  $remote_dayfile  \"$job_catalog\"  $local_dayfile"  >>  $job_to_send
          echo "[[ \"\$for_subjob_to_do\" != \"\" ]]  &&  eval \$for_subjob_to_do"  >>  $job_to_send
          echo "%%END%%"                           >>  $job_to_send

       elif [[ $remote_host = lccrayb || $remote_host = lccrayh ]]
       then
          echo "cat > scpjob.$identifier << %%END%%"        >>  $job_to_send
          echo "#!/bin/bash"                             >>  $job_to_send
          echo "#PBS -N job_protocol_transfer"           >>  $job_to_send
          echo "#PBS -l walltime=00:30:00"               >>  $job_to_send
          echo "#PBS -l nodes=1:ppn=1"                   >>  $job_to_send
          echo "#PBS -o \$HOME/job_queue/last_job_transfer_protocol"      >>  $job_to_send
          echo "#PBS -j oe"                              >>  $job_to_send
          echo " "                                       >>  $job_to_send
          echo "set -x"                                  >>  $job_to_send
          echo "batch_scp  $PORTOPT  -d  -w 10  -u $local_user $local_address  ${job_catalog}/$remote_dayfile  \"$job_catalog\"  $local_dayfile"  >>  $job_to_send
          echo "[[ \"\$for_subjob_to_do\" != \"\" ]]  &&  eval \$for_subjob_to_do"  >>  $job_to_send
          echo "%%END%%"                                 >>  $job_to_send

       elif [[ $remote_host = lcocean ]]
       then
          echo "cat > scpjob.${identifier}.tmp << %%END%%"                  >>  $job_to_send
          echo "#!/bin/bash"                                             >>  $job_to_send
          echo "SGEPREFIX -S /bin/bash"                                  >>  $job_to_send
          echo "SGEPREFIX -N transfer_$job_name"                         >>  $job_to_send
          echo "SGEPREFIX -cwd"                                          >>  $job_to_send
          echo "SGEPREFIX -j y"                                          >>  $job_to_send
          echo "SGEPREFIX -o ${local_host}_${job_name}_scpjob_$identifier"  >>  $job_to_send 
          echo " "                                                       >>  $job_to_send 
          echo "set -x"                                                  >>  $job_to_send 
          echo "export PALM_BIN=$PALM_BIN" | sed -e 's:'$HOME':$HOME:'   >>  $job_to_send
          echo "export PATH=\$PATH:\$PALM_BIN"                           >>  $job_to_send
          echo ""                                 >>  $job_to_send         
          echo "batch_scp  $PORTOPT  -d  -w 10  -u $local_user $local_address  ${job_catalog}/$remote_dayfile  \"$job_catalog\"  $local_dayfile"  >>  $job_to_send
          echo "[[ \"\$for_subjob_to_do\" != \"\" ]]  &&  eval \$for_subjob_to_do"  >>  $job_to_send
          echo "rm -f scpjob.${identifier}"                                 >>  $job_to_send         
          echo "%%END%%"                                                 >>  $job_to_send
          echo "sed -e 's/SGEPREFIX/#$/g' scpjob.${identifier}.tmp > scpjob.${identifier}" >>  $job_to_send         
          echo "rm -f scpjob.${identifier}.tmp"                             >>  $job_to_send         

       elif [[ $remote_host = lcflow ]]
       then
          echo "cat > scpjob.${identifier}.tmp << %%END%%"                  >>  $job_to_send
          echo "#!/bin/bash"                                             >>  $job_to_send
          echo "SGEPREFIX -S /bin/bash"                                  >>  $job_to_send
          echo "SGEPREFIX -N transfer_$job_name"                         >>  $job_to_send
          echo "SGEPREFIX -cwd"                                          >>  $job_to_send
          echo "SGEPREFIX -l h_rt=01:00:00"                              >>  $job_to_send
          echo "SGEPREFIX -l h_vmem=500M"                                >>  $job_to_send
          echo "SGEPREFIX -j y"                                          >>  $job_to_send
          echo "SGEPREFIX -o ${local_host}_${job_name}_scpjob_$identifier"  >>  $job_to_send 
          echo " "                                                       >>  $job_to_send 
          echo "set -x"                                                  >>  $job_to_send 
          echo "export PALM_BIN=$PALM_BIN" | sed -e 's:'$HOME':$HOME:'   >>  $job_to_send
          echo "export PATH=\$PATH:\$PALM_BIN"                           >>  $job_to_send
          echo ""                                 >>  $job_to_send         
          echo "batch_scp  $PORTOPT  -d  -w 10  -u $local_user $local_address  ${job_catalog}/$remote_dayfile  \"$job_catalog\"  $local_dayfile"  >>  $job_to_send
          echo "[[ \"\$for_subjob_to_do\" != \"\" ]]  &&  eval \$for_subjob_to_do"  >>  $job_to_send
          echo "rm -f scpjob.${identifier}"                                 >>  $job_to_send         
          echo "%%END%%"                                                 >>  $job_to_send
          echo "sed -e 's/SGEPREFIX/#$/g' scpjob.${identifier}.tmp > scpjob.${identifier}" >>  $job_to_send         
          echo "rm -f scpjob.${identifier}.tmp"                             >>  $job_to_send         
       elif [[ $remote_host = lcxe6 ]]
       then
          echo "cat > scpjob.${identifier}  << %%END%%"  >>  $job_to_send
          echo "#!/bin/ksh"                              >>  $job_to_send
          echo "#PBS -N job_protocol_transfer"           >>  $job_to_send
          echo "#PBS -l walltime=00:30:00"               >>  $job_to_send
          echo "#PBS -A $project_account"                >>  $job_to_send
          echo "#PBS -l mppwidth=1"                      >>  $job_to_send
          echo "#PBS -l mppnppn=1"                       >>  $job_to_send
          echo "#PBS -o \$HOME/job_queue/last_job_transfer_protocol"  >>  $job_to_send
          echo "#PBS -j oe"                              >>  $job_to_send
          echo " "                                       >>  $job_to_send
          echo "set -x"                                  >>  $job_to_send
          echo "batch_scp  $PORTOPT  -d  -w 10  -u $local_user $local_address  ${job_catalog}/$remote_dayfile  \"$job_catalog\"  $local_dayfile"  >>  $job_to_send
          echo "[[ \"\$for_subjob_to_do\" != \"\" ]]  &&  eval \$for_subjob_to_do"  >>  $job_to_send
          echo "%%END%%"                                 >>  $job_to_send
       else

          echo "cat > scpjob.$identifier << %%END%%"  >>  $job_to_send
          echo "# @\\\$-q $return_queue"           >>  $job_to_send
          echo "# @\\\$-l${qsubtime} 10"           >>  $job_to_send
          echo "# @\\\$-l${qsubmem} 10mb"          >>  $job_to_send
          if [[ $remote_host = t3ej2  ||  $remote_host = t3ej5  ||  $remote_host = t3es ]]
          then
             echo "# @\$-l mpp_p=0"                >>  $job_to_send
          fi
          echo '# @\$-lF 10mb'                     >>  $job_to_send
          echo '# @\$-o job_queue/last_job_transfer_protocol'    >>  $job_to_send
          echo '# @\\\$-eo'                          >>  $job_to_send
          echo " "                                 >>  $job_to_send
          if [[ $remote_host = t3ej2  ||  $remote_host = t3ej5 ]]
          then
             echo "set +vx"                        >>  $job_to_send
             echo ". .profile"                     >>  $job_to_send
          fi
          echo "set -x"                            >>  $job_to_send
          echo "batch_scp  $PORTOPT  -d  -w 10  -u $local_user $local_address  ${job_catalog}/$remote_dayfile  \"$job_catalog\"  $local_dayfile  >  /dev/null"  >>  $job_to_send
          echo "[[ \"\$for_subjob_to_do\" != \"\" ]]  &&  eval \$for_subjob_to_do"  >>  $job_to_send
          echo "%%END%%"                           >>  $job_to_send

       fi

       if [[ $(echo $remote_host | cut -c1-3) = ibm ]]
       then
          echo "llsubmit  scpjob.$identifier"      >>  $job_to_send
       elif [[ $remote_host = lcbullhh ]]
       then
          echo "sbatch  scpjob.$identifier"               >>  $job_to_send
       elif [[ $remote_host = lccrayb || $remote_host = lccrayh ]]
       then
          echo "msub -q $return_queue  scpjob.$identifier"               >>  $job_to_send
       elif [[ $remote_host = t3eb  ||  $remote_host = t3eh  ||  $remote_host = t3ej2  ||  $remote_host = t3ej5 ]]
       then
          echo "qsub -J n  scpjob.$identifier"     >>  $job_to_send
       elif [[ $remote_host = t3es ]]
       then
          echo "qsub -J n  -s /bin/ksh  scpjob.$identifier"     >>  $job_to_send
       elif [[ $remote_host = lckiaps ]]
       then
          echo "mv  scpjob.$identifier  $job_catalog"           >>  $job_to_send
          echo "ssh $SSH_PORTOPT ${remote_username}@${remote_address}  \"$submcom ${job_catalog}/scpjob.$identifier\" "  >>  $job_to_send
          echo "rm  ${job_catalog}/scpjob.$identifier"          >>  $job_to_send
       elif [[ $remote_host = lckyu* ]]
       then
          echo "scp $PORTOPT scpjob.$identifier  ${remote_username}@${remote_address}:job_queue"           >>  $job_to_send
          echo "ssh $SSH_PORTOPT ${remote_username}@${remote_address}  \"cd job_queue; $submcom scpjob.$identifier; rm scpjob.$identifier\" "  >>  $job_to_send
       elif [[ $remote_host = lcflow  ||  $remote_host = lcocean ]]
       then
          echo "mv  scpjob.$identifier  $job_catalog"           >>  $job_to_send
          echo "/usr/bin/ssh ${remote_username}@${remote_address}  \"$init_cmds $module_calls cd $job_catalog; $submcom scpjob.$identifier\" "  >>  $job_to_send
       else
          echo "$submcom  scpjob.$identifier"      >>  $job_to_send
       fi
       if [[ $remote_host != ibmku  &&  $remote_host != lckiaps ]]
       then
          echo "rm  scpjob.$identifier"            >>  $job_to_send
       fi
       if [[ $remote_host = nech ]]
       then
          echo "cd -"                           >>  $job_to_send
       fi
    else
#       echo "ftpcopy  -d  $local_address  ${job_catalog}/$remote_dayfile  \"$job_catalog\"  $local_dayfile"  >>  $job_to_send
       # ??? funktioniert das überhaupt noch ???
       echo "nohup  ftpcopy  -d  -w 15  $local_address  ${job_catalog}/$remote_dayfile  \"$job_catalog\"  $local_dayfile  >  /dev/null  &"  >>  $job_to_send
    fi
    echo "set -x"                               >>  $job_to_send
    echo "     ' exit"                          >>  $job_to_send
    echo "set -x"                               >>  $job_to_send
 fi


    # APPEND THE JOB-FILE (CREATE BY mrun) TO THE JOB-DIRECTIVES GENERATED ABOVE
 cat  $file_to_send  >>  $job_to_send

 if [[ $remote_host = ibm ]]
 then
    echo " "         >>  $job_to_send
    echo "exit"      >>  $job_to_send
 fi

    # REMOVE JOB-FILE
 if [[ $remote_host = lctit  ||  $remote_host = ibmku  ||  $remote_host = lcflow ]]
 then
    echo " "                               >>  $job_to_send
    echo "rm ~/job_queue/$job_on_remhost"  >>  $job_to_send
 fi


    # TRANSFER JOB TO THE TARGET HOST (JOB-DIRECTORY)
 if [[ $no_submit = false ]]
 then
    if [[ $remote_host != $local_host ]]
    then
       [[ $verify = true ]]  &&  printf "\n >>> transfering job to \"$remote_host\"..."
       if [[ $remote_host = ibms ]]    # ssh on ibms cannot handle "~/"
       then
          job_catalog_save=$job_catalog
          job_catalog=job_queue
       elif [[ $remote_host = nech ]]
       then
          job_catalog_save=$job_catalog
          job_catalog=/hpf/b/${remote_user}/job_queue
       fi
       if [[ $remote_host = nech ]]
       then
             # FILES CAN ONLY BE TRANSFERED VIA DKRZ'S ARCHIVE-SERVER
          scp  $PORTOPT  $job_to_send  ${remote_user}@136.172.44.205:${job_catalog}/$job_on_remhost
       else
          scp  $ssh_key  $PORTOPT  $job_to_send  ${remote_user}@${remote_address}:${job_catalog}/$job_on_remhost
       fi
       if [[ $? = 1 ]]
       then
          locat=scp; exit
       fi
       if [[ $remote_host = ibms ]]
       then
          job_catalog=$job_catalog_save
       fi
       [[ $verify = true ]]  &&  printf "\n >>> finished\n"
    else
       eval  job_catalog=$job_catalog
       cp  $job_to_send  ${job_catalog}/$job_on_remhost
    fi



       # START NQS- / LOADLEVELER-JOB
    if [[ $remote_host != $local_host ]]
    then
       [[ $verify = true ]]  &&  printf "\n >>> submitting job using \"qsub\"...\n"

       if [[ $remote_host = ibmku ]]
       then
          ssh  $SSH_PORTOPT $remote_address  -l $remote_user  "cd $job_catalog; $submcom $job_on_remhost"
       elif [[ $remote_host = lcflow ]]
       then
          /usr/bin/ssh  $SSH_PORTOPT $remote_address  -l $remote_user  "$init_cmds $module_calls cd $job_catalog; $submcom $job_on_remhost"
       else
          ssh  $ssh_key  $SSH_PORTOPT $remote_address  -l $remote_user  "cd $job_catalog; $submcom $job_on_remhost; rm $job_on_remhost"
       fi

       [[ $verify = true ]]  &&  printf " >>> o.k.\n"
    else
       cd  $job_catalog
       if [[ $(echo $local_host | cut -c1-3) = ibm  ||  $(echo $local_host | cut -c1-6) = lccray ]]
       then
          eval  $submcom  $job_on_remhost
       elif [[  $local_host = lctit  ||  $localhost = lcxe6  ||  $localhost = lck  || $localhost = lckordi ||  $localhost = lcyon || $localhost = lcsb  ||  $localhost = lckyu* ]]
       then
          chmod  u+x  $job_on_remhost
          eval  $submcom  $job_on_remhost
       elif [[ $local_host = nech ]]
       then
          if [[ $queue = default ]]
          then
             eval  $submcom  $job_on_remhost
          else
             eval  $submcom  -q $queue  $job_on_remhost
          fi
       elif [[ $local_host = lcbullhh ]]
       then
          if [[ $queue = default ]]
          then
             eval  $submcom  $job_on_remhost
          fi
       else
          qsub  $job_on_remhost
       fi

          # JOBFILE MUST NOT BE DELETED ON lctit/ibmku/lcflow. THIS WILL BE DONE
          # AT THE END OF THE JOB
       if [[ $local_host != lctit  &&  $local_host != ibmku  &&  $local_host != lcflow ]]
       then
          rm  $job_on_remhost
       fi
       cd  -  > /dev/null
    fi
 fi

    # FINAL ACTIONS
 if [[ $no_submit = false ]]
 then
    rm  -f $job_to_send
 fi
 [[ $verify = true ]]  &&  printf "\n\n *** SUBJOB finished \n\n"
