#!/usr/bin/ksh

###########################################################################
# 
#     create_mz_kpp_module 
#
#     create scalar code from .f90 sources created by KPP to be used in MECCA
#
#     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
# 
###########################################################################

set -eu


########################### User SetUp ####################################

export KPP_HOME=/home/modelle/kpp_for_palm/eclipse_neon_ws/kpp_for_palm/kpp
export KPP=$KPP_HOME/bin/kpp

BASE=`pwd`/kp4

########################## End User Setup ################################

WORK=tmp_kp4

# Default

PREFIX=kchem_kpp
OUTDIR=`pwd`/../palm/
DEFDIR=$BASE/def_small_strato
MODE="scalar"
VLEN=1
KEEP="NO"
DE_INDEX="NO"
DE_INDEX_FAST="NO"

export KPP_SOLVER=Rosenbrock

# get Command line option

while  getopts :d:ifkp:o:s:v:w:  c     # get options
do case $c in
      d)   DEFDIR=$OPTARG;;          # directory of definition files

      i)   DE_INDEX="YES";;          # if set, deindexing

      f)   DE_INDEX_FAST="YES";;     # if set, fast deindexing

      k)   KEEP="YES";;              # keep Working directory

      o)   OUTDIR=$OPTARG;;          # Output directory of Geneated Code

      p)   PREFIX=$OPTARG;;          # Name Prefix

      s)   KPP_SOLVER=$OPTARG;;      # Name Prefix

      v)   MODE="vector"
           VLEN=$OPTARG;;            # Set to vector Mode

      w)   WORK=$OPTARG;;            # Working directory

      \?)  print ${0##*/} "unknown option:" $OPTARG
           print "USAGE: ${0##*/} [ -d dir -e -k -o dir -p name -s solver -v length -w dir ] "
           exit 1;;
   esac
done
shift OPTIND-1

DEF_PREFIX=${PREFIX}.kpp

# Create or clean working directory

MY_PWD=`pwd`
mkdir -p $WORK
rm -rf $WORK/*
cd $WORK

# kpp dependend, may be changed

KPP_FILE_LIST="Initialize Integrator LinearAlgebra Jacobian Function Rates Util"


KPP_SUBROUTINE_LIST="Initialize"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST INTEGRATE Fun"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppSolve KppDecomp WLAMCH WLAMCH_ADD"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Jac_SP k_arr "
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Update_RCONST k_3rd "
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST initialize_kpp_ctrl error_output"

# if [[ $MODE = "vector" && $KPP_SOLVER = "ROS2" ]]
# then
#   cp $BASE/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90    # get vector Solver 
# else
# #  KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST FunTemplate JacTemplate Update_SUN "
#   KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
#   if [[ $MODE = "vector" ]]
#   then
#     cp $BASE/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90  # get vector Solver 
#   else
#     KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate JacTemplate Update_SUN"
#   fi
# fi
 if [[ $MODE = "vector" ]]
 then
   # get vector Solver 
   cp $BASE/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90
fi

# Interface ignore list
KPP_INTERFACE_IGNORE="WAXPY WCOPY"

case $KPP_SOLVER in
    ROS2) ;;

    Rosenbrock)   
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
    if [[ $MODE != "vector" ]]
    then
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate JacTemplate Update_SUN"
    fi;;

    rosenbrock_mz)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate JacTemplate Update_SUN";;

    rosenbrock)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate JacTemplate Update_SUN";;

    kpp_lsode)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppLsode DLSODE JAC_CHEM FUN_CHEM"
      KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE JAC_CHEM KppDecomp KppSolve";;

    kpp_radau5)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY FUN_CHEM JAC_CHEM SET2ZERO"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST RADAU5 Update_SUN"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppSolveCmplx KppDecompCmplx";;

    kpp_sdirk)
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST SDIRK JAC_CHEM SET2ZERO FUN_CHEM"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE Set2zero SET2ZERO FUN_CHEM";;

    kpp_seulex)
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST ATMSEULEX"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST SEULEX_ErrorMsg SEULEX_Integrator FUN_CHEM JAC_CHEM SEUL"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE SEULEX_Integrator SDIRK FUN_CHEM SEUL";;

   \?)  print "SORRY ONLY ROSENBROCK METHODS WORK AT THE MOMENT:" $KPP_SOLVER
        exit 1;;
esac
#mz-ak-20070509+

KPP_INCLUDE_LIST="Parameters Global JacobianSP Monitor"

#Get definition Files

cp $DEFDIR/*.eqn         .
cp $DEFDIR/*.spc         .
cp $DEFDIR/${PREFIX}.kpp     .

# Run kpp

$KPP $DEF_PREFIX

# Get templates for C++ program

cp $BASE/templates/module_header* .           # Use fixed Module_header
cp $BASE/templates/initialize_kpp_ctrl_template.f90 .  # CTRL kpp time stepping

# file with subroutine list for c++ program create_mz_kpp_module

for i in $KPP_FILE_LIST
do
  echo ${PREFIX}_${i} >> file_list
done
echo initialize_kpp_ctrl_template >> file_list

# file with subroutine list for c++ program create_mz_kpp_module

for i in $KPP_SUBROUTINE_LIST
do
  echo $i >> subroutine_list
done

# file with include list for c++ program create_mz_kpp_module

for i in $KPP_INCLUDE_LIST
do
  echo ${PREFIX}_${i} >> include_list
done

touch interface_ignore_list
for i in $KPP_INTERFACE_IGNORE
do
  echo $i >> interface_ignore_list
done

$BASE/bin/kp4.exe $PREFIX $MODE $VLEN $DE_INDEX $DE_INDEX_FAST


cp kk_kpp.f90    $OUTDIR/${PREFIX}.f90

echo " "
echo "Write kpp module -- > " $OUTDIR/${PREFIX}.f90

if [[ $KEEP = "NO" ]]
then
  cd  $MY_PWD
  rm -rf $WORK
fi
exit

