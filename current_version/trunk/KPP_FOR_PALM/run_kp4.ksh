#!/bin/ksh


echo "\$1 = " $1 

if [[ $1 == "clean" ]]
then
  cd kpp
  make distclean
  cd ..

  cd kp4
  make distclean                 
  cd ..
  exit
fi



export PATH=$PATH:`pwd`/kp4/bin

# build kpp

cd kpp
make
cd ..

# build kp4.exe

cd kp4
make install                    
cd ..

# run kp4 with default Setup
 
echo $PATH
kp4.ksh -d `pwd`/kp4/def_smog
