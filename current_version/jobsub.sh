#!/bin/bash



touch trunk/SOURCE/* 
mbuild -h lckeal -K parallel -v -u
mbuild -h lckeal -K parallel -v

mrun -B -v -d chem -DKPP_CHEM -h lckeal -K parallel -X 4 -T 4 -r "d3# 3d# pr#" -b -q magny -t 3600 -m 100 

# use option -B to keep the the simulation folder in the current temp folder
