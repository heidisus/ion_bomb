#!/bin/bash

T=$1  # Target T of surface

mpirun -np 4 lmp_mpi -in in.Cu_heat_surface > Cu_heat_surf/output/output_$T.tsv  # npt
mpirun -np 4 lmp_mpi -in in.Cu_reheat_surface > Cu_heat_surf/output/output_2_${T}.tsv  # nvt

