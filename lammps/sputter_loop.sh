#!/bin/bash

# A=$1  # deg
E=200  # eV
N=15  # 0-indexed, N=4 will produce 5 simulations
T=1358

# Loop over angles, simulate each angle 50 times
# for ((A=0; A<90; A+=5))
# do
#     for ((i=0; i<$N; i++))    
#     do
#         echo "Starting run $i, A: $A, E: $E"
#         mpirun -np 4 lmp_mpi -in input_scripts/in.sputter_loop > Cu_sputter_loops/output/output_sputter$i.tsv -var E $E -var A $A -var N $i
#     done

#     # Calculate number of atoms that left the box
#     echo "Running python script"
#     python3 pythonscripts/calc_atoms.py $A $E $N

# done

for ((A=0; A<80; A+=5))
do
    for ((i=0; i<$N; i++))    
    do
        echo "Starting run $i, A: $A, E: $E, T: $T"
        mpirun -np 4 lmp_mpi -in input_scripts/in.sputter_melted > Cu_heat_surf/output/output_sputter$i.tsv -var E $E -var A $A -var N $i -var T $T
    done

    # Calculate number of atoms that left the box
    echo "Running python script"
    python3 pythonscripts/calc_atoms.py $A $E $N $T

done
