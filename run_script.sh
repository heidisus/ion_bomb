#!/bin/bash

lmpenv

# Assign arguments to variables
TEMP=300  	                # K
TEMP_HOT=5000               # K
RADIUS=30  	                # Å
RUNTIME=5     	        # ps
FLUX=0.1  	                # part/ps/nm^2
ENERGY=100  	            # eV
FILENAME="lammps/input_scripts/in.hot_crater" # Lammps input file
DUMPFILE="testbighot.dump"
OUTPUT="output/testbighot.txt"

# REMEMBER TO CHECK DUMP FREQUENCY

# Run the simulation with provided arguments
# mpirun -np 4 python3 pythonscripts/sim_runners/run_crater_hot.py "$TEMP" "$RADIUS" "$RUNTIME" "$FLUX" "$ENERGY" "$INCRT" "$INCRD" "$INCRBUFF" "$FLUXDECR" "$FILENAME" "$DUMPFILE" > "$OUTPUT"

# No output file
mpirun -np 4 python3 pythonscripts/sim_runners/run_crater_hot.py "$TEMP" "$TEMP_HOT" "$RADIUS" "$RUNTIME" "$FLUX" "$ENERGY" "$FILENAME" "$DUMPFILE"
