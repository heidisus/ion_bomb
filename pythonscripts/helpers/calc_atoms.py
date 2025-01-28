import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import sys
from ovito.io import import_file


# Take in angle and E from the command line if given, set the except values to desired ones when not running through command line
try:
    angle = int(sys.argv[1])
    energy = int(sys.argv[2])
    repeats = int(sys.argv[3])
    temp = int(sys.argv[4])
except:
    angle = 75
    energy = 0
    repeats = 5
    temp = 'test'

cutoff_coordinate = 25

# Loop over output files
n_left = []
for i in range(repeats):
    filename = f"Cu_heat_surf/dump/sputter_{angle}_{i}.dump"
    pipeline = import_file(filename)  # Import the dumpfile into the pipeline
    n_frames = pipeline.source.num_frames  # Find the number of frames in the pipeline

    data0 = pipeline.compute(0)  # Compute values at the first frame
    n_at_0 = len(data0.particles.positions)  # Find number of atoms at the first frame

    data = pipeline.compute(n_frames-1)  # Compute values at the last frame
    n_at_end = len(data.particles.positions)  # Find number of atoms at the last frame

    positions = np.array(data.particles.positions)  # np.array of all particle positions in the order of the dump file
    
    # Mask the positions to return only those with z>20
    above_surf = np.where(positions[:, 2] > cutoff_coordinate, positions[:, 2], 0)  # Find the z-column for values above the cutoff_coordinate, set to 0 else
    n_above_surf = np.count_nonzero(above_surf)  # Count the number of atoms where z>cutoff_coordinate

    # Count the number of sputtered atoms, does not include the incident atoms
    atoms_left = n_at_0 - n_at_end + n_above_surf
    n_left.append(atoms_left)  


# Calculate sum of atoms that left the box
sum_left = sum(n_left)
print(n_left)
print(sum_left)

# Collect data to write into a file
cols = [angle, sum_left, energy, repeats]
sdata = [str(i) for i in cols]
empty = False

# Check if the header row has been added
with open(f'data/yielddata{temp}.csv', 'r') as f:
    if len(f.readlines()) == 0:
        empty = True

# Write the data        
with open(f'data/yielddata{temp}.csv', 'a') as f:
    if empty:
        f.write('Angle,Yield,Energy,N\n')
    datastr = ','.join(sdata)
    f.write(datastr + '\n')
