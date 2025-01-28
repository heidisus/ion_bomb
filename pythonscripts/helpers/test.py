import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import pandas as pd
import numpy as np
import ctypes
from scipy.stats import maxwell, lognorm
from scipy import constants as con
from random import randint
import cv2
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

# print(25*257/(-5.317140558637236-6.040140558637236)**2)
# print(100*108/5)

# nsteps = 1600
# print(f'run {int(nsteps*0.1):.0f}')
# def get_vel(arr):
#     return [rank,rank,rank]

# rank = MPI.COMM_WORLD.Get_rank()
# nprocs = MPI.COMM_WORLD.Get_size()


# coord_vectors = np.array([[1,2,3], [2,3,4], [3,4,5], [4,5,6],
#                           [1,2,3], [2,3,4], [3,4,5], [4,5,6],
#                           [1,2,3], [2,3,4], [3,4,5], [4,5,6],
#   #                        [1,2,3], [2,3,4], [3,4,5], [4,5,6],
#                           [1,2,3], [2,3,4], [3,4,5], [4,5,6],])  # Reshape array to 2D array [[x1, y1, z1], [x2, y2, z3] ...]
# vels_vectors = np.array([[1,2,3], [2,3,4], [3,4,5], [4,5,6],
#                          [1,2,3], [2,3,4], [3,4,5], [4,5,6],
#    #                      [1,2,3], [2,3,4], [3,4,5], [4,5,6],
#                          [1,2,3], [2,3,4], [3,4,5], [4,5,6],
#                          [1,2,3], [2,3,4], [3,4,5], [4,5,6],])       # Reshape array to 2D array [[vx1, vy1, vz1], [vx2, vy2, vz3] ...]

# atoms = len(coord_vectors)
# atoms_per_proc = int(atoms/nprocs)  # Number of atoms per process

# # Determine the chunk for this process, the remainder is added to the last process
# start = rank*atoms_per_proc
# end = (rank+1)*atoms_per_proc
# partial_coords = coord_vectors[start:end] if rank != nprocs-1 else coord_vectors[start:] 
# partial_vels = vels_vectors[start:end] if rank != nprocs-1 else vels_vectors[start:]        

# # Create a mask based on atom height, upper boundary is required to make sure sputtered atoms still move
# mask = np.logical_and(partial_coords[:, 2] > 3, partial_coords[:, 2] < 6)

# # Combine the vectors so get_vel can be applied along the axis
# combined_vectors = np.concatenate((partial_coords, partial_vels), axis=1)

# # Replace elements in the velocity vectors based on the coordinates
# partial_vels[mask] = np.apply_along_axis(get_vel, 1, combined_vectors[mask])

# # Gather all partial velocities back to the root process
# all_vels = MPI.COMM_WORLD.gather(partial_vels, root=0)

# print(rank)

# if rank == 0:
#     print(np.hstack(np.concatenate(all_vels)))


# df1 = pd.read_csv('copper/1358Edrift5fs.tsv', delim_whitespace=True)
# df2 = pd.read_csv('copper/1358Edrift10fs.tsv', delim_whitespace=True)
# df3 = pd.read_csv('copper/1358Edrift1fs.tsv', delim_whitespace=True)
# df4 = pd.read_csv('copper/1358Edrift2fs.tsv', delim_whitespace=True)
# df5 = pd.read_csv('copper/1358Edrift3fs.tsv', delim_whitespace=True)

# plt.plot(df1['Step'], df1['TotEng'], label='5 fs')
# plt.plot(df2['Step'], df2['TotEng'], label='10 fs')
# plt.plot(df3['Step'], df3['TotEng'], label='1 fs')
# plt.plot(df4['Step'], df4['TotEng'], label='2 fs')
# plt.plot(df5['Step'], df5['TotEng'], label='3 fs')

# plt.legend()
# plt.show()

import string

# a = 'Multi-scale model of vacuum arcing'
# print(a.title())

# data = pd.read_csv('data/csv_files/graphdata.csv')
# print(data.median())

# x_val = set()
# with open('copper/data/csv_files/1um.csv', 'r') as file:
#     data = file.readlines()
#     for line in data:
#         x, y = line.split(',')
#         if x in x_val:
#             print(x)
#         else:
#             x_val.add(x)


# data = pd.read_csv('copper/data/csv_files/cath_temp_30A_1um2.csv')
# # data = pd.read_csv('copper/data/csv_files/1um.csv')
# data['r'] = data['r']  # Convert the radius from m to Å

# plt.scatter(data['r'], data['T (K)'])
# # plt.loglog()
# plt.show()

# x = np.linspace(0, 100000, 50000)
# interpolator = CubicSpline(data['r'], data['T (K)'])
# print(interpolator(18848))
# plt.plot(x, interpolator(x))
# # plt.show()

with open('2908/crater_600_highflux.dump', 'r') as infile:
    # Open a new file to write the modified data
    with open('2908/fixed_dump.dump', 'w') as outfile:
        for line in infile:
            # Split the line into columns
            columns = line.split()
            
            # Check if the line is data (skip header lines)
            if len(columns) > 0 and columns[0].isdigit():
                # If the line has only 5 columns (id, type, xs, ys, zs), append zero velocities
                if len(columns) == 5:
                    line = line.strip() + ' 0.0 0.0 0.0\n'
                
            # Write the processed line to the new file
            outfile.write(line)