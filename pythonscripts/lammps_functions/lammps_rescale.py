"""
Module containing functions to rescale the velocities of atoms in a LAMMPS simulation based on interpolated data from a .csv file. 

Functions:
    calc_dist(x, y): Calculate 2D distance between (0,0) and a given point on the xy-plane.
    get_vel(data_vector): Return randomly rescaled velocity based on interpolated temperature based on the Maxwell-Boltzmann distribution and the x,y coordinates.
    rescale(lmp_ptr): Rescale the temperatures of particles using interpolated data. Requires continuous atoms ids to exist.

Usage in LAMMPS:
    Initialize the function in the LAMMPS input script. The function can be called as follows: 
    python rescale input 1 SELF format p file lammps_rescale.py
    fix pf all python/invoke <nsteps> end_of_step rescale

Note:
    The fix syntax needs a group to apply the fix to, however, since gather_atoms() and scatter_atoms() are used for rescaling, the fix will be applied to all 
    atoms between the lower_bound and upper_bound regardless of the given group. 
    The lower_bound and upper_bound must be defined as variables in the LAMMPS input script with the names 'heat_bot' and 'heat_top' respectively.
    or the values must be hardcoded in the rescale function.
"""

from lammps import lammps
import numpy as np
from scipy.interpolate import CubicSpline
import pandas as pd
from scipy.stats import maxwell
from scipy import constants as con
from mpi4py import MPI

# This block gets run when Lammps first reads through the file before the run begins
# Initiate the interpolator based on the datafile on all processes
data = pd.read_csv('data/csv_files/1um.csv') 
data['r'] = data['r']*0.007  # Rescale the radius to match the simulation box
interpolator = CubicSpline(data['r'], data['T (K)'])
print("interpolator initiated")

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

def calc_dist(x, y):
    """
    Calculate 2D distance between (0,0) and a given point on the xy-plane. 
    params: x, y
    """
    return np.sqrt(x**2 + y**2)


def get_vel(data_vector):
    """
    Return randomly rescaled velocity based on interpolated temperature based on the Maxwell-Boltzmann distribution and the x,y coordinates.
    params: 
        data_vector: vector of coordinates and velocities [c1, c2, c3, v1, v2, v3]
    returns: 
        rescaled velocity vector as an array [vx, vy, vz]
    """
    coordinates = data_vector[:3]
    velocities = data_vector[3:]

    x, y = coordinates[0], coordinates[1]
    dist = calc_dist(x, y)
    temp = interpolator(dist)   # Find the temperature at the given dist
    a = np.sqrt(con.Boltzmann * temp / (63.55*con.atomic_mass))  # [m/s]
    rv = maxwell.rvs(scale=a)   # Random variate from the distribution [m/s]
    v = rv/100                  # Conversion into [Å/ps]
    
    # Find the current velocity and scale the new velocity in the same direction
    v0 = np.sqrt(velocities[0]**2 + velocities[1]**2 + velocities[2]**2)
    if v0 == 0:
        return velocities

    vx = velocities[0]/v0 * v
    vy = velocities[1]/v0 * v
    vz = velocities[2]/v0 * v

    # Return the rescaled velocity vector
    return [vx, vy, vz]


def rescale(lmp_ptr):
    """
    Rescale the temperatures of particles using interpolated data. Requires continuous atoms ids to exist.
    params:
        lmp_ptr: Lammps pointer
    returns:
        None
    """
    # Currently scales the velocities of all atoms in a certain height range.
    # Depth profiling could be added, but would require a data set with a depth profile, and a different interpolator to suit it
    # If the function is applied to all atoms, it still needs a mask to avoid scaling the sputtered atoms. Where to set the height when the surf heats?

    lmp = lammps(ptr=lmp_ptr)             # Pointer to Lammps
    natoms = lmp.get_natoms()             # Number of atoms in the simulation
    atoms_per_proc = int(natoms/nprocs)   # Number of atoms per process

    try:
        # Extract the upper and lower bounds for the height range from Lammps
        upper_bound = lmp.extract_variable('heat_top', 'all', 0)  # Upper boundary for the height range
        lower_bound = lmp.extract_variable('heat_bot', 'all', 0)  # Lower boundary for the height range
    except:
        # Values to use if the variables are not defined
        upper_bound = 32                   
        lower_bound = 0

    coords = lmp.gather_atoms('x', 1, 3)  # Retrieve array of all atom positions
    vels = lmp.gather_atoms('v', 1, 3)    # Retrieve array of all atom velocities

    numpy_coords = np.ctypeslib.as_array(coords)                             # Form coordinates into a numpy array
    coord_vectors = np.reshape(numpy_coords, (int(len(numpy_coords)/3), 3))  # Reshape array to 2D array [[x1, y1, z1], [x2, y2, z3] ...]
    numpy_vels = np.ctypeslib.as_array(vels)                                 # Form velocities into a numpy array
    vels_vectors = np.reshape(numpy_vels, (int(len(numpy_vels)/3), 3))       # Reshape array to 2D array [[vx1, vy1, vz1], [vx2, vy2, vz3] ...]

    # Determine the chunk for this process, the remainder is added to the last process
    start = rank*atoms_per_proc
    end = (rank+1)*atoms_per_proc
    partial_coords = coord_vectors[start:end] if rank != nprocs-1 else coord_vectors[start:] 
    partial_vels = vels_vectors[start:end] if rank != nprocs-1 else vels_vectors[start:]        

    # Create a mask based on atom height, upper boundary is required to make sure sputtered atoms still move
    mask = np.logical_and(partial_coords[:, 2] > lower_bound, partial_coords[:, 2] < upper_bound)

    # Combine the coords and vels so get_vel() can be applied along the axis while getting both vectors
    combined_vectors = np.concatenate((partial_coords, partial_vels), axis=1)

    # Replace elements in the velocity vectors based on the coordinates
    try:
        partial_vels[mask] = np.apply_along_axis(get_vel, 1, combined_vectors[mask])
        # print(rank, "success")
    except ValueError:
        # A ValueError is raised when the mask is empty, meaning no atoms are in the height range for this process
        # print(rank, "failed")
        pass

    # Gather all partial velocities back to the root process
    all_vels = comm.gather(partial_vels, root=0)

    # Reshape velocities back into the original shape, change into ctypes and scatter back to all other processes
    if rank == 0:
        vels = np.hstack(np.concatenate(all_vels))
    else:
        vels = None

    vels = comm.bcast(vels, root=0)
    ctype_scaled_vels = np.ctypeslib.as_ctypes(vels)

    # Scatter the rescaled velocities back to Lammps
    lmp.scatter_atoms('v', 1, 3, ctype_scaled_vels)
