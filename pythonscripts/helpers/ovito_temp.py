"""Script to compute the temperature of the edge of a simulation box from a dump file using OVITO.
Parameters:
    filename: name of the dump file"
    dump_freq: dump frequency in ps
Usage:
    python3 ovito_temp.py filename dump_freq
"""

from ovito.io import import_file
import scipy.constants as cs
import numpy as np
import pandas as pd
import sys

# Compute the temperature of the edge
# Input file and its dump frequency can be given from the command line or by changing the default values:
filename = "hotcrater.dump"
dump_freq = 5  # Dump frequency in ps

try:
    filename = sys.argv[1]  # Name of the dump file
    dump_freq = int(sys.argv[2])  # Dump frequency in ps

except IndexError:
    print("No filename provided. Using default dump file.")


# Init pipeline from dump file
print('Importing dump file...')
pipeline = import_file(f"../dumpfiles/{filename}")
print('Dump file imported')

n_frames = pipeline.num_frames
temperatures = []  # Initialize list of edge temperatures 
edge_coordinates = [340, 350, 120]  # List of edge coordinates (x, y, z) (y-currently not used)

print(f'Computing temperature for {n_frames} frames...')

for frame in range(n_frames):
    data = pipeline.compute(frame)  # Compute values at the current frame

    # Get positions and velocities of all atoms
    atom_positions = np.array(data.particles.positions)
    atom_velocities = np.array(data.particles['Velocity Magnitude'])

    # Get indices of atoms in the edge region
    edge_indices = np.where((atom_positions[:, 0] > edge_coordinates[0]) & (atom_positions[:, 2] < edge_coordinates[2]))[0]
    edge_velocities = atom_velocities[edge_indices]  # Use indices to get corresponding velocities

    # Calculate the temperature of the edge atoms from average kinetic energy: ke = 1/2 * m*v^2
    boltzmann_constant = cs.Boltzmann  # Boltzmann constant in J/K
    edge_velocities = edge_velocities * 100  # Convert velocities from Å/ps to m/s
    atom_mass = 63.55 * cs.atomic_mass  # Atomic mass of Cu in kg

    # Compute kinetic energy for edge atoms
    kinetic_energy = 0.5 * atom_mass * np.sum(edge_velocities**2)

    # Calculate temperature
    temperature = (2/3) * (kinetic_energy / (len(edge_velocities) * boltzmann_constant))
    time = frame * dump_freq  # Frame number multiplied by the dump frequency to get time in ps
    # print(f"Temperature of edge atoms: {temperature} K")
    temperatures.append((time, temperature))

# Save the temperatures to a csv file
print('Saving temperatures to a csv file...')
df = pd.DataFrame(temperatures, columns=['Time (ps)', 'Temperature (K)'])
output_file = f"tables/temperatures_{filename[:-5]}.csv"
df.to_csv(output_file, index=False)