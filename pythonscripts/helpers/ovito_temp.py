from ovito.io import import_file
import scipy.constants as cs
import numpy as np

# Compute the temperature of the edge

# Init pipeline from dump file
print('Importing dump file...')
pipeline = import_file("../dumpfiles/olds/hotcrater.dump")
print('Dump file imported')

n_frames = pipeline.num_frames

# TODO: Add the temperatures into a list and add possibility to plot them
for frame in range(n_frames):
    data = pipeline.compute(frame)  # Compute values at the current frame

    atom_positions = np.array(data.particles.positions)
    atom_velocities = np.array(data.particles['Velocity Magnitude'])

    edge_indices = np.where((atom_positions[:, 0] > 340) & (atom_positions[:, 2] < 120))[0]  # Get indices of atoms with x > 340 and z < 120
    edge_velocities = atom_velocities[edge_indices]  # Use indices to get corresponding velocities

    # Calculate the temperature of the edge atoms
    # Temperature is proportional to the average kinetic energy: KE = 0.5 * m * v^2

    boltzmann_constant = cs.Boltzmann  # Boltzmann constant in J/K
    edge_velocities = edge_velocities * 100  # Convert velocities from Å/ps to m/s
    atom_mass = 63.55 * cs.atomic_mass  # Atomic mass of Cu in kg

    # Compute kinetic energy for edge atoms
    kinetic_energy = 0.5 * atom_mass * np.sum(edge_velocities**2)

    # Calculate temperature
    temperature = (2/3) * (kinetic_energy / (len(edge_velocities) * boltzmann_constant))

    print(f"Temperature of edge atoms: {temperature} K")
