"""
Create a simulation in which atoms are sputtered onto a surface at random positions and times over the center based on a given particle flux.
This version of the script uses an adaptive timestep, and sends an initial wave of ions to the surface.
Params: 
    temperature: temperature of the surface in K
    radius: radius of the showered area in Å
    runtime: total runtime in ps
    flux: particle flux in particles/ps/nm^2
    energy: incidence energy of the sputtered atoms in eV
    filename: name of the LAMMPS input file
Usage:
    mpirun -np nprocs run_wave.py temperature radius runtime flux energy filename
"""

from lammps import lammps
import numpy as np
from mpi4py import MPI
import sys

# Create the lammps instance based on an input script
lmp = lammps()

# Set initial run parameters
temperature = float(sys.argv[1])    # temperature of the of the surface in K
radius = int(sys.argv[2])           # radius of the showered area in Å
runtime = float(sys.argv[3])        # total runtime in ps
flux = float(sys.argv[4])           # particle flux in particles/ps/nm^2
energy = float(sys.argv[5])	        # incidence energy of the sputtered atoms in eV
input_file = sys.argv[6]            # name of the input file

lmp.command(f'variable r equal {radius}')
lmp.command(f'variable T equal {temperature}')
lmp.command(f'variable flux equal {flux}')

lmp.file(input_file)  # Input script must not have a run command, as the run is controlled through this python script
rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

# Set the seed for randomisation (ensure all processes create the same random numbers)
rng_seed = 214079644654894187569045748385860788528  # Random seed for the numpy generator, generated with SeedSequence.entropy
rng = np.random.default_rng(seed=rng_seed)

area = np.pi*(radius/10)**2  # area of the sputtered area in nm^2
sputter_time = runtime*0.75  # time in ps during which the sputtering occurs
flux_area = flux*area        # flux per nm^2 in particles/ps
nsteps = int(runtime*10000)  # this argument is passed to the run command as it requires a parameter, however, the actual number of steps taken is controlled by the halt command and the insertion times to use an adaptive timestep but also determine the total simulation time needed, rather than steps. 

v = np.sqrt(energy*2.0/63.55)*98.2269  # Total velocity of the sputtered atoms based on the incidence energy, 98 factor from scipy.constants
vx = 0   # x-component of the velocity
vy = 0   # y-component of the velocity
vz = -v  # z-component of the velocity

# Sample the arrival times of the sputtered atoms as a Poisson process 
current_time = 0
insert_times = []
wave_period = runtime/3
waves = 1

while current_time < sputter_time:
    # Create a list of (insertion time, number of particles to insert) tuples
    arrival_time = rng.exponential(scale=1/flux_area)
    current_time += arrival_time

    if current_time >= wave_period:
        waves += 1
        wave_period = waves * runtime/3
        insert_times.append((current_time, 10))
    else:
        insert_times.append((current_time, 1))

print(insert_times)
n_insertions = len(insert_times)  # number of insertions to be made

# Initialise counters
insertions_made = 0
next_insertion = insert_times[insertions_made][0]

# Function to insert particles
def insert_particles(n):
    global insertions_made

    for i in range(n):
        # Perform the insertion
        natoms = lmp.get_natoms()
        lmp.commands_list([f'reset_atoms id',  # Ensures that the correct atom is added to the group
                           f'create_atoms 1 random 1 {rng.integers(1000000, 2**31-1)} sputter overlap 1'])  # Create a new atom to be sputtered

        # Check if any atoms leave the box within the same timestep
        new_natoms = lmp.get_natoms()
        new_atom_id = new_natoms
        
        if new_natoms == natoms:
            # An atom has just left - the new atom's atom_id is natoms
            lmp.command(f'group newatom id {natoms}')
        else:
            # No atoms have left - the new atom's atom_id is new_natoms
            lmp.command(f'group newatom id {new_atom_id}')

        lmp.commands_list([f'velocity newatom set {vx} {vy} {vz}',      # Give the atom a downwards speed
                           f'group newatom delete',                     # Delete the group
                           f'variable inserted_atoms equal {insertions_made}'])  # TODO: Fix this to reflect the number of atoms, not the number of insertions

    insertions_made += 1


# TODO: Make sure the flux stays consistent
# Run the simulation
while insertions_made < n_insertions:
    # Set the insertion time
    next_insertion = insert_times[insertions_made][0]
    lmp.command(f'fix myhalt all halt 1 v_timee >= {next_insertion} error continue')

    # Run until we reach the next insertion point
    lmp.commands_list([f'run {nsteps}', 
                       f'unfix myhalt'])
    
    insert_particles(insert_times[insertions_made][1])


# Once all particles have been added, run until the end
lmp.command(f'fix myhalt all halt 100 v_timee >= {runtime} error continue')
lmp.command(f'run {nsteps}')

# lmp.command(f'write_data /scratch/djurabek/heilait/data.crater_{int(temperature)}_{flux}')

print("Proc %d out of %d procs has" % (rank,nprocs),lmp)

if rank == 0:
    print(f'Area: {area}')
    print(f'Inserted particles: {insertions_made}')
    print(f'Deposited energy per nm^2: {insertions_made*energy/area}')
    print(f'Flux: {insertions_made/sputter_time/area}')

MPI.Finalize()
