"""
Create a simulation in which atoms are sputtered onto a surface at random positions and times based on a given particle flux.
Usage:
    mpirun -np nprocs run_crater.py temperature radius nsteps timestep flux energy filename
"""

from lammps import lammps
import numpy as np
from random import randint, seed
from mpi4py import MPI
import sys

# Create the lammps instance based on an input script
lmp = lammps()

# Set initial run parameters
temperature = float(sys.argv[1])    # temperature of the of the surface in K
radius = int(sys.argv[2])           # radius of the showered area in Å
nsteps = int(sys.argv[3])           # number of timesteps in showering simulation
timestep = float(sys.argv[4])       # timestep in ps
flux = float(sys.argv[5])           # particle flux in particles/ps/nm^2
energy = float(sys.argv[6])	    # incidence energy of the sputtered atoms in eV
input_file = sys.argv[7]            # name of the input file

lmp.command(f'variable r equal {radius}')
lmp.command(f'variable T equal {temperature}')
lmp.command(f'variable flux equal {flux}')

area = np.pi*(radius/10)**2                  # area of the sputtered area in nm^2
runtime = nsteps*timestep                    # total runtime in ps
sputter_time = runtime*0.75                  # time in ps during which the sputtering occurs
flux_area = flux*area                        # flux per nm^2 in particles/ps

lmp.file("in.crater")  # Input script must not have a run command, as the run is controlled through this python script
rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

# Set the seed for randomisation
rng_seed = 6547
rng = np.random.default_rng(seed=rng_seed)

# Sample the arrival times of the sputtered atoms as a Poisson process 
current_time = 0
insert_times = []

while current_time < sputter_time:
    arrival_time = rng.exponential(scale=1/flux_area)  # ps
    current_time += arrival_time * 1000 # fs
    insert_times.append(current_time)

n_sputtered = len(insert_times)  # number of sputtered atoms

v = np.sqrt(energy*2.0/63.55)*98.2269  # Total velocity of the sputtered atoms based on the incidence energy, 98 factor from scipy.constants
vx = 0   # x-component of the velocity
vy = 0   # y-component of the velocity
vz = -v  # z-component of the velocity

# Set the region in which particles are created and the timestep in lammps
lmp.command(f'timestep {timestep}')

comm = MPI.COMM_WORLD

# Initialise counters
current_step = 0
particles_inserted = 0

while particles_inserted != n_sputtered:
    # Run until we reach the next insertion
    runsteps = insert_steps[particles_inserted] - current_step
    lmp.command(f'run {int(runsteps):d}')
    current_step += runsteps

    # Perform the insertion
    natoms = lmp.get_natoms()
    lmp.command(f'reset_atoms id')  # Ensures that the correct atom is added to the group
    lmp.command(f'create_atoms 1 random 1 {rng_seed*current_step} sputter overlap 1')  # Create a new atom to be sputtered

    # Check if any atoms leave the box within the same timestep
    new_natoms = lmp.get_natoms()
    new_atom_id = new_natoms
    
    if new_natoms == natoms:
        # An atom has just left - the new atom's atom_id is natoms
        lmp.command(f'group newatom id {natoms}')
    else:
        # No atoms have left - the new atom's atom_id is new_natoms
        lmp.command(f'group newatom id {new_atom_id}')

    lmp.command(f'velocity newatom set {vx} {vy} {vz}') # Give the atom a downwards speed
    lmp.command(f'group newatom delete')                # Delete the group
    particles_inserted += 1
    lmp.command(f'variable inserted_atoms equal {particles_inserted}')


# Once all particles have been added, run until the end
lmp.command(f'run {int(nsteps*0.75-current_step):d}')
current_step = nsteps

# Allow the particles to relax with no sputtering for 25% of the sputtered time
lmp.command(f'run {int(nsteps*0.25):d}')

lmp.command(f'write_data /scratch/djurabek/heilait/data.crater_{temperature}_{flux}')

print("Proc %d out of %d procs has" % (rank,nprocs),lmp)

if rank == 0:
    print(f'Area: {area}')
    print(f'Inserted particles: {particles_inserted}')
    print(f'Deposited energy per nm^2: {particles_inserted*energy/area}')

MPI.Finalize()
