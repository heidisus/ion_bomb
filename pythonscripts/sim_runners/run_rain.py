"""
Sputter atoms ontp a surface at random positions and times based on a given particle flux. Creates an "ion shower" effect over the entire surface.
Usage:
    mpirun -np nprocs python3 run_rain.py
"""
from lammps import lammps
import numpy as np
from random import randint, seed
from mpi4py import MPI
import run_rain_cmdargs

# Create the lammps instance based on an input script
lmp = lammps()
lmp.file("input_scripts/in.28x28x17_flux_py")  # Input script must not have a run command, as the run is controlled through this python script
rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

# Set the seed for randomisation
rng_seed = 6547
seed(rng_seed)

# Set initial run parameters
nsteps = 15000      # number of timesteps in full simulation
flux = 1e-1         # particle flux in particles/ps/nm^2
energy = 25         # incidence energy of the sputtered atoms in eV
timestep = 0.002    # timestep in ps

area = lmp.extract_variable('area')  # area of the surface in nm^2

# Calculate other parameters
n_sputtered = int(area * flux * (nsteps/(1/timestep)))  # number of sputtered atoms
insert_steps = [randint(1, nsteps*0.98) for i in range(n_sputtered)]  # generate the random timesteps at which atoms are created
insert_steps.sort()

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
    lmp.command(f'run {runsteps}')
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
lmp.command(f'run {nsteps-current_step}')
current_step = nsteps

# Add flux to the dump file name?
# lmp.command('write_data data/lammps_datafiles/data.halfway')
print("Proc %d out of %d procs has" % (rank,nprocs),lmp)

if rank == 0:
    print(f'Area: {area}')
    print(f'Inserted particles: {particles_inserted}')
    print(f'Deposited energy per nm^2: {particles_inserted*energy/area}')

MPI.Finalize()
