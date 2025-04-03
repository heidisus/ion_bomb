"""
Create a simulation in which atoms are sputtered onto a surface at random positions and times based on a given particle flux.
This version of the script uses an adaptive timestep, and increases the sputtering radius over time based on LAMMPS regions.
Params:
    temperature: temperature of the surface in K
    radius: radius of the showered area in Å
    runtime: total runtime in ps
    flux: particle flux in particles/ps/nm^2
    energy: incidence energy of the sputtered atoms in eV
    r_incr: crater increment distance in Å
    filename: name of the input file
    dump_file: name of the dump file
Usage:
    mpirun -np nprocs run_crater_expand.py temperature radius runtime flux energy r_incr filename dump_file
"""

from lammps import lammps
import numpy as np
from mpi4py import MPI
import sys

# Create the lammps instance based on an input script
lmp = lammps(cmdargs=["-sf", "opt"])

# Set initial run parameters
temperature = float(sys.argv[1])        # temperature of the of the surface in K
radius = int(sys.argv[2])               # initial radius of the showered area in Å
runtime = float(sys.argv[3])            # total runtime in ps
flux = float(sys.argv[4])               # particle flux in particles/ps/nm^2
energy = float(sys.argv[5])	            # incidence energy of the sputtered atoms in eV
r_incr = float(sys.argv[6])             # crater radius increase in Å
particle_threshold = int(sys.argv[7])   # number of particles in the region before it is expanded
input_file = sys.argv[8]                # name of the input file
dump_file = sys.argv[9]                 # name of the dump file

lmp.command(f'variable r equal {radius}')
lmp.command(f'variable T equal {temperature}')
lmp.command(f'variable flux equal {flux}')
lmp.command(f'variable dmpf string {dump_file}')

lmp.file(input_file)  # Input script must not have a run command, as the run is controlled through this python script

# Set the height of the sputtering region
box_dims = lmp.extract_box()  # Box dims in ([xlo, ylo, zlo], [xhi, yhi, zhi], xy, yz, xz)
zlo = box_dims[0][2] 
zhi = box_dims[1][2]  
cyl_bot = zhi - 15
cyl_top = zhi - 5

top = zhi - 50 + 1  # Top of the surface + 1 Å to not count the surface atoms
checkh = top + 5
check_xlo = radius + 5
check_xhi = radius + 15

# Create the region for sputtering 
lmp.commands_list([f'region sputter cylinder z 0.0 0.0 {radius} {cyl_bot} {cyl_top} units box',  
                   'group sputter region sputter',
                   f'region check block {check_xlo} {check_xhi} INF INF {top} {checkh} units box',
                   'group check dynamic all region check every 50'])


rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

# Set the seed for randomisation
rng_seed = 214079644654894187569045748385860788528  # Random seed for the numpy generator, generated with SeedSequence.entropy
rng = np.random.default_rng(seed=rng_seed)

area = np.pi*(radius/10)**2  # area of the sputtered area in nm^2
sputter_time = runtime   # time in ps during which the sputtering occurs
flux_area = flux*area    # flux per nm^2 in particles/ps
nsteps = int(runtime*10000)  # this argument is passed to the run command as it requires a parameter, however, the actual number of steps taken is controlled by the halt command and the insertion times

# Sample the arrival times of the sputtered atoms as a Poisson process
current_time = 0
insert_times = []

# TODO: Figure out whether N of ions stays consistent or grows with the sputtered area
while current_time < sputter_time:
    arrival_time = rng.exponential(scale=1/flux_area)
    current_time += arrival_time
    insert_times.append(current_time)

n_sputtered = len(insert_times)  # number of sputtered atoms

v = np.sqrt(energy*2.0/63.55)*98.2269  # Total velocity of the sputtered atoms based on the incidence energy, 98 factor from scipy.constants
vx = 0   # x-component of the velocity
vy = 0   # y-component of the velocity
vz = -v  # z-component of the velocity

comm = MPI.COMM_WORLD

# Initialise counters
current_step = 0
particles_inserted = 0
current_time = 0
next_insertion = insert_times[particles_inserted]

lmp.command('compute rimcount check count/type atom')

while particles_inserted != n_sputtered:
    # Set the insertion time
    next_insertion = insert_times[particles_inserted]
    lmp.command(f'fix myhalt all halt 1 v_timee >= {next_insertion} error continue')

    # Run until we reach the next insertion, the simulation will run until the next insertion point is reached
    lmp.commands_list([f'run {nsteps} post no',
                       f'unfix myhalt'])

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

    particles_inserted += 1
    lmp.commands_list([f'velocity newatom set {vx} {vy} {vz}',      # Give the atom a downwards speed
                       f'group newatom delete',                     # Delete the group
                       f'variable inserted_atoms equal {particles_inserted}'])
    
    # Create a new region with the new radius 
    compute_res = lmp.extract_compute('rimcount', 0, 1)
    check_particles = int(compute_res[0])

    if check_particles > particle_threshold:
        radius += r_incr
        lmp.commands_list(['region sputter delete',
                        f'region sputter cylinder z 0.0 0.0 {radius} {cyl_bot} {cyl_top} units box'])

        # Remove the old region and create a new one
        # TODO: Figure out if r_incr is a good value to increase the region by (might be too small) 
        check_xlo += r_incr
        check_xhi += r_incr  
        lmp.commands_list([f'uncompute rimcount',
                           'group check delete',
                           'region check delete',
                           f'region check block {check_xlo} {check_xhi} INF INF {top} {checkh} units box',
                           'group check dynamic all region check every 50',
                           'compute rimcount check count/type atom'])
        
        # TODO: Generate particles for the new ring here if using dynamic flux, and add to the insert_times list?

    

# Once all particles have been added, run until the end
lmp.command(f'fix myhalt all halt 100 v_timee >= {runtime} error continue')
lmp.command(f'run {nsteps}')


# TODO: Determine the area of the sputtered region / the flux
if rank == 0:
    print(f'Area: {area}')
    print(f'Inserted particles: {particles_inserted}')
    print(f'Deposited energy per nm^2: {particles_inserted*energy/area}')
    print(f'Flux: {particles_inserted/sputter_time/area}')


lmp.command(f'write_data {dump_file[:-5]}.data')
lmp.command(f'write_restart {dump_file[:-5]}_restart')

print("Proc %d out of %d procs has" % (rank,nprocs),lmp)

MPI.Finalize()
