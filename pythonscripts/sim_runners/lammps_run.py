from lammps import lammps, LMP_VAR_ATOM
from mpi4py import MPI

"""Script to test running LAMMPS from python"""

# Run parallel using mpirun -np N python3 pythonscripts/sim_runners/lammps_run.py

lmp = lammps()
# lmp.file("input_scripts/in.10nm_flux_1358")
lmp.file("lammps/relaxation_scripts/in.relaxation")
# lmp.file("input_scripts/in.heat_center")
# lmp.file("lammps/input_scripts/in.rescale_py")

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
print("Proc %d out of %d procs has" % (me,nprocs),lmp)

MPI.Finalize()
