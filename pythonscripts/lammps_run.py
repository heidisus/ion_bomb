from lammps import lammps, LMP_VAR_ATOM
from mpi4py import MPI

# Run parallel using mpirun -np N python3 lammps_run.py

lmp = lammps()
# lmp.file("input_scripts/in.10nm_flux_1358")
# lmp.file("input_scripts/in.28x28x17_relaxation")
# lmp.file("input_scripts/in.heat_center")
lmp.file("input_scripts/in.rescale_py")

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
print("Proc %d out of %d procs has" % (me,nprocs),lmp)

MPI.Finalize()
