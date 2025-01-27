from __future__ import print_function
import lammps
import ctypes
import traceback
import numpy as np
import pandas as pd


class LAMMPSFix(object):
    def __init__(self, ptr, group_name="all"):
        self.lmp = lammps.lammps(ptr=ptr)
        self.group_name = group_name


class LAMMPSFixMove(LAMMPSFix):
    def __init__(self, ptr, group_name="all"):
        super(LAMMPSFixMove, self).__init__(ptr, group_name)

    def init(self):
        pass

    def initial_integrate(self, vflag):
        pass

    def final_integrate(self):
        pass

    def initial_integrate_respa(self, vflag, ilevel, iloop):
        pass

    def final_integrate_respa(self, ilevel, iloop):
        pass

    def reset_dt(self):
        pass

class VRescale(LAMMPSFixMove):

    def __init__(self, ptr, group_name="all"):
        super().__init__(ptr, group_name)

    def init(self):
        self.data = pd.read_csv('data/cath_temp_250mA_10nm_sorted.csv')

    