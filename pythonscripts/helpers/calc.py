from math import sqrt
import scipy.constants as cs
import numpy as np
from scipy.stats import maxwell

# 1.8075 1.8075 0.0 0.0 1.8075 1.8075 1.8075 0.0 1.8075
# print(40/3.615)

# Calc sputtered atom energy
def v_to_E(vel):
    M = 63.55  # g/mol Cu
    m = M / cs.Avogadro  # g
    v = vel * cs.angstrom / 1e-12 # Å/ps -> m/s
    E_k = 1/2 * (m/1000) * v**2 # kgm^2/s^2
    # print(E_k)
    E_k_eV = E_k / cs.electron_volt  # eV
    return E_k_eV


def E_to_v(E):
    # E in eV
    # E = 1/2 m v**2 -> v = sqrt(2E/m)
    M = 63.55  # g/mol Cu
    m = M / cs.Avogadro  # g
    E_J = E * cs.electron_volt  # E_k in joules
    v = sqrt(2*E/(m/1000) / cs.angstrom**2 * 1e-12**2)  # m/s
    return v 


# print(v_to_E(245))
# print(E_to_v(32.9))
# for i in range(10):
#     print(0.1*i, end=" ")

# print(sqrt(cs.electron_volt/cs.atomic_mass)/100)
# print(53416-53312-114)

# temp = 300
# a = np.sqrt(cs.Boltzmann * temp / (63.55*cs.atomic_mass)) 
# # print(cs.Boltzmann)
# # print(cs.atomic_mass)
# rv = maxwell.rvs(scale=a)
# print(rv)


print(1/113129229.43*0.01)
print(2**31-1)

import numpy as np

print(np.random.SeedSequence().generate_state(1))
rng = np.random.default_rng(seed=214079644654894187569045748388860788528)

for i in range(10):
    print(rng.integers(1000000, 2**31-1))