import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read in data of n of atoms that have left the simulation box
data = pd.read_csv('data/yielddata1358.csv')

print(data)  # Print given data

N = data['N'][0]
E = data['Energy'][0]

# Remove incident atoms - Use only if the datafile has incident atoms included in the yield 
# data['Yield'] = data['Yield'].sub(N)

# Normalize the results
data['Yield'] = data['Yield'].truediv(N)

# Set negative values to 0 (values where less atoms left the box than were added)
data['Yield'] = data['Yield'].where(data['Yield'] > 0, 0)

print(data)  # Print modified data

# Plot the yield
data.plot('Angle', 'Yield')
plt.title(f'Sputtering yield {E} eV')
plt.xlabel('Angle (deg)')
plt.ylabel('Yield')
plt.show()