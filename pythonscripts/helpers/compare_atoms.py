import pandas as pd
import matplotlib.pyplot as plt

# value = 'Temp'
value = 'Atoms'
dir = 'output'

for i in ['crater', 'droplet', 'wave']:
    hf = pd.read_csv(f'{dir}/{i}.csv', delim_whitespace=True, skipfooter=0, engine='python')
    plt.plot(hf['Time'], hf[value], label=f'{i}')

plt.title(f'{value}')
plt.xlabel('Time (ps)')
plt.ylabel(value)
plt.grid()
plt.legend()
plt.show()

# for i in [300, 1000, 1358, 1680, 1818]:
#     lf = pd.read_csv(f'{dir}/{i}Klf.tsv', delim_whitespace=True, skipfooter=4, engine='python')
#     plt.plot(lf['Time'], lf[value], label=f'{i} K')

# plt.title(f'{value} at 10^24 flux')
# plt.xlabel('Time (ps)')
# plt.ylabel(f'{value}')
# plt.grid()
# plt.legend()
# plt.show()