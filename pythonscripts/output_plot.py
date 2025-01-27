import pandas as pd
import matplotlib.pyplot as plt

filename = 'copper/3008/hotcrater_1600_350ps.tsv'

# Determine the # of lines before output data, then read the output
rows_skipped = 0

# with open(filename, 'r') as input:
#     for line in input:
#         if line.split()[0] == 'Step':  # Break when the header row has been found
#             break
#         rows_skipped += 1
    
df = pd.read_csv(filename, delim_whitespace=True, skiprows=rows_skipped, engine='python', skipfooter=4)  # skipfooter=32,
print(df.head())
df['norm_atoms'] = df['Atoms']-df['Atoms'][0]
# Plot the output values
df.plot(x='Time', y='Temp', label='Atoms')

plt.draw()
plt.title('Change in number of atoms over time')
plt.xlabel('Time (ps)')
plt.ylabel('Change in atoms')
plt.show()
