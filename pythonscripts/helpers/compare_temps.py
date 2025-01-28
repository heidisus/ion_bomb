import pandas as pd
import matplotlib.pyplot as plt

# file1 = 'sflux_0_0_output.txt'
# file2 = 'sflux_low_0_0_output.txt'
# file3 = 'noflux_0_0_output.txt'

file1 = 'noflux_0_output.txt'
file2 = 'noflux_1_output.txt'
file3 = 'noflux_2_output.txt'

# Determine the # of lines before output data, then read the output
rows_skipped = 0

with open(file1, 'r') as input:
    for line in input:
        if line.split()[0] == 'Step':  # Break when the header row has been found
            break
        rows_skipped += 1
    
df1 = pd.read_csv(file1, delim_whitespace=True, skiprows=rows_skipped, skipfooter=32, engine='python')
df2 = pd.read_csv(file2, delim_whitespace=True, skiprows=rows_skipped, skipfooter=32, engine='python')
df3 = pd.read_csv(file3, delim_whitespace=True, skiprows=rows_skipped, skipfooter=32, engine='python')

# Plot the output values
df1.plot('Step', ['Temp', 'c_tsurf'])
df2.plot('Step', ['Temp', 'c_tsurf'])
df3.plot('Step', ['Temp', 'c_tsurf'])

print('No flux')
print(df1['Temp'].mean(), df1['c_tsurf'].mean())
print(df2['Temp'].mean(), df2['c_tsurf'].mean())
print(df3['Temp'].mean(), df3['c_tsurf'].mean())

plt.draw()
# plt.show()