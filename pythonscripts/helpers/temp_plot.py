import pandas as pd
import matplotlib.pyplot as plt

filename = 'tables/crater_300K_f26_600ps_100eV.tsv'

# Determine the # of lines before output data, then read the output
rows_skipped = 0

fig, ax = plt.subplots()

df = pd.read_csv(filename, delim_whitespace=True, skiprows=rows_skipped, engine='python', skipfooter=4)  # skipfooter=32,
print(df.head())
df['norm_atoms'] = df['Atoms']-df['Atoms'][0]

df.plot(x='Time', y=['Temp'], ax=ax, label=['Temp'])
# df.plot(x='Time', y='Temp', color='indigo', ax=ax)

filename = 'tables/crater_300K_f25_600ps_100eV.tsv'

# Determine the # of lines before output data, then read the output
rows_skipped = 0

    
# df = pd.read_csv(filename, delim_whitespace=True, skiprows=rows_skipped, engine='python', skipfooter=4)  # skipfooter=32,
# print(df.head())
# df['norm_atoms'] = df['Atoms']-df['Atoms'][0]

# df.plot(x='Time', y=['c_tbox'], ax=ax, label=['tbox'])


# Change legend font
plt.legend()
# plt.title('Temperature')
plt.xlabel('Time (ps)')
plt.ylabel('Temperature (K)')
plt.grid()
plt.show()
