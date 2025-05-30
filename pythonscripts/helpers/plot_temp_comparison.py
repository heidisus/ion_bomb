import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file1 = 'tables/hotcrater.tsv'
file2 = 'tables/temperatures_hotcrater.csv'

title = 'Heated center 1600K - 300K edge'
ylabel = 'T (K)'
xlabel = 'Time (ps)'

# Determine the # of lines after output data, then read the output
rows_skipped = 0

info = []
# Find number of lines to skip at the end of the file
with open(file1, 'r') as input:
    lines = input.readlines()
    for line_idx in range(len(lines)-1,0, -1):
        try:
            int(lines[line_idx].split()[0])
            break
        except:
            rows_skipped += 1
    info = lines[len(lines)-rows_skipped:]  # Add the data printed at the end of the run to the info list


fig, ax = plt.subplots()

df1 = pd.read_csv(file1, delim_whitespace=True, engine='python', skipfooter=rows_skipped) 

df1.plot(x='Time', y=['Temp'], label=['Entire box'], ax=ax)

df2 = pd.read_csv(file2)

df2.plot(x='Time (ps)', y='Temperature (K)', label='Edge (1 nm)', ax=ax)

plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.grid()

plt.legend()
plt.show()

