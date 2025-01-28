import pandas as pd
import numpy as np

# create points at which T is evaluated
r_arr = np.linspace(0, 200, 40)
# print(x_arr, y_arr)
df = pd.DataFrame({'r': r_arr, 'T (K)': [np.exp(i/2)/2e5 for i in range(40, 0, -1)]})

print(df)
df.to_csv('data/csv_files/T_data_file_test.csv', index=False)
