# FULLY WRITTEN BY GPT

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Example: Load your digitized data
data = np.loadtxt("e_distribution_plot_data_2.csv", delimiter=",", skiprows=1)
data = data[data[:, 0].argsort()]  # Sort the data by the x-values
energies = data[:, 0]  # x-values
densities = data[:, 1]  # y-values (PDF)

# Normalize the PDF
pdf = densities / np.trapz(densities, energies)

# Create the CDF
cdf = np.cumsum(pdf * np.diff(energies, prepend=0))

# Create an interpolation function for the inverse CDF
inverse_cdf = interp1d(cdf, energies, bounds_error=False, fill_value=(min(energies), max(energies)))

# Sampling from the distribution
n_samples = 100000
uniform_samples = np.random.uniform(0, 1, n_samples)
sampled_energies = inverse_cdf(uniform_samples)
sampled_energies = sampled_energies[sampled_energies <= 200]

# Plot the results
plt.figure()
plt.plot(energies, pdf, label="Normalized PDF")
plt.hist(sampled_energies, bins=50, density=True, alpha=0.5, label="Sampled Distribution")
plt.legend()
plt.xlabel("Energy")
plt.ylabel("Probability Density")
plt.title("Sampling from a Digitized Distribution")
plt.show()
