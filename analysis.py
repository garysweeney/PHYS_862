import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.optimize import curve_fit

# Fancy plots
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def read_data(element):

    data = np.genfromtxt("{}_data.dat".format(element))

    energy = np.array(data[:,1])
    flux = np.array(data[:,3])
    stat_uncert = np.array(data[:,4])

    return energy, flux, stat_uncert

boron_energy, boron_flux, boron_uncert = read_data("boron")
carbon_energy, carbon_flux, carbon_uncert = read_data("carbon")

# Compute boron-to-carbon ratio
ratios = boron_flux / carbon_flux

# Compute uncertainty
# dZ/Z = sqrt((dA/A)^2 + (dB/B)^2)
uncert = np.sqrt((boron_uncert/boron_flux)**2 + (carbon_uncert/carbon_flux)**2)
uncert = ratios * uncert

# Define power law function
def power_law(x, a, b):
    return a * np.power(x, b)

# Perform curve fitting
params, covariance = curve_fit(power_law, boron_energy, ratios, maxfev=2000)

# Extracting coefficients
a = params[0]
b = params[1]

# Plot ratio and fit
plt.figure()
plt.xlabel("Energy per nucleon (GeV/n)")
plt.ylabel("Boron-to-Carbon Ratio")
plt.scatter(boron_energy, ratios, color='black')
plt.errorbar(boron_energy, ratios, yerr=uncert, fmt='o', color='black', capsize=3)
plt.plot(boron_energy, power_law(boron_energy,a,b), color='red')
plt.xscale("log")
plt.yscale("log")
plt.show()

"""
plt.figure()
plt.xlabel("Energy per nucleon (GeV/n)")
plt.ylabel(r"Flux (GeV/n m$^{-2}$ sr$^{-1}$ s$^{-1}$)")
plt.scatter(boron_energy, boron_flux, color='blue')
plt.errorbar(boron_energy, boron_flux, yerr=boron_uncert, fmt='o', color='blue', capsize=3)
plt.scatter(carbon_energy, carbon_flux, color='red')
plt.errorbar(carbon_energy, carbon_flux, yerr=carbon_uncert, fmt='o', color='red', capsize=3)
plt.xscale("log")
plt.yscale("log")
plt.show()
"""