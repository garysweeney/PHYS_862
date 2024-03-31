import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.optimize import curve_fit

# Fancy plots
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# ========================================== PROBLEM 1.1 ==========================================
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
print(a, b)

"""
# Plot ratio and fit
plt.figure()
plt.xlabel("Energy per nuclei (GeV/n)", fontsize=14)
plt.ylabel("Boron-to-Carbon Ratio", fontsize=14)
plt.scatter(boron_energy, ratios, color='black')
plt.errorbar(boron_energy, ratios, yerr=uncert, fmt='o', color='black', capsize=3, label="Data")
plt.plot(boron_energy, power_law(boron_energy,a,b), color='red', label='Fit')
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale("log")
plt.yscale("log")
plt.show()
"""

# ========================================== PROBLEM 1.2 ==========================================
# Compute and plot Lambda_esc
# Formula is derived in latex doc
lambda_esc = [(16.21 * i)/(1 - 2.28 * i) for i in ratios]
lambda_uncert = uncert * (16.21 / ((1-2.28*ratios)**2))

ratio_fit = power_law(boron_energy,a,b)
lambda_fit = 16.21 * ratio_fit / (1 - 2.28*ratio_fit)
"""
# Plot lambda_esc vs energy per nuclei
# Plot ratio and fit
plt.figure()
plt.xlabel("Energy per nuclei (GeV/n)", fontsize=14)
plt.ylabel(r"$\lambda_{\rm{esc}}$ (g/cm$^2$)", fontsize=14)
plt.scatter(boron_energy, lambda_esc, color='black')
plt.errorbar(boron_energy, lambda_esc, yerr=lambda_uncert, fmt='o', color='black', capsize=3, label="Data")
plt.plot(boron_energy, lambda_fit, color='red', label='Fit')
plt.xscale("log")
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()
"""
# ========================================== PROBLEM 1.3 ==========================================
p_avg = 1.67e-24 #1 H/cm3 --> g/cm3
L_avg = [i / p_avg for i in lambda_esc] # cm
L_avg = [i/3.086e24 for i in L_avg]

L_uncert = [i / p_avg for i in lambda_uncert] 
L_uncert = [i / 3.086e24 for i in L_uncert]

L_fit = lambda_fit / p_avg
L_fit = [i/3.086e24 for i in L_fit]
"""
# Plot L_avg vs energy per nuclei
# Plot ratio and fit
plt.figure()
plt.xlabel("Energy per nuclei (GeV/n)", fontsize=14)
plt.ylabel(r"$\bar{L}$ (MPc)", fontsize=14)
plt.scatter(boron_energy, L_avg, color='black')
plt.errorbar(boron_energy, L_avg, yerr=L_uncert, fmt='o', color='black', capsize=3, label="Data")
plt.plot(boron_energy, L_fit, color='red', label='Fit')
plt.xscale("log")
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()
"""
# ========================================== PROBLEM 1.4 ==========================================
# Fit power law for observed boron flux vs energy
# Perform curve fitting
# Fit for energies up to 200 GeV/n
fit_energy = boron_energy[np.where(boron_energy < 200000)]
fit_flux = boron_flux[np.where(boron_energy < 20000)]
boron_params, boron_covariance = curve_fit(power_law, fit_energy, fit_flux, maxfev=2000)

# Extracting coefficients
a_boron = boron_params[0]
b_boron = boron_params[1]
print(a_boron, b_boron)
"""
# Plot boron and carbon flux
plt.figure()
plt.xlabel("Energy per nuclei (GeV/n)", fontsize=14)
plt.ylabel("Observed Flux", fontsize=14)
plt.scatter(boron_energy, boron_flux, color='black')
plt.errorbar(boron_energy, boron_flux, yerr=boron_uncert, fmt='o', color='black', capsize=3, label="Boron")

plt.plot(boron_energy, power_law(boron_energy, a_boron, b_boron), color='red')
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale("log")
plt.yscale("log")
plt.show()
"""

# ========================================== PROBLEM 1.5 ==========================================
energy = np.linspace(10, 1000000, 1000000001)
ratio = [0.46 * i**(-0.35) for i in energy]
lambda_carbon = [(16.28 * i) / (1 - 2.283 * i) for i in ratio]

print(np.where(lambda_carbon == 0.167))
plt.plot(energy, lambda_carbon)
plt.xscale("log")
plt.show()