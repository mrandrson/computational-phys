import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit

data = pd.read_csv('n_body_positions.csv')
data.columns = data.columns.str.strip()

central_mass = data[data['Mass'] == data['Mass'].max()]
central_x = central_mass['X'].iloc[0]
central_y = central_mass['Y'].iloc[0]

data['Radial_Distance'] = np.sqrt((data['X'] - central_x) ** 2 + (data['Y'] - central_y) ** 2)

particles = data['Particle'].unique()
time_steps = data['Step'].unique()



plt.figure(figsize=(10, 6))

mean_radius = []
massarr = []
for particle in particles:
    particle_data = data[data['Particle'] == particle]
    massarr.append(particle_data['Mass'].unique())
    mean_radius.append(np.mean(particle_data['Radial_Distance']))
    #plt.plot(particle_data['Step'], particle_data['Radial_Distance'], alpha=0.5, label=f'Particle {particle}')


mean_radius = np.array(mean_radius)/1.49e11
massarr = np.array(massarr)

hist_values, bin_edges = np.histogram(mean_radius[mean_radius > 0.25], bins='auto')
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

def fit_function(r, a, k):
    funcval = a*(1/r)**k
    return funcval 

popt, pcov = curve_fit(fit_function, bin_centers, hist_values, p0=[1, 1])

plt.hist(mean_radius[mean_radius > 0.25], bins = 'auto', histtype= 'step')
#plt.hist(massarr[massarr < 1e30], bins = 'auto', histtype= 'step', label = 'Mass')
fitexp = popt[1]
plt.plot(bin_centers, fit_function(bin_centers, *popt), 'r-', label=rf'Fit: ${popt[0]:.2f}\left(\frac{{AU}}{{r}}\right)^{{{fitexp:.2f}}}$')
# plt.plot(bin_centers, fit_function(bin_centers, *popt))
plt.xlabel(r'$r\ \left[AU\right]$')
plt.ylabel('Number of Particles')
plt.legend()
plt.show()
