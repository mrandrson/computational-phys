import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def gaussian(x, sigma, mu):
    return (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * pow((x - mu) / sigma, 2))

def error(args):
    sigma, mu, xvals, yvals = args
    return sum(gaussian(xvals.iloc[i], sigma, mu) - yvals.iloc[i] for i in range(len(xvals)))

def main():
    dat = pd.read_csv('results.csv')

    sigmavals = dat['Sigma'][::100]
    muvals = dat['Mu'][::100]
    yvals = dat['y'][::100]
    xvals = dat['x'][::100]

    sigma_range = np.linspace(sigmavals.min(), sigmavals.max(), 50)
    mu_range = np.linspace(muvals.min(), muvals.max(), 50)
    sigma_grid, mu_grid = np.meshgrid(sigma_range, mu_range)

    points = [(sigma_grid[i, j], mu_grid[i, j], xvals, yvals)
              for i in range(sigma_grid.shape[0]) for j in range(sigma_grid.shape[1])]

    with ProcessPoolExecutor(max_workers=16) as executor:
        error_values_flat = list(executor.map(error, points))

    error_values = np.array(error_values_flat).reshape(sigma_grid.shape)

    plt.figure(figsize=(10, 6))
    plt.contourf(sigma_grid, mu_grid, error_values, levels=50, cmap='viridis')
    plt.colorbar(label='Error')
    plt.xlabel('Sigma')
    plt.ylabel('Mu')
    plt.title('Error Landscape')
    plt.plot(sigmavals, muvals, color='white', marker='o', markersize=3, label='Sigma vs Mu')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()

