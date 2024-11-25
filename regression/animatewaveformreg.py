import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# Load the data from the results CSV file
data = pd.read_csv("results.csv")

# Extract unique iterations
iterations = data['Iteration'].unique()

# Set up the plot
fig, ax = plt.subplots()
scatter, = ax.plot([], [], 'o', label='Data')  # Scatter plot for data points
line, = ax.plot([], [], '-', label='Fit', color='red')  # Line for the fit
title = ax.text(0.5, 1.05, "", transform=ax.transAxes, ha="center", fontsize=12)
ax.legend()

# Initialize the plot
def init():
    ax.set_xlim(data['x'].min() - 1, data['x'].max() + 1)
    ax.set_ylim(data['y'].min() - 0.1, data['y'].max() + 0.1)
    scatter.set_data([], [])
    line.set_data([], [])
    title.set_text("")
    return scatter, line, title

# Update the plot for each frame
def update(frame):
    iter_data = data[data['Iteration'] == frame]

    x = iter_data['x']
    y = iter_data['y']
    scatter.set_data(x, y)

    sigma = iter_data['Sigma'].iloc[0]
    mu = iter_data['Mu'].iloc[0]

    # Compute the fit using the Gaussian model
    fit_x = np.linspace(x.min(), x.max(), 500)
    fit_y = (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((fit_x - mu) / sigma) ** 2)
    line.set_data(fit_x, fit_y)

    # Update the title with current iteration, sigma, and mu
    title.set_text(rf"Iteration: {frame}, $\sigma$: {sigma:.6f}, $\mu$: {mu:.6f}")
    return scatter, line, title

# Create the animation
ani = FuncAnimation(fig, update, frames=iterations, init_func=init, blit=False, repeat=False, interval=50)

# Save the animation as a video
# ani.save("fit_evolution.mp4", writer="ffmpeg", fps=20)

# Uncomment the line below to show the animation interactively
plt.show()

