import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

data = pd.read_csv("results.csv")

iterations = data['Iteration'].unique()

fig, ax = plt.subplots()
scatter, = ax.plot([], [], 'o', label='Data')
line, = ax.plot([], [], '-', label='Fit', color='red')
title = ax.text(0.5, 1.05, "", transform=ax.transAxes, ha="center", fontsize=12)
ax.legend()

def init():
    ax.set_xlim(data['x'].min() - 1, data['x'].max() - 80)
    ax.set_ylim(data['y'].min() - 0.1, data['y'].max() + 0.1)
    scatter.set_data([], [])
    line.set_data([], [])
    title.set_text("")
    return scatter, line, title

def update(frame):
    iter_data = data[data['Iteration'] == frame]

    print(f"Processing frame: {frame}")
    print(iter_data.head())

    x = iter_data['x']
    y = iter_data['y']
    scatter.set_data(x, y)

    sigma = iter_data['Sigma'].iloc[0]
    mu = iter_data['Mu'].iloc[0]

    print(f"Sigma: {sigma}, Mu: {mu}")

    fit_x = np.linspace(x.min(), x.max(), 500)
    fit_y = (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((fit_x - mu) / sigma) ** 2)
    line.set_data(fit_x, fit_y)

    title.set_text(rf"Iteration: {frame}, $\sigma$: {sigma:.6f}, $\mu$: {mu:.6f}")
    return scatter, line, title

ani = FuncAnimation(fig, update, frames=iterations, init_func=init, blit=False, repeat=False, interval = 1)

ani.save("fit_evolution.mp4", writer="ffmpeg", fps=1000)
# plt.show()  # Uncomment to display interactively

