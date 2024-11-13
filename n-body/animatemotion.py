import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

data = pd.read_csv('n_body_positions.csv')
print(data)
data.columns = data.columns.str.strip()

particles = data['Particle'].unique()
num_steps = data['Step'].max() + 1

fig, ax = plt.subplots()
fig.patch.set_facecolor('black')
ax.set_facecolor('black')
ax.set_xlim(data['X'].min() - 1, data['X'].max() + 1)
ax.set_ylim(data['Y'].min() - 1, data['Y'].max() + 1)
ax.set_title("N-Body Simulation")
ax.set_xlabel("X Position")
ax.set_ylabel("Y Position")

scatter = ax.scatter([], [], s=[], c='white', alpha=1)

central_mass = data[data['Mass'] == data['Mass'].max()]
central_mass_x = central_mass['X'].iloc[0]
central_mass_y = central_mass['Y'].iloc[0]

star_marker = ax.scatter([central_mass_x], [central_mass_y], s=100, color='orange', marker='*', label='Central Mass', alpha = 0.5)

def init():
    scatter.set_offsets([[0, 0]])
    return scatter, star_marker

def update(step):
    current_data = data[(data['Step'] == step) & (data['Mass'] < data['Mass'].max())]  # Exclude central mass
    
    positions = current_data[['X', 'Y']].values
    masses = current_data['Mass'].values
    sizes = 10 * (masses / masses.max())  # Scale sizes based on mass

    scatter.set_offsets(positions)
    scatter.set_sizes(sizes)

    return scatter, star_marker

ani = animation.FuncAnimation(fig, update, frames=num_steps, init_func=init, blit=True, interval=10)

# Optional: Save the animation as a .mp4 or .gif file
ani.save('n_body_simulation.mp4', writer='ffmpeg', fps=20)

# plt.show()

