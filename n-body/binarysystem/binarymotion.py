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
ax.tick_params(colors='black')

interval = 1.2e12

scatter = ax.scatter([], [], s=[], c='white', alpha=1)

sorted_masses = data[['Mass']].drop_duplicates().sort_values(by='Mass', ascending=False)
mass1, mass2 = sorted_masses['Mass'].iloc[0], sorted_masses['Mass'].iloc[1]

star1_data = data[(data['Mass'] == mass1) & (data['Step'] == 0)]
star2_data = data[(data['Mass'] == mass2) & (data['Step'] == 0)]
star1_x, star1_y = star1_data['X'].iloc[0], star1_data['Y'].iloc[0]
star2_x, star2_y = star2_data['X'].iloc[0], star2_data['Y'].iloc[0]

star1_marker = ax.scatter([star1_x], [star1_y], s=150, color='orange', marker='*', label='Star 1', alpha=0.8)
star2_marker = ax.scatter([star2_x], [star2_y], s=150, color='red', marker='*', label='Star 2', alpha=0.8)

def init():
    scatter.set_offsets([[0, 0]])
    return scatter, star1_marker, star2_marker

def compute_center_of_mass(step_data):
    total_mass = step_data['Mass'].sum()
    com_x = (step_data['X'] * step_data['Mass']).sum() / total_mass
    com_y = (step_data['Y'] * step_data['Mass']).sum() / total_mass
    return com_x, com_y

def update(step):
    current_data = data[data['Step'] == step]
    com_x, com_y = compute_center_of_mass(current_data)
    ax.set_xlim(com_x - interval, com_x + interval)
    ax.set_ylim(com_y - interval, com_y + interval)

    other_particles = current_data[current_data['Mass'] < mass2]
    positions = other_particles[['X', 'Y']].values
    masses = other_particles['Mass'].values
    sizes = 10 * (masses / masses.max())

    scatter.set_offsets(positions)
    scatter.set_sizes(sizes)

    star1_position = current_data[current_data['Mass'] == mass1][['X', 'Y']].values
    star2_position = current_data[current_data['Mass'] == mass2][['X', 'Y']].values
    star1_marker.set_offsets(star1_position)
    star2_marker.set_offsets(star2_position)

    return scatter, star1_marker, star2_marker

ani = animation.FuncAnimation(fig, update, frames=num_steps, init_func=init, blit=True, interval=10)
# ani.save('binary_simulation.mp4', writer='ffmpeg', fps=20)

plt.show()
