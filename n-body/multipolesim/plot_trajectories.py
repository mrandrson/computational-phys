import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.animation import FuncAnimation

def plot_trajectories(csv_file, output_image="trajectories.png"):
    dat = pd.read_csv(csv_file)
    
    x = lambda p: dat[dat['Index']==p]['x']
    y = lambda p: dat[dat['Index']==p]['y']

    for i in dat['Index'].unique():
        plt.plot(x(i), y(i))

    plt.show()

def animate_trajectories(csv_file, output_animation="trajectories.mp4", interval=2e12):
    data = pd.read_csv(csv_file)

    required_columns = {"Step", "Index", "x", "y", "mass"}
    if not required_columns.issubset(data.columns):
        raise ValueError(f"The file must contain the following columns: {required_columns}")

    data["Step"] = data["Step"].astype(int)
    data["Index"] = data["Index"].astype(int)
    data["mass"] = data["mass"].astype(float)

    num_steps = data["Step"].max() + 1

    unique_masses = data["mass"].unique()
    sorted_masses = np.sort(unique_masses)[::-1]
    if len(sorted_masses) >= 2:
        mass1, mass2 = sorted_masses[0], sorted_masses[0]
    else:
        mass1 = sorted_masses[0]
        mass2 = None

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.tick_params(colors='black')

    scatter = ax.scatter([], [], s=[], c='white', alpha=1)

    if mass2 is not None:
        star1_data = data[(data["mass"] == mass1) & (data["Step"] == 0)]
        star2_data = data[(data["mass"] == mass2) & (data["Step"] == 0)]
        if not star1_data.empty:
            star1_x, star1_y = star1_data["x"].iloc[0], star1_data["y"].iloc[0]
        else:
            star1_x, star1_y = 0, 0
        if not star2_data.empty:
            star2_x, star2_y = star2_data["x"].iloc[0], star2_data["y"].iloc[0]
        else:
            star2_x, star2_y = 0, 0

        star1_marker = ax.scatter([star1_x], [star1_y], s=200, color='orange', marker='*', label='Central Mass 1', alpha=0.8)
        star2_marker = ax.scatter([star2_x], [star2_y], s=200, color='red', marker = '*', label='Central Mass 2', alpha=0.8)

    else:
        central_data = data[(data["mass"] == mass1) & (data["Step"] == 0)]
        if not central_data.empty:
            central_x, central_y = central_data["x"].iloc[0], central_data["y"].iloc[0]
        else:
            central_x, central_y = 0, 0
        star1_marker = ax.scatter([central_x], [central_y], s=200, color='red', marker='*', label='Central Mass', alpha=0.8)
        star2_marker = None

    def init():
        scatter.set_offsets([[0, 0]])
        if mass2 is not None:
            return scatter, star1_marker, star2_marker
        else:
            return scatter, star1_marker

    global_max_interval = [0]

    def compute_center_of_mass(step_data):
        total_mass = step_data['mass'].sum()
        com_x = (step_data['x'] * step_data['mass']).sum() / total_mass
        com_y = (step_data['y'] * step_data['mass']).sum() / total_mass
        return com_x, com_y
    
    def update(step):
        current_data = data[data['Step'] == step]
        com_x, com_y = compute_center_of_mass(current_data)

        if mass2 is not None:
            central_masses = current_data[(current_data['mass'] == mass1) | (current_data['mass'] == mass2)]
        else:
            central_masses = current_data[current_data['mass'] == mass1]

        if not central_masses.empty:
            max_distance = np.sqrt((central_masses['x'] - com_x)**2 + (central_masses['y'] - com_y)**2).max()
            dynamic_interval = max(2 * max_distance, 1.0)
        else:
            dynamic_interval = 1.0

        global_max_interval[0] = max(global_max_interval[0], dynamic_interval)

        ax.set_xlim(com_x - global_max_interval[0], com_x + global_max_interval[0])
        ax.set_ylim(com_y - global_max_interval[0], com_y + global_max_interval[0])

        if mass2 is not None:
            orbiting_particles = current_data[(current_data['mass'] != mass1) & (current_data['mass'] != mass2)]
        else:
            orbiting_particles = current_data[current_data['mass'] != mass1]

        positions = orbiting_particles[['x', 'y']].values
        masses = orbiting_particles['mass'].values

        if len(masses) > 0:
            sizes = 10 * (masses / masses.max())
        else:
            sizes = []

        scatter.set_offsets(positions)
        scatter.set_sizes(sizes)

        if mass2 is not None:
            star1_pos = current_data[current_data['mass'] == mass1][['x', 'y']].values
            star2_pos = current_data[current_data['mass'] == mass2][['x', 'y']].values
            if len(star1_pos) > 0:
                star1_marker.set_offsets(star1_pos)
            if len(star2_pos) > 0:
                star2_marker.set_offsets(star2_pos)
            return scatter, star1_marker, star2_marker
        else:
            central_pos = current_data[current_data['mass'] == mass1][['x', 'y']].values
            if len(central_pos) > 0:
                star1_marker.set_offsets(central_pos)
            return scatter, star1_marker

    ani = animation.FuncAnimation(
        fig, update, frames=num_steps, init_func=init, blit=False, interval=1
    )
    #ani.save(output_animation, writer='ffmpeg', fps=50)
    plt.show()


csv_file = "particle_data.csv"
animate_trajectories(csv_file)
#plot_trajectories(csv_file)
