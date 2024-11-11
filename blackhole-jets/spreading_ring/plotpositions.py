import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Open the .h5 file and read particle positions for every timestep
with h5py.File("particle_positions.h5", "r") as file:
    # Retrieve all timesteps in the file
    steps = sorted(file.keys(), key=lambda x: int(x.split('_')[1]))  # Sort steps by timestep order

    # Get the initial number of particles from the first timestep
    initial_num_particles = file[steps[0]].shape[0]
    x_positions = []
    y_positions = []

    # Populate position arrays only if the particle count matches the initial count
    for step in steps:
        data = file[step][:]
        if data.shape[0] == initial_num_particles:
            x_positions.append(data[:, 0])  # x coordinates
            y_positions.append(data[:, 1])  # y coordinates

    # Convert lists to numpy arrays for easier indexing
    x_positions = np.array(x_positions)
    y_positions = np.array(y_positions)
    num_timesteps = x_positions.shape[0]

# Calculate the skip factor for subsampling
target_frames = 959789  # Approximate frames needed for a 1-minute video at 24 fps
skip_factor = max(1, num_timesteps // target_frames)

# Set up the figure and axis
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-5, 5)  # Adjust limits based on your simulation box size
ax.set_ylim(-5, 5)
ax.set_xlabel("x position")
ax.set_ylabel("y position")
ax.set_title("Orbital Plane Trajectories of Particles")
scat = ax.scatter([], [], s=1)

# Ensure equal scaling for realistic orbital paths
ax.axis('equal')

# Initialize function to set up the first frame
def init():
    scat.set_offsets(np.empty((0, 2)))
    return scat,

# Update function for each frame in the animation
def update(frame):
    actual_frame = frame * skip_factor
    x = x_positions[actual_frame]
    y = y_positions[actual_frame]
    positions = np.column_stack((x, y))
    scat.set_offsets(positions)
    return scat,

# Create the animation with subsampled frames
ani = FuncAnimation(fig, update, frames=num_timesteps // skip_factor, init_func=init, blit=True, interval=1000 / 24)

# Save the animation as an MP4 file
ani.save("orbital_plane_trajectories.mp4", writer="ffmpeg", dpi=300)

# Display the animation (optional, can be removed if only saving is needed)
plt.show()

