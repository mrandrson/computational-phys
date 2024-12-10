import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

def plot_leaf_cells(ax, leaf_cells_csv):
    if not os.path.exists(leaf_cells_csv):
        print(f"Error: {leaf_cells_csv} does not exist.")
        return
    data = pd.read_csv(leaf_cells_csv)
    for _, row in data.iterrows():
        rect = Rectangle((row['xmin'], row['ymin']),
                         row['xmax'] - row['xmin'],
                         row['ymax'] - row['ymin'],
                         edgecolor='lightgray', facecolor='none', lw=0.5)
        ax.add_patch(rect)

def plot_active_cells(ax, active_cells_csv):
    if not os.path.exists(active_cells_csv):
        print(f"Error: {active_cells_csv} does not exist.")
        return
    data = pd.read_csv(active_cells_csv)
    for _, row in data.iterrows():
        rect = Rectangle((row['cell_xmin'], row['cell_ymin']),
                         row['cell_xmax'] - row['cell_xmin'],
                         row['cell_ymax'] - row['cell_ymin'],
                         edgecolor='#9467bd', facecolor='none', lw=1)
        ax.add_patch(rect)

def plot_means(ax, means_csv):
    if not os.path.exists(means_csv):
        print(f"Error: {means_csv} does not exist.")
        return
    data = pd.read_csv(means_csv)
    if 'mass' not in data.columns:
        print(f"Error: 'mass' column not found in {means_csv}.")
        return
    masses = data['mass']
    sizes = (masses - masses.min()) / (masses.max() - masses.min()) * 100 + 20
    scatter = ax.scatter(data['mean_x'], data['mean_y'], s=sizes, color='#d62728', zorder=3, label='Cell Centers of Mass')

def plot_connections(ax, connections_csv):
    if not os.path.exists(connections_csv):
        print(f"Error: {connections_csv} does not exist.")
        return
    data = pd.read_csv(connections_csv)
    for _, row in data.iterrows():
        ax.plot([row['start_x'], row['end_x']], [row['start_y'], row['end_y']],
                '-', color='#2ca02c', lw=0.5, zorder=1)

def plot_query_position(ax, query_csv):
    if not os.path.exists(query_csv):
        print(f"Error: {query_csv} does not exist.")
        return
    data = pd.read_csv(query_csv)
    ax.scatter(data['query_x'], data['query_y'], color='blue', s=50, zorder=4, label='Query Position')


def main():
    fig, ax = plt.subplots(figsize=(10, 10))

    plot_leaf_cells(ax, 'quadtree_leaf_cells.csv')
    plot_active_cells(ax, 'quadtree_active_cells_with_particles.csv')
    plot_means(ax, 'quadtree_active_means.csv')
    plot_connections(ax, 'quadtree_active_connections.csv')
    plot_query_position(ax, 'quadtree_query_position.csv')

    ax.set_aspect('equal', adjustable='box')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())

    #plt.show()
    plt.savefig('0.1AUquadtree.png')

if __name__ == "__main__":
    main()

