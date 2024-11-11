from manim import *
import pandas as pd
import struct

def load_data(filename):
    format_str = 'ddd'
    data = []

    with open(filename, 'rb') as file:
        while True:
            bytes_chunk = file.read(struct.calcsize(format_str))
            if not bytes_chunk:
                break
            x, y, z = struct.unpack(format_str, bytes_chunk)
            data.append([x, y, z])

    df = pd.DataFrame(data, columns=['x', 'y', 'z'])
    return df

class DotMoving(ThreeDScene):
    def construct(self):
        df = load_data('output.bin')

        points = [np.array([row['x'], row['y'], row['z']]) for _, row in df.iterrows()]

        axes = ThreeDAxes()  # 3D axes
        self.set_camera_orientation(phi=75 * DEGREES, theta=-45 * DEGREES)  # Initial camera angle
        self.add(axes)

        dot = Dot3D(points[0], color=RED)  # Dot in 3D
        self.add(dot)

        origin = np.array([0, 0, 0])  # Origin

        line = Line3D(origin, points[0], color=BLUE)  # Line from origin to the first point in 3D
        self.add(line)

        # Variables to adjust camera orientation step
        total_points = len(points)
        initial_phi = 75 * DEGREES
        initial_theta = -45 * DEGREES

        for i, point in enumerate(points[1:], 1):
            new_line = Line3D(origin, point, color=BLUE)

            # Calculate camera rotation steps based on the number of points
            new_phi = initial_phi + (i / total_points) * 30 * DEGREES  # Slow rotation in phi
            new_theta = initial_theta + (i / total_points) * 60 * DEGREES  # Slow rotation in theta

            # Move dot and line
            self.play(dot.animate.move_to(point), Transform(line, new_line), run_time=0.1)

            # Update the camera orientation manually
            self.set_camera_orientation(phi=new_phi, theta=new_theta)

