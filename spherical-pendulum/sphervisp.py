from vispy import app, scene
from vispy.geometry import create_sphere
import numpy as np
import struct

with open('output.bin', 'rb') as f:
    data = f.read()

data_size = len(data)
steps = data_size // (3 * 8)

x, y, z = [], [], []
for i in range(steps):
    offset = i * 3 * 8
    x_val, y_val, z_val = struct.unpack('ddd', data[offset:offset + 24])
    x.append(x_val)
    y.append(y_val)
    z.append(z_val)

canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), bgcolor='white', show=True)
view = canvas.central_widget.add_view()
view.camera = scene.TurntableCamera(elevation=30, azimuth=30)

axes = scene.visuals.XYZAxis(parent=view.scene)

pendulum_line = scene.visuals.Line(
    pos=np.array([[0, 0, 0], [x[0], y[0], z[0]]]),
    color='black',  # Change line color to black to contrast with the white background
    parent=view.scene
)

sphere_data = create_sphere(20, 40, radius=0.05)
bob = scene.visuals.Mesh(meshdata=sphere_data, color='red', parent=view.scene)
bob.transform = scene.transforms.MatrixTransform()
bob.transform.translate((x[0], y[0], z[0]))

def update(event):
    global step_index
    if step_index >= steps:
        step_index = 0

    pendulum_line.set_data(pos=np.array([[0, 0, 0], [x[step_index], y[step_index], z[step_index]]]))

    bob.transform.reset()  # Reset transformation
    bob.transform.translate((x[step_index], y[step_index], z[step_index]))

    step_index += 1

step_index = 0
timer = app.Timer('auto', connect=update, start=True)

if __name__ == '__main__':
    app.run()

