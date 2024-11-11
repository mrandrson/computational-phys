import struct
import matplotlib.pyplot as plt

with open('output.bin', 'rb') as f:
    data = f.read()

data_size = len(data)
steps = data_size // (3 * 8)  # There are 'steps' number of double values for each of t, x, y

t = struct.unpack(f'{steps}d', data[:steps * 8])  # First part is t values
x = struct.unpack(f'{steps}d', data[steps * 8:2 * steps * 8])  # Second part is x values
y = struct.unpack(f'{steps}d', data[2 * steps * 8:])  # Third part is y values

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()

