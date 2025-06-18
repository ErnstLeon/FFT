import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Check for correct usage
if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} path/to/fft.dat")
    sys.exit(1)

# Get the input file path from command line argument
file_path = sys.argv[1]

# Load the data from the text file
data = np.loadtxt(file_path, skiprows=1)

# Extract columns for plotting
x1, y1 = data[:, 0], data[:, 1]
x2, y2 = data[:, 2], data[:, 3]

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Plot signal
axes[0].plot(x1, y1, marker='o', linestyle='-', color='b')
axes[0].set_title("signal")
axes[0].set_xlabel("time in s")
axes[0].set_ylabel("amplitude")
axes[0].set_xlim(0, 2)

# Plot FFT
axes[1].plot(x2, y2, marker='o', linestyle='None', color='r', label='Data Points')
axes[1].vlines(x2, 0, y2, colors='r', linestyles='dashed')
axes[1].set_title("FFT")
axes[1].set_xlabel("frequency")
axes[1].set_ylabel("amplitude")
axes[1].axhline(0, color='black', linewidth=0.8)
axes[1].set_xlim(-5, 5)

# Determine output path (same directory as input file)
output_dir = os.path.dirname(os.path.abspath(file_path))
output_path = os.path.join(output_dir, "fft.png")

# Save and close
plt.tight_layout()
plt.savefig(output_path)
plt.close()

print(f"Saved FFT plot to: {output_path}")