import numpy as np
import matplotlib.pyplot as plt

# Load the data from the text file
file_path = "./builddir/fft.dat"  # Replace with your file path
data = np.loadtxt(file_path, skiprows=1)

# Extract columns for plotting
x1, y1 = data[:, 0], data[:, 1]  # First and second columns
x2, y2 = data[:, 2], data[:, 3]  # Third and fourth columns

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # Two plots side by side

# Plot 2nd column vs. 1st column
axes[0].plot(x1, y1, marker='o', linestyle='-', color='b')
axes[0].set_title("signal")
axes[0].set_xlabel("time in s")
axes[0].set_ylabel("amplitude")
axes[0].set_xlim(0, 2)  # Restrict x-axis range to [0, 2]

# Plot 4th column vs. 3rd column
axes[1].plot(x2, y2, marker='o', linestyle='None', color='r', label='Data Points')
axes[1].vlines(x2, 0, y2, colors='r', linestyles='dashed')
axes[1].set_title("FFT")
axes[1].set_xlabel("frequency")
axes[1].set_ylabel("amplitude")
axes[1].axhline(0, color='black', linewidth=0.8)  # Optional: Add x-axis line
axes[1].set_xlim(-5, 5)  # Restrict x-axis range to [0, 2*pi]

# Adjust layout and show the plot
plt.tight_layout()
plt.savefig("fft.pdf")
