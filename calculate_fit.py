import sys 
import numpy as np
from mpmath import mp
import scipy.stats as st

# Set the desired precision for calculations 
mp.dps = 10 # dps = decimal places

# Check if arguments were provided 
if len(sys.argv) != 4:
    print("Usage: python calculate_fit.py <mean> <variance> <output_filepath>")
    sys.exit(1)

# Get parameters from command line 
# sys.argv[0] is the script name itself
# sys.argv[1] is the first argument, sys.argv[2] the second
mu = float(sys.argv[1])
sigma2 = float(sys.argv[2])
base_path = sys.argv[3]

# Compute k and theta (Method of Moments, rounded version)
term = 4 * sigma2 - mu**2
if term <= 0:
  print(f"Error: 4*sigma2 - mu**2 is not positive. Cannot calculate k and theta.")
  exit()

# Calculate k as a float first
k_float = (5 * mu**2) / term
#k_float = 6
# Apply conditional rounding
if k_float > 0.9:
    k = round(k_float)
    print(f"Calculated k ({k_float:.4f}) > 0.9, rounding to integer {k}.")
else:
    k = k_float
    print(f"Calculated k ({k_float:.4f}) <= 0.9, keeping as float.")

# Calculate theta using the final k
if k == 0:
    print(f"Error: Final k is 0. Cannot calculate theta.")
    exit()
theta = 1.5*mu / k

print(f"Received parameters: mu={mu}, sigma2={sigma2}")
print(f"Calculated parameters: k = {k}, theta = {theta}")


# Define functions g(q) and f(q) 
def g(q, k, theta):
    return (1.0 / float(mp.gamma(k)) / theta**k) * q**(k-1) * np.exp(-q/theta)
    

def f(q, k, theta):
    return float(mp.gammainc(k - 1, a=q/theta)) / (float(mp.gamma(k)) * theta)

# Generate the points for the curve 
# Determine a good range for based on the mean and variance
q_max = mu + 3 * np.sqrt(sigma2)
q_values = np.linspace(q_max / 500.0, q_max, 500) # Start from a small positive number

# Create an empty array to store the results
p_values = np.zeros_like(q_values)
g_values = np.zeros_like(q_values)
f_values = np.zeros_like(q_values)

# Loop through each q value and calculate P(q) individually
for i, q in enumerate(q_values):
    g_val = g(q, k, theta)
    f_val = f(q, k, theta)
    g_values[i] = g_val
    f_values[i] = f_val
    p_values[i] = (1.0/3.0) * g_val + (2.0/3.0) * f_val

# Save results
np.savetxt(f"{base_path}_P_fit.dat", np.column_stack((q_values, p_values)), header="q P(q)")
np.savetxt(f"{base_path}_g_fit.dat", np.column_stack((q_values, g_values)), header="q g(q)")
np.savetxt(f"{base_path}_f_fit.dat", np.column_stack((q_values, f_values)), header="q f(q)")

print(f"Analytical curve data saved to {base_path}_*.dat files") 
