# =============================================================================
# Calculates the optimal modes for a finite time interval (tau) for the
# barotropic model.
#
# Creator: Zhenyu You, Georgia Tech
# Contact: zyou2998@gmail.com
# =============================================================================

import numpy as np
import netCDF4 as nc
import os
from scipy.linalg import eig, lu_factor, lu_solve
from model_parameters import *

# --- 1. Setup File Paths ---
print(f"Calculating optimal modes for tau = {tau} days. Season: {season}")
damp_str = f"damp_o{damp_ocean}_l{damp_land}"
print(f"Model Config: Damping={damp_str}")

# Define input and output file paths
in_file = out_dir + f"normal_{season}_{damp_str}_{d_lat}_{lat_start}_{lat_end}.nc"
if loc_flag == 0:
    out_file = out_dir + f"optimal_{season}_tau{tau}_{damp_str}_{d_lat}_{lat_start}_{lat_end}.nc"
else:
    out_file = out_dir + f"optimal_{season}_local_{opt_region}_tau{tau}_{damp_str}_{d_lat}_{lat_start}_{lat_end}.nc"

# --- 2. Load Normal Mode Data ---
with nc.Dataset(in_file) as ds:
    sigma = ds["rate_real"][:] + 1j * ds["rate_imag"][:]
    P = ds["mode_real"][:] + 1j * ds["mode_imag"][:]
print("Normal mode data loaded.")

# --- 3. Define the Projection Operator for Local Optimization ---
eps = 0.01
if loc_flag == 1:
    print(f"Calculating local optimals for region: {opt_region}")
    M = np.full((Y - 4, X), eps)
    lat_interior = lat[2:Y-2]
    opt_lat_slice = slice(np.abs(lat_interior - opt_minlat).argmin(), np.abs(lat_interior - opt_maxlat).argmin() + 1)
    opt_lon_slice = slice(np.abs(lon - opt_minlon).argmin(), np.abs(lon - opt_maxlon).argmin() + 1)
    M[opt_lat_slice, opt_lon_slice] = 1.0
    D = M.flatten()
else:
    print("Calculating global optimals.")
    D = np.ones((N,))

# --- 4. Construct and Solve the Optimal Mode Eigenproblem ---
Gamma = np.diag(np.exp(sigma * tau * 86400)) # Propagator matrix
B_1 = np.conj(P.T) @ P                       # Norm matrix at t=0
B_0 = np.conj(P.T) @ (P * D[:, np.newaxis])   # Weighted norm matrix at t=0
B_tau = np.conj(Gamma.T) @ B_0 @ Gamma       # Weighted norm matrix at t=tau
print("Matrices for eigenproblem obtained.")

# Solve the generalized eigenvalue problem: B_tau * a = lambda * B_1 * a
lu, piv = lu_factor(B_1)
mat = lu_solve((lu, piv), B_tau)
print("Linear operator obtained.")
lamda, a = eig(mat)
print("Optimal modes obtained.")

# Sort by amplification factor (real part of lambda), descending
sorted_indices = np.argsort(np.real(lamda))[::-1]
D_sorted = np.real(lamda[sorted_indices])
a_sorted = a[:, sorted_indices]

# --- 5. Evolve Top Optimal Modes Forward in Time ---
stream = np.zeros((optimal_number, int(evo/td)+1, N), dtype=complex)
time_steps = np.arange(0, evo + td, td)

for ii in range(optimal_number):
    print("Optimal Mode:", ii+1)
    for it, time_val in enumerate(time_steps):
        stream[ii, it, :] = P @ (np.exp(sigma * time_val * 86400) * a_sorted[:, ii])

# --- 6. Save Results ---
if os.path.exists(out_file):
    os.remove(out_file)

with nc.Dataset(out_file, "w", format="NETCDF4") as ds:
    ds.description = f"Optimal modes and their time evolution for tau={tau} days."
    ds.createDimension("mode", optimal_number)
    ds.createDimension("time", len(time_steps))
    ds.createDimension("space", N)

    data_var = ds.createVariable("data", "f8", ("mode", "time", "space"))
    imag_var = ds.createVariable("imag", "f8", ("mode", "time", "space"))
    rate_var = ds.createVariable("rate", "f8", ("mode",))

    data_var.description = "Real part of the streamfunction evolution."
    imag_var.description = "Imaginary part of the streamfunction evolution."
    rate_var.description = "Amplification factor (lambda) over tau."

    data_var[:] = np.real(stream)
    imag_var[:] = np.imag(stream)
    rate_var[:] = D_sorted[0:optimal_number]

print(f"\nSuccessfully saved optimal modes to: {os.path.basename(out_file)}")