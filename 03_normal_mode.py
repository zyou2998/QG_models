# =============================================================================
# Calculates the normal modes and growth rates (eigenvalues) of the
# barotropic model by solving the eigenvalue problem L*phi = sigma*R*phi.
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
print("Calculating normal modes...")
damp_str = f"damp_o{damp_ocean}_l{damp_land}"
print(f"Model Config: {season}, Damping={damp_str}")

# Define input and output file paths
L_file = out_dir + f"LHS_{d_lat}_{lat_start}_{lat_end}.nc"
R_file = out_dir + f"RHS_{season}_{damp_str}_{d_lat}_{lat_start}_{lat_end}.nc"
out_file = out_dir + f"normal_{season}_{damp_str}_{d_lat}_{lat_start}_{lat_end}.nc"

# --- 2. Load L and R ---
with nc.Dataset(L_file) as ds:
    L = ds['L'][:]
with nc.Dataset(R_file) as ds:
    R = ds['R'][:]
print("LHS and RHS matrices loaded successfully.")

# --- 3. Solve the Eigenvalue Problem ---
lu, piv = lu_factor(L)
mat = lu_solve((lu, piv), R)
print("Linear operator obtained.")
sigma, phi = eig(mat)
print("Normal modes obtained.")

# --- 4. Sort Modes by Growth Rate ---
# Sort modes by the real part of sigma in descending order.
sorted_indices = np.argsort(np.real(sigma))[::-1]
sigma = sigma[sorted_indices]
phi = phi[:, sorted_indices]

# --- 5. Save Results ---
n_mode = N
if os.path.exists(out_file):
    os.remove(out_file)

with nc.Dataset(out_file, "w", format="NETCDF4") as ds:
    ds.description = "Normal modes and growth rates of the barotropic model."
    ds.createDimension("space", n_mode)
    ds.createDimension("mode", n_mode)

    mode_real = ds.createVariable("mode_real", 'f8', ("space", "mode"))
    mode_imag = ds.createVariable("mode_imag", 'f8', ("space", "mode"))
    rate_real = ds.createVariable("rate_real", 'f8', ("mode",))
    rate_imag = ds.createVariable("rate_imag", 'f8', ("mode",))

    mode_real.description = "Real part of eigenvectors (modes)."
    mode_imag.description = "Imaginary part of eigenvectors (modes)."
    rate_real.description = "Eigenvalue real part (growth rate, 1/s)."
    rate_imag.description = "Eigenvalue imaginary part (frequency, rad/s)."

    mode_real[:] = phi.real
    mode_imag[:] = phi.imag
    rate_real[:] = sigma.real
    rate_imag[:] = sigma.imag

print(f"\nSuccessfully saved normal modes to: {os.path.basename(out_file)}")