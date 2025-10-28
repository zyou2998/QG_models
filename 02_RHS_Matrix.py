# =============================================================================
# Constructs the RHS matrix operator, R, for the barotropic model.
# This matrix represents advection by the background flow plus dissipation.
#
# Creator: Zhenyu You, Georgia Tech
# Contact: zyou2998@gmail.com
# =============================================================================

import numpy as np
import netCDF4 as nc
from scipy.sparse import diags
import os
from toolbox import *
from model_parameters import *

# --- 1. Setup File Paths ---
print(f"Generating RHS matrix for {season}.")
damp_file = in_dir + f"damp_ocean_{damp_ocean}_land_{damp_land}.nc"
gph_file = in_dir + f"GP_{d_lat}_{lat_start}_{lat_end}_{season}.nc"
out_file = out_dir + f"RHS_{season}_damp_o{damp_ocean}_l{damp_land}_{d_lat}_{lat_start}_{lat_end}.nc"

# --- 2. Load Input Data ---
with nc.Dataset(gph_file) as ds:
    Z1 = ds["data"][0].flatten() # Background GPH (upper level)
print(f"Loaded background flow for {season}.")

# --- 3. Calculate Advection Operator Terms ---
f = Corio()
mlat = get_lat_all()
Vort1_y = reconstruct(Vort_y(Y_mat=Y) @ Z1)
Vort1_x = reconstruct(Vort_x(Y_mat=Y) @ Z1)
Z1_y = reconstruct(diff_y(Y_mat=Y) @ Z1)
Z1_x = reconstruct(diff_x(Y_mat=Y) @ Z1)

RA1 = -diags(1/(r*f*np.cos(mlat)) * Vort1_y) @ diff_x() # Advection of planetary vorticity by v'
RA2 = -diags(1/(r*f*np.cos(mlat)) * Z1_x) @ Vort_y()    # Advection of relative vorticity by u_bar
RA3 = -diags(2*ome/(r**2*f)) @ diff_x()                # Beta effect (advection of f by u')
RA4 =  diags(1/(r*f) * Vort1_x) @ diff_y()             # Advection of planetary vorticity by u'
RA5 =  diags(1/(r*f) * Z1_y) @ Vort_x()                # Advection of relative vorticity by v_bar

# --- 4. Assemble and Modify the RHS Matrix ---
R = (RA1 + RA2 + RA3 + RA4 + RA5).toarray()

# Add hyperdiffusion
vort = Vort()
lap = Laplacian()
hyper_diff = lap @ lap @ vort
R -= hype_coef * hyper_diff.toarray()
print(f"Hyperdiffusion of {hype_coef:.1e} added.")

# Add physical damping
if damp_flag == 1:
    with nc.Dataset(damp_file) as ds:
        eps = ds["damp"][:].flatten()
    
    damp_op = diags(eps) @ vort
    R -= damp_op.toarray()
    print(f"Spatially varying damping added (ocean={damp_ocean}, land={damp_land} days).")
else:
    print("damp_flag is 0. No physical damping added.")

# --- 5. Save the RHS Matrix ---
if os.path.exists(out_file):
    os.remove(out_file)

with nc.Dataset(out_file, "w", format="NETCDF4") as ds:
    ds.description = "RHS matrix operator for the barotropic model."
    ds.createDimension("row", N)
    ds.createDimension("col", N)

    data_var = ds.createVariable("R", np.float64, ("row", "col"))
    data_var[:] = R

print(f"\nSuccessfully saved RHS matrix to: {os.path.basename(out_file)}")