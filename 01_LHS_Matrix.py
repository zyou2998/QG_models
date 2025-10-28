# =============================================================================
# Constructs the LHS matrix operator, L, for the barotropic model.
# For the barotropic vorticity equation, L is the Laplacian operator.
#
# Creator: Zhenyu You, Georgia Tech
# Contact: zyou2998@gmail.com
# =============================================================================

import numpy as np
import netCDF4 as nc
import os
from toolbox import *
from model_parameters import *

# --- 1. Setup ---
print(f"Generating LHS matrix.")
out_file = out_dir + f"LHS_{d_lat}_{lat_start}_{lat_end}.nc"

# --- 2. Calculate the LHS Matrix ---
L = Vort().toarray()

# --- 3. Save the LHS Matrix to a NetCDF File ---
if os.path.exists(out_file):
    os.remove(out_file)

with nc.Dataset(out_file, "w", format="NETCDF4") as ds:
    ds.description = "LHS matrix operator (Laplacian) for the barotropic model."
    ds.createDimension("row", N)
    ds.createDimension("col", N)

    data_var = ds.createVariable("L", np.float64, ("row", "col"))
    data_var[:] = L

print(f"Successfully saved LHS matrix to: {os.path.basename(out_file)}")