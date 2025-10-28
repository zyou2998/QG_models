# =============================================================================
# Calculate and save damping coefficients for a barotropic model based on a
# land-sea mask from Sea Surface Temperature (SST) data.
#
# Creator: Zhenyu You, Georgia Tech
# Contact: zyou2998@gmail.com
# =============================================================================

import numpy as np
import netCDF4 as nc
from scipy.interpolate import RegularGridInterpolator
import os
from model_parameters import *

# --- 1. Load and Interpolate SST to Create a Land-Sea Mask ---
print("Loading and interpolating SST data...")
with nc.Dataset(in_dir + "sst.mon.ltm.1981-2010.nc", "r") as f:
    sst_raw = f.variables["sst"][0, ::-1, :]   # Flip latitude to be S-N
    lat_raw = f.variables["lat"][::-1]
    lon_raw = f.variables["lon"][:]

# Create an interpolator function from the raw SST grid
interpolator = RegularGridInterpolator(
    (lat_raw, lon_raw), sst_raw, method='linear', bounds_error=False, fill_value=0
)
# A small SST value is assumed to be land
sst_on_model_grid = interpolator((lat_grid, lon_grid))
is_land = sst_on_model_grid < 0.1

# --- 2. Calculate Damping Coefficients ---
print("Calculating damping coefficients...")
# Convert damping timescales (days) to damping rates (1/s)
s_per_day = 86400.0
rate_ocean = 1.0 / (damp_ocean * s_per_day)
rate_land = 1.0 / (damp_land * s_per_day)
# Use the land mask to assign damping values
damp = np.where(is_land, rate_land, rate_ocean)

# --- 3. Save Damping Field to NetCDF File ---

outfile = in_dir + f"damp_ocean_{damp_ocean}_land_{damp_land}.nc"
print(f"Saving damping field to: {os.path.basename(outfile)}")

if os.path.exists(outfile):
    os.remove(outfile)

with nc.Dataset(outfile, "w", format="NETCDF4") as ds:
    ds.title = "Damping Coefficient (1/s) for Barotropic Model"
    ds.createDimension("lat", Y - 4)
    ds.createDimension("lon", X)
    
    lat_var = ds.createVariable("lat", float, ("lat",))
    lon_var = ds.createVariable("lon", float, ("lon",))
    damp_var = ds.createVariable("damp", float, ("lat", "lon"))
    
    lat_var[:] = lat[2:Y-2]
    lon_var[:] = lon
    damp_var[:] = damp[2:Y-2, :]

print("\nDamping calculation complete.")