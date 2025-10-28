# =============================================================================
# Parameters for the one-layer QG (barotropic vorticity) model.
#
# Creator: Zhenyu You, Georgia Tech
# Contact: zyou2998@gmail.com
# Last Modified: 09/02/2025
# =============================================================================

import numpy as np

# --- I/O and Experiment Setup ---
model_dir = "../"
in_dir    = f"{model_dir}input_data/"
out_dir   = f"{model_dir}output_data/"
pic_dir   = f"{model_dir}pics/"

season = "JA"

# --- Time Integration & Optimization ---
# tau: optimiation time (days)
tau = 10.0
evo = 30.0          # Integration time
td = 0.5            # Time interval
optimal_number = 30 # Number of optimal modes

# --- Physical Parameters ---
damp_flag = 1
damp_ocean = 20.0
damp_land = 10.0
hype_coef = 1.2e16      # Hyperdiffusion coefficient (m^4/s)

# --- Domain and Grid ---
d_lat = 2.5
lat_start = 15.0
lat_end = 85.0
d_lon = d_lat

# --- Local Optimization ---
loc_flag = 1
opt_region = "NEUAS"
opt_minlat = 30.
opt_maxlat = 70.
opt_minlon = 0.
opt_maxlon = 150.

# --- Constants ---
r = 6.37122e6       # Earth radius (m)
ome = 2*np.pi/86400 # Earth angular velocity (rad/s)
g = 9.8             # Gravity (m/s^2)

# --- Derived Parameters ---
lat = np.arange(lat_start, lat_end + d_lat, d_lat)
lon = np.arange(0.0, 360.0, d_lon)
X = np.size(lon)
Y = np.size(lat)
N = (Y - 4) * X     # State vector size (removes 2-point N/S boundary)

# Grid for interpolation
lat_grid, lon_grid = np.meshgrid(lat, lon, indexing='ij')

# Grid spacing in radians
dy = np.deg2rad(d_lat)
dx = np.deg2rad(d_lon)
