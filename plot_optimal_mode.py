# =============================================================================
# Plots the time evolution of the leading optimal modes for the barotropic model.
# Each mode's evolution is saved as a separate page in a single PDF file.
#
# Creator: Zhenyu You, Georgia Tech
# Contact: zyou2998@gmail.com
# =============================================================================

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.backends.backend_pdf import PdfPages
import xarray as xr
from model_parameters import *
from cartopy.util import add_cyclic_point
import os

# --- 1. Setup ---
n_mode = 10  # Number of optimal modes to plot
n_time = 12  # Number of time steps to plot per mode
dt = 2       # Time interval between plots (days)

# Construct file paths
damp_str = f"damp_o{damp_ocean}_l{damp_land}"
if loc_flag == 0:
    base_name = f"optimal_{season}_tau{tau}_{damp_str}_{d_lat}_{lat_start}_{lat_end}"
else:
    base_name = f"optimal_{season}_local_{opt_region}_tau{tau}_{damp_str}_{d_lat}_{lat_start}_{lat_end}"

in_file = out_dir + base_name + ".nc"
in_u_file = in_dir + f"U_{d_lat}_{lat_start}_{lat_end}_{season}.nc"
out_pic = pic_dir + base_name + ".pdf"

# --- 2. Load Data ---
if ((n_time - 1) * dt > evo):
    print(f"Error: Requested plot time ({n_time*dt} days) exceeds total evolution time ({evo} days).")
    exit()

print(f"Loading optimal modes from: {os.path.basename(in_file)}")
with xr.open_dataset(in_file) as ds:
    data = ds["data"].values
    rate = ds["rate"].values
with xr.open_dataset(in_u_file) as ds:
    u_background = ds['u'][2:Y-2, :].values

# --- 3. Plot Modes and Save to PDF ---
cmap = plt.get_cmap('RdBu_r')
lat_interior = lat[2:Y-2]
u_cyclic, lon_cyclic_u = add_cyclic_point(u_background, coord=lon)
u_levels = np.arange(10, 35, 5)

with PdfPages(out_pic) as pdf:
    print(f"Generating PDF with {n_mode} modes...")
    for ii in range(n_mode):
        print(f"  - Plotting Mode {ii + 1}/{n_mode}")
        amp = rate[ii]
        mode_evolution = data[ii]
        
        fig = plt.figure(figsize=(10,15), dpi=300)
        
        for t in range(n_time):
            z = mode_evolution[int(t * dt / td)].reshape((Y - 4, X))
            z_cyclic, lon_cyclic = add_cyclic_point(z, coord=lon)
            
            # Set color scale individually for each panel
            max_s = 0.8 * np.max(np.abs(z)) + 1e-15
            shading_levels = np.linspace(-max_s, max_s, 11)
            
            # Create subplot
            num_rows = (n_time + 2) // 3
            ax = plt.subplot(num_rows, 3, t + 1, projection=ccrs.NorthPolarStereo(central_longitude=100))
            ax.coastlines()
            cf = ax.contourf(lon_cyclic, lat_interior, z_cyclic,
                             levels=shading_levels, cmap=cmap, extend='both',
                             transform=ccrs.PlateCarree())
            cs = ax.contour(lon_cyclic_u, lat_interior, u_cyclic,
                            levels=u_levels, colors='black', linewidths=1.0,
                            transform=ccrs.PlateCarree())
            
            ax.clabel(cs, fmt='%d', fontsize=8)
            ax.set_title(f'Day {t * dt}', fontsize=15)
            plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.03, shrink=1.0, format="%.1e")

        plt.tight_layout()
        fig.suptitle(f'Mode {ii + 1} amp factor: {amp:.4f}', fontsize=24, y=1.02)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

print(f"\nPDF generated successfully: {out_pic}")