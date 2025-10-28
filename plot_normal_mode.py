# =============================================================================
# Plots the spatial structure of the leading normal modes for the barotropic model.
# Each mode is saved as a separate page in a single PDF file.
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
n_mode = 50  # Number of modes to plot

# Construct file paths
damp_str = f"damp_o{damp_ocean}_l{damp_land}"
in_file = out_dir + f"normal_{season}_{damp_str}_{d_lat}_{lat_start}_{lat_end}.nc"
out_pic = pic_dir + f"normal_{season}_{damp_str}_{d_lat}_{lat_start}_{lat_end}.pdf"

# --- 2. Load Data ---
print(f"Loading normal modes from: {os.path.basename(in_file)}")
with xr.open_dataset(in_file) as ds:
    modes = ds["mode_real"].values
    growth_rates = ds["rate_real"].values

# Convert growth rate to e-folding time
e_folding_time = 1.0 / (growth_rates * 86400.0 + 1e-15)

# --- 3. Plot Modes and Save to PDF ---
cmap = plt.get_cmap('RdBu_r')
lat_interior = lat[2:Y-2]

with PdfPages(out_pic) as pdf:
    print(f"Generating PDF with {n_mode} modes...")
    for ii in range(n_mode):
        print(f"  - Plotting Mode {ii + 1}/{n_mode}")
        
        # Reshape the 1D mode vector into a 2D field
        mode = modes[:, ii].reshape((Y - 4, X))
        mode_cyclic, lon_cyclic = add_cyclic_point(mode, coord=lon)

        max_amp = np.max(np.abs(mode)) + 1e-15
        levels = np.linspace(-0.8 * max_amp, 0.8 * max_amp, 11)
        
        # Create figure
        fig = plt.figure(figsize=(10, 10), dpi=200)
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-100))
        ax.set_extent([0, 359.99, lat_start, lat_end], crs=ccrs.PlateCarree())
        
        cf = ax.contourf(lon_cyclic, lat_interior, mode_cyclic,
                         levels=levels, cmap=cmap, extend='both',
                         transform=ccrs.PlateCarree())
        
        ax.add_feature(cfeature.COASTLINE)
        ax.gridlines(linestyle='--', alpha=0.5)
        
        plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.05, shrink=0.7, label="Streamfunction Anomaly")
        ax.set_title(f'Mode {ii + 1}: e-folding time = {e_folding_time[ii]:.2f} days', fontsize=15)
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

print(f"\nPDF generated successfully: {out_pic}")