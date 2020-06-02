import matplotlib as plt
import pandas as pd
import matplotlib.pyplot as plt
import glob
from natsort import natsorted
from scipy import interpolate
import numpy as np
import seaborn as sns
from scipy.interpolate import griddata
from matplotlib import ticker, cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Useful links:
# https://stackoverflow.com/questions/42504987/python-interpolate-point-value-on-2d-grid/42505340
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html

# Get an ordered list of the luminosity tracks
lum_tracks = natsorted(glob.glob("Lum_m*.txt"))
print(lum_tracks)

fig = plt.figure(tight_layout=True, constrained_layout=False, figsize=[10, 9])
gs = fig.add_gridspec(nrows=1, ncols=1, wspace=0.15, hspace=0.15, left=0.07, right=0.999, top=0.99, bottom=0.08)
sns.set(style="ticks")
ax1 = fig.add_subplot(gs[0])


# Define data arrays for interpolation later on
xy_age_mass = []
z_lum       = []

# Fill the arrays
for lum_track in lum_tracks:

    # Read the specific data file
    luminosity_df = pd.read_csv(lum_track)

    # Cut the string to get the mass of the star
    star_mass = float(lum_track[-7:-4])

    # Read out age and luminosity
    age_list         = list(luminosity_df["age"]*1e+3)
    luminosity_list  = list(luminosity_df["lum"])

    mass_list        = np.ones(len(age_list))*star_mass

    # Fill the arrays
    zip_list = list(zip(age_list, mass_list))
    xy_age_mass.extend(zip_list)
    z_lum.extend(luminosity_list)


# Bring arrays in numpy shape
xy_age_mass = np.array(xy_age_mass)
z_lum = np.array(z_lum)

# Define the interpolation grids
grid_x, grid_y = np.mgrid[0.7:10000:100j, 0.1:1.4:100j]

# Interpolate
grid_z1 = griddata(xy_age_mass, z_lum, (grid_x, grid_y), method='linear', rescale=True)

# ---> CALCULATE SPECIFIC POINT
SPECIFIC_LUMINOSITY = griddata(xy_age_mass, z_lum, (100, 0.99), method='linear', rescale=True)
print(SPECIFIC_LUMINOSITY)


# Define color range for plotting
color_range = np.logspace(-4, 1, 255)

# Plot as contourf, becomes smoother the more often contourf is repeated
for i in range(0, 3):
    cont_linear = ax1.contourf(grid_x, grid_y, grid_z1, color_range, locator=ticker.LogLocator(), origin='lower', cmap="magma")

x_min = 0.7
x_max = 10000
y_min = 0.1
y_max = 1.4

ax1.set_xlim([x_min, x_max])
ax1.set_ylim([y_min, y_max])

ax1.set_xscale("log")

v_ticks  = [1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1]
# v_ticklabels = [ str(round(i,2)) for i in v_ticks ]

cbar = fig.colorbar(cont_linear, ax=ax1, shrink=0.9, ticks=v_ticks)
cbar.set_label(r'Bolometric luminosity ($L_{\star}/L_\odot$)')
# cbar.ax.tick_params(labelsize=12, labelrotation=0)
# cbar.ax.set_yticklabels(v_ticklabels)

plt.ylabel(r"Star mass ($M_{\star}/M_{\odot}$)")
plt.xlabel("Time (Myr)")

sns.despine()

plt.savefig("luminosity_grid.pdf")
