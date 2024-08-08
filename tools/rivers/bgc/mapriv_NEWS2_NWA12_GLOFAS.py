import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset

# name of netcdf file to be created
nc_file_name = './RiverNutrients_GlobalNEWS2_plusFe_Q100_GLOFAS_NWA12.nc'

# Parameters for the assignment algorithm
Q_min = 100  # minimum flow in m3 sec
plot_width = 3  # width of window (in degrees) for inspecting locations
min_dist = 2  # minimum distance (degrees) of the closest outflow point
inspect_map = 'n'  # flag enabling you to pause and inspect each river mapping as it is being done.

# Set the bio-availability of phosphorus and the fractionation of dissolved organic
frac_PP = 0.3
frac_ldon = 0.3
frac_sldon = 0.35
frac_srdon = 0.35
frac_ldop = 0.3
frac_sldop = 0.35
frac_srdop = 0.35
const_fed = 70.0e-6  # 40 nM dissolved iron concentration from De Baar and De Jong + 30nM Colloidal and nanoparticle flux

# GlobalNEWS2 data obtained from Emilio Mayorga
filename = 'GlobalNEWS2_RH2000Dataset-version1.0.xls'
basin = pd.read_excel(filename, sheet_name=1)
hydrology = pd.read_excel(filename, sheet_name=2)
loading = pd.read_excel(filename, sheet_name=3)

# Find all the river basins that empty into "land", e.g., lakes
ocean = basin['ocean']
land_index = ocean.apply(lambda x: 1 if x == 'Land' else 0)

river_names_all = basin['basinname']

# basin area in
area = basin['A']
lon_news_all = basin['mouth_lon']
lat_news_all = basin['mouth_lat']

# Loads in Mg yr-1 converted to moles per sec
DIN_load_all = loading['Ld_DIN'] * 1e6 / 14 / 86400 / 365
DIP_load_all = loading['Ld_DIP'] * 1e6 / 31 / 86400 / 365
DON_load_all = loading['Ld_DON'] * 1e6 / 14 / 86400 / 365
DOP_load_all = loading['Ld_DOP'] * 1e6 / 31 / 86400 / 365
Si_load_all = loading['Ld_DSi'] * 1e6 / 28.1 / 86400 / 365
PN_load_all = (loading['Ld_PN'] * 1e6 / 14 / 86400 / 365)
PP_load_all = (loading['Ld_PP'] * 1e6 / 31 / 86400 / 365) * frac_PP

# actual and natural discharge (convert from km3/yr to m3/sec)
# Used the actual hydrology to calculate concentrations
Qact_all = hydrology['Qact'] * 1e9 / (86400 * 365)
Qnat_all = hydrology['Qnat'] * 1e9 / (86400 * 365)


###########################################################################
# Filter for rivers in the region, set thresholds for minimum river size, #
# set parameters for plotting routines.                                   #
###########################################################################



import scipy.io as sio
import numpy as np
import xarray as xr

# Load runoff data from MATLAB .mat file
mat_data = sio.loadmat('./glofas_runoff_mean.mat')
runoff = mat_data['runoff']
lon_mod = mat_data['lon_mod']
lat_mod = mat_data['lat_mod']

# Convert runoff from kg m-2 sec-1 to m3 sec-1
area_mod = mat_data['area_mod']  # Assuming 'area_mod' is also loaded from the MATLAB file
Q_mod_monthly = np.zeros(runoff.shape)
for m in range(12):
    Q_mod_monthly[m, :, :] = runoff[m, :, :] * area_mod / 1000  # Conversion factor from kg to m3

Q_mod_ann = np.mean(Q_mod_monthly, axis=0)

# Clear variables to release memory
del runoff, area_mod

# Load grid data using xarray
grid_file = '/home/otel/Dropbox/trabalho_irado/Northeastern/202401_mom6_cobalt_compilation/run/tools_NWA25.COBALT/data/output/19930101.ocean_static.nc'
with xr.open_dataset(grid_file) as ds:
    depth = ds['deptho'].values.T  # Transpose to match MATLAB's permute function
    depth[np.isnan(depth)] = -1  # Replace NaN values with -1


###########################################################################
# Filter for rivers in the region, set thresholds for minimum river size, #
# set parameters for plotting routines.                                   #
###########################################################################



import numpy as np

# Use grid to filter rivers outside domain
lat_mod_max = np.max(lat_mod)
lat_mod_min = np.min(lat_mod)
lon_mod_max = np.max(lon_mod)
lon_mod_min = np.min(lon_mod)

in_region = np.where(
    (lon_news_all <= lon_mod_max) & (lon_news_all >= lon_mod_min) &
    (lat_news_all <= lat_mod_max) & (lat_news_all >= lat_mod_min) &
    (np.isfinite(Qact_all)) & (Qact_all > Q_min)
)[0]

# If you are using a high threshold, grab one smaller river to constrain
# Carribean Islands
if Q_min > 100:
    cuba_ind = np.where(river_names_all == 'GHAASBasin1808')[0]
    in_region = np.concatenate((in_region, cuba_ind))

num_rivers = len(in_region)

# Establish vectors of flow and nutrient loads for the NWA
Qact = Qact_all[in_region]
lon_news = lon_news_all[in_region]
lat_news = lat_news_all[in_region]
DIN_load = DIN_load_all[in_region]
DON_load = DON_load_all[in_region]
PN_load = PN_load_all[in_region]
DIP_load = DIP_load_all[in_region]
DOP_load = DOP_load_all[in_region]
PP_load = PP_load_all[in_region]
Si_load = Si_load_all[in_region]
river_names = river_names_all[in_region]

###########################################################################
# Following inspection of initial mapping, add any manual edits here to   #
# prevent anomalous extrapolations, etc.                                  #
###########################################################################

# see original script

###########################################################################
# END MANUAL EDITS                                                        #
###########################################################################

# Sort rivers by discharge
sort_ind = np.argsort(Qact.values)
Qact_sort = Qact.iloc[sort_ind]
lon_news_sort = lon_news.iloc[sort_ind]
lat_news_sort = lat_news.iloc[sort_ind]
DIN_load_sort = DIN_load.iloc[sort_ind]
DON_load_sort = DON_load.iloc[sort_ind]
PN_load_sort = PN_load.iloc[sort_ind]
DIP_load_sort = DIP_load.iloc[sort_ind]
DOP_load_sort = DOP_load.iloc[sort_ind]
PP_load_sort = PP_load.iloc[sort_ind]
Si_load_sort = Si_load.iloc[sort_ind]
river_names_sort = river_names.iloc[sort_ind]

# Total N and P load diagnostics
N_load_sort = DIN_load_sort + DON_load_sort + PN_load_sort
P_load_sort = DIP_load_sort + DOP_load_sort + PP_load_sort

# Calculate concentrations
# Loads are in moles N sec-1, Q in m3 s-1; conc in moles N m-3
DIN_conc_sort = DIN_load_sort / Qact_sort
DON_conc_sort = DON_load_sort / Qact_sort
DIP_conc_sort = DIP_load_sort / Qact_sort
DOP_conc_sort = DOP_load_sort / Qact_sort
PN_conc_sort = PN_load_sort / Qact_sort
PP_conc_sort = PP_load_sort / Qact_sort
Si_conc_sort = Si_load_sort / Qact_sort

# Initialize vectors to hold nutrient concentrations at each runoff point
aa = np.where(Q_mod_ann > 0)
Q_mod_vec = Q_mod_ann[aa]
din_conc_vec = np.zeros_like(Q_mod_vec)
don_conc_vec = np.zeros_like(Q_mod_vec)
pn_conc_vec = np.zeros_like(Q_mod_vec)
dip_conc_vec = np.zeros_like(Q_mod_vec)
dop_conc_vec = np.zeros_like(Q_mod_vec)
pp_conc_vec = np.zeros_like(Q_mod_vec)
si_conc_vec = np.zeros_like(Q_mod_vec)
lon_mod_runoff_vec = lon_mod[aa]
lat_mod_runoff_vec = lat_mod[aa]

###########################################################################
# Loop identifies points assigned to each river                           #
###########################################################################


from scipy.spatial.distance import cdist
from scipy.interpolate import NearestNDInterpolator
import matplotlib.pyplot as plt


for k in range(num_rivers):
    print(k + 1)
    dist = cdist([[lon_news_sort.iloc[k], lat_news_sort.iloc[k]]],
                 np.column_stack([lon_mod_runoff_vec, lat_mod_runoff_vec]))
    dist_sort_ind = np.argsort(dist.ravel())
    dist_sort = dist.ravel()[dist_sort_ind]
    
    if dist_sort[0] < min_dist:
        Q_sum1 = 0
        Q_sum2 = 0
        n = 0
        
        while Q_sum2 < Qact_sort.iloc[k]:
            Q_sum1 = Q_sum2
            n += 1
            Q_sum2 = Q_sum1 + Q_mod_vec[dist_sort_ind[n]]
        
        if abs(Q_sum1 - Qact_sort.iloc[k]) < abs(Q_sum2 - Qact_sort.iloc[k]):
            nrp = n - 1
            print(Q_sum1, Qact_sort.iloc[k])  # a quick check for comparable flow
        else:
            nrp = n
            print(Q_sum2, Qact_sort.iloc[k])  # a quick check for comparable flow
        
        # Update nutrient concentrations
        ind = dist_sort_ind[:nrp]
        din_conc_vec[ind] = DIN_conc_sort.iloc[k]
        don_conc_vec[ind] = DON_conc_sort.iloc[k]
        dip_conc_vec[ind] = DIP_conc_sort.iloc[k]
        dop_conc_vec[ind] = DOP_conc_sort.iloc[k]
        pn_conc_vec[ind] = PN_conc_sort.iloc[k]
        pp_conc_vec[ind] = PP_conc_sort.iloc[k]
        si_conc_vec[ind] = Si_conc_sort.iloc[k]
        
        if inspect_map == 'y':
            plt.figure()
            plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=3, c=np.log10(Q_mod_vec))
            plt.scatter(lon_mod_runoff_vec[ind], lat_mod_runoff_vec[ind], s=40, c=np.log10(Q_mod_vec[ind]), marker='o', edgecolors='k')
            plt.colorbar(label='log10(Q_mod_vec)')
            plt.plot(lon_news_sort.iloc[k], lat_news_sort.iloc[k], 'k.', markersize=20)
            plt.contour(lon_mod, lat_mod, depth, levels=[0], colors='k')
            plt.title(f'River number: {k+1}, Name: {river_names_sort.iloc[k]}')
            plt.xlabel('Longitude')
            plt.ylabel('Latitude')
            plt.show()
            
            # Provide diagnostics
            N_conc = DIN_conc_sort.iloc[k] + DON_conc_sort.iloc[k] + PN_conc_sort.iloc[k]
            P_conc = DIP_conc_sort.iloc[k] + DOP_conc_sort.iloc[k] + PP_conc_sort.iloc[k]
            print(f"River name: {river_names_sort.iloc[k]}")
            print(f"Coordinates: [{lon_news_sort.iloc[k]}, {lat_news_sort.iloc[k]}]")
            print(f"Total flow in m3 sec: {Qact_sort.iloc[k]}, {np.sum(Q_mod_vec[ind])}")
            print(f"N load in Gg per year: {N_load_sort.iloc[k]}, {np.sum(Q_mod_vec[ind]) * N_conc * 14 * 86400 * 365 / 1e9}")
            print(f"P load in Gg per year: {P_load_sort.iloc[k]}, {np.sum(Q_mod_vec[ind]) * P_conc * 31 * 86400 * 365 / 1e9}")
            print(f"N, P conc (mmoles m-3): [{DIN_conc_sort.iloc[k]}, {DON_conc_sort.iloc[k]}, {PN_conc_sort.iloc[k]}]")
            print(f"DIP, DOP, PP conc (mmoles m-3): [{DIP_conc_sort.iloc[k]}, {DOP_conc_sort.iloc[k]}, {PP_conc_sort.iloc[k]}]")
            print(f"Total N, Total P, Total N: Total P: {N_conc * 1e3}, {P_conc * 1e3}, {N_conc / P_conc}")
            print(f"DO:DI and P:DI ratios: {DON_conc_sort.iloc[k] / DIN_conc_sort.iloc[k]}, {PN_conc_sort.iloc[k] / DIN_conc_sort.iloc[k]}")
            print(f"DOP:DI and PP:DI ratios: {DOP_conc_sort.iloc[k] / DIP_conc_sort.iloc[k]}, {PP_conc_sort.iloc[k] / DIP_conc_sort.iloc[k]}")
            print(f"Silica concentration (mmoles m-3): {Si_conc_sort.iloc[k] * 1e3}")
            #input("Press Enter to continue...")
    else:
        if inspect_map == 'y':
            plt.figure()
            plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=3, c=np.log10(Q_mod_vec))
            plt.plot(lon_news_sort.iloc[k], lat_news_sort.iloc[k], 'k.', markersize=20)
            plt.colorbar(label='log10(Q_mod_vec)')
            plt.title(f'OUTSIDE: River number: {k+1}, Name: {river_names_sort.iloc[k]}')
            plt.xlabel('Longitude')
            plt.ylabel('Latitude')
            plt.show()


from scipy.interpolate import griddata

# Function to perform nearest neighbor interpolation
def fill_zeros(lon, lat, conc_vec):
    # Find indices where concentration values are zero
    zero_indices = np.where(conc_vec == 0)[0]
    non_zero_indices = np.where(conc_vec > 0)[0]
    
    # Perform nearest neighbor interpolation
    filled_values = griddata((lon[non_zero_indices], lat[non_zero_indices]), 
                             conc_vec[non_zero_indices], 
                             (lon[zero_indices], lat[zero_indices]), 
                             method='nearest')
    
    # Update concentration vector with filled values
    conc_vec[zero_indices] = filled_values
    
    return conc_vec

# Fill in zero values in concentration vectors
din_conc_vec = fill_zeros(lon_mod_runoff_vec, lat_mod_runoff_vec, din_conc_vec)
don_conc_vec = fill_zeros(lon_mod_runoff_vec, lat_mod_runoff_vec, don_conc_vec)
pn_conc_vec = fill_zeros(lon_mod_runoff_vec, lat_mod_runoff_vec, pn_conc_vec)
dip_conc_vec = fill_zeros(lon_mod_runoff_vec, lat_mod_runoff_vec, dip_conc_vec)
dop_conc_vec = fill_zeros(lon_mod_runoff_vec, lat_mod_runoff_vec, dop_conc_vec)
pp_conc_vec = fill_zeros(lon_mod_runoff_vec, lat_mod_runoff_vec, pp_conc_vec)
si_conc_vec = fill_zeros(lon_mod_runoff_vec, lat_mod_runoff_vec, si_conc_vec)

# Calculate total nutrient concentrations
totn_conc_vec = din_conc_vec + don_conc_vec + pn_conc_vec
totp_conc_vec = dip_conc_vec + dop_conc_vec + pp_conc_vec

# Calculate ratios of dissolved and particulate to inorganic
din_flux_vec = din_conc_vec * Q_mod_vec
dip_flux_vec = dip_conc_vec * Q_mod_vec
don_flux_vec = don_conc_vec * Q_mod_vec
dop_flux_vec = dop_conc_vec * Q_mod_vec
pn_flux_vec = pn_conc_vec * Q_mod_vec
pp_flux_vec = pp_conc_vec * Q_mod_vec

don_ratio = np.sum(don_flux_vec) / np.sum(din_flux_vec)
dop_ratio = np.sum(dop_flux_vec) / np.sum(dip_flux_vec)
pn_ratio = np.sum(pn_flux_vec) / np.sum(din_flux_vec)
pp_ratio = np.sum(pp_flux_vec) / np.sum(dip_flux_vec)


###########################################################################
# Produce plots to evaluate the mapping                                   #
###########################################################################



import numpy as np
import matplotlib.pyplot as plt

# Use total fluxes to scale dots, m3 sec-1 * moles m-3 = moles sec-1
totn_flux_vec = totn_conc_vec * Q_mod_vec
totp_flux_vec = totp_conc_vec * Q_mod_vec
totsi_flux_vec = si_conc_vec * Q_mod_vec

# Scale marker size with the total nitrogen flux
ms_vec = np.zeros_like(Q_mod_vec)
ms_vec[np.log10(Q_mod_vec) < 0] = 1
ms_vec[(np.log10(Q_mod_vec) > 0) & (np.log10(Q_mod_vec) < 1)] = 2.5
ms_vec[(np.log10(Q_mod_vec) > 1) & (np.log10(Q_mod_vec) < 2)] = 10
ms_vec[(np.log10(Q_mod_vec) > 2) & (np.log10(Q_mod_vec) < 3)] = 25
ms_vec[np.log10(Q_mod_vec) > 3] = 100

# Plotting
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

# Plot total nitrogen concentration
axs[0, 0].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=totn_conc_vec * 1e3, cmap='viridis', alpha=0.6)
axs[0, 0].set_title('Total Nitrogen Concentration (mmoles m$^{-3}$)')
axs[0, 0].set_xlabel('Longitude')
axs[0, 0].set_ylabel('Latitude')
axs[0, 0].grid(True)

# Plot total phosphorus concentration
axs[0, 1].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=totp_conc_vec * 1e3, cmap='viridis', alpha=0.6)
axs[0, 1].set_title('Total Phosphorus Concentration (mmoles m$^{-3}$)')
axs[0, 1].set_xlabel('Longitude')
axs[0, 1].set_ylabel('Latitude')
axs[0, 1].grid(True)

# Plot N:P ratio
axs[0, 2].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=totn_conc_vec / totp_conc_vec, cmap='viridis', alpha=0.6)
axs[0, 2].set_title('N:P Ratio')
axs[0, 2].set_xlabel('Longitude')
axs[0, 2].set_ylabel('Latitude')
axs[0, 2].grid(True)

# Plot DIN concentration
axs[1, 0].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=din_conc_vec * 1e3, cmap='viridis', alpha=0.6)
axs[1, 0].set_title('DIN Concentration (mmoles m$^{-3}$)')
axs[1, 0].set_xlabel('Longitude')
axs[1, 0].set_ylabel('Latitude')
axs[1, 0].grid(True)

# Plot DON concentration
axs[1, 1].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=don_conc_vec * 1e3, cmap='viridis', alpha=0.6)
axs[1, 1].set_title('DON Concentration (mmoles m$^{-3}$)')
axs[1, 1].set_xlabel('Longitude')
axs[1, 1].set_ylabel('Latitude')
axs[1, 1].grid(True)

# Plot PN concentration
axs[1, 2].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=pn_conc_vec * 1e3, cmap='viridis', alpha=0.6)
axs[1, 2].set_title('PN Concentration (mmoles m$^{-3}$)')
axs[1, 2].set_xlabel('Longitude')
axs[1, 2].set_ylabel('Latitude')
axs[1, 2].grid(True)

plt.tight_layout()
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Output total fluxes in gigagrams of N
total_annual_n = np.sum(totn_flux_vec) * 86400 * 365 * 14 / 1e12
total_annual_p = np.sum(totp_flux_vec) * 86400 * 365 * 31 / 1e12
total_annual_si = np.sum(totsi_flux_vec) * 86400 * 365 * 28.1 / 1e12

#input("Press any key to continue...")

# Initialize 2D concentration arrays
din_conc = np.zeros_like(lon_mod)
don_conc = np.zeros_like(lon_mod)
dip_conc = np.zeros_like(lon_mod)
dop_conc = np.zeros_like(lon_mod)
pn_conc = np.zeros_like(lon_mod)
pp_conc = np.zeros_like(lon_mod)
si_conc = np.zeros_like(lon_mod)

# Map concentration vectors onto 2D arrays
aa = np.where(Q_mod_ann > 0)
din_conc[aa] = din_conc_vec
don_conc[aa] = don_conc_vec
pn_conc[aa] = pn_conc_vec
dip_conc[aa] = dip_conc_vec
dop_conc[aa] = dop_conc_vec
pp_conc[aa] = pp_conc_vec
si_conc[aa] = si_conc_vec

NO3_CONC = din_conc.copy()
LDON_CONC = frac_ldon * don_conc.copy()
SLDON_CONC = frac_sldon * don_conc.copy()
SRDON_CONC = frac_srdon * don_conc.copy()
PO4_CONC = dip_conc.copy()
LDOP_CONC = frac_ldop * dop_conc.copy()
SLDOP_CONC = frac_sldop * dop_conc.copy()
SRDOP_CONC = frac_srdop * dop_conc.copy()
NDET_CONC = pn_conc.copy()
PDET_CONC = pp_conc.copy()
SI_CONC = si_conc.copy()

# Add iron concentrations
FED_CONC = NO3_CONC.copy()
FEDET_CONC = NO3_CONC.copy()
FED_CONC[FED_CONC > 0] = const_fed
FEDET_CONC[FEDET_CONC > 0] = 0.0

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define marker size
ms = 8
plt.close('all')
# Create figure 7
fig7 = plt.figure(7)
fig7.clf()

# Subplot 1: log10(NO3 CONC)
ax1 = fig7.add_subplot(3, 2, 1, projection='3d')
ax1.set_title('log10(NO3 CONC)')
sc = ax1.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(NO3_CONC.flatten()), s=ms, c=np.log10(NO3_CONC.flatten()), cmap='viridis')
ax1.set_zlabel('log10(NO3 CONC)')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')
ax1.set_zlim([-4, -1])
plt.colorbar(sc, ax=ax1)

# Subplot 2: log10(LDON CONC)
ax2 = fig7.add_subplot(3, 2, 2, projection='3d')
ax2.set_title('log10(LDON CONC)')
sc = ax2.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(LDON_CONC.flatten()), s=ms, c=np.log10(LDON_CONC.flatten()), cmap='viridis')
ax2.set_zlabel('log10(LDON CONC)')
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Latitude')
ax2.set_zlim([-4, -1])
plt.colorbar(sc, ax=ax2)

# Subplot 3: log10(SLDON CONC)
ax3 = fig7.add_subplot(3, 2, 3, projection='3d')
ax3.set_title('log10(SLDON CONC)')
sc = ax3.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SLDON_CONC.flatten()), s=ms, c=np.log10(SLDON_CONC.flatten()), cmap='viridis')
ax3.set_zlabel('log10(SLDON CONC)')
ax3.set_xlabel('Longitude')
ax3.set_ylabel('Latitude')
ax3.set_zlim([-4, -1])
plt.colorbar(sc, ax=ax3)

# Subplot 4: log10(SRDON CONC)
ax4 = fig7.add_subplot(3, 2, 4, projection='3d')
ax4.set_title('log10(SRDON CONC)')
sc = ax4.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SRDON_CONC.flatten()), s=ms, c=np.log10(SRDON_CONC.flatten()), cmap='viridis')
ax4.set_zlabel('log10(SRDON CONC)')
ax4.set_xlabel('Longitude')
ax4.set_ylabel('Latitude')
ax4.set_zlim([-4, -1])
plt.colorbar(sc, ax=ax4)

# Subplot 5: log10(NDET CONC)
ax5 = fig7.add_subplot(3, 2, 5, projection='3d')
ax5.set_title('log10(NDET CONC)')
sc = ax5.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(NDET_CONC.flatten()), s=ms, c=np.log10(NDET_CONC.flatten()), cmap='viridis')
ax5.set_zlabel('log10(NDET CONC)')
ax5.set_xlabel('Longitude')
ax5.set_ylabel('Latitude')
ax5.set_zlim([-4, -1])
plt.colorbar(sc, ax=ax5)

plt.show()

# Create figure 8
fig8 = plt.figure(8)
fig8.clf()

# Subplot 1: log10(PO4 CONC)
ax6 = fig8.add_subplot(3, 2, 1, projection='3d')
sc6 = ax6.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(PO4_CONC.flatten()), s=ms, c=np.log10(PO4_CONC.flatten()), cmap='viridis')
ax6.set_title('log10(PO4 CONC)')
ax6.set_xlabel('Longitude')
ax6.set_ylabel('Latitude')
ax6.set_zlabel('log10(PO4 CONC)')
ax6.set_zlim([-4, -2])
fig8.colorbar(sc6, ax=ax6)

# Subplot 2: log10(LDOP CONC)
ax7 = fig8.add_subplot(3, 2, 2, projection='3d')
sc7 = ax7.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(LDOP_CONC.flatten()), s=ms, c=np.log10(LDOP_CONC.flatten()), cmap='viridis')
ax7.set_title('log10(LDOP CONC)')
ax7.set_xlabel('Longitude')
ax7.set_ylabel('Latitude')
ax7.set_zlabel('log10(LDOP CONC)')
ax7.set_zlim([-4, -2])
fig8.colorbar(sc7, ax=ax7)

# Subplot 3: log10(SLDOP CONC)
ax8 = fig8.add_subplot(3, 2, 3, projection='3d')
sc8 = ax8.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SLDOP_CONC.flatten()), s=ms, c=np.log10(SLDOP_CONC.flatten()), cmap='viridis')
ax8.set_title('log10(SLDOP CONC)')
ax8.set_xlabel('Longitude')
ax8.set_ylabel('Latitude')
ax8.set_zlabel('log10(SLDOP CONC)')
ax8.set_zlim([-4, -2])
fig8.colorbar(sc8, ax=ax8)

# Subplot 4: log10(SRDOP CONC)
ax9 = fig8.add_subplot(3, 2, 4, projection='3d')
sc9 = ax9.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SRDOP_CONC.flatten()), s=ms, c=np.log10(SRDOP_CONC.flatten()), cmap='viridis')
ax9.set_title('log10(SRDOP CONC)')
ax9.set_xlabel('Longitude')
ax9.set_ylabel('Latitude')
ax9.set_zlabel('log10(SRDOP CONC)')
ax9.set_zlim([-4, -2])
fig8.colorbar(sc9, ax=ax9)

# Subplot 5: log10(PDET CONC)
ax10 = fig8.add_subplot(3, 2, 5, projection='3d')
sc10 = ax10.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(PDET_CONC.flatten()), s=ms, c=np.log10(PDET_CONC.flatten()), cmap='viridis')
ax10.set_title('log10(PDET CONC)')
ax10.set_xlabel('Longitude')
ax10.set_ylabel('Latitude')
ax10.set_zlabel('log10(PDET CONC)')
ax10.set_zlim([-4, -2])
fig8.colorbar(sc10, ax=ax10)

plt.show()

# Create figure 9
fig9 = plt.figure(9)
fig9.clf()

# Subplot 1: log10(FED CONC)
ax11 = fig9.add_subplot(3, 2, 1, projection='3d')
sc11 = ax11.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(FED_CONC.flatten()), s=ms, c=np.log10(FED_CONC.flatten()), cmap='viridis')
ax11.set_title('log10(FED CONC)')
ax11.set_xlabel('Longitude')
ax11.set_ylabel('Latitude')
ax11.set_zlabel('log10(FED CONC)')
ax11.set_zlim([-5, -3])
fig9.colorbar(sc11, ax=ax11)

# Subplot 2: log10(FEDET CONC)
ax12 = fig9.add_subplot(3, 2, 2, projection='3d')
sc12 = ax12.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(FEDET_CONC.flatten()), s=ms, c=np.log10(FEDET_CONC.flatten()), cmap='viridis')
ax12.set_title('log10(FEDET CONC)')
ax12.set_xlabel('Longitude')
ax12.set_ylabel('Latitude')
ax12.set_zlabel('log10(FEDET CONC)')
ax12.set_zlim([-5, -3])
fig9.colorbar(sc12, ax=ax12)

# Subplot 3: log10(SI CONC)
ax13 = fig9.add_subplot(3, 2, 3, projection='3d')
sc13 = ax13.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SI_CONC.flatten()), s=ms, c=np.log10(SI_CONC.flatten()), cmap='viridis')
ax13.set_title('log10(SI CONC)')
ax13.set_xlabel('Longitude')
ax13.set_ylabel('Latitude')
ax13.set_zlabel('log10(SI CONC)')
ax13.set_zlim([-3, 0])
fig9.colorbar(sc13, ax=ax13)

plt.show()

###########################################################################
# Save Files                                                              #
###########################################################################

import xarray as xr

# Define dimensions
time = 0
nlat = lat_mod.shape[0]
nlon = lat_mod.shape[1]

# Create xarray Dataset
ds = xr.Dataset()

# Create dimensions
ds['time'] = xr.DataArray([0], dims='time')
ds['y'] = xr.DataArray(np.arange(nlat), dims='y')
ds['x'] = xr.DataArray(np.arange(nlon), dims='x')

# Create variables
ds['lat'] = xr.DataArray(lat_mod, dims=('y', 'x'))
ds['lat'].attrs['units'] = 'degrees north'

ds['lon'] = xr.DataArray(lon_mod, dims=('y', 'x'))
ds['lon'].attrs['units'] = 'degrees east'

ds['NO3_CONC']   = xr.DataArray([NO3_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'DIN_CONC'})
ds['LDON_CONC']  = xr.DataArray([LDON_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.3*DON_CONC'})
ds['SLDON_CONC'] = xr.DataArray([SLDON_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.35*DON_CONC'})
ds['SRDON_CONC'] = xr.DataArray([SRDON_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.35*DON_CONC'})
ds['NDET_CONC']  = xr.DataArray([NDET_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '1.0*PN_CONC'})
ds['PO4_CONC']   = xr.DataArray([PO4_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'PO4_CONC'})
ds['LDOP_CONC']  = xr.DataArray([LDOP_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.3*DOP_CONC'})
ds['SLDOP_CONC'] = xr.DataArray([SLDOP_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.35*DOP_CONC'})
ds['SRDOP_CONC'] = xr.DataArray([SRDOP_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.35*DOP_CONC'})
ds['PDET_CONC']  = xr.DataArray([PDET_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.3*PP_CONC'})
ds['FED_CONC']   = xr.DataArray([FED_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'FED_CONC'})
ds['FEDET_CONC'] = xr.DataArray([FEDET_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'FEDET_CONC'})
ds['SI_CONC']    = xr.DataArray([SI_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'SI_CONC'})


import cftime
t = cftime.datetime(1993, 1, 1, calendar='365_day')
ds = ds.assign_coords(time=[t])
# Save to netCDF file

for d in ds:
    print(d)
    ds[d].values#*1e-6

ds.to_netcdf(nc_file_name,
        format='NETCDF4',
        engine='netcdf4',
        unlimited_dims='time'
)

