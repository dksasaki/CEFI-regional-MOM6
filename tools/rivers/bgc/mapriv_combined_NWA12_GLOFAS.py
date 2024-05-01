# Routine to map USGS nutrient data onto the MOM6 Northwest Atlantic (NWA) 
# grid. Run on matlab97 or above.

import numpy as np
import xarray as xr
import xarray_decorators as xrd
# name of netcdf file to be created
nc_file_name = './RiverNutrients_Integrated_NWA12_GLOFAS_RC4US1990to2022_2023_04_v2.nc'

# GLOBAL NEWS based map for filling in gaps
NEWS_file = './RiverNutrients_GlobalNEWS2_plusFe_Q100_GLOFAS_NWA12.nc'

# load in monthly world ocean T, S climatology for saturated oxygen calculation
woa_temp = xr.open_mfdataset('/home/otel/Dropbox/trabalho_irado/Northeastern/data/woa18/nwa25_regional/woa18_decav_t*_04.nc',
                           decode_times=False,
                           concat_dim='time',
                           combine='nested')


# Parameters for the assignment algorithm.
Q_min = 0  # minimum flow in m3 sec
plot_width = 4  # width of window (in degrees) for inspecting locations
min_dist = 2.0  # minimum distance (degrees) of the closest outflow point
max_dist = 2.0  # maximum distance (degrees) away that the algorithm looks for points for rivers that are in the domain
nutrient_option = 2  # option for deriving dissolved organic nutrients
inspect_map = 'n'  # flag enabling you to pause and inspect each river mapping as it is being done.
plot_progress = False

# set the bio-availability of phosphorus and the fractionation of dissolved organic
# PP is set to 30# based on Froelich; Partitioning of detritus between
frac_PP = 0.5
frac_ldon = 0.3
frac_sldon = 0.35
frac_srdon = 0.35
frac_ldop = 0.3
frac_sldop = 0.35
frac_srdop = 0.35
# 40 nM dissolved iron concentration from De Baar and De Jong + 30nM
# Colloidal and nanoparticle flux as reported in Canfield and Raiswell
const_fed = 70.0e-6

###########################################################################
# USGS data compiled by Fabian Gomez                                      #
###########################################################################


filename_chem = '/home/otel/Dropbox/trabalho_irado/Northeastern/data/RC4USCoast/0260455/3.3/data/0-data/mclim_19902022_chem.nc'


with xr.open_dataset(filename_chem) as ds:
    alk_monthly_RC4US = np.transpose(ds['alk'].values, (1, 0))
    dic_monthly_RC4US = np.transpose(ds['dic'].values, (1, 0))
    no3_monthly_RC4US = np.transpose(ds['no3'].values, (1, 0))
    nh4_monthly_RC4US = np.transpose(ds['nh4'].values, (1, 0))
    din_monthly_RC4US = no3_monthly_RC4US + nh4_monthly_RC4US
    dip_monthly_RC4US = np.transpose(ds['po4'].values, (1, 0))
    si_monthly_RC4US = np.transpose(ds['sio2'].values, (1, 0))

    o2_monthly_RC4US = np.transpose(ds['do'].values, (1, 0)) / 2.0
    don_monthly_RC4US = np.transpose(ds['don'].values, (1, 0))
    temp_monthly_RC4US = np.transpose(ds['temp'].values, (1, 0))

    if nutrient_option == 1:
        pn_monthly_RC4US = np.full_like(no3_monthly_RC4US, np.nan)
        dop_monthly_RC4US = np.full_like(no3_monthly_RC4US, np.nan)
        pp_monthly_RC4US = np.full_like(no3_monthly_RC4US, np.nan)
    elif nutrient_option == 2:
        tnf_monthly_RC4US = np.transpose(ds['tnf'].values, (1, 0))
        don2_monthly_RC4US = tnf_monthly_RC4US - din_monthly_RC4US
        tnu_monthly_RC4US = np.transpose(ds['tnu'].values, (1, 0))
        
        pn_monthly_RC4US = tnu_monthly_RC4US - tnf_monthly_RC4US
        pn_monthly_RC4US[pn_monthly_RC4US < 0] = np.nan
        
        tpf_monthly_RC4US = np.transpose(ds['tpf'].values, (1, 0))
        dop_monthly_RC4US = tpf_monthly_RC4US - dip_monthly_RC4US
        dop_monthly_RC4US[dop_monthly_RC4US < 0] = np.nan
        
        tpu_monthly_RC4US = np.transpose(ds['tpu'].values, (1, 0))
        pp_monthly_RC4US = (tpu_monthly_RC4US - tpf_monthly_RC4US) * frac_PP
        pp_monthly_RC4US[pp_monthly_RC4US < 0] = np.nan

dfe_monthly_RC4US = np.full_like(no3_monthly_RC4US, np.nan)
pfe_monthly_RC4US = np.full_like(no3_monthly_RC4US, np.nan)



filename_discharge = '/home/otel/Dropbox/trabalho_irado/Northeastern/data/RC4USCoast/0260455/3.3/data/0-data/mclim_19902022_disc.nc'

with xr.open_dataset(filename_discharge) as ds:
    Q_monthly_RC4US = np.transpose(ds['disc'].values, (1, 0))  # m^-3 sec^-1
    station_names_RC4US = ds['river_name'].values
    lon_stations_RC4US = ds['mouth_lon'].values
    lat_stations_RC4US = ds['mouth_lat'].values

Q_ann_RC4US = np.nanmean(Q_monthly_RC4US, axis=1)
dic_ann_RC4US = np.nanmean(dic_monthly_RC4US, axis=1)
alk_ann_RC4US = np.nanmean(alk_monthly_RC4US, axis=1)
no3_ann_RC4US = np.nanmean(no3_monthly_RC4US, axis=1)
nh4_ann_RC4US = np.nanmean(nh4_monthly_RC4US, axis=1)
o2_ann_RC4US = np.nanmean(o2_monthly_RC4US, axis=1)
dip_ann_RC4US = np.nanmean(dip_monthly_RC4US, axis=1)
si_ann_RC4US = np.nanmean(si_monthly_RC4US, axis=1)
din_ann_RC4US = no3_ann_RC4US + nh4_ann_RC4US
don_ann_RC4US = np.nanmean(don_monthly_RC4US, axis=1)

if nutrient_option == 1:
    don_ann_RC4US = np.full_like(lon_stations_RC4US, np.nan)
    pn_ann_RC4US = np.full_like(lon_stations_RC4US, np.nan)
    dop_ann_RC4US = np.full_like(lon_stations_RC4US, np.nan)
    pp_ann_RC4US = np.full_like(lon_stations_RC4US, np.nan)
elif nutrient_option == 2:
    don2_ann_RC4US = np.nanmean(don2_monthly_RC4US, axis=1)
    pn_ann_RC4US = np.nanmean(pn_monthly_RC4US, axis=1)
    dop_ann_RC4US = np.nanmean(dop_monthly_RC4US, axis=1)
    pp_ann_RC4US = np.nanmean(pp_monthly_RC4US, axis=1)

dfe_ann_RC4US = np.full_like(lon_stations_RC4US, np.nan)
pfe_ann_RC4US = np.full_like(lon_stations_RC4US, np.nan)

for n in range(len(lon_stations_RC4US)):
    if station_names_RC4US[n] == 'Susquehanna':
        lat_stations_RC4US[n] = 38.5
        lon_stations_RC4US[n] = -77.5
    elif station_names_RC4US[n] == 'Delaware':
        lat_stations_RC4US[n] = 39.5
        lon_stations_RC4US[n] = -75.5
    elif station_names_RC4US[n] == 'Potomac':
        lat_stations_RC4US[n] = 38.5
        lon_stations_RC4US[n] = -77.5
    elif station_names_RC4US[n] == 'Mississippi':
        lat_stations_RC4US[n] = 29.25
        lon_stations_RC4US[n] = -89.25
    elif station_names_RC4US[n] == 'Alabama':
        lat_stations_RC4US[n] = 30.5


# ###########################################################################
# # Gulf of Saint Lawrence Data from Diane Lavoie                           #
# ###########################################################################

# import pandas as pd

# # Load river data
# dic_GSL = np.loadtxt("/archive/cas/COBALT_EVAL/River_Data/Lavoie_GSL/river_DIC.dat")
# alk_GSL = np.loadtxt("/archive/cas/COBALT_EVAL/River_Data/Lavoie_GSL/river_Alk.dat")
# no3_GSL = np.loadtxt("/archive/cas/COBALT_EVAL/River_Data/Lavoie_GSL/river_nox.dat")

# # Load station data from CSV file
# file_name_GSL = '/archive/cas/COBALT_EVAL/River_Data/Lavoie_GSL/river_discharge_approx.csv'
# station_data_GSL = pd.read_csv(file_name_GSL)

# station_names_GSL = station_data_GSL.iloc[:, 3].values
# lon_stations_GSL = station_data_GSL.iloc[:, 1].values
# lat_stations_GSL = station_data_GSL.iloc[:, 2].values
# Q_ann_GSL = station_data_GSL.iloc[:, 4].values

# # Create arrays of monthly values for all river constituents.
# dic_monthly_GSL = dic_GSL[0:12, 2:]
# alk_monthly_GSL = alk_GSL[0:12, 2:]
# no3_monthly_GSL = no3_GSL[0:12, 2:]

# # Set NH4 to zero where not available
# nh4_monthly_GSL = np.zeros((42, lon_stations_GSL.shape[0]))

# # Calculate DIN
# din_monthly_GSL = no3_monthly_GSL + nh4_monthly_GSL

# # Initialize other monthly variables with NaN
# don_monthly_GSL = np.full((12, lon_stations_GSL.shape[0]), np.nan)
# pn_monthly_GSL = np.full((12, lon_stations_GSL.shape[0]), np.nan)
# dip_monthly_GSL = np.full((12, lon_stations_GSL.shape[0]), np.nan)
# dop_monthly_GSL = np.full((12, lon_stations_GSL.shape[0]), np.nan)
# pp_monthly_GSL = np.full((12, lon_stations_GSL.shape[0]), np.nan)
# dfe_monthly_GSL = np.full((12, lon_stations_GSL.shape[0]), np.nan)
# pfe_monthly_GSL = np.full((12, lon_stations_GSL.shape[0]), np.nan)
# si_monthly_GSL = np.full((12, lon_stations_GSL.shape[0]), np.nan)
# o2_monthly_GSL = np.full((12, lon_stations_GSL.shape[0]), np.nan)

# # Calculate annual means
# dic_ann_GSL = np.mean(dic_monthly_GSL, axis=0)
# alk_ann_GSL = np.mean(alk_monthly_GSL, axis=0)
# no3_ann_GSL = np.mean(no3_monthly_GSL, axis=0)
# nh4_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# din_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# don_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# pn_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# dip_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# dop_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# pp_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# dfe_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# pfe_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# si_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)
# o2_ann_GSL = np.full(lon_stations_GSL.shape[0], np.nan)

# # Update coordinates
# for n in range(lon_stations_GSL.shape[0]):
#     if station_names_GSL[n] == 'Saint John':
#         lat_stations_GSL[n] = 46.25
#         lon_stations_GSL[n] = -66.25



# ###########################################################################
# # Selected Large Rivers with Global NEWS nutrients + DIC/Alk (i.e.,       #
# # "NEWSplus").  These rivers fill in any large areas without direct obs   #
# # to prevent unrealistic extrapolation from much different areas.  For    #
# # NWA12, these are primary in Mexico, Central America and South America.  #
# ###########################################################################


# import pandas as pd
# import numpy as np

# # Read station data from CSV file
# file_name_extra = '/archive/cas/Regional_MOM6/NWA12/NWA12_ExtraRivers/stations_extra_NWA12.csv'
# station_data_extra = pd.read_csv(file_name_extra)

# # Extract data from DataFrame columns
# station_names_extra = station_data_extra.iloc[:, 0].values
# lon_stations_extra = station_data_extra.iloc[:, 2].values
# lat_stations_extra = station_data_extra.iloc[:, 3].values
# Q_ann_extra = station_data_extra.iloc[:, 4].values

# dic_ann_extra = station_data_extra.iloc[:, 5].values
# alk_ann_extra = station_data_extra.iloc[:, 6].values
# no3_ann_extra = station_data_extra.iloc[:, 7].values
# nh4_ann_extra = np.zeros_like(lon_stations_extra)
# din_ann_extra = station_data_extra.iloc[:, 7].values
# don_ann_extra = station_data_extra.iloc[:, 8].values
# pn_ann_extra = station_data_extra.iloc[:, 9].values
# dip_ann_extra = station_data_extra.iloc[:, 10].values
# dop_ann_extra = station_data_extra.iloc[:, 11].values
# pp_ann_extra = station_data_extra.iloc[:, 12].values
# dfe_ann_extra = np.full_like(lon_stations_extra, np.nan)
# pfe_ann_extra = np.full_like(lon_stations_extra, np.nan)
# si_ann_extra = station_data_extra.iloc[:, 13].values
# o2_ann_extra = np.full_like(lon_stations_extra, np.nan)

# # Initialize monthly variables with NaN
# dic_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# alk_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# no3_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# nh4_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# din_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# don_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# pn_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# dip_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# dop_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# pp_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# dfe_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# pfe_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# si_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)
# o2_monthly_extra = np.full((12, lon_stations_extra.shape[0]), np.nan)


import numpy as np

# Combine all annual station data
station_names_all = np.concatenate((station_names_RC4US,)) # station_names_GSL, station_names_extra))
lon_stations_all = np.concatenate((lon_stations_RC4US,)) # lon_stations_GSL, lon_stations_extra))
lat_stations_all = np.concatenate((lat_stations_RC4US,)) # lat_stations_GSL, lat_stations_extra))
Q_ann_all = np.concatenate((Q_ann_RC4US,)) # Q_ann_GSL, Q_ann_extra))

dic_ann_all = np.concatenate((dic_ann_RC4US,)) # dic_ann_GSL, dic_ann_extra))
alk_ann_all = np.concatenate((alk_ann_RC4US,)) # alk_ann_GSL, alk_ann_extra))
no3_ann_all = np.concatenate((no3_ann_RC4US,)) # no3_ann_GSL, no3_ann_extra))
nh4_ann_all = np.concatenate((nh4_ann_RC4US,)) # nh4_ann_GSL, nh4_ann_extra))
din_ann_all = np.concatenate((din_ann_RC4US,)) # din_ann_GSL, din_ann_extra))
don_ann_all = np.concatenate((don_ann_RC4US,)) # don_ann_GSL, don_ann_extra))
pn_ann_all = np.concatenate((pn_ann_RC4US,)) # pn_ann_GSL, pn_ann_extra))
dip_ann_all = np.concatenate((dip_ann_RC4US,)) # dip_ann_GSL, dip_ann_extra))
dop_ann_all = np.concatenate((dop_ann_RC4US,)) # dop_ann_GSL, dop_ann_extra))
pp_ann_all = np.concatenate((pp_ann_RC4US,)) # pp_ann_GSL, pp_ann_extra))
dfe_ann_all = np.concatenate((dfe_ann_RC4US,)) # dfe_ann_GSL, dfe_ann_extra))
pfe_ann_all = np.concatenate((pfe_ann_RC4US,)) # pfe_ann_GSL, pfe_ann_extra))
si_ann_all = np.concatenate((si_ann_RC4US,)) # si_ann_GSL, si_ann_extra))
o2_ann_all = np.concatenate((o2_ann_RC4US,)) # o2_ann_GSL, o2_ann_extra))

# Combine all monthly station data
dic_monthly_all = np.concatenate((dic_monthly_RC4US,)) # dic_monthly_GSL, dic_monthly_extra), axis=1)
alk_monthly_all = np.concatenate((alk_monthly_RC4US,)) # alk_monthly_GSL, alk_monthly_extra), axis=1)
no3_monthly_all = np.concatenate((no3_monthly_RC4US,)) # no3_monthly_GSL, no3_monthly_extra), axis=1)
nh4_monthly_all = np.concatenate((nh4_monthly_RC4US,)) # nh4_monthly_GSL, nh4_monthly_extra), axis=1)
din_monthly_all = np.concatenate((din_monthly_RC4US,)) # din_monthly_GSL, din_monthly_extra), axis=1)
don_monthly_all = np.concatenate((don_monthly_RC4US,)) # don_monthly_GSL, don_monthly_extra), axis=1)
pn_monthly_all = np.concatenate((pn_monthly_RC4US,)) # pn_monthly_GSL, pn_monthly_extra), axis=1)
dip_monthly_all = np.concatenate((dip_monthly_RC4US,)) # dip_monthly_GSL, dip_monthly_extra), axis=1)
dop_monthly_all = np.concatenate((dop_monthly_RC4US,)) # dop_monthly_GSL, dop_monthly_extra), axis=1)
pp_monthly_all = np.concatenate((pp_monthly_RC4US,)) # pp_monthly_GSL, pp_monthly_extra), axis=1)
dfe_monthly_all = np.concatenate((dfe_monthly_RC4US,)) # dfe_monthly_GSL, dfe_monthly_extra), axis=1)
pfe_monthly_all = np.concatenate((pfe_monthly_RC4US,)) # pfe_monthly_GSL, pfe_monthly_extra), axis=1)
si_monthly_all = np.concatenate((si_monthly_RC4US,)) # si_monthly_GSL, si_monthly_extra), axis=1)
o2_monthly_all = np.concatenate((o2_monthly_RC4US,)) # o2_monthly_GSL, o2_monthly_extra), axis=1)


###########################################################################
# Load in monthly climatology of river forcing from the regional grid.    #
# File contains:                                                          #
# runoff: monthly average runoff in kg m-2 sec-1                          #
# area: area of grid cell in m-2                                          #
# lon: longitude (0-360 degrees)                                          #
# lat: latitude                                                           #                                                                        #
###########################################################################

import scipy.io as sio
import numpy as np
import xarray as xr

# Load data from MATLAB .mat file
mat_contents = sio.loadmat('./glofas_runoff_mean.mat')
runoff = mat_contents['runoff']

# Convert runoff from kg m^-2 sec^-1 to m^3 sec^-1
area_mod = 1000  # Assuming area_mod is provided in square meters
Q_mod_monthly = np.zeros_like(runoff)
for m in range(12):  #TODO modify range(4) to range(12)
    Q_mod_monthly[m, :, :] = np.squeeze(runoff[m, :, :]) * area_mod / 1000

# Calculate annual mean runoff
Q_mod_ann = np.mean(Q_mod_monthly, axis=0)

# Clear unnecessary variables
del runoff, area_mod

# Load grid data/home/otel/Dropbox/trabalho_irado/Northeastern/202401_mom6_cobalt_compilation/run/NWA0.25.COBALT/INPUT/20210101.ocean_static.nc
grid_file = '/home/otel/Dropbox/trabalho_irado/Northeastern/202401_mom6_cobalt_compilation/run/tools_NWA25.COBALT/data/output/19930101.ocean_static.nc'
with xr.open_dataset(grid_file) as ds:
    depth = ds['deptho'].transpose()

# Replace NaN values in depth with -1
depth = depth.where(~np.isnan(depth), other=-1)


dsgrid = (xr.open_dataset(grid_file)
          .rename(geolon='lon', geolat='lat'))
woa_temp1 = woa_temp[['t_an']].roms.interpolate(dsgrid) 


###########################################################################
# Filter for rivers in the region, set thresholds for minimum river size, #
# set parameters for plotting routines.                                   #
###########################################################################

lat_mod = mat_contents['lat_mod']
lon_mod = mat_contents['lon_mod']



import numpy as np

# Use grid to filter rivers outside domain
lat_mod_max = np.max(lat_mod)
lat_mod_min = np.min(lat_mod)
lon_mod_max = np.max(lon_mod)
lon_mod_min = np.min(lon_mod)

in_region = np.where(
    (lon_stations_all <= lon_mod_max) & (lon_stations_all >= lon_mod_min) &
    (lat_stations_all <= lat_mod_max) & (lat_stations_all >= lat_mod_min) &
    np.isfinite(Q_ann_all) & (Q_ann_all > Q_min)
)[0]

num_rivers = len(in_region)

# Extract data for rivers within the domain
station_names_reg = station_names_all[in_region]
lon_stations_reg = lon_stations_all[in_region]
lat_stations_reg = lat_stations_all[in_region]
Q_ann_reg = Q_ann_all[in_region]
dic_ann_reg = dic_ann_all[in_region]
alk_ann_reg = alk_ann_all[in_region]
no3_ann_reg = no3_ann_all[in_region]
nh4_ann_reg = nh4_ann_all[in_region]
din_ann_reg = din_ann_all[in_region]
don_ann_reg = don_ann_all[in_region]
pn_ann_reg = pn_ann_all[in_region]
dip_ann_reg = dip_ann_all[in_region]
dop_ann_reg = dop_ann_all[in_region]
pp_ann_reg = pp_ann_all[in_region]
dfe_ann_reg = dfe_ann_all[in_region]
pfe_ann_reg = pfe_ann_all[in_region]
si_ann_reg = si_ann_all[in_region]
o2_ann_reg = o2_ann_all[in_region]


dic_monthly_reg = []
alk_monthly_reg = []
no3_monthly_reg = []
nh4_monthly_reg = []
din_monthly_reg = []
don_monthly_reg = []
pn_monthly_reg = []
dip_monthly_reg = []
dop_monthly_reg = []
pp_monthly_reg = []
dfe_monthly_reg = []
pfe_monthly_reg = []
si_monthly_reg = []
o2_monthly_reg = []

for m in range(12): #TODO range(12)
    dic_monthly_reg.append(dic_monthly_all.T[m, in_region])
    alk_monthly_reg.append(alk_monthly_all.T[m, in_region])
    no3_monthly_reg.append(no3_monthly_all.T[m, in_region])
    nh4_monthly_reg.append(nh4_monthly_all.T[m, in_region])
    din_monthly_reg.append(din_monthly_all.T[m, in_region])
    don_monthly_reg.append(don_monthly_all.T[m, in_region])
    pn_monthly_reg.append(pn_monthly_all.T[m, in_region])
    dip_monthly_reg.append(dip_monthly_all.T[m, in_region])
    dop_monthly_reg.append(dop_monthly_all.T[m, in_region])
    pp_monthly_reg.append(pp_monthly_all.T[m, in_region])
    dfe_monthly_reg.append(dfe_monthly_all.T[m, in_region])
    pfe_monthly_reg.append(pfe_monthly_all.T[m, in_region])
    si_monthly_reg.append(si_monthly_all.T[m, in_region])
    o2_monthly_reg.append(o2_monthly_all.T[m, in_region])


###########################################################################
# Assigning outflow points to rivers.                                     #
#  1. Assignment starts with the rivers with the smallest flow and works  #
#     to the largest, w/larger river characteristics taking precedence to #
#     ensure the most significant rivers are well represented.            #
#  2. The algorithm keeps choosing the closest points to each river mouth #
#     until the assigned flow is as close as possible to that observed    #
#  3. Once the outflow points are assigned using the mean flow values,    #
#     monthly concentrations are assigned to those points.                #
#  4. A simple "nearest neighbor" algorithm is used to fill in the gaps   #
###########################################################################

import numpy as np

# Sort rivers by discharge
sort_ind = np.argsort(Q_ann_reg)
station_names_sort = station_names_reg[sort_ind]
lon_stations_sort = lon_stations_reg[sort_ind]
lat_stations_sort = lat_stations_reg[sort_ind]
Q_ann_sort = Q_ann_reg[sort_ind]
dic_ann_sort = dic_ann_reg[sort_ind]
alk_ann_sort = alk_ann_reg[sort_ind]
no3_ann_sort = no3_ann_reg[sort_ind]
nh4_ann_sort = nh4_ann_reg[sort_ind]
din_ann_sort = din_ann_reg[sort_ind]
don_ann_sort = don_ann_reg[sort_ind]
pn_ann_sort = pn_ann_reg[sort_ind]
dip_ann_sort = dip_ann_reg[sort_ind]
dop_ann_sort = dop_ann_reg[sort_ind]
pp_ann_sort = pp_ann_reg[sort_ind]
dfe_ann_sort = dfe_ann_reg[sort_ind]
pfe_ann_sort = pfe_ann_reg[sort_ind]
si_ann_sort = si_ann_reg[sort_ind]
o2_ann_sort = o2_ann_reg[sort_ind]

dic_monthly_sort = []
alk_monthly_sort = []
no3_monthly_sort = []
nh4_monthly_sort = []
din_monthly_sort = []
don_monthly_sort = []
pn_monthly_sort = []
dip_monthly_sort = []
dop_monthly_sort = []
pp_monthly_sort = []
dfe_monthly_sort = []
pfe_monthly_sort = []
si_monthly_sort = []
o2_monthly_sort = []

for m in range(12):  #TODO range(12)
    dic_monthly_sort.append(dic_monthly_reg[m][sort_ind])
    alk_monthly_sort.append(alk_monthly_reg[m][sort_ind])
    no3_monthly_sort.append(no3_monthly_reg[m][sort_ind])
    nh4_monthly_sort.append(nh4_monthly_reg[m][sort_ind])
    din_monthly_sort.append(din_monthly_reg[m][sort_ind])
    don_monthly_sort.append(don_monthly_reg[m][sort_ind])
    pn_monthly_sort.append(pn_monthly_reg[m][sort_ind])
    dip_monthly_sort.append(dip_monthly_reg[m][sort_ind])
    dop_monthly_sort.append(dop_monthly_reg[m][sort_ind])
    pp_monthly_sort.append(pp_monthly_reg[m][sort_ind])
    dfe_monthly_sort.append(dfe_monthly_reg[m][sort_ind])
    pfe_monthly_sort.append(pfe_monthly_reg[m][sort_ind])
    si_monthly_sort.append(si_monthly_reg[m][sort_ind])
    o2_monthly_sort.append(o2_monthly_reg[m][sort_ind])

# Create vectors of values at the runoff points from the model grid
ind_ro = np.where(Q_mod_ann > 0)
Q_mod_vec = Q_mod_ann[ind_ro]
lon_mod_runoff_vec = lon_mod[ind_ro]
lat_mod_runoff_vec = lat_mod[ind_ro]
Q_mod_monthly_vecs = np.zeros((12, len(ind_ro[0])))  #TODO zeros(12...)
for m in range(12): #TODO range(12)
    Q_mod_monthly_vecs[m, :] = Q_mod_monthly[m,ind_ro[0],ind_ro[1]]

# Create a grid of saturated oxygen values using the world ocean atlas data
temp_woa_monthly_vecs = np.zeros((12, len(ind_ro[0]))) #TODO zeros((1,..))
o2sat_woa_monthly_vecs = np.zeros((12, len(ind_ro[0]))) #TODO zeros((1,..))

# Constants for o2 saturation calculation
a_0 = 2.00907
a_1 = 3.22014
a_2 = 4.05010
a_3 = 4.94457
a_4 = -2.56847e-1
a_5 = 3.88767
sal = 0
b_0 = -6.24523e-3
b_1 = -7.37614e-3
b_2 = -1.03410e-2
b_3 = -8.17083e-3
c_0 = -4.88682e-7


for m in range(12): #TODO range(12) woa should be the monthly one
    temp = woa_temp1['t_an'][m, :, :].squeeze()
    temp = temp.bfill('xh').ffill('yh').values
    temp[temp > 40] = 40
    temp[temp < 0] = 0
    temp_woa_monthly_vecs[m,:] = temp[m,ind_ro[0],ind_ro[1]]

    tt = 298.15 - temp_woa_monthly_vecs[m]
    tkb = 273.15 + temp_woa_monthly_vecs[m]
    ts = np.log(tt / tkb)
    ts2 = ts * ts
    ts3 = ts2 * ts
    ts4 = ts3 * ts
    ts5 = ts4 * ts

    o2sat_woa_monthly_vecs[m,:] = (1000.0 / 22391.6) * 1000 * np.exp(
        a_0 + a_1 * ts + a_2 * ts2 + a_3 * ts3 + a_4 * ts4 + a_5 * ts5 +
        (b_0 + b_1 * ts + b_2 * ts2 + b_3 * ts3 + c_0 * sal) * sal
    )
print(np.shape(o2sat_woa_monthly_vecs))



import xarray as xr



import xarray as xr

# Read the NetCDF file
NEWS_file = xr.open_dataset(NEWS_file)

# Read data from the NetCDF file
din_ann_NEWS = NEWS_file['NO3_CONC'].values
aa = np.where(din_ann_NEWS > 0)  # Finding indices where din_ann_NEWS > 0
temp1 = NEWS_file['LDON_CONC'].values
temp2 = NEWS_file['SLDON_CONC'].values
temp3 = NEWS_file['SRDON_CONC'].values
don_ann_NEWS = temp1 + temp2 + temp3
don_ratio_NEWS_vec = don_ann_NEWS[aa] / din_ann_NEWS[aa]

pn_ann_NEWS = NEWS_file['NDET_CONC'].values
pn_ratio_NEWS_vec = pn_ann_NEWS[aa] / din_ann_NEWS[aa]

dip_ann_NEWS = NEWS_file['PO4_CONC'].values
dip_ratio_NEWS_vec = dip_ann_NEWS[aa] / din_ann_NEWS[aa]
temp1 = NEWS_file['LDOP_CONC'].values
temp2 = NEWS_file['SLDOP_CONC'].values
temp3 = NEWS_file['SRDOP_CONC'].values
dop_ann_NEWS = temp1 + temp2 + temp3
dop_ratio_NEWS_vec = dop_ann_NEWS[aa] / din_ann_NEWS[aa]

pp_ann_NEWS = NEWS_file['PDET_CONC'].values
pp_ratio_NEWS_vec = pp_ann_NEWS[aa] / din_ann_NEWS[aa]

si_ann_NEWS = NEWS_file['SI_CONC'].values
si_ratio_NEWS_vec = si_ann_NEWS[aa] / din_ann_NEWS[aa]

# Close the NetCDF file
NEWS_file.close()


# Calculate ratios relative to DIN
aa = din_ann_NEWS > 0
don_ann_NEWS = temp1 + temp2 + temp3
don_ratio_NEWS_vec = don_ann_NEWS[aa] / din_ann_NEWS[aa]
dip_ratio_NEWS_vec = dip_ann_NEWS[aa] / din_ann_NEWS[aa]
dop_ann_NEWS = temp1 + temp2 + temp3
dop_ratio_NEWS_vec = dop_ann_NEWS[aa] / din_ann_NEWS[aa]
pp_ratio_NEWS_vec = pp_ann_NEWS[aa] / din_ann_NEWS[aa]
si_ratio_NEWS_vec = si_ann_NEWS[aa] / din_ann_NEWS[aa]

# Vectors to hold monthly values mapped onto model runoff points
#TODO zeros((12,...))
dic_mod_monthly_vecs = np.zeros((12, lon_mod_runoff_vec.size))
alk_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
no3_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
nh4_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
din_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
don_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
pn_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
dip_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
dop_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
pp_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
dfe_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
pfe_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
si_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))
o2_mod_monthly_vecs = np.zeros((12, lat_mod_runoff_vec.size))


###########################################################################
# Loop identifies points assigned to each river                           #
###########################################################################


import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


# Loop identifies points assigned to each river
for k in range(num_rivers):
# Assuming lon_stations_sort and lat_stations_sort are numpy arrays

    dist = cdist([[lon_stations_sort[k], lat_stations_sort[k]]],
                 np.column_stack([lon_mod_runoff_vec, lat_mod_runoff_vec]))[0]
    dist_sort_ind = np.argsort(dist)
    
    
    if dist[dist_sort_ind[0]] < min_dist:
        Q_sum1 = 0
        Q_sum2 = 0
        n = 0
        while Q_sum2 < Q_ann_sort[k] and dist[dist_sort_ind[n+1]] < max_dist:
            Q_sum1 = Q_sum2
            n += 1
            Q_sum2 = Q_sum1 + Q_mod_vec[dist_sort_ind[n]]
        if abs(Q_sum1 - Q_ann_sort[k]) < abs(Q_sum2 - Q_ann_sort[k]):
            nrp = n-1  # number of runoff points
        else:
            nrp = n

        # enter monthly concentration values into an array of monthly values
        # if no monthly value is available, use annuals.
        for m in range(12): #TODO range(12)
            if np.isnan(dic_monthly_sort[m][k]):
                if np.isfinite(dic_ann_sort[k]):
                    dic_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = dic_ann_sort[k]
                else:
                    dic_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = 0
            else:
                dic_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = dic_monthly_sort[m][k]

            if np.isnan(alk_monthly_sort[m][k]):
                if np.isfinite(dic_ann_sort[k]):
                    alk_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = alk_ann_sort[k]
                else:
                    alk_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = 0
            else:
                alk_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = alk_monthly_sort[m][k]

            # mapping assumes that DIN is defined for nutrient calculations since
            # ratios relative to DIN are used to fill in other components.  If
            # DIN is not defined, values are left at 0 and eventually filled with
            # a nearest neighbor filling (next section)
            if np.isfinite(din_monthly_sort[m][k]) or np.isfinite(din_ann_sort[k]):
                if np.isfinite(din_monthly_sort[m][k]):
                    din_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = din_monthly_sort[m][k]
                elif np.isfinite(din_ann_sort[k]):
                    din_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = din_ann_sort[k]

                if np.isfinite(no3_monthly_sort[m][k]):
                    no3_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = no3_monthly_sort[m][k]
                elif np.isfinite(no3_ann_sort[k]):
                    no3_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = no3_ann_sort[k]

                if np.isfinite(nh4_monthly_sort[m][k]):
                    nh4_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = nh4_monthly_sort[m][k]
                elif np.isfinite(nh4_ann_sort[k]):
                    nh4_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = nh4_ann_sort[k]

                if np.isnan(don_monthly_sort[m][k]) and np.isnan(don_ann_sort[k]):
                    don_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = \
                        din_mod_monthly_vecs[m, dist_sort_ind[:nrp]] * don_ratio_NEWS_vec[dist_sort_ind[:nrp]]
                elif np.isnan(don_monthly_sort[m][k]) and ~np.isnan(don_ann_sort[k]):
                    don_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = don_ann_sort[k]
                else:
                    don_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = don_monthly_sort[m][k]

                # Populate the monthly concentration values for various nutrients or use annual values if monthly data is not available
                if np.isnan(pn_monthly_sort[m][k]) and np.isnan(pn_ann_sort[k]):
                    pn_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = din_mod_monthly_vecs[m, dist_sort_ind[:nrp]] * pn_ratio_NEWS_vec[dist_sort_ind[:nrp]]
                elif np.isnan(pn_monthly_sort[m][k]) and ~np.isnan(pn_ann_sort[k]):
                    pn_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = pn_ann_sort[k]
                else:
                    pn_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = pn_monthly_sort[m][k]

                if np.isnan(dip_monthly_sort[m][k]) and np.isnan(dip_ann_sort[k]):
                    dip_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = din_mod_monthly_vecs[m, dist_sort_ind[:nrp]] * dip_ratio_NEWS_vec[dist_sort_ind[:nrp]]
                elif np.isnan(dip_monthly_sort[m][k]) and ~np.isnan(dip_ann_sort[k]):
                    dip_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = dip_ann_sort[k]
                else:
                    dip_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = dip_monthly_sort[m][k]

                if np.isnan(dop_monthly_sort[m][k]) and np.isnan(dop_ann_sort[k]):
                    dop_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = din_mod_monthly_vecs[m, dist_sort_ind[:nrp]] * dop_ratio_NEWS_vec[dist_sort_ind[:nrp]]
                elif np.isnan(dop_monthly_sort[m][k]) and ~np.isnan(dop_ann_sort[k]):
                    dop_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = dop_ann_sort[k]
                else:
                    dop_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = dop_monthly_sort[m][k]

                if np.isnan(pp_monthly_sort[m][k]) and np.isnan(pp_ann_sort[k]):
                    pp_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = din_mod_monthly_vecs[m, dist_sort_ind[:nrp]] * pp_ratio_NEWS_vec[dist_sort_ind[:nrp]]
                elif np.isnan(pp_monthly_sort[m][k]) and ~np.isnan(pp_ann_sort[k]):
                    pp_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = pp_ann_sort[k]
                else:
                    pp_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = pp_monthly_sort[m][k]

                if np.isnan(si_monthly_sort[m][k]) and np.isnan(si_ann_sort[k]):
                    si_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = din_mod_monthly_vecs[m, dist_sort_ind[:nrp]] * si_ratio_NEWS_vec[dist_sort_ind[:nrp]]
                elif np.isnan(si_monthly_sort[m][k]) and ~np.isnan(si_ann_sort[k]):
                    si_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = si_ann_sort[k]
                else:
                    si_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = si_monthly_sort[m][k]

                if np.isnan(o2_monthly_sort[m][k]):
                    o2_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = o2sat_woa_monthly_vecs[m, dist_sort_ind[:nrp]]
                else:
                    o2_mod_monthly_vecs[m, dist_sort_ind[:nrp]] = o2_monthly_sort[m][k]


        # plot to check location if inspect_map == 'y'.  The plot puts
        # open circles at each runoff location in the model grid and fills
        # those that are assigned to each river.  Note that some of the smaller
        # rivers may be replaced with larger ones as the fitting process
        # continues.
                            
        if inspect_map == 'y':
            plt.figure(1)
            plt.clf()
            ax = plt.axes(projection='3d')
            sc = ax.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, np.log10(Q_mod_vec), s=3, c=np.log10(Q_mod_vec))
            ax.scatter(lon_mod_runoff_vec[dist_sort_ind[0:nrp]], lat_mod_runoff_vec[dist_sort_ind[0:nrp]],
                    np.log10(Q_mod_vec[dist_sort_ind[0:nrp]]), s=25, c=np.log10(Q_mod_vec[dist_sort_ind[0:nrp]]), marker='o', cmap='viridis')
            ax.view_init(90, 0)
            ax.plot(lon_stations_sort[k], lat_stations_sort[k], 1e5, 'k.', markersize=20)
            cr =ax.contour(lon_mod, lat_mod, depth.T, [0], colors='k')
            # ax.set_xlim(lon_stations_sort[k] - plot_width/2, lon_stations_sort[k] + plot_width/2)
            # ax.set_ylim(lat_stations_sort[k] - plot_width/2, lat_stations_sort[k] + plot_width/2)
            # ax.set_zlim(np.log10(Q_mod_vec.min()), np.log10(Q_mod_vec.max()))
            ax.set_title('River number: {} Name: {}'.format(k, station_names_sort[k]))
            plt.colorbar(sc)
            
            # Check the values of each mapping
            N_check = np.mean(din_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            P_check = np.mean(dip_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            SI_check = np.mean(si_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            DIN_check = np.mean(din_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            DON_check = np.mean(don_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            PN_check = np.mean(pn_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            DIP_check = np.mean(dip_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            DOP_check = np.mean(dop_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            PP_check = np.mean(pp_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            SI_check = np.mean(si_mod_monthly_vecs[:, dist_sort_ind[0:nrp]])
            
            print(station_names_sort[k])
            ind = dist_sort_ind[0:nrp]
            print('total flow in m3 sec')
            print(Q_ann_sort[k], np.sum(Q_mod_vec[dist_sort_ind[0:nrp]]))
            print('N, P conc (mmoles m-3), DI, DO, P')
            print([DIN_check, DON_check, PN_check])
            print([DIP_check, DOP_check, PP_check])
            print('Total N, Total P, Total N: Total P')
            print(N_check, P_check, N_check / P_check)
            print('DO:DI and P:DI ratios')
            print([DON_check / DIN_check, PN_check / DIN_check])
            print([DOP_check / DIP_check, PP_check / DIP_check])
            print('silica concentration (mmoles m-3)')
            print(SI_check)
            
            plt.pause(0.1)

        # If river is outside the domain, skip all of the calculations above and
        # just plot for inspection/evaluation
        else:
            # This is for rivers that were outside of the domain
            if inspect_map == 'y':
                plt.figure(1)
                plt.clf()
                plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=3, c=np.log10(Q_mod_vec))
                plt.plot(lon_stations_sort[k], lat_stations_sort[k], 'k.', markersize=20)
                plt.axis([lon_stations_sort[k] - 10, lon_stations_sort[k] + 10,
                        lat_stations_sort[k] - 10, lat_stations_sort[k] + 10])
                plt.title('OUTSIDE: River number: {} Name: {}'.format(k, station_names_sort[k]))
                plt.colorbar()
                
                plt.pause(0.1)


###########################################################################
# nearest neighbor search to fill in any runoff points that were not      #
# assigned after the runoff mapping step                                  #
###########################################################################

def nearest_neighbor_fill(vec, lon, lat, m):
    vec_out = vec.copy()
    aa_dic = np.where(vec[m, :] == 0)[0]
    if len(aa_dic)>0:
        bb_dic = np.where(vec[m, :] > 0)[0]
        F_dic = interp.interp2d(lon[bb_dic], lat[bb_dic], vec[m, bb_dic], kind='linear')
        vec_out[m, aa_dic] = np.diag(F_dic(lon[aa_dic], lat[aa_dic]))
    return vec_out


from scipy import interpolate as interp


for m in range(12):
    dic_mod_monthly_vecs[m,:] = nearest_neighbor_fill(dic_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]
    alk_mod_monthly_vecs[m,:] = nearest_neighbor_fill(alk_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]
    no3_mod_monthly_vecs[m,:] = nearest_neighbor_fill(no3_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]

    nh4_mod_monthly_vecs[m,:] = nearest_neighbor_fill(nh4_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]

    din_mod_monthly_vecs[m,:] = nearest_neighbor_fill(din_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]

    pn_mod_monthly_vecs[m,:] = nearest_neighbor_fill(pn_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]

    dip_mod_monthly_vecs[m,:] = nearest_neighbor_fill(dip_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]

    dop_mod_monthly_vecs[m,:] = nearest_neighbor_fill(dop_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]

    pp_mod_monthly_vecs[m,:] = nearest_neighbor_fill(pp_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]

    si_mod_monthly_vecs[m,:] = nearest_neighbor_fill(si_mod_monthly_vecs,
                                                 lon_mod_runoff_vec,
                                                 lat_mod_runoff_vec,
                                                 m)[m]




# For o2sat, fill in any 0 values with saturated o2 at the world ocean atlas climatology
for m in range(12):
    aa = np.where(o2_mod_monthly_vecs[m, :] == 0)
    o2_mod_monthly_vecs[m, list(aa[0])] = o2sat_woa_monthly_vecs[m][list(aa[0])]

totn_mod_monthly_vecs = din_mod_monthly_vecs + don_mod_monthly_vecs + pn_mod_monthly_vecs
totp_mod_monthly_vecs = dip_mod_monthly_vecs + dop_mod_monthly_vecs + pp_mod_monthly_vecs

dicflux_mod_monthly_vecs = dic_mod_monthly_vecs * Q_mod_monthly_vecs
alkflux_mod_monthly_vecs = alk_mod_monthly_vecs * Q_mod_monthly_vecs
dinflux_mod_monthly_vecs = din_mod_monthly_vecs * Q_mod_monthly_vecs
no3flux_mod_monthly_vecs = no3_mod_monthly_vecs * Q_mod_monthly_vecs
nh4flux_mod_monthly_vecs = nh4_mod_monthly_vecs * Q_mod_monthly_vecs
dipflux_mod_monthly_vecs = dip_mod_monthly_vecs * Q_mod_monthly_vecs
donflux_mod_monthly_vecs = don_mod_monthly_vecs * Q_mod_monthly_vecs
dopflux_mod_monthly_vecs = dop_mod_monthly_vecs * Q_mod_monthly_vecs
pnflux_mod_monthly_vecs = pn_mod_monthly_vecs * Q_mod_monthly_vecs
ppflux_mod_monthly_vecs = pp_mod_monthly_vecs * Q_mod_monthly_vecs
siflux_mod_monthly_vecs = si_mod_monthly_vecs * Q_mod_monthly_vecs
totnflux_mod_monthly_vecs = totn_mod_monthly_vecs * Q_mod_monthly_vecs
totpflux_mod_monthly_vecs = totp_mod_monthly_vecs * Q_mod_monthly_vecs
o2flux_mod_monthly_vecs = o2_mod_monthly_vecs * Q_mod_monthly_vecs

dicflux_mod_ann_vec = np.mean(dicflux_mod_monthly_vecs, axis=0)
alkflux_mod_ann_vec = np.mean(alkflux_mod_monthly_vecs, axis=0)
dinflux_mod_ann_vec = np.mean(dinflux_mod_monthly_vecs, axis=0)
no3flux_mod_ann_vec = np.mean(no3flux_mod_monthly_vecs, axis=0)
nh4flux_mod_ann_vec = np.mean(nh4flux_mod_monthly_vecs, axis=0)
dipflux_mod_ann_vec = np.mean(dipflux_mod_monthly_vecs, axis=0)
donflux_mod_ann_vec = np.mean(donflux_mod_monthly_vecs, axis=0)
dopflux_mod_ann_vec = np.mean(dopflux_mod_monthly_vecs, axis=0)
pnflux_mod_ann_vec = np.mean(pnflux_mod_monthly_vecs, axis=0)
ppflux_mod_ann_vec = np.mean(ppflux_mod_monthly_vecs, axis=0)
siflux_mod_ann_vec = np.mean(siflux_mod_monthly_vecs, axis=0)
totnflux_mod_ann_vec = np.mean(totnflux_mod_monthly_vecs, axis=0)
totpflux_mod_ann_vec = np.mean(totpflux_mod_monthly_vecs, axis=0)
o2flux_mod_ann_vec = np.mean(o2flux_mod_monthly_vecs, axis=0)





import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Scale marker size with the freshwater flux
ms_vec = np.zeros_like(Q_mod_vec)
ms_vec[np.log10(Q_mod_vec) < 0] = 1
ms_vec[(np.log10(Q_mod_vec) > 0) & (np.log10(Q_mod_vec) < 1)] = 2.5
ms_vec[(np.log10(Q_mod_vec) > 1) & (np.log10(Q_mod_vec) < 2)] = 10
ms_vec[(np.log10(Q_mod_vec) > 2) & (np.log10(Q_mod_vec) < 3)] = 25
ms_vec[np.log10(Q_mod_vec) > 3] = 100

# DIC, Alk concentrations and DIC:Alk
if plot_progress:
    for m in range(12): #TODO range(12)
        plt.figure(figsize=(15, 5))

        plt.subplot(1, 3, 1, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=dic_mod_monthly_vecs[m, :], s=ms_vec,  cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'DIC, mmoles m-3; month = {m + 1}')

        plt.subplot(1, 3, 2, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=alk_mod_monthly_vecs[m, :], s=ms_vec,  cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'Alk, meq m-3; month = {m + 1}')

        plt.subplot(1, 3, 3, projection='3d')
        dic_alk_ratio = dic_mod_monthly_vecs[m, :] / alk_mod_monthly_vecs[m, :]
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=dic_alk_ratio, s=ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title('DIC:Alk ratio')

        plt.show()

    # Nitrogen Concentrations
    for m in range(12): #TODO range(12)
        plt.figure(figsize=(15, 10))

        plt.subplot(2, 2, 1, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=din_mod_monthly_vecs[m, :], s=ms_vec,  cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'DIN, mmoles m-3; month = {m + 1}')

        plt.subplot(2, 2, 2, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=no3_mod_monthly_vecs[m, :], s=ms_vec,  cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'no3, mmoles m-3; month = {m + 1}')

        plt.subplot(2, 2, 3, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=nh4_mod_monthly_vecs[m, :], s=ms_vec,  cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'nh4, mmoles m-3; month = {m + 1}')

        plt.subplot(2, 2, 4, projection='3d')
        no3_din_ratio = no3_mod_monthly_vecs[m, :] / din_mod_monthly_vecs[m, :]
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c= no3_din_ratio, s=ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title('NO3:DIN ratio')

        plt.show()

    for m in range(12): #TODO range(12)
        plt.figure(figsize=(15, 10))

        plt.subplot(2, 3, 1, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=din_mod_monthly_vecs[m, :], s=ms_vec,  cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'DIN, mmoles m-3; month = {m + 1}')

        plt.subplot(2, 3, 2, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=don_mod_monthly_vecs[m, :], s=ms_vec,  cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'DON, mmoles m-3; month = {m + 1}')

        plt.subplot(2, 3, 3, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=pn_mod_monthly_vecs[m, :],s= ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'PN, mmoles m-3; month = {m + 1}')

        plt.subplot(2, 3, 5, projection='3d')
        don_din_ratio = don_mod_monthly_vecs[m, :] / din_mod_monthly_vecs[m, :]
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c= don_din_ratio, s=ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title('DON:DIN ratio')

        plt.subplot(2, 3, 6, projection='3d')
        pn_din_ratio = pn_mod_monthly_vecs[m, :] / din_mod_monthly_vecs[m, :]
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c= pn_din_ratio, s=ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title('PN:DIN ratio')

        plt.show()

    # Phosphorus Concentrations
    for m in range(12): #TODO range(12)
        plt.figure(figsize=(15, 10))

        plt.subplot(2, 3, 1, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=dip_mod_monthly_vecs[m, :], s=ms_vec,  cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'DIP, mmoles m-3; month = {m + 1}')

        plt.subplot(2, 3, 2, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=dop_mod_monthly_vecs[m, :], s=ms_vec,  cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'DOP, mmoles m-3; month = {m + 1}')


        plt.subplot(2, 3, 3, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=pp_mod_monthly_vecs[m, :], s=ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'PP, mmoles m-3; month = {m + 1}')

        plt.subplot(2, 3, 5, projection='3d')
        dop_dip_ratio = dop_mod_monthly_vecs[m, :] / dip_mod_monthly_vecs[m, :]
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=dop_dip_ratio, s=ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title('DOP:DIP ratio')

        plt.subplot(2, 3, 6, projection='3d')
        pp_dip_ratio = pp_mod_monthly_vecs[m, :] / dip_mod_monthly_vecs[m, :]
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=pp_dip_ratio, s=ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title('PP:DIP ratio')

        plt.show()

    # Silica and Oxygen concentrations
    for m in range(12): #TODO range(12)
        plt.figure(figsize=(15, 10))

        plt.subplot(2, 1, 1, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=si_mod_monthly_vecs[m, :], s=ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'Si, mmoles m-3; month = {m + 1}')

        plt.subplot(2, 1, 2, projection='3d')
        plt.scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, c=o2_mod_monthly_vecs[m, :], s=ms_vec, cmap='viridis', marker='o')
        plt.colorbar()
        plt.title(f'o2, mmoles m-3; month = {m + 1}')

        plt.show()



# Initialize 2D concentration arrays; these are the ones read into MOM6 to
# specify the nutrient concentrations of river inputs.
import numpy as np

# Define concentration arrays
DIC_CONC = np.zeros((12, lon_mod.shape[0], lon_mod.shape[1]))
ALK_CONC = np.zeros_like(DIC_CONC)
NO3_CONC = np.zeros_like(DIC_CONC)
NH4_CONC = np.zeros_like(DIC_CONC)
LDON_CONC = np.zeros_like(DIC_CONC)
SLDON_CONC = np.zeros_like(DIC_CONC)
SRDON_CONC = np.zeros_like(DIC_CONC)
PO4_CONC = np.zeros_like(DIC_CONC)
LDOP_CONC = np.zeros_like(DIC_CONC)
SLDOP_CONC = np.zeros_like(DIC_CONC)
SRDOP_CONC = np.zeros_like(DIC_CONC)
NDET_CONC = np.zeros_like(DIC_CONC)
PDET_CONC = np.zeros_like(DIC_CONC)
SI_CONC = np.zeros_like(DIC_CONC)
O2_CONC = np.zeros_like(DIC_CONC)

# Map concentration vectors onto 2D arrays
temp = np.zeros_like(lon_mod)
for m in range(12):
    temp[ind_ro] = dic_mod_monthly_vecs[m][:]
    DIC_CONC[m, :, :] = temp.copy()
    temp[ind_ro] = alk_mod_monthly_vecs[m][:]
    ALK_CONC[m, :, :] = temp.copy()

    temp[ind_ro] = no3_mod_monthly_vecs[m][:]
    NO3_CONC[m, :, :] = temp.copy()
    temp[ind_ro] = nh4_mod_monthly_vecs[m][:]
    NH4_CONC[m, :, :] = temp.copy()
    temp[ind_ro] = don_mod_monthly_vecs[m][:]
    LDON_CONC[m, :, :] = frac_ldon * temp
    SLDON_CONC[m, :, :] = frac_sldon * temp
    SRDON_CONC[m, :, :] = frac_srdon * temp
    temp[ind_ro] = pn_mod_monthly_vecs[m][:]
    NDET_CONC[m, :, :] = temp.copy()

    temp[ind_ro] = dip_mod_monthly_vecs[m][:]
    PO4_CONC[m, :, :] = temp.copy()
    temp[ind_ro] = dop_mod_monthly_vecs[m][:]
    LDOP_CONC[m, :, :] = frac_ldop * temp
    SLDOP_CONC[m, :, :] = frac_sldop * temp
    SRDOP_CONC[m, :, :] = frac_srdop * temp
    temp[ind_ro] = pp_mod_monthly_vecs[m][:]
    PDET_CONC[m, :, :] = temp.copy()

    temp[ind_ro] = si_mod_monthly_vecs[m][:]
    SI_CONC[m, :, :] = temp.copy()
    temp[ind_ro] = o2_mod_monthly_vecs[m][:]
    O2_CONC[m, :, :] = temp.copy()

# Adjust units for consistency
DIC_CONC /= 1e3
ALK_CONC /= 1e3
NO3_CONC /= 1e3
NH4_CONC /= 1e3
LDON_CONC /= 1e3
SLDON_CONC /= 1e3
SRDON_CONC /= 1e3
NDET_CONC /= 1e3
PO4_CONC /= 1e3
LDOP_CONC /= 1e3
SLDOP_CONC /= 1e3
SRDOP_CONC /= 1e3
PDET_CONC /= 1e3
SI_CONC /= 1e3
O2_CONC /= 1e3

# Add iron concentrations
FED_CONC = NO3_CONC.copy()
FEDET_CONC = NO3_CONC.copy()
FED_CONC[FED_CONC > 0] = const_fed
FEDET_CONC[FEDET_CONC > 0] = 0.0

ms = 8  # Marker size


import matplotlib.pyplot as plt

if plot_progress:
    for m in range(12):
        fig, axs = plt.subplots(2, 1)
        fig.suptitle(f'Month {m+1}')
        scatter1 = axs[0].scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(DIC_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
        axs[0].set_title('log10(DIC CONC)')
        axs[0].set_xlabel('Longitude')
        axs[0].set_ylabel('Latitude')
        axs[0].grid(True)
        plt.colorbar(scatter1, ax=axs[0], orientation='vertical')
        
        scatter2 = axs[1].scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(ALK_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
        axs[1].set_title('log10(ALK CONC)')
        axs[1].set_xlabel('Longitude')
        axs[1].set_ylabel('Latitude')
        axs[1].grid(True)
        plt.colorbar(scatter2, ax=axs[1], orientation='vertical')
        
        plt.show()

    # Nitrogen
    for m in range(12):
        fig, axs = plt.subplots(3, 2)
        fig.suptitle(f'Month {m+1}')
        axs = axs.flatten()
        for i, ax in enumerate(axs):
            if i == 0:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(NO3_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(NO3 CONC)')
            elif i == 1:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(NH4_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(NH4 CONC)')
            elif i == 2:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(LDON_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(LDON CONC)')
            elif i == 3:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(SLDON_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(SLDON CONC)')
            elif i == 4:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(SRDON_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(SRDON CONC)')
            elif i == 5:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(NDET_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(NDET CONC)')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            ax.grid(True)
            plt.colorbar(scatter, ax=ax, orientation='vertical')
        plt.show()

    # Phosphorus
    for m in range(12):
        fig, axs = plt.subplots(3, 2)
        fig.suptitle(f'Month {m+1}')
        axs = axs.flatten()
        for i, ax in enumerate(axs):
            if i == 0:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(PO4_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(PO4 CONC)')
            elif i == 1:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(LDOP_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(LDOP CONC)')
            elif i == 2:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(SLDOP_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(SLDOP CONC)')
            elif i == 3:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(SRDOP_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(SRDOP CONC)')
            elif i == 4:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(PDET_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(PDET CONC)')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            ax.grid(True)
            plt.colorbar(scatter, ax=ax, orientation='vertical')
        plt.show()

    # Iron, Silica, Oxygen
    for m in range(12):
        fig, axs = plt.subplots(3, 2)
        fig.suptitle(f'Month {m+1}')
        axs = axs.flatten()
        for i, ax in enumerate(axs):
            if i == 0:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(FED_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(FED CONC)')
            elif i == 1:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(FEDET_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(FEDET CONC)')
            elif i == 2:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(SI_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(SI CONC)')
            elif i == 3:
                scatter = ax.scatter(lon_mod.flatten(), lat_mod.flatten(), c=np.log10(O2_CONC[m, :, :]).flatten(), s=ms, cmap='viridis', marker='o')
                ax.set_title('log10(O2 CONC)')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            ax.grid(True)
            plt.colorbar(scatter, ax=ax, orientation='vertical')
        plt.show()



###########################################################################
# Save Files                                                              #
###########################################################################
# option to save matlab file
# save River_DIC_ALK_RC4US_NWA ALK_CONC DIC_CONC

# Construct netcdf file following format used by nutrient input files to
# MOM6


import xarray as xr
import numpy as np
import datetime as dtt
import netCDF4

# Make reference dates for standard non-leap year
dates = np.array([[1993, 1, 16, 12, 0, 0],
                  [1993, 2, 15, 0, 0, 0],
                  [1993, 3, 16, 12, 0, 0],
                  [1993, 4, 16, 0, 0, 0],
                  [1993, 5, 16, 12, 0, 0],
                  [1993, 6, 16, 0, 0, 0],
                  [1993, 7, 16, 12, 0, 0],
                  [1993, 8, 16, 12, 0, 0],
                  [1993, 9, 16, 0, 0, 0],
                  [1993, 10, 16, 12, 0, 0],
                  [1993, 11, 16, 0, 0, 0],
                  [1993, 12, 16, 12, 0, 0]])

aux = [dtt.datetime(*i) for i in dates]
aux1 = [netCDF4.date2num(t, 'days since 1990-01-01T00:00:00') for t in aux]
time = np.array([netCDF4.num2date(t, 'days since 1990-01-01T00:00:00') for t in aux1])


nlat = lat_mod.shape[0]
nlon = lat_mod.shape[1]

# Create xarray dataset

variables = {
    'DIC_CONC': DIC_CONC,
    'ALK_CONC': ALK_CONC,
    'NO3_CONC': NO3_CONC,
    'NH4_CONC': NH4_CONC,
    'LDON_CONC': LDON_CONC,
    'SLDON_CONC': SLDON_CONC,
    'SRDON_CONC': SRDON_CONC,
    'NDET_CONC': NDET_CONC,
    'PO4_CONC': PO4_CONC,
    'LDOP_CONC': LDOP_CONC,
    'SLDOP_CONC': SLDOP_CONC,
    'SRDOP_CONC': SRDOP_CONC,
    'PDET_CONC': PDET_CONC,
    'FED_CONC': FED_CONC,
    'FEDET_CONC': FEDET_CONC,
    'O2_CONC': O2_CONC,
    'SI_CONC': SI_CONC
}

data_vars = {}
for var_name, var_data in variables.items():
    dims = ('time', 'y', 'x')
    data_vars[var_name] = ((dims), var_data)


coords = {
    'time': ('time', time),
    'lat': (('y', 'x'), lat_mod),
    'lon': (('y', 'x'), lon_mod)
}

ds = xr.Dataset(data_vars=data_vars, coords= coords)
# # Add time dimension
# ds['time'] = np.array(time)

# # Add latitude and longitude variables
# ds['lat'] = (('y', 'x'), lat_mod)
# ds['lon'] = (('y', 'x'), lon_mod)

# Add other variables


# Add attributes
ds.attrs['title'] = 'Nutrient Concentrations'
ds.attrs['history'] = 'Created using xarray'
# ds['time'].attrs['units'] = 'days since 1990-01-01T00:00:00'
# ds['time'].attrs['calendar'] = 'gregorian'
# Save to netCDF file
# ds.to_netcdf(nc_file_name)


# import cftime
# t = [cftime.datetime(1993, i+1, 1, calendar='365_day') for i in range(12)]
# ds = ds.assign_coords(time=t)
# Save to netCDF file
ds.to_netcdf(nc_file_name,
        format='NETCDF4',
        engine='netcdf4',
        encoding={'time': {'units': 'days since 1990-01-01T00:00:00', 'calendar': '365_day'}},
        unlimited_dims='time'
)
