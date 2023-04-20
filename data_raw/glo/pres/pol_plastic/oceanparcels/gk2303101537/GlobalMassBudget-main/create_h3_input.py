#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 10:10:26 2021

@author: kaandorp
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from h3.unstable import vect
import h3
# import pickle
from tools import from_pickle, pcolorhex, get_colors
import pandas as pd
import cartopy.crs as ccrs

def to_netcdf(output_filename,data,data_name,h3_index_ocean,time,explanation=''):
    '''
    All data is written to netcdf files to speed up computations
    '''
    dict_data = {}
    for name_, data_ in zip(data_name, data):
        dict_data[name_] = (( "time","h3_index"), data_ )
    dict_data['explanation'] = explanation
        
    ds = xr.Dataset(
        dict_data,
        coords={
            "h3_index": h3_index_ocean,
            "time": time,
        },
    )   
    ds.to_netcdf(output_filename)
    
    
h3_grid_res = 3
input_files = ['input_coastal_pop.nc','input_rivers.nc','input_fisheries_monthly.nc']

#3: ~40,000 cells, 60km
#2: ~6,000 cells, 160km
dict_coastlines = from_pickle('00_data_files/coastal_properties_h3_res%i_50m.pickle' % h3_grid_res)
h3_index_ocean = np.array([hex_ for hex_ in dict_coastlines.keys()], dtype=np.uint64)


#%% create coastal population h3 input

data = xr.open_dataset('00_data_files/' + input_files[0])
X,Y = np.meshgrid(data['lon'],data['lat'])
data_input_h3_pop = np.zeros([len(data['time']),len(h3_index_ocean)])

for i1 in range(len(data['time'])):
    mask = data['MPW_input'][i1,:,:] > 0
    lon_ = X[mask].astype(np.float64)
    lat_ = Y[mask].astype(np.float64)
    data_ = data['MPW_input'][i1,:,:].values[mask]
      
   
    h3_index = vect.geo_to_h3(lat_,lon_,h3_grid_res)
    
    assert(np.in1d(h3_index,h3_index_ocean).shape == np.in1d(h3_index,h3_index_ocean).sum()) #assert that every h3 cell is a valid cell
    
    # h3_index_unique = np.unique(h3_index)
    input_at_h3 = np.zeros(len(h3_index_ocean))

    df_index = pd.Series(index=h3_index_ocean,data=np.arange(len(h3_index_ocean)))

    np.add.at(input_at_h3,df_index[h3_index].values,data_)
    
    data_input_h3_pop[i1,:] = input_at_h3


mask = data_input_h3_pop[0,:]>0
colors = get_colors(np.log10(data_input_h3_pop[0,:][mask]),plt.cm.viridis)

plt.figure()
ax = plt.subplot(111, projection=ccrs.PlateCarree())
pcolorhex(ax,h3_index_ocean[mask],colors,draw_edges=False,alpha=1.)
ax.coastlines(resolution='110m')

to_netcdf('00_data_files/input_h3_coastal_pop.nc', [data_input_h3_pop], ['MPW_input'], h3_index_ocean, data['time'].values)


#%% create riverine input

data = xr.open_dataset('00_data_files/' + input_files[1])
X,Y = np.meshgrid(data['lon'],data['lat'])
data_input_h3_riv_mid = np.zeros([1,len(h3_index_ocean)])
data_input_h3_riv_lo = np.zeros([1,len(h3_index_ocean)])
data_input_h3_riv_hi = np.zeros([1,len(h3_index_ocean)])


mask = data['midpoint'] > 0
lon_ = X[mask].astype(np.float64)
lat_ = Y[mask].astype(np.float64)
data_mid = data['midpoint'].values[mask]
data_lo = data['68per_lower'].values[mask]
data_hi = data['68per_upper'].values[mask]
   
   
h3_index = vect.geo_to_h3(lat_,lon_,h3_grid_res)

assert(np.in1d(h3_index,h3_index_ocean).shape == np.in1d(h3_index,h3_index_ocean).sum()) #assert that every h3 cell is a valid cell

df_index = pd.Series(index=h3_index_ocean,data=np.arange(len(h3_index_ocean)))

np.add.at(data_input_h3_riv_mid[0,:],df_index[h3_index].values,data_mid)
np.add.at(data_input_h3_riv_lo[0,:],df_index[h3_index].values,data_lo)
np.add.at(data_input_h3_riv_hi[0,:],df_index[h3_index].values,data_hi)



mask = data_input_h3_riv_mid[0,:]>0
colors = get_colors(np.log10(data_input_h3_riv_mid[0,:][mask]),plt.cm.viridis)

plt.figure()
ax = plt.subplot(111, projection=ccrs.PlateCarree())
pcolorhex(ax,h3_index_ocean[mask],colors,draw_edges=False,alpha=1.)
ax.coastlines(resolution='110m')


to_netcdf('00_data_files/input_h3_rivers.nc', [data_input_h3_riv_mid,data_input_h3_riv_lo,data_input_h3_riv_hi], ['midpoint','68per_lower','68per_upper'], h3_index_ocean, np.array([np.datetime64('2015-01-01T00:00:00.000000000')]))


#%% create fisheries input

data = xr.open_dataset('00_data_files/' + input_files[2])
X,Y = np.meshgrid(data['lon'],data['lat'])
data_input_h3_fis = np.zeros([12,len(h3_index_ocean)])
# data_input_h3_fis_tot = np.zeros([12,len(h3_index_ocean)])

for i1 in range(len(data['time'])):
    mask = data['fishing_hours_monthly'][i1,:,:] > 0
    lon_ = X[mask].astype(np.float64)
    lat_ = Y[mask].astype(np.float64)
    data_ = data['fishing_hours_monthly'][i1,:,:].values[mask]
      
   
    h3_index = vect.geo_to_h3(lat_,lon_,h3_grid_res)
    
    assert(np.in1d(h3_index,h3_index_ocean).shape == np.in1d(h3_index,h3_index_ocean).sum()) #assert that every h3 cell is a valid cell
    
    # h3_index_unique = np.unique(h3_index)
    # input_at_h3 = np.zeros(len(h3_index_ocean))

    df_index = pd.Series(index=h3_index_ocean,data=np.arange(len(h3_index_ocean)))

    np.add.at(data_input_h3_fis[i1,:],df_index[h3_index].values,data_)
    
    # data_input_h3_fis[i1,:] = input_at_h3


mask = data_input_h3_fis[-1,:]>0
colors = get_colors(np.log10(data_input_h3_fis[-1,:][mask]),plt.cm.viridis)

plt.figure()
ax = plt.subplot(111, projection=ccrs.PlateCarree())
pcolorhex(ax,h3_index_ocean[mask],colors,draw_edges=False,alpha=1.)
ax.coastlines(resolution='110m')

to_netcdf('00_data_files/input_h3_fisheries_monthly.nc', [data_input_h3_fis], ['fishing_hours'], h3_index_ocean, data['time'].values)


#%% create fisheries input, total hrs

data = xr.open_dataset('00_data_files/' + 'input_fisheries.nc')
X,Y = np.meshgrid(data['lon'],data['lat'])
data_input_h3_fis = np.zeros([12,len(h3_index_ocean)])
# data_input_h3_fis_tot = np.zeros([12,len(h3_index_ocean)])

for i1 in range(len(data['time'])):
    mask = data['fishing_hours'][i1,:,:] > 0
    lon_ = X[mask].astype(np.float64)
    lat_ = Y[mask].astype(np.float64)
    data_ = data['fishing_hours'][i1,:,:].values[mask]
      
   
    h3_index = vect.geo_to_h3(lat_,lon_,h3_grid_res)
    
    assert(np.in1d(h3_index,h3_index_ocean).shape == np.in1d(h3_index,h3_index_ocean).sum()) #assert that every h3 cell is a valid cell
    
    # h3_index_unique = np.unique(h3_index)
    # input_at_h3 = np.zeros(len(h3_index_ocean))

    df_index = pd.Series(index=h3_index_ocean,data=np.arange(len(h3_index_ocean)))

    np.add.at(data_input_h3_fis[i1,:],df_index[h3_index].values,data_)
    
    # data_input_h3_fis[i1,:] = input_at_h3


mask = data_input_h3_fis[-1,:]>0
colors = get_colors(np.log10(data_input_h3_fis[-1,:][mask]),plt.cm.viridis)

plt.figure()
ax = plt.subplot(111, projection=ccrs.PlateCarree())
pcolorhex(ax,h3_index_ocean[mask],colors,draw_edges=False,alpha=1.)
ax.coastlines(resolution='110m')

to_netcdf('00_data_files/input_h3_fisheries_total_hrs.nc', [data_input_h3_fis], ['fishing_hours'], h3_index_ocean, data['time'].values)

