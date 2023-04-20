#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 14:55:14 2021

@author: kaandorp
"""
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import matplotlib.colors as colors
from datetime import datetime
import numpy.ma as ma
import cmocean

def to_netcdf(output_filename,data,data_name,lons,lats,time,explanation='',encode=False):
    '''
    All data is written to netcdf files to speed up computations
    '''
    dict_data = {}
    for name_, data_ in zip(data_name, data):
        dict_data[name_] = (( "time","lat", "lon"), data_ )
    dict_data['explanation'] = explanation
        
    ds = xr.Dataset(
        dict_data,
        coords={
            "lon": lons,
            "lat": lats,
            "time": time,
        },
    )   
    if encode:
        ds.to_netcdf(output_filename,encoding={name_: {'dtype': 'int16', 'scale_factor': 0.1, '_FillValue': -9999}})
    else:
        ds.to_netcdf(output_filename)
        
mask_land = xr.open_dataset('00_data_files/mask_land.nc')['mask_land']
mask_land_nan = mask_land.copy().values.astype(float)
mask_land_nan[~mask_land] = np.nan
lons,lats = mask_land['lon'].values, mask_land['lat'].values
fieldMesh_x,fieldMesh_y = np.meshgrid(lons,lats)
dlon = lons[1] - lons[0]
dlat = lats[1] - lats[0]
lons_edge = lons - .5*dlon
lats_edge = lats - .5*dlat
lons_edge = np.append(lons_edge,lons_edge[-1]+dlon)
lats_edge = np.append(lats_edge,lats_edge[-1]+dlat)
mesh_plot_x,mesh_plot_y = np.meshgrid(lons_edge,lats_edge)

# files_fisheries = sorted(glob.glob('/Users/kaandorp/Data/Fisheries/fleet*/**'))
files_fisheries = sorted(glob.glob('/Volumes/externe_SSD/kaandorp/Data/Fisheries/fleet*/**') )
                         
dlon_f,dlat_f = 0.01,0.01

date_array = pd.date_range(datetime(2012,1,1),datetime(2021,1,1),freq='MS')    
delta_date_array = (( np.array(date_array) - np.datetime64(datetime(1970,1,1)) ) ).astype('timedelta64[s]').astype(int)

input_fisheries = np.zeros([len(date_array)-1,fieldMesh_x.shape[0],fieldMesh_x.shape[1]])

for file_ in files_fisheries:
    
    date = datetime.fromisoformat(os.path.basename(file_).split('.')[0])
    if date.day == 1:
        print(date)
    delta_date = ((np.datetime64(date) - np.datetime64(datetime(1970,1,1))) ).astype('timedelta64[s]').astype(int)
    
    data_fisheries = pd.read_csv(file_)

    data_fishing_hours = data_fisheries[data_fisheries['fishing_hours'] > 0]

    i_lat = np.digitize(data_fishing_hours['cell_ll_lat']+.5*dlat_f, lats_edge)-1
    i_lon = (np.digitize(data_fishing_hours['cell_ll_lon']+.5*dlon_f, lons_edge)-1) % len(lons)
    i_time = (np.digitize(delta_date,delta_date_array)-1) * np.ones(len(i_lat),dtype=int)
    
    input_fisheries[i_time,i_lat,i_lon] += data_fishing_hours['fishing_hours'].values


plt.figure()
plt.pcolormesh(mesh_plot_x,mesh_plot_y,mask_land_nan)
plt.pcolormesh(mesh_plot_x,mesh_plot_y,input_fisheries[0,:,:],norm=colors.LogNorm(vmin=1e-3,vmax=1e2))

input_fisheries[:,mask_land] = 0

SAVE=False
if SAVE:
    to_netcdf('00_data_files/input_fisheries.nc',[input_fisheries],['fishing_hours'],lons,lats,date_array[:-1],explanation='Total fishing hours per month, calculated from the Global Fishing Watch dataset')

# data = xr.open_dataset('00_data_files/input_fisheries.nc')
# input_fisheries = data['fishing_hours'].values


n_years = 9
hours_per_month = input_fisheries.sum(axis=(1,2))
hours_per_year = np.array([np.sum(hours_per_month[i1*12:i1*12+12]) for i1 in range(n_years)])
input_fisheries_normalized = np.copy(input_fisheries)

for i1 in range(n_years):
    input_fisheries_normalized[i1*12:i1*12+12,:,:] /= hours_per_year[i1]

# calculate the monthly fishing intensity: 
    #every year is normalized to have a sum of 1
    #afterwards, the fishing hours per month for each year are summed up (2012 is left out)
    #the monthly fishing hours are normalized to a sum of 1
input_fisheries_monthly = np.zeros([12,fieldMesh_x.shape[0],fieldMesh_x.shape[1]])
for i1 in range(12):
    input_fisheries_monthly[i1,:,:] = input_fisheries_normalized[i1::12,:,:][1:,:,:].sum(axis=0) #1:,:,: as 2012 is left out
input_fisheries_monthly /= input_fisheries_monthly.max()


if SAVE:
    to_netcdf('00_data_files/input_fisheries_monthly.nc',[input_fisheries_monthly],['fishing_hours_monthly'],lons,lats,
              pd.date_range(datetime(2012,1,1),datetime(2012,12,1),freq='MS'),explanation='Normalized fishing hours per month, calculated from the Global Fishing Watch dataset over 2012-2020',
              encode=False)


fig,ax = plt.subplots(1)
ax.plot(pd.date_range(datetime(2012,1,1),datetime(2012,12,1),freq='MS').month,input_fisheries_monthly.sum(axis=(1,2)),'ko-')
ax.set_xlabel('Month')
ax.set_ylabel('Fishing intensity [-]')

fig,ax = plt.subplots(1)
ax.plot(np.arange(2012,2012+n_years),hours_per_year,'ko-')
ax.set_xticks(np.arange(2012,2012+n_years))
ax.set_xlabel('year')
ax.set_ylabel('Total fishing hours per year')

fig,ax = plt.subplots(1)
ax.plot(date_array[:-1],input_fisheries.sum(axis=(1,2)))
ax.set_xlabel('date')
ax.set_ylabel('Total fishing hours per month')

cmap = plt.cm.viridis
cmap.set_under('k')
cmap.set_bad('k')
str_m = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
for i1 in range(12):
    print(i1)
    fig,ax = plt.subplots(1,dpi=120,figsize=(10,5))
    # ax.pcolormesh(mesh_plot_x,mesh_plot_y,~mask_land,cmap=cmocean.cm.deep_r,vmin=.95,vmax=1.8)
    cplt = ax.pcolormesh(mesh_plot_x,mesh_plot_y,input_fisheries_monthly[i1,:,:],cmap=cmap,norm=colors.LogNorm(vmin=1e-5,vmax=1e0))
    ax.pcolormesh(mesh_plot_x,mesh_plot_y,mask_land_nan,cmap=cmocean.cm.turbid,vmin=.95,vmax=1.8)
    cbar = plt.colorbar(cplt,extend='min')
    cbar.set_label('Fishing intensity [-]')
    ax.set_title('%s' % str_m[i1])
    fig.savefig('00_figures/%4.4i.png' % i1)
    plt.close(fig)