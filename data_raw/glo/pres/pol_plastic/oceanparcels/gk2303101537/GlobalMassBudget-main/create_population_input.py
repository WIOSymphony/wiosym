#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 16:02:54 2021

@author: kaandorp
"""

import numpy as np
import os
import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import xarray as xr
import pandas as pd
from scipy.interpolate import griddata
import shapely
import cartopy.io.shapereader as shpreader
import matplotlib.colors as colors
from datetime import datetime

def to_netcdf(output_filename,data,data_name,lons,lats,time,explanation=''):
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
    ds.to_netcdf(output_filename)
    


data_mpw = xr.open_dataset('00_data_files/mpw_gridded.nc')

lons = data_mpw['lon']
lats = data_mpw['lat']
fieldMesh_x,fieldMesh_y = np.meshgrid(lons,lats)
dlon = lons[1] - lons[0]
dlat = lats[1] - lats[0]
lons_plot = lons - .5*dlon
lats_plot = lats - .5*dlat
lons_plot = np.append(lons_plot,lons_plot[-1]+dlon)
lats_plot = np.append(lats_plot,lats_plot[-1]+dlat)

mesh_plot_x,mesh_plot_y = np.meshgrid(lons_plot,lats_plot)

mpw_gridded = data_mpw['mpw_gridded_filled']

pop_density_data = xr.open_dataset('/Users/kaandorp/Data/populationDensity/gpw_v4_e_atotpopbt_dens_2pt5_min.nc')

coastal_distance = xr.open_dataset('00_data_files/coastal_distance.nc')
mask_land = xr.open_dataset('00_data_files/mask_land.nc')['mask_land']
mask_50km = (coastal_distance['total_distance'] < 50) & (mask_land)

mask_land_nan = mask_land.copy().values.astype(float)
mask_land_nan[~mask_land] = np.nan

# relative closest coastal points
index_closest_u = coastal_distance['index_closest_u']
index_closest_v = coastal_distance['index_closest_v']

# 2 matrices containing the absolute (instead of relative) indices of the closest coastal points
index_mat_u,index_mat_v = np.meshgrid(np.arange(len(lons),dtype=float),np.arange(len(lats),dtype=float))

index_mat_u[mask_50km] += index_closest_u.values[mask_50km]
index_mat_v[mask_50km] += index_closest_v.values[mask_50km]

index_mat_u = (index_mat_u % len(lons)) #periodic BC's
index_mat_v = (index_mat_v % len(lats))


input_matrix = np.zeros([5,fieldMesh_x.shape[0],fieldMesh_x.shape[1]])

MPW_matrix = np.zeros([5,fieldMesh_x.shape[0],fieldMesh_x.shape[1]])
popden_matrix = np.zeros([5,fieldMesh_x.shape[0],fieldMesh_x.shape[1]])
#area of each grid cell. km's in lon is multiplied by km's in lat, scaled by the reduction in latitude
grid_area = ((dlon*1.11e2) * (dlat*1.11e2)).values * np.ones(fieldMesh_x.shape) * np.cos(fieldMesh_y * (np.pi/180))

for i1 in range(5):

    print(i1)
    pop_density = pop_density_data['Population Density, v4.10 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][i1,:,:]
    X_pop,Y_pop = np.meshgrid(pop_density['longitude'],pop_density['latitude'])
  
    print('Interpolating data...')
    pop_density_i = griddata((X_pop.ravel(),Y_pop.ravel()),pop_density.values.ravel(),(fieldMesh_x,fieldMesh_y),method='nearest')
    
    #               kg/pop/day    pop/km2         km2 -> kg/day per grid cell
    MPW_density = mpw_gridded * pop_density_i * grid_area
    
    MPW_matrix[i1,:,:] = MPW_density
    popden_matrix[i1,:,:] = pop_density_i
    
    MPW_density_50km = MPW_density.values.copy()
    values_50km_zero = ( (np.isnan(MPW_density_50km) | (MPW_density_50km==0)) & mask_50km)
    
    #filling in missing data using nearest neighbor
    MPW_density_50km[values_50km_zero] = griddata((fieldMesh_x[~values_50km_zero],fieldMesh_y[~values_50km_zero]),MPW_density_50km[~values_50km_zero],(fieldMesh_x[values_50km_zero],fieldMesh_y[values_50km_zero]),method='nearest')
    
    MPW_density_50km[~mask_50km] = 0
    MPW_density_50km[np.isnan(MPW_density_50km)] = 0
    # mask_nan_50km = (MPW_density_50km[mask_50km] == np.nan)

    # if mask_nan_50km.sum() > 0:
    #     print(mask_nan_50km.sum())
        
    plt.figure()
    plt.pcolormesh(mesh_plot_x,mesh_plot_y,MPW_density_50km,norm=colors.LogNorm(vmin=1e-5,vmax=1e3))
        
    index_array_closest_u = index_mat_u[mask_50km].astype(int)
    index_array_closest_v = index_mat_v[mask_50km].astype(int)
    
    np.add.at(input_matrix[i1,:,:],(index_array_closest_v,index_array_closest_u),MPW_density_50km[mask_50km])
    
    plt.figure()
    plt.pcolormesh(mesh_plot_x,mesh_plot_y,mask_land_nan)
    plt.pcolormesh(mesh_plot_x,mesh_plot_y,input_matrix[i1,:,:],norm=colors.LogNorm(vmin=1e-5,vmax=1e3))
    
    print(input_matrix[i1,:,:].sum())
    
date_array = [datetime(2000,1,1),datetime(2005,1,1),datetime(2010,1,1),datetime(2015,1,1),datetime(2020,1,1)] #.date_range(datetime(2000,1,1),datetime(2020,1,1),periods=5)    

to_netcdf('00_data_files/input_coastal_pop.nc',[input_matrix],['MPW_input'],lons,lats,date_array,explanation='Plastic input from 50km coastal population (assuming 100% MPW enters the sea)')

to_netcdf('00_data_files/MPW_popden_gridded.nc',[MPW_matrix,popden_matrix],['MPW_input','population_density'],lons,lats,date_array,explanation='Gridded MPW data [kg/day] per gridcell, and population density data [pop/km2], converted to the same grid')


