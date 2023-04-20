#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 14:39:52 2021

@author: kaandorp
"""
import numpy as np
import os
import glob
# from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from datetime import timedelta,datetime
import xarray as xr
import pandas as pd
# from scipy.interpolate import griddata
# import shapely
# import cartopy.io.shapereader as shpreader
# import shapefile
# import csv
import matplotlib
# from time import sleep
# from parcels import rng as random
# import math
# from parcels.tools.converters import Geographic, GeographicPolar 
import cmocean
from argparse import ArgumentParser
# import parcels.rng as ParcelsRandom
from tools import NEMO_select_section 

def to_netcdf(output_filename,data,data_name,lons,lats,explanation=''):
    '''
    All data is written to netcdf files to speed up computations
    '''
    
    dict_data = {}
    for data_, name_ in zip(data,data_name):
        dict_data[name_] = (["x","y"], data_)
    dict_data['explanation'] = explanation
        
    dict_coords = {}
    dict_coords['lon'] = (["x","y"], lons)
    dict_coords['lat'] = (["x","y"], lats)   
  
    ds = xr.Dataset(
        dict_data,
        coords=dict_coords,
    )   
    ds.to_netcdf(output_filename)


def get_mask_land(field,outfile='./tmp_mask_land'):
    '''
    Mask with true on land cells, false on ocean cells
    '''
    if os.path.exists(outfile):
        ds = xr.open_dataset(outfile)
        mask_land = np.array(ds['mask_land'],dtype=bool)    
    else:
        mask_land = np.isnan(field)
        to_netcdf(outfile,[mask_land],['mask_land'],lons,lats,explanation='land mask')
    return mask_land


# def get_mask_landborder(mask_land,outfile='./tmp_mask_landborder'):
#     '''
#     calculate the landborder mask. With landborder cells, we mean cells on land, adjacent to water
#     '''
#     if os.path.exists(outfile):
#         ds = xr.open_dataset(outfile)
#         mask_landborder = np.array(ds['mask_landborder'],dtype=bool)    
#     else:
#         #check the upper,lower,left & right neighbor: if one of these is an ocean cell, set to landborder
#         mask_landborder = mask_land & (np.roll(~mask_land,1,axis=0) | np.roll(~mask_land,-1,axis=0) | 
#                                        np.roll(~mask_land,1,axis=1) | np.roll(~mask_land,-1,axis=1))
#         to_netcdf(outfile,[mask_landborder],['mask_landborder'],lons,lats,explanation='landborder mask (land cells adjacent to water')
#     return mask_landborder 


# def get_mask_coast(mask_land,outfile='./tmp_mask_coast'):
#     '''
#     calculate the coast mask. With coastal cells, we mean cells in the water, adjacent to land
#     '''
#     if os.path.exists(outfile):
#         ds = xr.open_dataset(outfile)
#         mask_coast = np.array(ds['mask_coast'],dtype=bool)    
#     else:
#         #check the upper,lower,left & right neighbor: if one of these is an ocean cell, set to landborder
#         mask_coast = ~mask_land & (np.roll(mask_land,1,axis=0) | np.roll(mask_land,-1,axis=0) | 
#                                        np.roll(mask_land,1,axis=1) | np.roll(mask_land,-1,axis=1))
#         to_netcdf(outfile,[mask_coast],['mask_coast'],lons,lats,explanation='coast mask (ocean cells adjacent to land')
#     return mask_coast


# def extend_landborder(mask_landborder,mask_land,n_iter=1):
#     '''
#     calculate the landborder mask. With landborder cells, we mean cells on land, adjacent to water
#     '''
#     mask_landborder_extended = mask_landborder.copy()
#     for i1 in range(n_iter):
#         mask_landborder_extended = mask_land & (np.roll(mask_landborder_extended,1,axis=0) | np.roll(mask_landborder_extended,-1,axis=0) | 
#                                            np.roll(mask_landborder_extended,1,axis=1) | np.roll(mask_landborder_extended,-1,axis=1))
#     return mask_landborder_extended


# def get_land_current(mask_landborder,mask_coast,mask_land,n_iter=1,outfile='./tmp_land_current'):
#     '''
#     Calculate current which pushes particles back from the land into the sea.
#     The scripts adds 1's and -1's to locations which are on the landborder, and adjacent to sea
#     depending on whether the sea is on the left/right/upper/lower cell
    
#     Scipt can be iterated for some safety margin in case particles jump more than one cell
#     '''
#     if os.path.exists(outfile):
#         ds = xr.open_dataset(outfile)
#         land_current_u = np.array(ds['land_current_u'],dtype=float)    
#         land_current_v = np.array(ds['land_current_v'],dtype=float)    
#     else:
#         land_current_u = np.zeros([len(lats),len(lons)])
#         land_current_v = np.zeros([len(lats),len(lons)])
        
#         # checks on which side the sea is located, e.g. True when there is a sea cell on the left
#         # in which case a u-velocity of -1 is added
#         mask_left = mask_landborder & np.roll(mask_coast,1,axis=1)
#         mask_right = mask_landborder & np.roll(mask_coast,-1,axis=1)
#         mask_lower = mask_landborder & np.roll(mask_coast,1,axis=0)
#         mask_upper = mask_landborder & np.roll(mask_coast,-1,axis=0)
    
#         land_current_u[mask_left] -= 1.
#         land_current_u[mask_right] += 1.
#         land_current_v[mask_lower] -= 1.
#         land_current_v[mask_upper] += 1.
    
#         if n_iter > 1:
#             mask_landborder_ = mask_landborder.copy()
#             for i1 in range(1,n_iter):
#                 # we will iterate though land cells which are adjacent to the previous land border
#                 # this means that the new cells should be land, be adjacent to the previous mask_landborder
#                 #(mask_landborder_), and not already been included previously as a landborder cell
#                 new_landborder = (mask_land 
#                                   & (np.roll(mask_landborder_,1,axis=0) | np.roll(mask_landborder_,-1,axis=0) | 
#                                                np.roll(mask_landborder_,1,axis=1) | np.roll(mask_landborder_,-1,axis=1))
#                                   & ~mask_landborder_)
               
#                 mask_left = new_landborder & np.roll(mask_landborder_,1,axis=1)
#                 mask_right = new_landborder & np.roll(mask_landborder_,-1,axis=1)
#                 mask_lower = new_landborder & np.roll(mask_landborder_,1,axis=0)
#                 mask_upper = new_landborder & np.roll(mask_landborder_,-1,axis=0)            
    
#                 land_current_u[mask_left] -= 1.
#                 land_current_u[mask_right] += 1.
#                 land_current_v[mask_lower] -= 1.
#                 land_current_v[mask_upper] += 1.
    
#                 mask_landborder_ = (mask_landborder_ | new_landborder)
    
#         land_current_mag = np.sqrt(land_current_u**2 + land_current_v**2)
#         land_current_u /= land_current_mag
#         land_current_v /= land_current_mag
        
#         to_netcdf(outfile,[land_current_u,land_current_v],['land_current_u','land_current_v'],lons,lats,explanation='land current, pusing particles on land back to the sea, magnitude of 1')
        
#     return land_current_u, land_current_v


def get_coastal_nodes(landmask):
    """Function that detects the coastal nodes, i.e. the ocean nodes directly
    next to land. Computes the Laplacian of landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the coastal nodes, the coastal nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap -= 4*landmask
    coastal = np.ma.masked_array(landmask, mask_lap > 0)
    coastal = coastal.mask.astype('int')

    return coastal


def get_shore_nodes(landmask):
    """Function that detects the shore nodes, i.e. the land nodes directly
    next to the ocean. Computes the Laplacian of landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the shore nodes, the shore nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap -= 4*landmask
    shore = np.ma.masked_array(landmask, mask_lap < 0)
    shore = shore.mask.astype('int')

    return shore


def get_coastal_nodes_diagonal(landmask):
    """Function that detects the coastal nodes, i.e. the ocean nodes where 
    one of the 8 nearest nodes is land. Computes the Laplacian of landmask
    and the Laplacian of the 45 degree rotated landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the coastal nodes, the coastal nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap += np.roll(landmask, (-1,1), axis=(0,1)) + np.roll(landmask, (1, 1), axis=(0,1))
    mask_lap += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (1, -1), axis=(0,1))
    mask_lap -= 8*landmask
    coastal = np.ma.masked_array(landmask, mask_lap > 0)
    coastal = coastal.mask.astype('int')
    
    return coastal
    

def get_shore_nodes_diagonal(landmask):
    """Function that detects the shore nodes, i.e. the land nodes where 
    one of the 8 nearest nodes is ocean. Computes the Laplacian of landmask 
    and the Laplacian of the 45 degree rotated landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the shore nodes, the shore nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap += np.roll(landmask, (-1,1), axis=(0,1)) + np.roll(landmask, (1, 1), axis=(0,1))
    mask_lap += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (1, -1), axis=(0,1))
    mask_lap -= 8*landmask
    shore = np.ma.masked_array(landmask, mask_lap < 0)
    shore = shore.mask.astype('int')

    return shore


def create_displacement_field(landmask, double_cell=False, outfile='./tmp_mask_coast'):
    """Function that creates a displacement field 1 m/s away from the shore.

    - landmask: the land mask dUilt using `make_landmask`.
    - double_cell: Boolean for determining if you want a double cell.
      Default set to False.

    Output: two 2D arrays, one for each camponent of the velocity.
    """
    
    if os.path.exists(outfile):
        ds = xr.open_dataset(outfile)
        v_x = np.array(ds['land_current_u'],dtype=float)    
        v_y = np.array(ds['land_current_v'],dtype=float)   
    
    else:
        shore = get_shore_nodes(landmask)
        shore_d = get_shore_nodes_diagonal(landmask) # bordering ocean directly and diagonally
        shore_c = shore_d - shore                    # corner nodes that only border ocean diagonally
        
        Ly = np.roll(landmask, -1, axis=0) - np.roll(landmask, 1, axis=0) # Simple derivative
        Lx = np.roll(landmask, -1, axis=1) - np.roll(landmask, 1, axis=1)
        
        Ly_c = np.roll(landmask, -1, axis=0) - np.roll(landmask, 1, axis=0)
        Ly_c += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (-1,1), axis=(0,1)) # Include y-component of diagonal neighbours
        Ly_c += - np.roll(landmask, (1,-1), axis=(0,1)) - np.roll(landmask, (1,1), axis=(0,1))
        
        Lx_c = np.roll(landmask, -1, axis=1) - np.roll(landmask, 1, axis=1)
        Lx_c += np.roll(landmask, (-1,-1), axis=(1,0)) + np.roll(landmask, (-1,1), axis=(1,0)) # Include x-component of diagonal neighbours
        Lx_c += - np.roll(landmask, (1,-1), axis=(1,0)) - np.roll(landmask, (1,1), axis=(1,0))
        
        v_x = -Lx*(shore)
        v_y = -Ly*(shore)
        
        v_x_c = -Lx_c*(shore_c)
        v_y_c = -Ly_c*(shore_c)
        
        v_x = v_x + v_x_c
        v_y = v_y + v_y_c
    
        magnitude = np.sqrt(v_y**2 + v_x**2)
        # the coastal nodes between land create a problem. Magnitude there is zero
        # I force it to be 1 to avoid problems when normalizing.
        ny, nx = np.where(magnitude == 0)
        magnitude[ny, nx] = 1
    
        v_x = v_x/magnitude
        v_y = v_y/magnitude

        to_netcdf(outfile,[v_x,v_y],['land_current_u','land_current_v'],lons,lats,explanation='land current, pusing particles on land back to the sea, magnitude of 1')
    return v_x, v_y


#%%
            
data_T = xr.open_dataset('/Volumes/externe_SSD/kaandorp/Data/NEMO-MEDUSA/ORCA0083-N006/means/ORCA0083-N06_20100125d05T.nc')

data_grid = xr.open_dataset('/Volumes/externe_SSD/kaandorp/Data/NEMO-MEDUSA/ORCA0083-N006/domain/coordinates.nc',decode_times=False)


sst = data_T['sst'].values[0,:,:]
lons = data_grid['nav_lon'][:,:].values
lats = data_grid['nav_lat'][:,:].values
# lat[1494,:] = 0

# lons = data_currents['longitude'].values
# lats = data_currents['latitude'].values
# dlon = lons[1] - lons[0]
# dlat = lats[1] - lats[0]
# fieldMesh_x,fieldMesh_y = np.meshgrid(lons,lats)


mask_land = get_mask_land(sst,outfile='00_data_files/mask_land_NEMO0083.nc')

v_x, v_y = create_displacement_field(mask_land*1, outfile='00_data_files/land_current_NEMO0083.nc')

# mask_landborder = get_mask_landborder(mask_land,outfile='00_data_files/mask_landborder.nc')

# mask_coast = get_mask_coast(mask_land,outfile='00_data_files/mask_coast.nc')

# mask_landborder_extended = extend_landborder(mask_landborder, mask_land, 3)

# land_current_u, land_current_v = get_land_current(mask_landborder,mask_coast,mask_land,n_iter=3,outfile='00_data_files/land_current.nc')
 
# coastal_distance, index_coast_dlon, index_coast_dlat = get_coastal_distance(mask_landborder,mask_coast,mask_land,n_iter=30,outfile='00_data_files/coastal_distance.nc')



extent = [-10,20,30,40]
X_,Y_,mask_land_ = NEMO_select_section(extent,lons,lats,mask_land)
X_,Y_,v_x_ = NEMO_select_section(extent,lons,lats,v_x)
X_,Y_,v_y_ = NEMO_select_section(extent,lons,lats,v_y)

plt.figure()
plt.contourf(X_,Y_,mask_land_)
plt.quiver(X_,Y_,v_x_,v_y_,scale=10)

