# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import xarray as xr
import geopandas as gpd


def to_netcdf(output_filename,data,data_name,lons,lats,explanation=''):
    '''
    All data is written to netcdf files to speed up computations
    '''
    dict_data = {}
    for name_, data_ in zip(data_name, data):
        dict_data[name_] = (( "lat", "lon"), data_ )
    dict_data['explanation'] = explanation
        
    ds = xr.Dataset(
        dict_data,
        coords={
            "lon": lons,
            "lat": lats,
        },
    )   
    ds.to_netcdf(output_filename)


def find_closest_coast(lon, lat):
   
    distance = np.sqrt( ( (lon_coast - lon)*np.cos(lat * (np.pi/180)))**2 + (lat_coast- lat)**2 )

    return indices_lon[mask_coast][np.argmin(distance)], indices_lat[mask_coast][np.argmin(distance)]


river_data = '/Users/kaandorp/Data/PlasticData/PlasticRiverInputs_Meijer/Meijer2021_midpoint_emissions/'

data_mask_coast = xr.open_dataset('00_data_files/mask_coast.nc')
mask_coast = data_mask_coast['mask_coast'].values

lons = data_mask_coast['lon']
lats = data_mask_coast['lat']
fieldMesh_x,fieldMesh_y = np.meshgrid(lons,lats)


data = gpd.read_file(river_data)
lon_river = data.geometry.x
lat_river = data.geometry.y
output_river = data['dots_exten'].values #in tonnes per year

input_matrix_rivers = np.zeros(fieldMesh_x.shape)

lon_coast = fieldMesh_x[mask_coast]
lat_coast = fieldMesh_y[mask_coast]
indices_lon,indices_lat = np.meshgrid(np.arange(len(lons)),np.arange(len(lats)))


for i1, (lon_,lat_,output_) in enumerate(zip(lon_river,lat_river,output_river)):
    
    i_lon, i_lat = find_closest_coast(lon_,lat_)    
    
    input_matrix_rivers[i_lat,i_lon] += output_
    
    if i1 % 100 == 0:
        print('%i/%i' % (i1,len(lon_river)))
        
        
to_netcdf('00_data_files/input_rivers.nc',[input_matrix_rivers,input_matrix_rivers/4,input_matrix_rivers*4],['midpoint','68per_lower','68per_upper'],lons,lats,explanation='River input estimates Meijer et al. (2021), midpoints and lower/upper estimates based on 68 percent C.I.')        
    

