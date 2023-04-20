#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 16:02:54 2021

@author: kaandorp
"""
import numpy as np
import xarray as xr
import pandas as pd
from scipy.interpolate import griddata
import shapely
import cartopy.io.shapereader as shpreader


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
    
def list_contains_string(string_,list_):
    result = False
    for element in list_:
        if string_ in element:
            result = True
    return result

def list_contains_string_index(string1,string2,list_):
    result = None
    for i1,element in enumerate(list_):
        if string1 in element or string2 in element:
            result = i1
    return result    

shpfilename = shpreader.natural_earth(resolution='50m',
                                      category='cultural',
                                      name='admin_0_countries')

data_mpw = pd.read_excel('00_data_files/plastic_rivers2sea_v4.xlsx',sheet_name='Mismanaged Plastic Waste')

data_mask_land = xr.open_dataset('00_data_files/mask_land.nc')

mask_land = data_mask_land['mask_land'].values
lons = data_mask_land['lon'].values
lats = data_mask_land['lat'].values
fieldMesh_x,fieldMesh_y = np.meshgrid(lons,lats)


all_country_names = list(data_mpw['Country'])


mpw_mat = np.zeros(fieldMesh_x.shape)
country_names = []
mpw_ = []


reader = shpreader.Reader(shpfilename)
countries = reader.records()
for country in countries:

    country_name = country.attributes['NAME_LONG']
    country_name_2 = country.attributes['NAME']
    
    geom = country.geometry
    
    if list_contains_string(country_name_2,all_country_names) or list_contains_string(country_name,all_country_names):
       
        country_names.append(country_name)

        mpw = np.nan
        
        i1 = list_contains_string_index(country_name,country_name_2,all_country_names)
        country_mpw = data_mpw.loc[i1,'Country']
        mpw = data_mpw.iloc[i1,6]
       
        print(country_mpw,mpw,i1)
         
        if np.isnan(mpw):
            print(mpw)
        
        # create initial matrix containing mpw
        for i1 in range(len(lons)):
            for i2 in range(len(lats)):
                xy_point = shapely.geometry.Point(fieldMesh_x[i2,i1],fieldMesh_y[i2,i1]) 
                
                if geom.contains(xy_point):
                    mpw_mat[i2,i1] = mpw

    else:
        print('Not found:')
        print(country_name_2,country_name)
    
                



#create refined matrix using the given landmask   
nanmask = np.isnan(mpw_mat)
mpw_mat_final = np.zeros(fieldMesh_x.shape)

for i1 in range(len(lons)):
    for i2 in range(len(lats)):
        
        if mask_land[i2,i1]:
            
            fillVal = mpw_mat[i2,i1]
            
            if np.isnan(fillVal):
                fillVal = griddata((fieldMesh_x[~nanmask],fieldMesh_y[~nanmask]),mpw_mat[~nanmask],(fieldMesh_x[i2,i1],fieldMesh_y[i2,i1]))
                
            mpw_mat_final[i2,i1] = fillVal


to_netcdf('00_data_files/mpw_gridded.nc',[mpw_mat,mpw_mat_final],['mpw_gridded','mpw_gridded_filled'],lons,lats,explanation='Gridded MPW data based on Jambeck et al. (2015), raw data, and data where nan values are interpolated')        
