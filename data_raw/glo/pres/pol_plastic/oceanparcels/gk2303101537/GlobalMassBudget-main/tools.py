#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 11:40:25 2021

@author: kaandorp
"""

import numpy as np
import h3
import cartopy.crs as ccrs
import pickle
import matplotlib.pyplot as plt

def pcolorhex(ax,hexagons,colors,draw_edges=True,fill_polygons=True,transform=ccrs.PlateCarree(),extent=None,dict_h3_location=None,**kwargs):
  
    # if extent:
    #     mask = (dict_h3_location['lon']>extent[0]) & (dict_h3_location['lon']<extent[1]) & (dict_h3_location['lat']>extent[2]) & (dict_h3_location['lat']<extent[3])
    #     hexagons = hexagons[mask]
    #     colors = colors[mask]
        
    for i1,hex_ in enumerate(hexagons):
        x = np.array(h3.h3_to_geo_boundary(hex(hex_)))[:,1]
        y = np.array(h3.h3_to_geo_boundary(hex(hex_)))[:,0]
        
        x_hexagon = np.append(x,x[0]) 
        y_hexagon = np.append(y,y[0])
         
        if x_hexagon.max() - x_hexagon.min() > 25:
            x_hexagon[x_hexagon < 0] += 360
            
        if draw_edges:
            ax.plot(x_hexagon, y_hexagon, 'k-',transform=transform,linewidth=.2)
        
        if fill_polygons:
            ax.fill(x_hexagon,y_hexagon,color=colors[i1],transform=transform,**kwargs)    


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))        

def get_colors_tableau():
    return ['#006BA4', '#FF800E', '#ABABAB', '#595959',
                     '#5F9ED1', '#C85200', '#898989', '#A2C8EC', '#FFBC79', '#CFCFCF']


def to_pickle(item,filename):
    outfile = open(filename,'wb')
    pickle.dump(item,outfile)
    outfile.close()   
    
  
def from_pickle(filename):
    with open(filename, 'rb') as f:
        item = pickle.load(f)
    return item


def to_cell_edges(vals):
    assert(np.all(vals[1:]-vals[:-1] == vals[1]-vals[0] ))    
    dval = vals[1]-vals[0]
    vals_edges = .5*(vals[1:]+vals[:-1])
    vals_edges = np.append(vals_edges[0]-dval,vals_edges)
    vals_edges = np.append(vals_edges,vals_edges[-1]+dval)
    
    return vals_edges


def NEMO_select_section(extent,lon,lat,val):
    
    lon_mean = .5*(extent[0]+extent[1])
    lat_mean = .5*(extent[2]+extent[3])
    
    i_min = np.unravel_index( (np.sqrt( (lon-lon_mean)**2 + (lat-lat_mean)**2   )).argmin() , lon.shape)
    
    i_lon_s = np.where(  (lon[i_min[0],:] > extent[0]) & (lon[i_min[0],:] < extent[1]) )[0][0]
    i_lon_e = np.where(  (lon[i_min[0],:] > extent[0]) & (lon[i_min[0],:] < extent[1]) )[0][-1]

    i_lat_s = np.where(  (lat[:,i_min[1]] > extent[2]) & (lat[:,i_min[1]] < extent[3]) )[0][0]
    i_lat_e = np.where(  (lat[:,i_min[1]] > extent[2]) & (lat[:,i_min[1]] < extent[3]) )[0][-1]
    
    return lon[i_lat_s:i_lat_e,i_lon_s:i_lon_e], lat[i_lat_s:i_lat_e,i_lon_s:i_lon_e], val[i_lat_s:i_lat_e,i_lon_s:i_lon_e]
    
    

    