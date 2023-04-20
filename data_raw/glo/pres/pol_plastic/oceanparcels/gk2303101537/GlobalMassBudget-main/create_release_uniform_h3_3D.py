#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 14:30:13 2021

@author: kaandorp
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from h3.unstable import vect
import h3
import cartopy.crs as ccrs
from scipy.interpolate import griddata
import os

def to_netcdf(output_filename,data,data_name):
    '''
    Create particle netcdf file
    '''

    dict_data = {}
    for data_,data_name_ in zip(data,data_name):
        dict_data[data_name_] = (["traj","obs"], data_)
        

    dict_attrs = {'feature_type' : 'trajectory',
                  'Conventions' : 'CF-1.6/CF-1.7',
                  'ncei_template_version':  'NCEI_NetCDF_Trajectory_Template_v2.0',
                  'parcels_version':        '2.2.3.dev7+g8df4758b',
                  'parcels_mesh':           'spherical'}
  
    ds = xr.Dataset(
        dict_data,
        attrs=dict_attrs,
    )   
    ds.to_netcdf(output_filename)
    
def pcolormesh_NEMO(extent,lon,lat,val,**kwargs):
    
    lon_mean = .5*(extent[0]+extent[1])
    lat_mean = .5*(extent[2]+extent[3])
    
    i_min = np.unravel_index( (np.sqrt( (lon-lon_mean)**2 + (lat-lat_mean)**2   )).argmin() , lon.shape)
    
    i_lon_s = np.where(  (lon[i_min[0],:] > extent[0]) & (lon[i_min[0],:] < extent[1]) )[0][0]
    i_lon_e = np.where(  (lon[i_min[0],:] > extent[0]) & (lon[i_min[0],:] < extent[1]) )[0][-1]

    i_lat_s = np.where(  (lat[:,i_min[1]] > extent[2]) & (lat[:,i_min[1]] < extent[3]) )[0][0]
    i_lat_e = np.where(  (lat[:,i_min[1]] > extent[2]) & (lat[:,i_min[1]] < extent[3]) )[0][-1]
    
    ax.pcolormesh(lon[i_lat_s:i_lat_e,i_lon_s:i_lon_e], lat[i_lat_s:i_lat_e,i_lon_s:i_lon_e], val[i_lat_s:i_lat_e,i_lon_s:i_lon_e], **kwargs)
    
#%% 
if os.environ['USER'] == 'kaandorp': # desktop
    dir_data = '/Volumes/externe_SSD/kaandorp/Data/'
elif os.environ['USER'] == 'kaand004': #gemini
    dir_data = '/data/oceanparcels/input_data'



#%% Create release points on the h3 mesh, resolution e.g. 4 (22km -> 200k surface particles) res 3 (60km) 

res_parent = 2
res_release = 4

test_release = False
if test_release:
    extent = (-15,-5,35,55)
    extent_release = ()

geoJson1 = {'type': 'Polygon', 'coordinates': [[[90,-180],[90,0],[-90,0],[-90,-180]]]}
geoJson2 = {'type': 'Polygon', 'coordinates': [[[90,0],[90,180],[-90,180],[-90,0]]]}

hexagons = list(h3.polyfill(geoJson1, res_parent)) + list(h3.polyfill(geoJson2, res_parent))

# fig = plt.figure(figsize=(10,8),dpi=120)  
# ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
# ax.coastlines(resolution='110m')

# for hex_ in hexagons:
#     x = np.array(h3.h3_to_geo_boundary(hex_))[:,1]
#     y = np.array(h3.h3_to_geo_boundary(hex_))[:,0]
    
    # if (np.all(x>0) or np.all(x<0) or np.all(np.abs(x)<100)) and (np.all(y>0) or np.all(y<0) or np.all(np.abs(y)<100)):
    #     ax.plot(x,y,transform=ccrs.PlateCarree())

# ax.set_extent((-180,180,-90,90))


hexagons_release = list(h3.polyfill(geoJson1, res_release)) + list(h3.polyfill(geoJson2, res_release))
hexagons_release_centre = np.array([h3.h3_to_geo(hex_) for hex_ in hexagons_release])

release_lon = hexagons_release_centre[:,1]
release_lat = hexagons_release_centre[:,0]

if test_release:
    mask = (release_lon > extent[0]) & (release_lon < extent[1]) & (release_lat > extent[2]) & (release_lat < extent[3])
    release_lon = release_lon[mask]
    release_lat = release_lat[mask]
# ax.plot(release_lon,release_lat,'ko',markersize=3)


#%% Select the particles in the water by looking at the land mask

# a mask with True values on land cells, and False on ocean cells
data_mask_land = xr.open_dataset('00_data_files/mask_land_NEMO0083.nc')
mask_land = data_mask_land['mask_land'].values
mask_land_lon = data_mask_land['lon'].values
mask_land_lat = data_mask_land['lat'].values

# Interpolate the release points onto this True/False mask using nearest neighbor. 
# This mask can now be used to filter out points on land (by using ~land_val_release)
land_val_release = griddata((mask_land_lon.ravel(),mask_land_lat.ravel()), mask_land.ravel(), (release_lon,release_lat), method='nearest')


# fig = plt.figure(figsize=(10,8),dpi=120)  
# ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
# ax.coastlines(resolution='110m')
# ax.plot(release_lon[~land_val_release], release_lat[~land_val_release], 'ko',markersize=1)
# ax.set_extent((-180,180,-90,90))


extent = (-30,20,30,60)

fig = plt.figure(figsize=(10,8),dpi=120)  
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
pcolormesh_NEMO(extent,mask_land_lon,mask_land_lat,mask_land,transform=ccrs.PlateCarree(),alpha=.4)
ax.plot(release_lon[~land_val_release], release_lat[~land_val_release], 'ro',markersize=1,transform=ccrs.PlateCarree())
ax.set_extent(extent)


#%% Write to netcdf release particlefile
# def calculate_particle_z(z_layers,n_per_layer):
#     p_z = np.array([])
#     for z_upper,z_lower in zip(z_layers[:-1],z_layers[1:]):
#         asdads
release_type = '2D'

if release_type == '3D':
    data_bathymetry = xr.open_dataset(os.path.join(dir_data, 'NEMO-MEDUSA/ORCA0083-N006/domain/bathymetry_ORCA12_V3.3.nc') )
    bathymetry = data_bathymetry['Bathymetry']
    
    # data_u_example = xr.open_dataset(os.path.join(dir_data, 'NEMO-MEDUSA/ORCA0083-N006/means/ORCA0083-N06_20100130d05U.nc') )
    # z_max = np.inf#500 # maximum depth to be simulated
    # mask = (data_u_example['depthu'].values < z_max)
    # z_levels_release = data_u_example['depthu'][mask][::3].values
    # z_levels_MOi = np.array([4.940254e-01, 1.541375e+00, 2.645669e+00, 3.819495e+00, 5.078224e+00,
    #        6.440614e+00, 7.929560e+00, 9.572997e+00, 1.140500e+01, 1.346714e+01,
    #        1.581007e+01, 1.849556e+01, 2.159882e+01, 2.521141e+01, 2.944473e+01,
    #        3.443415e+01, 4.034405e+01, 4.737369e+01, 5.576429e+01, 6.580727e+01,
    #        7.785385e+01, 9.232607e+01, 1.097293e+02, 1.306660e+02, 1.558507e+02,
    #        1.861256e+02, 2.224752e+02, 2.660403e+02, 3.181274e+02, 3.802130e+02,
    #        4.539377e+02, 5.410889e+02, 6.435668e+02, 7.633331e+02, 9.023393e+02,
    #        1.062440e+03, 1.245291e+03, 1.452251e+03, 1.684284e+03, 1.941893e+03,
    #        2.225078e+03, 2.533336e+03, 2.865703e+03, 3.220820e+03, 3.597032e+03,
    #        3.992484e+03, 4.405224e+03, 4.833291e+03, 5.274784e+03, 5.727917e+03])
    # z_levels_release = z_levels_MOi[::2]
    
    z_start = 0.5
    z_end = 5000
    n_steps = 4
    d_log = (np.log10(z_end) - np.log10(z_start))
    n_per_level = 3
    
    d_log_per_level = d_log / (n_per_level*n_steps)
    
    log_start = np.log10(z_start) + .5*d_log_per_level
    log_end = np.log10(z_end) - .5*d_log_per_level
    
    z_levels_release = np.logspace(log_start,log_end,n_steps*n_per_level)
    
    plt.figure()
    plt.semilogx(np.logspace(np.log10(z_start),np.log10(z_end),n_steps+1),np.zeros(5),'o')
    plt.semilogx(z_levels_release,np.zeros(len(z_levels_release)),'ro')

elif release_type == '2D':
    z_levels_release = np.array([0.5])

date_release = np.datetime64('2010-01-05T12:00:00')

n_particles = len(release_lon[~land_val_release])


p_lon = release_lon[~land_val_release].reshape([n_particles,1])
p_lat = release_lat[~land_val_release].reshape([n_particles,1])

bath_interpolated = griddata((data_bathymetry['nav_lon'].values.ravel(),data_bathymetry['nav_lat'].values.ravel()),
                             data_bathymetry['Bathymetry'].values.ravel(),(p_lon,p_lat),method='nearest')

p_lon_tiled = np.tile(p_lon,len(z_levels_release))
p_lat_tiled = np.tile(p_lat,len(z_levels_release))
p_z_tiled = np.ones(p_lon_tiled.shape)*z_levels_release

mask_keep = np.less(p_z_tiled,bath_interpolated.reshape(len(bath_interpolated),1))

n_particles_2 = int(mask_keep.sum())
p_lon_keep = p_lon_tiled[mask_keep].reshape([n_particles_2,1])
p_lat_keep = p_lat_tiled[mask_keep].reshape([n_particles_2,1])
p_z_keep = p_z_tiled[mask_keep].reshape([n_particles_2,1])

# p_lon = release_lon[~land_val_release].reshape([1,n_particles]).T
# p_lat = release_lat[~land_val_release].reshape([1,n_particles]).T
# p_z = np.zeros(len(release_lon)).reshape([1,n_particles]).T
# p_time = np.array([date_release for i in range(len(release_lon))]).reshape([1,len(release_lon)]).T
# trajectory = np.arange(len(release_lon)).reshape([1,len(release_lon)]).T
# to_netcdf('00_release_files/h3_uniform_%i_%.10s.nc' % (len(release_lon),date_release),[trajectory,p_time,p_lat,p_lon,p_z],['trajectory','time','lat','lon','z'])


# p_z = np.zeros(n_particles).reshape([n_particles,1])



p_time = np.array([date_release for i in range(n_particles_2)]).reshape([n_particles_2,1])
trajectory = np.arange(n_particles_2).reshape([n_particles_2,1])
to_netcdf('00_release_files/h%i_uniform_%s_%i.nc' % (res_release,release_type,n_particles_2),[trajectory,p_time,p_lat_keep,p_lon_keep,p_z_keep],['trajectory','time','lat','lon','z'])


#%% plot
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import cmocean
from tools import NEMO_select_section



extent = (-20,00,40,50)
cmap = cmocean.cm.deep_r
newcmap = cmocean.tools.crop_by_percent(cmap, 5, which='max')
newcmap.set_over(color='darkgoldenrod')
# fig = plt.figure(figsize=(10,8),dpi=120)  
# ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
# ax.coastlines()
# ax.set_extent(extent)

mask_plot = (p_lon_keep > extent[0]) & (p_lon_keep < extent[1]) & (p_lat_keep > extent[2]) & (p_lat_keep < extent[3])
lons_plot = p_lon_keep[mask_plot]
lats_plot = p_lat_keep[mask_plot]
z_plot = p_z_keep[mask_plot]

X,Y,Z = NEMO_select_section(extent,data_bathymetry['nav_lon'],data_bathymetry['nav_lat'],bathymetry)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X,Y,-Z,cmap=newcmap,vmax=-0.1,zorder=0)
ax.scatter(lons_plot,lats_plot,-z_plot,s=10)
