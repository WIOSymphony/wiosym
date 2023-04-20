#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:08:11 2021
v5: Implemented adjusted surface transport, permanent fouling
v6: remove suface transport for now, remove unnecessary parameters
@author: kaandorp
"""
import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from datetime import datetime, timedelta
from h3.unstable import vect
from h3 import h3_to_geo, geo_to_h3
from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize
import cartopy.crs as ccrs
from tools import from_pickle, to_pickle, pcolorhex
from utils import get_settling_velocity
# from tools import pcolorhex
from cartopy import feature
import matplotlib.gridspec as gridspec
from scipy import special
from time import time
import cmocean
from scipy import sparse
import gc
import socket
import multiprocessing
from multiprocessing import set_start_method
from scipy.interpolate import interp1d
from argparse import ArgumentParser
from tools import get_colors_tableau


def Gaussian(x,mu,sigma):
    return np.exp(-.5*((x-mu)/sigma)**2)

def draw_coastlines(ax,h3_indices,dict_coastlines,colors='r',**kwargs):
    for i1,index_ in enumerate(h3_indices):
        list_lon = dict_coastlines[index_]['lon_coast']
        list_lat = dict_coastlines[index_]['lat_coast']
        
        for lon_,lat_ in zip(list_lon,list_lat):
            if isinstance(colors,np.ndarray):
                ax.plot(lon_,lat_,'-',color=colors[i1],**kwargs)
            else:
                ax.plot(lon_,lat_,'-',color=colors,**kwargs)
                
def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))

    
def NB_model(k_arr,i_f,p,d_N=3):
    '''
    The fragmentation model by Charalambous (2015)
    k_arr: the size classes
    i_f: the fragmentation 'index' (denoted by f in the paper)
    p: fragmentation probability
    d_N: spatial dimensionality
    '''
    
    pmf_m = (special.gamma(k_arr+i_f) / (special.gamma(k_arr+1)*special.gamma(i_f)))*(p**k_arr)*(1-p)**i_f
    pmf_N = 2**(d_N*k_arr) * pmf_m
    
    return pmf_m,pmf_N


def geo_to_h3_depth(lat,lon,depth,res):
    h3_index = vect.geo_to_h3(lat,lon,res)
    
    if len(depth_layers) > 2:
        
        for i1,(d_upper,d_lower,index_remove) in enumerate(zip(depth_layers[:-1],depth_layers[1:],index_remove_depth)):
            mask = (depth >= d_upper) & (depth < d_lower)
            if i1 == 0: #surface: don't change the h3 index
                pass
            elif i1 == len(depth_layers) - 2: #below the last depth layer: move particles to index 1 (removal due to sinking too deep)
                h3_index[mask] = 1
                print('%i particles removed below depth %f meters' % (mask.sum(),d_upper))
            else:
                h3_index[mask] -= index_remove #at other depth layers: remove the appropriate depth index
            print('n particles in depth layer: %i' % mask.sum())

    return h3_index


def get_fragmentation_fraction(lambda_f_week = 2e-4,units='m'):
    l_0 = 160
    k_arr = np.arange(10)
    l_arr = l_0 / (2**k_arr)

    pmf_m_0,pmf_N_0 = NB_model(k_arr,1,.4,d_N=2.5)

    lambda_f = lambda_f_week*4 #f/week to f/month
    
    m_next, N_next = NB_model(k_arr, lambda_f, .4, d_N=2.5)

    macro_to_macro = 0
    macro_to_micro = 0
    for i1, (m_, l_) in enumerate(zip(pmf_m_0, l_arr)):
        if i1 == 0:
            i_start = None
            i_end = None
        else:
            i_start = i1
            i_end = -i1
            
        m_arr_ = m_next[:i_end]
        l_arr_ = l_arr[i_start:]    
        
        mask_macro = l_arr_ > 5
        
        macro_to_macro += m_*(m_arr_[mask_macro].sum())
        macro_to_micro += m_*(m_arr_[~mask_macro].sum())
        
    return macro_to_macro,macro_to_micro


def h3_cell_volume():
    
    h3_index_upper = h3_index_base + np.uint64(1e17)
    
    dict_h3_normalize = pd.Series(index=map_h3_to_mat.index,dtype=np.float32)

    for i1 in range(len(depth_layers)-1):
        
        mask_select = (map_h3_to_mat.index < (h3_index_upper - index_remove_depth[i1])) & (map_h3_to_mat.index >= (h3_index_upper - index_remove_depth[i1+1]))
        
        h3_orig = map_h3_to_mat[mask_select].index + index_remove_depth[i1]
        
        areas = np.array([dict_coastlines[h3_orig[i1]]['area'] for i1 in range(len(h3_orig))]) *1e6 #m2
        
        if i1 == len(depth_layers)-2: #in the lowest layer, calculate difference between local bathymetry and upper depth level (otherwise this results in np.inf)
            bath_ = df_bath[h3_orig].values
            height = bath_ - depth_layers[i1]
        else:
            height = depth_layers[i1+1] - depth_layers[i1]

        # dict_h3_area[h3_orig] = areas
        # dict_h3_volume[h3_orig] = areas*height
        
        dict_h3_normalize[map_h3_to_mat[mask_select].index] = areas*height
     
    #coastlines
    mask_select = map_h3_to_mat.index > h3_index_upper
    h3_orig = map_h3_to_mat[mask_select].index - np.uint64(1e17)
    lengths = np.array([dict_coastlines[h3_orig[i1]]['length_coast_D1.27'] for i1 in range(len(h3_orig))]) *1e3
    dict_h3_normalize[map_h3_to_mat[mask_select].index] = lengths

    return dict_h3_normalize

def h3_cell_location():
 
    h3_index_upper = h3_index_base + np.uint64(1e17)
    
    df_h3_location = pd.DataFrame(index=map_h3_to_mat.index,columns=('lon','lat'),dtype=np.float32)

    for i1 in range(len(depth_layers)-1):
        
        mask_select = (map_h3_to_mat.index < (h3_index_upper - index_remove_depth[i1])) & (map_h3_to_mat.index >= (h3_index_upper - index_remove_depth[i1+1]))
        
        h3_orig = map_h3_to_mat[mask_select].index + index_remove_depth[i1]
        
        lons_ = []
        lats_ = []
        for i2 in range(len(h3_orig)):
            (lat,lon) = h3_to_geo(hex(h3_orig[i2]))
            lats_.append(lat)
            lons_.append(lon)
            
        
        df_h3_location.loc[map_h3_to_mat[mask_select].index,'lon'] = np.array(lons_)
        df_h3_location.loc[map_h3_to_mat[mask_select].index,'lat'] = np.array(lats_)
        
     
    #coastlines
    mask_select = map_h3_to_mat.index > h3_index_upper
    h3_orig = map_h3_to_mat[mask_select].index - np.uint64(1e17)

    lons_ = []
    lats_ = []
    for i2 in range(len(h3_orig)):
        (lat,lon) = h3_to_geo(hex(h3_orig[i2]))
        lats_.append(lat)
        lons_.append(lon)
    df_h3_location.loc[map_h3_to_mat[mask_select].index,'lon'] = np.array(lons_)
    df_h3_location.loc[map_h3_to_mat[mask_select].index,'lat'] = np.array(lats_)
    
    return df_h3_location      
    

# def plot_dist(ax,vmax=1.,vmin=1e-16,cmap=plt.cm.tab10):
    
#     for i1,index_ in enumerate(mat_index_plot_dist.values):
        
#         dist = tracer_matrix_mass[index_,:] 
#         # n_points = (dist > 0).sum()
    
#         dist_valid = dist[dist>0]
#         # k_valid = k_arr[dist>0]
#         # l_arr = l0 / (2**k_arr)
#         l_valid = l_arr[dist>0]
       
        
#         if vmax==None:
#             vmax = dist.max()*2
#         if vmin==None:
#             vmin = max(1e-12,dist.min())    
    
#         for l_, n_ in zip(l_valid,dist_valid):
#             ax.loglog(l_,n_,'o',color=cmap(i1) )
#             ax.loglog([l_,l_],[vmin,n_],'-',color=cmap(i1) )

#     ax.set_xlim(0.8*l_arr.min(),1.2*l_arr.max())
#     ax.set_ylim(vmin,vmax)


def create_release_vector(type_='r',mu_log10_release=-.7,sigma_log10_release=.15):
   
    weights = Gaussian(np.log10(l_arr),mu_log10_release,sigma_log10_release)
    weights = weights / weights.sum()
    
    if type_ == 'r':
        file_release = os.path.join(data_folder,'00_data_files/input_h3_rivers.nc')
        data_release = xr.open_dataset(file_release)

        i_valid_h3 = np.in1d(data_release['h3_index'].values, map_mat_to_h3.values) #release locations present in transition matrix
        h3_release = data_release['h3_index'].values[i_valid_h3]
        tracer_release = (data_release['midpoint'].values[0,i_valid_h3] / 12) * 1e6 #tonnes/year to grams/month

        release_per_month = {}
        for i1 in range(12):
            release_per_month[i1] = np.zeros([n_total,len(l_arr)])

            for i2 in range(len(l_arr)):
                np.add.at(release_per_month[i1][:,i2], map_h3_to_mat[h3_release].values, tracer_release*weights[i2])    

    elif type_ == 'f':
        file_release = os.path.join(data_folder,'00_data_files/input_h3_fisheries_monthly.nc')
        data_release = xr.open_dataset(file_release)

        i_valid_h3 = np.in1d(data_release['h3_index'].values, map_mat_to_h3.values) #release locations present in transition matrix
        h3_release = data_release['h3_index'].values[i_valid_h3]

        release_per_month = {}
        total_release = 0
        # release_per_month = np.zeros([12,n_total])
        for i1 in range(12):
            release_per_month[i1] = np.zeros([n_total,len(l_arr)])
            tracer_release = data_release['fishing_hours'].values[i1,i_valid_h3]     

            for i2 in range(len(l_arr)):
                np.add.at(release_per_month[i1][:,i2], map_h3_to_mat[h3_release].values, tracer_release*weights[i2])

            total_release += release_per_month[i1].sum()
            # np.add.at(release_per_month[i1,:], map_h3_to_mat[h3_release].values, tracer_release)

        for i1 in range(12):
            #normalize to have a mean monthly input of 1
            release_per_month[i1] /= (total_release / 12)

    elif type_ == 'p':
        file_release = os.path.join(data_folder,'00_data_files/input_h3_coastal_pop.nc')
        data_release = xr.open_dataset(file_release)

        i_valid_h3 = np.in1d(data_release['h3_index'].values, map_mat_to_h3.values) #release locations present in transition matrix
        h3_release = data_release['h3_index'].values[i_valid_h3]
        tracer_release = data_release['MPW_input'].values[3,i_valid_h3] #choose 2015?

        release_per_month = {}
        total_release = 0
        for i1 in range(12):
            release_per_month[i1] = np.zeros([n_total,len(l_arr)])

            for i2 in range(len(l_arr)):
                np.add.at(release_per_month[i1][:,i2], map_h3_to_mat[h3_release].values, tracer_release*weights[i2])    
            total_release += release_per_month[i1].sum()

        for i1 in range(12):
            #normalize to have a mean monthly input of 1
            release_per_month[i1] /= (total_release / 12)
   
    return release_per_month


def calculate_rise_velocity_Dietrich(l,rho_p=1010,rho_sw=1025):
    w_b_arr = []
    # d_star_arr = []
    for r_tot in l:
        # rho_sw = 1029
        g = 9.81
        sw_kin_visc = 1e-6
        
        dn = r_tot  # equivalent spherical diameter [m], calculated from Dietrich (1982) from A = pi/4 * dn**2
        delta_rho = (rho_p - rho_sw) / rho_sw  # normalised difference in density between total plastic+bf and seawater[-]
        dstar = (abs(rho_p - rho_sw) * g * dn ** 3.) / (rho_sw * sw_kin_visc ** 2.)  # [-]
    
        if dstar > 5e9:
            w_star = 265000
        elif dstar < 0.05:
            w_star = (dstar ** 2.) * 1.71E-4
        else:
            w_star = 10. ** (-3.76715 + (1.92944 * np.log10(dstar)) - (0.09815 * np.log10(dstar) ** 2.) - (
                        0.00575 * np.log10(dstar) ** 3.) + (0.00056 * np.log10(dstar) ** 4.))
        
        if delta_rho > 0:  # sinks
            w_b = -(g * sw_kin_visc * w_star * delta_rho) ** (1. / 3.)
        else:  # rises
            a_del_rho = delta_rho * -1.
            w_b = (g * sw_kin_visc * w_star * a_del_rho) ** (1. / 3.)  # m s-1    
    
        w_b_arr.append(w_b)
        # d_star_arr.append(dstar)
        
    return np.array(w_b_arr)

def print_time(t1,str_,div_by=1.):
    units = ''
    if div_by == 1:
        units = 'seconds'
    elif div_by == 60:
        units = 'minutes'
    print('time %s: %f %s' % (str_,((time()-t1)/div_by),units))
    return time()
#%%
gc.enable()
# def plot_h3_cells(h3_cells):
#     pcolorhex(ax,hexagons,colors

def get_available_data(filename,map_h3_to_mat,plot_missing=False,index_add=0):
    data = pd.read_csv(filename,index_col=0,parse_dates=['eventDate'])
    date_cols = [col for col in data.columns if 'Date' in col]
    for col in date_cols:
        data[col] = pd.to_datetime(data[col])
    
    mask_present_in_mat = ( np.isin(data['h%i'%h3_grid_res]+index_add,map_h3_to_mat.index) & ~(data['h%i'%h3_grid_res] == 0) )
    data_return = data.loc[mask_present_in_mat,:].copy()
    print('%f percent of data is located, n = %i' % ((mask_present_in_mat.sum()/len(mask_present_in_mat))*100,mask_present_in_mat.sum()) )
    if plot_missing:
        from cartopy.feature import LAND
        fig = plt.figure(figsize=(13, 6))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.add_feature(LAND, zorder=0,edgecolor='black')
        h3_plot = np.unique(data.loc[~mask_present_in_mat,'h%i'%h3_grid_res])
        colors = [plt.cm.tab10(1)]*len(h3_plot)
        pcolorhex(ax,h3_plot,colors,alpha=.7)
        ax.plot(data.loc[~mask_present_in_mat,'decimalLongitude'],data.loc[~mask_present_in_mat,'decimalLatitude'],'ko')
    
    if 'l_min' in data_return.keys():
        data_return.loc[np.isnan(data_return.loc[:,'l_min']),'l_min'] = 0
    if 'l_max' in data_return.keys():
        data_return.loc[np.isnan(data_return.loc[:,'l_max']),'l_max'] = np.inf
    return data_return, mask_present_in_mat

if 'mbp13' in socket.gethostname(): # desktop
    home_folder = '.'
    data_folder = '.'
    n_sizes = 2
elif 'node0' in socket.gethostname() or 'lorenz' in socket.gethostname(): #Lorenz
    home_folder = '/nethome/kaand004/Git_repositories/Global_Analysis_Mikael'    
    data_folder = '/nethome/kaand004/Git_repositories/Global_Analysis_Mikael' 
    n_sizes = 15
elif 'science-bs' in socket.gethostname(): #gemini
    home_folder = '/nethome/kaand004/Git_repositories/Global_Analysis_Mikael'    
    data_folder = '/scratch/kaand004' 
    n_sizes = 15
    
folder_figure = os.path.join(home_folder,'00_figures/17_inversion_01/')
if not os.path.exists(folder_figure):
    os.mkdir(folder_figure)

#3: ~40,000 cells, 60km
#2: ~6,000 cells, 160km
h3_grid_res = 3
coastline_res = '50m'
use_Hinata = True
use_Stokes = False
depth_layers = np.array([0,5,50,500,np.inf])
particles_l = [0.0001,0.0004,0.0016,0.0064,0.0256,0.1024]
months = [1,2,3,4,5,6,7,8,9,10,11]
years = [2019]
year_start = 1980
dt_days = 30
l0 = 1.6384 
# l0 = 0.2048

release = 'mix' #simple, 'r','f','p', or 'mix'
mass_count_ratio = 500. #grams/unit, 0.4096 meter objects
k_arr = np.arange(0,n_sizes) #11: 0.4mm to 0.4m
l_arr = l0 / (2**k_arr)

str_info = 'Year_' + '-'.join('%2.2i'%s_ for s_ in years) + '_Stokes%s_'%use_Stokes

#-----------------parameters------------------------
h3_index_base = np.uint64(5e17)
index_add_coast = np.uint64(1e17)     #indices for coastal cells are the ocean surface cells + index_add_coast
index_remove_depth = np.array([int(0),int(1e17),int(2e17),int(3e17),int(4e17)])

# the dict in which all coastline data is stored
dict_coastlines = from_pickle(os.path.join(data_folder,'00_data_files/coastal_properties_h3_res%i_%s_Dc.pickle' % (h3_grid_res,coastline_res)))
h3_index_coastline_data = np.array([hex_ for hex_ in dict_coastlines.keys() if 'length_coast' in dict_coastlines[hex_].keys() ],dtype=np.uint64)

#TODO: monthly oceanic transitions below
filename_dict_transition = os.path.join(data_folder,'00_transition_files/transitions_h%i_%s_monthly_' % (h3_grid_res,coastline_res) + 'year_' + '-'.join('%2.2i'%s_ for s_ in years) + '_merge_months_False.pickle')


P = from_pickle(filename_dict_transition)
map_mat_to_h3 = P['map_mat_to_h3']
map_h3_to_mat = P['map_h3_to_mat']

# a map from the matrix entry to removal (-1), surface (0), and deeper layers (1,2,3,4) to make computations easier
# TODO: can be made variable for different depth layers
# map_mat_to_domain = pd.Series(index=map_mat_to_h3.index,dtype=int)
# map_mat_to_domain[map_mat_to_h3.values < 1e17] = -1
# map_mat_to_domain[map_mat_to_h3.values > 6e17] = 0
# map_mat_to_domain[(map_mat_to_h3.values < 6e17) & (map_mat_to_h3.values > 5e17)] = 1
# map_mat_to_domain[(map_mat_to_h3.values < 5e17) & (map_mat_to_h3.values > 4e17)] = 2
# map_mat_to_domain[(map_mat_to_h3.values < 4e17) & (map_mat_to_h3.values > 3e17)] = 3
# map_mat_to_domain[(map_mat_to_h3.values < 3e17) & (map_mat_to_h3.values > 2e17)] = 4

# not all coastline sections are present in the transition matrix: take the intersection
h3_index_coast = np.intersect1d(h3_index_coastline_data, map_mat_to_h3)
length_coast = np.array([dict_coastlines[hex_]['length_coast_D1.27'] for hex_ in h3_index_coast])
map_h3_to_coast_length = pd.Series(index=h3_index_coast, data=length_coast )

dict_variance = {'n_surf':0.28,'m_surf':0.40,'n_beach':0.2,'m_beach':0.14,'fish_frac':0.008,
                 'MDMAP':0.21,'OSPAR_100m':0.29,'OSPAR_1km':0.32,
                 'surf_nm-3_Egger':0.08,'surf_gm-3_Egger':0.09,
                 'depth_nm-3_Egger':0.04,'depth_gm-3_Egger':0.04,}

data_trawl_m, _ = get_available_data(os.path.join(data_folder,'00_data_files/plastic_data/converted_surface_m.csv'),map_h3_to_mat)

data_trawl_n, _ = get_available_data(os.path.join(data_folder,'00_data_files/plastic_data/converted_surface_n.csv'),map_h3_to_mat)

data_macro_n, _ = get_available_data(os.path.join(data_folder,'00_data_files/plastic_data/converted_surface_macro.csv'),map_h3_to_mat)

data_beach_inst, _ = get_available_data(os.path.join(data_folder,'00_data_files/plastic_data/converted_beach_inst.csv'),map_h3_to_mat,index_add=index_add_coast)

data_beach_OSPAR, _ = get_available_data(os.path.join(data_folder,'00_data_files/plastic_data/converted_agg_MDMAP_OSPAR.csv'),map_h3_to_mat,index_add=index_add_coast)

data_beach_mean, _ = get_available_data(os.path.join(data_folder,'00_data_files/plastic_data/converted_beach_mean.csv'),map_h3_to_mat,index_add=index_add_coast)
data_beach_mean['modelled_n'] = 0
data_beach_mean['modelled_sum'] = 0
data_beach_mean['modelled_sum2'] = 0

# mask_present_in_mat = np.isin(data_vanSebille2015['h%i'%h3_grid_res],map_h3_to_mat.index)
# data_vanSebille2015 = data_vanSebille2015.loc[mask_present_in_mat,:] 

# data_vanSebille2015 = pd.read_csv(os.path.join(data_folder,'00_data_files/plastic_data/converted_vanSebille2015.csv'),index_col=0)
# # remove measurements for which the h3 cell is not present in transition matrix
# # TODO: find closest existing cell? 
# mask_present_in_mat = np.isin(data_vanSebille2015['h%i'%h3_grid_res],map_h3_to_mat.index)
# data_vanSebille2015 = data_vanSebille2015.loc[mask_present_in_mat,:] 

# data_AtlantECO = pd.read_csv(os.path.join(data_folder,'00_data_files/plastic_data/converted_AtlantECO.csv'),index_col=0)
# mask_present_in_mat = np.isin(data_AtlantECO['h%i'%h3_grid_res],map_h3_to_mat.index)
# data_AtlantECO = data_AtlantECO.loc[mask_present_in_mat,:] 

# data_BBCT = pd.read_csv(os.path.join(data_folder,'00_data_files/plastic_data/converted_beach_BBCT.csv'),index_col=0)
# mask_present_in_mat = np.isin(data_BBCT['h%i'%h3_grid_res]+index_add_coast,map_h3_to_mat.index)
# data_BBCT = data_BBCT.loc[mask_present_in_mat,:] 

# data_OSPAR = pd.read_csv(os.path.join(data_folder,'00_data_files/plastic_data/converted_beach_ospar.csv'),index_col=0,parse_dates=['eventDate','eventDate_end'])
# mask_present_in_mat = np.isin(data_OSPAR['h%i'%h3_grid_res]+index_add_coast,map_h3_to_mat.index)
# data_OSPAR = data_OSPAR.loc[mask_present_in_mat,:] 
# data_OSPAR['modelled_n'] = 0
# data_OSPAR['modelled_sum'] = 0
# data_OSPAR['modelled_sum2'] = 0

data_Egger_surf = pd.read_csv(os.path.join(data_folder,'00_data_files/plastic_data/converted_Egger2020_ERL_surface.csv'),index_col=0)
data_Egger_surf = pd.concat([data_Egger_surf,
                             pd.read_csv(os.path.join(data_folder,'00_data_files/plastic_data/converted_Egger2020_SR_surface.csv'),index_col=0)])
mask_present_in_mat = np.isin(data_Egger_surf['h%i'%h3_grid_res],map_h3_to_mat.index)
data_Egger_surf = data_Egger_surf.loc[mask_present_in_mat,:] 

data_Egger_depth = pd.read_csv(os.path.join(data_folder,'00_data_files/plastic_data/converted_Egger2020SR_depth_frag.csv'),index_col=0)
data_Egger_depth['h%i'%h3_grid_res] = data_Egger_depth['h%i'%h3_grid_res] - index_remove_depth[np.digitize(data_Egger_depth['Depth'],depth_layers)-1]
mask_present_in_mat = np.isin(data_Egger_depth['h%i'%h3_grid_res],map_h3_to_mat.index)
data_Egger_depth = data_Egger_depth.loc[mask_present_in_mat,:] 

data_Zhao_depth = pd.read_csv(os.path.join(data_folder,'00_data_files/plastic_data/converted_Zhao2022_depth.csv'),index_col=0)
data_Zhao_depth['h%i'%h3_grid_res] = data_Zhao_depth['h%i'%h3_grid_res] - index_remove_depth[np.digitize(data_Zhao_depth['Depth'],depth_layers)-1]
mask_present_in_mat = np.isin(data_Zhao_depth['h%i'%h3_grid_res],map_h3_to_mat.index)
data_Zhao_depth = data_Zhao_depth.loc[mask_present_in_mat,:] 


df_bath = from_pickle(os.path.join(data_folder,'00_data_files/bathymetry_h%i.pickle' % h3_grid_res))
df_pp = from_pickle(os.path.join(data_folder,'00_data_files/mean_pp_h%i_2019.pickle' % h3_grid_res)).astype(float)
df_pp = df_pp / np.quantile(df_pp.values[~np.isnan(np.array(df_pp.values,dtype=float))],.99) #normalize for outliers
df_pp[df_pp>1] = 1.
df_pp[df_pp<0] = 0.

dict_h3_normalize = h3_cell_volume()
dict_h3_location = h3_cell_location()

n_total = len(map_h3_to_mat)

lon_plot_dist = np.array([-140.,-112.,-140.,-100.])
lat_plot_dist = np.array([20.,23.,50.,15.])
depth_plot_dist = np.array([0,0,0,0])


mask_neuston_net = (l_arr > 3e-4) & (l_arr < 0.2)
mask_beach_cleanup = l_arr > 2.5e-2
# mask_Egger = [(l_arr>=4e-4)&(l_arr<1.5e-3),(l_arr>1.5e-3)&(l_arr<5.0e-3),
#               (l_arr>5e-3)&(l_arr<1.5e-2),(l_arr>1.5e-2)&(l_arr<5.2e-2)]
mask_Egger = [(l_arr>=5.6e-4)&(l_arr<1.13e-3),(l_arr>1.13e-3)&(l_arr<4.5e-3),
              (l_arr>4.5e-3)&(l_arr<18e-3),(l_arr>18e-3)&(l_arr<36e-3)]

mask_Egger_sizes = [[5e-4,1.5e-3],[1.5e-3,5e-3],[5e-3,1.5e-2],[1.5e-2,5e-2]]
Egger_corr_n = [1.5029822975240785,
 0.6564421139154354,
 0.7997358204039169,
 1.7005002991087623]
Egger_corr_m = [1.6803211760254935,
 0.9670585583407109,
 0.7636187765408871,
 1.8350271441960349]

mask_Zhao = (l_arr >= 1e-4) & (l_arr < 5e-3)

    
parameter_names = [r'$S_{riv}$',r'$S_{pop}$',r'$S_{fis}$',
                   r'$\tau_{beach}$',r'$l_{coast,beach,max}$',r'$a_{resus.}$',r'$b_{resus.}$',
                   r'$f_{fouled}$',r'$f_{nb}$',r'$f_{pf}$',
                   r'$p_{land,remove}$',
                   r'$\lambda_{f}$',r'$p_f$',r'$d_N$',
                   r'$a_{L_{eff}}$',r'$b_{L_{eff}}$',
                  r'$\mu_{release,log10}$','PW_factor']

log_parameters = np.array([True,True,True,
                           True,False,False,False,
                           False,False,True,
                           True,
                           True,False,False,
                           False,False,
                          False,False])

# mu_parameters = np.array([-1.0,0.7,0.0,
#                           np.log10(100),125,260,7.1,
#                           .2,.2,.2,
#                           np.log10(0.01),
#                           -3.0,0.4,2.18,
#                           .64,-1.,
#                          -.4,0.05])

# sigma_parameters = np.array([0.16,0.08,0.22,
#                              .25,30,1.,0.1,
#                              .1,.1,.1,
#                              .3,
#                              0.16,.001,.1,
#                              0.06,0.15,
#                             .13,0.01]) #sigma is given for the log-transformed parameters 

# mu_parameters = np.array([-9.89405350e-01,  6.98970004e-01,  3.57155523e-02,  
#                           2.00000000e+00, 4.70465152e+02,  2.60000000e+02,  7.10000000e+00,  
#                           2.50000000e-01, 2.50000000e-01,  -8.92197672e-01, 
#                           -2.69929574e+00, 
#                           -2.66549661e+00, 4.00000000e-01,  2.18718327e+00,  
#                           6.37325580e-01, -9.48615760e-01,
#                           -3.92022000e-01,  5.00000000e-02])

# sigma_parameters = np.array([4.10814466e-01, 7.04270806e-02, 1.77824465e-01, 
#                              2.00686664e-01, 1.24409219e+02, 2.60000000e+00, 7.10000000e-02, 
#                              8.10000000e-02, 8.10000000e-02, 2.93321072e-01, 
#                              3.33333333e-01, 
#                              1.43804468e-01, 1.00000000e-03, 4.60050728e-02, 
#                              8.74484991e-02, 2.26095201e-01,
#                              1.50947667e-01, 1.66666667e-02])


#2 sigma
mu_parameters = np.array([-9.89405350e-01,  6.98970004e-01,  3.57155523e-02,  
                          2.00000000e+00,  4.70465152e+02,  2.60000000e+02,  7.10000000e+00,  
                          2.50000000e-01, 2.50000000e-01, -8.92197672e-01, 
                          -2.69929574e+00, 
                          -2.66549661e+00, 4.20000000e-01,  2.18718327e+00,  
                          6.37325580e-01, -9.48615760e-01,
       -3.92022000e-01,  2.50000000e-02])

sigma_parameters = np.array([4.10814466e-01, 7.04270806e-02, 1.77824465e-01, 2.00686664e-01,
       1.24409219e+02, 2.60000000e+00, 7.10000000e-02, 8.10000000e-02,
       8.10000000e-02, 2.93321072e-01, 3.33333333e-01, 1.43804468e-01,
       1.00000000e-02, 4.60050728e-02, 8.74484991e-02, 2.26095201e-01,
       1.50947667e-01, 8.33333333e-03])


cap_0_1 = np.array([False,False,False,
                    False,False,False,False,
                    True,True,True,
                    True,
                    False,False,False,
                    False,False,
                   False,False])

def forward_model(parameters,final_run=False):
                
    def draw_map(vec_,extent,gs,col_,base_title='k0, %3.3f-%3.3f m depth',stop_earlier=0,log_scale=True):
    
        h3_index_upper = h3_index_base + np.uint64(1e17)
    
    
        for i1 in range(len(depth_layers)-1-stop_earlier):
    
    
            ax_ = plt.subplot(gs[i1,col_], projection=ccrs.PlateCarree())
    
            mask_select = (map_h3_to_mat.index < (h3_index_upper - index_remove_depth[i1])) & (map_h3_to_mat.index >= (h3_index_upper - index_remove_depth[i1+1]))
    
            h3_orig = map_h3_to_mat[mask_select].index + index_remove_depth[i1]
            vec_select = vec_[mask_select]
    
            if not log_scale:
                colors = get_colors(vec_select,cmap,
                                    vmin=vmin,vmax=vmax)
                h3_plot = h3_orig
            else:
                colors = get_colors(np.log10(vec_select[vec_select>0]),cmap,
                                    vmin=np.log10(vmin),vmax=np.log10(vmax))
                h3_plot = h3_orig[vec_select > 0]
    
            mask = (dict_h3_location.loc[h3_plot,'lon']>extent[0]) & (dict_h3_location.loc[h3_plot,'lon']<extent[1]) & (dict_h3_location.loc[h3_plot,'lat']>extent[2]) & (dict_h3_location.loc[h3_plot,'lat']<extent[3])
            h3_plot = h3_plot[mask]
            colors = colors[mask,:]        
    
            pcolorhex(ax_,h3_plot,colors,draw_edges=False,alpha=1.)
    
            ax_.add_feature(feature.NaturalEarthFeature('physical','land',coastline_res),facecolor='grey',zorder=1000)
            if i1 == 0:
                mask_select = map_h3_to_mat.index > h3_index_upper # coastlines
                h3_orig = map_h3_to_mat[mask_select].index - np.uint64(1e17)
                vec_select = vec_[mask_select]
    
                if not log_scale:
                    colors = get_colors(vec_select,cmap2,
                                        vmin=vmin2,vmax=vmax2)
                    h3_plot = h3_orig
                else:
                    colors = get_colors(np.log10(vec_select[vec_select>0]),cmap2,
                                        vmin=np.log10(vmin2),vmax=np.log10(vmax2))
                    h3_plot = h3_orig[vec_select > 0]
    
                mask = (dict_h3_location.loc[h3_plot,'lon']>extent[0]) & (dict_h3_location.loc[h3_plot,'lon']<extent[1]) & (dict_h3_location.loc[h3_plot,'lat']>extent[2]) & (dict_h3_location.loc[h3_plot,'lat']<extent[3])
                h3_plot = h3_plot[mask]
                colors = colors[mask,:]  
    
                # colors = get_colors(np.log10(vec_select[vec_select>0]),cmap2,
                #                     vmin=np.log10(vmin2),vmax=np.log10(vmax2))
                # h3_plot = h3_orig[vec_select > 0]
                if vec_select.sum() > 0:
                    draw_coastlines(ax_,h3_plot,dict_coastlines,transform=ccrs.PlateCarree(),linewidth=7)    
                    draw_coastlines(ax_,h3_plot,dict_coastlines,colors=colors,transform=ccrs.PlateCarree(),linewidth=6)   
    
            ax_.set_extent((extent[0]+1,extent[1]-1,extent[2]+1,extent[3]-1),ccrs.PlateCarree()) #reduce the margins to hide white cells
            ax_.set_title(base_title % (depth_layers[i1],depth_layers[i1+1]))
    
            
    t1 = time()

    
    def to_scaled(p_unscaled,mu_params_,sigma_params_):
        p_scaled = (p_unscaled - mu_params_) / sigma_params_
        return p_scaled
    
    def to_unscaled(p_scaled,mu_params_,sigma_params_):
        p_unscaled = (p_scaled*sigma_params_) + mu_params_
        return p_unscaled
    
    def lin_max(length,l_max,max_val=1.):
        if length.size > 1:
            p = (1/l_max)*length
            p[length > l_max] = max_val
        elif length.size == 1:
            p = min((1/l_max)*length,max_val)
        else:
            p = None
        return p
    #-------------get back the parameters
    parameters_unscaled = to_unscaled(parameters,mu_parameters,sigma_parameters)
    
    parameters_unscaled[log_parameters] =10**(parameters_unscaled[log_parameters])
    parameters_unscaled[cap_0_1 & (parameters_unscaled < 0)] = 0.01
    parameters_unscaled[cap_0_1 & (parameters_unscaled > 1)] = .99
    
    for name_,p_ in zip(parameter_names,parameters_unscaled):
        print('%s: %4.4e' % (name_,p_))
    
    (r_riv,r_pop,r_fis,
     tau_b_,l_max,a_resus,b_resus,
     frac_fouled,frac_nb,frac_pf,
     p_land_remove,
     frag_if_month,p_f,d_N,
     a_L_Leff,b_L_Leff,
     mu_log10_release,PW_factor) = parameters_unscaled
    
    sigma_L_eff = .13
    pf_mode = 'permanent_fouled' # permanent_fouled_removal500
    
    # unrealisic scenario: per time step more than the available tracer is removed (shouldnt happen)
    # if 1-p_ocean_max_remove <= 0:
    #     print('max removal ocean outside of bounds (%f), setting to .99' % (p_ocean_max_remove))
    #     p_ocean_max_remove = .99
    if 1-p_land_remove <= 0:
        print('max removal land outside of bounds (%f), setting to .95' % (p_land_remove))
        p_land_remove = .95
    if b_resus < 0:
        b_resus = 0
    if l_max < 0:
        l_max = .001
    # fractions of fouled and neutrally bouyant particles shouldn't make up more than 100%
    if frac_fouled + frac_nb + frac_pf >= 1:
        print('adjusting fouled/nb/pf fractions to 98% in total')
        sum_ = frac_fouled + frac_nb + frac_pf
        frac_fouled = (frac_fouled / sum_) - 0.01
        frac_nb = (frac_nb / sum_) - 0.01
        frac_pf = (frac_pf / sum_) - 0.01
        
    #t1 = print_time(t1,'init')
    
    # p_ocean_min_remain = (1-p_ocean_max_remove)
    p_land_remain = (1-p_land_remove)
    
    release_vector = np.array([r_riv,r_pop,r_fis]) 
    
    def calc_mass_per_item(l_arr,d_N):
        l_calib = 0.0064 #from calculate_priors.py
        m_calib = 0.0161
        mass_per_item = np.zeros(len(l_arr))
        mask_small = l_arr <= l_calib
        mask_large = l_arr > l_calib
        
        mass_per_item[mask_small] =  ((l_arr[mask_small] / l_calib)**2.187)*m_calib
        mass_per_item[mask_large] = ((l_arr[mask_large] / l_calib)**d_N)*m_calib
        
        return mass_per_item
        
    # mass_per_item = ((l_arr / 0.0064)**d_N)*0.02039 #from test_count_mass_ratio.py
    mass_per_item = calc_mass_per_item(l_arr,d_N)
    
    tau_b = tau_b_ * np.ones(len(l_arr)) #days
    if not use_Hinata:
        tau_r = 25 * np.ones(len(l_arr)) #days
    else:
        w_b = calculate_rise_velocity_Dietrich(l_arr)
        tau_r = a_resus*w_b + b_resus
    
    p_frag_m,p_frag_N = NB_model(k_arr,frag_if_month,p_f,d_N=d_N)
    
    p_beach = 1-np.exp(-dt_days/tau_b)
    p_resus = 1-np.exp(-dt_days/tau_r)
        # sigma_log10_wb = 0.2
    
    if release == 'rivers' or release == 'r':
        vec_start = create_release_vector(type_='r')  
    elif release == 'coastlines' or release == 'p':
        vec_start_r = create_release_vector(type_='r')
        mass_calibrate = vec_start_r.sum() / 12
        vec_start = create_release_vector(type_='p')  #1/month
        vec_start *= release_vector[1]*mass_calibrate
    elif release == 'fisheries' or release == 'f':
        vec_start_r = create_release_vector(type_='r')
        mass_calibrate = vec_start_r.sum() / 12
        vec_start = create_release_vector(type_='f')  #1/month
        vec_start *= release_vector[2]*mass_calibrate
    
    elif release == 'mix': #mixed scenario
        vec_start_r = create_release_vector(type_='r',mu_log10_release=mu_log10_release,sigma_log10_release=.3) #grams/month
        vec_start_p = create_release_vector(type_='p',mu_log10_release=mu_log10_release,sigma_log10_release=.3)
        vec_start_f = create_release_vector(type_='f',mu_log10_release=mu_log10_release,sigma_log10_release=.3)
    
        mass_calibrate = vec_start_r[0].sum()  #grams/month
    
        vec_start = {}
        vec_start_fN = {}
        vec_start_N = {}
        for i2 in range(12):        
            vec_start_r[i2] *= release_vector[0]
            vec_start_p[i2] *= (release_vector[1]*release_vector[0]*mass_calibrate)
            vec_start_f[i2] *= (release_vector[2]*release_vector[0]*mass_calibrate)
    
            vec_start[i2] = vec_start_r[i2] + vec_start_p[i2] + vec_start_f[i2] #each column: in grams per month
            vec_start_fN[i2] = vec_start_f[i2] / mass_per_item
            vec_start_N[i2] = vec_start[i2] / mass_per_item
    
        print('Releasing %3.3f kilotons of plastics from rivers per year' % (np.array([vec_start_r[i].sum() for i in range(12)]).sum()/1e9))
        print('Compare to: 1000 kilotons (Meijer), 6 kilotons (Weiss)')
        print('Total: %3.3f kilotons per year from rivers, pop, fisheries' % (np.array([vec_start[i].sum() for i in range(12)]).sum()/1e9))
    
    #t1 = print_time(t1,'sources init.')
    #-------------------------set-up matrix-------------------
    i_r_ll,i_c_ll,_ = sparse.find(P['land_land'])
    cols_land = np.unique(i_c_ll)
    P_land_removed = coo_matrix( (np.ones(len(cols_land)),(2*np.ones(len(cols_land)),cols_land)), 
                                shape=(n_total,n_total))# * (1-p_land_remain)
    
    i_r_o,i_c_o,_ = sparse.find(P['water_water']+P['coast_water'])
    # cols_ocean = np.unique(i_c_o)
    
    i_r_ww,i_c_ww,_ = sparse.find(P['water_water'])
    i_r_cw,i_c_cw,_ = sparse.find(P['coast_water'])
    i_r_cl,i_c_cl,_ = sparse.find(P['coast_land'])
    
    #t1 = print_time(t1,'mat init')
    # -------------------sinking parameterization-------------------------
    # pp_water_water = P['pp_water_water']
    # pp_coast_water = P['pp_coast_water']
    # pp_coast_land = P['pp_coast_land']
    # P_ocean_removed = P['P_ocean_removed']

        
    # -------------------beaching parameterization-------------------
    p_coast_cl = lin_max(map_h3_to_coast_length[map_mat_to_h3[i_c_cl]],l_max)
    p_coast_cw = lin_max(map_h3_to_coast_length[map_mat_to_h3[i_c_cw]],l_max)
    
    P_cascade_m = {}
    P_cascade_m[0] = P['cascade_water'] + P['cascade_land']*p_frag_m[0]
    for i3 in range(1,len(k_arr)):
        P_cascade_m[i3] = P['cascade_land']*p_frag_m[i3]
    
    P_cascade_N = {}
    P_cascade_N[0] = P['cascade_water'] + P['cascade_land']*p_frag_N[0]
    for i3 in range(1,len(k_arr)):
        P_cascade_N[i3] = P['cascade_land']*p_frag_N[i3]
    
    #t1 = print_time(t1,'casc')
    
    P_total_l = {}
    for month_ in np.arange(12):
        P_ocean_l = {}
        P_total_l[month_] = {}
        available_lengths = np.sort(np.array(list(P['ocean'][month_]['nonfouled'].keys())))
    
    
        def l_to_leff(l_,n_pad=5):
            # for the rise velocity, take a linear combination for different particle sizes (i.e. w_b's)
            # as to try to solve the problem of too little mixing. A Gaussian function is used to weigh the samples
            n_pad = 5
            assert(available_lengths[1]/available_lengths[0] == 4.)
            left_ = np.sort(np.array([available_lengths[0] / (4**n) for n in range(1,n_pad+1)]))
            right_ = np.sort(np.array([available_lengths[-1] * (4**n) for n in range(1,n_pad+1)]))
            l_arr_padded = np.append(left_,np.append(available_lengths,right_))
    
            mu_L_eff = a_L_Leff*np.log10(l_) + b_L_Leff
            weights_Leff_padded = Gaussian(np.log10(l_arr_padded),mu_L_eff,sigma_L_eff)
            weights_Leff = weights_Leff_padded[n_pad:-n_pad].copy()
            weights_Leff[0] += weights_Leff_padded[:n_pad].sum()
            weights_Leff[-1] += weights_Leff_padded[-n_pad:].sum()
            weights_Leff /= weights_Leff.sum() 
            l_use = available_lengths[weights_Leff > 0]
            weight_use = weights_Leff[weights_Leff > 0]
    
            if not len(weight_use)>0:
                #exception where the the spread gets too small -> no weights
                l_use = np.array([l_])
                weight_use = np.array([1])
    
            return l_use, weight_use
    
    
        P_ocean_nb = P['ocean'][month_]['neutral'][1e-5] #the 1e-5 is the arbitrary small size denotation that was used in this case
        for i1,l_ in enumerate(l_arr):
    
            l_use,weight_use = l_to_leff(l_)
    
            #------------------------------------------------------------------------
            # non-fouled bouyant particles: create a combination of rise velocities
            P_ocean_nonfouled = P['ocean'][month_]['nonfouled'][l_use[0]]*weight_use[0]
            if len(weight_use) > 1:
                for i2 in range(1,len(weight_use)):
                    P_ocean_nonfouled += P['ocean'][month_]['nonfouled'][l_use[i2]]*weight_use[i2]
            #t1 = print_time(t1,'non-fouled')
            #------------------------------------------------------------------------
            # fouled bouyant particles: take the closest fouled case for transport
            i_right = np.searchsorted(available_lengths,l_)
            if i_right == 0:
                i_right = 1
            if i_right == len(available_lengths): #above max avail. length: assume same behavior for macroplastics
                l_max = available_lengths.max()
                P_ocean_pf = P['ocean'][month_][pf_mode][l_max]
                P_ocean_fouled = P['ocean'][month_]['fouled'][l_max] 
            else:
                l_left = available_lengths[i_right-1]
                l_right = available_lengths[i_right]
                w_2 = (np.log10(l_) - np.log10(l_left)) / (np.log10(l_right) - np.log10(l_left)) #0 at 0.0004, 1 at 0.102
                assert(0<=w_2<=1)
                P_ocean_fouled =  P['ocean'][month_]['fouled'][l_left]*(1-w_2) + P['ocean'][month_]['fouled'][l_right]*(w_2)
                P_ocean_pf = P['ocean'][month_][pf_mode][l_left]*(1-w_2) + P['ocean'][month_][pf_mode][l_right]*(w_2)
            #t1 = print_time(t1,'fouled')
            #------------------------------------------------------------------------
            # A linear combination is made of fouled (oscill.), non-fouled, neutrally bouyant, and permanently fouled particles
            P_ocean_l[l_] = P_ocean_nonfouled*(1-frac_fouled-frac_nb-frac_pf) + P_ocean_fouled*frac_fouled + P_ocean_nb*frac_nb + P_ocean_pf*frac_pf
    
            #t1 = print_time(t1,'combin.')
            #------------------------------------------------------------------------
            # Surface transport can be adjusted, using bi-linear interpolation with wind and Stokes
            # frac_Stokes = lin_max(l_,l_max_Stokes,max_val=.99)
            # frac_wind = lin_max(l_,l_max_wind,max_val=.99)
            # # frac_Stokes = .99
            # # frac_wind = .99
            # P_ocean_l[l_] = (P_ocean_l[l_]*(1-frac_Stokes) + P['ocean'][month_]['2D_Stokes'][1.0]*(frac_Stokes) )*(1 - frac_wind) + \
            #                 (P['ocean'][month_]['2D_wind'][1.0]*(1-frac_Stokes) + P['ocean'][month_]['2D_Stokes_wind'][1.0]*(frac_Stokes) )*(frac_wind)
    
            # # needs to be normalized again to retain mass balance in the deeper water layer matrix columns
            # P_ocean_l[l_] = normalize(P_ocean_l[l_], norm='l1', axis=0)
            
            #t1 = print_time(t1,'surf+norm')
            
            p_beach_per_cell_cl = p_beach[i1]*p_coast_cl
            p_beach_per_cell_cw = 1-(p_beach[i1]*p_coast_cw)
            P_coast_land_beach = coo_matrix( (p_beach_per_cell_cl,(i_r_cl,i_c_cl)), shape=(n_total,n_total))
            P_coast_water_beach = coo_matrix( (p_beach_per_cell_cw,(i_r_cw,i_c_cw)), shape=(n_total,n_total))
            # P_coast_water_beach = P['coast_water']*(1-p_beach[i1]) 
            # P_coast_land_beach = P['coast_land']*p_beach[i1]
            P_land_coast_beach = P['land_coast']*p_resus[i1]
            P_land_land_beach = P['land_land']*(1-p_resus[i1])
            
            #t1 = print_time(t1,'beach')
            
            #cells on land have permanent sinks
            P_land_coast_sinks = P_land_coast_beach.multiply(p_land_remain)
            P_land_land_sinks = P_land_land_beach.multiply(p_land_remain)
    
            # sinking due to primary productivity
            # A linear map is made from min pp (0 -> 1), to max pp (1 -> p_ocean_min_remain), such that no particles sink when there's no pp, and the prescribed max amount at pp=1
            # P_water_water_remain = (p_ocean_min_remain - 1)*pp_water_water[month_] + P['water_water'] 
            # P_coast_water_remain = (p_ocean_min_remain - 1)*pp_coast_water[month_] + P['coast_water'] 
            # P_coast_land_remain = (p_ocean_min_remain - 1)*pp_coast_land[month_] + P['coast_land'] 

            P_water_water_remain = P['water_water'] 
            P_coast_water_remain = P['coast_water'] 
            P_coast_land_remain = P['coast_land'] 
            #t1 = print_time(t1,'pp')
            
            #TO DO: size-dependent bottom hitting can be added here as well for fouled/non-fouled particles
            P_water_water = P_water_water_remain
            P_coast_water = P_coast_water_remain.multiply(P_coast_water_beach) #probabilities after fouling and beaching
            P_coast_land = P_coast_land_remain.multiply(P_coast_land_beach)
    
    
            # P_total_l[month_][l_] = (P_ocean_l[l_].multiply(P_water_water + P_coast_water)) + P_coast_land + \
            # P_land_coast_sinks + P_land_land_sinks + P_land_removed*(1-p_land_remain) + P_ocean_removed[month_]*(1-p_ocean_min_remain) + P['removal']
            
            P_total_l[month_][l_] = (P_ocean_l[l_].multiply(P_water_water + P_coast_water)) + P_coast_land + \
            P_land_coast_sinks + P_land_land_sinks + P_land_removed*(1-p_land_remain) + P['removal']
            
            #t1 = print_time(t1,'total')
            colsum = np.array(np.sum(P_total_l[month_][l_],axis=0)).ravel()
            leakiness = 1 - (colsum.sum()/len(colsum))
            if leakiness > 0.0001: #check that mass conservation is not violated too much
                print('leakiness bigger than 0.01 percent (m:%i l:%4.4f): %f' % (month_,l_,leakiness))
    
    #-------------------------Iterate-------------------
    
    # df_aggregated = pd.DataFrame(columns=['modelled_n','modelled_sum','modelled_sum2'])
    df_results = pd.DataFrame(columns=['Year','Month',
                                       'h3_cell_id','Modelled','Measured',
                                       'Units','var_log10','l_min','l_max',
                                       'ParentEventID','index_parent','Depth','detection_limit'])
    
    extent = (-179,179,-80,80)
    
    
    plot_map = False
    do_plot_dist = False
    plot_tracer_amount = False
    plot_per = 120#120*2#12*10
    i_start = 0
    n_rows = 1
    n_cols = 1
    n_steps = (2020-1980)*12#12*11
    
    
    vmin = 1e-9
    vmax = 1e-3
    cmap = plt.cm.viridis
    vmin2 = 1e-5 #beaches
    vmax2 = 1e1
    cmap2 = cmocean.cm.turbid_r
    
    tracer_matrix_mass = np.zeros([n_total,len(k_arr)],dtype=np.float32)
    tracer_matrix_num = np.zeros([n_total,len(k_arr)],dtype=np.float32)
    tracer_matrix_num_f = np.zeros([n_total,len(k_arr)],dtype=np.float32)
    
    if plot_map:
        fig = plt.figure(figsize=(20,9))
        gs = gridspec.GridSpec(1, 4, bottom=.05, top=.95, wspace=.25, width_ratios=[30,30,1,1])
        ax_1 = plt.subplot(gs[0,0])
        ax_3 = plt.subplot(gs[0,2])
        cplt1 = ax_1.scatter(np.random.random(1000),np.random.random(1000),c=np.linspace(np.log10(vmin),np.log10(vmax),1000),cmap=cmap)
        cbar1 = plt.colorbar(cplt1,cax=ax_3)
        cbar1.set_label('Tracer concentration [log$_{10}$(g m$^{-3}$)]')
        ax_2 = plt.subplot(gs[0,1])
        ax_4 = plt.subplot(gs[0,3])
        cplt2 = ax_2.scatter(np.random.random(1000),np.random.random(1000),c=np.linspace(np.log10(vmin2),np.log10(vmax2),1000),cmap=cmap2)
        cbar2 = plt.colorbar(cplt2,cax=ax_4)
        cbar2.set_label('Tracer concentration [log$_{10}$(g m$^{-1}$)]')
        fig.savefig(folder_figure+'colorbar.png')
    
        fig = plt.figure(figsize=(25, 9))
        gs = gridspec.GridSpec(n_rows, n_cols, bottom=.05, top=.95, wspace=.1)
        for i4,l_ in enumerate(l_arr[0:n_cols]):
            str_size = 'year %i' % (0)
            # draw_map(tracer_matrix_mass.sum(axis=1)/dict_h3_normalize.values,extent,gs,i4,base_title=str_size+', %3.0f-%3.0f m depth',stop_earlier=(len(depth_layers)-n_rows-1) )
            draw_map(tracer_matrix_mass[:,i4]/dict_h3_normalize.values,extent,gs,i4,base_title=str_size+', %3.0f-%3.0f m depth',stop_earlier=(len(depth_layers)-n_rows-1) )
        if do_plot_dist:
            ax_ = plt.subplot(gs[0,0], projection=ccrs.PlateCarree())
            for i1,(lon_,lat_) in enumerate(zip(lon_plot_dist,lat_plot_dist)):
                ax_.plot(lon_,lat_,'x',markersize=10,color=plt.cm.tab10(i1))
            
        fig.savefig(folder_figure+'%3.3i.png' % (0))
        plt.close('all')
        
    gc.collect()
    # gc.disable()
        
    for i1 in range(n_steps): # time loop
        t0 = time()
        month_ = i1 % 12
        
        year_current = year_start + i1//12
        month_current = (i1 % 12) + 1 # from 0-11 to 1-12
        date_current = datetime(year_current,month_current,1) 
        
        year_from_2015 = year_current - 2015
        frac_PW = np.exp(PW_factor*year_from_2015)
        
        # if release == 'simple' and i1 == 0:
        #     tracer_matrix_mass[:,0] += vec_start
        # elif release == 'mix':
        #     tracer_matrix_mass[:,0] += vec_start[month_,:]
        #     tracer_matrix_num[:,0] += vec_start_N[month_,:]
        #     tracer_matrix_num_f[:,0] += vec_start_fN[month_,:]        
        # else:
        #     tracer_matrix_mass[:,0] += vec_start[month_,:]
        if release == 'simple' and i1 == 0:
            tracer_matrix_mass += vec_start[month_]
        elif release == 'mix':
            tracer_matrix_mass += vec_start[month_]*frac_PW
            tracer_matrix_num += vec_start_N[month_]*frac_PW
            tracer_matrix_num_f += vec_start_fN[month_]*frac_PW       
        else:
            tracer_matrix_mass += vec_start[month_]      
        # t01 = time()
        
        for i2,l_ in enumerate(l_arr): # particle size loop
            
            # t02 = time()
            # first, calculate how much material is not fragmented
            tracer_matrix_mass[:,i2] = P_cascade_m[0].dot(tracer_matrix_mass[:,i2])
            tracer_matrix_num[:,i2] = P_cascade_N[0].dot(tracer_matrix_num[:,i2])
            tracer_matrix_num_f[:,i2] = P_cascade_N[0].dot(tracer_matrix_num_f[:,i2])
            
            # t002 = time()
            # then, let some of the particles cascade into the next size classes.
            # start at size class i2+1
            for i3 in range(i2+1,len(k_arr)):
                tracer_matrix_mass[:,i3] += P_cascade_m[i3-i2].dot(tracer_matrix_mass[:,i2]) #i3-i2 is used here, since this number should run from 1 upwards
                tracer_matrix_num[:,i3] += P_cascade_N[i3-i2].dot(tracer_matrix_num[:,i2]) #i3-i2 is used here, since this number should run from 1 upwards
                tracer_matrix_num_f[:,i3] += P_cascade_N[i3-i2].dot(tracer_matrix_num_f[:,i2]) #i3-i2 is used here, since this number should run from 1 upwards
     
            # t003 = time()
            # then, advect the tracer
            tracer_matrix_mass[:,i2] = P_total_l[month_][l_].dot( tracer_matrix_mass[:,i2] )
            tracer_matrix_num[:,i2] = P_total_l[month_][l_].dot( tracer_matrix_num[:,i2] )
            tracer_matrix_num_f[:,i2] = P_total_l[month_][l_].dot( tracer_matrix_num_f[:,i2] )
            
            # t004 = time()
        mass_total = tracer_matrix_mass.sum()
        mass_removed = (tracer_matrix_mass[0:2,:].sum(axis=1) / mass_total)*100
        mass_bottom = (tracer_matrix_mass[2,:].sum() / mass_total)*100
    
        tracer_matrix_mass[0:3,:] = 0.
        tracer_matrix_num[0:3,:] = 0.
        tracer_matrix_num_f[0:3,:] = 0.
               
            
            # t005 = time()
            # print('-------')
            # print(t002-t02,t003-t002,t004-t003,t005-t004)
            # print(time()-t02)
            
        t1 = time()
        # print(t01-t0,t1-t02)
        # print('-------------------------------------')
    
        # start = 01 jan
        # end = 31 jan
        # compare all measurements in jan to 31 jan
        

        
        """
        df_aggregated = pd.DataFrame(columns=['modelled_n','modelled_sum','modelled_sum2'])
        df_results = pd.DataFrame(columns=['Year','Month',
                                           'h3_cell_id','Modelled','Measured',
                                           'Units','l_min','l_max','ParentEventID','detection_limit'])
        """
        mask_select = (data_trawl_m['Year'] == year_current) & (data_trawl_m['Month'] == month_current)
        if mask_select.sum() > 0 and mask_neuston_net.sum() > 0:
            # print('Adding %i results' % mask_select.sum())
            h3_selection = data_trawl_m.loc[mask_select,'h%i'%h3_grid_res]
            tracer_total_mass = tracer_matrix_mass[map_h3_to_mat[h3_selection]][:,mask_neuston_net].sum(axis=1)
            tracer_concentration_mass = tracer_total_mass / dict_h3_normalize[h3_selection]
            tracer_concentration_mass_meas = data_trawl_m['measurementValue_vol_5m'][mask_select]
            l_min = l_arr[mask_neuston_net][-1]
            l_max = l_arr[mask_neuston_net][0]
            
            df_tmp = pd.DataFrame({'Year': year_current,
                                   'Month': month_current,
                                   'h3_cell_id': h3_selection.values,
                                   'Modelled': tracer_concentration_mass,
                                   'Measured': tracer_concentration_mass_meas.values,
                                   'Units': 'g m-3',
                                   'var_log10': dict_variance['m_surf'],
                                   'l_min': l_min,
                                   'l_max': l_max,
                                   'ParentEventID': data_trawl_m['ParentEventID'][mask_select].values,
                                   'index_parent':data_trawl_m.index[mask_select].values,
                                   'Depth': data_trawl_m['Depth'][mask_select].values,
                                   'detection_limit': 0})
            df_results = pd.concat((df_results,df_tmp), ignore_index=True)
            
        mask_select = (data_trawl_n['Year'] == year_current) & (data_trawl_n['Month'] == month_current)
        if mask_select.sum() > 0 and mask_neuston_net.sum() > 0:
            # print('Adding %i results' % mask_select.sum())
            h3_selection = data_trawl_n.loc[mask_select,'h%i'%h3_grid_res]
            tracer_total_num = tracer_matrix_num[map_h3_to_mat[h3_selection]][:,mask_neuston_net].sum(axis=1)
            tracer_concentration_num = tracer_total_num / dict_h3_normalize[h3_selection]        
            tracer_concentration_num_meas = data_trawl_n['measurementValue_vol_5m'][mask_select]
            l_min = l_arr[mask_neuston_net][-1]
            l_max = l_arr[mask_neuston_net][0]
            
            df_tmp = pd.DataFrame({'Year': year_current,
                                   'Month': month_current,
                                   'h3_cell_id': h3_selection.values,
                                   'Modelled': tracer_concentration_num,
                                   'Measured': tracer_concentration_num_meas.values,
                                   'Units': 'n m-3',
                                   'var_log10': dict_variance['n_surf'],
                                   'l_min': l_min,
                                   'l_max': l_max,
                                   'ParentEventID': data_trawl_n['ParentEventID'][mask_select].values,
                                   'index_parent':data_trawl_n.index[mask_select].values,
                                   'Depth': data_trawl_n['Depth'][mask_select].values,
                                   'detection_limit': 0})
            df_results = pd.concat((df_results,df_tmp), ignore_index=True)
                
            
        mask_select = (data_macro_n['Year'] == year_current) & (data_macro_n['Month'] == month_current)
        if mask_select.sum() > 0:
            size_bnds = np.unique(np.array([data_macro_n.loc[mask_select,'l_min'],data_macro_n.loc[mask_select,'l_max']]),axis=1)
        
            for size_min,size_max in zip(size_bnds[0,:],size_bnds[1,:]):
                mask_size_class = (l_arr >= size_min) & (l_arr <= size_max)
                mask_size_measurement = data_macro_n['l_min'] == size_min
                    
                if mask_size_class.sum() > 0:
                    h3_selection = data_macro_n.loc[(mask_select&mask_size_measurement),'h%i'%h3_grid_res] 
        
                    tracer_total = tracer_matrix_num[map_h3_to_mat[h3_selection]][:,mask_size_class].sum(axis=1)
                    tracer_concentration = tracer_total / dict_h3_normalize[h3_selection]        
                    tracer_concentration_meas = data_macro_n['measurementValue_vol_5m'][(mask_select&mask_size_measurement)]
        
                    df_tmp = pd.DataFrame({'Year': year_current,
                                           'Month': month_current,
                                           'h3_cell_id': h3_selection.values,
                                           'Modelled': tracer_concentration,
                                           'Measured': tracer_concentration_meas.values,
                                           'Units': 'n m-3',
                                           'var_log10': dict_variance['n_surf'],
                                           'l_min': size_min,
                                           'l_max': size_max,
                                           'ParentEventID': data_macro_n['ParentEventID'][(mask_select&mask_size_measurement)].values,
                                           'index_parent':data_macro_n.index[(mask_select&mask_size_measurement)].values,
                                           'Depth': data_macro_n['Depth'][(mask_select&mask_size_measurement)].values,
                                           'detection_limit': 0})
                    df_results = pd.concat((df_results,df_tmp), ignore_index=True)     
                   
                
        mask_select = (data_beach_inst['Year'] == year_current) & (data_beach_inst['Month'] == month_current)
        if mask_select.sum() > 0 and mask_beach_cleanup.sum() > 0:
            
            for unit_ in np.unique(data_beach_inst.loc[mask_select,'measurementUnit'].values):
                mask_units = data_beach_inst['measurementUnit'] == unit_
                # print('Adding %i results' % mask_select.sum())
                h3_selection = data_beach_inst.loc[(mask_select & mask_units),'h%i'%h3_grid_res] + index_add_coast
                
                if unit_ == 'g m-1':
                    tracer_mat = tracer_matrix_mass
                    label_var = 'm_beach'
                elif unit_ == 'n m-1':
                    tracer_mat = tracer_matrix_num
                    label_var = 'n_beach'
                else:
                    raise RuntimeError('unknown units %s' % unit_)
                
                tracer_total = tracer_mat[map_h3_to_mat[h3_selection]][:,mask_beach_cleanup].sum(axis=1)
                tracer_concentration = tracer_total / dict_h3_normalize[h3_selection]        
                tracer_concentration_meas = data_beach_inst['measurementValue'][(mask_select & mask_units)]
                l_min = l_arr[mask_beach_cleanup][-1]
                l_max = l_arr[mask_beach_cleanup][0]
                
                df_tmp = pd.DataFrame({'Year': year_current,
                                       'Month': month_current,
                                       'h3_cell_id': h3_selection.values,
                                       'Modelled': tracer_concentration,
                                       'Measured': tracer_concentration_meas.values,
                                       'Units': unit_,
                                       'var_log10': dict_variance[label_var],
                                       'l_min': l_min,
                                       'l_max': l_max,
                                       'ParentEventID': data_beach_inst['ParentEventID'][(mask_select & mask_units)].values,
                                       'index_parent':data_beach_inst.index[(mask_select & mask_units)].values,
                                       'Depth': data_beach_inst['Depth'][(mask_select & mask_units)].values,
                                       'detection_limit': 0})
                df_results = pd.concat((df_results,df_tmp), ignore_index=True)
            
            
        def add_OSPAR_results(year_current,month_current,key_,df_results):
            
            mask_select = (data_beach_OSPAR['Year'] == year_current) & (data_beach_OSPAR['Month'] == month_current) & (data_beach_OSPAR['ParentEventID'] == key_)
            
            if mask_select.sum() > 0 and mask_beach_cleanup.sum() > 0:
                    
                h3_selection = data_beach_OSPAR.loc[mask_select,'h%i'%h3_grid_res] + index_add_coast
                size_ = np.unique(data_beach_OSPAR.loc[mask_select,'l_min'])
                assert(len(size_) == 1)
                
                mask_size_class = l_arr > size_
                # TO DO: remove this. This is set such that macroplastic classes can still be compared
                # if mask_size_class.sum() == 0:
                #     mask_size_class[0] = True
                    
                tracer_total_num = tracer_matrix_num[map_h3_to_mat[h3_selection]][:,mask_size_class].sum(axis=1)
                tracer_total_fis = tracer_matrix_num_f[map_h3_to_mat[h3_selection]][:,mask_size_class].sum(axis=1)
                
                tracer_concentration_num = tracer_total_num / dict_h3_normalize[h3_selection]        
                tracer_concentration_num_meas = data_beach_OSPAR['measurementValue'][mask_select]
                tracer_fraction_fis = tracer_total_fis / tracer_total_num 
                mask_invalid = (np.isnan(tracer_fraction_fis) | np.isinf(tracer_fraction_fis))
                tracer_fraction_fis[mask_invalid] = 0
                tracer_fraction_fis_meas = data_beach_OSPAR['fishing_fraction'][mask_select]
                
                l_min = l_arr[mask_size_class][-1]
                l_max = l_arr[mask_size_class][0]
                
                df_tmp = pd.DataFrame({'Year': year_current,
                                       'Month': month_current,
                                       'h3_cell_id': h3_selection.values,
                                       'Modelled': tracer_concentration_num,
                                       'Measured': tracer_concentration_num_meas.values,
                                       'Units': 'n m-1',
                                       'var_log10': dict_variance[key_],
                                       'l_min': l_min,
                                       'l_max': l_max,
                                       'ParentEventID': data_beach_OSPAR['ParentEventID'][mask_select].values,
                                       'index_parent':data_beach_OSPAR.index[mask_select].values,
                                       'Depth': data_beach_OSPAR['Depth'][mask_select].values,
                                       'detection_limit': 0})
                df_results = pd.concat((df_results,df_tmp), ignore_index=True)    
                
                #fishing fraction
                df_tmp = pd.DataFrame({'Year': year_current,
                                       'Month': month_current,
                                       'h3_cell_id': h3_selection.values,
                                       'Modelled': 10**(tracer_fraction_fis), #10** since the log is taken in cost fn, want to compare standard fractions
                                       'Measured': 10**(tracer_fraction_fis_meas.values),
                                       'Units': 'n n-1',
                                       'var_log10': dict_variance['fish_frac'],
                                       'l_min': l_min,
                                       'l_max': l_max,
                                       'ParentEventID': data_beach_OSPAR['ParentEventID'][mask_select].values,
                                       'index_parent':data_beach_OSPAR.index[mask_select].values,
                                       'Depth': data_beach_OSPAR['Depth'][mask_select].values,
                                       'detection_limit': 0})
                df_results = pd.concat((df_results,df_tmp), ignore_index=True)  
            return df_results
        
        df_results = add_OSPAR_results(year_current,month_current,'MDMAP',df_results)
        # df_results = add_OSPAR_results(year_current,month_current,'OSPAR_100m',df_results)
        df_results = add_OSPAR_results(year_current,month_current,'OSPAR_1km',df_results)
                                
            
        def add_beach_mean_results(date_current,data_beach_mean):        
                
            mask_select = (date_current >= data_beach_mean['eventDate']) & (date_current < data_beach_mean['eventDate_end'])
            if mask_select.sum() > 0 and mask_beach_cleanup.sum() > 0: #mask_beach_cleanup: perhaps not necessay since l_min/max are given below
                    
                size_bnds = np.unique(np.array([data_beach_mean.loc[mask_select,'l_min'],data_beach_mean.loc[mask_select,'l_max']]),axis=1)
                
                for size_min,size_max in zip(size_bnds[0,:],size_bnds[1,:]):
                    mask_size_class = (l_arr >= size_min) & (l_arr <= size_max)
                    
                    for unit_ in np.unique(data_beach_mean.loc[mask_select,'measurementUnit'].values):
                        mask_units = data_beach_mean['measurementUnit'] == unit_
                        h3_selection = data_beach_mean.loc[(mask_select & mask_units),'h%i'%h3_grid_res] + index_add_coast    
                
                        if unit_ == 'g m-1':
                            tracer_mat = tracer_matrix_mass
                        elif unit_ == 'n m-1':
                            tracer_mat = tracer_matrix_num
                        else:
                            raise RuntimeError('unknown units %s' % unit_)
                        
                        tracer_total = tracer_mat[map_h3_to_mat[h3_selection]][:,mask_size_class].sum(axis=1)
                        tracer_concentration = tracer_total / dict_h3_normalize[h3_selection]        
                
                        data_beach_mean.loc[(mask_select & mask_units),'modelled_n'] += 1
                        data_beach_mean.loc[(mask_select & mask_units),'modelled_sum'] += tracer_concentration.values
                        data_beach_mean.loc[(mask_select & mask_units),'modelled_sum2'] += tracer_concentration.values**2 
                        
        add_beach_mean_results(date_current,data_beach_mean)
        
            
        mask_select = (data_Egger_surf['Year'] == year_current) & (data_Egger_surf['Month'] == month_current)
        if mask_select.sum() > 0:
            # print('Adding %i results' % mask_select.sum())
            for i3,mask_size_range in enumerate(mask_Egger):
                if mask_size_range.sum() > 0:                 
                    size_range = mask_Egger_sizes[i3]
                    corr_fac_n = Egger_corr_n[i3]
                    corr_fac_m = Egger_corr_m[i3]
                    
                    col_name = 'measurementValue_vol_%im_%5.5fm_%5.5fm' % (depth_layers[1],size_range[0],size_range[1])
                    
                    mask_units = data_Egger_surf['measurementUnit_vol'] == 'nm-3'
                    h3_selection = data_Egger_surf.loc[(mask_select&mask_units),'h%i'%h3_grid_res]
                    tracer_total_num = tracer_matrix_num[map_h3_to_mat[h3_selection]][:,mask_size_range].sum(axis=1)
                    tracer_concentration_num = (tracer_total_num / dict_h3_normalize[h3_selection])*corr_fac_n        
                    tracer_concentration_num_meas = data_Egger_surf[col_name][(mask_select&mask_units)]
                    l_min = l_arr[mask_size_range][-1]
                    l_max = l_arr[mask_size_range][0]
                    
                    df_tmp = pd.DataFrame({'Year': year_current,
                                           'Month': month_current,
                                           'h3_cell_id': h3_selection.values,
                                           'Modelled': tracer_concentration_num,
                                           'Measured': tracer_concentration_num_meas.values,
                                           'Units': 'n m-3',
                                           'var_log10': dict_variance['surf_nm-3_Egger'],
                                           'l_min': l_min,
                                           'l_max': l_max,
                                           'ParentEventID': data_Egger_surf['ParentEventID'][(mask_select&mask_units)].values,
                                           'index_parent':data_Egger_surf.index[(mask_select & mask_units)].values,
                                           'Depth': data_Egger_surf['Depth'][(mask_select & mask_units)].values,
                                           'detection_limit': 0})
                    df_results = pd.concat((df_results,df_tmp), ignore_index=True)            
    
                    mask_units = data_Egger_surf['measurementUnit_vol'] == 'gm-3'
                    h3_selection = data_Egger_surf.loc[(mask_select&mask_units),'h%i'%h3_grid_res]
                    tracer_total_mass = tracer_matrix_mass[map_h3_to_mat[h3_selection]][:,mask_size_range].sum(axis=1)
                    tracer_concentration_mass = (tracer_total_mass / dict_h3_normalize[h3_selection])*corr_fac_m        
                    tracer_concentration_mass_meas = data_Egger_surf[col_name][(mask_select&mask_units)]
                    l_min = l_arr[mask_size_range][-1]
                    l_max = l_arr[mask_size_range][0]
                    
                    df_tmp = pd.DataFrame({'Year': year_current,
                                           'Month': month_current,
                                           'h3_cell_id': h3_selection.values,
                                           'Modelled': tracer_concentration_mass,
                                           'Measured': tracer_concentration_mass_meas.values,
                                           'Units': 'g m-3',
                                           'var_log10': dict_variance['surf_gm-3_Egger'],
                                           'l_min': l_min,
                                           'l_max': l_max,
                                           'ParentEventID': data_Egger_surf['ParentEventID'][(mask_select&mask_units)].values,
                                           'index_parent':data_Egger_surf.index[(mask_select & mask_units)].values,
                                           'Depth': data_Egger_surf['Depth'][(mask_select & mask_units)].values,
                                           'detection_limit': 0})
                    df_results = pd.concat((df_results,df_tmp), ignore_index=True)    
    
        mask_select = (data_Egger_depth['Year'] == year_current) & (data_Egger_depth['Month'] == month_current)
        if mask_select.sum() > 0:
            # print('Adding %i results' % mask_select.sum())
            for i3,mask_size_range in enumerate(mask_Egger):
                if mask_size_range.sum() > 0: 
                    size_range = mask_Egger_sizes[i3]
                    corr_fac_n = Egger_corr_n[i3]
                    corr_fac_m = Egger_corr_m[i3]
                    
                    col_name = 'measurementValue_vol_%5.5fm_%5.5fm' % (size_range[0],size_range[1])
                    
                    mask_units = data_Egger_depth['measurementUnit_vol'] == 'nm-3'
                    h3_selection = data_Egger_depth.loc[(mask_select&mask_units),'h%i'%h3_grid_res]
                    tracer_total_num = tracer_matrix_num[map_h3_to_mat[h3_selection]][:,mask_size_range].sum(axis=1)
                    tracer_concentration_num = (tracer_total_num / dict_h3_normalize[h3_selection])*corr_fac_n        
                    tracer_concentration_num_meas = data_Egger_depth[col_name][(mask_select&mask_units)]
                    tracer_detection_limit = data_Egger_depth[col_name+'_detlim.'][(mask_select&mask_units)]
                    l_min = l_arr[mask_size_range][-1]
                    l_max = l_arr[mask_size_range][0]
                    
                    df_tmp = pd.DataFrame({'Year': year_current,
                                           'Month': month_current,
                                           'h3_cell_id': h3_selection.values,
                                           'Modelled': tracer_concentration_num,
                                           'Measured': tracer_concentration_num_meas.values,
                                           'Units': 'n m-3',
                                           'var_log10': dict_variance['surf_nm-3_Egger'],
                                           'l_min': l_min,
                                           'l_max': l_max,
                                           'ParentEventID': data_Egger_depth['ParentEventID'][(mask_select&mask_units)].values,
                                           'index_parent':data_Egger_depth.index[(mask_select & mask_units)].values,
                                           'Depth': data_Egger_depth['Depth'][(mask_select & mask_units)].values,  
                                           'detection_limit': tracer_detection_limit.values})
                    df_results = pd.concat((df_results,df_tmp), ignore_index=True)            
    
                    mask_units = data_Egger_depth['measurementUnit_vol'] == 'gm-3'
                    h3_selection = data_Egger_depth.loc[(mask_select&mask_units),'h%i'%h3_grid_res]
                    tracer_total_mass = tracer_matrix_mass[map_h3_to_mat[h3_selection]][:,mask_size_range].sum(axis=1)
                    tracer_concentration_mass = (tracer_total_mass / dict_h3_normalize[h3_selection])*corr_fac_m        
                    tracer_concentration_mass_meas = data_Egger_depth[col_name][(mask_select&mask_units)]
                    tracer_detection_limit = data_Egger_depth[col_name+'_detlim.'][(mask_select&mask_units)]
                    l_min = l_arr[mask_size_range][-1]
                    l_max = l_arr[mask_size_range][0]
                    
                    df_tmp = pd.DataFrame({'Year': year_current,
                                           'Month': month_current,
                                           'h3_cell_id': h3_selection.values,
                                           'Modelled': tracer_concentration_mass,
                                           'Measured': tracer_concentration_mass_meas.values,
                                           'Units': 'g m-3',
                                           'var_log10': dict_variance['surf_gm-3_Egger'],
                                           'l_min': l_min,
                                           'l_max': l_max,
                                           'ParentEventID': data_Egger_depth['ParentEventID'][(mask_select&mask_units)].values,
                                           'index_parent':data_Egger_depth.index[(mask_select & mask_units)].values,
                                           'Depth': data_Egger_depth['Depth'][(mask_select & mask_units)].values,  
                                           'detection_limit': tracer_detection_limit.values})
                    df_results = pd.concat((df_results,df_tmp), ignore_index=True)    
    
        mask_select = (data_Zhao_depth['Year'] == year_current) & (data_Zhao_depth['Month'] == month_current)
        if mask_select.sum() > 0 and mask_Zhao.sum()>0:
            # print('Adding %i results' % mask_select.sum())
            col_name = 'measurementValue_vol_0.00030m_0.00500m'
            
            mask_units = data_Zhao_depth['measurementUnit_vol'] == 'nm-3'
            h3_selection = data_Zhao_depth.loc[(mask_select&mask_units),'h%i'%h3_grid_res]
            tracer_total_num = tracer_matrix_num[map_h3_to_mat[h3_selection]][:,mask_Zhao].sum(axis=1)
            tracer_concentration_num = tracer_total_num / dict_h3_normalize[h3_selection]        
            tracer_concentration_num_meas = data_Zhao_depth[col_name][(mask_select&mask_units)]
            tracer_detection_limit = data_Zhao_depth[col_name+'_detlim.'][(mask_select&mask_units)]
            l_min = l_arr[mask_Zhao][-1]
            l_max = l_arr[mask_Zhao][0]
            
            df_tmp = pd.DataFrame({'Year': year_current,
                                   'Month': month_current,
                                   'h3_cell_id': h3_selection.values,
                                   'Modelled': tracer_concentration_num,
                                   'Measured': tracer_concentration_num_meas.values,
                                   'Units': 'n m-3',
                                   'var_log10': dict_variance['n_surf'],
                                   'l_min': l_min,
                                   'l_max': l_max,
                                   'ParentEventID': data_Zhao_depth['ParentEventID'][(mask_select&mask_units)].values,
                                   'index_parent':data_Zhao_depth.index[(mask_select & mask_units)].values,
                                   'Depth': data_Zhao_depth['Depth'][(mask_select & mask_units)].values,  
                                   'detection_limit': tracer_detection_limit.values})
            df_results = pd.concat((df_results,df_tmp), ignore_index=True)            
    
            mask_units = data_Zhao_depth['measurementUnit_vol'] == 'gm-3'
            h3_selection = data_Zhao_depth.loc[(mask_select&mask_units),'h%i'%h3_grid_res]
            tracer_total_mass = tracer_matrix_mass[map_h3_to_mat[h3_selection]][:,mask_Zhao].sum(axis=1)
            tracer_concentration_mass = tracer_total_mass / dict_h3_normalize[h3_selection]        
            tracer_concentration_mass_meas = data_Zhao_depth[col_name][(mask_select&mask_units)]
            tracer_detection_limit = data_Zhao_depth[col_name+'_detlim.'][(mask_select&mask_units)]
            l_min = l_arr[mask_Zhao][-1]
            l_max = l_arr[mask_Zhao][0]
            
            df_tmp = pd.DataFrame({'Year': year_current,
                                   'Month': month_current,
                                   'h3_cell_id': h3_selection.values,
                                   'Modelled': tracer_concentration_mass,
                                   'Measured': tracer_concentration_mass_meas.values,
                                   'Units': 'g m-3',
                                   'var_log10': dict_variance['m_surf'],
                                   'l_min': l_min,
                                   'l_max': l_max,
                                   'ParentEventID': data_Zhao_depth['ParentEventID'][(mask_select&mask_units)].values,
                                   'index_parent':data_Zhao_depth.index[(mask_select & mask_units)].values,
                                   'Depth': data_Zhao_depth['Depth'][(mask_select & mask_units)].values,  
                                   'detection_limit': tracer_detection_limit.values})
            df_results = pd.concat((df_results,df_tmp), ignore_index=True)    
            
        t2 = time()
        if i1 % 40 == 0:
            print('%i/%i, t=%2.2e sec (matmul), t=%2.2e sec (pandas), %f tracer' % (i1,n_steps,(t1-t0),(t2-t1),tracer_matrix_mass.sum()))
            print('removed masses: %f, %f percent' % (mass_removed[0],mass_removed[1]))
            print('mass hit bottom: %f percent' % mass_bottom)
            
        if i1 % plot_per == 0 and plot_map and i1 >= i_start:
            fig = plt.figure(figsize=(25, 9))
            gs = gridspec.GridSpec(n_rows, n_cols, bottom=.05, top=.95, wspace=.1)
            for i4,l_ in enumerate(l_arr[0:n_cols]):
                str_size = 'year %i' % (i1/12)
                draw_map(tracer_matrix_mass[:,i4]/dict_h3_normalize.values,extent,gs,i4,base_title=str_size+', %3.0f-%3.0f m depth',stop_earlier=(len(depth_layers)-n_rows-1) )
            
            if plot_tracer_amount:
                fig.text(0.54,0.9,'%3.3f tracer'%tracer_matrix_mass.sum(),fontdict={'color':'red',
                                                                           'size':13})
            fig.savefig(folder_figure+'%3.3i.png' % (i1+1))
            plt.close('all')
      
        
    # for aggregated values, calculate the mean/std/CV    
    data_beach_mean['modelledValue_mean'] = data_beach_mean['modelled_sum'] / data_beach_mean['modelled_n']
    data_beach_mean['modelledValue_std'] = np.sqrt((data_beach_mean['modelled_sum2'] / data_beach_mean['modelled_n']) - data_beach_mean['modelledValue_mean']**2)
    data_beach_mean['modelledValue_CV'] = data_beach_mean['modelledValue_std']/data_beach_mean['modelledValue_mean']
        
    df_tmp = data_beach_mean.rename(columns={'Year_start':'Year','Month_start':'Month','h3':'h3_cell_id',
                                        'modelledValue_mean':'Modelled','measurementValue':'Measured',
                                        'measurementUnit':'Units'})
    # df_tmp['Units'] = 'n m-1'
    # df_tmp['var_log10'] = dict_variance['n_beach']
    df_tmp.loc[df_tmp['Units']=='n m-1','var_log10'] = dict_variance['n_beach']
    df_tmp.loc[df_tmp['Units']=='g m-1','var_log10'] = dict_variance['m_beach']
    df_tmp['l_min'] = l_arr[mask_beach_cleanup][-1]
    df_tmp['l_max'] = l_arr[mask_beach_cleanup][0]
    df_tmp['h3_cell_id'] += index_add_coast
    df_tmp['index_parent'] = df_tmp.index.values
    df_tmp['Depth'] = 0
    df_tmp['detection_limit'] = 0
    df_tmp =df_tmp[['Year','Month','h3_cell_id','Modelled','Measured','Units','var_log10','l_min','l_max','ParentEventID','index_parent','Depth','detection_limit']]
    
    df_results = pd.concat((df_results,df_tmp), ignore_index=True)
    
    #-------------------------save tracer results-------------------
    
    file_tracer = os.path.join(data_folder,'00_tracer_files/tracer_info_findiff_inversion_nsteps%i_'%n_steps + str_info + str(datetime.now())[0:10] + '_' + str(datetime.now())[11:] + '.npy')
    dict_save = {}
    dict_save['tracer_matrix_mass'] = tracer_matrix_mass
    dict_save['tracer_matrix_num'] = tracer_matrix_num
    dict_save['tracer_matrix_num_f'] = tracer_matrix_num_f
    
    dict_save['df_results'] = df_results
    to_pickle(dict_save,file_tracer)

    if final_run:
        file_ = os.path.join(data_folder,'00_assimilation_files/tracer_info_final_nsteps%i_'%n_steps + str_info  + str(datetime.now())[0:10] + '_' + str(datetime.now())[11:] + '.pickle')
        to_pickle(dict_save,file_)
    
    return df_results


#%%

os.system("taskset -p 0xffffff %d" % os.getpid())

def to_unscaled(p_scaled,mu_params_,sigma_params_):
    p_unscaled = (p_scaled*sigma_params_) + mu_params_
    return p_unscaled
    

def bound_parameters(p_scaled,mu_params_,sigma_params_):
    
    parameters_unscaled = to_unscaled(p_scaled,mu_params_,sigma_params_)
    
    parameters_unscaled[log_parameters] =10**(parameters_unscaled[log_parameters])
    
    capping_necessary = False
    mask_near_lower = (cap_0_1 & (parameters_unscaled <= 0.0) ) | (p_scaled < -9)
    mask_near_upper = cap_0_1 & (parameters_unscaled >= 1.) | (p_scaled > 9)
    
    # print(parameters_unscaled)
    
    if mask_near_lower.sum() > 0 or mask_near_upper.sum() > 0:
        capping_necessary = True
    
        if mask_near_lower.sum() > 0:
            # print(mask_near_lower)
            p_scaled[mask_near_lower] *= .99
            # print('capping %s' % np.array(parameter_names)[mask_near_lower][0])
            # print(p_scaled[mask_near_lower])
        if mask_near_upper.sum() > 0:
            # print(mask_near_lower)
            p_scaled[mask_near_upper] *= .99
            # print('capping %s' % np.array(parameter_names)[mask_near_upper][0])
            # print(p_scaled[mask_near_upper])
    return p_scaled, capping_necessary


def iter_bound_parameters(p_scaled,mu_params_,sigma_params_):
    iter_necessary = True
    while iter_necessary:
        p_scaled,iter_necessary = bound_parameters(p_scaled,mu_params_,sigma_params_)
    return p_scaled

    
    
    
if __name__ == '__main__':
    p = ArgumentParser(description="""Parcels runs to construct global transition matrices""")
    p.add_argument('-n_i', '--n_i', default=4, type=int, help='amount of iterations')
    p.add_argument('-n_p', '--n_p', default=11, type=int, help='amount of parallel runs')

    args = p.parse_args()

    n_iter = int(args.n_i)
    n_parallel = int(args.n_p)
    
    parameters_start = np.zeros(len(mu_parameters))

    df_prior = forward_model(parameters_start)
    
    mask_usable_data = (~np.isnan(df_prior['Measured']).values & (df_prior['Measured'] > 0).values & ~np.isnan(df_prior['Modelled']).values & (df_prior['Modelled'] > 0).values)
    min_val_fill = 1e-15    
    n_data = len(df_prior.loc[mask_usable_data,:]) #n_data
    
    data_measured = np.log10(df_prior.loc[mask_usable_data,'Measured'].values)
    mask_detection_limit = df_prior.loc[mask_usable_data,'detection_limit'].values.astype(bool)
    mask_fis_frac = (df_prior.loc[mask_usable_data,'Units'] == 'n n-1')
    
    variance_data = df_prior.loc[mask_usable_data,'var_log10'].values
    

    
    def run_parallel(param_f,data_f,dict_ESMDA,i1=0,startup=False):
        pool = multiprocessing.Pool(processes=n_parallel)
        output_parallel = pool.map(forward_model, param_f.T)
    
        df_mean = None
        for i in range(n_ensemble):
            # df_member = forward_model(param_f[:,i])
            df_member = output_parallel[i]
            if i == 0:
                df_mean = df_member.copy() # save the best estim. for plotting
            modelled_vals = df_member.loc[mask_usable_data,'Modelled'].values
            mask_0 = (modelled_vals == 0)
            if mask_0.sum() > 0:
                print('Zeros found in modelled values, setting to minval %e' % min_val_fill)
                modelled_vals[mask_0] = min_val_fill
            data_f[:,i] = np.log10(modelled_vals)  
    
            dict_ESMDA['data']['log10_fore_%i_%i' % (i1,i)] = data_f[:,i]
            dict_ESMDA['parameters']['fore_%i_%i' % (i1,i)] = param_f[:,i]
        pool.close()   
        return data_f, dict_ESMDA, df_mean
        
    n_iter = 6
    n_parallel = 11
    
    
    time_ = time()
    time_tot = time()
    
    # parameters_start = np.array([ -0.26500744, -3,  -1.01769122,  -2.78810885,
    #          3,  -3,   2.04645507,   1.87197762,
    #         -1.79971577,  -1,   2.76901717,  3,
    #          1.00191666,   1.2814201 ,  -2.26288741,   0.06911989,
    #          0.41531513,  -0.43232664])
    
    # parameters_start = np.array([-0.85951897, -2.66709247,  2.2389081 , -0.92924097,  1.74702634,
    #        -5, -3.99194948, -0.12170348, -0.52174837, -1.78690392,
    #         5, 5,  0.77339775, -2.92964171,  0.7727559 ,
    #         4.39988365, -0.47139919,  0.68114616])
    
    PLOT = True
    eps = 0.01
    mu = .8
    mu_end = .1
    n_param = len(mu_parameters) #mLength
    n_ensemble = n_param + 1
    m_n = parameters_start.copy()
    
    n_sigma_param = 3 #confidence interval: in how many std is the difference between upper/mid estimates given?
    C_D = np.diag(variance_data)
    C_M = ((1/n_sigma_param)**2)*np.eye(n_param)
    
    dict_ESMDA = {}
    dict_ESMDA['param_names'] = parameter_names 
    dict_ESMDA['param_mu'] = mu_parameters
    dict_ESMDA['param_sigma'] = sigma_parameters
    dict_ESMDA['param_log'] = log_parameters
    dict_ESMDA['data'] = df_prior.loc[mask_usable_data,:].copy()
    dict_ESMDA['data']['log10_meas'] = data_measured
    dict_ESMDA['parameters'] = {} 
        
    f_mu = (mu_end/mu)**(1/n_iter)
    
    data_f = np.zeros([n_data, n_ensemble])
    for i1 in range(n_iter+1):
        
        param_f = np.zeros([n_param, n_ensemble]) #Prior ensemble
        param_f[:,0] = m_n
        
        if i1 < n_iter:
    
            for i2 in range(n_param):
                perturbation = np.zeros(n_param)
                perturbation[i2] += eps
                param_f[:,i2+1] = m_n + perturbation
    
            data_f,dict_ESMDA,csv_res = run_parallel(param_f,data_f,dict_ESMDA) 
    
            data_central = data_f[:,0]
            dg_dm = np.zeros([n_data,n_param])
            for i2 in range(n_param):
                dg_dm[:,i2] = (data_f[:,i2+1] - data_central) / eps #forward difference
    
            gc.collect()

            print('Inverting matrix...')
            t1 = np.linalg.inv(np.dot(dg_dm.T,np.dot(np.linalg.inv(C_D),dg_dm)) + np.linalg.inv(C_M))
            t2 = np.dot(dg_dm.T,np.dot(np.linalg.inv(C_D),(data_measured - data_f[:,0]))) 
            t3 = np.dot(np.linalg.inv(C_M),(param_f[:,0] - np.zeros(len(param_f[:,0]))))
            print('Inversion done')
    
            m_np1 = m_n - mu*np.dot(t1,(t2+t3))
            m_n = m_np1
            mu = mu*f_mu
            print(m_n)   
            
            m_n = iter_bound_parameters(m_n,mu_parameters,sigma_parameters)
            print(m_n)
        
        else:
            csv_res = forward_model(m_n)
        
    to_pickle(dict_ESMDA, os.path.join(data_folder,'00_assimilation_files/%i_iter_%i_members_'%(n_iter,n_ensemble) + str_info  + str(datetime.now())[0:10] + '_' + str(datetime.now())[11:16] + '.pickle'))
    csv_res.to_csv(os.path.join(data_folder,'00_assimilation_files/analysis_%i_iter_%i_members_'%(n_iter,n_ensemble) + str_info  + str(datetime.now())[0:10] + '_' + str(datetime.now())[11:16] + '.csv'))
    time_tot = print_time(time_tot,'total time',div_by=60)
    
    