#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:08:11 2021
v2: adding 2D modes as well (where cells in the depth go to the surface)
@author: kaandorp
"""

# from parcels import FieldSet, ParticleSet, AdvectionRK4_3D, ErrorCode, ParticleFile, Field, JITParticle, AdvectionRK4
import numpy as np
import os
from glob import glob
# import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from datetime import datetime, timedelta
from h3.unstable import vect
from h3 import geo_to_h3
from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize
# import cartopy.crs as ccrs
from tools import from_pickle, to_pickle, pcolorhex
# from tools import pcolorhex
# from cartopy import feature
# import matplotlib.gridspec as gridspec
# from scipy import special
# from time import time
# import cmocean
from scipy import sparse


def geo_to_h3_depth(lat,lon,depth,res,remove_last=False):
    h3_index = vect.geo_to_h3(lat,lon,res)
    
    if len(depth_layers) > 2:
        
        for i1,(d_upper,d_lower,index_remove) in enumerate(zip(depth_layers[:-1],depth_layers[1:],index_remove_depth)):
            mask = (depth >= d_upper) & (depth < d_lower)
            if i1 == 0: #surface: don't change the h3 index
                pass
            elif i1 == (len(depth_layers) - 2) and remove_last: #below the last depth layer: move particles to index 1 (removal due to sinking too deep)
                h3_index[mask] = 1
                print('%i particles removed below depth %f meters' % (mask.sum(),d_upper))
            else:
                h3_index[mask] -= index_remove #at other depth layers: remove the appropriate depth index
            print('n particles in depth layer (%f-%f): %i' % (d_upper,d_lower,mask.sum()) )

    return h3_index


def get_transition_matrices_2D(filename, map_h3_to_mat=None, map_mat_to_h3=None, dt_days=30, removal_dict=None):
    """
    For 2D transport, it is assumed all particles at depth move at the surface, where they then experience transport
    """
    # transport at the ocean is resolved:
    data_t1 = xr.open_dataset(filename)

    t0 = data_t1['time'][0,0].values
    t1 = pd.Timestamp(t0) + timedelta(days=dt_days)
    
    mask_complete = (data_t1['time'][:,-1] == t1)
    
    i_bottom = np.unique(np.where(data_t1['hit_bottom'] == 1)[0])
    mask_hit_bottom = np.zeros(len(mask_complete),dtype=bool)
    mask_hit_bottom[i_bottom] = True

    mask_oob_arctic = (~mask_complete) & ( (data_t1['lat'][:,0] < -75) | (data_t1['lat'][:,0] > 85) ) & (~mask_hit_bottom)
    
    h3_index_t0 = geo_to_h3_depth(np.double(data_t1['lat'][:,0].values), np.double(data_t1['lon'][:,0].values), data_t1['z'][:,0].values, h3_grid_res)
    
    h3_index_t1 = geo_to_h3_depth(np.double(data_t1['lat'][:,-1].values), np.double(data_t1['lon'][:,-1].values), data_t1['z'][:,-1].values, h3_grid_res)
    
    # set elements which ended up in cells not present in the original data to zero (i.e. removal)
    mask_valid = np.in1d(h3_index_t1,map_mat_to_h3.values) # which indices are present in the original indices
    h3_index_t1[~mask_valid] = np.uint(0)
        
    # add separate index for arctic out-of-bounds particles
    h3_index_t1[mask_oob_arctic] = np.uint(1)
    
    index_0 = map_h3_to_mat[h3_index_t0].values
    index_1 = map_h3_to_mat[h3_index_t1].values
    n_total = len(map_h3_to_mat)
    
    transitions = np.ones(len(index_0))
    W = coo_matrix((transitions, (index_1, index_0)), shape=(n_total, n_total))
    # W[0,0] = 1 #removed particles remain removed
    
    #baseline transition matrix for transport in the ocean
    P_ocean = normalize(W, norm='l1', axis=0)#.tolil()
    
    ##----------this block was used to copy surface transport to deeper layers, left out now since it does not make sense after all
    # for d_ in depth_layers[1:-1]:
    #     h3_index_t0_d = geo_to_h3_depth(np.double(data_t1['lat'][:,0].values), np.double(data_t1['lon'][:,0].values), d_*np.ones(len(data_t1['lon'][:,0].values))+0.01, h3_grid_res)

    #     # not all depth cells exist in the transition matrix (limited bathymetry)
    #     mask_use = np.isin(h3_index_t0_d,map_h3_to_mat.index)

    #     index_0_d = map_h3_to_mat[h3_index_t0_d[mask_use]].values
    #     index_1_d = map_h3_to_mat[h3_index_t1[mask_use]].values
        
    #     transitions = np.ones(len(index_0_d))
    #     W = coo_matrix((transitions, (index_1_d, index_0_d)), shape=(n_total, n_total))

    #     #transition matrix for transport in the ocean, particles at depth go to the surface
    #     P_ocean_d = normalize(W, norm='l1', axis=0)#.tolil()
        
    #     P_ocean = P_ocean + P_ocean_d    
        
    mask_valid = np.in1d(h3_index_coast, map_mat_to_h3.values) # which indices are present in the original indices
    h3_index_coast_water = h3_index_coast[mask_valid]
    
    mask_from_coast = np.in1d(index_0,map_h3_to_mat[h3_index_coast_water]) #check which starting points start @ coastal water
    indices_unique_fromcoast = np.unique(np.array([index_0[mask_from_coast],index_1[mask_from_coast]]), axis=1) #for the coo_matrix we only want the unique indices
    indices_unique_fromocean = np.unique(np.array([index_0[~mask_from_coast],index_1[~mask_from_coast]]), axis=1)
    
    P_coast_water = coo_matrix((np.ones(indices_unique_fromcoast.shape[1]), 
                                (indices_unique_fromcoast[1,:],indices_unique_fromcoast[0,:])), shape = (n_total,n_total))
    P_water_water = coo_matrix((np.ones(indices_unique_fromocean.shape[1]), 
                                (indices_unique_fromocean[1,:],indices_unique_fromocean[0,:])), shape = (n_total,n_total))
    
    
    # non-oceanic transport is only calculated the first time, as this is the same for all 'hard parameters'
    P_coast_land = None 
    P_land_coast = None 
    P_land_land = None
    P_removal = None 
    P_cascade_land = None
    P_cascade_water = None
    vec_p_bottom = np.zeros(len(map_h3_to_mat))
    
    return (P_ocean, P_coast_water, P_water_water, P_coast_land, P_land_coast, P_land_land, 
            P_removal, P_cascade_land, P_cascade_water, vec_p_bottom, map_h3_to_mat, map_mat_to_h3, dt_days)
    
def get_transition_matrices(filename, map_h3_to_mat=None, map_mat_to_h3=None, dt_days=30, removal_dict=None):
    """
    Get the transition matrix for the ocean gridcells based on the parcels simulations
    """
    data_t1 = xr.open_dataset(filename)

    t0 = data_t1['time'][0,0].values
    t1 = pd.Timestamp(t0) + timedelta(days=dt_days)
    
    mask_complete = (data_t1['time'][:,-1] == t1)
    
    i_bottom = np.unique(np.where(data_t1['hit_bottom'] == 1)[0])
    mask_hit_bottom = np.zeros(len(mask_complete),dtype=bool)
    mask_hit_bottom[i_bottom] = True

    mask_oob_arctic = (~mask_complete) & ( (data_t1['lat'][:,0] < -75) | (data_t1['lat'][:,0] > 85) ) & (~mask_hit_bottom)
    
    # mask_unaccounted = (~mask_complete) & (~mask_oob_arctic) & (~mask_hit_bottom)
    
    
    h3_index_t0 = geo_to_h3_depth(np.double(data_t1['lat'][:,0].values), np.double(data_t1['lon'][:,0].values), data_t1['z'][:,0].values, h3_grid_res)
    
    h3_index_t1 = geo_to_h3_depth(np.double(data_t1['lat'][:,-1].values), np.double(data_t1['lon'][:,-1].values), data_t1['z'][:,-1].values, h3_grid_res)
    
    # set elements which ended up in cells not present in the original data to zero (i.e. removal)
    mask_valid = np.in1d(h3_index_t1,h3_index_t0) # which indices are present in the original indices
    h3_index_t1[~mask_valid] = np.uint(0)
        
    # add separate index for arctic out-of-bounds particles
    h3_index_t1[mask_oob_arctic] = np.uint(1)
    
    # separate index for particles hitting the ocean floor
    if removal_dict:
        
        if removal_dict['remove_bottom_hits'] == True:
            print('bottom hits are moved to matrix element 2')
            print('Total: %i/%i particles hit the bottom' % (mask_hit_bottom.sum(),len(mask_hit_bottom)))
            h3_index_t1[mask_hit_bottom] = np.uint(2)
    
        if removal_dict['remove_under_500'] == True:
            print('Particles under 500 are also moved to matrix element 2')
            # no primary productivity at deepest layer, making the permanent sinking in this layer artificially low, 
            # therefore assume for the special permanent fouling scenario that particles starting here sink to the bottom
            mask_start_below_500 = (data_t1['z'][:,0] > depth_layers[-2] )
            h3_index_t1[mask_start_below_500] = np.uint(2)
    
    # dt_days = float((data_t1['time'][0,:][-1] - data_t1['time'][0,:][0]) / np.timedelta64(1,'D'))
    
    # if no map exists from h3 indices to matrix indices, initialize this. This is the first step that needs to be done.
    # all unique h3 indices are determined here, so make sure particles are released in every h3 cell at the start of this analysis
    if map_h3_to_mat is None:
        
        # indices_h3: hex format from the uber h3 grid
        indices_h3,counts = np.unique(h3_index_t0, return_counts=True) #TODO: think about what cells to take in consideration...
        # for now only cells are taken into consideration which are present for which release points are available, which makes sense
        # as otherwise there would be no transport defined for them?
        
        # h3_index_water = indices_h3.copy()
        

        # -----------coastal cells-------------
        
        # dict_coastlines = from_pickle('00_data_files/coastal_properties_h3_res%i_110m.pickle' % h3_grid_res)
        # h3_index_coast = np.array([hex_ for hex_ in dict_coastlines.keys() if 'length_coast' in dict_coastlines[hex_].keys() ],dtype=np.uint64)
        
        mask_valid = np.in1d(h3_index_coast, h3_index_t0) # which indices are present in the original indices
        
        h3_index_coast_water = h3_index_coast[mask_valid]
        
        h3_index_coast_land = h3_index_coast_water + index_add_coast 
        
        indices_h3 = np.append(indices_h3, h3_index_coast_land)
        
        
        #------------deep cells--------------
        if depth_layers[1] == np.inf:
            pass
        else:
            print('Running transition matrix in 3D mode')
        
        
        #---------------a map from index_h3 to index_mat. A zero is added for removal of particles----------------
        n_total = len(indices_h3) + 3
        indices_h3 = np.append(np.uint([0,1,2]),indices_h3)
        # indices_mat: running index from 0 to n_total
        indices_mat = np.arange(0,n_total,dtype=np.uint(64))
        
        map_h3_to_mat = pd.Series(index=indices_h3,data=indices_mat)
        map_mat_to_h3 = pd.Series(index=map_h3_to_mat.values, data=map_h3_to_mat.index)

        # map h3_index to mat_index for storage in transition matrix
        index_0 = map_h3_to_mat[h3_index_t0].values
        index_1 = map_h3_to_mat[h3_index_t1].values        
    
        #########################################################
        # The rest of this part initializes the connectivity matrices between land/ocean
        # these matrices only contain 0/1's, these are later multiplied by the parameters defining the beaching/resus. probabilities
        # the same goes for the fragmentation connectivity matrices
        #----------------add coastal stuff-----------------------

        # from coastal (ocean) cells to any other ocean (coast/open ocean) cell: multiply 
        # with (1-beach) probability
        mask_from_coast = np.in1d(index_0,map_h3_to_mat[h3_index_coast_water]) #check which starting points start @ coastal water
        indices_unique_fromcoast = np.unique(np.array([index_0[mask_from_coast],index_1[mask_from_coast]]), axis=1) #for the coo_matrix we only want the unique indices
        indices_unique_fromocean = np.unique(np.array([index_0[~mask_from_coast],index_1[~mask_from_coast]]), axis=1)
        
        P_coast_water = coo_matrix((np.ones(indices_unique_fromcoast.shape[1]), 
                                    (indices_unique_fromcoast[1,:],indices_unique_fromcoast[0,:])), shape = (n_total,n_total))
        P_water_water = coo_matrix((np.ones(indices_unique_fromocean.shape[1]), 
                                    (indices_unique_fromocean[1,:],indices_unique_fromocean[0,:])), shape = (n_total,n_total))
        
        # from coastal water to land: entries get p_beach
        P_coast_land = coo_matrix((np.ones(len(h3_index_coast_water)), 
                                    (map_h3_to_mat[h3_index_coast_land],map_h3_to_mat[h3_index_coast_water])), shape = (n_total,n_total))
        
        
        # from coastal land to coastal water: add resuspension probabilities
        P_land_coast = coo_matrix((np.ones(len(h3_index_coast_water)), 
                                    (map_h3_to_mat[h3_index_coast_water],map_h3_to_mat[h3_index_coast_land])), shape = (n_total,n_total))
        
        # remaining on land: multiply with probability of not resuspending (1-p_resus)
        P_land_land = coo_matrix((np.ones(len(h3_index_coast_water)), 
                                    (map_h3_to_mat[h3_index_coast_land],map_h3_to_mat[h3_index_coast_land])), shape = (n_total,n_total))
        
        # removed particles remain removed
        P_removal = coo_matrix((np.array([1,1,1]),(np.array([0,1,2]),np.array([0,1,2]))),shape=(n_total,n_total))
        
        # fragmentation cascade matrices
        P_cascade_water = coo_matrix((np.ones(n_total),(indices_mat,indices_mat)), shape=(n_total,n_total))
        
        P_cascade_land = P_land_land.copy()
        
        P_cascade_water -= P_cascade_land
        

    # all unique h3 indices are already determined. Put the particle locations in the existing indices    
    else:
    
        mask_valid = np.in1d(h3_index_t1,map_mat_to_h3.values) # which indices are present in the original indices
        h3_index_t1[~mask_valid] = 0
        # map h3_index to mat_index for storage in transition matrix
        index_0 = map_h3_to_mat[h3_index_t0].values
        index_1 = map_h3_to_mat[h3_index_t1].values
        n_total = len(map_h3_to_mat)
  
        mask_valid = np.in1d(h3_index_coast, map_mat_to_h3.values) # which indices are present in the original indices
        h3_index_coast_water = h3_index_coast[mask_valid]
        
        mask_from_coast = np.in1d(index_0,map_h3_to_mat[h3_index_coast_water]) #check which starting points start @ coastal water
        indices_unique_fromcoast = np.unique(np.array([index_0[mask_from_coast],index_1[mask_from_coast]]), axis=1) #for the coo_matrix we only want the unique indices
        indices_unique_fromocean = np.unique(np.array([index_0[~mask_from_coast],index_1[~mask_from_coast]]), axis=1)
        
        P_coast_water = coo_matrix((np.ones(indices_unique_fromcoast.shape[1]), 
                                    (indices_unique_fromcoast[1,:],indices_unique_fromcoast[0,:])), shape = (n_total,n_total))
        P_water_water = coo_matrix((np.ones(indices_unique_fromocean.shape[1]), 
                                    (indices_unique_fromocean[1,:],indices_unique_fromocean[0,:])), shape = (n_total,n_total))
        
        
        # non-oceanic transport is only calculated the first time, as this is the same for all 'hard parameters'
        P_coast_land = None 
        P_land_coast = None 
        P_land_land = None
        P_removal = None 
        P_cascade_land = None
        P_cascade_water = None
        
    transitions = np.ones(len(index_0))
    
    W = coo_matrix((transitions, (index_1, index_0)), shape=(n_total, n_total))
    # W[0,0] = 1 #removed particles remain removed
    
    #baseline transition matrix for transport in the ocean
    P_ocean = normalize(W, norm='l1', axis=0)#.tolil()
    
    
    # calculate bottom hit probabilities
    index_add_bottom_hit = map_h3_to_mat[h3_index_t0[mask_hit_bottom]]
    index_add_n_start = map_h3_to_mat[h3_index_t0]
    
    vec_n_start = np.zeros(len(map_h3_to_mat))
    vec_n_bottom = np.zeros(len(map_h3_to_mat))
    
    np.add.at(vec_n_start, index_add_n_start, 1)
    np.add.at(vec_n_bottom, index_add_bottom_hit, 1)
    
    vec_p_bottom = vec_n_bottom / vec_n_start
    vec_p_bottom[np.isnan(vec_p_bottom)] = 0
        
    
    
    return (P_ocean, P_coast_water, P_water_water, P_coast_land, P_land_coast, P_land_land, 
            P_removal, P_cascade_land, P_cascade_water, vec_p_bottom, map_h3_to_mat, map_mat_to_h3, dt_days)


def merge_connection_matrices(list_):
    """
    Merges two connection matrices: these are the matrices containing 0/1's, 
    specifying which cells are connected to which
    """
    mat_tot = None
    for i1, mat_ in enumerate(list_):
        if i1 == 0:
            mat_tot = mat_
        else:
            mat_tot = mat_tot + mat_
    # mat_tot = mat1 + mat2
    mat_tot[mat_tot>1] = 1
    return mat_tot


def glob_defined_date_old(pattern,years=np.arange(2004,2005),months=np.arange(1,13),days=np.arange(1,32)):
    
    files = []
    for yr_ in years:
        for mo_ in months:
            for da_ in days:
                files_ = glob(pattern % (yr_,mo_,da_))
                if files_:
                    files.extend(glob(pattern % (yr_,mo_,da_)))
    return files


def glob_defined_date(pattern,years=np.arange(2004,2005),months=np.arange(1,13),days=np.arange(1,32)):
    
    files = []
    for yr_,mo_,da_ in zip(years,months,days):
        files_ = glob(pattern % (yr_,mo_,da_))
        if files_:
            files.extend(glob(pattern % (yr_,mo_,da_)))
    return files


def add_transition_matrices(string_,map_h3_to_mat=None, map_mat_to_h3=None, dt_days=None,
                            years=np.arange(2004,2005),months=np.arange(1,13),days=np.arange(1,32),mode='3D'):
    
    files = glob_defined_date(string_,years=years,months=months,days=days)
    # files = glob(string_)
    print('merging %s' % files)
    
    P_ocean = None
    P_coast_water = None
    P_water_water = None
    P_coast_land  = None
    P_land_coast = None
    P_land_land = None 
    P_removal = None
    P_cascade_land = None
    P_cascade_water = None
    P_bottom_hit = None
    # map_h3_to_mat = None
    # map_mat_to_h3 = None
    # dt_days = None

    removal_dict = None
    if 'permanent' in mode:
        removal_dict = {}
        removal_dict['remove_bottom_hits'] = True
    
        if 'removal500' in mode:
            removal_dict['remove_under_500'] = True
        else:
            removal_dict['remove_under_500'] = False
            
    
    if '2D' in mode:
        fn_get_transition_mat = lambda a,b,c,d,e: get_transition_matrices_2D(a,b,c,d,e)
    else:
        fn_get_transition_mat = lambda a,b,c,d,e: get_transition_matrices(a,b,c,d,e)
    
    for i1,file_ in enumerate(files):
        
        if i1 == 0:
            (P_ocean, P_coast_water, P_water_water, P_coast_land, P_land_coast, P_land_land, 
             P_removal, P_cascade_land, P_cascade_water, P_bottom_hit, map_h3_to_mat, map_mat_to_h3, dt_days) = fn_get_transition_mat(file_, map_h3_to_mat, map_mat_to_h3, 
                                                                                                                                      dt_days,removal_dict)
        else:
            (P_ocean_, P_coast_water_, P_water_water_, _, _, _, 
             _, _, _, P_bottom_hit_, _, _, dt_days_) = fn_get_transition_mat(file_, map_h3_to_mat, map_mat_to_h3, 
                                                                             dt_days,removal_dict)
            
            assert (dt_days_ == dt_days), 'different timesteps: {} {}, file {}'.format(dt_days,dt_days_,file_)
            
            P_ocean = P_ocean + P_ocean_
            P_bottom_hit = P_bottom_hit + P_bottom_hit_
            
            P_coast_water = merge_connection_matrices([P_coast_water,P_coast_water_])
            P_water_water = merge_connection_matrices([P_water_water,P_water_water_])
            
    P_ocean = P_ocean / len(files)
    #TODO: normalize P_ocean once again w.r.t. column sums (sometimes they do not add up to 1 after summing?)
    P_bottom_hit = P_bottom_hit / len(files)
       
    return (P_ocean, P_coast_water, P_water_water, P_coast_land, P_land_coast, P_land_land, 
            P_removal, P_cascade_land, P_cascade_water, P_bottom_hit, map_h3_to_mat, map_mat_to_h3, dt_days)
 
if os.environ['USER'] == 'kaandorp': # desktop
    folder_runfiles = '/Volumes/externe_SSD/kaandorp/Data/Global_project_data/'
    home_folder = '.'
elif os.environ['USER'] == 'kaand004': #gemini
    folder_runfiles = '/storage/shared/oceanparcels/output_data/data_Mikael/00_run_files/'
    home_folder = '/nethome/kaand004/Git_repositories/Global_Analysis_Mikael'    
    

#3: ~40,000 cells, 60km
#2: ~6,000 cells, 160km
h3_grid_res = 3
coastline_res = '50m'
depth_layers = np.array([0,5,50,500,np.inf])

h3_index_base = np.uint64(5e17)
index_add_coast = np.uint64(1e17)     #indices for coastal cells are the ocean surface cells + index_add_coast
index_remove_depth = np.array([int(0),int(1e17),int(2e17),int(3e17),int(4e17)])

dict_coastlines = from_pickle(os.path.join(home_folder,'00_data_files/coastal_properties_h3_res%i_%s.pickle' % (h3_grid_res,coastline_res)))
h3_index_coast = np.array([hex_ for hex_ in dict_coastlines.keys() if 'length_coast' in dict_coastlines[hex_].keys() ],dtype=np.uint64)
    
#%% Calculate transition matrices

def disconnect_cell(P,hex_disconnect):
    P = P.tolil()
    for i1 in range(len(depth_layers)-2):
        hex_ = np.int64(hex_disconnect - index_remove_depth[i1])
        
        if hex_ in map_h3_to_mat.index: #doesnt need to be true for shallow bathymetry
        
            P[map_h3_to_mat[hex_],:] = 0
            P[:,map_h3_to_mat[hex_]] = 0
        
    hex_ = np.int64(hex_disconnect) + np.int64(index_add_coast)
    P[map_h3_to_mat[hex_],:] = 0
    P[:,map_h3_to_mat[hex_]] = 0    
    
    return P.tocsc()
    

hex_med = int(geo_to_h3(35.7,-5.6,h3_grid_res),16)
hex_pan = int(geo_to_h3(9.0,-79.6,h3_grid_res),16)

# months = [np.arange(1,12)]
# days = [np.arange(1,32)]
years = [2019]
merge_months = False #if merge months, all data is merged into 1 transition matrix


str_read = {}
str_read['fouled'] = os.path.join(folder_runfiles,'output_MOi_h4_uniform_3D_2301043_%4.4i-%2.2i-%2.2iT12:00:00_30days_MixingTrue_FoulingTrue_nbTrue_StokesFalse_l{:.6f}_*.nc')
str_read['nonfouled'] = os.path.join(folder_runfiles,'output_MOi_h4_uniform_3D_2301043_%4.4i-%2.2i-%2.2iT12:00:00_30days_MixingTrue_FoulingFalse_nbFalse_StokesFalse_l{:.6f}_*.nc')
str_read['neutral'] = os.path.join(folder_runfiles,'output_MOi_h4_uniform_3D_2301043_%4.4i-%2.2i-%2.2iT12:00:00_30days_MixingTrue_FoulingFalse_nbFalse_StokesFalse_l{:.6f}_*.nc')
str_read['permanent_fouled'] = os.path.join(folder_runfiles,'output_MOi_h4_uniform_3D_2301043_%4.4i-%2.2i-%2.2iT12:00:00_30days_MixingTrue_FoulingFalse_nbTrue_PermFoulingTrue_StokesFalse_windage0.00_l{:.6f}_*.nc')
str_read['permanent_fouled_removal500'] = os.path.join(folder_runfiles,'output_MOi_h4_uniform_3D_2301043_%4.4i-%2.2i-%2.2iT12:00:00_30days_MixingTrue_FoulingFalse_nbTrue_PermFoulingTrue_StokesFalse_windage0.00_l{:.6f}_*.nc')
str_read['2D_baseline'] = os.path.join(folder_runfiles,'output_MOi_h4_uniform_2D_205599_%4.4i-%2.2i-%2.2iT12:00:00_30days_MixingFalse_FoulingFalse_nbFalse_StokesFalse_windage0.00_l{:.6f}_*.nc')
str_read['2D_Stokes'] = os.path.join(folder_runfiles,'output_MOi_h4_uniform_2D_205599_%4.4i-%2.2i-%2.2iT12:00:00_30days_MixingFalse_FoulingFalse_nbFalse_StokesTrue_windage0.00_l{:.6f}_*.nc')
str_read['2D_Stokes_wind'] = os.path.join(folder_runfiles,'output_MOi_h4_uniform_2D_205599_%4.4i-%2.2i-%2.2iT12:00:00_30days_MixingFalse_FoulingFalse_nbFalse_StokesTrue_windage0.02_l{:.6f}_*.nc')
str_read['2D_wind'] = os.path.join(folder_runfiles,'output_MOi_h4_uniform_2D_205599_%4.4i-%2.2i-%2.2iT12:00:00_30days_MixingFalse_FoulingFalse_nbFalse_StokesFalse_windage0.02_l{:.6f}_*.nc')

particles_l = {}
particles_l['fouled'] = [0.0001,0.0004,0.0016,0.0064,0.0256,0.1024]
particles_l['nonfouled'] = [0.0001,0.0004,0.0016,0.0064,0.0256,0.1024]
particles_l['neutral'] = [0.00001]
particles_l['permanent_fouled'] = [0.0001,0.0004,0.0016,0.0064,0.0256,0.1024]
particles_l['permanent_fouled_removal500'] = [0.0001,0.0004,0.0016,0.0064,0.0256,0.1024]
particles_l['2D_baseline'] = [1.0]
particles_l['2D_Stokes'] = [1.0]
particles_l['2D_Stokes_wind'] = [1.0]
particles_l['2D_wind'] = [1.0]




if merge_months: 
    months_s = [[]]
    days_s = [[]]
    years_s = [[]]
    for year_ in years:
        dates = pd.date_range('%i-01-01-12' % year_,periods=12,freq='30d')
        for i1,date_ in enumerate(dates):
            months_s[0].append(date_.month)
            days_s[0].append(date_.day)
            years_s[0].append(date_.year)      
else: #keep transition matrix monthly. In this case, create lists in which the starting month and day for each year to be merged are given
    months_s = [[] for i in range(12)]
    days_s = [[] for i in range(12)]
    years_s = [[] for i in range(12)]
    for year_ in years:
        dates = pd.date_range('%i-01-01-12' % year_,periods=12,freq='30d')
        for i1,date_ in enumerate(dates):
            months_s[i1].append(date_.month)
            days_s[i1].append(date_.day)
            years_s[i1].append(date_.year)  

# particles_l = [0.0001,0.0004]#,0.0016,0.0064,0.0256,0.1020]
# months_s = [1]
# years = [2019]

filename_dict_transition_ = '00_transition_files/transitions_h%i_%s_monthly_year_'%(h3_grid_res,coastline_res)+ '-'.join('%2.2i'%s_ for s_ in years) + '_merge_months_%s'%merge_months + '.pickle'
filename_dict_transition = os.path.join(home_folder,filename_dict_transition_)

map_h3_to_mat = None
map_mat_to_h3 = None
dt_days = 30

P_dict = {}
P_ocean = {}
P_bottom_hit = {}
list_P_water_water = []
list_P_coast_water = []

for i_month in range(len(months_s)):
    P_ocean[i_month] = {}
    P_bottom_hit[i_month] = {}
    
    #make sure to start with a case where all cells have been resolved! (i.e. not one of the surface runs, to initialize the map_h3_to_mat)
    for mode_ in ['fouled','nonfouled','neutral','permanent_fouled','permanent_fouled_removal500','2D_baseline','2D_Stokes','2D_Stokes_wind','2D_wind']: 
        P_ocean[i_month][mode_] = {}
        P_bottom_hit[i_month][mode_] = {}
        
        for l_ in particles_l[mode_]:
            string_ = str_read[mode_].format(l_)
    
    
            if map_h3_to_mat is None:       
                (P_ocean_, P_coast_water, P_water_water, P_dict['coast_land'], P_dict['land_coast'], P_dict['land_land'], 
                  P_dict['removal'], P_dict['cascade_land'], P_dict['cascade_water'], P_bottom_, map_h3_to_mat, map_mat_to_h3, dt_days) = add_transition_matrices(string_, dt_days=dt_days, years=years_s[i_month], months=months_s[i_month], days=days_s[i_month], mode=mode_)
            else:
                (P_ocean_, P_coast_water, P_water_water, _, _, _, 
                  _, _, _, P_bottom_, _, _, dt_days_) = add_transition_matrices(string_, map_h3_to_mat, map_mat_to_h3, dt_days=dt_days, years=years_s[i_month], months=months_s[i_month], days=days_s[i_month], mode=mode_)
            
                assert(dt_days == dt_days_) #check the matrices have the same transition time
    
            # P_ocean[mode_][l_] = P_ocean_
            P_ocean[i_month][mode_][l_] = disconnect_cell(P_ocean_,hex_med)
            P_ocean[i_month][mode_][l_] = disconnect_cell(P_ocean_,hex_pan)
            P_bottom_hit[i_month][mode_][l_] = P_bottom_
            
            list_P_water_water.append(P_water_water)
            list_P_coast_water.append(P_coast_water)
    

P_coast_water = merge_connection_matrices(list_P_coast_water) 
P_water_water = merge_connection_matrices(list_P_water_water)

P_dict['coast_water'] = P_coast_water
P_dict['water_water'] = P_water_water
P_dict['map_h3_to_mat'] = map_h3_to_mat
P_dict['map_mat_to_h3'] = map_mat_to_h3
P_dict['ocean'] = P_ocean
P_dict['bottom_hit'] = P_bottom_hit
P_dict['depth_layers'] = depth_layers

#%% primary productivity calculations
n_total = len(map_h3_to_mat)
df_pp = from_pickle(os.path.join(home_folder,'00_data_files/mean_pp_h%i_2019.pickle' % h3_grid_res)).astype(float)
df_pp = df_pp / np.quantile(df_pp.values[~np.isnan(np.array(df_pp.values,dtype=float))],.99) #normalize for outliers
df_pp[df_pp>1] = 1.
df_pp[df_pp<0] = 0.

i_r_o,i_c_o,_ = sparse.find(P_dict['water_water']+P_dict['coast_water'])
cols_ocean = np.unique(i_c_o)

i_r_ww,i_c_ww,_ = sparse.find(P_dict['water_water'])
i_r_cw,i_c_cw,_ = sparse.find(P_dict['coast_water'])
i_r_cl,i_c_cl,_ = sparse.find(P_dict['coast_land'])

pp_water_water = {}
pp_coast_water = {}
pp_coast_land = {}
P_ocean_removed = {}

for i1 in range(12): # for every month

    pp_water_water[i1] = coo_matrix( (df_pp.loc[map_mat_to_h3[i_c_ww],i1].values,(i_r_ww,i_c_ww)), shape=(n_total,n_total))
    pp_coast_water[i1] = coo_matrix( (df_pp.loc[map_mat_to_h3[i_c_cw],i1].values,(i_r_cw,i_c_cw)), shape=(n_total,n_total))
    pp_coast_land[i1] = coo_matrix( (df_pp.loc[map_mat_to_h3[i_c_cl],i1].values,(i_r_cl,i_c_cl)), shape=(n_total,n_total))

    P_ocean_removed[i1] = coo_matrix( (df_pp.loc[map_mat_to_h3[cols_ocean],i1].values,(1*np.ones(len(cols_ocean)),cols_ocean)), 
                                     shape=(n_total,n_total))
    
P_dict['pp_water_water'] = pp_water_water
P_dict['pp_coast_water'] = pp_coast_water
P_dict['pp_coast_land'] = pp_coast_land
P_dict['P_ocean_removed'] = P_ocean_removed


to_pickle(P_dict,filename_dict_transition)

