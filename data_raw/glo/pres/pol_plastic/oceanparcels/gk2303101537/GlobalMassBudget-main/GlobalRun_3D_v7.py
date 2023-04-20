#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 17:14:05 2021
Script to run the parcels simulation from a given start date to create a set of transition matrices
v2: mixing and rise velocities added
v3: biofouling added
v4: with chunking
v5: add grazing
v6: MOi data
v7: option for Stokes/wind forcing
@author: kaandorp
"""

from parcels import FieldSet, ParticleSet, AdvectionRK4_3D, ErrorCode, ParticleFile, Field, \
    JITParticle, AdvectionRK4, ParcelsRandom
# from parcels.application_kernels.TEOSseawaterdensity import 
from parcels.tools.converters import Geographic, GeographicPolar 
import numpy as np
import os
from glob import glob
# import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from datetime import datetime, timedelta
# import math
from argparse import ArgumentParser
from kernels_v3 import periodicBC, delete_particle, delete_particle_interp, markov_0_mixing, \
    PolyTEOS10_bsq, plastic_particle, MOi_biofouling, MOi_permanent_biofouling, initialize_neutral_bouyancy, remove_at_bounds,\
        settling_velocity, delete_at_depth, Stokes_drift, windage_drift, unbeaching, neutral_buoyancy
from utils import getclosest_ij, select_files

import dask
dask.config.set({"array.slicing.split_large_chunks": False})

import warnings
warnings.simplefilter("ignore", category=xr.SerializationWarning)

if __name__=="__main__":
    p = ArgumentParser(description="""Parcels runs to construct global transition matrices""")
    p.add_argument('-mode', '--mode', default='3D', help='2D (surface) or 3D mode')
    p.add_argument('-mixing', '--mixing', default=1, type=int, help='Use the KPP mixing model')
    p.add_argument('-biofouling', '--biofouling', default=0, type=int, help='Use the Lobelle (2021) biofouling model')
    p.add_argument('-permanent_fouling', '--permanent_fouling', default=0, type=int, help='Use the Lobelle (2021) biofouling model without the respiration')
    p.add_argument('-fouled_nb', '--fouled_nb', default=0, type=int, help='Initialize fouled particles to be neutrally bouyant')
    p.add_argument('-stokes', '--stokes', default=0, type=int, help='Include Stokes drift')
    p.add_argument('-windage', '--windage', default=0.0, type=float, help='Fraction of windage to add')
    p.add_argument('-nb', '--nb', default=0, type=int, help='Neutrally buoyant particles')
    
    p.add_argument('-particle_l', '--particle_l', default=0.001, type=float, help='Particle size in meters')
    p.add_argument('-rho_bf', '--rho_bf', default=1388., type=float, help='Biofilm density')
    p.add_argument('-rho_pl', '--rho_pl', default=1010., type=float, help='Plastic density')

    p.add_argument('-date_start', '--date_start', default='2019-01-01-12', type=str, help='Starting date of the transition matrix')    
    p.add_argument('-index_transition', '--index_transition', default=0, type=int, help='Which transition matrix to simulate')
    p.add_argument('-dt_transition', '--dt_transition', default=30, type=float, help='Amount of days for each transition matrix')
    p.add_argument('-dt_write', '--dt_write', default=30, type=float, help='Write particle data every n days')

    p.add_argument('-filename_release', '--filename_release', default='h3_uniform_4212', type=str, help='file with release positions')
    p.add_argument('-test_run', '--test_run', default=0, type=int, help='Try limited domain test run')
    
    p.add_argument('-use_chunking', '--use_chunking', default=0, type=int, help='Use chunking (1) or not (0)' )
    p.add_argument('-chunks_latlon', '--chunks_latlon', default=128, type=int, help='Chunking size in the horizontal')
    p.add_argument('-chunks_d', '--chunks_d', default=25, type=int, help='Chunking size in the vertical')
    
    args = p.parse_args()
    mode = args.mode
    do_mixing = bool(args.mixing)
    do_biofouling = bool(args.biofouling)
    do_permanent_fouling = bool(args.permanent_fouling)
    fouled_nb = bool(args.fouled_nb)
    do_Stokes = bool(args.stokes)
    do_nb = bool(args.nb)
    windage = float(args.windage)
    if windage < 1e-5:
        windage = 0
        print('No windage added...')
        
    particle_l = args.particle_l
    rho_bf = args.rho_bf
    rho_pl = args.rho_pl
    
    date_start = pd.Timestamp(args.date_start)
    index_transition = args.index_transition
    dt_transition = args.dt_transition
    dt_write = args.dt_write   
  
    filename_release = args.filename_release
    test_run = bool(args.test_run)
    
    do_chunking = bool(args.use_chunking)
    chunks_latlon = args.chunks_latlon
    chunks_d = args.chunks_d
    chunks_t = 1
    
    verbose_delete = 0
    
    date_current = date_start + index_transition*timedelta(days = dt_transition)

    dt_mins = 20
    
    if mode == '2D':
        assert(do_mixing==False)
        assert(do_biofouling==False)
        assert(do_nb==False)
        assert(do_permanent_fouling==False)
        
    print('------------Settings overview------------')
    print('Mixing: %s' % do_mixing)
    print('Fouling (oscill.): %s' % do_biofouling)
    print('Fouling (perm.): %s' % do_permanent_fouling)
    print('Initialized fouling: %s' % fouled_nb)
    print('Particle length: %f' % particle_l)
    print('Density plastic: %f' % rho_pl)
    print('Density biofilm: %f' % rho_bf)
    print('Date: %s' % date_current)
    print('Dt transition: %f' % dt_transition)
    print('Mode: %s' % mode)
    print('Stokes drift: %s' % do_Stokes)
    print('windage factor: %f' % windage)
    print('Release file: %s' % filename_release)
    print('Test run: %s' % test_run)
    print('Use chunking: %i' % do_chunking)
    print('Chunking lat/lon: %i' % chunks_latlon)
    print('Chunking depth: %i' % chunks_d)

    print('-----------------------------------------')
    
    
    folder_run = '00_run_files'
    
    if os.environ['USER'] == 'kaandorp': # desktop
        dir_home = '/Users/kaandorp/Git_repositories/Global_Analysis_Mikael/'
        dir_input_data = '/Volumes/externe_SSD/kaandorp/Data/'
        dir_input_data2 = '/Volumes/externe_SSD/kaandorp/Data/'
        dir_output_data = '/Users/kaandorp/Git_repositories/Global_Analysis_Mikael/'
    #     dir_write = '/Users/kaandorp/Git_repositories/Global_Analysis_Mikael/'
    elif os.environ['USER'] == 'kaand004': #gemini
        dir_home = '/storage/home/kaand004/Git_repositories/Global_Analysis_Mikael/'
        dir_input_data = '/storage/shared/oceanparcels/input_data/'
        dir_input_data2 = '/storage/shared/oceanparcels/output_data/data_Mikael/'        
        dir_output_data = '/storage/shared/oceanparcels/output_data/data_Mikael/'
    
    elif os.environ['USER'] == 'mikaelk': #snellius
        dir_home = '/home/mikaelk/Git_repositories/Global_Analysis_Mikael/'
        dir_input_data = '/projects/0/topios/hydrodynamic_data/'
        dir_input_data2 = '/projects/0/topios/hydrodynamic_data/'
        dir_output_data = '/scratch-shared/mikaelk/output_data/'
        # dir_output_data = '/home/mikaelk/Git_repositories/Global_Analysis_Mikael/' # write to the back-upped home folder (200GB)
                
    dir_output_files = os.path.join(dir_output_data,folder_run)
    if not os.path.exists(dir_output_files):
        print('Creating dir %s' % dir_output_files)
        os.makedirs(dir_output_files)
    
    file_release = os.path.join(dir_output_data,'00_release_files/%s.nc' % filename_release)
    data_release = xr.open_dataset(file_release)
    n_particles = len(data_release['lon'])
    
    dirread = os.path.join(dir_input_data, 'MOi/psy4v3r1/')
    dirread_bgc = os.path.join(dir_input_data, 'MOi/biomer4v2r1/')
    dirread_mesh_12th = os.path.join(dir_input_data, 'MOi/domain_ORCA0083-N006/')
    dirread_mesh_4th = os.path.join(dir_input_data, 'MOi/domain_ORCA025-N006/')
    dirread_Stokes = os.path.join(dir_input_data2, 'ERA5/waves/')
    dirread_wind = os.path.join(dir_input_data2, 'ERA5/wind/')

    # get files for 2 years of data (max transition matrix duration: 1 year)
    ufiles = select_files(dirread,'psy4v3r1-daily_U_%4i*.nc',date_current,dt_transition,dt_margin=3)
    vfiles = select_files(dirread,'psy4v3r1-daily_V_%4i*.nc',date_current,dt_transition,dt_margin=3)
    wfiles = select_files(dirread,'psy4v3r1-daily_W_%4i*.nc',date_current,dt_transition,dt_margin=3)   
    tfiles = select_files(dirread,'psy4v3r1-daily_T_%4i*.nc',date_current,dt_transition,dt_margin=3)    
    sfiles = select_files(dirread,'psy4v3r1-daily_S_%4i*.nc',date_current,dt_transition,dt_margin=3)   
    kzfiles = select_files(dirread,'psy4v3r1-daily_KZ_%4i*.nc',date_current,dt_transition,dt_margin=3)   
    ppfiles = select_files(dirread_bgc,'biomer4v2r1-weekly_nppv_%4i*.nc',date_current,dt_transition,dt_margin=8)
    phy1files = select_files(dirread_bgc,'biomer4v2r1-weekly_phy_%4i*.nc',date_current,dt_transition,dt_margin=8)
    phy2files = select_files(dirread_bgc,'biomer4v2r1-weekly_phy2_%4i*.nc',date_current,dt_transition,dt_margin=8)
    wavesfiles = select_files(dirread_Stokes,'ERA5_global_waves_monthly_converted_%4i*.nc',date_current,dt_transition,dt_margin=32)
    windfiles = select_files(dirread_wind,'ERA5_global_wind_monthly_converted_%4i*.nc',date_current,dt_transition,dt_margin=32)

    assert(len(ufiles) == len(kzfiles))
    assert(len(ppfiles) == len(phy2files))
    
    mesh_mask = dirread_mesh_12th+'coordinates.nc'      
    mesh_mask_4th = dirread_mesh_4th+'mesh_hgr_PSY4V3_deg.nc'     


    if test_run:
        test_data = xr.open_dataset(ufiles[0])
        latvals = test_data['nav_lat'].values
        lonvals = test_data['nav_lon'].values
        minlat = 40
        maxlat = 50
        minlon = -20
        maxlon = 0
        iy_min, ix_min = getclosest_ij(latvals, lonvals, minlat, minlon)
        iy_max, ix_max = getclosest_ij(latvals, lonvals, maxlat, maxlon)
        iy_min -= 1
        ix_min -= 1
        iy_max += 1
        ix_max += 1
        indices = {'lat': range(iy_min, iy_max), 'lon': range(ix_min, ix_max), 'depth': range(0,33)} #depth : range(0,2000)
        
        mask_release = (data_release['lon'] > minlon) & (data_release['lon'] < maxlon) & (data_release['lat'] > minlat) & (data_release['lat'] < maxlat)
        
        data_release = data_release.where(mask_release,drop=True)
        n_particles = len(data_release['lon'])
    elif mode=='2D':
        indices = {'depth': range(0,1)}
    else:
        indices = {'depth': range(0,46)} # almost all ocean is covered, up to 4000m
    
    
    filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles}, #'depth': wfiles,
                 'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
                 'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles},
                'cons_temperature': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': tfiles},
                'abs_salinity': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': sfiles}}
    
    variables = {'U': 'vozocrtx',
                 'V': 'vomecrty',
                 'W': 'vovecrtz',
                'cons_temperature': 'votemper',
                'abs_salinity': 'vosaline'}
    
    dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}, #time_centered
                  'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                  'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                 'cons_temperature': {'lon': 'glamf', 'lat': 'gphif','depth': 'depthw', 'time': 'time_counter'},
                 'abs_salinity': {'lon': 'glamf', 'lat': 'gphif','depth': 'depthw', 'time': 'time_counter'}}
    
    
    if do_mixing:
        filenames['mixing_kz'] = {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': kzfiles}

        variables['mixing_kz'] = 'votkeavt'

        dimensions['mixing_kz'] = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}

    if do_biofouling or do_permanent_fouling:
        filenames_bio = {'pp_phyto': {'lon': mesh_mask_4th, 'lat': mesh_mask_4th, 'depth': wfiles[0], 'data': ppfiles},
                        'bio_nanophy': {'lon': mesh_mask_4th, 'lat': mesh_mask_4th, 'depth': wfiles[0], 'data': phy1files},
                        'bio_diatom': {'lon': mesh_mask_4th, 'lat': mesh_mask_4th, 'depth': wfiles[0], 'data': phy2files}}

        variables_bio = {'pp_phyto': 'nppv',
                        'bio_nanophy': 'phy',
                        'bio_diatom': 'phy2'}

        dimensions_bio = {'pp_phyto': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                         'bio_nanophy': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                         'bio_diatom': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}}

    if do_Stokes:
        filenames_S = {'Stokes_U': wavesfiles, #Cannot be U for codegenerator!!
                       'Stokes_V': wavesfiles,
                      'wave_Tp': wavesfiles}

        variables_S = {'Stokes_U': 'ust',
                       'Stokes_V': 'vst',
                      'wave_Tp': 'pp1d'}

        dimensions_S = {'lat': 'lat',
                        'lon': 'lon',
                        'time': 'time'} 

    if windage > 0:
        filenames_wind = {'wind_U': windfiles, #Cannot be U for codegenerator!!
                       'wind_V': windfiles}

        variables_wind = {'wind_U': 'u10',
                       'wind_V': 'v10'}

        dimensions_wind = {'lat': 'lat',
                        'lon': 'lon',
                        'time': 'time'} 

    if do_Stokes or windage > 0:
        unbeachfiles = os.path.join(dir_output_data, '00_data_files/land_current_NEMO0083.nc')
        filenames_unbeach = {'unbeach_U': unbeachfiles, 
                           'unbeach_V': unbeachfiles}

        variables_unbeach = {'unbeach_U': 'land_current_u',
                           'unbeach_V': 'land_current_v'}   

        dimensions_unbeach = {'lat': 'lat',
                            'lon': 'lon'}

        
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices=indices, allow_time_extrapolation=False)
    fieldset.add_constant('G', 9.81) #gravitational constant

    
    if do_mixing:
        fieldset.add_constant('z_start',0.5)
        fieldset.add_constant('rho_bf', rho_bf)                   # density of biofilm [g m-3]
        fieldset.add_constant('V_a', 2.0E-16)                    # Volume of 1 algal cell [m3]  

        variables = ('bathymetry', 'Bathymetry')
        dimensions = {'lon': 'nav_lon', 'lat': 'nav_lat'}
        bathymetry_field = Field.from_netcdf(os.path.join(dirread_mesh_12th,'bathymetry_ORCA12_V3.3.nc'), variables, dimensions)
        fieldset.add_field(bathymetry_field) 
    
    if do_biofouling or do_permanent_fouling:
        
        if not do_mixing:
            raise RuntimeError('Mixing must be enabled for biofouling model to work')
      
        fieldset.add_constant('collision_eff', 1.)
        fieldset.add_constant('K', 1.0306E-13 / (86400. ** 2.))  # Boltzmann constant [m2 kg d-2 K-1] now [s-2] (=1.3804E-23)
        fieldset.add_constant('R20', 0.1 / 86400.)               # respiration rate, now [s-1]

        fieldset.add_constant('Q10', 2.13)                         # temperature coefficient respiration [-]
        fieldset.add_constant('Gamma', 1.728E5 / 86400.)         # shear [d-1], now [s-1]
        fieldset.add_constant('Wt_C', 12.)                    # atomic weight of nitrogen

        bio_fieldset = FieldSet.from_nemo(filenames_bio,variables_bio,dimensions_bio)

        fieldset.add_field(bio_fieldset.pp_phyto) #phytoplankton primary productivity 
        fieldset.add_field(bio_fieldset.bio_nanophy) #nanopyhtoplankton concentration [mmol C m-3]
        fieldset.add_field(bio_fieldset.bio_diatom) #diatom concentration [mmol C m-3]   

        
    if do_Stokes:
        fieldset_Stokes = FieldSet.from_netcdf(filenames_S, variables_S, dimensions_S, mesh='spherical')
        fieldset_Stokes.Stokes_U.units = GeographicPolar()
        fieldset_Stokes.Stokes_V.units = Geographic()
        fieldset_Stokes.add_periodic_halo(zonal=True)
        
        fieldset.add_field(fieldset_Stokes.Stokes_U)
        fieldset.add_field(fieldset_Stokes.Stokes_V)
        fieldset.add_field(fieldset_Stokes.wave_Tp)

    if windage > 0:
        fieldset_wind = FieldSet.from_netcdf(filenames_wind, variables_wind, dimensions_wind, mesh='spherical')
        fieldset_wind.wind_U.units = GeographicPolar()
        fieldset_wind.wind_V.units = Geographic()
        fieldset_wind.wind_U.set_scaling_factor(windage)
        fieldset_wind.wind_V.set_scaling_factor(windage)
        
        fieldset_wind.add_periodic_halo(zonal=True)
        
        fieldset.add_field(fieldset_wind.wind_U)
        fieldset.add_field(fieldset_wind.wind_V)

    if do_Stokes or windage > 0:
        fieldset_unbeach = FieldSet.from_netcdf(filenames_unbeach, variables_unbeach, dimensions_unbeach, mesh='spherical')
        fieldset_unbeach.unbeach_U.units = GeographicPolar()
        fieldset_unbeach.unbeach_V.units = Geographic()

        fieldset.add_field(fieldset_unbeach.unbeach_U)
        fieldset.add_field(fieldset_unbeach.unbeach_V)
    
    if test_run:
        print('Running test run...')
        fieldset.add_constant('lat_min', minlat)
        fieldset.add_constant('lat_max', maxlat)
        fieldset.add_constant('lon_min', minlon)
        fieldset.add_constant('lon_max', maxlon)

        
    fieldset.add_constant('verbose_delete',verbose_delete)
    
    print('-------------------fieldset created-------------------')
    #%% 

    filename_out = 'output_MOi_' + filename_release + '_%s_%idays_Mixing%s_Fouling%s_nb%s_' % (str(date_current.to_datetime64())[0:19],dt_transition,do_mixing,do_biofouling,fouled_nb)
    if do_permanent_fouling:
        filename_out += 'PermFoulingTrue_'
    filename_out += 'Stokes%s_windage%2.2f_l%f_test%s_chunking%i' % (do_Stokes,windage,particle_l,test_run,do_chunking)

    pset = ParticleSet.from_list(fieldset, plastic_particle, 
                                 lon=data_release['lon'].values,
                                 lat=data_release['lat'].values,
                                 time=date_current.to_datetime64(),
                                 depth=data_release['z'].values,
                                 rho_pl=rho_pl*np.ones(n_particles),
                                 l_pl=particle_l*np.ones(n_particles))
    
    
    #------------------------------ Initialization of variables if necessary ------------------------------
    kernels_init = None
    if fouled_nb: #use biofouling, and initialize a biofilm making the particles neutrally buoyant
        print('-------------------Initializing fouling + settling velocity...-------------------')
        if test_run:
            kernels_init = pset.Kernel(remove_at_bounds)+pset.Kernel(PolyTEOS10_bsq)+pset.Kernel(initialize_neutral_bouyancy)+pset.Kernel(settling_velocity)
        else:
            kernels_init = pset.Kernel(PolyTEOS10_bsq)+pset.Kernel(initialize_neutral_bouyancy)+pset.Kernel(settling_velocity)
                   
    elif do_mixing:
        print('-------------------Initializing settling velocity...-------------------')
        if test_run:
            kernels_init = pset.Kernel(remove_at_bounds)+pset.Kernel(PolyTEOS10_bsq)+pset.Kernel(settling_velocity)
        else:
            kernels_init = pset.Kernel(PolyTEOS10_bsq)+pset.Kernel(settling_velocity)
    elif do_nb:
        kernels_init = pset.Kernel(neutral_buoyancy)
        
    else:
        print('-------------------No initialization...-------------------')
    
    if kernels_init:
        pset.execute(kernels_init, runtime=timedelta(days=0), dt=timedelta(minutes=0), verbose_progress=True, 
             recovery={ErrorCode.ErrorOutOfBounds: delete_particle, ErrorCode.ErrorInterpolation: delete_particle_interp})
        print('-------------------Init finished-------------------')

    #------------------------------ Main kernels ------------------------------
    pfile= ParticleFile(os.path.join(dir_output_files,'%s.nc' % filename_out), pset, outputdt=timedelta(days=dt_write))
    
    kernels = pset.Kernel(PolyTEOS10_bsq) + pset.Kernel(delete_at_depth)
    
    if mode =='3D':
        kernels += AdvectionRK4_3D
    else:
        kernels += AdvectionRK4
  
    if test_run:
        kernels += remove_at_bounds   
    else:
        kernels += pset.Kernel(periodicBC)
    
   
    if do_mixing:
        kernels += markov_0_mixing
    if do_biofouling:
        kernels += MOi_biofouling
    if do_permanent_fouling:
        kernels += MOi_permanent_biofouling
    if do_Stokes:
        kernels += Stokes_drift
    if windage > 0:
        kernels += windage_drift
    if do_Stokes or windage > 0:
        kernels += unbeaching
        
    print('-------------------Ready to execute main integration-------------------')      
    pset.execute(kernels, runtime=timedelta(days=dt_transition), dt=timedelta(minutes=dt_mins), output_file=pfile,
                 verbose_progress=True, recovery={ErrorCode.ErrorOutOfBounds: delete_particle, ErrorCode.ErrorInterpolation: delete_particle_interp})
    pfile.close()
    
    print('--------Output written to:--------')
    print(os.path.join(dir_output_files,'%s.nc' % filename_out))