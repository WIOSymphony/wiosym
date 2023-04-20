#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 16:31:18 2021
v2: grazing added to biofouling kernel
@author: kaandorp
"""
from parcels import ParcelsRandom, Variable, JITParticle
import math
import numpy as np


class plastic_particle(JITParticle):
    # properties necessary for vertical mixing
    sw_density = Variable('sw_density', dtype=np.float32, to_write=False)
    sw_surface_density = Variable('sw_surface_density', dtype=np.float32, to_write=False)

    l_pl = Variable('l_pl', dtype=np.float32, to_write=False)
    rho_pl = Variable('rho_pl', dtype=np.float32, to_write=False)

    v_s = Variable('v_s', dtype=np.float64, to_write=False)
   
    # biofouling related properties
    a = Variable('a', dtype=np.float32, to_write=False)
    
    
#     s_unbeached = Variable('s_unbeached', dtype=np.float32, to_write=True)
#     s_unbeach_u = Variable('s_unbeach_u', dtype=np.float32, to_write=True)
#     s_unbeach_v = Variable('s_unbeach_v', dtype=np.float32, to_write=True)

#     s_stokes = Variable('s_stokes', dtype=np.float32, to_write=True)
#     s_stokes_u = Variable('s_stokes_u', dtype=np.float32, to_write=True)
#     s_stokes_v = Variable('s_stokes_v', dtype=np.float32, to_write=True)
    
#     s_interp_lon = Variable('s_interp_lon', dtype=np.float32, to_write=True)
#     s_interp_lat = Variable('s_interp_lat', dtype=np.float32, to_write=True)
    
#     s_oob_lon = Variable('s_oob_lon', dtype=np.float32, to_write=True)
#     s_oob_lat = Variable('s_oob_lat', dtype=np.float32, to_write=True)

    hit_bottom = Variable('hit_bottom', dtype=np.int32, to_write=True, initial=0)
    below_500 = Variable('below_500', dtype=np.int32, to_write=True, initial=0)

    
def periodicBC(particle, fieldset, time):
    if particle.lon <= -180.:
        particle.lon += 360.
    elif particle.lon >= 180.:
        particle.lon -= 360.

        
# def periodicBC(particle, fieldset, time):
#     if particle.lon < fieldset.halo_west:
#         particle.lon += fieldset.halo_east - fieldset.halo_west
#     elif particle.lon > fieldset.halo_east:
#         particle.lon -= fieldset.halo_east - fieldset.halo_west


def delete_at_depth(particle, fieldset, time):
    if particle.depth > 500:
        particle.below_500 = 1
        # particle.delete()
        
def delete_particle(particle, fieldset, time):
    """Kernel for deleting particles if they are out of bounds."""
    if fieldset.verbose_delete == 1:
        print('particle is deleted out of bounds at lon = ' + str(particle.lon) + ', lat =' + str(
            particle.lat) + ', depth =' + str(particle.depth))

    # s_oob_lon = particle.lon
    # s_oob_lat = particle.lat
    particle.delete()

def remove_at_bounds(particle, fieldset, time):
    """Kernel for deleting particles if they are out of bounds in the small test run."""
    flag_ = False
    if particle.lat < fieldset.lat_min:
        particle.delete()
        flag_ = True
    if particle.lat > fieldset.lat_max:
        particle.delete()        
        flag_ = True
    if particle.lon < fieldset.lon_min:
        particle.delete()
        flag_ = True
    if particle.lon > fieldset.lon_max:
        particle.delete()
        flag_ = True
    if flag_:
        print('particle deleted at bounds (test run)') 
        
def delete_particle_interp(particle, fieldset, time):
    """Kernel for deleting particles if they are out of bounds."""
    if fieldset.verbose_delete == 1:
        print('particle is deleted due to an interpolation error at lon = ' + str(particle.lon) + ', lat =' + str(
            particle.lat) + ', depth =' + str(particle.depth))
    
    # s_interp_lon = particle.lon
    # s_interp_lat = particle.lat
    particle.delete()
  
    
def initialize_neutral_bouyancy(particle, fieldset, time):
    rho_sw_ = particle.sw_density
    rho_bf_ = fieldset.rho_bf
    rho_pl_ = particle.rho_pl

    r_pl = 0.5*particle.l_pl
    
    theta_pl = 4. * math.pi * r_pl ** 2.  # surface area of plastic particle [m2]

    vol_pl_ = (4. / 3.) * math.pi * r_pl ** 3.
    vol_a_ = fieldset.V_a
    
    #now, set vol_pl*rho_pl + vol_bf*rho_bf == vol_tot        * rho_sw,
    #           tmp_a         tmp_b  tmp_c    (tmp_d + tmp_b) * tmp_e      
    # such that the particles become neutrally bouyant
    vol_bf_ = (vol_pl_*rho_sw_ - vol_pl_*rho_pl_) / (rho_bf_ - rho_sw_)

    # avoid negative algae concentrations when the plastic density is already lower than sea water (might e.g. occur in low saline waters such as Baltic)
    a_neutral_bouyancy = vol_bf_ / (vol_a_ * theta_pl)
    if a_neutral_bouyancy < 0:
        a_neutral_bouyancy = 0
        
    particle.a = a_neutral_bouyancy


def settling_velocity(particle, fieldset, time):
    """
    Calculate settling velocity based on plastic properties (l_pl, rho_pl),
    biofilm properties (amount of algae a, rho_bf),
    and seawater properties (rho_sw, sw_(kin)_visc )
    """
    g = fieldset.G  # gravitational acceleration [m s-2]
    v_a = fieldset.V_a  # Volume of 1 algal cell [m-3]
    temp = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]
    rho_sw = particle.sw_density  # [kg m-3]
    
    mu_w = 4.2844E-5 + (1 / ((0.157 * (temp + 64.993) ** 2) - 91.296))
    A = 1.541 + 1.998E-2 * temp - 9.52E-5 * temp ** 2
    B = 7.974 - 7.561E-2 * temp + 4.724E-4 * temp ** 2
    S_sw = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon] / 1000
    sw_visc = mu_w * (1 + A * S_sw + B * S_sw ** 2)
    sw_kin_visc = sw_visc / particle.sw_density

    r_pl = 0.5*particle.l_pl

    # ------ Volumes -----
    v_pl = (4. / 3.) * math.pi * r_pl ** 3.  # volume of plastic [m3]
    theta_pl = 4. * math.pi * r_pl ** 2.  # surface area of plastic particle [m2]

    v_bfa = (v_a * particle.a) * theta_pl  # volume of living biofilm [m3]
    v_tot = v_bfa + v_pl  # volume of total (biofilm + plastic) [m3]
    r_tot = ((v_tot * (3. / (4. * math.pi))) ** (1. / 3.)) # total radius [m]
    
    t_bf = r_tot - r_pl  # biofilm thickness [m]

    rho_tot = (v_pl * particle.rho_pl + v_bfa * fieldset.rho_bf) / v_tot  # total density [kg m-3]
#     delta_rho = (rho_tot - rho_sw) / rho_sw
#     particle.s_delta_rho = delta_rho
    
#     # Use the equations in Poulain et al. (2019) to calculate the rise velocity of a sphere
#     Re = (r_tot*math.fabs(particle.v_s))/sw_kin_visc
#     if Re == 0:
#         C_d = 30.
#     else:
#         C_d = (12./Re + (6./(1.+math.sqrt(2.*Re))) + 0.4)
#     RHS = (8/3)*r_tot*delta_rho*g
    
#     if delta_rho > 0: #particle heavier than water -> sinking -> positive settling velocity
#         vs = (RHS/C_d)**0.5
#     else: #particle lighter than water -> rising -> negative settling velocity (depth decreases)
#         vs = - ((-RHS/C_d)**0.5 )

    dn = 2. * (r_tot)  # equivalent spherical diameter [m], calculated from Dietrich (1982) from A = pi/4 * dn**2
    delta_rho = (rho_tot - rho_sw) / rho_sw  # normalised difference in density between total plastic+bf and seawater[-]
    dstar = (math.fabs(rho_tot - rho_sw) * g * dn ** 3.) / (rho_sw * sw_kin_visc ** 2.)  # [-]

    if dstar > 5e9:
        w_star = 265000
    elif dstar < 0.05:
        w_star = (dstar ** 2.) * 1.71E-4
    else:
        w_star = 10. ** (-3.76715 + (1.92944 * math.log10(dstar)) - (0.09815 * math.log10(dstar) ** 2.) - (
                    0.00575 * math.log10(dstar) ** 3.) + (0.00056 * math.log10(dstar) ** 4.))

    # ------ Settling velocity of particle -----
    if delta_rho > 0:  # sinks
        vs_new = (g * sw_kin_visc * w_star * delta_rho) ** (1. / 3.)
    else:  # rises
        a_del_rho = delta_rho * -1.
        vs_new = -1. * (g * sw_kin_visc * w_star * a_del_rho) ** (1. / 3.)  # m s-1
        
    particle.v_s = vs_new   
    # particle.s_cd = C_d
    # particle.s_RHS = RHS    


def MOi_biofouling(particle, fieldset, time):
    """
    Kernel to compute the vertical velocity (Vs) of particles due to changes in ambient algal concentrations, growth and death of attached algae based on Kooi et al. 2017 
    model settling velocity and MEDUSA 2.0 biofilm dynamics, including modelling of the 3D mesozooplankton grazing of diatoms
    """
    # ------ Constants and algal properties -----
    g = fieldset.G  # gravitational acceleration [m s-2]
    k = fieldset.K  # Boltzmann constant [m2 kg d-2 K-1] now [s-2] (=1.3804E-23)

    v_a = fieldset.V_a  # Volume of 1 algal cell [m-3]
    r20 = fieldset.R20  # respiration rate [s-1]
    q10 = fieldset.Q10  # temperature coefficient respiration [-]
    gamma = fieldset.Gamma  # shear rate [s-1]

    # ------ sample fields ------
    temp = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]
    mu_w = 4.2844E-5 + (1 / ((0.157 * (temp + 64.993) ** 2) - 91.296))
    A = 1.541 + 1.998E-2 * temp - 9.52E-5 * temp ** 2
    B = 7.974 - 7.561E-2 * temp + 4.724E-4 * temp ** 2
    S_sw = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon] / 1000
    sw_visc = mu_w * (1 + A * S_sw + B * S_sw ** 2)
    sw_kin_visc = sw_visc / particle.sw_density
    rho_sw = particle.sw_density  # [kg m-3]
    vs_init = particle.v_s  # settling velocity [m s-1]

    med_C_cell = 2726e-9 # mg of C per cell
    wt_C_ = fieldset.Wt_C # grams C per mol of C 
    
    
    ###---------------Growth--------------------###
    mmol_conc_diatom = fieldset.bio_diatom[time, particle.depth, particle.lat, particle.lon]
    mmol_conc_nondiat = fieldset.bio_nanophy[time, particle.depth, particle.lat, particle.lon] #nanophytoplankton
    
    no_conc_diatom = mmol_conc_diatom * (wt_C_ / med_C_cell) # conversion from [mmol C m-3] to [mg C m-3] to [no. m-3]
    no_conc_nondiat = mmol_conc_nondiat * (wt_C_ / med_C_cell)
    no_conc_total = no_conc_diatom + no_conc_nondiat

    if no_conc_total < 0:
        no_conc_total = 0.
        print('negative diat/non-diat. concentration')
        
    pp_phyto_ = fieldset.pp_phyto[time, particle.depth, particle.lat, particle.lon] # mg C /m3/day
    
    pp_per_cell = pp_phyto_ / no_conc_total # primary productivity per cell, in mg C / cell / day
    
    pp_ncell_per_cell = pp_per_cell * (1 / med_C_cell) #primary productivity in terms of amount of cells per cell, in cells / cell / day
    
    if pp_ncell_per_cell < 0:
        mu_a = 0
    elif pp_ncell_per_cell > 1.85:
        mu_a = 1.85 / 86400. # maximum growth rate
    else:
        mu_a = pp_ncell_per_cell / 86400. #d-1 to s-1
    
    a_growth = mu_a * particle.a #productivity in amount of cells/m2/s
    

    ### ------------Collisions-------------- ###
    r_pl = 0.5*particle.l_pl

    # ------ Volumes -----
    v_pl = (4. / 3.) * math.pi * r_pl ** 3.  # volume of plastic [m3]
    theta_pl = 4. * math.pi * r_pl ** 2.  # surface area of plastic particle [m2]
    r_a = ((3. / 4.) * (v_a / math.pi)) ** (1. / 3.)  # radius of an algal cell [m]

    v_bfa = (v_a * particle.a) * theta_pl  # volume of living biofilm [m3]
    v_tot = v_bfa + v_pl  # volume of total (biofilm + plastic) [m3]
    # t_bf = ((v_tot * (3. / (4. * math.pi))) ** (1. / 3.)) - r_pl  # biofilm thickness [m]
    r_tot = ((v_tot * (3. / (4. * math.pi))) ** (1. / 3.)) # total radius [m]
    t_bf = r_tot - r_pl  # biofilm thickness [m]
  
    # ------ Diffusivity -----
    r_tot = r_pl + t_bf  # total radius [m]
    rho_tot = (v_pl * particle.rho_pl + v_bfa * fieldset.rho_bf) / v_tot  # total density [kg m-3]
    # theta_tot = 4. * math.pi * r_tot ** 2.  # surface area of total [m2]
    d_pl = k * (temp + 273.16) / (6. * math.pi * sw_visc * r_tot)  # diffusivity of plastic particle [m2 s-1]
    d_a = k * (temp + 273.16) / (6. * math.pi * sw_visc * r_a)  # diffusivity of algal cells [m2 s-1]

    # ------ Encounter rates -----
    beta_abrown = 4. * math.pi * (d_pl + d_a) * (r_tot + r_a)  # Brownian motion [m3 s-1]
    beta_ashear = 1.3 * gamma * ((r_tot + r_a) ** 3.)  # advective shear [m3 s-1]
    beta_aset = (1. / 2.) * math.pi * r_tot ** 2. * math.fabs(vs_init)  # differential settling [m3 s-1]
    beta_a = beta_abrown + beta_ashear + beta_aset  # collision rate [m3 s-1]

    # ------ Attached algal growth (Eq. 11 in Kooi et al. 2017) -----
    a_coll = (beta_a * no_conc_diatom) / theta_pl * fieldset.collision_eff  # [no. m-2 s-1] collisions with diatoms


    ### ------------Respiration-------------- ###
    a_resp = (q10 ** ((temp - 20.) / 10.)) * r20 * particle.a  # [no. m-2 s-1] respiration


    a_graze = 0
    
    #--------------------Total-----------------------------
    particle.a += (a_coll + a_growth - a_resp - a_graze) * particle.dt

    dn = 2. * (r_tot)  # equivalent spherical diameter [m], calculated from Dietrich (1982) from A = pi/4 * dn**2
    delta_rho = (rho_tot - rho_sw) / rho_sw  # normalised difference in density between total plastic+bf and seawater[-]
    dstar = (math.fabs(rho_tot - rho_sw) * g * dn ** 3.) / (rho_sw * sw_kin_visc ** 2.)  # [-]

    if dstar > 5e9:
        w_star = 265000
    elif dstar < 0.05:
        w_star = (dstar ** 2.) * 1.71E-4
    else:
        w_star = 10. ** (-3.76715 + (1.92944 * math.log10(dstar)) - (0.09815 * math.log10(dstar) ** 2.) - (
                    0.00575 * math.log10(dstar) ** 3.) + (0.00056 * math.log10(dstar) ** 4.))

    # ------ Settling velocity of particle -----
    if delta_rho > 0:  # sinks
        vs_new = (g * sw_kin_visc * w_star * delta_rho) ** (1. / 3.)
    else:  # rises
        a_del_rho = delta_rho * -1.
        vs_new = -1. * (g * sw_kin_visc * w_star * a_del_rho) ** (1. / 3.)  # m s-1
        
    
    particle.v_s = vs_new    
    

def MOi_permanent_biofouling(particle, fieldset, time):
    """
    Kernel to compute how particles might sink down permanently (respiration is turned off to turn off the oscillations)
    """
    # ------ Constants and algal properties -----
    g = fieldset.G  # gravitational acceleration [m s-2]
    k = fieldset.K  # Boltzmann constant [m2 kg d-2 K-1] now [s-2] (=1.3804E-23)

    v_a = fieldset.V_a  # Volume of 1 algal cell [m-3]
    # r20 = fieldset.R20  # respiration rate [s-1]
    # q10 = fieldset.Q10  # temperature coefficient respiration [-]
    gamma = fieldset.Gamma  # shear rate [s-1]

    # ------ sample fields ------
    temp = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]
    mu_w = 4.2844E-5 + (1 / ((0.157 * (temp + 64.993) ** 2) - 91.296))
    A = 1.541 + 1.998E-2 * temp - 9.52E-5 * temp ** 2
    B = 7.974 - 7.561E-2 * temp + 4.724E-4 * temp ** 2
    S_sw = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon] / 1000
    sw_visc = mu_w * (1 + A * S_sw + B * S_sw ** 2)
    sw_kin_visc = sw_visc / particle.sw_density
    rho_sw = particle.sw_density  # [kg m-3]
    vs_init = particle.v_s  # settling velocity [m s-1]

    med_C_cell = 2726e-9 # mg of C per cell
    wt_C_ = fieldset.Wt_C # grams C per mol of C 
    
    
    ###---------------Growth--------------------###
    mmol_conc_diatom = fieldset.bio_diatom[time, particle.depth, particle.lat, particle.lon]
    mmol_conc_nondiat = fieldset.bio_nanophy[time, particle.depth, particle.lat, particle.lon] #nanophytoplankton
    
    no_conc_diatom = mmol_conc_diatom * (wt_C_ / med_C_cell) # conversion from [mmol C m-3] to [mg C m-3] to [no. m-3]
    no_conc_nondiat = mmol_conc_nondiat * (wt_C_ / med_C_cell)
    no_conc_total = no_conc_diatom + no_conc_nondiat

    if no_conc_total < 0:
        no_conc_total = 0.
        print('negative diat/non-diat. concentration')
        
    pp_phyto_ = fieldset.pp_phyto[time, particle.depth, particle.lat, particle.lon] # mg C /m3/day
    
    pp_per_cell = pp_phyto_ / no_conc_total # primary productivity per cell, in mg C / cell / day
    
    pp_ncell_per_cell = pp_per_cell * (1 / med_C_cell) #primary productivity in terms of amount of cells per cell, in cells / cell / day
    
    if pp_ncell_per_cell < 0:
        mu_a = 0
    elif pp_ncell_per_cell > 1.85:
        mu_a = 1.85 / 86400. # maximum growth rate
    else:
        mu_a = pp_ncell_per_cell / 86400. #d-1 to s-1
    
    a_growth = mu_a * particle.a #productivity in amount of cells/m2/s
    

    ### ------------Collisions-------------- ###
    r_pl = 0.5*particle.l_pl

    # ------ Volumes -----
    v_pl = (4. / 3.) * math.pi * r_pl ** 3.  # volume of plastic [m3]
    theta_pl = 4. * math.pi * r_pl ** 2.  # surface area of plastic particle [m2]
    r_a = ((3. / 4.) * (v_a / math.pi)) ** (1. / 3.)  # radius of an algal cell [m]

    v_bfa = (v_a * particle.a) * theta_pl  # volume of living biofilm [m3]
    v_tot = v_bfa + v_pl  # volume of total (biofilm + plastic) [m3]
    # t_bf = ((v_tot * (3. / (4. * math.pi))) ** (1. / 3.)) - r_pl  # biofilm thickness [m]
    r_tot = ((v_tot * (3. / (4. * math.pi))) ** (1. / 3.)) # total radius [m]
    t_bf = r_tot - r_pl  # biofilm thickness [m]
  
    # ------ Diffusivity -----
    r_tot = r_pl + t_bf  # total radius [m]
    rho_tot = (v_pl * particle.rho_pl + v_bfa * fieldset.rho_bf) / v_tot  # total density [kg m-3]
    # theta_tot = 4. * math.pi * r_tot ** 2.  # surface area of total [m2]
    d_pl = k * (temp + 273.16) / (6. * math.pi * sw_visc * r_tot)  # diffusivity of plastic particle [m2 s-1]
    d_a = k * (temp + 273.16) / (6. * math.pi * sw_visc * r_a)  # diffusivity of algal cells [m2 s-1]

    # ------ Encounter rates -----
    beta_abrown = 4. * math.pi * (d_pl + d_a) * (r_tot + r_a)  # Brownian motion [m3 s-1]
    beta_ashear = 1.3 * gamma * ((r_tot + r_a) ** 3.)  # advective shear [m3 s-1]
    beta_aset = (1. / 2.) * math.pi * r_tot ** 2. * math.fabs(vs_init)  # differential settling [m3 s-1]
    beta_a = beta_abrown + beta_ashear + beta_aset  # collision rate [m3 s-1]

    # ------ Attached algal growth (Eq. 11 in Kooi et al. 2017) -----
    a_coll = (beta_a * no_conc_diatom) / theta_pl * fieldset.collision_eff  # [no. m-2 s-1] collisions with diatoms


    ### ------------Respiration-------------- ###
    # a_resp = (q10 ** ((temp - 20.) / 10.)) * r20 * particle.a  # [no. m-2 s-1] respiration
    a_resp = 0

    a_graze = 0
    
    #--------------------Total-----------------------------
    particle.a += (a_coll + a_growth - a_resp - a_graze) * particle.dt

    dn = 2. * (r_tot)  # equivalent spherical diameter [m], calculated from Dietrich (1982) from A = pi/4 * dn**2
    delta_rho = (rho_tot - rho_sw) / rho_sw  # normalised difference in density between total plastic+bf and seawater[-]
    dstar = (math.fabs(rho_tot - rho_sw) * g * dn ** 3.) / (rho_sw * sw_kin_visc ** 2.)  # [-]

    if dstar > 5e9:
        w_star = 265000
    elif dstar < 0.05:
        w_star = (dstar ** 2.) * 1.71E-4
    else:
        w_star = 10. ** (-3.76715 + (1.92944 * math.log10(dstar)) - (0.09815 * math.log10(dstar) ** 2.) - (
                    0.00575 * math.log10(dstar) ** 3.) + (0.00056 * math.log10(dstar) ** 4.))

    # ------ Settling velocity of particle -----
    if delta_rho > 0:  # sinks
        vs_new = (g * sw_kin_visc * w_star * delta_rho) ** (1. / 3.)
    else:  # rises
        a_del_rho = delta_rho * -1.
        vs_new = -1. * (g * sw_kin_visc * w_star * a_del_rho) ** (1. / 3.)  # m s-1
        
    
    particle.v_s = vs_new    
    
    
    
def neutral_buoyancy(particle, fieldset, time):
    particle.v_s = 0.


def markov_0_mixing(particle, fieldset, time):
    """
    simple markov-0 kernel for vertical mixing, where the deterministic component
    is determined using forward-difference with a given delta_z
    """
    
    delta_z = 1.
    
    kz = fieldset.mixing_kz[time, particle.depth, particle.lat, particle.lon]
    kz_delta = fieldset.mixing_kz[time, particle.depth+delta_z, particle.lat, particle.lon]
    
    dkz_dz = (kz_delta - kz) / delta_z
    
    # According to Ross & Sharples (2004), first the deterministic part of equation 1
    dz_deterministic = dkz_dz * particle.dt

    # The random walk component
    dz_random = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3) * math.sqrt(2 * kz)
    
    # rise velocity component
    dz_wb = particle.v_s * particle.dt

    potential = particle.depth + dz_deterministic + dz_random  + dz_wb
    bathymetry_local = fieldset.bathymetry[time, fieldset.z_start, particle.lat, particle.lon]
    
    if potential < fieldset.z_start:
        particle.depth = fieldset.z_start
    elif potential > bathymetry_local:
        particle.depth += 0 #TO DO: keep particles 'beached'
        particle.hit_bottom = 1
    elif particle.depth > 100 and potential > (bathymetry_local*0.99): # for deeper particles; since bathymetry can be quite rough (and is interpolated linearly) look at the 99% value instead
        particle.depth += 0
        particle.hit_bottom = 1
    elif potential > 3900:
        particle.depth += 0
        particle.hit_bottom = 1
    else:
        particle.depth = potential        
        
        
def Stokes_drift(particle, fieldset, time):
    """
    Stokes drift kernel, using the Breivik (2016) approach assuming
    a Phillips wave spectrum to determine the Stokes drift in depth
    """        
    #U/V components Stokes drift
    U_stokes_ = fieldset.Stokes_U[time, particle.depth, particle.lat, particle.lon]
    V_stokes_ = fieldset.Stokes_V[time, particle.depth, particle.lat, particle.lon]
    
    #peak wave period
    T_p = fieldset.wave_Tp[time, particle.depth, particle.lat, particle.lon]
    
    if T_p > 1e-14:
        #peak wave frequency and wavelength
        omega_p = 2*math.pi / T_p
        k_p = omega_p**2 / fieldset.G

        kp_z_2 = 2*k_p*particle.depth

        decay = math.exp(-kp_z_2) - math.sqrt(math.pi*kp_z_2)*math.erfc(math.sqrt(kp_z_2))

        U_s_z = U_stokes_*decay
        V_s_z = V_stokes_*decay


        dlon = U_s_z * particle.dt
        dlat = V_s_z * particle.dt
        
        # particle.s_stokes = 1
        # particle.s_stokes_u = U_s_z
        # particle.s_stokes_v = V_s_z

    else:
        dlon = 0.
        dlat = 0.

        
    particle.lon += dlon
    particle.lat += dlat

def windage_drift(particle, fieldset, time):
    """
    Simple windage addition
    """
    U_wind_ = fieldset.wind_U[time, particle.depth, particle.lat, particle.lon]
    V_wind_ = fieldset.wind_V[time, particle.depth, particle.lat, particle.lon]
      
    dlon = U_wind_ * particle.dt
    dlat = V_wind_ * particle.dt

    particle.lon += dlon
    particle.lat += dlat
    
def unbeaching(particle, fieldset, time):
    (vel_u, vel_v, vel_w) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
    
    if math.fabs(vel_u) < 1e-14 and math.fabs(vel_v) < 1e-14:
        U_ub = fieldset.unbeach_U[time, particle.depth, particle.lat, particle.lon]
        V_ub = fieldset.unbeach_V[time, particle.depth, particle.lat, particle.lon]

        dlon = U_ub * particle.dt
        dlat = V_ub * particle.dt  

        particle.lon += dlon
        particle.lat += dlat

#         print('unbeaching')
#         print(U_ub)
#         print(V_ub)
        
#         if math.fabs(dlon) > 0.01 or math.fabs(dlat) > 0.01:
#             print('---------large unbeaching value-------------')
#             print(dlon)
#             print(dlat)        

        # particle.s_unbeached = 1
        # particle.s_unbeach_u = U_ub
        # particle.s_unbeach_v = V_ub       
    # else:
    #     particle.s_unbeached = 0
    #     particle.s_unbeach_u = 0
    #     particle.s_unbeach_v = 0
            

def PolyTEOS10_bsq(particle, fieldset, time):
    '''
    # calculates density based on the polyTEOS10-bsq algorithm from Appendix A.2 of
    # https://www.sciencedirect.com/science/article/pii/S1463500315000566
    # requires fieldset.abs_salinity and fieldset.cons_temperature Fields in the fieldset
    # and a particle.density Variable in the ParticleSet
    #
    # References:
    #  Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014: Accurate
    #   polynomial expressions for the density and specific volume of
    #   seawater using the TEOS-10 standard. Ocean Modelling.
    #  McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003:
    #   Accurate and computationally efficient algorithms for potential
    #   temperature and density of seawater.  Journal of Atmospheric and
    #   Oceanic Technology, 20, 730-741.
    '''

    Z = - particle.depth  # note: use negative depths!
    SA = fieldset.abs_salinity[time, particle.depth, particle.lat, particle.lon]
    CT = fieldset.cons_temperature[time, particle.depth, particle.lat, particle.lon]

    SAu = 40 * 35.16504 / 35
    CTu = 40
    Zu = 1e4
    deltaS = 32
    R000 = 8.0189615746e+02
    R100 = 8.6672408165e+02
    R200 = -1.7864682637e+03
    R300 = 2.0375295546e+03
    R400 = -1.2849161071e+03
    R500 = 4.3227585684e+02
    R600 = -6.0579916612e+01
    R010 = 2.6010145068e+01
    R110 = -6.5281885265e+01
    R210 = 8.1770425108e+01
    R310 = -5.6888046321e+01
    R410 = 1.7681814114e+01
    R510 = -1.9193502195e+00
    R020 = -3.7074170417e+01
    R120 = 6.1548258127e+01
    R220 = -6.0362551501e+01
    R320 = 2.9130021253e+01
    R420 = -5.4723692739e+00
    R030 = 2.1661789529e+01
    R130 = -3.3449108469e+01
    R230 = 1.9717078466e+01
    R330 = -3.1742946532e+00
    R040 = -8.3627885467e+00
    R140 = 1.1311538584e+01
    R240 = -5.3563304045e+00
    R050 = 5.4048723791e-01
    R150 = 4.8169980163e-01
    R060 = -1.9083568888e-01
    R001 = 1.9681925209e+01
    R101 = -4.2549998214e+01
    R201 = 5.0774768218e+01
    R301 = -3.0938076334e+01
    R401 = 6.6051753097e+00
    R011 = -1.3336301113e+01
    R111 = -4.4870114575e+00
    R211 = 5.0042598061e+00
    R311 = -6.5399043664e-01
    R021 = 6.7080479603e+00
    R121 = 3.5063081279e+00
    R221 = -1.8795372996e+00
    R031 = -2.4649669534e+00
    R131 = -5.5077101279e-01
    R041 = 5.5927935970e-01
    R002 = 2.0660924175e+00
    R102 = -4.9527603989e+00
    R202 = 2.5019633244e+00
    R012 = 2.0564311499e+00
    R112 = -2.1311365518e-01
    R022 = -1.2419983026e+00
    R003 = -2.3342758797e-02
    R103 = -1.8507636718e-02
    R013 = 3.7969820455e-01
    ss = math.sqrt((SA + deltaS) / SAu)
    tt = CT / CTu
    zz = -Z / Zu
    rz3 = R013 * tt + R103 * ss + R003
    rz2 = (R022 * tt + R112 * ss + R012) * tt + (R202 * ss + R102) * ss + R002
    rz1 = (((R041 * tt + R131 * ss + R031) * tt + (R221 * ss + R121) * ss + R021) * tt + ((R311 * ss + R211) * ss + R111) * ss + R011) * tt + (((R401 * ss + R301) * ss + R201) * ss + R101) * ss + R001
    rz0 = (((((R060 * tt + R150 * ss + R050) * tt + (R240 * ss + R140) * ss + R040) * tt + ((R330 * ss + R230) * ss + R130) * ss + R030) * tt + (((R420 * ss + R320) * ss + R220) * ss + R120) * ss + R020) * tt + ((((R510 * ss + R410) * ss + R310) * ss + R210) * ss + R110) * ss + R010) * tt + (((((R600 * ss + R500) * ss + R400) * ss + R300) * ss + R200) * ss + R100) * ss + R000
    particle.sw_density = ((rz3 * zz + rz2) * zz + rz1) * zz + rz0
    particle.sw_surface_density = rz0
    
    
