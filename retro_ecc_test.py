# retro_ecc_test.py

import numpy as np
import scipy

from inputs import ReadInputs
from setup import setupdiskblackholes
from physics.eccentricity import orbital_ecc

def main():

    rng = np.random.default_rng(1234)
    gamma = -3.0

    fname = "inputs/model_choice.txt"
    
    mass_smbh, trap_radius, disk_outer_radius, alpha, n_bh, mode_mbh_init, max_initial_bh_mass, \
         mbh_powerlaw_index, mu_spin_distribution, sigma_spin_distribution, \
             spin_torque_condition, frac_Eddington_ratio, max_initial_eccentricity, \
                 timestep, number_of_timesteps, disk_model_radius_array, disk_inner_radius,\
                     disk_outer_radius, surface_density_array, aspect_ratio_array, retro, feedback, capture_time, outer_capture_radius, crit_ecc, \
                        r_nsc_out, M_nsc, r_nsc_crit, nbh_nstar_ratio, mbh_mstar_ratio, nsc_index_inner, nsc_index_outer, h_disk_average\
                        = ReadInputs.ReadInputs_ini(fname)

    # create surface density & aspect ratio functions from input arrays
    surf_dens_func_log = scipy.interpolate.UnivariateSpline(disk_model_radius_array, np.log(surface_density_array))
    surf_dens_func = lambda x, f=surf_dens_func_log: np.exp(f(x))

    aspect_ratio_func_log = scipy.interpolate.UnivariateSpline(disk_model_radius_array, np.log(aspect_ratio_array))
    aspect_ratio_func = lambda x, f=aspect_ratio_func_log: np.exp(f(x))

    #Set up number of BH in disk
    n_bh = setupdiskblackholes.setup_disk_nbh(M_nsc,nbh_nstar_ratio,mbh_mstar_ratio,r_nsc_out,nsc_index_outer,mass_smbh,disk_outer_radius,h_disk_average,r_nsc_crit,nsc_index_inner)

    # generate initial BH parameter arrays
    #print("Generate initial BH parameter arrays")
    bh_initial_locations = setupdiskblackholes.setup_disk_blackholes_location(rng, n_bh, disk_outer_radius)
    bh_initial_masses = setupdiskblackholes.setup_disk_blackholes_masses(rng, n_bh, mode_mbh_init, max_initial_bh_mass, mbh_powerlaw_index)
    bh_initial_spins = setupdiskblackholes.setup_disk_blackholes_spins(rng, n_bh, mu_spin_distribution, sigma_spin_distribution)
    #bh_initial_spin_angles = setupdiskblackholes.setup_disk_blackholes_spin_angles(rng, n_bh, bh_initial_spins)
    bh_initial_orb_ang_mom = setupdiskblackholes.setup_disk_blackholes_orb_ang_mom(rng, n_bh)
    bh_initial_orb_ecc = setupdiskblackholes.setup_disk_blackholes_eccentricity_uniform(rng, n_bh)
    #bh_initial_orb_incl = setupdiskblackholes.setup_disk_blackholes_inclination(rng, n_bh)
    #bh_initial_generations = np.ones((n_bh,),dtype=int)

    # assign functions to variable names (continuity issue)
    # Disk surface density (in kg/m^2) is a function of radius, where radius is in r_g
    disk_surface_density = surf_dens_func
    # and disk aspect ratio is also a function of radius, where radius is in r_g
    disk_aspect_ratio = aspect_ratio_func
    # Housekeeping: Set up time
    initial_time = 0.0
    final_time = timestep*number_of_timesteps

    # Find prograde BH orbiters. Identify BH with orb. ang mom =+1
    bh_orb_ang_mom_indices = np.array(bh_initial_orb_ang_mom)
    #prograde_orb_ang_mom_indices = np.where(bh_orb_ang_mom_indices == 1)
    retrograde_orb_ang_mom_indices = np.where(bh_orb_ang_mom_indices == -1)
    retrograde_bh_locations = bh_initial_locations[retrograde_orb_ang_mom_indices]
    #prograde_bh_locations = bh_initial_locations[prograde_orb_ang_mom_indices]
    retrograde_bh_masses = bh_initial_masses[retrograde_orb_ang_mom_indices]
    retrograde_bh_orb_ecc = bh_initial_orb_ecc[retrograde_orb_ang_mom_indices]

    print("Retrograde orbital eccentricity before")
    print(retrograde_bh_orb_ecc)
    retrograde_bh_orb_ecc = orbital_ecc.delta_retro_ecc(mass_smbh, retrograde_bh_locations, retrograde_bh_masses, retrograde_bh_orb_ecc, disk_surface_density, disk_aspect_ratio, timestep, gamma)
    print("Retrograde orbital eccentricity after")
    print(retrograde_bh_orb_ecc)
    

if __name__ == "__main__":
    main()
