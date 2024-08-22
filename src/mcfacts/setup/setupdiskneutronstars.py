import numpy as np


def setup_disk_neutronstars_location(rng, n_ns, disk_outer_radius):
    #Return an array of NS locations distributed randomly uniformly in disk
    integer_nns = int(n_ns)
    ns_initial_locations = disk_outer_radius*rng.random(integer_nns)
    return ns_initial_locations

def setup_prior_neutronstars_indices(rng, prograde_n_ns, prior_ns_locations):
    #Return an array of indices which allow us to read prior NS properties & replace prograde NS with these.
    integer_nns = int(prograde_n_ns)
    len_prior_locations = (prior_ns_locations.size)-1
    ns_indices = np.rint(len_prior_locations*rng.random(integer_nns))
    return ns_indices

def setup_disk_neutronstars_masses(n_ns):
    #Return an array of NS initial masses for a given powerlaw index and max mass
    integer_nns = int(n_ns)
    ns_initial_masses = 1.4 * np.ones(integer_nns)
    return ns_initial_masses


def setup_disk_neutronstars_spins(rng, n_ns, mu_spin_distribution, sigma_spin_distribution):
    #Return an array of NS initial spin magnitudes for a given mode and sigma of a distribution
    integer_nns = int(n_ns)
    ns_initial_spins = rng.normal(mu_spin_distribution, sigma_spin_distribution, integer_nns)
    return ns_initial_spins


def setup_disk_neutronstars_spin_angles(rng, n_ns, ns_initial_spins):
    #Return an array of NS initial spin angles (in radians).
    #Positive (negative) spin magnitudes have spin angles [0,1.57]([1.5701,3.14])rads
    #All NS spin angles drawn from [0,1.57]rads and +1.57rads to negative spin indices
    integer_nns = int(n_ns)
    ns_initial_spin_indices = np.array(ns_initial_spins)
    negative_spin_indices = np.where(ns_initial_spin_indices < 0.)
    ns_initial_spin_angles = rng.uniform(0.,1.57,integer_nns)
    ns_initial_spin_angles[negative_spin_indices] = ns_initial_spin_angles[negative_spin_indices] + 1.57
    return ns_initial_spin_angles


def setup_disk_neutronstars_orb_ang_mom(rng, n_ns):
    #Return an array of NS initial orbital angular momentum.
    #Assume either fully prograde (+1) or retrograde (-1)
    integer_nns = int(n_ns)
    random_uniform_number = rng.random((integer_nns,))
    ns_initial_orb_ang_mom = (2.0*np.around(random_uniform_number)) - 1.0
    return ns_initial_orb_ang_mom

def setup_disk_neutronstars_eccentricity_thermal(rng, n_ns):
    # Return an array of NS orbital eccentricities
    # For a thermal initial distribution of eccentricities, select from a uniform distribution in e^2.
    # Thus (e=0.7)^2 is 0.49 (half the eccentricities are <0.7). 
    # And (e=0.9)^2=0.81 (about 1/5th eccentricities are >0.9)
    # So rnd= draw from a uniform [0,1] distribution, allows ecc=sqrt(rnd) for thermal distribution.
    # Thermal distribution in limit of equipartition of energy after multiple dynamical encounters
    integer_nns = int(n_ns)
    random_uniform_number = rng.random((integer_nns,))
    ns_initial_orb_ecc = np.sqrt(random_uniform_number)
    return ns_initial_orb_ecc

def setup_disk_neutronstars_eccentricity_uniform(rng, n_ns):
    # Return an array of NS orbital eccentricities
    # For a uniform initial distribution of eccentricities, select from a uniform distribution in e.
    # Thus half the eccentricities are <0.5
    # And about 1/10th eccentricities are >0.9
    # So rnd = draw from a uniform [0,1] distribution, allows ecc = rnd for uniform distribution
    # Most real clusters/binaries lie between thermal & uniform (e.g. Geller et al. 2019, ApJ, 872, 165)
    integer_nns = int(n_ns)
    random_uniform_number = rng.random((integer_nns,))
    ns_initial_orb_ecc = random_uniform_number
    return ns_initial_orb_ecc

def setup_disk_neutronstars_eccentricity_uniform_modified(rng, mod_factor, n_ns):
    # Return an array of NS orbital eccentricities
    # For a uniform initial distribution of eccentricities, select from a uniform distribution in e.
    # Thus half the eccentricities are <0.5
    # And about 1/10th eccentricities are >0.9
    # So rnd = draw from a uniform [0,1] distribution, allows ecc = rnd for uniform distribution
    # Most real clusters/binaries lie between thermal & uniform (e.g. Geller et al. 2019, ApJ, 872, 165)
    integer_nns = int(n_ns)
    random_uniform_number = rng.random((integer_nns,))
    ns_initial_orb_ecc = mod_factor*random_uniform_number
    return ns_initial_orb_ecc

def setup_disk_neutronstars_inclination(rng, n_ns):
    # Return an array of NS orbital inclinations
    # Return an initial distribution of inclination angles that are 0.0
    #
    # To do: initialize inclinations so random draw with i <h (so will need to input ns_locations and disk_aspect_ratio)
    # and then damp inclination.
    # To do: calculate v_kick for each merger and then the (i,e) orbital elements for the newly merged NS. 
    # Then damp (i,e) as appropriate
    integer_nns = int(n_ns)
    # For now, inclinations are zeros
    ns_initial_orb_incl = np.zeros((integer_nns,),dtype = float)
    return ns_initial_orb_incl

def setup_disk_neutronstars_circularized(rng, n_ns,crit_ecc):
    # Return an array of NS orbital inclinations
    # Return an initial distribution of inclination angles that are 0.0
    #
    # To do: initialize inclinations so random draw with i <h (so will need to input ns_locations and disk_aspect_ratio)
    # and then damp inclination.
    # To do: calculate v_kick for each merger and then the (i,e) orbital elements for the newly merged ns. 
    # Then damp (i,e) as appropriate
    integer_nns = int(n_ns)
    # For now, inclinations are zeros
    #ns_initial_orb_ecc = crit_ecc*np.ones((integer_nns,),dtype = float)
    #Try zero eccentricities
    ns_initial_orb_ecc = crit_ecc*np.zeros((integer_nns,),dtype = float)
    return ns_initial_orb_ecc

def setup_disk_nns(M_nsc,nbh_nns_ratio,mns_mstar_ratio,r_nsc_out,nsc_index_outer,mass_smbh,disk_outer_radius,h_disk_average,r_nsc_crit,nsc_index_inner):
    # Return the integer number of NS in the AGN disk as calculated from NSC inputs assuming isotropic distribution of NSC orbits
    # To do: Calculate when R_disk_outer is not equal to the r_nsc_crit
    # To do: Calculate when disky NSC population of NS in plane/out of plane.
    # Housekeeping:
    # Convert outer disk radius in r_g to units of pc. 1r_g =1AU (M_smbh/10^8Msun) and 1pc =2e5AU =2e5 r_g(M/10^8Msun)^-1
    pc_dist = 2.e5*((mass_smbh/1.e8)**(-1.0))
    critical_disk_radius_pc = disk_outer_radius/pc_dist
    #Total average mass of NS in NSC
    M_ns_nsc = M_nsc * nbh_nns_ratio * mns_mstar_ratio
    #print("M_ns_nsc",M_ns_nsc)
    #Total number of NS in NSC
    N_ns_nsc = M_ns_nsc / mns_mstar_ratio
    #print("N_ns_nsc",N_ns_nsc)
    #Relative volumes:
    #   of central 1 pc^3 to size of NSC
    relative_volumes_at1pc = (1.0/r_nsc_out)**(3.0)
    #   of r_nsc_crit^3 to size of NSC
    relative_volumes_at_r_nsc_crit = (r_nsc_crit/r_nsc_out)**(3.0)
    #print(relative_volumes_at1pc)
    #Total number of NS 
    #   at R<1pc (should be about 10^4 for Milky Way parameters; 3x10^7Msun, 5pc, r^-5/2 in outskirts)
    N_ns_nsc_pc = N_ns_nsc * relative_volumes_at1pc * (1.0/r_nsc_out)**(-nsc_index_outer)
    #   at r_nsc_crit
    N_ns_nsc_crit = N_ns_nsc * relative_volumes_at_r_nsc_crit * (r_nsc_crit/r_nsc_out)**(-nsc_index_outer)
    #print("Normalized N_ns at 1pc",N_ns_nsc_pc)
    
    #Calculate Total number of NS in volume R < disk_outer_radius, assuming disk_outer_radius<=1pc.
    
    if critical_disk_radius_pc >= r_nsc_crit:
        relative_volumes_at_disk_outer_radius = (critical_disk_radius_pc/1.0)**(3.0)
        Nns_disk_volume = N_ns_nsc_pc * relative_volumes_at_disk_outer_radius * ((critical_disk_radius_pc/1.0)**(-nsc_index_outer))          
    else:
        relative_volumes_at_disk_outer_radius = (critical_disk_radius_pc/r_nsc_crit)**(3.0)
        Nns_disk_volume = N_ns_nsc_crit * relative_volumes_at_disk_outer_radius * ((critical_disk_radius_pc/r_nsc_crit)**(-nsc_index_inner))
     
    # Total number of NS in disk
    Nns_disk_total = np.rint(Nns_disk_volume * h_disk_average)
    #print("Nns_disk_total",Nns_disk_total)  
    return np.int64(Nns_disk_total)

