import numpy as np

def orbital_ecc_damping(mass_smbh, prograde_bh_locations, prograde_bh_masses, disk_surf_model, disk_aspect_ratio_model, bh_orb_ecc, timestep, crit_ecc):
    """"Return array of BH orbital eccentricities damped according to a prescription
    Use Tanaka & Ward (2004)  t_damp = M^3/2 h^4 / (2^1/2 m Sigma a^1/2 G )
     where M is the central mass, h is the disk aspect ratio (H/a), m is the orbiter mass, 
     Sigma is the disk surface density, a is the semi-major axis, G is the universal gravitational constant.
     From McKernan & Ford (2023) eqn 4. we can parameterize t_damp as 
        t_damp ~ 0.1Myr (q/10^-7)^-1 (h/0.03)^4 (Sigma/10^5 kg m^-2)^-1 (a/10^4r_g)^-1/2

     For eccentricity e<2h
        e(t)=e0*exp(-t/t_damp)......(1)
        So 
            in 0.1 damping time, e(t_damp)=0.90*e0
            in 1 damping time,  e(t_damp)=0.37*e0
            in 2 damping times, e(t_damp)=0.135*e0
            in 3 damping times, e(t_damp)=0.05*e0

    For now assume e<2h condition. To do: Add condition below (if ecc>2h..)

     For eccentricity e>2h eqn. 9 in McKernan & Ford (2023), based on Horn et al. (2012) the scaling time is now t_ecc.
        t_ecc = (t_damp/0.78)*[1 - (0.14*(e/h)^2) + (0.06*(e/h)^3)] ......(2)
        which in the limit of e>0.1 for most disk models becomes
        t_ecc ~ (t_damp/0.78)*[1 + (0.06*(e/h)^3)] 

        Parameters
    ----------
    mass_smbh : float
        mass of supermassive black hole in units of solar masses
    prograde_bh_locations : float array
        locations of prograde singleton BH at start of timestep in units of gravitational radii (r_g=GM_SMBH/c^2)
    prograde_bh_masses : float array
        mass of prograde singleton BH at start of timestep in units of solar masses
    disk_surf_model : function
        returns AGN gas disk surface density in kg/m^2 given a distance from the SMBH in r_g
        can accept a simple float (constant), but this is deprecated
    disk_aspect_ratio_model : function
        returns AGN gas disk aspect ratio given a distance from the SMBH in r_g
        can accept a simple float (constant), but this is deprecated
    bh_orb_ecc : float array
        orbital eccentricity of singleton BH     
    timestep : float
        size of timestep in years
        
        Returns
    -------
    bh_new_orb_ecc : float array
        updated orbital eccentricities damped by AGN gas
        """
    # get surface density function, or deal with it if only a float
    if isinstance(disk_surf_model, float):
        disk_surface_density = disk_surf_model
    else:
        disk_surface_density = disk_surf_model(prograde_bh_locations)
    # ditto for aspect ratio
    if isinstance(disk_aspect_ratio_model, float):
        disk_aspect_ratio = disk_aspect_ratio_model
    else:
        disk_aspect_ratio = disk_aspect_ratio_model(prograde_bh_locations)
    #Set up new_bh_orb_ecc
    new_bh_orb_ecc=np.empty_like(bh_orb_ecc)
    
    #Calculate & normalize all the parameters above in t_damp
    # E.g. normalize q=bh_mass/smbh_mass to 10^-7 
    mass_ratio = prograde_bh_masses/mass_smbh
    
    normalized_mass_ratio = mass_ratio/10**(-7)
    normalized_bh_locations = prograde_bh_locations/1.e4
    normalized_disk_surf_model = disk_surface_density/1.e5
    normalized_aspect_ratio = disk_aspect_ratio/0.03

    #Assume all incoming eccentricities are prograde (for now)
    prograde_bh_orb_ecc = bh_orb_ecc

    #Calculate (e/h) ratio for all prograde BH for use in eqn. 2 above
    e_h_ratio = prograde_bh_orb_ecc/disk_aspect_ratio
    
    # Modest orb eccentricities: e < 2h (experience simple exponential damping): mask entries > 2*aspect_ratio; only show BH with e<2h
    prograde_bh_modest_ecc = np.ma.masked_where(prograde_bh_orb_ecc > 2.0*disk_aspect_ratio_model(prograde_bh_locations),prograde_bh_orb_ecc)
    # Large orb eccentricities: e > 2h (experience more complicated damping)
    prograde_bh_large_ecc = np.ma.masked_where(prograde_bh_orb_ecc < 2.0*disk_aspect_ratio_model(prograde_bh_locations),prograde_bh_orb_ecc)
    #Indices of orb eccentricities where e<2h
    modest_ecc_prograde_indices = np.ma.nonzero(prograde_bh_modest_ecc)
    #Indices of orb eccentricities where e>2h
    large_ecc_prograde_indices = np.ma.nonzero(prograde_bh_large_ecc)
    #print('modest ecc indices', modest_ecc_prograde_indices)
    #print('large ecc indices', large_ecc_prograde_indices)
    #Calculate the 1-d array of damping times at all locations since we need t_damp for both modest & large ecc (see eqns above)
    t_damp =1.e5*(1.0/normalized_mass_ratio)*(normalized_aspect_ratio**4)*(1.0/normalized_disk_surf_model)*(1.0/np.sqrt(normalized_bh_locations))
    #timescale ratio for modest ecc damping
    modest_timescale_ratio = timestep/t_damp
    #timescale for large ecc damping from eqn. 2 above
    t_ecc = (t_damp/0.78)*(1 - (0.14*(e_h_ratio)**(2.0)) + (0.06*(e_h_ratio)**(3.0)))
    large_timescale_ratio = timestep/t_ecc
    #print("t_damp",timestep/t_damp)
    #print("t_ecc",timestep/t_ecc)
    
    #print("timescale_ratio",timescale_ratio)
    new_bh_orb_ecc[modest_ecc_prograde_indices] = bh_orb_ecc[modest_ecc_prograde_indices]*np.exp(-modest_timescale_ratio[modest_ecc_prograde_indices])
    new_bh_orb_ecc[large_ecc_prograde_indices] = bh_orb_ecc[large_ecc_prograde_indices]*np.exp(-large_timescale_ratio[large_ecc_prograde_indices])
    new_bh_orb_ecc = np.where(new_bh_orb_ecc<crit_ecc, crit_ecc,new_bh_orb_ecc)
    #print("Old ecc, New ecc",bh_orb_ecc,new_bh_orb_ecc)
    return new_bh_orb_ecc
