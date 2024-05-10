import numpy as np
import scipy

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


def dimnless_energy_integral(gamma, ecc):
    """Numerically computes the dimensionless energy integral required for eccentricity pumping
    for retrograde orbiters. Uses simple trapezoid rule. Eqn 19 from Secunda et al. 2021
    (2021ApJ...908L..27S). Derivation assumes AGN disk midplane density is powerlaw in radius
    (see below definition of gamma). This is close enough to true for SG and TQM, and should
    be true in general... but maybe you have a really strange disk model, so, be aware.

    Parameters
    ----------
    gamma : float
        power-law index for radial distribution of AGN disk midplane density (ie rho(r) \propto r^gamma)
    ecc : float array
        current eccentricities of orbiters

    Returns
    -------
    dless_energy_integral : float array
        value of the dimensionless energy integral (Eqn 19) from Secunda et al. 2021 for each ecc
    """
    # set up output array
    dless_energy_integral = np.zeros_like(ecc)
    # compute value of prefix to integral
    integral_prefix = 1.0/(2.0*np.pi*np.sqrt(1.0-ecc**2))
    # set resolution and range of integration, which goes from 0-2pi around the orbit
    d_phi = 0.01
    phi = np.arange(0.0, 2.0*np.pi, d_phi)
    # for each eccentricity (yes, has to be a for loop I think...)
    for j in range(len(ecc)):
        # set u substitution variable, per Eqn 21
        u_sub_var = (1.0 + ecc[j]*np.cos(phi))/(1.0 - ecc[j]**2)
        # compute denominator of the integrand (note same as for ang mom integral, eqn 20)
        integrand_denom = pow((-1.0 + 3.0*u_sub_var + 2.0*np.sqrt(1.0-ecc[j]**2)*pow(u_sub_var, 1.5)), 1.5)
        # compute numerator of the integrand
        integrand_numerator = (-1.0 + 3.0*u_sub_var + 2.0*np.sqrt(1.0-ecc[j]**2)*pow(u_sub_var, 1.5))*pow(u_sub_var,(-gamma-2.0))
    
        #trapezoid rule integration
        for i in range(len(phi)-1):
            trapezoid = integral_prefix[j]*0.5*((integrand_numerator[i]/integrand_denom[i]) + (integrand_numerator[i+1]/integrand_denom[i+1]))*d_phi
            dless_energy_integral[j] = dless_energy_integral[j] + trapezoid
    
    return dless_energy_integral

def dimnless_angmom_integral(gamma, ecc):
    """Numerically computes the dimensionless angular momentum integral required for eccentricity pumping
    for retrograde orbiters. Uses simple trapezoid rule. Eqn 20 from Secunda et al. 2021
    (2021ApJ...908L..27S). Derivation assumes AGN disk midplane density is powerlaw in radius
    (see below definition of gamma). This is close enough to true for SG and TQM, and should
    be true in general... but maybe you have a really strange disk model, so, be aware.

    Parameters
    ----------
    gamma : float
        power-law index for radial distribution of AGN disk midplane density (ie rho(r) \propto r^gamma)
    ecc : float array
        current eccentricities of orbiters

    Returns
    -------
    dless_angmom_integral : float array
        value of the dimensionless angular momentum integral (Eqn 20) from Secunda et al. 2021 for each ecc
    """
    # set up output array
    dless_angmom_integral = np.zeros_like(ecc)
    # compute value of prefix to integral
    integral_prefix = 1.0/(2.0*np.pi*np.sqrt(1.0-ecc**2))
    # set resolution and range of integration, which goes from 0-2pi around the orbit
    d_phi = 0.01
    phi = np.arange(0.0, 2.0*np.pi, d_phi)
    # for each eccentricity (yes, has to be a for loop I think...)
    for j in range(len(ecc)):
        # set u substitution variable, per Eqn 21
        u_sub_var = (1.0 + ecc[j]*np.cos(phi))/(1.0 - ecc[j]**2)
        # compute denominator of the integrand (note same as for energy integral, eqn 19)
        integrand_denom = pow((-1.0 + 3.0*u_sub_var + 2.0*np.sqrt(1.0-ecc[j]**2)*pow(u_sub_var, 1.5)), 1.5)
        # compute numerator of the integrand
        integrand_numerator = (-np.sqrt(1.0-ecc[j]**2) - pow(u_sub_var, -0.5))*pow(u_sub_var,(-gamma-2.0))

        #trapezoid rule integration
        for i in range(len(phi)-1):
            trapezoid = integral_prefix[j]*0.5*((integrand_numerator[i]/integrand_denom[i]) + (integrand_numerator[i+1]/integrand_denom[i+1]))*d_phi
            dless_angmom_integral[j] = dless_angmom_integral[j] + trapezoid

    return dless_angmom_integral

def delta_retro_ecc(mass_smbh, retrograde_bh_locations, retrograde_bh_masses, retrograde_bh_orb_ecc, disk_surf_model, disk_aspect_ratio_model, timestep, gamma):
    """Computes the changed orbital eccentricities of retrograde orbiters fully embedded in an
    AGN disk, computed approximately via Eqn 18 of Secunda et al. 2021 (2021ApJ...908L..27S), which
    gives de^2/dt averaged over one orbit. The derivation assumes the AGN disk midplane density is 
    powerlaw in radius ie rho(r) \propto r^(gamma). We make a number of possibly stupid approximations
    here (in order from least concerning to most concerning, I think):

    -powerlaw density function (this is ok in outskirts of SG disk where gamma=-3)
    -numerically integrate some dimensionless integrals using a trapezoid rule (cf utility functions)
    -assume the midplane density is the disk surface density/(aspect ratio * semi-major axis)
    -assume the semi-major axis does not change over the whole timestep--OK this is probably an issue
    -assume the rate of change of eccentricity per orbit is the same for the whole timestep--along with this
    -neglect eccentricity damping (and semi-major axis change) due to GW 
    (this is only important for small pericenter passages tho)

    Parameters
    ----------
    mass_smbh : float
        mass of the supermassive black hole in units of solar masses
    retrograde_bh_locations : float array
        semi-major axes of all retrograde orbiting embedded black holes in units of r_g(=G*mass_smbh/c^2)
    retrograde_bh_masses : float array
        masses of all retrograde orbiting embedded black holes in units of solar masses
    retrograde_bh_orb_ecc : float array
        eccentricities of all retrograde orbiting embedded black holes
    disk_surf_model : function
        returns AGN gas disk surface density in kg/m^2 given a distance from the SMBH in r_g
    disk_aspect_ratio_model : function
        returns AGN gas disk aspect ratio given a distance from the SMBH in r_g
    timestep : float
        size of timestep in years
    gamma : float
        powerlaw index of AGN disk midplane surface density function, ie rho(r) \propto r^(gamma)

    Returns
    -------
    retrograde_bh_orb_ecc : float array
        eccentricities of all retrograde orbiting embedded black holes after one timestep
    """
    # compute the period of each orbiter in years, given semi-major axis in r_g and mass_smbh in Msun
    period = 2.0*np.pi*(scipy.constants.G*2.0e30)/scipy.constants.c**3 \
        *(np.sqrt(retrograde_bh_locations**3)*mass_smbh)/3.15e7
    num_orbits = timestep/period
    
    lnLambda = 1.0 # this is probably a reasonable value most of the time
    # compute mass density in kg/m^3 from surface density (in kg/m^2) and scale 
    # height in meters(=aspect ratio*semi-major axis in r_g * r_g)
    rho = disk_surf_model(retrograde_bh_locations)/ \
        (disk_aspect_ratio_model(retrograde_bh_locations)*retrograde_bh_locations* \
         scipy.constants.G*mass_smbh*2.0e30/scipy.constants.c**2)
    # compute prefactor function f(a) Eqn 13 from Secunda et al 2021 (in output in SI--convert semi-major axes from r_g)
    prefactor = 4.0*np.pi*lnLambda*(scipy.constants.G*retrograde_bh_masses*2.0e30)**2*rho* \
        np.sqrt(retrograde_bh_locations/scipy.constants.c**2)

    # now compute de^2 per period (from Eqn 18)--computed in SI as above--convert from r_g and Msun
    # note that we are multiplying by -period because t_final-t_initial is actually -ive period
    # this also yields a positive definite value for d(ecc^2) as expected
    # I think the above is wrong! changing it below!
    # and note also that we have to convert the period back from years to seconds here 
    delta_ecc_sq_per_period = -(prefactor*2.0*retrograde_bh_locations/(retrograde_bh_masses*2.0e30*scipy.constants.c**2)*(1.0-retrograde_bh_orb_ecc**2)* \
                                    (dimnless_energy_integral(gamma, retrograde_bh_orb_ecc) + \
                                     dimnless_angmom_integral(gamma, retrograde_bh_orb_ecc)/np.sqrt(1.0-retrograde_bh_orb_ecc**2)))*(period*3.15e7)
    delta_ecc_sq_per_ts = delta_ecc_sq_per_period*num_orbits
    # delta ecc^2 = ecc_final^2 - ecc_initial^2 = delta_ecc_sq_per_ts
    # solve for ecc_final:
    retrograde_bh_orb_ecc = np.sqrt(delta_ecc_sq_per_ts + retrograde_bh_orb_ecc**2)

    return retrograde_bh_orb_ecc

def delta_retro_smaj_axis(mass_smbh, retrograde_bh_locations, retrograde_bh_masses, retrograde_bh_orb_ecc, disk_surf_model, disk_aspect_ratio_model, timestep, gamma):
    # compute the period of each orbiter in years, given semi-major axis in r_g and mass_smbh in Msun
    period = 2.0*np.pi*(scipy.constants.G*2.0e30)/scipy.constants.c**3 \
        *(np.sqrt(retrograde_bh_locations**3)*mass_smbh)/3.15e7
    num_orbits = timestep/period
    
    lnLambda = 1.0 # this is probably a reasonable value most of the time
    # compute mass density in kg/m^3 from surface density (in kg/m^2) and scale 
    # height in meters(=aspect ratio*semi-major axis in r_g * r_g)
    rho = disk_surf_model(retrograde_bh_locations)/ \
        (disk_aspect_ratio_model(retrograde_bh_locations)*retrograde_bh_locations* \
         scipy.constants.G*mass_smbh*2.0e30/scipy.constants.c**2)
    # compute prefactor function f(a) Eqn 13 from Secunda et al 2021 (in output in SI--convert semi-major axes from r_g)
    prefactor = 4.0*np.pi*lnLambda*(scipy.constants.G*retrograde_bh_masses*2.0e30)**2*rho* \
        np.sqrt(retrograde_bh_locations/scipy.constants.c**2)

    # now compute da per period (from Eqn 17)--computed in R_g--convert from r_g and Msun to get SI, then multiply by R_g
    # note that we are multiplying by -period because t_final-t_initial is actually -ive period
    # this also yields a negative definite value for da as expected
    # I think the above is wrong! changing it below!
    # and note also that we have to convert the period back from years to seconds here 
    delta_smaj_axis_per_period = -prefactor*retrograde_bh_locations*(2.0*retrograde_bh_locations/(retrograde_bh_masses*2.0e30*scipy.constants.c**2))*dimnless_energy_integral(gamma, retrograde_bh_orb_ecc)*(period*3.15e7)
    delta_smaj_axis_per_ts = delta_smaj_axis_per_period*num_orbits

    # check signs, maybe not ok...
    retrograde_bh_locations = retrograde_bh_locations + delta_smaj_axis_per_ts

    return retrograde_bh_locations

def new_fucking_plan():
    