import numpy as np

#NS feedback function. Saavik says feedback shouldn't be super important, hence why the whole thing is an array of zeros, but this may need to be reconsidered later.

def feedback_hankla(prograde_ns_locations, disk_surf_model, frac_Eddington_ratio, alpha):
    """_summary_
    This feedback model uses Eqn. 28 in Hankla, Jiang & Armitage (2020)
    which yields the ratio of heating torque to migration torque.
    Heating torque is directed outwards. 
    So, Ratio <1, slows the inward migration of an object. Ratio>1 sends the object migrating outwards.
    The direction & magnitude of migration (effected by feedback) will be executed in type1.py.

    The ratio of torque due to heating to Type 1 migration torque is calculated as
    R   = Gamma_heat/Gamma_mig 
        ~ 0.07 (speed of light/ Keplerian vel.)(Eddington ratio)(1/optical depth)(1/alpha)^3/2
    where Eddington ratio can be >=1 or <1 as needed,
    optical depth (tau) = Sigma* kappa
    alpha = disk viscosity parameter (e.g. alpha = 0.01 in Sirko & Goodman 2003)
    kappa = 10^0.76 cm^2 g^-1=5.75 cm^2/g = 0.575 m^2/kg for most of Sirko & Goodman disk model (see Fig. 1 & sec 2)
    but e.g. electron scattering opacity is 0.4 cm^2/g
    So tau = Sigma*0.575 where Sigma is in kg/m^2.
    Since v_kep = c/sqrt(a(r_g)) then
    R   ~ 0.07 (a(r_g))^{1/2}(Edd_ratio) (1/tau) (1/alpha)^3/2
    So if assume a=10^3r_g, Sigma=7.e6kg/m^2, alpha=0.01, tau=0.575*Sigma (SG03 disk model), Edd_ratio=1, 
    R   ~5.5e-4 (a/10^3r_g)^(1/2) (Sigma/7.e6) v.small modification to in-migration at a=10^3r_g
        ~0.243 (R/10^4r_g)^(1/2) (Sigma/5.e5)  comparable.
        >1 (a/2x10^4r_g)^(1/2)(Sigma/) migration is *outward* at >=20,000r_g in SG03
        >10 (a/7x10^4r_g)^(1/2)(Sigma/) migration outwards starts to runaway in SG03

    TO DO: Need alpha as an input for disk model (alpha=0.01 is SG03 default)
    TO (MAYBE) DO: kappa default as an input? Or kappa table? Or kappa user set?
    
    Parameters
    ----------
    
    prograde_ns_locations : float array
        locations of prograde singleton NS at start of timestep in units of gravitational radii (r_g=GM_SMBH/c^2)
    disk_surf_model : function
        returns AGN gas disk surface density in kg/m^2 given a distance from the SMBH in r_g
        can accept a simple float (constant), but this is deprecated

    Returns
    -------
    ratio_feedback_to_mig : float array
        ratio of feedback torque to migration torque for each entry in prograde_ns_locations
    """
    # get surface density function, or deal with it if only a float
    if isinstance(disk_surf_model, float):
        disk_surface_density = disk_surf_model
    else:
        disk_surface_density = disk_surf_model(prograde_ns_locations)

    #Calculate ratio
    #
    #Define kappa (or set up a function to call). 
    #kappa = 10^0.76 cm^2/g = 10^(0.76) (10^-2m)^2/10^-3kg=10^(0.76-1)=10^(-0.24) m^2/kg to match units of Sigma
    kappa = 10**(-0.24)
    #Define alpha parameter for disk in Readinputs.py
    #alpha = 0.01

    Ratio_feedback_ns_migration_torque = np.zeros(len(prograde_ns_locations))

    #print((1/kappa),((alpha)**(-1.5)),frac_Eddington_ratio)
    #print("Ratio", Ratio_feedback_migration_torque) 
    #print("NS locations", prograde_ns_locations) 

    return Ratio_feedback_ns_migration_torque  
