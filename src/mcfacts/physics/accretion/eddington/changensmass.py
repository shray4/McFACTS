import numpy as np


def change_mass(prograde_ns_masses, frac_Eddington_ratio, mass_growth_Edd_rate, timestep):
    """Given initial neutron star masses at start of timestep, add mass according to
        chosen NS mass accretion prescription

    Parameters
    ----------
    prograde_ns_masses : float array
        initial masses of neutron stars in prograde orbits around SMBH, units of solar masses
    frac_Eddington_ratio : float
        user chosen input set by input file; Accretion rate of fully embedded neutron star
        in units of Eddington accretion rate. 1.0=embedded NS accreting at Eddington.
        Super-Eddington accretion rates are permitted.
    mass_growth_Edd_rate : float
        fractional rate of mass growth AT Eddington accretion rate per year (2.3e-8)
    timestep : float
        length of timestep in units of years

    Returns
    -------
    ns_new_masses : float array
        masses of neutron stars after accreting at prescribed rate for one timestep
    """
    # Mass grows exponentially for length of timestep:
    ns_new_masses = prograde_ns_masses*np.exp(mass_growth_Edd_rate*frac_Eddington_ratio*timestep)

    return ns_new_masses
