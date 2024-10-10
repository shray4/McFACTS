"""
Module for EMRI specific calculations.
"""
import numpy as np
from astropy import units as u

from mcfacts.physics.gw import gw_strain_freq


def evolve_emri_gw(blackholes_inner_disk, timestep_duration_yr, old_gw_freq, smbh_mass, agn_redshift):
    """
    This function evaluates the EMRI gravitational wave frequency and strain at the end of each timestep_duration_yr

    Parameters
    ----------
    blackholes_inner_disk : AGNBlackHole
        black holes in the inner disk
    timestep_duration_yr : float
        timestep in years
    old_gw_freq : numpy array
        previous GW frequency in Hz
    smbh_mass : float
        mass of the SMBH in Msun
    agn_redshift : float
        redshift of the AGN
    """

    old_gw_freq = old_gw_freq * u.Hz

    # If number of EMRIs has grown since last timestep_duration_yr, add a new component to old_gw_freq to carry out dnu/dt calculation
    while (blackholes_inner_disk.num < len(old_gw_freq)):
        old_gw_freq = np.delete(old_gw_freq, 0)
    while blackholes_inner_disk.num > len(old_gw_freq):
        old_gw_freq = np.append(old_gw_freq, (9.e-7) * u.Hz)

    char_strain, nu_gw = gw_strain_freq(mass_1=smbh_mass,
                                        mass_2=blackholes_inner_disk.mass,
                                        obj_sep=blackholes_inner_disk.orb_a,
                                        timestep_duration_yr=timestep_duration_yr,
                                        old_gw_freq=old_gw_freq,
                                        smbh_mass=smbh_mass,
                                        agn_redshift=agn_redshift,
                                        flag_include_old_gw_freq=1)

    return (char_strain, nu_gw)
