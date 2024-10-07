"""
Module for calculating change of mass, spin magnitude, and spin angle due to accretion.
"""

import numpy as np


def change_star_mass(disk_star_pro_masses,
                     disk_star_eddington_ratio,
                     mdisk_star_eddington_mass_growth_rate,
                     timestep_duration_yr):
    """Given initial star masses at start of timestep, add mass according to
        chosen stellar mass accretion prescription

    Parameters
    ----------
    disk_star_pro_masses : float array
        initial masses of stars in prograde orbits around SMBH, units of solar masses
    disk_star_eddington_ratio : float
        user chosen input set by input file; Accretion rate of fully embedded stars
          in units of Eddington accretion rate. 1.0=embedded star accreting at Eddington.
        Super-Eddington accretion rates are permitted.
    mdisk_star_eddington_mass_growth_rate : float
        fractional rate of mass growth AT Eddington accretion rate per year (2.3e-8)
    timestep_duration_yr : float
        length of timestep in units of years

    Returns
    -------
    disk_star_pro_new_masses : float array
        masses of stars after accreting at prescribed rate for one timestep
    """
    # Mass grows exponentially for length of timestep:
    disk_star_pro_new_masses = disk_star_pro_masses*np.exp(mdisk_star_eddington_mass_growth_rate*disk_star_eddington_ratio*timestep_duration_yr)

    return disk_star_pro_new_masses


def change_bh_mass(disk_bh_pro_masses, disk_bh_eddington_ratio, disk_bh_eddington_mass_growth_rate, timestep_duration_yr):
    """Given initial black hole masses at start of timestep, add mass according to
        chosen BH mass accretion prescription

    Parameters
    ----------
    disk_bh_pro_masses : float array
        initial masses of black holes in prograde orbits around SMBH, units of solar masses
    disk_bh_eddington_ratio : float
        user chosen input set by input file; Accretion rate of fully embedded stellar mass
        black hole in units of Eddington accretion rate. 1.0=embedded BH accreting at Eddington.
        Super-Eddington accretion rates are permitted.
    mdisk_bh_eddington_mass_growth_rate : float
        fractional rate of mass growth AT Eddington accretion rate per year (2.3e-8)
    timestep_duration_yr : float
        length of timestep in units of years

    Returns
    -------
    disk_bh_pro_new_masses : float array
        masses of black holes after accreting at prescribed rate for one timestep
    """
    # Mass grows exponentially for length of timestep:
    disk_bh_pro_new_masses = disk_bh_pro_masses*np.exp(disk_bh_eddington_mass_growth_rate*disk_bh_eddington_ratio*timestep_duration_yr)

    return disk_bh_pro_new_masses


def change_star_spin_magnitudes(disk_star_pro_spins,
                                disk_star_eddington_ratio,
                                disk_star_torque_condition,
                                timestep_duration_yr,
                                disk_star_pro_orbs_ecc,
                                disk_star_pro_orbs_ecc_crit):
    """
    Update the spin magnitude of the embedded stars based on their accreted
    mass in this timestep.
    Right now this is just a copy of the BH function. No changes.

    Parameters
    ----------
    disk_star_pro_spins : float array
        initial spins of black holes in prograde orbits around SMBH
    disk_star_eddington_ratio : float
        user chosen input set by input file; Accretion rate of fully embedded stellar mass
        black hole in units of Eddington accretion rate. 1.0=embedded BH accreting at Eddington.
        Super-Eddington accretion rates are permitted.
    disk_star_torque_condition : float
        user chosen input set by input file; fraction of initial mass required to be
        accreted before BH spin is torqued fully into alignment with the AGN disk.
        We don't know for sure but Bogdanovic et al. says between 0.01=1% and 0.1=10%
        is what is required
    timestep_duration_yr : float
        length of timestep in units of years
    disk_star_pro_orbs_ecc : float array
        orbital eccentricity of BH in prograde orbits around SMBH
    disk_star_pro_orbs_ecc_crit : float
        critical value of orbital eccentricity below which prograde accretion (& migration & binary formation) occurs
    Returns
    -------
    disk_star_pro_spins_new : float array
        spin magnitudes of black holes after accreting at prescribed rate for one timestep
    """
    # A retrograde BH a=-1 will spin down to a=0 when it accretes a factor sqrt(3/2)=1.22 in mass (Bardeen 1970).
    # Since M_edd/t = 2.3 e-8 M0/yr or 2.3e-4M0/10kyr then M(t)=M0*exp((M_edd/t)*f_edd*time)
    # so M(t)~1.2=M0*exp(0.2) so in 10^7yr, spin should go a=-1 to a=0. Or delta a ~ 10^-3 every 10^4yr.

    disk_bh_eddington_ratio_normalized = disk_star_eddington_ratio/1.0
    timestep_duration_yr_normalized = timestep_duration_yr/1.e4
    disk_star_torque_condition_normalized = disk_star_torque_condition/0.1

    spin_iteration = (1.e-3*disk_bh_eddington_ratio_normalized*disk_star_torque_condition_normalized*timestep_duration_yr_normalized)

    disk_star_pro_spins_new = disk_star_pro_spins
    # Singleton star with orb ecc > disk_star_pro_orbs_ecc_crit will spin down b/c accrete retrograde
    disk_star_pro_spin_down = np.ma.masked_where(disk_star_pro_orbs_ecc <= disk_star_pro_orbs_ecc_crit, disk_star_pro_orbs_ecc)
    # Singleton star with orb ecc < disk_star_pro_orbs_ecc_crit will spin up b/c accrete prograde
    disk_star_pro_spin_up = np.ma.masked_where(disk_star_pro_orbs_ecc > disk_star_pro_orbs_ecc_crit, disk_star_pro_orbs_ecc)
    # Indices of singleton star with orb ecc > disk_star_pro_orbs_ecc_crit
    disk_star_pro_indices_spin_down = np.ma.nonzero(disk_star_pro_spin_down)
    # Indices of singleton star with orb ecc < disk_star_pro_orbs_ecc_crit
    disk_star_pro_indices_spin_up = np.ma.nonzero(disk_star_pro_spin_up)
    # disk_star_pro_spins_new[prograde_orb_ang_mom_indices]=disk_star_pro_spins_new[prograde_orb_ang_mom_indices]+(4.4e-3*disk_bh_eddington_ratio_normalized*disk_star_torque_condition_normalized*timestep_duration_yr_normalized)
    disk_star_pro_spins_new[disk_star_pro_indices_spin_up] = disk_star_pro_spins[disk_star_pro_indices_spin_up] + spin_iteration
    # Spin down stars with orb ecc > disk_star_pro_orbs_ecc_crit
    disk_star_pro_spins_new[disk_star_pro_indices_spin_down] = disk_star_pro_spins[disk_star_pro_indices_spin_down] - spin_iteration

    # Housekeeping: max possible spins, do not spin above or below these values.
    disk_star_pro_spin_max = 0.98
    disk_star_pro_spin_min = -0.98

    for i in range(len(disk_star_pro_spins)):
        if (disk_star_pro_spins_new[i] < disk_star_pro_spin_min):
            disk_star_pro_spins_new[i] = disk_star_pro_spin_min

        if (disk_star_pro_spins_new[i] > disk_star_pro_spin_max):
            disk_star_pro_spins_new[i] = disk_star_pro_spin_max

    return disk_star_pro_spins_new


def change_star_spin_angles(prograde_stars_spin_angles,
                            disk_star_eddington_ratio,
                            disk_star_torque_condition,
                            disk_star_spin_minimum_resolution,
                            timestep_duration_yr,
                            disk_star_pro_orbs_ecc,
                            disk_star_pro_orbs_ecc_crit):
    """
    Update the spin angles of the embedded black holes based on their
    accreted mass in this timestep. Right now this is just a copy
    of the BH function, no changes.

    Parameters
    ----------
    prograde_stars_spin_angles : float array
        spin angles of stars in prograde orbits around SMBH
    disk_star_eddington_ratio : float
        user chosen input set by input file; Accretion rate of fully embedded stellar mass
        black hole in units of Eddington accretion rate. 1.0=embedded BH accreting at Eddington.
        Super-Eddington accretion rates are permitted.
    disk_star_torque_condition : float
        user chosen input set by input file; fraction of initial mass required to be
        accreted before BH spin is torqued fully into alignment with the AGN disk.
        We don't know for sure but Bogdanovic et al. says between 0.01=1% and 0.1=10%
        is what is required
    disk_star_spin_minimum_resolution : float
        user chosen input set by input file; minimum resolution of spin change followed by code.
    timestep_duration_yr : float
        length of timestep in units of years
    prograde_bh_orb ecc : float array
        orbital eccentricity of stars in prograde orbist around SMBH
    disk_star_pro_orbs_ecc_crit : float
        critical eccentricity of stars below which prograde accretion & spin
        will torque into disk alignment else retrograde accretion
    Returns
    -------
    disk_star_spin_angles_new : float array
        spin magnitudes of stars after accreting at prescribed rate for one timestep
    """
    # Calculate change in spin angle due to accretion during timestep
    disk_bh_eddington_ratio_normalized = disk_star_eddington_ratio/1.0
    timestep_duration_yr_normalized = timestep_duration_yr/1.e4
    disk_star_torque_condition_normalized = disk_star_torque_condition/0.1

    spin_torque_iteration = (6.98e-3*disk_bh_eddington_ratio_normalized*disk_star_torque_condition_normalized*timestep_duration_yr_normalized)

    # Assume same angles as before to start
    disk_star_spin_angles_new = prograde_stars_spin_angles
    # Singleton star with orb ecc > disk_star_pro_orbs_ecc_crit will spin down b/c accrete retrograde
    disk_star_pro_spin_down = np.ma.masked_where(disk_star_pro_orbs_ecc <= disk_star_pro_orbs_ecc_crit, disk_star_pro_orbs_ecc)
    # Singleton star with orb ecc < disk_star_pro_orbs_ecc_crit will spin up b/c accrete prograde
    disk_star_pro_spin_up = np.ma.masked_where(disk_star_pro_orbs_ecc > disk_star_pro_orbs_ecc_crit, disk_star_pro_orbs_ecc)
    # Indices of singleton star with orb ecc > disk_star_pro_orbs_ecc_crit
    disk_star_pro_indices_spin_down = np.ma.nonzero(disk_star_pro_spin_down)
    # Indices of singleton star with orb ecc < disk_star_pro_orbs_ecc_crit
    disk_star_pro_indices_spin_up = np.ma.nonzero(disk_star_pro_spin_up)

    # Spin up stars are torqued towards zero (ie alignment with disk,
    # so decrease mag of spin angle)
    disk_star_spin_angles_new[disk_star_pro_indices_spin_up] = prograde_stars_spin_angles[disk_star_pro_indices_spin_up] - spin_torque_iteration
    # Spin down stars with orb ecc > disk_star_pro_orbs_ecc_crit are torqued toward anti-alignment with disk, incr mag of spin angle.
    disk_star_spin_angles_new[disk_star_pro_indices_spin_down] = prograde_stars_spin_angles[disk_star_pro_indices_spin_down] + spin_torque_iteration

    # TO DO: Include a condition to keep spin angle at or close to zero once it gets there
    # Housekeeping
    # Max bh spin angle in rads (pi rads = anti-alignment)
    disk_star_spin_angle_max = 3.10
    disk_star_spin_angles_new[disk_star_spin_angles_new < disk_star_spin_minimum_resolution] = 0.0
    disk_star_spin_angles_new[disk_star_spin_angles_new > disk_star_spin_angle_max] = disk_star_spin_angle_max
    return disk_star_spin_angles_new


def change_bh_spin_magnitudes(disk_bh_pro_spins,
                              disk_bh_eddington_ratio,
                              disk_bh_torque_condition,
                              timestep_duration_yr,
                              disk_bh_pro_orbs_ecc,
                              disk_bh_pro_orbs_ecc_crit):
    """
    Update the spin magnitude of the embedded black holes based on their
    accreted mass in this timestep.

    Parameters
    ----------
    disk_bh_pro_spins : float array
        initial spins of black holes in prograde orbits around SMBH
    disk_bh_eddington_ratio : float
        user chosen input set by input file; Accretion rate of fully
        embedded stellar mass black hole in units of Eddington accretion rate.
        1.0=embedded BH accreting at Eddington.
        Super-Eddington accretion rates are permitted.
    disk_bh_torque_condition : float
        user chosen input set by input file; fraction of initial mass required
        to be accreted before BH spin is torqued fully into alignment with
        the AGN disk. We don't know for sure but Bogdanovic et al. says
        between 0.01=1% and 0.1=10% is what is required
    timestep_duration_yr : float
        length of timestep in units of years
    disk_bh_pro_orbs_ecc : float array
        orbital eccentricity of BH in prograde orbits around SMBH
    disk_bh_pro_orbs_ecc_crit : float
        critical value of orbital eccentricity below which prograde accretion
        (& migration & binary formation) occurs
    Returns
    -------
    disk_bh_pro_spins_new : float array
        spin magnitudes of black holes after accreting at prescribed rate for
        one timestep
    """
    # A retrograde BH a=-1 will spin down to a=0 when it accretes a factor sqrt(3/2)=1.22 in mass (Bardeen 1970).
    # Since M_edd/t = 2.3 e-8 M0/yr or 2.3e-4M0/10kyr then M(t)=M0*exp((M_edd/t)*f_edd*time)
    # so M(t)~1.2=M0*exp(0.2) so in 10^7yr, spin should go a=-1 to a=0. Or delta a ~ 10^-3 every 10^4yr.

    normalized_Eddington_ratio = disk_bh_eddington_ratio/1.0
    normalized_timestep = timestep_duration_yr/1.e4
    normalized_spin_torque_condition = disk_bh_torque_condition/0.1

    # Magnitude of spin iteration per normalized timestep
    spin_iteration = (1.e-3*normalized_Eddington_ratio*normalized_spin_torque_condition*normalized_timestep)

    disk_bh_pro_spins_new = disk_bh_pro_spins
    # Singleton BH with orb ecc > disk_bh_pro_orbs_ecc_crit will spin down b/c accrete retrograde
    prograde_bh_spin_down = np.ma.masked_where(disk_bh_pro_orbs_ecc <= disk_bh_pro_orbs_ecc_crit, disk_bh_pro_orbs_ecc)
    # Singleton BH with orb ecc < disk_bh_pro_orbs_ecc_crit will spin up b/c accrete prograde
    prograde_bh_spin_up = np.ma.masked_where(disk_bh_pro_orbs_ecc > disk_bh_pro_orbs_ecc_crit, disk_bh_pro_orbs_ecc)
    # Indices of singleton BH with orb ecc > disk_bh_pro_orbs_ecc_crit
    indices_bh_spin_down = np.ma.nonzero(prograde_bh_spin_down)
    # Indices of singleton BH with orb ecc < disk_bh_pro_orbs_ecc_crit
    indices_bh_spin_up = np.ma.nonzero(prograde_bh_spin_up)
    # disk_bh_pro_spins_new[prograde_orb_ang_mom_indices]=disk_bh_pro_spins_new[prograde_orb_ang_mom_indices]+(4.4e-3*normalized_Eddington_ratio*normalized_spin_torque_condition*normalized_timestep)
    disk_bh_pro_spins_new[indices_bh_spin_up] = disk_bh_pro_spins[indices_bh_spin_up] + spin_iteration
    # Spin down BH with orb ecc > disk_bh_pro_orbs_ecc_crit
    disk_bh_pro_spins_new[indices_bh_spin_down] = disk_bh_pro_spins[indices_bh_spin_down] - spin_iteration
    # Housekeeping: Max possible spins. Do not spin above or below these values
    disk_bh_pro_spin_max = 0.98
    disk_bh_pro_spin_min = -0.98

    for i in range(len(disk_bh_pro_spins)):
        if (disk_bh_pro_spins_new[i] < disk_bh_pro_spin_min):
            disk_bh_pro_spins_new[i] = disk_bh_pro_spin_min

        if (disk_bh_pro_spins_new[i] > disk_bh_pro_spin_max):
            disk_bh_pro_spins_new[i] = disk_bh_pro_spin_max

    return disk_bh_pro_spins_new


def change_bh_spin_angles(disk_bh_pro_spin_angles,
                          disk_bh_eddington_ratio,
                          disk_bh_torque_condition,
                          disk_bh_spin_minimum_resolution,
                          timestep_duration_yr,
                          disk_bh_pro_orbs_ecc,
                          disk_bh_pro_orbs_ecc_crit):
    """Update the spin magnitude of the embedded black holes based on their accreted mass
        in this timestep.

    Parameters
    ----------
    disk_bh_pro_spin_angles : float array
        initial spin angles of black holes in prograde orbits around SMBH
    disk_bh_eddington_ratio : float
        user chosen input set by input file; Accretion rate of fully embedded stellar mass
        black hole in units of Eddington accretion rate. 1.0=embedded BH accreting at Eddington.
        Super-Eddington accretion rates are permitted.
    disk_bh_torque_condition : float
        user chosen input set by input file; fraction of initial mass required to be
        accreted before BH spin is torqued fully into alignment with the AGN disk.
        We don't know for sure but Bogdanovic et al. says between 0.01=1% and 0.1=10%
        is what is required
    disk_bh_spin_minimum_resolution : float
        user chosen input set by input file; minimum resolution of spin change followed by code.
    timestep_duration_yr : float
        length of timestep in units of years
    disk_bh_orbs_ecc : float array
        orbital eccentricity of BH in prograde orbits around SMBH
    disk_bh_orbs_ecc_crit : float
        critical value of orbital eccentricity below which prograde accretion (& migration & binary formation) occurs
    Returns
    -------
    disk_bh_pro_spin_new : float array
        spin magnitudes of black holes after accreting at prescribed rate for one timestep
    """

    # Calculate change in spin angle due to accretion during timestep
    normalized_Eddington_ratio = disk_bh_eddington_ratio/1.0
    normalized_timestep = timestep_duration_yr/1.e4
    normalized_spin_torque_condition = disk_bh_torque_condition/0.1

    spin_torque_iteration = (6.98e-3*normalized_Eddington_ratio*normalized_spin_torque_condition*normalized_timestep)

    # Assume same angles as before to start
    disk_bh_spin_angles_new = disk_bh_pro_spin_angles
    # Singleton BH with orb ecc > disk_bh_pro_orbs_ecc_crit will spin down b/c accrete retrograde
    prograde_bh_spin_down = np.ma.masked_where(disk_bh_pro_orbs_ecc <= disk_bh_pro_orbs_ecc_crit, disk_bh_pro_orbs_ecc)
    # Singleton BH with orb ecc < disk_bh_pro_orbs_ecc_crit will spin up b/c accrete prograde
    prograde_bh_spin_up = np.ma.masked_where(disk_bh_pro_orbs_ecc > disk_bh_pro_orbs_ecc_crit, disk_bh_pro_orbs_ecc)
    # Indices of singleton BH with orb ecc > disk_bh_pro_orbs_ecc_crit
    indices_bh_spin_down = np.ma.nonzero(prograde_bh_spin_down)
    # Indices of singleton BH with orb ecc < disk_bh_pro_orbs_ecc_crit
    indices_bh_spin_up = np.ma.nonzero(prograde_bh_spin_up)

    # Spin up BH are torqued towards zero (ie alignment with disk, so decrease mag of spin angle)
    disk_bh_spin_angles_new[indices_bh_spin_up] = disk_bh_pro_spin_angles[indices_bh_spin_up] - spin_torque_iteration
    # Spin down BH with orb ecc > disk_bh_pro_orbs_ecc_crit are torqued toward anti-alignment with disk, incr mag of spin angle.
    disk_bh_spin_angles_new[indices_bh_spin_down] = disk_bh_pro_spin_angles[indices_bh_spin_down] + spin_torque_iteration
    # print(disk_bh_spin_angles_new[indices_bh_spin_down])

    # Housekeeping
    # Max bh spin angle in rads (pi rads = anti-alignment). Do not grow bh spin angle < 0 or > bh_max_spin_angle
    bh_max_spin_angle = 3.10
    disk_bh_spin_angles_new[disk_bh_spin_angles_new < disk_bh_spin_minimum_resolution] = 0.0
    disk_bh_spin_angles_new[disk_bh_spin_angles_new > bh_max_spin_angle] = bh_max_spin_angle

    return disk_bh_spin_angles_new