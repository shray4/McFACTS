"""
Module for calculating the gw strain and freq of a binary.
"""
import numpy as np
from astropy import units as u, constants as const
from astropy.units import cds


def gw_strain_freq(mass_1, mass_2, obj_sep, timestep_duration_yr, old_gw_freq, smbh_mass, agn_redshift, flag_include_old_gw_freq=1):
    """
    This function takes in two masses, their separation, the previous frequency, and the redshift and
    calculates the new GW strain (unitless) and frequency (Hz).

    Note from Saavik about hardcoding strain_factor to 4e3 if nu_gw > 1e-6:
        basically we are implicitly assuming if the frequency is low enough the source is monochromatic
        in LISA over the course of 1yr, so that's where those values come from... and we do need to make
        a decision about that... and that's an ok decision for now. But if someone were to be considering
        a different observatory they might not like that decision?

    Parameters
    ----------
    mass_1 : numpy array
        mass of object 1
    mass_2 : numpy array
        mass of object 2
    obj_sep : numpy array
        separation between both objects, in R_g of the SMBH = GM_smbh/c^2
    timestep_duration_yr : float, or -1 if not given
        current timestep in years
    old_gw_freq : numpy array, or -1 if not given
        previous GW frequency
    smbh_mass : float
        mass of the SMBH in Msun
    agn_redshift : float
        redshift of the SMBH
    flag_include_old_gw_freq : boolean
        flag indicating if old_gw_freq should be included in calculations
        if not, we use the hardcoded value (see note above)
        0 if no, 1 if yes

    Returns
    -------
    char_strain : numpy array
        dimensionless characteristic strain
    nu_gw : numpy array
        GW frequency in Hz
    """

    redshift_d_obs_dict = {0.1: 421*u.Mpc,
                           0.5: 1909*u.Mpc}

    timestep_units = (timestep_duration_yr*u.year).to(u.second)

    # 1rg =1AU=1.5e11m for 1e8Msun
    rg = 1.5e11*(smbh_mass/1.e8)*u.meter
    mass_1 = (mass_1 * cds.Msun).to(u.kg)
    mass_2 = (mass_2 * cds.Msun).to(u.kg)
    mass_total = mass_1 + mass_2
    bin_sep = obj_sep * rg

    mass_chirp = np.power(mass_1 * mass_2, 3./5.) / np.power(mass_total, 1./5.)
    rg_chirp = ((const.G * mass_chirp) / np.power(const.c, 2)).to(u.meter)

    # If separation is less than rg_chirp then cap separation at rg_chirp.
    bin_sep[bin_sep < rg_chirp] = rg_chirp[bin_sep < rg_chirp]

    nu_gw = (1.0/np.pi) * np.sqrt(mass_total * const.G / np.power(bin_sep, 3))
    nu_gw = nu_gw.to(u.Hz)

    # For local distances, approx d=cz/H0 = 3e8m/s(z)/70km/s/Mpc =3.e8 (z)/7e4 Mpc =428 Mpc
    # From Ned Wright's calculator https://www.astro.ucla.edu/~wright/CosmoCalc.html
    # (z=0.1)=421Mpc. (z=0.5)=1909 Mpc
    d_obs = redshift_d_obs_dict[agn_redshift].to(u.meter)
    strain = (4/d_obs) * rg_chirp * np.power(np.pi * nu_gw * rg_chirp / const.c, 2./3.)

    # But power builds up in band over multiple cycles!
    # So characteristic strain amplitude measured by e.g. LISA is given by h_char^2 = N/8*h_0^2 where N is number of cycles per year & divide by 8 to average over viewing angles
    strain_factor = np.ones(len(nu_gw))

    # char amplitude = strain_factor*h0
    #                = sqrt(N/8)*h_0 and N=freq*1yr for approx const. freq. sources over ~~yr.
    strain_factor[nu_gw < (1e-6)*u.Hz] = np.sqrt(nu_gw[nu_gw < (1e-6)*u.Hz]*np.pi*(1e7)/8)

    # For a source changing rapidly over 1 yr, N~freq^2/ (dfreq/dt).
    # char amplitude = strain_factor*h0
    #                = sqrt(freq^2/(dfreq/dt)/8)*h0
    if (flag_include_old_gw_freq == 1):
        delta_nu = np.abs(old_gw_freq - nu_gw)
        delta_nu_delta_timestep = delta_nu/timestep_units
        nu_squared = (nu_gw*nu_gw)
        strain_factor[nu_gw > (1e-6)*u.Hz] = np.sqrt((nu_squared[nu_gw > (1e-6)*u.Hz] / delta_nu_delta_timestep[nu_gw > (1e-6)*u.Hz])/8.)
    # Condition from evolve_gw
    elif (flag_include_old_gw_freq == 0):
        strain_factor[nu_gw > (1e-6)*u.Hz] = np.full(np.sum(nu_gw > (1e-6)*u.Hz), 4.e3)
    char_strain = strain_factor*strain

    return (char_strain.value, nu_gw.value)


def evolve_gw(blackholes_binary, smbh_mass, agn_redshift):

    char_strain, nu_gw = gw_strain_freq(mass_1=blackholes_binary.mass_1,
                                        mass_2=blackholes_binary.mass_2,
                                        obj_sep=blackholes_binary.bin_sep,
                                        timestep_duration_yr=-1,
                                        old_gw_freq=-1,
                                        smbh_mass=smbh_mass,
                                        agn_redshift=agn_redshift,
                                        flag_include_old_gw_freq=0)

    # Update binaries
    blackholes_binary.gw_freq = nu_gw
    blackholes_binary.gw_strain = char_strain

    return (blackholes_binary)


def bbh_gw_params(blackholes_binary, bh_binary_id_num_gw, smbh_mass, timestep_duration_yr, old_bbh_freq, agn_redshift):
    """
    This function evaluates the binary gravitational wave frequency and strain at the end of each timestep_duration_yr.
    Set up binary GW frequency nu_gw = 1/pi *sqrt(GM_bin/a_bin^3). Set up binary strain of GW
    h = (4/d_obs) *(GM_chirp/c^2)*(pi*nu_gw*GM_chirp/c^3)^(2/3)
    where m_chirp =(M_1 M_2)^(3/5) /(M_bin)^(1/5)


    Parameters
    ----------
    blackholes_binary : AGNBinaryBlackHole
        Full binary array.
    bh_binary_id_num_gw : numpy array
        ID numbers of binaries with separations below min_bbh_gw_separation
    smbh_mass : float
        Mass of the SMBH
    timestep_duration_yr : float
        timestep in years
    old_bbh_freq : numpy array
        Previous gw_freq
    agn_redshift : float
        Redshift of the AGN, used to set d_obs

    Returns
    -------
    char_strain : numpy array
        characteristic strain, unitless
    nu_gw : numpy array
        GW frequency in Hz
    """

    num_tracked = bh_binary_id_num_gw.size

    old_bbh_freq = old_bbh_freq * u.Hz

    while (num_tracked > len(old_bbh_freq)):
        old_bbh_freq = np.append(old_bbh_freq, (9.e-7) * u.Hz)

    while (num_tracked < len(old_bbh_freq)):
        old_bbh_freq = np.delete(old_bbh_freq, 0)

    char_strain, nu_gw = gw_strain_freq(mass_1=blackholes_binary.at_id_num(bh_binary_id_num_gw, "mass_1"),
                                        mass_2=blackholes_binary.at_id_num(bh_binary_id_num_gw, "mass_2"),
                                        obj_sep=blackholes_binary.at_id_num(bh_binary_id_num_gw, "bin_sep"),
                                        timestep_duration_yr=timestep_duration_yr,
                                        old_gw_freq=old_bbh_freq,
                                        smbh_mass=smbh_mass,
                                        agn_redshift=agn_redshift,
                                        flag_include_old_gw_freq=1)

    return (char_strain, nu_gw)