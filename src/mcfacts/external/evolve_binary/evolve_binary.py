import numpy as np
import quaternionic
from sxs import WaveformModes

import fit_modeler

import juliacall

PostNewtonian = juliacall.newmodule("PN")
PostNewtonian.seval("using PostNewtonian")


def log_q_chi_hat_chi_a(q, chi_1, chi_2):
    """
    Convert q, chi_1, and chi_2 to parameters
    expected by the NR remnant surrogate.

    Parameters
    ----------
    q : float
        mass ratio > 1
    chi_1 : numpy array
        dimensionless inertial spin of heavier object
    chi_2 : numpy array
        dimensionless inertial spin of lighter object
    """
    chi_eff = (q * chi_1 + chi_2) / (1 + q)
    chi_a = 0.5 * (chi_1 - chi_2)

    eta = q / (1 + q) ** 2

    chi_hat = (chi_eff[2] - 38 * eta * (chi_1[2] + chi_2[2]) / 113) / (
        1 - 76 * eta / 113
    )

    return np.log(q), *[chi_eff[0], chi_eff[1], chi_hat], *chi_a


def compute_orbital_omega(mass_1, mass_2, chi_1, chi_2, r, pn_order=3.5):
    """
    Compute the orbital frequency omega using
    Eqs. (373) and (589) of https://arxiv.org/pdf/1310.1528.
    Assumes that z is the direction of orbital angular momentum.
    NS terms >= 3PN order excluded because of UV length scale gauge constant.

    Parameters
    ----------
    mass_1 : float
        mass of heavier object
    mass_2 : float
        mass of lighter object
    chi_1 : numpy array
        dimensionless inertial spin of heavier object
    chi_2 : numpy array
        dimensionless inertial spin of lighter object
    r : float
        binary separation
    pn_order : float
        PN order to use when computing omega.

    Returns
    -------
    omega : float
        orbital frequency
    """
    m = mass_1 + mass_2
    delta = mass_1 - mass_2
    nu = mass_1 * mass_2 / m**2
    gamma = m / r
    S_l = mass_1**2 * chi_1[2] + mass_2**2 * chi_2[2]
    Sigma_l = -(
        mass_1 / m * mass_1**2 * chi_1[2] - mass_2 / m * mass_2**2 * chi_2[2]
    )

    corrections = 1
    if pn_order >= 1:
        corrections += (-3 + nu) * gamma
    if pn_order >= 1.5:
        corrections += gamma ** (3 / 2) / m**2 * (-5 * S_l - 3 * delta / m * Sigma_l)
    if pn_order >= 2:
        corrections += (6 + 41 / 4 * nu + nu**2) * gamma**2
    if pn_order >= 2.5:
        corrections += (
            gamma ** (5 / 2)
            / m**2
            * (
                (45 / 2 - 27 / 2 * nu) * S_l
                + delta / m * (27 / 2 - 13 / 2 * nu) * Sigma_l
            )
            * gamma
        )
    if pn_order >= 3.5:
        corrections += (
            gamma ** (7 / 2)
            / m**2
            * (
                (-495 / 8 - 561 / 8 * nu - 51 / 8 * nu**2) * S_l
                + delta / m * (-297 / 8 - 341 / 8 * nu - 21 / 8 * nu**2) * Sigma_l
            )
            * gamma**2
        )

    omega_sq = m / r**3 * corrections

    omega = np.sqrt(omega_sq)

    return omega


def compute_precessional_omega(mass_1, mass_2, chi_1, chi_2, r, omega, pn_order=5.5):
    """
    Compute the precessional frequency omega using
    Eq. (4.4) of https://arxiv.org/pdf/1212.5520.
    Assumes that z is the direction of orbital angular momentum.

    Parameters
    ----------
    mass_1 : float
        mass of heavier object
    mass_2 : float
        mass of lighter object
    chi_1 : numpy array
        dimensionless inertial spin of heavier object
    chi_2 : numpy array
        dimensionless inertial spin of lighter object
    r : float
        binary separation
    omega : float
        orbital angular momementum
    pn_order : float
        PN order to use when computing omega.

    Returns
    -------
    prec_omega : float
        precessional frequency
    """
    m = mass_1 + mass_2
    delta = mass_1 - mass_2
    nu = mass_1 * mass_2 / m**2
    S_n = mass_1**2 * chi_1[0] + mass_2**2 * chi_2[0]
    Sigma_n = -(
        mass_1 / m * mass_1**2 * chi_1[0] - mass_2 / m * mass_2**2 * chi_2[0]
    )

    x = (m * omega) ** (2 / 3)

    corrections = 0
    if PN_order >= 3.5:
        corrections += 7 * S_n + 3 * delta / m * Sigma_n
    if PN_order >= 4.5:
        corrections += (
            (-10 - 29 / 3 * nu) * S_n + delta / m * (-6 - 9 / 2 * nu) * Sigma_n
        ) * x
    if PN_order >= 5.5:
        corrections += (
            (3 / 2 + 59 / 4 * nu + 52 / 9 * nu**2) * S_n
            + delta / m * (3 / 2 + 73 / 8 * nu + 17 / 6 * nu**2) * Sigma_n
        ) * x**2

    a_l = x ** (7 / 2) / m**3 * corrections

    prec_omega = a_l / (r * omega)

    return prec_omega


def compute_frame_rotation(h, idx):
    """
    Compute the frame rotation that maps the system to have it's
    - maximal emission direction along the z-axis;
    - (2,2) phase = 0;
    - Re[(2,1)] > 0.
    This should correspond to the orbital angular momentum
    in the positive z-direction and the heavier BH on the positve x-axis.

    Parameters
    ----------
    h : sxs.WaveformModes
        strain waveform
    idx : int
        index of data to use to fix the frame

    Returns
    -------
    frame_rotation : quaternionic.array
        rotation to the coorbital frame at h.t[idx].
    """
    # maximal emission reset
    omega = h.angular_velocity[idx]
    omega_as_q = quaternionic.array.from_vector_part(omega).normalized
    maximal_emission_rotation = (
        quaternionic.array(1, 0, 0, 0) - omega_as_q * quaternionic.z
    ).normalized

    # (2,2) phase reset
    phase = (
        -np.angle(h.copy().rotate(maximal_emission_rotation).data[idx, h.index(2, 2)])
        / 2
    )
    phase_rotation = quaternionic.array.from_axis_angle(
        phase * np.array([0, 0, 1])
    ).normalized

    partial_rotation = maximal_emission_rotation * phase_rotation

    # (2,1) real part reset
    pi_freedom_rotation = quaternionic.array.from_axis_angle(
        0.0 * np.array([0, 0, 1])
    ).normalized
    if h.copy().rotate(partial_rotation).data[idx, h.index(2, 1)] < 0:
        pi_freedom_rotation = quaternionic.array.from_axis_angle(
            np.pi * np.array([0, 0, 1])
        ).normalized

    return partial_rotation * pi_freedom_rotation


def evolve_binary_w_PN(
    mass_1, mass_2, chi_1, chi_2, bin_sep, waveform_pn_order=None, verbose=False
):
    """
    Evolve a binary to merger using PN equations.

    Parameters
    ----------
    mass_1 : float
        mass of heavier object
    mass_2 : float
        mass of lighter object
    chi_1 : numpy array
        dimensionless inertial spin of heavier object
    chi_2 : numpy array
        dimensionless inertial spin of lighter object
    bin_sep : float
        binary separation
    waveform_pn_order : float
        PN order to use when computing omega
    verbose : bool
        Whether or not to print useful information

    Returns
    -------
    mass_1 : float
        mass of heavier object 100M before merger
    mass_2 : float
        mass of lighter object 100M before merger
    chi_1 : numpy array
        dimensionless inertial spin of heavier object 100M before merger
    chi_2 : numpy array
        dimensionless inertial spin of lighter object 100M before merger
    frame_rotation : quaternionic.array
        rotation to the coorbital frame 100M before merger
    """
    if verbose:
        print("Evolving binary to merger using PN equations...")

    omega_i = compute_orbital_omega(mass_1, mass_2, chi_1, chi_2, bin_sep)

    inspiral = PostNewtonian.orbital_evolution(mass_1, mass_2, chi_1, chi_2, omega_i)
    values = PostNewtonian.stack(inspiral.u)

    waveform_pn_order = waveform_pn_order or PostNewtonian.typemax(PostNewtonian.Int)
    w = (
        PostNewtonian.inertial_waveform(
            inspiral, ell_min=2, ell_max=4, PNOrder=waveform_pn_order
        )
        .to_numpy()
        .T
    )

    h = WaveformModes(
        w, time=inspiral.t, modes_axis=1, ell_min=2, ell_max=4, spin_weight=-2
    )

    idx_surrogate_matching = np.argmin(abs(h.t - (h.t[-1] - 100)))

    mass_1, mass_2, chi_1, chi_2 = (
        values[0, :].to_numpy()[idx_surrogate_matching],
        values[1, :].to_numpy()[idx_surrogate_matching],
        values[2:5, :].to_numpy().T[idx_surrogate_matching],
        values[5:8, :].to_numpy().T[idx_surrogate_matching],
    )

    frame_rotation = compute_frame_rotation(h, idx_surrogate_matching)

    return mass_1, mass_2, chi_1, chi_2, frame_rotation


def evolve_binary_w_surrogate(mass_1, mass_2, chi_1, chi_2, surrogate, verbose=False):
    """
    Evolve a binary through merger using a NR surrogate.
    Note that the frame expected by the surrogate is the coorbital frame, i.e.,
    the orbital angular momentum is in the positive z-axis and the
    heavier black hole is on the positive x-axis. Consequently,
    the remnant parameters are returned in this frame.

    Parameters
    ----------
    mass_1 : float
        mass of heavier object
    mass_2 : float
        mass of lighter object
    chi_1 : numpy array
        dimensionless inertial spin of heavier object
    chi_2 : numpy array
        dimensionless inertial spin of lighter object
    surrogate : fit_modeler.GPRFitters
        surrogate model of remnant parameters
    verbose : bool
        Whether or not to print useful information

    Returns
    -------
    M_f : float
        remnant mass
    chi_f : numpy array
        remnant dimensionless inertial spin
    v_f : numpy array
        remnant velocity
    """
    if verbose:
        print("Evolving binary through merger with surrogate...")

    q = mass_1 / mass_2

    input_values_for_surrogate = log_q_chi_hat_chi_a(q, chi_1, chi_2)
    remnant_parameters = surrogate(input_values_for_surrogate)

    M_f = (mass_1 + mass_2) * remnant_parameters[0]
    chi_f = remnant_parameters[1:4]
    v_f = remnant_parameters[4:7]

    return M_f, chi_f, v_f


def evolve_binary(
    mass_1,
    mass_2,
    spin_1,
    spin_2,
    spin_angle_1,
    spin_angle_2,
    phi_12,
    bin_sep,
    bin_orb_inc,
    surrogate,
    verbose=False,
):
    """
    Evolve a binary through merger.

    Parameters
    ----------
    mass_1 : float
        mass of heavier object
    mass_2 : float
        mass of lighter object
    spin_1 : float
        dimensionful inertial spin magnitude of heavier object
    spin_2 : float
        dimensionful inertial spin magnitude of lighter object
    spin_1_angle : float
        angle between spin of heavier object and positive z-axis
    spin_2_angle : float
        angle between spin of lighter object and positive z-axis
    phi_12 : float
        angle between spin vectors in the x/y-plane
    bin_sep : float
        binary separation
    bin_orb_inc : numpy array
        orbital inclination of the binary wrt SMBH
    surrogate : fit_modeler.GPRFitters
        surrogate model of remnant parameters
    verbose : bool
        Whether or not to print useful information

    Returns
    -------
    M_f : float
        remnant mass
    spin_f : numpy array
        remnant dimensionful inertial spin
    v_f : numpy array
        remnant velocity
    """
    if verbose:
        print("Evolving binary...")

    # Compute dimensionless spin parameters
    chi_1 = (
        spin_1
        / mass_1**2
        * np.array(
            [
                np.cos(phi_12) * np.sin(np.pi / 2 - spin_angle_1),
                np.sin(phi_12) * np.sin(np.pi / 2 - spin_angle_1),
                np.cos(np.pi / 2 - spin_angle_1),
            ]
        )
    )
    chi_2 = (
        spin_2
        / mass_2**2
        * np.array(
            [
                np.cos(phi_12) * np.sin(np.pi / 2 - spin_angle_2),
                np.sin(phi_12) * np.sin(np.pi / 2 - spin_angle_2),
                np.cos(np.pi / 2 - spin_angle_2),
            ]
        )
    )

    # Flip labels, if necessary
    if mass_1 < mass_2:
        if verbose:
            print("mass_1 < mass_2! Flipping black hole labels.")

        mass_1_copy = mass_1
        mass_1 = mass_2
        mass_2 = mass_1_copy

        chi_1_copy = chi_1
        chi_1 = chi_2
        chi_2 = chi_1_copy

    # Map to the frame expected by PN;
    # see https://moble.github.io/PostNewtonian.jl/dev/internals/fundamental_variables/#PostNewtonian.FundamentalVariables.R,
    # i.e., Newtonian orbital angular momentum in the positive z-direction
    # and the heavier black hole on the positive x-axis

    # Evolve binary to merger using PN equations
    mass_1, mass_2, chi_1, chi_2, frame_rotation = evolve_binary_w_PN(
        mass_1, mass_2, chi_1, chi_2, bin_sep, verbose
    )

    # Map to the frame expected by the NR surrogate
    if np.linalg.norm(chi_1) > 1e-4:
        chi_1_prime = frame_rotation.rotate(chi_1)
    else:
        chi_1_prime = chi_1

    if np.linalg.norm(chi_1) > 1e-4:
        chi_2_prime = frame_rotation.rotate(chi_2)
    else:
        chi_2_prime = chi_2

    # Evolve binary through merger using the NR surrogate
    M_f, chi_f_prime, v_f_prime = evolve_binary_w_surrogate(
        mass_1, mass_2, chi_1_prime, chi_2_prime, surrogate, verbose
    )

    # Map back to the frame of the binary
    spin_f = M_f**2 * frame_rotation.inverse.rotate(chi_f_prime)
    v_f = frame_rotation.inverse.rotate(v_f_prime)

    return M_f, spin_f, v_f
