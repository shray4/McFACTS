import juliacall
import numpy as np
import fit_modeler
import evolve_binary

if __name__ == "__main__":
    mass_1 = 0.8
    mass_2 = 0.2
    spin_1 = mass_1**2 * 0.7
    spin_2 = mass_2**2 * 0.2
    spin_angle_1 = np.pi / 3
    spin_angle_2 = np.pi / 2
    phi_12 = np.pi / 4
    bin_sep = 10
    bin_orb_inc = [0, 0, 0]

    surrogate = fit_modeler.GPRFitters.read_from_file(f"surrogate.joblib")

    M_f, spin_f, v_f = evolve_binary.evolve_binary(
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
        True,
    )

    print("M_f = ", M_f)
    print("spin_f = ", spin_f)
    print("v_f = ", v_f)
