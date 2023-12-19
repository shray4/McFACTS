import numpy as np

def add_to_binary_array2(rng, bin_array, bh_locations, bh_masses, bh_spins, bh_spin_angles, bh_gens, close_encounters, bindex, retro, verbose=False,):
    """This is where we add new binaries. We take the locations, masses, spins, spin angles and generations
    from the relevant singletons, found in hillsphere.binary_check, and sort those parameters into bin_array. 
    We then ADD additional parameters relevant only for binaries, including semi-major axis of the binary,
    semi-major axis of the orbit of the center of mass of the binary around the SMBH, a flag to permit or 
    suppress retrograde binaries, eventually eccentricity and inclination. There is also a verbose flag 
    that is by default set to False, to suppress overly detailed output.

    Parameters
    ----------
    bin_array : [17, bindex] mixed array
        binary black hole array, multi-dimensional;
            [0,j]: float
            location of object 1 at formation of binary (distance from SMBH in r_g)
            [1,j]: float
            location of object 2 at formation of binary (distance from SMBH in r_g)
            [2,j]: float
            mass of obj 1 at time t in units of solar masses
            [3,j]: float
            mass of obj 2 at time t in units of solar masses
            [4,j]: float
            dimensionless spin magnitude of obj 1 at time t
            [5,j]: float
            dimensionless spin magnitude of obj 2 at time t
            [6,j]: float
            spin angle of obj 1 wrt disk gas in radians at time t
            [7,j]: float
            spin angle of obj 2 wrt disk gas in radians at time t
            [8,j]: float
            binary semi-major axis in units of r_g of SMBH
            [9,j]: float
            binary center of mass location wrt SMBH in r_g
            [10,j]: float
            time to merger through GW alone (not set here)
            [11,j]: int
            merger flag = -2 if merging this timestep, else = 0 (not set here)
            [12,j]: float
            time of merger if binary has already merged (not set here)
            [13,j]: ?
            not set here
            [14,j]: int
            generation of obj 1 (1=natal black hole, no prior mergers)
            [15,j]: int
            generation of obj 2 (1=natal black hole, no prior mergers)
            [16,j]: int
            binary angular momentum switch +1/-1 for pro/retrograde
    bh_locations : float array
        locations of prograde singleton BH at start of timestep in units of gravitational radii (r_g=GM_SMBH/c^2)
    bh_masses : float array
        mass of prograde singleton BH at start of timestep in units of solar masses
    bh_spins : float array
        dimensionless spin magnitude of prograde singleton BH at start of timestep
    bh_spin_angles : float array
       spin angle of prograde singleton BH wrt the gas disk angular momentum in radians??? at start of time step
    bh_gens : int array
        generation of prograde singleton BH (1=no previous merger history)
    close_encounters : [2,N] int array
        array of indices corresponding to locations in prograde_bh_locations, prograde_bh_masses,
        prograde_bh_spins, prograde_bh_spin_angles, and prograde_bh_generations which corresponds
        to binaries that form in this timestep. it has a length of the number of binaries to form (N)
        and a width of 2.
    bindex : int
        counter for length of bin_array before adding new binaries
    retro : int
        switch to inhibit retrograde binaries.
        retro = 0 turns all retrograde BBH at formation into prograde BBH. 
        retro = 1 keeps retrograde BBH once they form. Eventually will turn this into
        something physically motivated rather than a brute force switch.
        default : 0
    verbose : bool, optional
        how much output do you want? set True for debugging, by default False

    Returns
    -------
    bin_array : [17, bindex+N] float array
        as for input, but updated (thus longer), to include newly formed binaries
    """
    # find number of new binaries based on indices from hillsphere.binary_check 
    num_new_bins = np.shape(close_encounters)[1]

    # If there are new binaries, actually form them!
    if num_new_bins > 0:
        # send close encounter indices to new array
        array_of_indices = close_encounters
        print("Close encounters ", np.shape(close_encounters)[1], array_of_indices)
        bincount = 0
        # for all the new binaries that need to be created
        for j in range(bindex, bindex + num_new_bins):
            # for each member of the binary
            for i in range(0,2):
                # pick the [0,N] or [1,N] index for member 1 and member 2 of binary
                thing1 = array_of_indices[i][bincount]
                bin_array[i,j] = bh_locations[thing1]
                bin_array[i+2,j] = bh_masses[thing1]
                bin_array[i+4,j] = bh_spins[thing1]
                bin_array[i+6,j] = bh_spin_angles[thing1]
                # For new binary create initial binary semi-major axis
                temp_loc_1 = bin_array[0,j]
                temp_loc_2 = bin_array[1,j]
                temp_bin_separation = np.abs(temp_loc_1 - temp_loc_2)
                bin_array[8,j] = temp_bin_separation
                # Binary c.o.m.= location_1 + separation*M_2/(M_1+M_2)
                temp_mass_1 = bin_array[2,j]
                temp_mass_2 = bin_array[3,j]
                temp_bin_mass = temp_mass_1 + temp_mass_2
                bin_array[9,j] = temp_loc_1 + (temp_bin_separation*temp_mass_2/temp_bin_mass)
                # Set up binary member generations
                bin_array[i+14,j] = bh_gens[thing1]
                # Set up bin orb. ang. mom. (randomly +1 (pro) or -1(retrograde))
                # random number
                random_uniform_number = rng.random()
                bh_initial_orb_ang_mom = (2.0*np.around(random_uniform_number)) - 1.0
                # If retro switch is zero, turn all retro BBH at formation into prograde.
                if retro == 0:
                    bh_initial_orb_ang_mom = np.fabs(bh_initial_orb_ang_mom)
                bin_array[16,j] = bh_initial_orb_ang_mom
                print("Random uniform number =", random_uniform_number )
                print("New orb ang mom =", bh_initial_orb_ang_mom)
            bincount = bincount + 1
        if verbose:
            print("New Binary")
            print(bin_array)
        

    return bin_array
