def ReadInputs():
    """This function reads your input choices from a file named
    model_choice.txt, and returns the chosen variables for manipulation
    by main.    

    Required input formats and units are given in the model_choice.txt file.

    See below for full output list, including units & formats

    Example
    -------
    To run, ensure a model_choice.txt is in the same directory and type:

        $ python ReadInputs.py

    Notes
    -----
        This should (will) be updated to allow a user to specify an arbitrary
        input file name.

    Function will tell you what it is doing via print statements along the way.

    Attributes
    ----------
    Output variables:
    mass_smbh : float
        Mass of the supermassive black hole (M_sun)
    trap_radius : float
        Radius of migration trap in gravitational radii (r_g = G*mass_smbh/c^2)
        Should be set to zero if disk model has no trap
    n_bh : int
        Number of stellar mass black holes embedded in disk 
    mode_mbh_init : float
        Initial mass distribution for stellar bh is assumed to be Pareto
        with high mass cutoff--mode of initial mass dist (M_sun)
    max_initial_bh_mass : float
        Initial mass distribution for stellar bh is assumed to be Pareto
        with high mass cutoff--mass of cutoff (M_sun)
    mbh_powerlaw_index : float
        Initial mass distribution for stellar bh is assumed to be Pareto
        with high mass cutoff--powerlaw index for Pareto dist
    mu_spin_distribution : float
        Initial spin distribution for stellar bh is assumed to be Gaussian
        --mean of spin dist
    sigma_spin_distribution : float
        Initial spin distribution for stellar bh is assumed to be Gaussian
        --standard deviation of spin dist
    spin_torque_condition : float
        fraction of initial mass required to be accreted before BH spin is 
        torqued fully into alignment with the AGN disk. We don't know for sure but 
        Bogdanovic et al. says between 0.01=1% and 0.1=10% is what is required.
    frac_Eddington_ratio : float
        assumed accretion rate onto stellar bh from disk gas, in units of Eddington
        accretion rate
    max_initial_eccentricity : float
        assuming initially flat eccentricity distribution among single orbiters around SMBH
        out to max_initial_eccentricity. Eventually this will become smarter.
    timestep : float
        How long is your timestep in years?
    number_of_timesteps : int
        How many timesteps are you taking (timestep*number_of_timesteps = disk_lifetime)
    disk_model_radius_array : float array
        The radii along which your disk model is defined in units of r_g (=G*mass_smbh/c^2)
        drawn from modelname_surface_density.txt
    disk_inner_radius : float
        0th element of disk_model_radius_array (units of r_g)
    disk_outer_radius : float
        final element of disk_model_radius_array (units of r_g)
    surface_density_array : float array
        Surface density corresponding to radii in disk_model_radius_array (units of kg/m^2)
        Yes, it's in SI not cgs. Get over it. Kisses.
        drawn from modelname_surface_density.txt
    aspect_ratio_array : float array
        Aspect ratio corresponding to radii in disk_model_radius_array
        drawn from modelname_aspect_ratio.txt

    """

    # ReadInputs.py
    # This reads and parses the inputs from model_choice.txt
    # should eventually change to specify infilename generically
    # Even assuming this works, need to pass info out...
    # also note: some items in B's original inputs.txt are NOT INPUTS

    # create a dictionary to store numerical variables in
    input_variables = {}

    # open the main input file for reading
    model_inputs = open("model_choice.txt", 'r')

    # go through the file line by line
    for line in model_inputs:
        line = line.strip()
        # If it is NOT a comment line
        if (line.startswith('#') == 0):
            # split the line between the variable name and its value
            varname, varvalue = line.split("=")
            # remove whitespace
            varname = varname.strip()
            varvalue = varvalue.strip()
            # if the variable is the one that's a filename for the disk model, deal with it
            if (varname == 'disk_model_name'):
                disk_model_name = varvalue.strip("'")
            # or if it's the number of black holes, typecast to int
            elif (varname == 'n_bh'):
                input_variables[varname] = int(varvalue)
            # or the number of timesteps, typecast to int
            elif (varname == 'number_of_timesteps'):
                input_variables[varname] = int(varvalue)
            # otherwise, typecast to float and stick it in the dictionary
            else:
                input_variables[varname] = float(varvalue)
 
    # close the file
    model_inputs.close()

    print("I read the main input file and closed it")

    # Recast the inputs from the dictionary lookup to actual variable names
    #   !!!is there a better (automated?) way to do this?
    mass_smbh = input_variables['mass_smbh']
    trap_radius = input_variables['trap_radius']
    n_bh = input_variables['n_bh']
    mode_mbh_init = input_variables['mode_mbh_init']
    max_initial_bh_mass = input_variables['max_initial_bh_mass']
    mbh_powerlaw_index = input_variables['mbh_powerlaw_index']
    mu_spin_distribution = input_variables['mu_spin_distribution']
    sigma_spin_distribution = input_variables['sigma_spin_distribution']
    spin_torque_condition = input_variables['spin_torque_condition']
    frac_Eddington_ratio = input_variables['frac_Eddington_ratio']
    max_initial_eccentricity = input_variables['max_initial_eccentricity']
    timestep = input_variables['timestep']
    number_of_timesteps = input_variables['number_of_timesteps']

    print("I put your variables where they belong")

    # open the disk model surface density file and read it in
    # Note format is assumed to be comments with #
    #   density in SI in first column
    #   radius in r_g in second column
    #   infile = model_surface_density.txt, where model is user choice
    infile_suffix = '_surface_density.txt'
    infile = disk_model_name+infile_suffix
    surface_density_file = open(infile, 'r')
    density_list = []
    radius_list = []
    for line in surface_density_file:
        line = line.strip()
        # If it is NOT a comment line
        if (line.startswith('#') == 0):
            columns = line.split()
            density_list.append(float(columns[0]))
            radius_list.append(float(columns[1]))
    # close file
    surface_density_file.close()

    # re-cast from lists to arrays
    surface_density_array = np.array(density_list)
    disk_model_radius_array = np.array(radius_list)

    # open the disk model aspect ratio file and read it in
    # Note format is assumed to be comments with #
    #   aspect ratio in first column
    #   radius in r_g in second column must be identical to surface density file
    #       (radius is actually ignored in this file!)
    #   filename = model_aspect_ratio.txt, where model is user choice
    infile_suffix = '_aspect_ratio.txt'
    infile = disk_model_name+infile_suffix
    aspect_ratio_file = open(infile, 'r')
    aspect_ratio_list = []
    for line in aspect_ratio_file:
        line = line.strip()
        # If it is NOT a comment line
        if (line.startswith('#') == 0):
            columns = line.split()
            aspect_ratio_list.append(float(columns[0]))
    # close file
    aspect_ratio_file.close()

    # re-cast from lists to arrays
    aspect_ratio_array = np.array(aspect_ratio_list)

    # Housekeeping from input variables
    disk_outer_radius = disk_model_radius_array[-1]
    disk_inner_radius = disk_model_radius_array[0]

    print("I read and digested your disk model")

    print("Sending variables back")

    return mass_smbh, trap_radius, n_bh, mode_mbh_init, max_initial_bh_mass, \
        mbh_powerlaw_index, mu_spin_distribution, sigma_spin_distribution, \
            spin_torque_condition, frac_Eddington_ratio, max_initial_eccentricity, \
                timestep, number_of_timesteps, disk_model_radius_array, disk_inner_radius,\
                    disk_outer_radius, surface_density_array, aspect_ratio_array
