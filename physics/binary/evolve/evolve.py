import numpy as np
import scipy


def change_bin_mass(bin_array, frac_Eddington_ratio, mass_growth_Edd_rate, timestep, integer_nbinprop, bin_index):
    """_summary_

    Parameters
    ----------
    bin_array : _type_
        _description_
    frac_Eddington_ratio : _type_
        _description_
    mass_growth_Edd_rate : _type_
        _description_
    timestep : _type_
        _description_
    integer_nbinprop : _type_
        _description_
    bin_index : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    #Return new updated mass array due to accretion for prograde orbiting BH after timestep
    #Extract the binary locations and masses
    bindex = int(bin_index)
    # Run over active binaries (j is jth binary; i is the ith property of the jth binary, e.g. mass1,mass 2 etc)
    
    for j in range(0, bindex): 
            if bin_array[11,j] < 0:
                #do nothing -merger happened!
                print()
            else:            
                #for i in range(0, integer_nbinprop):
                temp_bh_mass_1 = bin_array[2,j] 
                temp_bh_mass_2 = bin_array[3,j]
                mass_growth_factor = np.exp(mass_growth_Edd_rate*frac_Eddington_ratio*timestep)
                new_bh_mass_1 = temp_bh_mass_1*mass_growth_factor
                new_bh_mass_2 = temp_bh_mass_2*mass_growth_factor
                #Update new bh masses in bin_array
                bin_array[2,j] = new_bh_mass_1
                bin_array[3,j] = new_bh_mass_2

    return bin_array


def change_bin_spin_magnitudes(bin_array, frac_Eddington_ratio, spin_torque_condition, timestep, integer_nbinprop, bin_index):
    #def change_spin_magnitudes(bh_spins,prograde_orb_ang_mom_indices,frac_Eddington_ratio,spin_torque_condition,mass_growth_Edd_rate,timestep):
    #bh_new_spins=bh_spins
    normalized_Eddington_ratio = frac_Eddington_ratio/1.0
    normalized_timestep = timestep/1.e4
    normalized_spin_torque_condition = spin_torque_condition/0.1
    #Extract the binary locations and spin magnitudes
    #max allowed spin
    max_allowed_spin=0.98
    bindex = int(bin_index)
    # Run over active binaries (j is jth binary; i is the ith property of the jth binary, e.g. mass1,mass 2 etc)
    
    for j in range(0, bindex):
            if bin_array[11,j] < 0:
                #do nothing -merger happened!
                print()
            else:
                #for i in range(0, integer_nbinprop):
                temp_bh_spin_1 = bin_array[4,j] 
                temp_bh_spin_2 = bin_array[5,j]
                spin_change_factor = 4.4e-3*normalized_Eddington_ratio*normalized_spin_torque_condition*normalized_timestep
                #print("Spin change factor", spin_change_factor)
                new_bh_spin_1 = temp_bh_spin_1 + spin_change_factor
                #print("Old spin1, new spin1 =",temp_bh_spin_1, new_bh_spin_1)
                new_bh_spin_2 = temp_bh_spin_2 + spin_change_factor
                if new_bh_spin_1 > max_allowed_spin:
                    new_bh_spin_1 = max_allowed_spin
                if new_bh_spin_2 > max_allowed_spin:
                    new_bh_spin_2 = max_allowed_spin
                #Update new bh masses in bin_array
                bin_array[4,j] = new_bh_spin_1
                bin_array[5,j] = new_bh_spin_2

    return bin_array


def change_bin_spin_angles(bin_array, frac_Eddington_ratio, spin_torque_condition, spin_minimum_resolution, timestep, integer_nbinprop, bin_index):
    #Calculate change in spin angle due to accretion during timestep
    normalized_Eddington_ratio = frac_Eddington_ratio/1.0
    normalized_timestep = timestep/1.e4
    normalized_spin_torque_condition = spin_torque_condition/0.1

    #Extract the binary locations and spin magnitudes
    bindex = int(bin_index)
    # Run over active binaries (j is jth binary; i is the ith property of the jth binary, e.g. mass1,mass 2 etc)
    
    for j in range(0, bindex):
            if bin_array[11,j] < 0:
                #do nothing -merger happened!
                print()
            else:
            #for i in range(0, integer_nbinprop):
                temp_bh_spin_angle_1 = bin_array[6,j] 
                temp_bh_spin_angle_2 = bin_array[7,j]
                #bh_new_spin_angles[prograde_orb_ang_mom_indices]=bh_new_spin_angles[prograde_orb_ang_mom_indices]-(6.98e-3*normalized_Eddington_ratio*normalized_spin_torque_condition*normalized_timestep)
                spin_angle_change_factor = (6.98e-3*normalized_Eddington_ratio*normalized_spin_torque_condition*normalized_timestep)
                new_bh_spin_angle_1 = temp_bh_spin_angle_1 - spin_angle_change_factor
                new_bh_spin_angle_2 = temp_bh_spin_angle_2 - spin_angle_change_factor
                if new_bh_spin_angle_1 < spin_minimum_resolution:
                    new_bh_spin_angle_1 = 0.0
                if new_bh_spin_angle_2 < spin_minimum_resolution:
                    new_bh_spin_angle_2 = 0.0
                bin_array[6,j] = new_bh_spin_angle_1
                bin_array[7,j] = new_bh_spin_angle_2

    return bin_array

def com_feedback_hankla(bin_array, disk_surf_model, frac_Eddington_ratio, alpha):
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
    
    bin_array : float array
        binary array. Row 9 is center of mass of binary BH at start of timestep in units of gravitational radii (r_g=GM_SMBH/c^2)
    disk_surf_model : function
        returns AGN gas disk surface density in kg/m^2 given a distance from the SMBH in r_g
        can accept a simple float (constant), but this is deprecated

    Returns
    -------
    ratio_feedback_to_mig : float array
        ratio of feedback torque to migration torque for each entry in prograde_bh_locations
    """
    #Extract the binary locations and masses
    #bindex = int(bin_index)
    # Run over active binaries (j is jth binary; i is the ith property of the jth binary, e.g. mass1,mass 2 etc)
    

    temp_bin_com_locations = bin_array[9,:]
    
    print("Bin com locations",temp_bin_com_locations)

    # get surface density function, or deal with it if only a float
    if isinstance(disk_surf_model, float):
        disk_surface_density = disk_surf_model
    else:
        disk_surface_density = disk_surf_model(temp_bin_com_locations)

    #Calculate ratio
    #
    #Define kappa (or set up a function to call). 
    #kappa = 10^0.76 cm^2/g = 10^(0.76) (10^-2m)^2/10^-3kg=10^(0.76-1)=10^(-0.24) m^2/kg to match units of Sigma
    kappa = 10**(-0.24)
    #Define alpha parameter for disk in Readinputs.py
    #alpha = 0.01

    Ratio_feedback_migration_torque_bin_com = 0.07 *(1/kappa)* ((alpha)**(-1.5))*frac_Eddington_ratio*np.sqrt(temp_bin_com_locations)/disk_surface_density

    print("ratio", Ratio_feedback_migration_torque_bin_com)
    #print((1/kappa),((alpha)**(-1.5)),frac_Eddington_ratio)
    #print("Ratio", Ratio_feedback_migration_torque) 
    #print("BH locations", prograde_bh_locations) 

    return Ratio_feedback_migration_torque_bin_com  


def com_migration(bin_array, disk_surf_model, disk_aspect_ratio_model, timestep, integer_nbinprop, bin_index):
    """_summary_

    Parameters
    ----------
    bin_array : _type_
        _description_
    disk_surf_model : _type_
        _description_
    disk_aspect_ratio_model : _type_
        _description_
    timestep : _type_
        _description_
    integer_nbinprop : _type_
        _description_
    bin_index : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    # !!!This should be re-written to do computations NOT in the loop 
    # (especially the surface density & aspect ratio)
    # Saavik has some worries about the functional forms working correctly
    # because it didn't actually fail when I didn't update them (I did singleton
    # migration first, tested, then did binary evolution, but it didn't fail
    # the first time... so that's a bit worrisome from a physics perspective)
    #Return updated locations of binary center of mass (com) and location 1,2 
    # based on Type 1 migration prescription
    #sg_norm is a normalization factor for the Sirko & Goodman (2003) disk model
    #38Myrs=3.8e7yrs is the time for a 5Msun BH to undergo Type I migration to
    #the SMBH from 10^4r_g in that model.
    sg_norm = 3.8e7
    #scaled mass= BH mass/lower bound mass (e.g. 5Msun, upper end of lower mass gap)
    scaled_mass = 5.0
    #scaled_aspect=disk_aspect ratio scaled to 0.02 as a fiducial value.
    scaled_aspect = 0.02
    #for test fixed disk aspect ratio
    #disk_aspect_ratio = 0.03
    #scaled location= BH location scaled to 10^4r_g
    scaled_location = 1.e4
    #scaled sigma= Disk surface density scaled to 10^5kg/m^2
    scaled_sigma = 1.e5
    
    #Extract the binary locations and masses
    bindex = int(bin_index)
    # Run over active binaries (j is jth binary; i is the ith property of the jth binary, e.g. mass1,mass 2 etc)
    
    #If binary has been removed, it's all zeroes. So find first non-zero binary in array.
    #loc1_bins=bin_array[:,0]
    #live_loc1_bins=np.count_nonzero(loc1_bins) 
    #if live_loc1_bins > 0:
    #    bin_indices = np.where(live_loc1_bins > 0.0)

    for j in range(0,bindex):
            if bin_array[11,j] < 0:
                #do nothing -merger happened!
                print()
            else:
            #for i in range(0, integer_nbinprop):
                temp_bh_loc_1 = bin_array[0,j]
                temp_bh_loc_2 = bin_array[1,j]
                temp_bh_mass_1 = bin_array[2,j] 
                temp_bh_mass_2 = bin_array[3,j]
                temp_bin_sep = bin_array[8,j]
                temp_bin_com = bin_array[9,j]
                bin_mass = temp_bh_mass_1 + temp_bh_mass_2
                #Normalize the com location and BH masses 
                if isinstance(disk_surf_model, float):
                    disk_surface_density = disk_surf_model
                else:
                    disk_surface_density = disk_surf_model(temp_bin_com)
                if isinstance(disk_aspect_ratio_model, float):
                    disk_aspect_ratio = disk_aspect_ratio_model
                else:
                    disk_aspect_ratio = disk_aspect_ratio_model(temp_bin_com)
                normalized_com = temp_bin_com/scaled_location
                normalized_bin_mass = bin_mass/scaled_mass
                normalized_com_sqrt = np.sqrt(normalized_com)
                # Can normalize the aspect ratio and sigma to these scales when we
                # implement the 1d disk model (interpolate over SG03)
                normalized_sigma = disk_surface_density/scaled_sigma
                normalized_aspect_ratio = disk_aspect_ratio/scaled_aspect
                normalized_aspect_ratio_squared = np.square(normalized_aspect_ratio)
                #So our fiducial timescale should now be 38Myrs as calcd below
                dt_mig = sg_norm*(normalized_aspect_ratio_squared)/((normalized_bin_mass)*(normalized_com_sqrt)*(normalized_sigma))
                #Effective fractional time of migration is timestep/dt_mig
                fractional_migration_timestep = timestep/dt_mig
                #Migration distance is location of bin com * fractional_migration_timestep
                migration_distance = temp_bin_com*fractional_migration_timestep
                new_bin_com = temp_bin_com-migration_distance
                new_bh_loc1 = new_bin_com -(temp_bin_sep*temp_bh_mass_2/bin_mass)
                new_bh_loc2 = new_bh_loc1 + temp_bin_sep
                #Write new values of R1,R2,com to bin_array
                bin_array[0,j] = new_bh_loc1
                bin_array[1,j] = new_bh_loc2
                bin_array[9,j] = new_bin_com

    return bin_array

def bin_migration(mass_smbh, bin_array, disk_surf_model, disk_aspect_ratio_model, timestep):
    """This function calculates how far the center of mass of a binary migrates in an AGN gas disk in a time
    of length timestep, assuming a gas disk surface density and aspect ratio profile, for
    objects of specified masses and starting locations, and returns their new locations
    after migration over one timestep.

    This function replaces com_migration which did not include smbh mass and was unreadable.

    Parameters
    ----------
    mass_smbh : float
        mass of supermassive black hole in units of solar masses
    bin_array : 2d float array (?)
        Bane of my existence, giant pain in the ass. All the binary parameters hacked into one.
    disk_surf_model : function
        returns AGN gas disk surface density in kg/m^2 given a distance from the SMBH in r_g
        can accept a simple float (constant), but this is deprecated
    disk_aspect_ratio_model : function
        returns AGN gas disk aspect ratio given a distance from the SMBH in r_g
        can accept a simple float (constant), but this is deprecated
    timestep : float
        size of timestep in years


    Returns
    -------
    bin_array : 2d float array (?)
        Bane of my existence, giant pain in the ass. All the binary parameters hacked into one.
    """
    
    # get locations of center of mass of binary from bin_array
    bin_com = bin_array[9,:]
    # get masses of each binary by adding their component masses
    bin_mass = bin_array[2,:] + bin_array[3,:]
    # get surface density function, or deal with it if only a float
    if isinstance(disk_surf_model, float):
        disk_surface_density = disk_surf_model
    else:
        disk_surface_density = disk_surf_model(bin_com)
    # ditto for aspect ratio
    if isinstance(disk_aspect_ratio_model, float):
        disk_aspect_ratio = disk_aspect_ratio_model
    else:
        disk_aspect_ratio = disk_aspect_ratio_model(bin_com) 

    # compute migration timescale for each binary in seconds
    # eqn from Paardekooper 2014, rewritten for R in terms of r_g of SMBH = GM_SMBH/c^2
    # tau = (pi/2) h^2/(q_d*q) * (1/Omega)
    #   where h is aspect ratio, q is m/M_SMBH, q_d = pi R^2 disk_surface_density/M_SMBH
    #   and Omega is the Keplerian orbital frequency around the SMBH
    # here mass_smbh/prograde_bh_masses are both in M_sun, so units cancel
    # c, G and disk_surface_density in SI units
    tau_mig = ((disk_aspect_ratio**2)* scipy.constants.c/(2.0*scipy.constants.G) * (mass_smbh/bin_mass) / disk_surface_density) / np.sqrt(bin_com)
    # ratio of timestep to tau_mig (timestep in years so convert)
    dt = timestep * scipy.constants.year / tau_mig
    # migration distance is original locations times fraction of tau_mig elapsed
    migration_distance = bin_com * dt
    # new locations are original ones - distance traveled
    bh_new_locations = bin_com - migration_distance
    # send locations back to bin_array and DONE!
    bin_array[9,:] = bh_new_locations

    return bin_array
