import numpy as np

#NS hill sphere function. Is currently copy-pasted from the BH function, with only changes to variable names.
#Given that its only use is in the NS-specific binary formation function,
#it can probably be deleted as soon as the combined binary function is gotten working.

def calculate_hill_sphere(prograde_ns_locations, prograde_ns_masses, mass_smbh):
    #Return the Hill sphere radius (R_Hill) for an array of prograde NS where
    # R_Hill=a(q/3)^1/3 where a=semi-major axis, q=m_ns/M_SMBH
    
    ns_smbh_mass_ratio = prograde_ns_masses/(3.0*mass_smbh)
    mass_ratio_factor = (ns_smbh_mass_ratio)**(1./3.)
    ns_hill_sphere = prograde_ns_locations*mass_ratio_factor
    #Return the NS Hill sphere radii for all orbiters. Prograde should have much larger Hill sphere
    return ns_hill_sphere

def encounter_test(prograde_ns_locations, ns_hill_sphere):
    #Using Hill sphere size and NS locations see if there are encounters within the Hill sphere
    # return indices of NS involved.

    # First sort the prograde ns locations in order from inner disk to outer disk
    sorted_ns_locations = np.sort(prograde_ns_locations)
    #Returns the indices of the original array in order, to get the sorted array
    sorted_ns_location_indices = np.argsort(prograde_ns_locations)
   
    #Find the appropriate (sorted) Hill sphere radii
    sorted_hill_spheres = ns_hill_sphere[sorted_ns_location_indices]

    # Find the distances between [r1,r2,r3,r4,..] as [r2-r1,r3-r2,r4-r3,..]=[delta1,delta2,delta3..]
    separations = np.diff(sorted_ns_locations)

    # Note that separations are -1 of length of ns_locations
    # Take 1st location off locations array
    sorted_ns_locations_minus_first_element = sorted_ns_locations[1:len(sorted_ns_locations)]
    #Take last location off locations array
    sorted_ns_locations_minus_last_element = sorted_ns_locations[0:len(sorted_ns_locations)-1]

    # Separations are -1 of length of Hill_sphere array
    # Take 1st Hill sphere off Hill sphere array
    sorted_hill_spheres_minus_first = sorted_hill_spheres[1:len(sorted_hill_spheres)]
    # Take last Hill sphere off Hill sphere array
    sorted_hill_spheres_minus_last = sorted_hill_spheres[0:len(sorted_hill_spheres)-1]

    # Compare the Hill sphere distance for each NS with separation to neighbor NS
    #so we compare e.g. r2-r1 vs R_H1, r3-r2 vs R_H2
    comparison_distance_inwards = separations-sorted_hill_spheres_minus_last
    # and e.g. compare r2-r1 vs R_H2, r3-r2 vs R_H3
    comparison_distance_outwards = separations-sorted_hill_spheres_minus_first

    index_in = np.where(comparison_distance_inwards < 0)
    if isinstance(index_in,tuple):
        index_in = index_in[0]
    # E.g say r3-r2 <R_H2 then we'll want the info for the NS at r2 & r3. (i,i+1)
    index_out = np.where(comparison_distance_outwards < 0)
    if isinstance(index_out,tuple):
        index_out = index_out[0]
    #E.g. say r3-r2 <R_H3 then we'll want the info for the NS at r3 and r2 (i,i-1)
    length_index_in = len(index_in)
    length_index_out = len(index_out)

    new_indx_in = list(range(2*len(index_in)))
    new_indx_out = list(range(2*len(index_out)))
    #new_indx_in_bin_array = np.array[2*len(index_in),2]
    #temp_new_bin_array = np.ndarray[2*length_index_in,2]

    for ind in range(length_index_in):
        temp_index = index_in[ind]
        new_indx_in[2*ind] = temp_index
        new_indx_in[(2*ind)+1] = temp_index+1
    #    temp_new_bin_array = 
    #print("new_indx_in",new_indx_in)
    temp_new_bin_in_array=np.reshape(new_indx_in,(len(index_in),2))
    #print("ordered as bins",temp_new_bin_in_array)
    
    test_remove_bin = new_indx_in
    
    # Dynamics Here! Potential Double-binary or triple interaction!
    # For now! 
    # 0. Construct array of pairs of indices.
    # 1. Select binaries based on distance
    #     For binaries [i-1,i], [i,i+1] compare distance r(i)-r(i-1) to r(i+1)-r(i) and select smallest.
    #     Remove larger distane pair. E.g. if [i+1,i] is smaller, remove [i-1,i]
    #    But want to calculate 
    # 2. Fractional R_Hill for [i-1,i] vs [i,i+1]. Smaller fractional Hill radius wins
    # 3. TO DO Write a module for DiLaurentii+22 or Rowan+22 or LANL+22 phase space encounter and apply to all encounters 
    #          over timestep (10kyrs; assume random phase & number of encounters during timestep; pick randomly 
    #          from phase plots.) Also look at LANL group papers on binding energy of encounter.
    # 4. Ideally, consider the triple-dynamics encounter 
    
    # Search for repeat indices in binary array 
    unique_element,unique_index,unique_ct = np.unique(temp_new_bin_in_array,return_inverse = True,return_counts = True)
    repeats = unique_ct > 1
    
    repeat_indices = unique_element[repeats]
    
    duplicate_indices,dup_index,dup_ct = np.unique(unique_index, return_inverse = True, return_counts = True)
    rep_idx = dup_ct >1
    dupl_indices = duplicate_indices[rep_idx]
    #print("indices of repeats",dupl_indices)
    

    #return all those elements that only occur once
    new_array = unique_element[unique_ct == 1]
    #print("new array",new_array)
    #return the index where the element occurs first e.g. if element 59 occurs at index 9,10 in the array and
    #element 60 occurs at index 11,12, then index returned is =[...8,9,11,13,...]
    unique_values,indices_list = np.unique(temp_new_bin_in_array, return_index=True)
    #unique_values,repeat_indices_list = np.unique(temp_new_bin_in_array, )
    #print("indices list",indices_list)   

    
    #print("repeat_indices",repeat_indices)
    #test_remove_bin=np.delete(test_remove_bin,repeat_indices)
    #print("test remove bin",test_remove_bin)

    
    #Outer bin (i,i+1) here
    temp_dist_outer_bin_in = comparison_distance_inwards[repeat_indices]
    temp_dist_outer_bin_out = comparison_distance_outwards[repeat_indices]
    smallest_sep_outer_bin = np.fmin(temp_dist_outer_bin_in,temp_dist_outer_bin_out)
    #print("outer bin",smallest_sep_outer_bin)
    #Inner bin (i-1,i) here
    temp_dist_in_minus1 = comparison_distance_inwards[repeat_indices-1]
    temp_dist_out_minus1 = comparison_distance_outwards[repeat_indices-1]
    temp_dist_inner_bin = np.fmin(temp_dist_in_minus1,temp_dist_out_minus1)
    #print("inner bin",temp_dist_inner_bin)
    #print("repeat_indices",repeat_indices)
    #These are the indices of binaries to be deleted. 
    #element_to_be_kept = np.where(smallest_sep_outer_bin < temp_dist_inner_bin, repeat_indices, repeat_indices + 1)
    #element_to_be_deleted = np.where(smallest_sep_outer_bin > temp_dist_inner_bin,repeat_indices,repeat_indices + 1)

    new_test_remove_bin = test_remove_bin
    #print("test_remove_bin",new_test_remove_bin)
    #set up array of indices to be removed
    idx_array=np.zeros(len(2*dupl_indices))

    for ind in range(len(dupl_indices)):
        idx = dupl_indices[ind] - (ind)
        value=2*ind
        # If outer bin to be kept    
        if smallest_sep_outer_bin[ind] < temp_dist_inner_bin[ind]:
            #delete inner bin 
            remove_indxs=[idx - 1,idx]
            
            #idx_array[value]=idx - 1
            #idx_array[value+1]=idx
            new_test_remove_bin = np.delete(new_test_remove_bin,remove_indxs)
            #print("test remove",new_test_remove_bin)
        else:
            #delete outer bin
            remove_indxs=[idx+1,idx+2]

            #idx_array[value]=idx+1
            #idx_array[value]=idx+2

            new_test_remove_bin = np.delete(new_test_remove_bin,remove_indxs)
            #print("test remove",new_test_remove_bin)
    
    #print("idx_array",idx_array)

    # Say [10,11,15] Not 9-10;so 10-11 then 13:14.  
    # New=[0:9,11:14,15:end]
    #while ind in range(2*dupl_indices):
    #    new_array=test_remove_bin[1:ind(i),ind(i+2):ind(i+3)-1,ind(i)]
    
    #CHANGE THIS going forward: For removing repeats. Compute mutual Hill sphere (R_H=((M1+M_2)/M_SMBH)^{1/3}) and
    # see how small a fraction of the mutual Hill sphere, the binary separation *would* be.
    #E.g. (r2-R1)=0.5R_H(m1,m2) vs (r3-r2)=0.2R_H(m3,m2) form the tighter binary (r3,r2)

    #Compare nearest neighbour separations for NS in repeat_indices
    # E.g. separations =[r2-r1,r3-r2,r3-4]. If NS at r2 repeats in binary array 
    # e.g. bin_array=[[r1,r2] [r2,r3]..] then compare separations[r2-r1] to separations[r3-r2]
    # If r2-r1 < r3-r2 then make [r1,r2] the binary and remove [r2,r3]
    
    #if separations[repeat_indices] < separations[repeat_indices + 1]:
    #    temp_new_bin_in_array = np.delete(temp_new_bin_in_array,repeat_indices + 1)
    #    print("temp_bin_array",temp_new_bin_in_array)

#    for j in range(2*len(index_in)):    
#        temp_bin_in = new_indx_in_bin_array[temp_index[]]
   

#    for j in range(2*len(index_in)):    
#        temp_bin_in = new_indx_in_bin_array[temp_index[]]
   
    for ind in range(length_index_out):
        temp_index = index_out[ind]
        new_indx_out[2*ind] = temp_index
        new_indx_out[(2*ind)+1] = temp_index+1

    #print("new_indx_out",new_indx_out)    
    temp_new_bin_out_array=np.reshape(new_indx_out,(len(index_out),2))
    #print("ordered as bins",temp_new_bin_out_array)

    new_indxs = new_indx_in+new_indx_out
    #rindx = np.sort(new_indxs)
    result = np.asarray(new_test_remove_bin)
    #result = np.asarray(new_indx_in)
    sorted_in_result = np.sort(result)

    new_result = np.asarray(new_indx_out)
    sorted_out_result = np.sort(new_result)
    #print("sorted in result",sorted_in_result)
    #print("sorted out result",sorted_out_result)
    # Concatenate the two lists, and remove duplicates
    final_bin_indices = np.array(list(set(list(sorted_in_result) + list(sorted_out_result))))
    sorted_final_bin_indices = np.sort(final_bin_indices)
    #print("total final bin indices",sorted_final_bin_indices)
    #print("check if sorted_in & sorted_out arrays are the same")
    check = np.array_equiv(sorted_in_result, sorted_out_result)
    #print(check)
    # Return the indices of those elements in separation array <0
    # (ie NS are closer than 1 R_Hill)
    # In inwards case, r_i+1 -r_i <R_H_i, so relevant NS indices are i,i+1
   
    # In outwards case, r_i - r_i-1 <R_H_i so relevant NS indices are i,i-1
    final_1d_indx_array = sorted_in_result.flatten()
    sorted_final_1d_indx_array = np.sort(final_1d_indx_array)

    #if len(sorted_final_1d_indx_array) > 0:
    #     print("Binary", sorted_final_1d_indx_array)

    return sorted_final_1d_indx_array

def binary_check(prograde_ns_locations, prograde_ns_masses, mass_smbh, prograde_ns_orb_ecc, e_crit):
    """Determines which prograde NS will form binaries in this timestep. Takes as inputs
    the singleton NS locations & masses, and checks if their separations are less than
    the mutual Hill sphere of any 2 adjacent NS. If this is the case, determine the
    smallest separation pairs (in units of their mutual Hill sphere) to form a set of
    actual binaries (this module does handle cases where 3 or more bodies *might* form
    some set of binaries which would be mutually exclusive; however it does not handle
    or even flag the implied triple system dynamics). Returns a 2xN array of the relevant
    binary indices, for further handling to form actual binaries & assign additional
    parameters (e.g. angular momentum of the binary).

    Parameters
    ----------
    prograde_ns_locations : float array
        locations of prograde singleton NS at start of timestep in units of 
        gravitational radii (r_g=GM_SMBH/c^2)
    prograde_ns_masses : float array
        initial masses of ns in prograde orbits around SMBH in units of solar masses
    mass_smbh : float
        mass of supermassive black hole in units of solar masses
    prograde_ns_orb_ecc : float array
        Orbital ecc of singleton NS after damping during timestep
    e_crit : float
        Critical eccentricity allowing bin formation and migration    
    Returns
    -------
    all_binary_indices : [2,N] int array
        array of indices corresponding to locations in prograde_ns_locations, prograde_ns_masses,
        prograde_ns_spins, prograde_ns_spin_angles, and prograde_ns_generations which corresponds
        to binaries that form in this timestep. it has a length of the number of binaries to form (N)
        and a width of 2.
    """

    #print('ns locations',prograde_ns_locations)
    # First sort the prograde ns locations in order from inner disk to outer disk
    sorted_ns_locations = np.sort(prograde_ns_locations)
    #print('sorted ns locations',sorted_ns_locations)
    # Returns the indices of the original array in order, to get the sorted array
    sorted_ns_location_indices = np.argsort(prograde_ns_locations)
    # Returns the indices of the orb ecc array, to get sorted array of orb ecc
    #print('sorted ns location indices',sorted_ns_location_indices)
    #print('ns orb ecc',prograde_ns_orb_ecc)
    sorted_ns_ecc_array = np.empty_like(prograde_ns_orb_ecc)
    sorted_ns_ecc_array = prograde_ns_orb_ecc[np.argsort(prograde_ns_locations)]
    #print('orb ecc of sorted ns', sorted_ns_ecc_array)
    # Find the distances between [r1,r2,r3,r4,..] as [r2-r1,r3-r2,r4-r3,..]=[delta1,delta2,delta3..]
    # Note length of separations is 1 less than prograde_ns_locations
    separations = np.diff(sorted_ns_locations)
    # Now compute mutual hill spheres of all possible binaries
    # same length as separations
    R_Hill_possible_binaries = (sorted_ns_locations[:-1] + separations/2.0) * \
        pow(((prograde_ns_masses[sorted_ns_location_indices[:-1]] + \
              prograde_ns_masses[sorted_ns_location_indices[1:]]) / \
                (mass_smbh * 3.0)), (1.0/3.0))
    # compare separations to mutual Hill spheres - negative values mean possible binary formation
    minimum_formation_criteria = separations - R_Hill_possible_binaries
    # collect indices of possible real binaries (where separation is less than mutual Hill sphere)
    index_formation_criteria = np.where(minimum_formation_criteria < 0)
    # Now deal with sequences: compute separation/R_Hill for all
    sequences_to_test = (separations[index_formation_criteria])/(R_Hill_possible_binaries[index_formation_criteria])
    #print(sequences_to_test)
    # sort sep/R_Hill for all 'binaries' that need checking & store indices
    sorted_sequences = np.sort(sequences_to_test)
    #print(sorted_sequences)
    sorted_sequences_indices = np.argsort(sequences_to_test)
    #print(sorted_sequences_indices)
    # the smallest sep/R_Hill should always form a binary, so
    checked_binary_index = np.array([sorted_sequences_indices[0]])
    #print(checked_binary_index)
    for i in range(len(sorted_sequences)): 
        # if we haven't already counted it
        if (sorted_sequences_indices[i] not in checked_binary_index):
            # and it isn't the implicit partner of something we've already counted
            if (sorted_sequences_indices[i] not in checked_binary_index+1):
                # and the implicit partner of this thing isn't already counted
                if (sorted_sequences_indices[i]+1 not in checked_binary_index):
                    # and the implicit partner of this thing isn't already an implicit partner we've counted
                    if (sorted_sequences_indices[i]+1 not in checked_binary_index+1):
                        # then you can count it as a real binary
                        checked_binary_index = np.append(checked_binary_index, sorted_sequences_indices[i])

    #print("THIS IS SAAVIK'S OUTPUT!!")
    # create array of all real binaries
    # BUT what are we returning? need indices of original arrays
    # that go to singleton vs binary assignments--actual binary formation *should* happen elsewhere!
    # I have the indices of the sorted_ns_locations array that correspond to actual binaries
    # these are checked_binary_index, checked_binary_index+1
    all_binary_indices = np.array([sorted_ns_location_indices[checked_binary_index], sorted_ns_location_indices[checked_binary_index+1]])
    #print(np.shape(all_binary_indices))
    #print(np.shape(all_binary_indices)[1])
    
    #HERE is where we check that both NS have damped orbital eccentricity (e<=0.01). Otherwise do not form binary
    # TO DO: Make new bin formation condition at modest ecc & encounters. 
    # E.g. Remember in a timestep of 10^4yrs, there are 10^4 orbits at R=10^3r_g
    #Singleton NS with orb ecc > e_crit
    prograde_ns_not_form_bins = np.ma.masked_where(prograde_ns_orb_ecc <= e_crit, prograde_ns_orb_ecc)
    #Singleton NS with orb ecc < e_crit
    prograde_ns_can_form_bins = np.ma.masked_where(prograde_ns_orb_ecc >e_crit, prograde_ns_orb_ecc)
    #Indices of singleton NS with orb ecc > e_crit
    indices_ns_not_form_bins = np.ma.nonzero(prograde_ns_not_form_bins) 
    #Indices of singleton NS with orb ecc < e_crit
    indices_ns_can_form_bins = np.ma.nonzero(prograde_ns_can_form_bins)
    

    #print('ORIGINAL INDICES NS ALLOWED FORM BINS',np.array(indices_ns_can_form_bins[0]))
    allowed_to_form_bins = np.array(indices_ns_can_form_bins[0])
    #print('allowed ns locs',prograde_ns_locations[allowed_to_form_bins])
    sorted_allowed_ns_loc = np.sort(prograde_ns_locations[allowed_to_form_bins])
    #print('sorted allowed ns locs',sorted_allowed_ns_loc)
    #print('allowed ns eccs',prograde_ns_orb_ecc[allowed_to_form_bins])
    #print('allowed ns masses',prograde_ns_masses[allowed_to_form_bins])
    allowed_separations = np.diff(sorted_allowed_ns_loc)
    # Now compute mutual hill spheres of all possible binaries
    # same length as separations
    R_Hill_allowed_bin_test = (sorted_allowed_ns_loc[:-1] + allowed_separations/2.0) * \
        pow(((prograde_ns_masses[allowed_to_form_bins[:-1]] + \
              prograde_ns_masses[allowed_to_form_bins[1:]]) / \
                (mass_smbh * 3.0)), (1.0/3.0))
    # compare separations to mutual Hill spheres - negative values mean possible binary formation
    allowed_min_form_criteria = allowed_separations - R_Hill_allowed_bin_test
    print('criteria',allowed_min_form_criteria)
    # collect indices of possible real binaries (where separation is less than mutual Hill sphere)
    allowed_indx_form_criteria = np.where(allowed_min_form_criteria < 0)
    print(allowed_indx_form_criteria)
    allowed_idx_crit = allowed_indx_form_criteria[0]
    print(allowed_idx_crit)

    if np.size(allowed_indx_form_criteria) >0: 
        #If multiple negative results in criteria
        item1 = np.empty(len(allowed_idx_crit))
        item2 = np.empty(len(allowed_idx_crit))
        item1_idx = np.empty(len(allowed_idx_crit))
        item2_idx = np.empty(len(allowed_idx_crit))
        idx1 = np.empty(len(allowed_idx_crit))
        idx2 = np.empty(len(allowed_idx_crit))
        for i in range(len(allowed_idx_crit)):
            item1[i] = sorted_allowed_ns_loc[allowed_idx_crit[i]]
            item2[i] = sorted_allowed_ns_loc[allowed_idx_crit[i]+1]         
            print(item1[0],item2[0])
            item1_idx = np.where(prograde_ns_locations == item1[i])
            item2_idx = np.where(prograde_ns_locations == item2[i])
            print(item1_idx[0],item2_idx[0])
            idx1 = item1_idx[0] 
            idx2 = item2_idx[0]
            print(idx1,idx2)

        #for j in range(len(idx1)):
            #print(idx1[j],idx2[j]) 
            #print(prograde_ns_locations[idx1[j]],prograde_ns_locations[idx2[j]],prograde_ns_orb_ecc[idx1[j]],prograde_ns_orb_ecc[idx2[j]])
    
    # Now deal with sequences: compute separation/R_Hill for all
    allowed_sequences_to_test = (allowed_separations[allowed_indx_form_criteria])/(R_Hill_allowed_bin_test[allowed_indx_form_criteria])
    #print('allowed seqs to test', allowed_sequences_to_test)
    # sort sep/R_Hill for all 'binaries' that need checking & store indices
    sorted_allowed_sequences = np.sort(allowed_sequences_to_test)
    #print(sorted_allowed_sequences)
    sorted_allowed_sequences_indices = np.argsort(allowed_sequences_to_test)
    #print(sorted_allowed_sequences_indices)
    # the smallest sep/R_Hill should always form a binary, so
    if np.count_nonzero(sorted_allowed_sequences_indices) > 0:
        allowed_checked_binary_index = np.array([sorted_allowed_sequences_indices[0]],dtype = int)
        #print(allowed_checked_binary_index)
    
    for i in range(len(sorted_allowed_sequences)): 
        allowed_checked_binary_index = np.array([sorted_allowed_sequences_indices[0]],dtype = int)
        # if we haven't already counted it
        if (sorted_allowed_sequences_indices[i] not in allowed_checked_binary_index):
            # and it isn't the implicit partner of something we've already counted
            if (sorted_allowed_sequences_indices[i] not in allowed_checked_binary_index+1):
                # and the implicit partner of this thing isn't already counted
                if (sorted_allowed_sequences_indices[i]+1 not in allowed_checked_binary_index):
                    # and the implicit partner of this thing isn't already an implicit partner we've counted
                    if (sorted_allowed_sequences_indices[i]+1 not in allowed_checked_binary_index+1):
                        # then you can count it as a real binary
                        allowed_checked_binary_index = np.append(allowed_checked_binary_index, sorted_allowed_sequences_indices[i])

    
    #Check if any of the checked binary index array are in the array of indices that can form bins
    #check_overlap = np.isin(checked_binary_index,allowed_to_form_bins)
    #final_overlap = np.array(check_overlap.nonzero()[0])
    #print('final overlap',final_overlap)
    #print('sorted ns locs',sorted_ns_locations[final_overlap])
    #print('orb eccs',sorted_ns_ecc_array[final_overlap])
    #Check for adjacent integers in final overlap
    #Take the difference of integers in final overlap
    #diffs_final_overlap = np.diff(final_overlap)
    #print(diffs_final_overlap)
    #possible_bins = np.where(diffs_final_overlap == 1)
    #print(possible_bins)
    #print('final checked binary index ',checked_binary_index)
    
    # create array of all real binaries
    # BUT what are we returning? need indices of original arrays
    # that go to singleton vs binary assignments--actual binary formation *should* happen elsewhere!
    # I have the indices of the sorted_ns_locations array that correspond to actual binaries
    # these are checked_binary_index, checked_binary_index+1
    if np.count_nonzero(sorted_allowed_sequences) > 0:
        final_binary_indices = np.array([idx1,idx2])
#        final_binary_indices = np.array([sorted_ns_location_indices[allowed_checked_binary_index], sorted_ns_location_indices[allowed_checked_binary_index+1]])
        print(np.shape(final_binary_indices))
        print(np.shape(final_binary_indices)[1])
        #print(prograde_ns_locations[allowed_checked_binary_index],prograde_ns_locations[allowed_checked_binary_index+1])
    else: 
        final_binary_indices=np.empty_like(allowed_sequences_to_test)
    #return all_binary_indices
    return final_binary_indices

def binary_check2(prograde_ns_locations, prograde_ns_masses, mass_smbh, prograde_ns_orb_ecc, e_crit):
    """Determines which prograde NS will form binaries in this timestep. Takes as inputs
    the singleton NS locations,masses & orbital eccentricities, and takes the candidate binary population from 
    NS with orbital eccentricities damped to <e_crit.
    Among this damped population, checks if their separations are less than
    the mutual Hill sphere of any 2 adjacent NS. If this is the case, determine the
    smallest separation pairs (in units of their mutual Hill sphere) to form a set of
    actual binaries (this module does handle cases where 3 or more bodies *might* form
    some set of binaries which would be mutually exclusive; however it does not handle
    or even flag the implied triple system dynamics). Returns a 2xN array of the relevant
    binary indices, for further handling to form actual binaries & assign additional
    parameters (e.g. angular momentum of the binary).

    Parameters
    ----------
    prograde_ns_locations : float array
        locations of prograde singleton NS at start of timestep in units of 
        gravitational radii (r_g=GM_SMBH/c^2)
    prograde_ns_masses : float array
        initial masses of ns in prograde orbits around SMBH in units of solar masses
    mass_smbh : float
        mass of supermassive black hole in units of solar masses
    prograde_ns_orb_ecc : float array
        Orbital ecc of singleton NS after damping during timestep
    e_crit : float
        Critical eccentricity allowing bin formation and migration    
    Returns
    -------
    all_binary_indices : [2,N] int array
        array of indices corresponding to locations in prograde_ns_locations, prograde_ns_masses,
        prograde_ns_spins, prograde_ns_spin_angles, prograde_ns_orb_ecc and prograde_ns_generations which corresponds
        to binaries that form in this timestep. it has a length of the number of binaries to form (N)
        and a width of 2.
    """

    #First check for NS with sufficiently damped orbital eccentricity (e<=e_crit (usually 0.01)). 
    #This population is the sub-set of prograde NS from which we CAN form binaries.
    
    #Singleton NS with orb ecc < e_crit (candidates for binary formation)
    prograde_ns_can_form_bins = np.ma.masked_where(prograde_ns_orb_ecc >e_crit, prograde_ns_orb_ecc)    
    indices_ns_can_form_bins = np.ma.nonzero(prograde_ns_can_form_bins)
    # Indices of those candidates for binary formation
    allowed_to_form_bins = np.array(indices_ns_can_form_bins[0])
    #Sort the location of the candidates
    sorted_ns_locations = np.sort(prograde_ns_locations[allowed_to_form_bins])
    #Sort the indices of all singleton NS (the superset)
    sorted_ns_location_indices_superset = np.argsort(prograde_ns_locations)
    #Set the condition for membership in candidate array to be searched/tested
    condition = np.isin(sorted_ns_location_indices_superset, allowed_to_form_bins)
    #Here is the subset of indices that can be tested for binarity
    subset = np.extract(condition,sorted_ns_location_indices_superset)
    
    # Find the distances between [r1,r2,r3,r4,..] as [r2-r1,r3-r2,r4-r3,..]=[delta1,delta2,delta3..]
    # Note length of separations is 1 less than prograde_ns_locations
    
    #This is the set of separations between the sorted candidate NS
    separations = np.diff(sorted_ns_locations)
 
    # Now compute mutual hill spheres of all possible candidate binaries if can test    
    if len(separations) > 0:
        R_Hill_possible_binaries = (sorted_ns_locations[:-1] + separations/2.0) * \
            pow(((prograde_ns_masses[subset[:-1]] + \
                  prograde_ns_masses[subset[1:]]) / \
                    (mass_smbh * 3.0)), (1.0/3.0))
        # compare separations to mutual Hill spheres - negative values mean possible binary formation
        minimum_formation_criteria = separations - R_Hill_possible_binaries
        #print("sep -R_hill", minimum_formation_criteria)
        # collect indices of possible real binaries (where separation is less than mutual Hill sphere)
        index_formation_criteria = np.where(minimum_formation_criteria < 0)
        
        #Here's the index of the array of candidates
        test_idx = index_formation_criteria[0]
        #print("test_idx", test_idx)
        
        #If we actually have any candidates this time step
        if np.size(test_idx) >0:
            #print("subset(test_idx)",subset[test_idx][0],subset[test_idx+1][0])
            #Start with real index (from full singleton array) of 1st candidate binary component (implicit + 1 partner since separations are ordered )
            bin_indices = np.array([subset[test_idx[0]],subset[test_idx[0]+1]])
            #If only 1 binary this timestep, return this binary!
            all_binary_indices = np.array([subset[test_idx],subset[test_idx+1]])
            
            for i in range(len(test_idx)):
                #If more than 1 binary
                if i >0:
                    # append nth binary indices formed this timestep
                    bin_indices = np.append(bin_indices,[subset[test_idx[i]],subset[test_idx[i]+1]])
                    #print("2XBin", bin_indices)
                    
                    #Check to see if repeat binaries among the set of binaries formed (e.g. (1,2)(2,3) )
                    #If repeats, only form a binary from the pair with smallest fractional Hill sphere separation

                    # Compute separation/R_Hill for all
                    sequences_to_test = (separations[test_idx])/(R_Hill_possible_binaries[test_idx])
                    #print(sequences_to_test)
                    # sort sep/R_Hill for all 'binaries' that need checking & store indices
                    sorted_sequences = np.sort(sequences_to_test)
                    #print(sorted_sequences)
                    #Sort the indices for the test
                    sorted_sequences_indices = np.argsort(sequences_to_test)
                    #print(sorted_sequences_indices)

                    # Assume the smallest sep/R_Hill should form a binary, so
                    if len(sorted_sequences) > 0:
                        #Index of smallest sorted fractional Hill radius binary so far
                        checked_binary_index = np.array([test_idx[sorted_sequences_indices[0]]])
                    else:
                        checked_binary_index = []    
                    #print(checked_binary_index)
                    for j in range(len(sorted_sequences)): 
                        # if we haven't already counted it
                        if (test_idx[sorted_sequences_indices[j]] not in checked_binary_index):
                            # and it isn't the implicit partner of something we've already counted
                            if (test_idx[sorted_sequences_indices[j]] not in checked_binary_index+1):
                                # and the implicit partner of this thing isn't already counted
                                if (test_idx[sorted_sequences_indices[j]]+1 not in checked_binary_index):
                                    # and the implicit partner of this thing isn't already an implicit partner we've counted
                                    if (test_idx[sorted_sequences_indices[j]]+1 not in checked_binary_index+1):
                                        # then you can count it as a real binary
                                        checked_binary_index = np.append(checked_binary_index, test_idx[sorted_sequences_indices[j]])
                                        #print("checked_binary_index",checked_binary_index)
                    all_binary_indices = np.array([subset[checked_binary_index],subset[checked_binary_index+1]])

            
            

        else:
            #No binaries from candidates this time step
            all_binary_indices = []
        
    else:
        #No candidate for binarity testing yet
        all_binary_indices = []
    
    return all_binary_indices
        
    