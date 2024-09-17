"""Unit test for test_tgw.py"""
import numpy as np
import pytest
import scipy

from mcfacts.physics.binary.formation.add_new_binary import add_to_binary_array2
from mcfacts.utility import unit_test_variable_provider as provider
from mcfacts.utility.unit_test_variable_provider import InputParameterSet

INPUT_PARAMETERS = "disk_bins_bhbh, disk_bh_pro_orbs_a, disk_bh_pro_masses, disk_bh_pro_spins, \
                       disk_bh_pro_spin_angles, disk_bh_pro_gens, disk_bin_bhbh_pro_indices, bindex,\
                       fraction_bin_retro, smbh_mass"


# @pytest.mark.parametrize("disk_bh_pro_orbs_a", provider.INPUT_PARAMETERS["bh_orbital_semi_major_axis_inner"][InputParameterSet.SINGLETON])
# @pytest.mark.parametrize("disk_bh_pro_masses", provider.INPUT_PARAMETERS["bh_masses"][InputParameterSet.SINGLETON])
# @pytest.mark.parametrize("disk_bh_pro_spins", provider.INPUT_PARAMETERS["bh_spins"][InputParameterSet.SINGLETON])
# @pytest.mark.parametrize("disk_bh_pro_spin_angles", provider.INPUT_PARAMETERS["bh_spin_angles"][InputParameterSet.SINGLETON])
# @pytest.mark.parametrize("disk_bh_pro_gens", [0, 1, 2, 3])
# @pytest.mark.parametrize("disk_bin_bhbh_pro_indices", [0])
# @pytest.mark.parametrize("bindex", [0])
# @pytest.mark.parametrize("fraction_bin_retro", [0, .1, .5, .6, 1])
# @pytest.mark.parametrize("smbh_mass", provider.INPUT_PARAMETERS["smbh_mass"][InputParameterSet.BASE])
def test_add_to_binary_array2(disk_bh_pro_orbs_a, disk_bh_pro_masses, disk_bh_pro_spins,
        disk_bh_pro_spin_angles, disk_bh_pro_gens, disk_bin_bhbh_pro_indices, bindex, fraction_bin_retro, smbh_mass):
    """test function"""

    return_value = add_to_binary_array2([], disk_bh_pro_orbs_a, disk_bh_pro_masses, disk_bh_pro_spins,
                                        disk_bh_pro_spin_angles, disk_bh_pro_gens, disk_bin_bhbh_pro_indices, bindex,
                                        fraction_bin_retro, smbh_mass)

    # base_value = test([], disk_bh_pro_orbs_a, disk_bh_pro_masses, disk_bh_pro_spins,
    #                                     disk_bh_pro_spin_angles, disk_bh_pro_gens, disk_bin_bhbh_pro_indices, bindex,
    #                                     fraction_bin_retro, smbh_mass)

    assert return_value == base_value