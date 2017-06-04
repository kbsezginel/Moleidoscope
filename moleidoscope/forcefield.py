# Moleidoscope Force Field Functions
# Date: March 2017
# Author: Kutay B. Sezginel
import os
import math
import xlrd
import numpy as np


def read_ff_parameters(excel_file_path, ff_selection='uff'):
    """
    Read force field parameters from an excel file according to force field selection
    """
    # Read Excel File
    force_field_data = xlrd.open_workbook(excel_file_path)
    # Read columns to acquire force field parameters
    atom_names = force_field_data.sheets()[0].col_values(0)[2:]
    uff_sigma = force_field_data.sheets()[0].col_values(1)[2:]
    uff_epsilon = force_field_data.sheets()[0].col_values(2)[2:]
    dre_sigma = force_field_data.sheets()[0].col_values(3)[2:]
    dre_epsilon = force_field_data.sheets()[0].col_values(4)[2:]

    uff = {'atom': atom_names, 'sigma': uff_sigma, 'epsilon': uff_epsilon}
    dre = {'atom': atom_names, 'sigma': dre_sigma, 'epsilon': dre_epsilon}

    if ff_selection == 'uff':
        return uff
    if ff_selection == 'dre':
        return dre
    else:
        print('No such force field')


def get_ff_par(atom_name, ff_parameters):
    """ Get sigma and epsilon values for given atom name and force field parameters dictionary. """
    atom_index = ff_parameters['atom'].index(atom_name)
    sigma = ff_parameters['sigma'][atom_index]
    epsilon = ff_parameters['epsilon'][atom_index]
    return sigma, epsilon


def lb_mix(sigma_1, sigma_2, epsilon_1, epsilon_2):
    """ Calculate Lorentz-Barthelot mixing rules """
    sigma_mix = (sigma_1 + sigma_2) / 2
    epsilon_mix = math.sqrt(epsilon_1 * epsilon_2)
    return sigma_mix, epsilon_mix


def lennard_jones(r, sig, eps):
    """
    Calculate Lennard Jones potential for given distance, sigma, and epsilon values.
    Energy unit: (kB)
    """
    return 4 * eps * ((sig / r)**12 - (sig / r)**6)
