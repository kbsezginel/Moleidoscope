# Moleidoscope hostdesigner methods
# Author: Kutay B Sezginel
# Date: February 2017
import os
import numpy as np


def read_library(library_path, rename='N'):
    """
    Read HostDesigner linker library.
    """
    with open(library_path, 'r') as library_file:
        lib_lines = library_file.readlines()

    connectivity_index = []
    connectivity = []
    coordinate_index = []
    number_of_atoms = []
    linker_index = []
    linker_names = []
    num_of_linkers = 1
    for line_index, line in enumerate(lib_lines):
        if 'LINK' in line:
            linker_names.append(lib_lines[line_index].split()[-1])
            linker_index.append(num_of_linkers)
            connectivity_index.append(line_index + 1)
            coordinate_index.append(line_index + 6)
            num_atom = float(lib_lines[line_index + 5].split()[0])
            number_of_atoms.append(num_atom)
            num_of_linkers += 1

    atom_names = []
    atom_coors = []
    for num_index, num in enumerate(number_of_atoms):
        coor_index = coordinate_index[num_index]
        connect_index = connectivity_index[num_index]
        atom_names.append([])
        atom_coors.append([])
        for i in range(int(num)):
            atom_name = lib_lines[coor_index + i].split()[1]
            if atom_name == 'X':
                atom_name = rename
            x = float(lib_lines[coor_index + i].split()[2])
            y = float(lib_lines[coor_index + i].split()[3])
            z = float(lib_lines[coor_index + i].split()[4])
            atom_names[num_index].append(atom_name)
            atom_coors[num_index].append([x, y, z])
        # Add connectivity information
        connectivity.append(read_connnection(lib_lines, connect_index, atom_coors[num_index]))

    library = {'atom_names': atom_names,
               'atom_coors': atom_coors,
               'number_of_atoms': number_of_atoms,
               'linker_index': linker_index,
               'coordinate_index': coordinate_index,
               'connectivity_index': connectivity_index,
               'connectivity': connectivity,
               'linker_names': linker_names}

    return library


def read_connnection(library_lines, connectivity_index, atom_coors):
    """ Read connectivity information for given linker from library """
    connectivity = dict(dummy_dist=0, carbon_dist=0, angle1=0, angle2=0, dihedral=0, atoms=[], coors=[])
    atom1_index = int(library_lines[connectivity_index].split()[0])
    atom2_index = int(library_lines[connectivity_index + 1].split()[0])
    connectivity['atoms'] = [atom1_index, atom2_index]
    connectivity['carbon_dist'] = float(library_lines[connectivity_index + 3].split()[0])
    connectivity['angle1'] = float(library_lines[connectivity_index + 3].split()[1])
    connectivity['angle2'] = float(library_lines[connectivity_index + 3].split()[2])
    connectivity['dihedral'] = float(library_lines[connectivity_index + 3].split()[3])
    try:
        connect1_coor = np.array(atom_coors[atom1_index - 1])
        connect2_coor = np.array(atom_coors[atom2_index - 1])
        connectivity['dummy_dist'] = np.linalg.norm(connect1_coor - connect2_coor)
        connectivity['coors'] = [connect1_coor, connect2_coor]
    except:
        connectivity['dummy_dist'] = None
        connectivity['coors'] = None
    return connectivity
