# Moleidoscope hostdesigner methods
# Author: Kutay B Sezginel
# Date: February 2017
import os


def read_library(library_path):
    """
    Read HostDesigner linker library.
    """
    with open(library_path, 'r') as library_file:
        lib_lines = library_file.readlines()

    connectivity_index = []
    coordinate_index = []
    number_of_atoms = []
    linker_index = []
    linker_names = []
    num_of_linkers = 0
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
                atom_name = 'O'
            x = float(lib_lines[coor_index + i].split()[2])
            y = float(lib_lines[coor_index + i].split()[3])
            z = float(lib_lines[coor_index + i].split()[4])
            atom_names[num_index].append(atom_name)
            atom_coors[num_index].append([x, y, z])

    library = {'atom_names': atom_names,
               'atom_coors': atom_coors,
               'number_of_atoms': number_of_atoms,
               'linker_index': linker_index,
               'coordinate_index': coordinate_index,
               'connectivity_index': connectivity_index,
               'linker_names': linker_names}

    return library
