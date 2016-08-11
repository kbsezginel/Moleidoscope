# Linker class for molecular kaleidoscope
# Date: August 2016
# Authors: Kutay B. Sezginel and Yan Gui
import os

main_dir = os.getcwd()
library_path = os.path.join(main_dir, 'LIBRARY')
library_file = open(library_path, 'r')
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
atom_coors_symmetry = []
for num_index, num in enumerate(number_of_atoms):
    coor_index = coordinate_index[num_index]
    connect_index = connectivity_index[num_index]
    atom_names.append([])
    atom_coors.append([])
    atom_coors_symmetry.append([])
    for i in range(int(num)):
        atom_name = lib_lines[coor_index + i].split()[1]
        if atom_name == 'X':
            atom_name = 'O'
        x = float(lib_lines[coor_index + i].split()[2])
        y = float(lib_lines[coor_index + i].split()[3])
        z = float(lib_lines[coor_index + i].split()[4])
        x1 = float(x * (-1))
        atom_names[num_index].append(atom_name)
        atom_coors[num_index].append([x, y, z])
        atom_coors_symmetry[num_index].append([x1, y, z])


class Linker:
    """
    Linker class.
    """
    def __init__(self, linker_index):
        linker_info = self.read_linker(linker_index)
        self.index = linker_index
        self.name = linker_info['name']
        self.atom_names = linker_info['atom_names']
        self.atom_coors = linker_info['atom_coors']
        self.num_of_atoms = linker_info['num_of_atoms']
        self.connectivity = linker_info['connectivity']

    def read_linker(self, linker_index):
        linker_info = {'atom_names': [], 'atom_coors': [], 'num_of_atoms': 0, 'connectivity': []}
        linker_info['num_of_atoms'] = number_of_atoms[linker_index]
        linker_info['atom_names'] = atom_names[linker_index]
        linker_info['atom_coors'] = atom_coors[linker_index]
        linker_info['name'] = linker_names[linker_index]
        return linker_info

    def mirror(self, mirror_info):
        if mirror_info == [1, 1, 1]:
            mirror_operator = [-1, -1, -1]
        else:
            mirror_operator = [1, 1, 1]
            for axis_index, axis in enumerate(mirror_info):
                if axis == 0:
                    mirror_operator[axis_index] = -1
        mirror_coordinates = []
        mirror_names = []
        for coor, name in zip(self.atom_coors, self.atom_names):
            mirror_coor = []
            for axis_index, mirror_axis in enumerate(mirror_operator):
                mirror_coor.append(coor[axis_index] * mirror_axis)
            mirror_names.append(name)
            mirror_coordinates.append(mirror_coor)

        mirror_linker = Linker(self.index)
        mirror_linker.atom_coors = mirror_coordinates
        mirror_linker.atom_names = mirror_names
        return mirror_linker

    def translate(self, translation_info):
        translation_coordinates = []
        for coor in self.atom_coors:
            translation_coor = []
            for axis_index, translation_amount in enumerate(translation_info):
                translation_coor.append(coor[axis_index] + translation_amount)
            translation_coordinates.append(translation_coor)
        self.atom_coors = translation_coordinates

    def join(self, other_linker):
        joined_linker = Linker(self.index)
        joined_linker.atom_coors = self.atom_coors + other_linker.atom_coors
        joined_linker.atom_names = self.atom_names + other_linker.atom_names
        joined_linker.num_of_atoms = len(self.atom_names) + len(other_linker.atom_names)
        if self.name == other_linker.name:
            joined_linker.name = self.name + 'JOINED'
        else:
            joined_linker.name = self.name + '_' + other_linker.name
        return joined_linker

    def export(self):
        export_index = os.getcwd()
        linker_path = os.path.join(export_index, self.name + '.xyz')
        linker_file = open(linker_path, 'w')
        linker_file.write(str(len(self.atom_coors)) + '\n')
        linker_file.write(self.name + '\n')
        for atom, coor in zip(self.atom_names, self.atom_coors):
            coor_line = atom + ' ' + str(coor[0]) + ' '
            coor_line += str(coor[1]) + ' ' + str(coor[2]) + '\n'
            linker_file.write(coor_line)
        linker_file.close()
