# Linker class for molecular kaleidoscope
# Date: August 2016
# Authors: Kutay B. Sezginel and Yan Gui
import os
import math
import nglview
import numpy as np
from moleidoscope.geo.coor import Coor
from moleidoscope.geo.quaternion import Quaternion
from moleidoscope.mirror import Mirror
from moleidoscope.hd import read_library
from moleidoscope.visualize import write_pdb


main_dir = os.getcwd()
library_path = os.path.join(main_dir, 'LIBRARY')
linker_export_dir = os.path.join(main_dir, 'doc', 'tmp')
library = read_library(library_path)


class Linker:
    """
    Linker class.
    """
    def __init__(self, linker_index=None):
        if linker_index is not None:
            linker_info = self.read_linker(linker_index)
            self.index = linker_index
            self.name = linker_info['name']
            self.atom_names = linker_info['atom_names']
            self.atom_coors = linker_info['atom_coors']
            self.num_of_atoms = linker_info['num_of_atoms']
            self.connectivity = linker_info['connectivity']
        else:
            self.name = ''
            self.atom_names = []
            self.atom_coors = []

    def __repr__(self):
        return "<Linker object %s with:%s atoms>" % (self.name, len(self.atom_coors))

    def copy(self):
        copy_linker = Linker()
        copy_linker.atom_coors = []
        copy_linker.atom_names = []
        for name, coor in zip(self.atom_names, self.atom_coors):
            copy_linker.atom_names.append(name)
            copy_linker.atom_coors.append(coor)
        copy_linker.num_of_atoms = len(copy_linker.atom_names)
        copy_linker.name = self.name
        return copy_linker

    def read_linker(self, linker_index):
        linker_info = {'atom_names': [], 'atom_coors': [], 'num_of_atoms': 0, 'connectivity': []}
        linker_info['num_of_atoms'] = library['number_of_atoms'][linker_index + 1]
        linker_info['atom_names'] = library['atom_names'][linker_index + 1]
        linker_info['atom_coors'] = library['atom_coors'][linker_index + 1]
        linker_info['name'] = library['linker_names'][linker_index + 1]
        linker_info['connectivity'] = library['connectivity_index'][linker_index + 1]
        return linker_info

    def reflect(self, mirror_plane, translate=None):
        """ Reflect the linker off a plane.
            - translate: Translate the linker after reflection
            - amount: Translation amount
              (default is 0 which means no additional translation after reflection)
        """
        if isinstance(mirror_plane, list) and len(mirror_plane) == 3:
            p1, p2, p3 = mirror_plane
            m = Mirror(p1, p2, p3)
        elif isinstance(mirror_plane, str):
            m = Mirror(mirror_plane)
        else:
            m = mirror_plane

        mirror_coordinates = []
        mirror_names = []
        for coor, name in zip(self.atom_coors, self.atom_names):
            mirror_coor = m.symmetry(coor)
            mirror_coordinates.append(mirror_coor)
            mirror_names.append(name)

        mirror_linker = Linker()
        mirror_linker.name = self.name + '_M'
        mirror_linker.atom_coors = mirror_coordinates
        mirror_linker.atom_names = mirror_names

        if translate is not None:
            normal_vector = [m.a, m.b, m.c]
            translation_vector = [i * translate for i in normal_vector]
            # translation_vector = [(i * amount for i, j in zip(self.atom_coors[0], mirror_coordinates[0])]
            mirror_linker.translate(translation_vector)

        return mirror_linker

    def translate(self, translation_info):
        translation_coordinates = []
        for coor in self.atom_coors:
            translation_coor = []
            for axis_index, translation_amount in enumerate(translation_info):
                translation_coor.append(coor[axis_index] + translation_amount)
            translation_coordinates.append(translation_coor)
        self.atom_coors = translation_coordinates

    def rotate(self, rotation_info):
        Q = Quaternion([0, 1, 1, 1])
        rotated_coors = []
        for coor in self.atom_coors:
            x, y, z = Q.rotation(coor, rotation_info[2], rotation_info[1], rotation_info[0]).xyz()
            rotated_coors.append([x, y, z])

        rotated_linker = self.copy()
        rotated_linker.name = self.name + '_R' + str(round(math.degrees(rotation_info[0])))
        rotated_linker.atom_coors = rotated_coors
        return rotated_linker

    def rotoreflect(self, rotation_info, mirror_plane, translate=None):
        rotated_linker = self.rotate(rotation_info)
        reflected_linker = rotated_linker.reflect(mirror_plane, translate=translate)
        return reflected_linker

    def get_center(self):
        xsum = 0
        ysum = 0
        zsum = 0
        for coor in self.atom_coors:
            xsum += coor[0]
            ysum += coor[1]
            zsum += coor[2]
        num_of_atoms = len(self.atom_coors)
        return [xsum / num_of_atoms, ysum / num_of_atoms, zsum / num_of_atoms]

    def center(self, coor=[0, 0, 0], mirror=None):
        if mirror is None:
            center_coor = self.get_center()
            center_vector = [i - j for i, j in zip(coor, center_coor)]
            new_coors = []
            for atom in self.atom_coors:
                new_coors.append([i + j for i, j in zip(atom, center_vector)])
                self.atom_coors = new_coors
        else:
            if type(mirror) is list:
                mir, scale = mirror
            else:
                mir = mirror
                scale = 0
            linker_center = np.array(self.get_center())
            mirror_center = np.array(mir.get_center())
            mirror_vector = np.array([mir.a, mir.b, mir.c]) * scale
            translation = mirror_center - linker_center + mirror_vector
            self.translate(translation)

    def join(self, *args):
        joined_linker = self.copy()
        for other_linker in args:
            joined_linker.atom_coors += other_linker.atom_coors
            joined_linker.atom_names += other_linker.atom_names
            joined_linker.num_of_atoms += len(other_linker.atom_names)
            if joined_linker.name == other_linker.name:
                joined_linker.name += 'JOINED'
            else:
                joined_linker.name += '_' + other_linker.name
        return joined_linker

    def view(self):
        linker_path = self.export()
        return nglview.show_structure_file(linker_path)

    def export(self, file_name=None, export_dir=None):
        if file_name is None:
            file_name = self.name
        if export_dir is None:
            export_dir = linker_export_dir
        linker_path = os.path.join(export_dir, file_name + '.pdb')
        write_pdb(linker_path, self.atom_names, self.atom_coors)
        return linker_path
