# Linker class for molecular kaleidoscope
# Date: August 2016
# Authors: Kutay B. Sezginel and Yan Gui
import os
import math
import nglview
import numpy as np
from moleidoscope.geo.quaternion import Quaternion
from moleidoscope.mirror import Mirror
from moleidoscope.hd import read_library
from moleidoscope.output import write_pdb


hd_dir = os.environ['HD_DIR']
library_path = os.path.join(hd_dir, 'LIBRARY')
if os.path.exists(library_path):
    hd_lib = read_library(library_path)
else:
    print('HostDesigner directory not found in environment variables!')


class Linker:
    """
    Linker class.
    """
    def __init__(self, linker_index=None, host=None):
        if linker_index is not None:
            self.read_linker(linker_index)
        elif host is not None:
            self.read_host(host)
        else:
            self.name = ''
            self.atom_names = []
            self.atom_coors = []

    def __repr__(self):
        return "<Linker object %s with:%s atoms>" % (self.name, len(self.atom_coors))

    def copy(self):
        """ Copy linker object """
        linker_copy = Linker()
        variables = [i for i in vars(self)]
        for var in variables:
            setattr(linker_copy, var, getattr(self, var))
        return linker_copy

    def read_linker(self, linker_index, library=hd_lib):
        """ Read linker information into object from the library """
        self.index = linker_index
        self.name = library['linker_names'][linker_index - 1]
        self.atom_names = library['atom_names'][linker_index - 1]
        self.atom_coors = library['atom_coors'][linker_index - 1]
        self.num_of_atoms = library['number_of_atoms'][linker_index - 1]
        self.connectivity = library['connectivity'][linker_index - 1]

    def read_host(self, host):
        """ Read host object into linker object """
        self.name = host.name
        self.atom_coors = host.atom_coors
        self.atom_names = host.atom_names
        self.host = host

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

    def translate(self, translation_vector):
        """ Translate linker using given vector """
        translation = np.array(translation_vector)
        translation_coordinates = []
        for coor in self.atom_coors:
            translation_coordinates.append(np.array(coor) + translation)
        self.atom_coors = translation_coordinates

    def rotate(self, angle, axis):
        """ Rotate linker with given angle and axis """
        Q = Quaternion([0, 1, 1, 1])
        rotated_coors = []
        for coor in self.atom_coors:
            x, y, z = Q.rotation(coor, [0, 0, 0], axis, angle).xyz()
            rotated_coors.append([x, y, z])

        rotated_linker = self.copy()
        rotated_linker.name = '%s_R' % self.name
        rotated_linker.atom_coors = rotated_coors
        return rotated_linker

    def rotoreflect(self, angle, axis, mirror_plane, translate=None):
        """ Improper rotation with given angle, axis and reflection plane """
        rotated_linker = self.rotate(angle, axis)
        reflected_linker = rotated_linker.reflect(mirror_plane, translate=translate)
        return reflected_linker

    def get_center(self):
        """ Get center coordinates of the linker """
        xsum, ysum, zsum = 0, 0, 0
        for coor in self.atom_coors:
            xsum += coor[0]
            ysum += coor[1]
            zsum += coor[2]
        num_of_atoms = len(self.atom_coors)
        return [xsum / num_of_atoms, ysum / num_of_atoms, zsum / num_of_atoms]

    def center(self, coor=[0, 0, 0], mirror=None):
        """ Move linker to given coordinates using it's center """
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
        """ Join multiple linker objects into single linker object """
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

    def save(self, file_name=None, save_dir=None):
        """ Save linker object """
        if file_name is None:
            file_name = self.name
        if save_dir is None:
            save_dir = os.getcwd()
        linker_path = os.path.join(save_dir, file_name + '.pdb')
        with open(linker_path, 'w') as pdb:
            write_pdb(pdb, self.atom_names, self.atom_coors, header=self.name)
        return linker_path
