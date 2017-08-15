# Date: August 2016
# Author: Kutay B. Sezginel
"""
Polyhedra object to create polyhedral molecules
"""
import os
import copy
import math
import yaml
import numpy as np
from random import randint
from moleidoscope.geo.quaternion import Quaternion
from moleidoscope.geo.vector import align
from moleidoscope.forcefield import get_ff_par, lennard_jones, lb_mix, read_ff_parameters
from moleidoscope.linker import Linker
from moleidoscope.output import save


class Polyhedra:
    """ Polyhedra object."""
    def __init__(self, lib=None, name=None, atom='C'):
        if name is not None:
            self.path = os.path.join(lib, '%s.yaml' % name)
            self.load(self.path, atom='C')
        self.name = name
        self.lib = lib

    def load(self, polyhedra_path, atom='C'):
        """ Load polyhedra yaml file """
        ph = yaml.load(open(polyhedra_path, 'r'))
        self.vertices = ph['vertices']
        self.edges = ph['edges']
        self.faces = ph['faces']
        self.size = ph['size']
        self.symbol = ph['symbol']
        self.atom_coors = self.vertices
        self.atom_names = [atom] * len(self.atom_coors)

    def skeleton(self, atom='C'):
        """ Return new polyhedra only with vertices """
        skeleton = Polyhedra()
        skeleton.atom_coors = self.vertices
        skeleton.atom_names = [atom] * len(skeleton.atom_coors)
        return skeleton

    def resize(self, size):
        """ Resize the polyhedra to given size """
        scaling_coefficient = size / self.size
        scaled_vertices = []
        for v in self.vertices:
            scaled_vertices.append([i * scaling_coefficient for i in v])
        self.vertices = scaled_vertices
        self.size = size

    def get_edge_vectors(self, norm=False):
        """ Calculate polyhedra edge vectors"""
        edge_vectors = []
        for edge in self.edges:
            v1_index, v2_index = edge
            v1 = np.array(self.vertices[v1_index])
            v2 = np.array(self.vertices[v2_index])
            vec = v2 - v1
            if norm:
                edge_vectors.append(vec / np.linalg.norm(vec))
            else:
                edge_vectors.append(vec)
        self.edge_vectors = edge_vectors

    def rotate_edge(self, edge, angle):
        """ Rotate selected edge of the polyhedra """
        linker = self.edge_linkers[edge]
        linker = linker.rotate(angle, self.edge_vectors[edge])
        linker.center(self.edge_centers[edge])
        self.edge_linkers[edge] = linker
        self.update()

    def build(self, linker, scale='auto', metal=None, bond_length=1.5):
        """ Build polyhedra to generate coordinates """
        if scale is 'auto':
            scale = linker.length + bond_length * 2
        self.linker = linker               # Might not be wise for memory
        self.resize(scale)
        self.get_edge_vectors(norm=True)   # Calculate normalized edge vectors

        self.edge_coors = []
        self.edge_linkers = []
        self.edge_centers = []
        for i, vec in enumerate(self.edge_vectors):
            rotation_axis, angle = align(linker.vector, vec)   # Get rotation axis and angle to align linker
            dest_index_1, dest_index_2 = self.edges[i]    # Get linker destination (mid point of edge)
            dest = (np.array(self.vertices[dest_index_1]) + np.array(self.vertices[dest_index_2])) / 2
            self.edge_centers.append(dest)
            aligned_linker = linker.rotate(angle, rotation_axis)
            aligned_linker.center(dest)
            self.edge_linkers.append(aligned_linker)
            self.edge_coors.append(aligned_linker.atom_coors)

        if metal is not None:
            self.metal = metal
        self.update()

    def add_metal(self, metal='Pd'):
        """ Add metal atoms to vertices """
        bond_length = 1.5
        self.metal = metal
        for v in self.vertices:
            self.atom_names += [metal]
            self.atom_coors += [v]

    def update(self):
        """ Update coordinates and atom names for each linker and metal if exists """
        self.atom_coors = []
        self.atom_names = []
        for l in self.edge_linkers:
            self.atom_coors += l.atom_coors
            self.atom_names += l.atom_names
        if hasattr(self, 'metal'):
            self.add_metal(metal=self.metal)

    def get_force_field(self, ff_path, ff_selection='uff'):
        """ Get force field parameters """
        ff_parameters = read_ff_parameters(ff_path, ff_selection)
        self.ff = dict(type=ff_selection, atom_names=[], sigma=[], epsilon=[])
        for atom_name in self.atom_names:
            self.ff['atom_names'].append(atom_name)
            sig, eps = get_ff_par(atom_name, ff_parameters)
            self.ff['sigma'].append(sig)
            self.ff['epsilon'].append(eps)

    def get_energy(self):
        """ Calculate Lennard-Jones energy for structure """
        min_dist = 1E-5
        self.energy = 0
        for i_1, (name_1, coor_1) in enumerate(zip(self.atom_names, self.atom_coors)):
            for i_2, (name_2, coor_2) in enumerate(zip(self.atom_names, self.atom_coors)):
                if i_2 > i_1:
                    dist = np.linalg.norm(np.array(coor_2) - np.array(coor_1))
                    dist = max(dist, min_dist)  # To avoid getting very large energy values

                    sig_1 = self.ff['sigma'][i_1]
                    sig_2 = self.ff['sigma'][i_2]
                    eps_1 = self.ff['epsilon'][i_1]
                    eps_2 = self.ff['epsilon'][i_2]

                    sig_mix, eps_mix = lb_mix(sig_1, sig_2, eps_1, eps_2)
                    self.energy += lennard_jones(dist, sig_mix, eps_mix)
        return self.energy

    def copy(self):
        """ Return deepcopy of polyhedra """
        return copy.deepcopy(self)

    def save(self, file_format='yaml', save_dir=None, file_name=None, setup=None):
        """ Save polyhedra object """
        if file_name is None:
            file_name = self.name
        if save_dir is None:
            save_dir = os.getcwd()
        save(self, file_format=file_format, file_name=file_name, save_dir=save_dir, setup=setup)

    def get_coordination_vectors(self):
        """ Get coordination vectors for metal atoms """
        coordination = self.symbol[1]
        self.coordination_vectors = []
        for vertice in self.vertices:
            dist_list = []
            for coor in self.atom_coors:
                v = np.array(vertice)
                c = np.array(coor)
                dist = np.linalg.norm(c - v)
                if dist > 1E-6:
                    dist_list.append(dist)
            sorted_dist = sorted(dist_list)
            coord_vec = sorted_dist[:coordination]
            self.coordination_vectors.append([coord_vec])

    def relax_edges(self, angle=15, scan_limit=180, verbose=False):
        """ Rotate each edge and select the configuration with min energy """
        inc = int(scan_limit / angle)
        rot_angles = [math.radians(i * angle) for i in range(1, inc)]
        energies = []
        configurations = []
        for a in rot_angles:
            new_poly = self.copy()
            for i, e in enumerate(new_poly.edges):
                new_poly.rotate_edge(i, a)
            energies.append(new_poly.get_energy())
            configurations.append(new_poly)
            print('Angle: %.1f | Energy: %.1e' % (math.degrees(a), new_poly.energy)) if verbose else None
        min_energy = sorted(energies)[0]
        min_index = energies.index(min_energy)
        min_poly = configurations[min_index]
        print('Selected %.1f rotation' % math.degrees(rot_angles[min_index])) if verbose else None
        return min_poly
