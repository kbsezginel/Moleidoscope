# Date: May 2017
# Authors: Kutay B. Sezginel
"""
File input methods (formats: xyz)
"""
import os
import numpy as np


def read_xyz(xyz_path):
    """
    Reads xyz file and return dictionary as:
        - name:       Name of the structure
        - n_atoms:    Number of atoms
        - atom_names: Atoms names as list
        - Atom_coors: Atom coordinates as list of 3D lists
    """
    with open(xyz_path, 'r') as xyz_file:
        xyz_lines = xyz_file.readlines()
    n_atoms = int(xyz_lines[0].strip())
    name = xyz_lines[1].strip()
    atom_names = []
    atom_coors = []
    for line in xyz_lines[2:]:
        atom = line.split()[0]
        x, y, z = line.split()[1:4]
        atom_names.append(atom)
        atom_coors.append(np.array([float(x), float(y), float(z)]))
    return dict(name=name, n_atoms=n_atoms,
                atom_names=atom_names,
                atom_coors=atom_coors)
