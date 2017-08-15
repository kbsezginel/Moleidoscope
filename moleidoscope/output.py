# Date: March 2017
# Authors: Kutay B. Sezginel
"""
File output methods (formats: pdb / xyz / yaml / orca)
"""
import os
import yaml


def save(molecule, file_name='mol', file_format='yaml', save_dir=None, setup=None):
    """ Save object to selected format """
    file_path = os.path.join(save_dir, '%s.%s' % (file_name, file_format))
    with open(file_path, 'w') as file_object:
        if file_format is 'yaml':
            yaml.dump(molecule, file_object)
        elif file_format is 'xyz':
            write_xyz(file_object, molecule.atom_names, molecule.atom_coors, header=file_name)
        elif file_format is 'pdb':
            write_pdb(file_object, molecule.atom_names, molecule.atom_coors, header=file_name)
        elif file_format is 'orca':
            write_orca(file_object, molecule.atom_names, molecule.atom_coors, header=file_name, setup=setup)
    return file_path


def write_pdb(pdb_file, names, coors, header='mol'):
    """ Write given atomic coordinates to file object in pdb format """
    pdb_file.write('HEADER    ' + header + '\n')
    format = 'HETATM%5d%3s  MOL     1     %8.3f%8.3f%8.3f  1.00  0.00          %2s\n'
    for atom_index, (atom_name, atom_coor) in enumerate(zip(names, coors), start=1):
        x, y, z = atom_coor
        pdb_file.write(format % (atom_index, atom_name, x, y, z, atom_name.rjust(2)))
    pdb_file.write('END\n')
    pdb_file.flush()


def write_xyz(xyz_file, names, coors, header='mol'):
    """ Write given atomic coordinates to file object in xyz format """
    xyz_file.write(str(len(coors)) + '\n')
    xyz_file.write(header + '\n')
    format = '%s %.4f %.4f %.4f\n'
    for atom, coor in zip(names, coors):
        xyz_file.write(format % (atom, coor[0], coor[1], coor[2]))
    xyz_file.flush()


def write_orca(orca_file, names, coors, header='mol', setup=None):
    """ Write given coordinates to file object in orca input format """
    with open(orca_setup_path, 'r') as orca_setup:
        orca_setup_lines = orca_setup.readlines()
    for line in orca_setup_lines:
        orca_file.write(line)
    orca_file.write('\n')
    format = '%s %.4f %.4f %.4f\n'
    for atom, coor in zip(names, coors):
        xyz_file.write(format % (atom, coor[0], coor[1], coor[2]))
    orca_file.write('\n*')
    orca_file.flush()
