# Moleidoscope visualization methods
# Author: Kutay B Sezginel
# Date: February 2017
import os
import math
import nglview

molecule_export_dir = os.path.join(os.getcwd(), 'doc', 'tmp')


def show(*args, camera='perspective', move='auto', div=5, distance=(-10, -10), axis=0, rename='F'):
    """
    Show given structures using nglview
        - camera: 'perspective' / 'orthographic'
        - move: separate multiple structures equally distant from each other
        - distance: separation distance
        - axis: separation direction
        - rename: rename dummy atoms to given atom name for visualization
    """
    if move is 'auto':
        translation_vectors = arrange_structure_positions(len(args), div=div, distance=distance, rename='F')
    elif move is 'single':
        translation_vectors = axis_translation(len(args), distance=distance[0], axis=axis)
    else:
        translation_vectors = [[0, 0, 0]] * len(args)
    atom_names = []
    atom_coors = []
    for molecule, vec in zip(args, translation_vectors):
        atom_names += molecule.atom_names
        atom_coors += translate(molecule.atom_coors, vector=vec)
    pdb_path = os.path.join(molecule_export_dir, 'view.pdb')
    export_pdb(pdb_path, atom_names, atom_coors, rename=rename)
    view = nglview.show_structure_file(pdb_path)
    view.camera = camera
    return view


def arrange_structure_positions(n_structures, div=5, distance=(10, 10), rename='F'):
    """
    Arrange structure positions according to number of structures given.
    """
    n_structures_lateral = div
    split = math.ceil(n_structures / n_structures_lateral)
    vertical_positions = axis_translation(split, distance=distance[1], axis=1)
    translation_vectors = []
    for s in range(split):
        if s < 2:
            n_split = math.ceil(n_structures / split)
        else:
            n_split = math.floor(n_structures / split)
        new_vectors = axis_translation(n_split, distance=distance[0], axis=0)
        new_vectors = translate(new_vectors, vector=vertical_positions[s])
        translation_vectors += new_vectors
    return translation_vectors


def translate(atom_coors, vector=[-10, 0, 0]):
    """ Translate given coordinates with given vector """
    translated_coors = []
    x, y, z = vector
    for coor in atom_coors:
        new_coor = [coor[0] + x, coor[1] + y, coor[2] + z]
        translated_coors.append(new_coor)
    return translated_coors


def axis_translation(n_structures, distance=-10, axis=0):
    """
    Automatically adjust structure positions equally distant from each other in given axis
        - distance: distance between each structure
        - axis: axis selection for translation (0: x-axis, 1: y-axis, 2: z-axis)
    """
    translation_vectors = []
    lim = (n_structures - 1) * distance / 2
    for i in range(n_structures):
        vec = [0, 0, 0]
        vec[axis] = -lim + i * distance
        translation_vectors.append(vec)
    return translation_vectors


def export_pdb(pdb_path, names, coors, rename=False):
    """ Export given atomic coordinates to pdb format """
    replace_list = ['Du', 'DU', 'X']
    if len(rename) > 1:
        replace_list += [rename[0]]
        rename = rename[1]
    structure_name = os.path.splitext(os.path.basename(pdb_path))[0]
    with open(pdb_path, 'w') as pdb_file:
        pdb_file.write('HEADER    ' + structure_name + '\n')
        format = 'HETATM%5d%3s  MOL     1     %8.3f%8.3f%8.3f  1.00  0.00          %2s\n'
        atom_index = 1
        for atom_name, atom_coor in zip(names, coors):
            if rename is not False:
                if atom_name in replace_list:
                    atom_name = rename
            x, y, z = atom_coor
            pdb_file.write(format % (atom_index, atom_name, x, y, z, atom_name.rjust(2)))
            atom_index += 1
        pdb_file.write('END\n')


def write_pdb(pdb_path, names, coors):
    structure_name = os.path.splitext(os.path.basename(pdb_path))[0]
    with open(pdb_path, 'w') as pdb_file:
        pdb_file.write('HEADER    ' + structure_name + '\n')
        format = 'HETATM%5d%3s  MOL     1     %8.3f%8.3f%8.3f  1.00  0.00          %2s\n'
        atom_index = 1
        for atom_name, atom_coor in zip(names, coors):
            x, y, z = atom_coor
            pdb_file.write(format % (atom_index, atom_name, x, y, z, atom_name.rjust(2)))
            atom_index += 1
        pdb_file.write('END\n')
