# Author: Kutay B. Sezginel
# Date: August 2017
"""
Animation of molecules
"""
import os
import mdtraj
import nglview
import tempfile
from .output import write_pdb


def animate(frames, gui=False, delete=True,):
    """
    Creates nglview widget for given list of molecule files (frames).
    """
    T = mdtraj.load(frames, top=frames[0])
    view = nglview.show_mdtraj(T, gui=gui)
    view.add_ball_and_stick()
    if delete:
        for frame in frames:
            os.remove(frame)
    return view


def rotate(molecule, angle=2, axis=[1, 0, 0], n_frames=100, temp_dir=None):
    """
    Rotate molecules in given axis, angle increment and number of steps.
    """
    angle = np.deg2rad(angle)
    frames = []
    for frame in range(1, n_frames + 1):
        mol = molecule.rotate(angle * frame, axis)
        temp_pdb_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False)
        write_pdb(temp_pdb_file, mol.atom_names, mol.atom_coors)
        frames.append(temp_pdb_file.name)
    return frames
