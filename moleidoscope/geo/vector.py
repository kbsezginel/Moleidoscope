# Date: March 2017
# Author: Kutay B. Sezginel
"""
3D vector operations
"""
import math
import numpy as np


def align(v1, v2, norm=True):
    """ Calculates the rotation axis and angle to align v1 with v2. """
    if norm:
        v1 = np.array(v1) / np.linalg.norm(v1)
        v2 = np.array(v2) / np.linalg.norm(v2)
    rotation_axis = np.cross(v1, v2)
    d = np.dot(v1, v2)
    angle = np.arccos(d)
    if math.isnan(angle):
        angle = 0
    return rotation_axis, angle


def find_closest(target, coor_list):
    """ Find closest coordinate to a set of coordinates """
    target = np.array(target)
    distances = []
    for coor in coor_list:
        c = np.array(coor)
        distances.append(np.linalg.norm(target - c))
    sorted_distances = sorted(distances)
    min_d = sorted_distances[0]
    dist_index = distances.index(min_d)
    return coor_list[dist_index]
