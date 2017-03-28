# 3D vector operations
# Date: March 2017
# Author: Kutay B. Sezginel
import math
import numpy as np


def align(v1, v2):
    """ Calculates the rotation axis and angle to align v1 with v2. """
    rotation_axis = np.cross(v1, v2)
    d = np.dot(v1, v2)
    angle = np.arccos(d)
    if math.isnan(angle):
        angle = 0
    return rotation_axis, angle
