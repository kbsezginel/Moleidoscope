# Date: August 2016
# Authors: Kutay B. Sezginel
"""
Mirror object to create prisms and perform reflections
"""
import numpy as np
from moleidoscope.geo.quaternion import Quaternion


class Mirror:
    """ Mirror class."""
    def __init__(self, *args, size=10):
        if len([*args]) == 3:
            p1 = np.array(args[0])
            p2 = np.array(args[1])
            p3 = np.array(args[2])
            p1, p2, p3 = p1 * size, p2 * size, p3 * size    # Change the size of the mirror plane
            self.name = 'p1: ' + str(p1) + ' p2: ' + str(p2) + ' p3: ' + str(p3)
        elif str(*args) == 'xy':
            p1 = np.array([0, 0, 0])
            p2 = np.array([size, 0, 0])
            p3 = np.array([size, size, 0])
            self.name = str(*args)
        elif str(*args) == 'xz':
            p1 = np.array([0, 0, 0])
            p2 = np.array([size, 0, 0])
            p3 = np.array([size, 0, size])
            self.name = str(*args)
        elif str(*args) == 'yz':
            p1 = np.array([0, 0, 0])
            p2 = np.array([0, 0, size])
            p3 = np.array([0, size, size])
            self.name = str(*args)
        # Source: http://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/
        # These two vectors are in the plane
        self.v1 = p3 - p1
        self.v2 = p2 - p1
        # the cross product is a vector normal to the plane
        cp = np.cross(self.v1, self.v2)
        self.a, self.b, self.c = cp
        # This evaluates a * x3 + b * y3 + c * z3 which equals d
        self.d = np.dot(cp, p3)
        self.p1, self.p2, self.p3 = p1, p2, p3

    def symmetry(self, coor):
        """ Get symmetrical points through the mirror. """
        s0 = (self.a * coor[0] + self.b * coor[1] + self.c * coor[2] - self.d)
        s0 /= (self.a**2 + self.b**2 + self.c**2)

        x = coor[0] - 2 * s0 * self.a
        y = coor[1] - 2 * s0 * self.b
        z = coor[2] - 2 * s0 * self.c

        return [x, y, z]

    def rotate(self, axis, angle, size=1):
        """ Rotate mirror around an axis.
            - axis -> ex: [[1, 0, 0], [0, 0, 0]]
            - angle -> ex: math/pi / 2 (must be in radians)
            - size -> default is 1 which means same size as before
        """
        q = Quaternion([0, 1, 1, 1])
        p1r = q.rotation(self.p1, axis[0], axis[1], angle).xyz()
        p2r = q.rotation(self.p2, axis[0], axis[1], angle).xyz()
        p3r = q.rotation(self.p3, axis[0], axis[1], angle).xyz()
        new_mirror = Mirror(p1r, p2r, p3r, size=size)
        return new_mirror

    def coordinates(self, length, grid_size=20):
        """ Get meshgrid coordinates for the mirror plane. """
        nx, ny = (grid_size, grid_size)
        x = np.linspace(0, length, nx)
        y = np.linspace(0, length, ny)
        xv, yv = np.meshgrid(x, y)

    def grid_plane(self, grid_size):
        """ Three points P1(x1, y1, z1), P2(x2, y2, z2), P3(x3, y3, z3)
            Make sure this is the same with the plane defined using 3 points!!!
                  ________
              P2 /   /    / P3
             Pn /---/----/
            P1 /___/____/
                  Pm
        """
        p1, p2, p3 = self.p1, self.p2, self.p3
        grid_plane_coors = []
        for m in range(grid_size + 1):
            for n in range(grid_size + 1):
                x = (grid_size - n) / grid_size * p1[0] + (n - m) / grid_size * p2[0] + m / grid_size * p3[0]
                y = (grid_size - n) / grid_size * p1[1] + (n - m) / grid_size * p2[1] + m / grid_size * p3[1]
                z = (grid_size - n) / grid_size * p1[2] + (n - m) / grid_size * p2[2] + m / grid_size * p3[2]
                grid_plane_coors.append([x, y, z])

        return grid_plane_coors

    def to_linker(self, grid_size=10, show_edges=False, atom_type='N'):
        """ Convert mirror class to linker """
        mirror_coors = self.grid_plane(grid_size)
        mirror_names = [atom_type] * len(mirror_coors)
        if show_edges:
            mirror_coors += [self.p1, self.p2, self.p3]
            mirror_names += ['S', 'S', 'S']
        from moleidoscope.linker import Linker
        l_mirror = Linker()
        l_mirror.name = self.name
        l_mirror.atom_coors = mirror_coors
        l_mirror.atom_names = mirror_names
        return l_mirror

    def get_center(self):
        return self.p1 + (self.p3 - self.p1) / 2

    def scale(self, size=5):
        self.p1 *= size
        self.p2 *= size
        self.p3 *= size
