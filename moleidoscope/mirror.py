import numpy as np


class Mirror:
    """ Mirror class."""
    def __init__(self, *args, size=10):
        if len([*args]) == 3:
            p1 = np.array(args[0])
            p2 = np.array(args[1])
            p3 = np.array(args[2])
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
        v1 = p3 - p1
        v2 = p2 - p1
        # the cross product is a vector normal to the plane
        cp = np.cross(v1, v2)
        self.a, self.b, self.c = cp
        # This evaluates a * x3 + b * y3 + c * z3 which equals d
        self.d = np.dot(cp, p3)
        self.p1, self.p2, self.p3 = p1 * size, p2 * size, p3 * size

    def symmetry(self, coor):
        """ Get symmetrical points through the mirror. """
        s0 = (self.a * coor[0] + self.b * coor[1] + self.c * coor[2] + self.d)
        s0 /= (self.a**2 + self.b**2 + self.c**2)

        x = coor[0] - 2 * s0 * self.a
        y = coor[1] - 2 * s0 * self.b
        z = coor[2] - 2 * s0 * self.c

        return [x, y, z]

    def coordinates(self, length, grid_size=20):
        """ Get meshgrid coordinates for the mirror plane. """
        nx, ny = (grid_size, grid_size)
        x = np.linspace(0, length, nx)
        y = np.linspace(0, length, ny)
        xv, yv = np.meshgrid(x, y)

    def grid_plane(self, grid_size):
        """ Three points P1(x1, y1, z1), P2(x2, y2, z2), P3(x3, y3, z3)
                  ________
              P2 /   /    / P3
             Pn /---/----/
            P1 /___/____/
                  Pm
        """
        p1, p2, p3 = self.p1, self.p2, self.p3
        grid_plane_coors = []

        for m in range(grid_size):
            for n in range(grid_size):
                x = (grid_size - n) / grid_size * p1[0] + (n - m) / grid_size * p2[0] + m / grid_size * p3[0]
                y = (grid_size - n) / grid_size * p1[1] + (n - m) / grid_size * p2[1] + m / grid_size * p3[1]
                z = (grid_size - n) / grid_size * p1[2] + (n - m) / grid_size * p2[2] + m / grid_size * p3[2]
                grid_plane_coors.append([x, y, z])

        return grid_plane_coors

    def to_linker(self, l_axis, grid_size=10, atom_type='N'):
        mirror_coors = self.grid_plane(grid_size)
        mirror_names = [atom_type] * len(mirror_coors)
        l_axis.name = self.name
        l_axis.atom_coors = mirror_coors
        l_axis.atom_names = mirror_names
        return l_axis

    def scale(self, size=5):
        self.p1 *= size
        self.p2 *= size
        self.p3 *= size
