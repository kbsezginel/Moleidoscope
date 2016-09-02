import math
import numpy as np


class Coor(object):
    """
    Coor class for holding 3D space coordinates.
    """
    def __init__(self, input):
        self.x = input[0]
        self.y = input[1]
        self.z = input[2]

    def __repr__(self):
        return "<Coordinate object x:%s y:%s z:%s>" % (self.x, self.y, self.z)

    def __str__(self):
        return "%s %s %s" % (round(self.x, 4), round(self.y, 4), round(self.z, 4))

    def __add__(self, coor2):
        return Coor([self.x + coor2.x, self.y + coor2.y, self.z + coor2.z])

    def __sub__(self, coor2):
        return Coor([self.x - coor2.x, self.y - coor2.y, self.z - coor2.z])

    def dist(self, coor2):
        """
        Calculates distance between this coordinate and another given coordinate.

        Example usage:
         >>> coor1 = Coor([1, 2, 3])
         >>> coor2 = Coor([-3, 4, -2])
         >>> coor1.dist(coor2) -> 6.708203932499369
        """
        dx = self.x - coor2.x
        dy = self.y - coor2.y
        dz = self.z - coor2.z

        return math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    def xyz(self):
        """
        Return a list containing x, y, z coordinates.

        Example usage:
         >>> coor1 = Coor([1, 2, 3])
         >>> coor1.xyz() -> [1, 2, 3]
        """
        return [self.x, self.y, self.z]

    def frac(self, uc_size, uc_angle, frac_ucv):
        """
        Converts cartesian coordinates to fractional coordinates.
        *The fractional unit cell volume is calculated each time. Instead can be given as input.

        Example usage:
         >>> coor1 = Coor([-11, 22, 33])
         >>> coor1.frac([26, 26, 26], [90, 90, 90]) -> <Coordinate object x:-0.423 y:0.846 z:1.269>
        """
        alp = math.radians(uc_angle[0])
        bet = math.radians(uc_angle[1])
        gam = math.radians(uc_angle[2])

        v = frac_ucv
        # v = 1 - math.cos(alp) ** 2 - math.cos(bet) ** 2
        # v += - math.cos(gam) ** 2 + 2 * math.cos(alp) * math.cos(bet) * math.cos(gam)
        # v = math.sqrt(v)

        a = uc_size[0]
        b = uc_size[1]
        c = uc_size[2]

        x = self.x
        y = self.y
        z = self.z

        x_frac = 1 / a * x
        x_frac += - math.cos(gam) / (a * math.sin(gam)) * y
        x_frac += (math.cos(alp) * math.cos(gam) - math.cos(bet)) / (a * v * math.sin(gam)) * z

        y_frac = 1 / (b * math.sin(gam)) * y
        y_frac += (math.cos(bet) * math.cos(gam) - math.cos(alp)) / (b * v * math.sin(gam)) * z

        z_frac = math.sin(gam) / (c * v) * z

        return Coor([x_frac, y_frac, z_frac])

    def car(self, uc_size, uc_angle, frac_ucv):
        """
        Converts fractional coordinates to cartesian coordinates.
        Takes unit cell size, angles and fractional unit cell volume as input.

        Example usage:
         >>> coor1 = Coor([0.23, 2.21, -1.33])
         >>> coor1.car([26, 26, 26], [90, 90, 90]) -> <Coordinate object x:5.98 y:57.46 z:-34.58>
        """
        alp = math.radians(uc_angle[0])
        bet = math.radians(uc_angle[1])
        gam = math.radians(uc_angle[2])

        v = frac_ucv
        # v = 1 - math.cos(alp) ** 2 - math.cos(bet) ** 2
        # v += - math.cos(gam) ** 2 + 2 * math.cos(alp) * math.cos(bet) * math.cos(gam)
        # v = math.sqrt(v)

        a = uc_size[0]
        b = uc_size[1]
        c = uc_size[2]

        x_frac = self.x
        y_frac = self.y
        z_frac = self.z

        x = a * x_frac
        x += b * math.cos(gam) * y_frac
        x += c * math.cos(bet) * z_frac

        y = b * math.sin(gam) * y_frac
        y += c * (math.cos(alp) - math.cos(bet) * math.cos(gam)) / math.sin(gam) * z_frac

        z = c * v / math.sin(gam) * z_frac

        return Coor([x, y, z])

    def pbc(self, uc_size, uc_angle, frac_ucv):
        """
        Apply perodic boundary conditions to given cartesian coordinates and unit cell parameters.

        Example usage:
         >>> coor1 = Coor([0.23, 2.21, -1.33])
         >>> coor1.pbc([26, 26, 26], [90, 90, 90]) -> <Coordinate object x:0.23 y:0.21 z:0.67>
        """
        frac_coor = self.frac(uc_size, uc_angle, frac_ucv)
        frac_pbc_coor = frac_coor.frac_pbc()
        car_pbc_coor = frac_pbc_coor.car(uc_size, uc_angle, frac_ucv)

        return car_pbc_coor

    def frac_pbc(self):
        """
        Apply perodic boundary conditions to given fractional coordinates.

        Example usage:
         >>> coor1 = Coor([0.23, 2.21, -1.33])
         >>> coor1.frac_pbc() -> <Coordinate object x:0.23 y:0.21 z:0.67>
        """
        pbc_x = self.x - math.floor(self.x)
        pbc_y = self.y - math.floor(self.y)
        pbc_z = self.z - math.floor(self.z)

        return Coor([pbc_x, pbc_y, pbc_z])


class Mirror:
    """ Mirror class."""
    def __init__(self, *args, size=10):
        if len([*args]) == 3:
            p1 = np.array(p1)
            p2 = np.array(p2)
            p3 = np.array(p3)
        elif str(*args) == 'xy':
            p1 = np.array([0, 0, 0])
            p2 = np.array([size, 0, 0])
            p3 = np.array([size, size, 0])
        elif str(*args) == 'xz':
            p1 = np.array([0, 0, 0])
            p2 = np.array([size, 0, 0])
            p3 = np.array([size, 0, size])
        elif str(*args) == 'yz':
            p1 = np.array([0, 0, 0])
            p2 = np.array([0, 0, size])
            p3 = np.array([0, size, size])
        # This can fail for non-string inputs
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
        self.p1, self.p2, self.p3 = p1, p2, p3

    def symmetry(self, coor):
        """ Get symmetrical points through the mirror. """
        s0 = (self.a * coor.x + self.b * coor.y + self.c * coor.z + self.d)
        s0 /= (self.a**2 + self.b**2 + self.c**2)

        x = coor.x - 2 * s0 * self.a
        y = coor.y - 2 * s0 * self.b
        z = coor.z - 2 * s0 * self.c

        return Coor([x, y, z])

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
