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
