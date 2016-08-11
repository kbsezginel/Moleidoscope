# Python library for 3 dimensional geometric operations
# Quaternion class:
# - Multiplication  - Division  - Inverse   - Rotation  - Output List <Q.xyz()>
#  >>> Q = Quaternion([0, 1, 1, 1])
# Coor class:
#  - Addition  - Subtraction - Distance  - To fractional - To cartesian
#  - Periodic boundary conditions
#  >>> coor = Coor([1, 2, 3])
# Date: June 2016
# Author: Kutay B. Sezginel
import math


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


class Quaternion(object):
    """
    Quaternion class for quaternion operations and 3D rotations.
    """
    def __init__(self, input):
        self.w = input[0]
        self.x = input[1]
        self.y = input[2]
        self.z = input[3]

    def __repr__(self):
        return "<Quaternion object w:%s x:%s y:%s z:%s>" % (self.w, self.x, self.y, self.z)

    def __str__(self):
        return "x:%s y:%s z:%s" % (self.x, self.y, self.z)

    def xyz(self):
        """
        Returns x, y, z values of the quaternion in list format.
        """
        return [self.x, self.y, self.z]

    def __mul__(self, quat2):
        """
        Multiply quaternion by another.

        Example usage::
          >>> q1 = Quaternion([1, 2, 3, 4])
          >>> q2 = Quaternion([2, 3, 4, 5])
          >>> q1 * q2 -> <Quaternion object w:-36 x:6 y:12 z:12>
        """

        q1 = self
        q2 = quat2

        w3 = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
        x3 = q1.x * q2.w + q1.w * q2.x - q1.z * q2.y + q1.y * q2.z
        y3 = q1.y * q2.w + q1.z * q2.x + q1.w * q2.y - q1.x * q2.z
        z3 = q1.z * q2.w - q1.y * q2.x + q1.x * q2.y + q1.w * q2.z

        return Quaternion([w3, x3, y3, z3])

    def __truediv__(self, quat2):
        """
        Divide one quaternion by another. Performs the operation as q1 * inverse q2.

        Example usage::
          >>> q1 = Quaternion([1, 2, 3, 4])
          >>> q2 = Quaternion([2, 3, 4, 5])
          >>> q1 / q2 -> <Quaternion object w:0.7407 x:0.0370 y:0.0 z:0.0741>
        """
        return self * quat2.inv()

    def inv(self):
        """
        Returns the inverse of the quaternion as a new quaternion.

        """
        norm = self.w ** 2 + self.x ** 2 + self.y ** 2 + self.z ** 2

        return Quaternion([self.w / norm, -self.x / norm, -self.y / norm, -self.z / norm])

    def rotation(self, rotation_point, axis_point1, axis_point2, rotation_angle):
        """
        Rotation of a point around an axis defined by two points in 3D space.
        Rotation angle needs to be given in radians.

        Example usage::
         >>> Q = Quaternion([0, 1, 1, 1])
         >>> Q = Q.rotation(Q.xyz(), [-2, 4, 6.1], [0.3, 1.2, -0.76], math.pi/6)
         >>> [2.1192250600275795, 2.2773560513200133, 5.890236840657188]
        """
        i = axis_point2[0] - axis_point1[0]
        j = axis_point2[1] - axis_point1[1]
        k = axis_point2[2] - axis_point1[2]
        length = math.sqrt(i**2 + j**2 + k**2)
        i = i / length
        j = j / length
        k = k / length
        qp_w = 0
        qp_x = rotation_point[0] - axis_point2[0]
        qp_y = rotation_point[1] - axis_point2[1]
        qp_z = rotation_point[2] - axis_point2[2]
        Q_point = Quaternion([qp_w, qp_x, qp_y, qp_z])

        qr_w = math.cos(rotation_angle / 2.0)
        qr_x = math.sin(rotation_angle / 2.0) * i
        qr_y = math.sin(rotation_angle / 2.0) * j
        qr_z = math.sin(rotation_angle / 2.0) * k
        Q_rot = Quaternion([qr_w, qr_x, qr_y, qr_z])

        Quat = (Q_rot * Q_point) * Q_rot.inv()
        Quat.x = Quat.x + axis_point2[0]
        Quat.y = Quat.y + axis_point2[1]
        Quat.z = Quat.z + axis_point2[2]

        return Quat

    def coor(self):
        """
        Converts Quaternion object to Coor object.
        """
        return Coor(self.xyz())
