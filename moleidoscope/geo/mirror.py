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


class Mirror:
    """ Mirror class."""
    def __init__(self, input):
        self.a = input[0]
        self.b = input[1]
        self.c = input[2]
        self.d = input[3]
