# Author: Kutay B Sezginel
# Date: August 2017
"""
Geometric shapes library.
In elementary geometry, a polytope is a geometric object with "flat" sides.
"""
import os
import yaml


lib_dir = os.path.abspath(os.path.dirname(__file__))


def read_yaml(file_name):
    with open(os.path.join(lib_dir, file_name), 'r') as f:
        return yaml.load(f)


polytope = {}
polytope['cube'] = read_yaml('cube.yaml')
polytope['triangle'] = read_yaml('triangle.yaml')
polytope['octahedron'] = read_yaml('octahedron.yaml')
