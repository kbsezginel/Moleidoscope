import os
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="moleidoscope",
    version="0.1",
    description="Molecular structure generator",
    author="Kutay B. Sezginel",
    author_email="kbs37@pitt.edu",
    install_requires=requirements,
    dependency_links=['http://github.com/kbsezginel/HostDesigner/tarball/master#egg=package-1.0'],
    packages=find_packages(),
    include_package_data=True
)

try:
    hd_dir = os.environ['HD_DIR']
    library_path = os.path.join(hd_dir, 'LIBRARY')
    if not os.path.exists(library_path):
        print('Linker library not found in HostDesigner directory!')
        print('Please add LIBRARY to %s' % hd_dir)
    else:
        print('Using LIBRARY in %s' % library_path)
except:
    print('HostDesigner directory not found in environment variables!')
    print('Please add HostDesigner directory (including LIBRARY file) as HD_DIR variable')
