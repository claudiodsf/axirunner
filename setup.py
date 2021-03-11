# -*- coding: utf-8 -*-
"""setup.py: setuptools control."""
from setuptools import setup

import inspect
import os
import sys

# Import the version string.
path = os.path.join(os.path.abspath(os.path.dirname(inspect.getfile(
    inspect.currentframe()))), 'axirunner')
sys.path.insert(0, path)
from version import get_git_version


with open('README.md', 'rb') as f:
    long_descr = f.read().decode('utf-8')


setup(
    name='axirunner',
    packages=['axirunner', 'axirunner.configobj'],
    include_package_data=True,
    entry_points={
        'console_scripts': ['axirunner = axirunner.axirunner:main']
        },
    version=get_git_version(),
    description='Make synthetic seismograms using Axitra. The easy way™.',
    long_description=long_descr,
    long_description_content_type='text/markdown',
    author='Claudio Satriano',
    author_email='satriano@ipgp.fr',
    url='http://www.ipgp.fr/~satriano',
    license='CeCILL Free Software License Agreement, Version 2.1',
    platforms='OS Independent',
    classifiers=[
            'Development Status :: Development Status :: 2 - Pre-Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre '
                'License, version 2.1 (CeCILL-2.1)',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics'],
    install_requires=['obspy>=1.2.0', 'scipy>=0.17', 'proj'],
    python_requires='>3'
    )
