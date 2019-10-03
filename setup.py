#!/usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup

setup(name='topoflow',
      version='3.5',
      description='d8-based, spatial hydrologic model',
      author='Scott D. Peckham',
      author_email='Scott.Peckham@colorado.edu',
      license='MIT',
      url='https://github.com/peckhams/topoflow',
      # url='http://csdms.colorado.edu/wiki/Model:TopoFlow',
      packages=['topoflow',
                'topoflow.components',
                'topoflow.components.tests',
                'topoflow.examples',
                'topoflow.framework',
                'topoflow.framework.tests',
                'topoflow.gui',
                'topoflow.utils',
                'topoflow.utils.tests'], 
      install_requires=['numpy', 'scipy', 'h5py', 'netCDF4', 'cfunits'],
     )
