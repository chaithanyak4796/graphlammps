### This file installs the graphlammps library ####
from setuptools import find_packages, setup

setup(
      name        = 'graphlammps',
      packages    = find_packages(),
      version     = '0.1.0',
      description = 'Python library to pre-process and post-process LAMMPS simulations of graphitic systems',
      author      = 'Chaithanya Kondur',
      license     = 'MIT',
      )