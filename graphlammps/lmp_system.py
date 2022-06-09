#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import random
import sys

# from graphlammps import create, atom
# from graphlammps.params import mass_C, mass_O, mass_N, kb, NA

class lmp_system:
    """
        This class contains the parameters of the lammps system.
        It will be used both pre and post processing.
    """
    
    def __init__(self, name='general'):
        self.name = name
        self.atoms_list = []

    