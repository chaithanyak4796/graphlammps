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
        self.name           = name
        self.atoms_list     = []
        self.num_atom_types = 0
        self.num_atoms      = 0
        self.ortho_box      = True

    def get_num_atom_types(self):
        atype = []
        for at in self.atoms_list:
            atype.append(at.type)
        self.num_atom_types = max(atype)
    
    def renumber_atoms(self):
        idx = 1
        for at in self.atoms_list:
            at.idx = idx
            idx   += 1
        self.num_atoms = len(self.atoms_list)
        self.get_num_atom_types()

    def get_atom(self, idx):
        for at in self.atoms_list:
            if(at.idx == idx):
                return at
        sys.exit(f"Atom with atom id {idx} not found in the atoms_list.")