#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import random
import sys
from typing import List

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

    def get_atom_idx(self, idx: int):
        if idx in self.idx_map:
            return self.idx_map[idx]
        else:
            raise Exception(f'Error in system.get_atom_idx(); No atom index with {idx} found.')
    
    def get_atom(self, idx : int):
        return self.atoms_list[self.get_atom_idx(idx)]
    
    
    def remove_atoms(self, idx : List[int]):
        for i in sorted(idx, reverse=True):
            self.atoms_list.pop(i)