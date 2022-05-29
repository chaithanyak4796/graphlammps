#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class atom:
    """
        This is the class for any atom
    """
    def __init__(self, mass=12.0107, idx=1, a_type=1, a_name='C'):
        self.pos  = np.zeros((3,), dtype=float)
        self.vel  = np.zeros_like(self.pos)
        self.type = a_type
        self.name = a_name
        self.idx  = idx
        self.mass = mass
        self.q    = 0.0
        self.corner  = False
        self.layer   = 0
        self.lattice = True

