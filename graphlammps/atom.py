#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from graphlammps import params

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

    def sample_vel_Temp(self, Temp, IP='reaxFF-CHO'):
        m   = self.mass/1000/params.NA
        sig = (params.kb*Temp/m)**0.5
        
        if('reaxFF' in IP): unit_corr = 1E-5
        else              : unit_corr = 1E-2
        
        v = np.zeros_like(self.vel)
        v = np.random.normal(0, sig, v.shape).reshape(v.shape)
        self.vel = np.copy(v*unit_corr)
        