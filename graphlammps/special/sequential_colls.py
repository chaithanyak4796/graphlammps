#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 16:25:48 2022

@author: ckondur
"""

import numpy as np
import sys

class product:
    def __init__(self, name):
        self.name      = name
        self.count_old = 0
        self.count_new = 0
        
        self.timestep_prod = []
        
def parse_species(spec_info, spec_count):
    
    for i in range(len(spec_info)):
        name  = spec_info[i][0]
        count = spec_info[i][1]
   
        if 'H' in name:   # lattice
            if 'O' in name:
                c = name.split('O')
                # print (c)
                if (c[-1] == ''): count = 1
                else            : count = int(c[-1])
                name  = 'O_ads'
        
        if name in spec_count:
            spec_count[name] = count