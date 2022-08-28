#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 15:08:49 2022

@author: chaithanya
"""
import numpy as np
import sys


class molecule:
    def __init__(self, name):
        self.name = name
        self.timestep_prod = 0
        self.mol_id = 0
        self.atom_id = []
