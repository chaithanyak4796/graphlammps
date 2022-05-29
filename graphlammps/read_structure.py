#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 11:34:43 2022

@author: ckondur
"""

import numpy as np
import sys
import subprocess

from graphlammps import lmp_system, atom

class read_structure:
    """ This class contains the functions to read the various formats of lammps structure outputs.
        Currently supported : .xyz, .dump
    
    """
    def __init__(self, fname, cols):
        self.num_timesteps = 1
        self.fname = fname
        self.cols  = cols
        self.default_cols = ['id', 'type', 'element', 'q', 'x', 'y', 'z']
        
        if(self.cols != self.default_cols):
            sys.exit("The column mapping provided is currently not supported.")
        
        if(".dump" in self.fname):
            self.read_next_timestep = self.read_time_step_dump
            self.init_read_dump()
        else:
            sys.exit(" Invalid file extension. ")
            
    def init_read_dump(self):
        try:
            self.fr = open(self.fname, "r")
        except:
            print(" Unable to open file : ", self.fname)
            sys.exit(1)
        
        # Find the total number of time steps in the dump file
        cmd = ["grep", "TIMESTEP", self.fname]
        res = subprocess.run(cmd, stdout = subprocess.PIPE)
        self.num_timesteps = len(res.stdout)//len('ITEM: TIMESTEP\n')
        
        print("Total number of timesteps in file = ", self.num_timesteps)
        
    def read_time_step_dump(self):
        """ This function reads the next time step in the dump file and returns a lmp_system object """
        
        system = lmp_system.lmp_system()
        
        line  = self.fr.readline()  # ITEM: TIMESTEP
        line  = self.fr.readline()
        system.timestep = int(line)
        
        line = self.fr.readline() # ITEM: NUMBER OF ATOMS
        line = self.fr.readline()
        system.num_atoms = int(line)
        
        line   = self.fr.readline().split() # ITEM: BOX BOUNDS
        line_1 = self.fr.readline().split()
        line_2 = self.fr.readline().split()
        line_3 = self.fr.readline().split()
        
        if(line[3] == 'xy'): # non-ortho box
            system.xyz = np.zeros((3,), dtype=float)
            
            system.xlo_bound = float(line_1[0])
            system.xhi_bound = float(line_1[1])
            system.xyz[0]    = float(line_1[2])
            
            system.ylo_bound = float(line_2[0])
            system.yhi_bound = float(line_2[1])
            system.xyz[1]    = float(line_2[2])
            
            system.zlo_bound = float(line_3[0])
            system.zhi_bound = float(line_3[1])
            system.xyz[2]    = float(line_3[2])
            
        else: # ortho_box
            system.xlo, system.xhi = float(line_1)
            system.ylo, system.yhi = float(line_2)
            system.zlo, system.zhi = float(line_3)
            
            
        line = self.fr.readline().split() # ITEM: ATOMS
        line = line[2:]
        
        if (line != self.cols):
            sys.exit("Coloumns in file do not match the given mapping")
    
        for i in range(system.num_atoms):
            system.atoms_list.append(atom.atom())
            line = self.fr.readline().split()
            at = system.atoms_list[-1]
            
            at.idx    = int(line[0])
            at.type   = int(line[1])
            at.name   = line[2]
            at.q      = float(line[3])
            at.pos[0] = float(line[4])
            at.pos[1] = float(line[5])
            at.pos[2] = float(line[6])
        
        return system
    
            
            