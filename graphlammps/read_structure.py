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

default_cols = ['id', 'type', 'element', 'q', 'x', 'y', 'z']
#----------------------------------------------------------------------------------------------------------------#
class read_structure:
    """ This class contains the functions to read the various formats of lammps structure outputs.
        Currently supported : .xyz, .dump
    
    """
    def __init__(self, fname, cols=default_cols, dt=0.1E-3):
        self.num_timesteps = 1
        self.fname = fname
        self.cols  = cols
        self.default_cols = default_cols
        self.dt    = dt
        
        if(self.cols != self.default_cols):
            sys.exit("The column mapping provided is currently not supported.")
        
        if(".dump" in self.fname):
            self.read_next_timestep   = self.read_dump_next_timestep
            self.read_struc_timestep  = self.read_dump_timestep
            self.init_read_dump()
        else:
            sys.exit(" Invalid file extension. ")
            
    def init_read_dump(self):
        try:
            self.fr = open(self.fname, "r")
        except:
            print(" Unable to open file : ", self.fname)
            sys.exit(1)
        
        # # Find the total number of time steps in the dump file
        # cmd = ["grep", "TIMESTEP", self.fname]
        # res = subprocess.run(cmd, stdout = subprocess.PIPE)
        # self.num_timesteps = len(res.stdout)//len('ITEM: TIMESTEP\n')
        
        # print("Total number of timesteps in file = ", self.num_timesteps)
        
        
        # # Create a mapping between the timestep and the line numbers
        # cmd = ["grep", "-n", "Timestep", self.fname]
        # res = subprocess.run(cmd, stdout = subprocess.PIPE).stdout.decode()
        # res = res.split('\n')[:-1]
        # # print(res)
        # self.num_timesteps = len(res)
        
        # self.line_map = np.zeros((self.num_timesteps,2), dtype=int)
        
        # for i in range(self.num_timesteps):
        #     line     = res[i].split()
        #     line_no  = line[0].split(':')
        #     line_no  = int(line_no[0])
        #     timestep = int(line[-1])
            
        #     self.line_map[i][0], self.line_map[i][1] = timestep, line_no
        
        # Create the position of file pointers necessary to jump around
        self.line_offset = []
        offset = 0
        
        with open(self.fname) as fr:
            for line in fr:
                offset_beg = offset
                offset += len(line)
        
                if("TIMESTEP" in line.split()):
                    line = fr.readline()
                    offset += len(line)
                    ts = int(line)
        
                    self.line_offset.append([ts, offset_beg])
                    
                    line = fr.readline()
                    offset += len(line)
                    
                    line = fr.readline()
                    offset += len(line)
                    natoms = int(line)
                    
                    for i in range(natoms+5):
                        line = fr.readline()
                        offset += len(line)

        self.fr.seek(0)
        self.line_offset   = np.array(self.line_offset)
        self.num_timesteps = len(self.line_offset)
    
    def read_dump_next_timestep(self):
        system = self.read_dump_info()
        return system
        
    def read_dump_timestep(self, step):
        # Check if the step is valid
        idx = np.where(self.line_offset[:,0] == step)[0]
        
        if (len(idx) == 0):
            sys.exit("Time step not found in bonds file.")
        else:
            idx = idx[0]
            
        # Get the file pointer to move to the appropriate location
        self.fr.seek(self.line_offset[idx][1])
        
        # Read the dump information 
        system = self.read_dump_info()
        return system
        
    def read_dump_info(self):
        """ This function reads the next time step in the dump file and returns a lmp_system object """
        
        system = lmp_system.lmp_system()

        line  = self.fr.readline()  # ITEM: TIMESTEP
        if(len(line) == 0):
            self.fr.close()
            sys.exit("Reached the end of the file. Safely exiting.")
            
        line  = self.fr.readline().split()
        system.timestep = int(line[0])
        system.time     = system.timestep * self.dt
        
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
    
    def close(self):
        self.fr.close()
            
#----------------------------------------------------------------------------------------------------------------#

