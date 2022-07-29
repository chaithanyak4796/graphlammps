#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 11:34:43 2022

@author: ckondur
"""

import numpy as np
import sys
import subprocess
import os

from graphlammps import lmp_system, atom

cols_vel_no  = ['id', 'type', 'element', 'q', 'x', 'y', 'z']
cols_vel_yes = ['id', 'type', 'element', 'q', 'x', 'y', 'z', 'vx', 'vy', 'vz']
use_cache    = True # Use cache or not?
warn_cache   = True # Print the cache warnings or not? 
#----------------------------------------------------------------------------------------------------------------#
class read_structure:
    """ This class contains the functions to read the various formats of lammps structure outputs.
        Currently supported : .xyz, .dump
    
    """
    def __init__(self, fname, dt=0.1E-3, read_vel=False):
        self.num_timesteps = 1
        self.fname = fname
        self.dt    = dt
        
        self.read_vel = read_vel 
        if(self.read_vel):
            self.cols = [cols_vel_yes]
        else:
            self.cols = [cols_vel_no, cols_vel_yes]
            
        self.use_cache  = use_cache     # Use cache or not?
        self.warn_cache = warn_cache    # Print the cache warnings or not?
        
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
        
        # Create the position of file pointers necessary to jump around
        self.line_offset = []
        offset = 0
        
        if(self.use_cache):
            cache_struc_fname = self.fname + ".cache" 
            cache_exist       = os.path.isfile(cache_struc_fname)
            
            if not cache_exist:
                write_cache = True
                if(self.warn_cache): print("Dump cache does not exist. Will write new cache.")
            else:
                write_cache = False
                fc = open(cache_struc_fname, "r+")
                file_fname = fc.readline()
                file_mtime = fc.readline()
                mtime = "%s"%(os.path.getmtime(self.fname))
                if (file_fname.strip() != self.fname.strip() or file_mtime.strip() != mtime.strip() ):
                    write_cache = True
                    if(self.warn_cache): print(" Dump cache exists, but is outdated. Will write new cache.")
                    # print(file_fname.strip(),self.fname.strip())
                    # print(file_mtime.strip(),mtime.strip())
                else:
                    if(self.warn_cache): print("Dump cache exists and is up to date. Will read the cache.")
                    self.num_timesteps = int(fc.readline())
                    for i in range(self.num_timesteps):
                        line = fc.readline().split()
                        self.line_offset.append([int(line[0]), int(line[1])])
                fc.close()
        
        if ((self.use_cache == False) or (self.use_cache and write_cache)):
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
                            
        if self.use_cache and write_cache:
            if(self.warn_cache): print("Writing cache for the dump file.")
            fc = open(cache_struc_fname, "w")
            fc.write("%s\n"%(self.fname))
            mtime = os.path.getmtime(self.fname)
            fc.write("%s\n"%(mtime))
            fc.write("%d\n"%(len(self.line_offset)))
            for i in range(len(self.line_offset)):
                fc.write("%d %d\n"%(self.line_offset[i][0], self.line_offset[i][1]))
            fc.close()

        self.fr.seek(0)
        self.line_offset   = np.array(self.line_offset)
        self.num_timesteps = len(self.line_offset)
    
    def read_dump_next_timestep(self):
        system = self.__read_dump_info()
        return system
        
    def read_dump_timestep(self, step):
        # Check if the step is valid
        idx = np.where(self.line_offset[:,0] == step)[0]
        
        if (len(idx) == 0):
            sys.exit("Time step not found in structure file.")
        else:
            idx = idx[0]
            
        # Get the file pointer to move to the appropriate location
        self.fr.seek(self.line_offset[idx][1])
        
        # Read the dump information 
        system = self.__read_dump_info()
        return system
        
    def __read_dump_info(self):
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
            system.xlo, system.xhi = [float(l) for l in line_1]
            system.ylo, system.yhi = [float(l) for l in line_2]
            system.zlo, system.zhi = [float(l) for l in line_3]
            system.box = np.array([system.xlo, system.xhi, system.ylo, system.yhi, system.zlo, system.zhi])
            
        line = self.fr.readline().split() # ITEM: ATOMS
        line = line[2:]
        
        if (line not in self.cols):
            sys.exit("Coloumns in file do not match the given mapping")
        
        if not self.read_vel:
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
        else:
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
                at.vel[0] = float(line[7])
                at.vel[1] = float(line[8])
                at.vel[2] = float(line[9])
                
        return system
    
    def close(self):
        self.fr.close()
            
#----------------------------------------------------------------------------------------------------------------#

