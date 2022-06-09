#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 15:51:56 2022

@author: ckondur
"""

import numpy as np
import sys
import subprocess

class bond_info:
    def __init__(self):
        self.id    = 1
        self.type  = 1
        self.nb    = 1
        self.id_nb = []
        self.mol   = 0
        self.bo_nb = []
        self.abo   = 0
        self.nlp   = 0
        self.q     = 0

class bonds:
    def __init__(self, fname, dt = 0.1E-3):
        self.fname = fname
        self.dt    = dt
        
        try:
            self.fr = open(fname, "r")
        except:
            sys.exit(f"Unable to open file {self.fname}. \nExiting.")
        
        # Create a mapping between the timestep and the line numbers
        cmd = ["grep", "-n", "Timestep", self.fname]
        res = subprocess.run(cmd, stdout = subprocess.PIPE).stdout.decode()
        res = res.split('\n')[:-1]
        # print(res)
        self.num_timesteps = len(res)
        
        self.line_map = np.zeros((self.num_timesteps,2), dtype=int)
        
        for i in range(self.num_timesteps):
            line     = res[i].split()
            line_no  = line[0].split(':')
            line_no  = int(line_no[0])
            timestep = int(line[-1])
            
            self.line_map[i][0], self.line_map[i][1] = timestep, line_no
        
        # Create the position of file pointers necessary to jump around
        self.line_offset = []
        offset = 0
        
        for line in self.fr:
            self.line_offset.append(offset)
            offset += len(line)
        self.fr.seek(0)
        
    def read_bonds_timestep(self, step):
        # Check if the index is valid
        idx = np.where(self.line_map[:,0] == step)[0]
        
        if (len(idx) == 0):
            sys.exit("Time step not found in bonds file.")
        else:
            idx = idx[0]
           
        # Get the file pointer to move to the appropriate location
        line_no = self.line_map[idx][1] - 1
        self.fr.seek(self.line_offset[line_no])
        
        # Read the bonds information 
        self.read_bonds_info()
        
    def read_bonds_info(self):
        """ This function assumes the file pointer is at the begining of a timestep and reads the bonds information """
        line = self.fr.readline()
        # print(line)
        if (len(line) == 0):
            self.fr.close()
            raise Exception(f" Reached end of file : {self.fname}")
        else:
            self.timestep = int(line.split()[-1])
            self.time     = self.timestep * self.dt
            
            line = self.fr.readline()
            line = self.fr.readline().split()
            
            self.num_atoms = int(line[-1])
            
            for _ in range(4): 
                line = self.fr.readline()
            
            self.bonds_list = []
            
            for i in range(self.num_atoms):
                bond = bond_info()
                line = self.fr.readline().split()
                bond.id    = int(line[0])
                bond.type  = int(line[1])
                bond.nb    = int(line[2])
                bond.id_nb = np.zeros((bond.nb,), dtype=int)
                # print(bond.id_nb)
                for j in range(bond.nb):
                    bond.id_nb[j] = int(line[j+3])
                bond.mol   = int(line[3+bond.nb])
                bond.bo_nb = np.zeros((bond.nb,), dtype=float)
                for j in range(bond.nb):
                    bond.bo_nb[j] = float(line[j+4+bond.nb])
                    
                bond.abo   = float(line[-3])
                bond.nlp   = float(line[-2])
                bond.q      = float(line[-1])
                
                self.bonds_list.append(bond)
                
            line = self.fr.readline()
            
    def read_next_timestep(self):
        self.read_bonds_info()
        
    def close(self):
        self.fr.close()
        
    def get_O2_atoms_id(self, O_type=3):
        num_atoms = len(self.bonds_list)
        O_atoms   = []
        O2_idx    = []
        for i in range(num_atoms):
            b = self.bonds_list[i]
            if(b.type == O_type and b.nb == 1):
                O_atoms.append(i)
        
        for i in range(len(O_atoms)):
            Oa = self.bonds_list[int(O_atoms[i])]
            for j in range(i+1, len(O_atoms)):
                Ob = self.bonds_list[int(O_atoms[j])]
                if (Oa.id in Ob.id_nb):
                    temp = [Oa.id, Ob.id]
                    O2_idx.append([min(temp), max(temp)])
        
        # print("%d O2 molecules found at timestep %d"%(len(O2_idx), self.timestep))
        return O2_idx