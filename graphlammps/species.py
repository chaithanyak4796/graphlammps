#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 16:08:13 2022

@author: ckondur
"""

import numpy as np
import sys

name_list_default = ['O_ads', 'O', 'O2', 'CO', 'CO2']

class species:
    
    def __init__(self, fname, dt=0.1E-3):
        self.fname = fname
        self.dt    = dt
        
        try:
            self.fr = open(fname, "r")
        except:
            sys.exit(f"Unable to open file {self.fname}. \nExiting.")
    
        self.initialize_species_count()   
        
        self.line_offset = []
        offset = 0
        
        with open(self.fname) as fr:
            for line in fr:
                offset_beg = offset
                if("Timestep" in line.split()):
                    offset += len(line)
                    
                    line    = fr.readline()
                    offset += len(line)
                    
                    ts = int(line.split()[0])
                    self.line_offset.append([ts, offset_beg])
                else:
                    sys.exit("Error while mapping the species file. Timestep not found at the expected location in the species file !!!!!")
        
        self.fr.seek(0)
        self.line_offset   = np.array(self.line_offset)
        self.num_timesteps = len(self.line_offset)
        # self.spec_dict     = {}
        
    def initialize_species_count(self, name_list = name_list_default):
        """ name_list is the list of species names to be considered """
        self.time = 0.0
        self.num_species = len(name_list)

        self.spec_count = {}
        for i in name_list:
            self.spec_count[i] = 0
    
        self.adsorb = True    # Do you want to track adsorbed oxygen as a distinct species as oppposed to graphene-oxide
        
    
    def read_next_timestep(self):
        self.__read_species_info()
    
    def read_species_timestep(self, step):
        
        # Check if the step is valid
        idx = np.where(self.line_offset[:,0] == step)[0]
        
        if (len(idx) == 0):
            sys.exit("Time step not found in bonds file.")
        else:
            idx = idx[0]
            
        # Get the file pointer to move to the appropriate location
        self.fr.seek(self.line_offset[idx][1])
        
        # Read the bonds information 
        self.__read_species_info()
    
    def __read_species_info(self):
        line = self.fr.readline()
        if (len(line) == 0):
            self.fr.close()
            raise Exception(f" Reached end of file : {self.fname}")
            # sys.exit(f"Reached end of file {self.fname}. \nExiting.")
        else:
            line_head = line.split()
            line_info = self.fr.readline().split()
            
            self.timestep = int(line_info[0])
            self.time     = self.timestep * self.dt
            
            self.num_mols = int(line_info[1])
            self.num_spec = int(line_info[2])
            
            self.spec_info = []
            
            for i in range(self.num_spec):
                spec_name  = line_head[i+4]
                spec_count = int(line_info[i+3])
                self.spec_info.append([spec_name, spec_count])
    
    def close(self):
        self.fr.close()
            
    def parse_species(self):
        """ This function parses the line read from the species file """
        for i in range(len(self.spec_info)):
            name  = self.spec_info[i][0]
            count = self.spec_info[i][1]
            
            # print(name, count)
            
            if 'H' in name:   # lattice
                if 'O' in name:
                    c = name.split('O')
                    # print (c)
                    if (c[-1] == ''): count = 1
                    else            : count = int(c[-1])
                    name  = 'O_ads'
            
            if name in self.spec_count:
                self.spec_count[name] = count