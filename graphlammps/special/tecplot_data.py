#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:44:06 2023

@author: chaithanya
"""
import numpy as np
import sys

class tecplot_data:
    """ This is a lightweight class primarily intended to write 2D numpy arrays in the .dat format that Tecplot can read natively """
    def __init__(self):
        self.data   = np.zeros((2,1), dtype=float)
        self.header_var  = ''
        self.header_zone = ''
        self.default_format = '%8.4E'
        self.zone = 'Data'
        
    def update_data(self, data, var_list, zone = "Data"):
        self.data = data
        self.var_names  = []
        self.var_format = []
        self.zone = zone
        
        if isinstance(var_list[0], list):
            for l in var_list:
                self.var_names.append(l[0])
                self.var_format.append(l[1])
        else:
            self.var_names = var_list
            self.var_format = [self.default_format for i in range(len(var_list))]
            
            
        
        if(len(self.data.shape) > 2):
            sys.exit("tecplot_data class in graphlammps currently only supports 1D or 2D data.")
        
        self.num_rows = self.data.shape[0]
        if(len(self.data.shape) == 2):
            self.num_cols = self.data.shape[1]
        else:
            self.num_cols = self.num_rows
            self.num_rows = 1
        
        
        if(self.num_cols != len(self.var_names)):
            sys.exit("The number of cols in data do not match the number of variables in the list")
        
        # Header for Variables names
        self.header_var = 'VARIABLES = '
        for i in range(self.num_cols):
            var = self.var_names[i]
            self.header_var += f'\"{var}\", '
        self.header_var = self.header_var[:-2] + ' \n'
        
        # Header for zone info
        self.header_zone = f'ZONE T=\"{self.zone}\", I={self.num_rows}, F=POINT\n'
        
        
    def write_data_file(self, fname='sample.dat', write_mode="w"):
        self.fname = fname
        if (write_mode != "w" and write_mode != "a"):
            print("Warning: Invalid write mode. Rewriting the file")
            write_mode = "w"
        with open(self.fname, write_mode) as ftec:
            ftec.write(self.header_var)
            ftec.write(self.header_zone)
            
            for i in range(self.num_rows):
                for j in range(self.num_cols):
                    ftec.write(f"{self.var_format[j]} "%(self.data[i][j]))
                ftec.write("\n")
    