#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 12:03:12 2023

@author: chaithanya
"""
import numpy as np
from graphlammps.special import tecplot_data

Tgas_list = [1000, 1200, 1500, 1800, 2000]
# Tgas_list = [1000]
Temp_surf = 1200

tecplot_fname = 'Prob_%d-test.dat'%(Temp_surf)

mech_list = ['ER', 'ER-long', 'Oxidation', 'Chemisorbed']

# 2D version with individual data formatting
var_list = [['T<sub>gas</sub> [K]', '%6.2f'],
            ['<greek>g</greek><sub>ER</sub>', '%8.4E'],
            ['<greek>s</greek><sub>ER</sub>', '%8.4E'],
            ['<greek>g</greek><sub>ER-long</sub>', '%8.4E'],
            ['<greek>s</greek><sub>ER-long</sub>', '%8.4E'],
            ['<greek>g</greek><sub>oxy</sub>', '%8.4E'],
            ['<greek>s</greek><sub>oxy</sub>', '%8.4E'],
            ['<greek>g</greek><sub>ads</sub>', '%8.4E'],
            ['<greek>s</greek><sub>ads</sub>', '%8.4E']]

# 1D version with default data formatting
# var_list = ['T<sub>gas</sub> [K]',
#             '<greek>g</greek><sub>ER</sub>',
#             '<greek>s</greek><sub>ER</sub>',
#             '<greek>g</greek><sub>ER-long</sub>',
#             '<greek>s</greek><sub>ER-long</sub>',
#             '<greek>g</greek><sub>oxy</sub>',
#             '<greek>s</greek><sub>oxy</sub>',
#             '<greek>g</greek><sub>ads</sub>',
#             '<greek>s</greek><sub>ads</sub>']

tec_data = np.zeros((len(Tgas_list), len(var_list)), dtype=float)

for i,Temp_gas in enumerate(Tgas_list):
    prob_fname = './Prob_%d_%d.npy'%(Temp_surf, Temp_gas)
    Prob = np.load(prob_fname, allow_pickle=True).item()
    
    tec_data[i][0] = Temp_gas
    tec_data[i][1] = Prob['ER'][0]
    tec_data[i][2] = Prob['ER'][1]
    tec_data[i][3] = Prob['ER-long'][0]
    tec_data[i][4] = Prob['ER-long'][1]
    tec_data[i][5] = Prob['Oxidation'][0]
    tec_data[i][6] = Prob['Oxidation'][1]
    tec_data[i][7] = Prob['Chemisorbed'][0]
    tec_data[i][8] = Prob['Chemisorbed'][1]
    
mytec = tecplot_data.tecplot_data()
mytec.update_data(tec_data, var_list)
mytec.write_data_file(tecplot_fname)
