#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 09:50:30 2022

@author: chaithanya
"""

import sys
from graphlammps import molecule

class parse_bonds():
    def __init__(self, bonds):
        self.bonds = bonds
        self.O_type = [3]
        self.C_type = [1,2]

    def idenitfy_bond(self, at_id):
        bonds = self.bonds
        
        mol = molecule.molecule(name='None')
        mol.atom_id = [at_id]
            
        b0 = bonds.get_bond(at_id)
        # print(b0.id, b0.id_nb)
        
        if(b0.type in self.C_type):
            sys.exit('C atom is the chosen atom')
        
        bonds.get_neighbor_info(at_id)
        
        if(b0.nb == 2):
            if (b0.type_nb[0] in self.C_type and b0.type_nb[1] in self.C_type):
                mol.name    = 'Bridge site O atom'
                mol.atom_id = [b0.id]
            else:
                # Find the corner O atom
                for i in range(2):
                    b1 = bonds.get_bond(b0.id_nb[i])
                    if(b1.type in self.O_type and b1.nb == 1):
                        b0 = b1
                        break
                if(b1.type not in self.O_type or b1.nb != 1):
                    LH_cand=False
                    for i in range(2):
                        if(b0.type_nb[i] in self.O_type):
                            b1 = bonds.get_bond(b0.id_nb[i])
                       
                            bonds.get_neighbor_info(b1.id)
                            if(b1.nb > 1):
                               mol.name = 'LH candidate'
                               mol.atom_id = [b0.id, b1.id]
                               LH_cand = True
                    if not LH_cand:
                        sys.exit(f'Error: Corner O atom not found\n Timestep = {bonds.timestep}\n Atom_id = {at_id}')

        
        # Now that corner atom is found, identify the molecule
        bonds.get_neighbor_info(b0.id)
        # print(b0.id, b0.id_nb)
        
        if(b0.nb == 0):
            mol.name    = "Gas-phase O atom"
            mol.atom_id = [b0.id]
            
        elif(b0.nb == 1):
            b1 = bonds.get_bond(b0.id_nb[0])
            bonds.get_neighbor_info(b1.id)
    
            # O-O
            if(b1.type in self.O_type):
                if(b1.nb == 1):
                    mol.name    = "Gas-phase O2 molecule"
                    mol.atom_id = [b0.id, b1.id]
                    
                elif(b1.nb == 2):
                    b2 = bonds.get_bond(b1.id_nb[0])
                    if(b0.id == b2.id):
                        b2 = bonds.get_bond(b1.id_nb[1])    
    
                    if(b2.type in self.C_type):
                        mol.name    = "O2 adsorbed"
                        mol.atom_id = [b0.id, b1.id]
                        
                    elif(b2.type in self.O_type and b2.nb == 1):
                        mol.name    = 'Gas-phase O3'
                        mol.atom_id = [b0.id, b1.id, b2.id]
            # O-C
            elif(b1.type in self.C_type):
                if(b1.nb == 1):
                    mol.name    = "Gas-phase CO molecule"
                    mol.atom_id = [b0.id, b1.id] 
                    
                elif(b1.nb == 2):
                    b2 = bonds.get_bond(b1.id_nb[0])
                    if(b0.id == b2.id):
                        b2 = bonds.get_bond(b1.id_nb[1])
                  
                    if(b2.type in self.O_type):
                        mol.name    = "Gas-phase CO2 molecule"
                        mol.atom_id = [b0.id, b1.id, b2.id]
                        
                    else:
                        mol.name    = "CO adsorbed"
                        mol.atom_id = [b1.id, b0.id]
                        
                elif(b1.nb == 3):
                    bonds.get_neighbor_info(b1.id)
                    count = [0, 0]
                    for i in range(3):
                        if(b1.type_nb[i] in self.O_type):
                            b2 = bonds.get_bond(b1.id_nb[i])
                            if(b2.nb == 1):
                                count[0] += 1
                        elif(b1.type_nb[i] in self.C_type):
                            count[1] += 1
                    if(count[0] == 2 and count[1] == 1):
                        mol.name    = "CO2 adsorbed"
                        for i in range(3):
                            b2 = bonds.get_bond(b1.id_nb[i])
                            if(b2.type in self.O_type and b2.id != b0.id): break
                        
                        mol.atom_id = [b0.id, b1.id, b2.id]
                        
                    elif(count[0] == 3):
                        mol.name = 'Gas-phase CO3'
                        mol.atom_id = [b1.id, b1.id_nb[0], b1.id_nb[1], b1.id_nb[2]]
                    else:
                        mol.name    = "Top site O atom"
                        mol.atom_id = [b0.id]
                elif(b1.nb == 4):
                    bonds.get_neighbor_info(b1.id)
                    count_C = 0
                    for i in range(4):
                        if(b1.type_nb[i] in self.C_type):
                            count_C += 1
    
                    mol.name    = "Top site O atom"
                    mol.atom_id = [b0.id]
        return mol
    


    def process_timestep(self, ts):
        bonds = self.bonds   
        mols_count = {}
        
        # ts = bonds.line_offset[0][0]
        bonds.read_bonds_timestep(ts)
        # print(ts)
        
        mols_list  = []
        iden_at_id = []
        
        iden_names = []
        
        num_O_atoms = 0
        for b in bonds.bonds_list:
            if b.type in self.O_type:
                num_O_atoms += 1
                
                if(b.id not in iden_at_id):
                    mol = self.idenitfy_bond(b.id)
                    mol.mol_id = len(mols_list)
                    mols_list.append(mol)
                    mol.atom_id = sorted(mol.atom_id)
                    iden_names.append(mol.name)
                    
                    for i in mol.atom_id:
                        iden_at_id.append(i)
                
        mols_count = {}
        for name in iden_names:
            mols_count[name] = []
            
        for mol in mols_list:
            mols_count[mol.name].append(mol)
        
        # print("")  
        # for name in mols_count:
        #     print(name, " : ", len(mols_count[name]))
        
        return mols_count
    
    def print_count(self, count):
        """ Print the count of the molecules """
        print(" ")
        for mech in count:
            print(mech, len(count[mech]))
        print(" ")
        
    def get_molecule(self, count, idx):
        """ Return the molecule object which contains the atom """
        for mech in count:
            for mol in count[mech]:
                if (idx in mol.atom_id):
                    return mol
        
        mol = molecule.molecule(name='None')
        if(self.bonds.get_bond(idx).type not in self.O_type):   
            print("Given atom is not an O atom.")
        else:
            print(" Atom not found in any molecule")
        return mol