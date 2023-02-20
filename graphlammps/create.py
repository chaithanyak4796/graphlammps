#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sys
from graphlammps import params, lmp_system, atom
from graphlammps.params import mass_C, mass_O, mass_N, kb, NA

class create:
    """
        This is a class that provides various functions to create the initial structure
    """
    def __init__(self):
        self.is_graphene  = True   # Is the material graphene or graphite?
        self.lat_a        = 1.42   # Lattice constant a (A)
        self.lat_c        = 3.35   # Lattice constant c (A)
        self.num_rep      = np.array([10,10,2],dtype=int)         # Number of unit cells in each direction
        self.vac_height   = np.array([10.0, 120.0], dtype=float)  # Vaccuum height along the z direction [bottom, top]
        self.ortho_box    = True  # Using an orthogonal (True) or hexagonal (False) basis
        
        self.use_lat_file = True   # Get the lattice constants from file
        self.Lat_file     = params.Lat_const_file 
        self.Lat_Temp     = 300    # Lattice Temperature in K
        
        # Marking means assigning the atoms a different atom_type
        self.mark_corner     = True    # Do you want to mark the corners? Useful to ensure the lattice does not displace during simulations
        self.mark_all_layers = False   # Do you want to mark each graphene layers individually?
        self.mark_top_layer  = True    # Mark top layer alone seperately
        
        self.struc_fname_prefix = "./graphite_slab"
        self.write_lmp = True
        self.write_xyz = True
        
        
        self.available_IPs = ['GAP-20', 'COMB3', 'reaxFF', 'reaxFF-CHO', 'reaxFF-CHON', 'LCBOP']
        
    def create_system(self, IP='reaxFF-CHO'):
        """ This function creates the lammps system """
        
        system = lmp_system.lmp_system()
        
        self.IP = IP
        
        # Estimate the number of lattice atoms
        if self.is_graphene:
            self.num_rep[2] = 1
            system.num_atoms    = 2 * np.product(self.num_rep)
        else:
            system.num_atoms    = 4 * np.product(self.num_rep)
        
        # Determine the lattice parameters
        if (self.use_lat_file):
            if not (self.IP in self.available_IPs):
                sys.exit(f" Invlaid IP given.\n Available IPs are  {self.available_IPs}")
            self.lat_a , self.lat_c = self.get_lat_const_from_file()
            
        self.lat_a *= 3**0.5
        self.lat_c *= 2.0
        
        
        # Calculate the box dimensions
        self.create_box(system)
        
        # Populate the lattice
        pos             = np.zeros((3,), dtype=float)
        vel             = np.zeros_like(pos)
        base_lattice    = np.zeros_like(pos)
        
        num_rep     = np.copy(self.num_rep)
        lattice_vec = system.lattice_vec
        
        idx = 1
        
        if not self.ortho_box:   # Non-orthogonal box
            for i in range(num_rep[0]):
                for j in range(num_rep[1]):
                    for k in range(num_rep[2]):
                        base_lattice[0] = i*lattice_vec[0][0] + j*lattice_vec[1][0] + k*lattice_vec[2][0]
                        base_lattice[1] = i*lattice_vec[0][1] + j*lattice_vec[1][1] + k*lattice_vec[2][1]
                        base_lattice[2] = i*lattice_vec[0][2] + j*lattice_vec[1][2] + k*lattice_vec[2][2]
                        
                        pos = np.copy(base_lattice) 
                        
                        atom_type = 2*k + 1
                        system.atoms_list.append(atom.atom(mass_C, idx, atom_type,'C'))
                        system.atoms_list[-1].pos     = np.copy(pos)
                        system.atoms_list[-1].vel     = np.copy(vel)
                        system.atoms_list[-1].layer   = 2*k
                        system.atoms_list[-1].lattice = True
                                 
                        pos[0] = base_lattice[0] + 0.5*self.lat_a
                        pos[1] = base_lattice[1] + (0.5/3**0.5) * self.lat_a
                        pos[2] = base_lattice[2] + 0
                        
     
                        system.atoms_list.append(atom.atom(mass_C, idx+1, atom_type,'C'))
                        system.atoms_list[-1].pos     = np.copy(pos)
                        system.atoms_list[-1].vel     = np.copy(vel)
                        system.atoms_list[-1].layer   = 2*k 
                        system.atoms_list[-1].lattice = True
                        
                        if self.mark_corner:
                            if  (i == 0             and j == 0           ): system.atoms_list[-2].corner = True
                            elif(i == 0             and j == num_rep[1]-1): system.atoms_list[-1].corner = True
                            elif(i == num_rep[0] -1 and j == 0           ): system.atoms_list[-2].corner = True
                            elif(i == num_rep[0] -1 and j == num_rep[1]-1): system.atoms_list[-1].corner = True
                            
                        if not self.is_graphene:
                            pos[0] = base_lattice[0] + 0.5*self.lat_a
                            pos[1] = base_lattice[1] + (0.5/3.0**0.5) * self.lat_a
                            pos[2] = base_lattice[2] + 0.5*self.lat_c
                            
                            atom_type = 2*k + 2
                            system.atoms_list.append(atom.atom(mass_C, idx+2, atom_type,'C'))
                            system.atoms_list[-1].pos     = np.copy(pos)
                            system.atoms_list[-1].vel     = np.copy(vel)
                            system.atoms_list[-1].layer   = 2*k + 1
                            system.atoms_list[-1].lattice = True
                        
                            pos[0] = base_lattice[0] + 1.0*self.lat_a
                            pos[1] = base_lattice[1] + (1.0/3.0**0.5) * self.lat_a
                            pos[2] = base_lattice[2] + 0.5*self.lat_c
    
                            system.atoms_list.append(atom.atom(mass_C, idx+3,atom_type,'C'))
                            system.atoms_list[-1].pos     = np.copy(pos)
                            system.atoms_list[-1].vel     = np.copy(vel)
                            system.atoms_list[-1].layer   = 2*k + 1
                            system.atoms_list[-1].lattice = True
                            
                            if self.mark_corner:
                                if  (i == 0             and j == 0           ): system.atoms_list[-2].corner = True
                                elif(i == 0             and j == num_rep[1]-1): system.atoms_list[-1].corner = True
                                elif(i == num_rep[0] -1 and j == 0           ): system.atoms_list[-2].corner = True
                                elif(i == num_rep[0] -1 and j == num_rep[1]-1): system.atoms_list[-1].corner = True
            
                            idx += 4
                        else:
                            idx += 2
                            
                        z_max = np.copy(pos[2])
        
        else:   # Orthogonal box
            for i in range(num_rep[0]):
                for j in range(num_rep[1]):
                    for k in range(num_rep[2]):
                        base_lattice[0] = i*lattice_vec[0][0] + j*lattice_vec[1][0] + k*lattice_vec[2][0]
                        base_lattice[1] = i*lattice_vec[0][1] + j*lattice_vec[1][1] + k*lattice_vec[2][1]
                        base_lattice[2] = i*lattice_vec[0][2] + j*lattice_vec[1][2] + k*lattice_vec[2][2]
                        
                        height_mod = k%2 
                        
                        pos[0] = base_lattice[0] + (3/4 + height_mod) * self.lat_a/3**0.5
                        pos[1] = base_lattice[1] + (1/4*3**0.5) * self.lat_a/3**0.5
                        pos[2] = base_lattice[2] + self.lat_c/4
                        
                        atom_type = k + 1
                        system.atoms_list.append(atom.atom(mass_C, idx, atom_type,'C'))
                        system.atoms_list[-1].pos     = np.copy(pos)
                        system.atoms_list[-1].vel     = np.copy(vel)
                        system.atoms_list[-1].layer   = k
                        system.atoms_list[-1].lattice = True
                                 
                        pos[0] = base_lattice[0] + (5/4 + height_mod) * self.lat_a/3**0.5
                        pos[1] = base_lattice[1] + (3/4*3**0.5) * self.lat_a/3**0.5
                        pos[2] = base_lattice[2] + self.lat_c/4
                        
                        system.atoms_list.append(atom.atom(mass_C, idx+1, atom_type,'C'))
                        system.atoms_list[-1].pos     = np.copy(pos)
                        system.atoms_list[-1].vel     = np.copy(vel)
                        system.atoms_list[-1].layer   = k 
                        system.atoms_list[-1].lattice = True
                        
                        pos[0] = base_lattice[0] + (9/4 - 2*height_mod) * self.lat_a/3**0.5
                        pos[1] = base_lattice[1] + (3/4*3**0.5) * self.lat_a/3**0.5
                        pos[2] = base_lattice[2] + self.lat_c/4
                        
                        system.atoms_list.append(atom.atom(mass_C, idx+2, atom_type,'C'))
                        system.atoms_list[-1].pos     = np.copy(pos)
                        system.atoms_list[-1].vel     = np.copy(vel)
                        system.atoms_list[-1].layer   = k 
                        system.atoms_list[-1].lattice = True
                        
                        pos[0] = base_lattice[0] + (11/4 - 2*height_mod) * self.lat_a/3**0.5
                        pos[1] = base_lattice[1] + (1/4*3**0.5) * self.lat_a/3**0.5
                        pos[2] = base_lattice[2] + self.lat_c/4
                        
                        system.atoms_list.append(atom.atom(mass_C, idx+3, atom_type,'C'))
                        system.atoms_list[-1].pos     = np.copy(pos)
                        system.atoms_list[-1].vel     = np.copy(vel)
                        system.atoms_list[-1].layer   = k 
                        system.atoms_list[-1].lattice = True
                        
                        if self.mark_corner:
                            if(height_mod == 0):
                                if  (i == 0             and j == 0           ): system.atoms_list[-4].corner = True
                                elif(i == 0             and j == num_rep[1]-1): system.atoms_list[-3].corner = True
                                elif(i == num_rep[0] -1 and j == 0           ): system.atoms_list[-1].corner = True
                                elif(i == num_rep[0] -1 and j == num_rep[1]-1): system.atoms_list[-2].corner = True
                            else:
                                if  (i == 0             and j == 0           ): system.atoms_list[-1].corner = True
                                elif(i == 0             and j == num_rep[1]-1): system.atoms_list[-2].corner = True
                                elif(i == num_rep[0] -1 and j == 0           ): system.atoms_list[-4].corner = True
                                elif(i == num_rep[0] -1 and j == num_rep[1]-1): system.atoms_list[-3].corner = True
                            
                        idx += 4  
                        z_max = np.copy(pos[2])
                        
        num_atoms   = len(system.atoms_list)
        
        # Adjust the z coordinates so that the top layer has z = 0
        system.box[4] -= z_max
        system.box[5] -= z_max
        for i in range(num_atoms):
            system.atoms_list[i].pos[2] -= z_max
        
        ########### Various types of marking the atoms ####################
        # Adjust the marking of the layers
        if (self.mark_all_layers == False or self.mark_top_layer == True) :
            for i in range(num_atoms):
                system.atoms_list[i].type = 1
                if (self.mark_top_layer and system.atoms_list[i].pos[2] == 0): 
                    system.atoms_list[i].type = 2
            
            if self.mark_top_layer: 
                atom_type = 2
            else:
                atom_type = 1
                
        # Special case of marking for graphene
        if self.is_graphene:
            for i in range(num_atoms):
                system.atoms_list[i].type = 1
            atom_type = 1
        
        # Increase the atom_type of the corner atoms
        if self.mark_corner:
            atom_type += 1
            for at in system.atoms_list:
                if at.corner: 
                    at.type = atom_type
                    at.name = 'N'
            
        ##################################################################
        system.num_atom_types = atom_type
        
        # Write the structure files
        # self.write_struc(system)
    
        return system
    
    def generate_velocities_lattice(self, system):
        """ This function generates the velocities of each layer. They are randomly drawn from a Gaussian at lattice temperature."""
        if self.is_graphene:
            num_layers = 1 * self.num_rep[2]
        else:
            num_layers = 2 * self.num_rep[2]
            
        if self.ortho_box:
            num_layers = 1 * self.num_rep[2]
            
        layer_count = np.zeros((num_layers,), dtype=int)

        for at in system.atoms_list:
            if(at.lattice):
                layer_count[at.layer] += 1
        
        if self.mark_corner:
            layer_count -= 4
        
        # Draw the velocities from Gaussian dist
        m   = mass_C/1000/NA
        sig = (kb*self.Lat_Temp/m)**0.5
        
        self.vels = []
        for i in range(num_layers):
            v = np.zeros((layer_count[i], 3), dtype=float)
            v = np.random.normal(0, sig, v.shape).reshape(v.shape)
            
            # Zeroing out the com velocities of each layer
            vcom = 0
            vcom = np.sum(v, axis=0)/layer_count[i]
            v   -= vcom
            
            self.vels.append(v)
        
        if('reaxFF' in self.IP): self.unit_corr = 1E-5   # real  units
        else                   : self.unit_corr = 1E-2   # metal units
        
        # Scaling to ensure correct temperature
        vel2 = 0
        for i in range(num_layers):
            vel2 +=  sum(self.vels[i]**2)
        
        T         = m*sum(vel2)/(3*np.sum(layer_count)*kb)
        for i in range(num_layers):
            self.vels[i] *= (self.Lat_Temp/T)**0.5 * self.unit_corr 
        
        #### Assign the velocities to the lattice #######
        for i in range(len(system.atoms_list)):
            at = system.atoms_list[i]
            if at.lattice:
                if ( (self.mark_corner and at.corner==False) or (self.mark_corner == False) ):
                    layer_count[at.layer] -= 1
                    at.vel = np.copy(self.vels[at.layer][layer_count[at.layer]])
                else:
                    at.vel *= 0.0
                    
        ### First compute the com vel and lattice temperatures #######
        com_vel, tot_mass = self.compute_com_lattice(system)  
        print(" before com_vel = ", com_vel) 
        vel2     = 0
        if self.mark_corner:
            for at in system.atoms_list:
                if (at.lattice and at.corner == False):
                    at.vel -= com_vel
                    vel2 += sum(at.vel**2)
        else:
            for at in system.atoms_list:
                if (at.lattice):
                    at.vel -= com_vel
                    vel2 += sum(at.vel**2)
                    
        com_vel, tot_mass = self.compute_com_lattice(system)             
        T_calc = m * vel2/self.unit_corr**2 /(3*tot_mass/mass_C*kb)
        print(" after com_vel = ", com_vel)   
        print(" Calculated Temp [K] = ", T_calc )
        
    def compute_com_lattice(self, system):
        com_vel = np.zeros((3,), dtype = float)
        tot_mass = 0
        
        if self.mark_corner:
            for at in system.atoms_list:
                if (at.lattice and at.corner == False):
                    com_vel  += at.vel * at.mass
                    tot_mass += at.mass
        else:
            for at in system.atoms_list:
                if (at.lattice):
                    com_vel  += at.vel * at.mass
                    tot_mass += at.mass
        
        com_vel /= tot_mass
        
        return com_vel, tot_mass
                
    def write_struc(self, system, write_lmp=True, write_xyz=True):
        """ This function writes the positions of the atoms in the .lmp and .xyz formats"""
        num_atoms = len(system.atoms_list)
        if(system.num_atom_types <= 0): system.get_num_atom_types()

        if write_lmp:
            self.lmp_fname = self.struc_fname_prefix + ".lmp"
            flmp  = open(self.lmp_fname, "w")
            
            flmp.write(" Hexagonal lattice structure of dimension %d x %d x %d with a = %.4f and c = %.4f\n"%(self.num_rep[0], self.num_rep[1], self.num_rep[2], self.lat_a, self.lat_c))
            flmp.write("%d atoms\n"%(num_atoms))
            flmp.write("%d atom types\n\n"%(system.num_atom_types))
            flmp.write("%6.8E %6.8E xlo xhi\n"%(system.box[0], system.box[1]))
            flmp.write("%6.8E %6.8E ylo yhi\n"%(system.box[2], system.box[3]))
            flmp.write("%6.8E %6.8E zlo zhi\n\n"%(system.box[4], system.box[5]))
            if not self.ortho_box:
                flmp.write("%6.8E  %6.8E  %6.8E  xy xz yz\n\n"%(system.xyz[0], system.xyz[1], system.xyz[2]))
            flmp.write("Atoms\n\n")
            for i in range(num_atoms):
                at = system.atoms_list[i]
                flmp.write("%5d  %2d  %6.4E  %6.8E  %6.8E  %6.8E\n"%(at.idx, at.type, at.q, at.pos[0], at.pos[1], at.pos[2]))
            
            flmp.close()
        
        if write_xyz:
            self.xyz_fname = self.struc_fname_prefix + ".xyz"
            fxyz  = open(self.xyz_fname,"w")
            
            fxyz.write("%d \n"%(num_atoms))
            fxyz.write("Atoms\n")
            for i in range(num_atoms):
                at = system.atoms_list[i]
                fxyz.write("%s %6.8f %6.8f %6.8f\n"%(at.name, at.pos[0], at.pos[1], at.pos[2]))
            fxyz.close()

    def write_velocities(self, system):
        """ This function appends the velocities of all the atoms to the lmp file """
        
        try:
            flmp = open(self.lmp_fname, "a")
        except:
            sys.exit("No lmp file exists to write velocities. Write the atom positions first. ")
            
        flmp.write("\nVelocities\n\n")
        for i in range(len(system.atoms_list)):
            at = system.atoms_list[i]
            flmp.write("%5d  %6.8E  %6.8E  %6.8E\n"%(at.idx, at.vel[0], at.vel[1], at.vel[2]))

        flmp.close()
 
    def create_box(self, system):
        """ This function creates the LAMMPS simulation box """
        self.vac_height = abs(self.vac_height)
        
        if self.ortho_box:
            system.lattice_vec = np.zeros((3,3))
            system.lattice_vec[0][0] = 3**0.5 * self.lat_a
            system.lattice_vec[1][1] = self.lat_a
            system.lattice_vec[2][2] = 0.5 * self.lat_c
            
            system.box = np.zeros((6,))
            system.box[1] = 3**0.5 * self.lat_a * self.num_rep[0]
            system.box[3] = self.lat_a * self.num_rep[1]
            system.box[4] = -1 * self.vac_height[0]
            system.box[5] =  1 * self.vac_height[1] + 0.5 * self.lat_c * self.num_rep[2]
            
            system.lattice_pos = np.zeros((4,3))
            system.lattice_pos[0] = np.array([3/4/3**0.5  * self.lat_a , 0.25 * self.lat_a,  0.25 * self.lat_c])
            system.lattice_pos[1] = np.array([5/4/3**0.5  * self.lat_a , 0.75 * self.lat_a,  0.25 * self.lat_c])
            system.lattice_pos[2] = np.array([9/4/3**0.5  * self.lat_a , 0.75 * self.lat_a,  0.25 * self.lat_c])
            system.lattice_pos[3] = np.array([11/4/3**0.5 * self.lat_a , 0.25 * self.lat_a,  0.25 * self.lat_c])
            
        else:
            system.lattice_vec       = np.zeros((3,3))
            system.lattice_vec[0][0] = self.lat_a
            system.lattice_vec[1][0] = 0.5 * self.lat_a
            system.lattice_vec[1][1] = 0.5 * 3**0.5 * self.lat_a
            system.lattice_vec[2][2] = self.lat_c
            
            system.xyz    = np.zeros((3,))
            system.xyz[0] = 0.5 * self.lat_a * self.num_rep[1] - 0.00001
            
            system.box    = np.zeros((6,))
            system.box[1] = self.lat_a * self.num_rep[0]
            system.box[3] = 0.5 * 3**0.5 * self.lat_a * self.num_rep[1]
            system.box[4] = -1 * self.vac_height[0]
            system.box[5] =  1 * self.vac_height[1] + self.lat_c * self.num_rep[2] 
            
            system.lattice_pos    = np.zeros((4,3))
            system.lattice_pos[1] = np.array([0.5*self.lat_a, (0.5/3**0.5)*self.lat_a, 0])
            system.lattice_pos[2] = np.array([0.5*self.lat_a, (0.5/3**0.5)*self.lat_a, 0.5*self.lat_c])
            system.lattice_pos[3] = np.array([1.0*self.lat_a, (1.0/3**0.5)*self.lat_a, 0.5*self.lat_c])
            
            
    def get_lat_const_from_file(self):
        """ This function reads the lattice constant from the file """
        # Reading the Lattice file 
        fname = self.Lat_file
        print("Extracting the lattice params from file : ", fname)
        pot_name = []
        Lat_const = []
        num_pot = 0
        fr = open(fname,'r')
        line = fr.readline().strip("\n")
        while(line != 'End'):
            if(line [0] == '*'):
                pot = line.strip(" *\n")
                # print (pot)
                pot_name.append(pot)
                num_pot += 1
                Lat_const.append([''])
                line = fr.readline().strip("\n") 
                while(True):
                    line = fr.readline().strip("\n") 
                    if(len(line) > 0): 
                        # print(line)
                        temp = [float(j) for j in line.split()]
                        Lat_const[num_pot-1].append(temp)
                    else:
                        break
            
            line = fr.readline().strip("\n")
    
        pot_name = np.array(pot_name)
        for i in range(num_pot) :
            Lat_const[i].pop(0)
            Lat_const[i] = np.array(Lat_const[i])
        
        # Extracting the lattice params
        idx_pot = np.where(pot_name == self.IP)[0]
        if(len(idx_pot) == 0):
            print(" Error : Potential not found.")
            sys.exit()
        idx_pot = idx_pot[0]
        idx_Temp = np.where(Lat_const[idx_pot][:,0] == self.Lat_Temp)[0]
        if(len(idx_Temp) > 0):
            lat_param = Lat_const[idx_pot][idx_Temp][0]
            print("Lattice params found = ",lat_param)
        else:
            
            print("\n Warning : Temp not found in the potential.")
            
            lat_param    = np.copy(Lat_const[idx_pot][0])
            lat_param[0] = self.Lat_Temp
            lat_param[1] = np.interp(self.Lat_Temp, Lat_const[idx_pot][:,0], Lat_const[idx_pot][:,1])
            lat_param[2] = np.interp(self.Lat_Temp, Lat_const[idx_pot][:,0], Lat_const[idx_pot][:,2])
            print("Lattice params interpolated = ",lat_param)
            
        return lat_param[1], lat_param[2]

