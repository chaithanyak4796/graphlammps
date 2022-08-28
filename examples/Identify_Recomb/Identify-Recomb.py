#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:22:56 2022

@author: chaithanya
"""
import numpy as np
import matplotlib.pyplot as plt
from graphlammps import read_structure, bond, species
from graphlammps.special import sequential_colls
import time
import sys

En = 1.0
dt = 1e-4  # ps
z_GP = 5.0

if len(sys.argv) == 2:
    En = float(sys.argv[1])
# else:
#     z_GP = 5.0

write_reaction_count = True
reaction_count_fname = f"Count.dat"

spec_name_list = ["O_ads", "O", "O2", "CO", "CO2"]
prod_name_list = ["O2", "CO", "CO2"]


def get_mol_idx(mols_list, atom_id):
    iden = -1000  # If mol is not found
    for i in range(len(mols_list)):
        if atom_id[0] == mols_list[i].idx[0] and atom_id[1] == mols_list[i].idx[1]:
            iden = i
            break

    return iden


def identify(En):
    t_start = time.time()
    print("Energy [eV] = ", En)
    Dir = "../" + str(En) + "eV/"

    fname_dump = Dir + "coords.dump"
    fname_spec = Dir + "species.out"
    fname_bond = Dir + "bonds.out"

    # Initialize the read_structure class
    read_dump = read_structure.read_structure(fname_dump, dt=dt)
    max_timesteps = read_dump.num_timesteps - 21
    # max_timesteps = 25000

    # Initialize the bonds class
    bonds = bond.bonds(fname_bond, dt=dt)

    product_list = {}
    for prod in prod_name_list:
        product_list[prod] = []

    O2_global = []
    O2_mols = []
    short_O2_mols = []

    # Identify the indices of the O2 molecules at every nskip timesteps, and create the list of O2 molecule product classes.
    print(
        "\nCoarsely parsing the bonds file to find instances of O2 molecules formation."
    )

    nskip = 5  # Skip this many timesteps in the coarse sweep
    n_check_mol = 20  # Number of further time steps to see
    n_cutoff = 10  # It should be a molecule for at least this many time steps

    for i in range(0, max_timesteps, nskip):
        ts = read_dump.line_offset[i + 1][0]
        if i % 5000 == 0:
            print("Time [ps] = %.2f" % (ts * dt - 0.01))
        bonds.read_bonds_timestep(ts)
        O2_idx = bonds.get_O2_atoms_id()

        product_list["O2"].append([ts, O2_idx])

        for x in O2_idx:
            if x not in O2_global:

                molecule = sequential_colls.product("O2")
                molecule.idx = np.copy(x)
                molecule.ts = ts

                idx_ts = np.where(ts == bonds.line_offset[:, 0])[0][0]

                # Go back nskip timesteps and adjust the timestamp when it was created
                idx_beg = idx_ts - nskip
                j = idx_beg
                while j < idx_ts:
                    ts_curr = bonds.line_offset[j][0]
                    bonds.read_bonds_timestep(ts_curr)
                    temp = bonds.get_O2_atoms_id()
                    if molecule.idx.tolist() in temp:
                        molecule.ts = ts_curr
                        break
                    j += 1

                # In the next n_check_mol timesteps, figure out for how many timesteps, this remains a molecule
                count_mol = 0
                bonds.read_bonds_timestep(ts_curr)
                for j in range(n_check_mol):
                    bonds.read_next_timestep()
                    temp = bonds.get_O2_atoms_id()
                    if molecule.idx.tolist() in temp:
                        count_mol += 1

                molecule.count_mol = count_mol

                # Compare count_mol to cutoff to see if it remains a molecule.
                if count_mol < n_cutoff:
                    is_mol = False
                    molecule.mechanism = "Short"
                    short_O2_mols.append(molecule)
                else:
                    is_mol = True
                    molecule.mechanism = "None"
                    O2_global.append(x)
                    O2_mols.append(molecule)

    #%%

    print("Idenitfying the reaction mechanism of %d O2 molecules\n" % (len(O2_mols)))

    # Go through each molecule and idenitfy the reaction mechanism for each O2 molecule
    n_check_trap = 100  # 1.0 ps
    cutoff_exch = 50  # 0.5 ps
    cutoff_trap = 10  # 0.1 ps
    type_C = 1
    type_O = 3

    for mol in O2_mols:
        ts_curr = mol.ts

        # First check if it is a gas-phase reaction
        system = read_dump.read_dump_timestep(ts_curr)
        Oa = system.atoms_list[
            mol.idx[0] - 1
        ]  # Assuming here that the dump coords is sorted
        Ob = system.atoms_list[
            mol.idx[1] - 1
        ]  # Assuming here that the dump coords is sorted

        if Oa.pos[2] > z_GP and Ob.pos[2] > z_GP:
            mol.mechanism = "GP"
        else:  # It is not a gas-phase reaction. Proceed
            ## Find the amount of tim both the O atoms are bonded to various atoms
            trap_1, trap_2 = 0.0, 0.0
            trap_C = 0
            ex = 0

            ts_curr = mol.ts
            idx_curr = np.where(ts_curr == bonds.line_offset[:, 0])[0][0]
            idx = idx_curr
            while idx > idx_curr - n_check_trap:
                ts = bonds.line_offset[idx][0]
                bonds.read_bonds_timestep(ts)
                b = bonds.bonds_list

                # Amount of time the 1st O atom is trapped
                num_nb1 = b[mol.idx[0] - 1].nb
                if num_nb1 > 0:
                    trap_1 += 1

                # Amount of time the 1st O atom is bonded to C
                for c in b[mol.idx[0] - 1].id_nb:
                    if b[c - 1].type == type_C:
                        trap_C += 1

                # Amount of time the 2nd O atom is trapped
                if len(b) > mol.idx[1] - 1:
                    num_nb2 = b[mol.idx[1] - 1].nb

                    for c in b[mol.idx[1] - 1].id_nb:
                        if b[c - 1].type == type_C:
                            trap_C += 1
                else:
                    num_nb2 = 0

                if num_nb2 > 0:
                    trap_2 += 1

                # Exchange reactions
                c = mol.idx[0] - 1
                if b[c].nb == 1:
                    d = b[c].id_nb[0] - 1
                    if b[d].type == type_O and b[d].id != mol.idx[1]:
                        ex += 1

                c = mol.idx[1] - 1
                if len(b) > mol.idx[1] - 1:
                    if b[c].nb == 1:
                        d = b[c].id_nb[0] - 1
                        if b[d].type == type_O and b[d].id != mol.idx[0]:
                            ex += 1

                idx -= 1

            min_trap = min(trap_1, trap_2)
            max_trap = max(trap_1, trap_2)

            if ex > cutoff_exch:
                mol.mechanism = "Exch"
            elif max_trap < 10:
                mol.mechanism = "Phys"
            else:
                if trap_C < cutoff_trap:
                    mol.mechanism = "Phys"
                elif max_trap == n_check_trap and min_trap == n_check_trap:
                    mol.mechanism = "LH"
                else:
                    mol.mechanism = "ER"

    # Close the necessary files
    bonds.close()
    read_dump.close()

    # Create the reaction_mechanism dictionary
    mech_list = ["ER", "LH", "Phys", "GP", "Exch"]
    reaction_mechanism = {}
    for mech in mech_list:
        reaction_mechanism[mech] = []
    reaction_mechanism["None"] = []
    reaction_mechanism["Short"] = []

    for mol in O2_mols:
        if mol.mechanism in mech_list:
            reaction_mechanism[mol.mechanism].append(mol)
        else:
            reaction_mechanism["None"].append(mol)

    for mol in short_O2_mols:
        reaction_mechanism["Short"].append(mol)

    for mech in reaction_mechanism:
        print(f"Number of {mech} = ", len(reaction_mechanism[mech]))

    # Write mechanism files
    pref = Dir + "Mechanism."
    for mech in reaction_mechanism:
        fname = pref + mech
        fw = open(fname, "w")

        for mol in reaction_mechanism[mech]:
            fw.write(
                "%10d  %4d  %4d  %3d\n"
                % (mol.ts, mol.idx[0], mol.idx[1], mol.count_mol)
            )
        fw.close()

    t_end = time.time()
    print("Time taken [s] = ", t_end - t_start, "\n")

    return reaction_mechanism


#%%
if __name__ == "__main__":
    En_list = [0.1, 1.0, 5.0, 8.0, 10.0]
    # En_list = [5.0, 10.0]
    reaction = {}

    for En in En_list:
        reaction[En] = identify(En)

    if write_reaction_count:
        fw = open(reaction_count_fname, "w")
        fw.write("En[eV] ER    LH  Phys GP   Exch None  Short\n")
        for En in En_list:
            fw.write(
                "%4.2f  %3d  %3d  %3d  %3d  %3d   %3d    %3d\n"
                % (
                    En,
                    len(reaction[En]["ER"]),
                    len(reaction[En]["LH"]),
                    len(reaction[En]["Phys"]),
                    len(reaction[En]["GP"]),
                    len(reaction[En]["Exch"]),
                    len(reaction[En]["None"]),
                    len(reaction[En]["Short"]),
                )
            )
        fw.close()
