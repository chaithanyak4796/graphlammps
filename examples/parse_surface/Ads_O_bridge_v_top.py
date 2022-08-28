#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 12:22:18 2022

@author: chaithanya
"""
import numpy as np
import matplotlib.pyplot as plt
from graphlammps import bond
from graphlammps.special import parse_bonds

plt.close("all")

bond_fname = "bonds.out"
dt = 0.1e-3

bonds = bond.bonds(bond_fname, dt)
parse_bonds = parse_bonds.parse_bonds(bonds)

num_steps = len(bonds.line_offset[:, 0])

time = []
num_bridge = []
num_top = []


for i in range(0, num_steps, 100):
    ts = bonds.line_offset[i][0]

    time.append(ts * dt)
    num_bridge.append(-1)
    num_top.append(-1)

    mols_count = parse_bonds.process_timestep(ts)

    nb = nt = 0

    if "Bridge site O atom" in mols_count:
        nb = len(mols_count["Bridge site O atom"])
    if "Top site O atom" in mols_count:
        nt = len(mols_count["Top site O atom"])

    if nb + nt > 0:
        num_bridge[-1] = nb / (nb + nt)
        num_top[-1] = nt / (nb + nt)


plt.figure(1)
plt.plot(time, num_bridge, "o", label="Bridge")
plt.plot(time, num_top, "o", label="Top")
plt.ylim([0, 1])
plt.legend()
plt.xlabel("Time [ps]")
plt.ylabel("Fraction of different sites")
