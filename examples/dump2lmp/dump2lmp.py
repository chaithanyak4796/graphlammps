import numpy as np
from graphlammps import read_structure, create

read_structure.use_cache = False

fname_pref = 'coords'
fname_dump = fname_pref + '.dump'

# Reading the dump
read_dump = read_structure.read_structure(fname_dump, read_vel=True)
system    = read_dump.read_dump_timestep(read_dump.line_offset[0][0])

# Creating an instance of the create class
cr = create.create()
cr.struc_fname_prefix  = fname_pref
cr.ortho_box           = True

# Manipulating the number of atom types in the lmp file. 
system.get_num_atom_types()   # The number of atom types is not updated during read_dump_timestep(). This has to be initialized before we can manipulate.
system.num_atom_types += 5    # If no initialization is done, the number of atom types currently present is evaluated during cr.write_struc() and can not be manipulated. 

cr.write_struc(system, write_xyz=False)
cr.write_velocities(system)
