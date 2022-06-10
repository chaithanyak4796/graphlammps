import numpy as np
from graphlammps import create, lmp_system

cr = create.create()

cr.is_graphene  = False
cr.use_lat_file = False  # Important (related to temp)

cr.lat_a = 1.45
cr.lat_c = 3.45

cr.num_rep      = np.array([10,10,2],dtype=int)
cr.vac_height   = np.array([10.0, 120.0], dtype=float)  # Vaccuum height along the z direction [bottom, top]

cr.mark_corner     = False    # Do you want to mark the corners? Useful to ensure the lattice does not displace during simulations
cr.mark_all_layers = False   # Do you want to mark each graphene layers individually?
cr.mark_top_layer  = False    # Mark top layer alone seperately

cr.struc_fname_prefix = "./graphite_slab"
cr.write_lmp = True
cr.write_xyz = True


system = cr.create_system()

cr.write_struc(system, cr.write_lmp, cr.write_xyz)
