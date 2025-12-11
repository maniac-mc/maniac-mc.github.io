# LAMMPS parameter file 

mass 1 1.008    # HO (CH4O)
mass 2 1.008    # HC (CH4O)
mass 3 12.011  # CH (CH4O)
mass 4 15.9994 # OH (CH4O)

# pair_coeff
pair_coeff 1 1 0.0000 1.000     # HO-HO (CH4O)
pair_coeff 2 2 0.0300 2.500     # HC-HC (CH4O)
pair_coeff 3 3 0.066 3.500    # CH-CH (CH4O)
pair_coeff 4 4 0.170 3.120    # OH-OH (CH4O)

# bond_coeff
bond_coeff 1 0.0 1.09             # HC-CH (CH4O)
bond_coeff 2 0.0 1.41             # CH-OH (CH4O)
bond_coeff 3 0.0 0.945           # HO-OH (CH4O)

# angle_coeff
angle_coeff 1 0.000 108.6         # CH-OH-HO (CH4O)
angle_coeff 2 0.000 109.5         # HC-CH-OH (CH4O)
angle_coeff 3 0.000 107.8         # HC-CH-HC (CH4O)

# dihedral_coeff
dihedral_coeff 1 0.0 0.0 0.0 0.0  # HC-CH-OH-HO (CH4O)
