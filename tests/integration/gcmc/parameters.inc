# LAMMPS parameter file 

mass 1 12.011   # C1 (Zif8)
mass 2 12.011   # C2 (Zif8)
mass 3 12.011   # C3 (Zif8)
mass 4 1.008    # H2 (Zif8)
mass 5 1.008    # H3 (Zif8)
mass 6 14.007   # N  (Zif8)
mass 7 65.38    # Zn (Zif8)
mass 8 1.008    # HO (CH4O)
mass 9 1.008    # HC (CH4O)
mass 10 12.011  # CH (CH4O)
mass 11 15.9994 # OH (CH4O)

# pair_coeff
pair_coeff 1 1 0.0860 3.400     # C1-C1 (Zif8)
pair_coeff 2 2 0.0860 3.400     # C2-C2 (Zif8)
pair_coeff 3 3 0.1094 3.400     # C3-C3 (Zif8)
pair_coeff 4 4 0.0150 2.511     # H2-H2 (Zif8)
pair_coeff 5 5 0.0157 2.650     # H3-H3 (Zif8)
pair_coeff 6 6 0.1700 3.250     # N-N  (Zif8)
pair_coeff 7 7 0.0125 1.960     # Zn-Zn (Zif8)
pair_coeff 8 8 0.0000 1.000     # HO-HO (CH4O)
pair_coeff 9 9 0.0300 2.500     # HC-HC (CH4O)
pair_coeff 10 10 0.066 3.500    # CH-CH (CH4O)
pair_coeff 11 11 0.170 3.120    # OH-OH (CH4O)

# bond_coeff
bond_coeff 1 346.543 1.490          # C1-C3 (Zif8)
bond_coeff 2 488.000 1.335          # C1-N (Zif8)
bond_coeff 3 540.249 1.350          # C2-C2 (Zif8)
bond_coeff 4 367.000 1.080          # C2-H2 (Zif8)
bond_coeff 5 440.221 1.370          # C2-N (Zif8)
bond_coeff 6 340.000 1.090          # C3-H3 (Zif8)
bond_coeff 7 78.500 2.011           # N-Zn (Zif8)
bond_coeff 8 340.0 1.09             # HC-CH (CH4O)
bond_coeff 9 320.0 1.41             # CH-OH (CH4O)
bond_coeff 10 553.0 0.945           # HO-OH (CH4O)

# angle_coeff
angle_coeff 1 48.088 109.32         # C1-C3-H3 (Zif8)
angle_coeff 2 71.254 105.27         # C1-N-C2 (Zif8)
angle_coeff 3 48.680 128.33         # C1-N-Zn (Zif8)
angle_coeff 4 49.451 125.67         # C2-C2-H2 (Zif8)
angle_coeff 5 76.750 108.65         # C2-C2-N (Zif8)
angle_coeff 6 32.477 126.40         # C2-N-Zn (Zif8)
angle_coeff 7 66.015 123.92         # C3-C1-N (Zif8)
angle_coeff 8 49.954 125.68         # H2-C2-N (Zif8)
angle_coeff 9 35.000 109.50         # H3-C3-H3 (Zif8)
angle_coeff 10 75.484 112.16        # N-C1-N (Zif8)
angle_coeff 11 35.240 109.48        # N-Zn-N (Zif8)
angle_coeff 12 55.000 108.6         # CH-OH-HO (CH4O)
angle_coeff 13 35.000 109.5         # HC-CH-OH (CH4O)
angle_coeff 14 33.000 107.8         # HC-CH-HC (CH4O)

# dihedral_coeff

dihedral_coeff * charmm 0.000 0 0 0.0
dihedral_coeff * opls 0.0 0.0 0.0 0.0

dihedral_coeff 1 charmm 0.000 0 0 0.0      # C1-N-Zn-N (Zif8)
dihedral_coeff 2 charmm 2.325 2 180 0.0    # C2-C2-N-C1 (Zif8)
dihedral_coeff 3 charmm 2.325 2 180 0.0    # C2-C2-N-Zn (Zif8)
dihedral_coeff 4 charmm 0.000 0 0 0.0      # C2-N-Zn-N (Zif8)
dihedral_coeff 5 charmm 5.000 2 180 0.0    # C3-C1-N-C2 (Zif8)
dihedral_coeff 6 charmm 5.000 2 180 0.0    # C3-C1-N-Zn (Zif8)
dihedral_coeff 7 charmm 5.150 2 180 0.0    # H2-C2-C2-H2 (Zif8)
dihedral_coeff 8 charmm 5.150 2 180 0.0    # H2-C2-C2-N (Zif8)
dihedral_coeff 9 charmm 2.325 2 180 0.0    # H2-C2-N-C1 (Zif8)
dihedral_coeff 10 charmm 2.325 2 180 0.0   # H2-C2-N-Zn (Zif8)
dihedral_coeff 11 charmm 0.000 0 0 0.0     # H3-C3-C1-N (Zif8)
dihedral_coeff 12 charmm 5.000 2 180 0.0   # N-C1-N-C2 (Zif8)
dihedral_coeff 13 charmm 5.000 2 180 0.0   # N-C1-N-Zn (Zif8)
dihedral_coeff 14 charmm 5.150 2 180 0.0   # N-C2-C2-N (Zif8)
dihedral_coeff 15 opls 0.0 0.0 0.45 0.0  # HC-CH-OH-HO (CH4O)

# improper_coeff
improper_coeff 1 0.0 0.0            # C1-C3-H3-H3 (Zif8)
improper_coeff 2 2.0 180.0          # C1-N-C2-Zn (Zif8)
improper_coeff 3 2.0 180.0          # C2-C2-H2-N (Zif8)
improper_coeff 4 2.0 180.0          # C3-C1-N-N (Zif8)
improper_coeff 5 0.0 0.0            # H3-C3-H3-H3 (Zif8)
improper_coeff 6 0.0 0.0            # N-Zn-N-N (Zif8)
improper_coeff 7 2.0 180.0          # Zn-N-C1-C2 (Zif8)
