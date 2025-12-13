import numpy as np

# Lennard-Jones parameters for pair 1-2
epsilon = 0.070711      # kcal/mol
sigma = 2.449490        # Angstrom
q1 = 0.5 # e
q2 = -0.5 # e

# Atomic positions (Angstrom)
r1 = np.array([5.0, 0.0, 0.0])
r2 = np.array([-5.0, 0.0, 0.0])

# Distance
r = np.linalg.norm(r1 - r2)

# Lennard-Jones 12-6 potential
sr6 = (sigma / r) ** 6
sr12 = sr6 ** 2

E_LJ = 4.0 * epsilon * (sr12 - sr6)

# Output
print(f"E vdwl = {E_LJ:.12f} kcal/mol")

coulomb_const = 332.06371  # kcal·Å / (mol·e^2)
E_coul = coulomb_const * q1 * q2 / r

print(f"Coulomb energy = {E_coul:.12f} kcal/mol")

print(f"Total energy = {E_coul + E_LJ:.12f} kcal/mol")
