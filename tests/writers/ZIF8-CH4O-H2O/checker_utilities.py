import re

def read_lammps_data(filename):
    """
    Reads a LAMMPS data file and extracts atoms, bonds, angles, dihedrals, and impropers.
    Ignores Velocities and other non-essential sections.
    """
    atoms = {}
    bonds = []
    angles = []
    dihedrals = []
    impropers = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Section flags
    in_atoms = in_bonds = in_angles = in_dihedrals = in_impropers = False

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('!'):
            continue

        # Section headers
        if line.startswith("Atoms"):
            in_atoms = True
            in_bonds = in_angles = in_dihedrals = in_impropers = False
            continue
        elif line.startswith("Bonds"):
            in_bonds = True
            in_atoms = in_angles = in_dihedrals = in_impropers = False
            continue
        elif line.startswith("Angles"):
            in_angles = True
            in_atoms = in_bonds = in_dihedrals = in_impropers = False
            continue
        elif line.startswith("Dihedrals"):
            in_dihedrals = True
            in_atoms = in_bonds = in_angles = in_impropers = False
            continue
        elif line.startswith("Impropers"):
            in_impropers = True
            in_atoms = in_bonds = in_angles = in_dihedrals = False
            continue
        elif line.startswith("Velocities"):
            # Ignore velocities
            in_atoms = in_bonds = in_angles = in_dihedrals = in_impropers = False
            continue
        elif re.match(r"^(Masses|PairCoeffs|BondCoeffs|AngleCoeffs|DihedralCoeffs|ImproperCoeffs)", line):
            # Other sections to ignore
            in_atoms = in_bonds = in_angles = in_dihedrals = in_impropers = False
            continue

        # Parse sections
        if in_atoms:
            parts = line.split()
            atom_id = int(parts[0])
            atom_type = int(parts[2])
            atoms[atom_id] = atom_type

        elif in_bonds:
            parts = line.split()
            bond_type = int(parts[1])
            a1, a2 = int(parts[2]), int(parts[3])
            bonds.append((bond_type, a1, a2))

        elif in_angles:
            parts = line.split()
            angle_type = int(parts[1])
            a1, a2, a3 = int(parts[2]), int(parts[3]), int(parts[4])
            angles.append((angle_type, a1, a2, a3))

        elif in_dihedrals:
            parts = line.split()
            dihedral_type = int(parts[1])
            a1, a2, a3, a4 = int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])
            dihedrals.append((dihedral_type, a1, a2, a3, a4))

        elif in_impropers:
            parts = line.split()
            improper_type = int(parts[1])
            a1, a2, a3, a4 = int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])
            impropers.append((improper_type, a1, a2, a3, a4))

    return atoms, bonds, angles, dihedrals, impropers

def check_bonds(atoms, bonds, allowed_bonds):

    cpt_bond = 0

    for bond_type, a1, a2 in bonds:
        t1, t2 = atoms[a1], atoms[a2]

        # Ensure sorted order so (1,7) = (7,1)
        pair = tuple(sorted((t1, t2)))
        pair_rev = (t2, t1)  # bonds are undirected

        if bond_type in allowed_bonds:
            if pair not in allowed_bonds[bond_type] and pair_rev not in allowed_bonds[bond_type]:
                # Stop immediately if bond is invalid
                raise ValueError(f"Inconsistent bond found: "
                                 f"bond_type={bond_type}, atoms=({a1}-{a2}), "
                                 f"atom_types=({t1}-{t2}), allowed={allowed_bonds[bond_type]}")

        cpt_bond += 1

    print(cpt_bond, "bonds found")

def check_angles(atoms, angles, allowed_angles):
    """
    Stop immediately if any angle is inconsistent with allowed_angles.
    Prints a clean, readable message.
    """

    cpt_angle = 0

    for angle_type, a1, a2, a3 in angles:
        t1, t2, t3 = atoms[a1], atoms[a2], atoms[a3]
        triplet = (t1, t2, t3)
        triplet_rev = (t3, t2, t1)

        if angle_type in allowed_angles:
            if triplet not in allowed_angles[angle_type] and triplet_rev not in allowed_angles[angle_type]:
                allowed_str = ", ".join([f"({x},{y},{z})" for x,y,z in allowed_angles[angle_type]])
                raise ValueError(
                    f"\n Inconsistent angle found!\n"
                    f"  Angle type: {angle_type}\n"
                    f"  Atoms: {a1}-{a2}-{a3}\n"
                    f"  Atom types: {t1}-{t2}-{t3}\n"
                    f"  Allowed types: {allowed_str}\n"
                )
            
        cpt_angle += 1

    print(cpt_angle, "angles found")

def check_dihedrals(atoms, dihedrals, allowed_dihedrals):
    """
    Stop immediately if any dihedral is inconsistent with allowed_dihedrals.
    Dihedrals: 4 atoms, check exact and reversed (first/last swapped).
    """

    cpt_dihedral = 0

    for dih_type, a1, a2, a3, a4 in dihedrals:
        t1, t2, t3, t4 = atoms[a1], atoms[a2], atoms[a3], atoms[a4]
        quad = (t1, t2, t3, t4)
        quad_rev = (t4, t3, t2, t1)  # reverse symmetry

        if dih_type in allowed_dihedrals:
            if quad not in allowed_dihedrals[dih_type] and quad_rev not in allowed_dihedrals[dih_type]:
                allowed_str = ", ".join([f"({w},{x},{y},{z})" for w,x,y,z in allowed_dihedrals[dih_type]])
                raise ValueError(
                    f"\n Inconsistent dihedral found!\n"
                    f"  Dihedral type: {dih_type}\n"
                    f"  Atoms: {a1}-{a2}-{a3}-{a4}\n"
                    f"  Atom types: {t1}-{t2}-{t3}-{t4}\n"
                    f"  Allowed types: {allowed_str}\n"
                )
            
        cpt_dihedral += 1

    print(cpt_dihedral, "dihedrals found")

def check_impropers(atoms, impropers, allowed_impropers):
    """
    Stop immediately if any improper is inconsistent with allowed_impropers.
    Atom order is ignored (multiset comparison).
    """

    cpt_improper = 0

    for imp_type, a1, a2, a3, a4 in impropers:
        # Atom types
        t1, t2, t3, t4 = atoms[a1], atoms[a2], atoms[a3], atoms[a4]
        quad_set = frozenset([t1, t2, t3, t4])  # unordered comparison

        if imp_type in allowed_impropers:
            # Make sets for all allowed patterns
            allowed_sets = [frozenset(p) for p in allowed_impropers[imp_type]]
            if quad_set not in allowed_sets:
                allowed_str = ", ".join([f"{{{','.join(map(str,p))}}}" for p in allowed_impropers[imp_type]])
                raise ValueError(
                    f"\n Inconsistent improper found!\n"
                    f"  Improper type: {imp_type}\n"
                    f"  Atoms: {a1}-{a2}-{a3}-{a4}\n"
                    f"  Atom types: {t1}-{t2}-{t3}-{t4}\n"
                    f"  Allowed (unordered) sets: {allowed_str}\n"
                )
        
        cpt_improper += 1

    print(cpt_improper, "impropers found")

