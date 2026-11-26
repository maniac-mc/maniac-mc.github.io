#!/usr/bin/env python3
"""LAMMPS-compatibility validator for data files written by MANIAC.

Usage:
    python3 check_output.py outputs/topology.data

Checks that the file is a well-formed LAMMPS data file (atom_style full):
  - header counts match the actual number of entries in each section
    (a missing header line counts as zero)
  - no duplicate lines inside Masses / Atoms
  - box bounds exist, are finite and ordered (lo < hi), optional tilt finite
  - atom ids unique and positive; atom types in [1, natomtypes];
    mol ids non-negative; coordinates finite
  - one positive finite mass per atom type
  - bond/angle/dihedral/improper entries reference existing atoms and
    valid type ids

Exit code 0 = OK, 1 = violation (details on stderr).
No reference file is needed; this only validates the output shape.
"""
import argparse
import math
import re
import sys

DATA_SECTIONS = {'Masses', 'Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers'}
HEADER_RE = re.compile(r'^(\d+)\s+(.+)$')
BOUNDS_RE = re.compile(r'^(\S+)\s+(\S+)\s+(xlo xhi|ylo yhi|zlo zhi)$')
TILT_RE = re.compile(r'^(\S+)\s+(\S+)\s+(\S+)\s+xy xz yz$')
NUM_RE = re.compile(r'^[-+]?[0-9.]')


def parse_data_file(path):
    d = dict(header={}, bounds={}, tilt=None, masses={}, n_mass_lines=0,
             atoms={}, n_atom_lines=0,
             bonds=[], angles=[], dihedrals=[], impropers=[])
    section = None
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line[0] in '!#':
                continue
            m = BOUNDS_RE.match(line)
            if m:
                d['bounds'][m.group(3)] = (float(m.group(1)), float(m.group(2)))
                continue
            m = TILT_RE.match(line)
            if m:
                d['tilt'] = tuple(float(m.group(i)) for i in (1, 2, 3))
                continue
            if section is None:
                m = HEADER_RE.match(line)
                if m:
                    d['header'][m.group(2).strip()] = int(m.group(1))
                    continue
            if not NUM_RE.match(line):
                # section header (e.g. "Atoms # full", "Pair Coeffs # ...")
                key = line.split('#')[0].strip()
                section = key if key in DATA_SECTIONS else None
                continue
            p = line.split()
            if section == 'Masses':
                d['n_mass_lines'] += 1
                d['masses'][int(p[0])] = float(p[1])
            elif section == 'Atoms':
                d['n_atom_lines'] += 1
                d['atoms'][int(p[0])] = dict(mol=int(p[1]), type=int(p[2]),
                                             q=float(p[3]),
                                             xyz=tuple(float(v) for v in p[4:7]))
            elif section == 'Bonds':
                d['bonds'].append(tuple(int(v) for v in p))
            elif section == 'Angles':
                d['angles'].append(tuple(int(v) for v in p))
            elif section == 'Dihedrals':
                d['dihedrals'].append(tuple(int(v) for v in p))
            elif section == 'Impropers':
                d['impropers'].append(tuple(int(v) for v in p))
            # other sections (Velocities, coeffs): ignored
    return d


def check_structure(d):
    errs = []
    h = d['header']
    ntypes = h.get('atom types', 0)

    # A missing header line means zero (some writers omit "0 dihedrals" lines)
    for key, attr in [('atoms', 'atoms'), ('bonds', 'bonds'), ('angles', 'angles'),
                      ('dihedrals', 'dihedrals'), ('impropers', 'impropers')]:
        if h.get(key, 0) != len(d[attr]):
            errs.append(f"header declares {h.get(key, 0)} {key} but section has {len(d[attr])}")
    if h.get('atom types', 0) != len(d['masses']):
        errs.append(f"header declares {h.get('atom types', 0)} atom types "
                    f"but Masses has {len(d['masses'])}")
    if d['n_mass_lines'] != len(d['masses']):
        dup = d['n_mass_lines'] - len(d['masses'])
        errs.append(f"Masses section has {d['n_mass_lines']} lines but only "
                    f"{len(d['masses'])} unique atom types ({dup} duplicate entries)")

    for ax in ('xlo xhi', 'ylo yhi', 'zlo zhi'):
        if ax not in d['bounds']:
            errs.append(f"missing box bounds line '{ax}'")
        else:
            lo, hi = d['bounds'][ax]
            if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
                errs.append(f"bad bounds '{ax}': lo={lo}, hi={hi}")
    if d['tilt'] is not None and not all(math.isfinite(v) for v in d['tilt']):
        errs.append(f"non-finite tilt factors: {d['tilt']}")

    ids = list(d['atoms'])
    if d['n_atom_lines'] != len(ids):
        errs.append(f"Atoms section has {d['n_atom_lines']} lines but "
                    f"{len(ids)} unique atom ids (duplicate atom ids present)")
    for aid, a in d['atoms'].items():
        if aid < 1:
            errs.append(f"atom id {aid} is not positive")
            break
        if not 1 <= a['type'] <= ntypes:
            errs.append(f"atom {aid} has type {a['type']} outside [1, {ntypes}]")
            break
        if a['mol'] < 0:
            errs.append(f"atom {aid} has negative mol id {a['mol']}")
            break
        if not all(math.isfinite(v) for v in a['xyz']):
            errs.append(f"atom {aid} has non-finite coordinates")
            break

    for t, m in d['masses'].items():
        if not 1 <= t <= ntypes or not math.isfinite(m) or m <= 0:
            errs.append(f"bad mass for atom type {t}: {m}")
            break

    for name, entries, nref, ntypes_key in [
            ('bond', d['bonds'], 2, 'bond types'),
            ('angle', d['angles'], 3, 'angle types'),
            ('dihedral', d['dihedrals'], 4, 'dihedral types'),
            ('improper', d['impropers'], 4, 'improper types')]:
        ntk = h.get(ntypes_key, 0)
        for e in entries:
            if not 1 <= e[1] <= ntk:
                errs.append(f"{name} {e[0]} has type {e[1]} outside [1, {ntk}]")
                break
            bad = [ref for ref in e[2:2 + nref] if ref not in d['atoms']]
            if bad:
                errs.append(f"{name} {e[0]} references missing atom id {bad[0]}")
                break
    return errs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('output', help='LAMMPS data file written by MANIAC (e.g. outputs/topology.data)')
    args = ap.parse_args()

    d = parse_data_file(args.output)
    errs = check_structure(d)

    if errs:
        for e in errs:
            print(f"FAIL {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()