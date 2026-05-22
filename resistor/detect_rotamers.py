"""Open-source replacement for Maestro Step 5 (rotamer/dihedral picker).

Auto-detects rotatable dihedrals in a small-molecule ligand and writes them
in OSPREY's `.rot` format (Figure 3 in the STAR Protocol).

Algorithm
---------
1. RDKit reads the protonated ligand PDB.
2. For each rotatable bond (default RotatableBondSmarts: any non-ring single
   bond between two heavy atoms with > 1 heavy neighbor), pick the two
   "pivot" atoms b–c and the most-connected neighbors a (of b) and d (of c)
   to form the dihedral a–b–c–d.
3. Use the input PDB coordinates to compute the current dihedral angle; the
   `.rot` file lists each dihedral on its own line, then a final line with
   one rotamer's angles (rounded to whole degrees).
4. Output format mirrors Figure 3 exactly so it can be dropped into
   `extra_rotamers` in the YAML.

Limitations vs Maestro
----------------------
- RDKit's default rotatable-bond SMARTS may include some bonds Maestro
  would exclude (amide-like). We additionally filter out:
    * amide C–N (rotation chemically restricted)
    * triple bonds, double bonds (already filtered by SMARTS, just safety)
    * bonds to terminal H/D
- For ATP analogs / multi-ring pharmacophores you should sanity-check the
  output by visualizing in PyMOL: every red bond in `pymol> show sticks` +
  the dihedral list should be 1:1.
- We pick *one* current-conformer rotamer. Maestro's workflow can include
  scanning multiple rotamers per dihedral; for OSPREY's continuous-flex
  voxel search a single seed is usually enough.

Usage:
    python detect_rotamers.py --pdb amppnp.h.pdb --resname ANP --out amppnp.rot
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


# Default SMARTS: any single bond, both atoms in non-ring or non-aromatic env,
# both atoms have >1 heavy neighbor (i.e. non-terminal). Same convention as
# RDKit's NumRotatableBonds + SMILESLipinski.
DEFAULT_ROTSMARTS = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"

# Excluded patterns (chemically restricted bonds OpenBabel/Maestro skip):
EXCLUDE_SMARTS = [
    "[NX3;H2,H1;!$(NC=O)]-[CX3]=O",   # amide N-C(=O)
    "[NX3]-[CX3](=[OX1])",            # amide N-C
    "[#6]=[#6]",                      # any C=C (safety; SMARTS already excludes)
]


def load_mol(pdb: Path):
    from rdkit import Chem
    mol = Chem.MolFromPDBFile(str(pdb), removeHs=False, sanitize=False)
    if mol is None:
        raise RuntimeError(f"RDKit could not parse {pdb}")
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    return mol


def get_pdb_atom_names(mol) -> list[str]:
    """Return PDB-style atom names per atom index."""
    names = []
    for at in mol.GetAtoms():
        info = at.GetPDBResidueInfo()
        if info is not None:
            names.append(info.GetName().strip())
        else:
            # Fallback: element + index (rarely hit if PDB is well-formed)
            names.append(f"{at.GetSymbol()}{at.GetIdx()}")
    return names


def find_rotatable_dihedrals(mol) -> list[tuple[int, int, int, int]]:
    from rdkit import Chem
    patt = Chem.MolFromSmarts(DEFAULT_ROTSMARTS)
    bond_atom_pairs = mol.GetSubstructMatches(patt)

    excl_patts = [Chem.MolFromSmarts(s) for s in EXCLUDE_SMARTS]
    excl_bonds: set[frozenset[int]] = set()
    for ep in excl_patts:
        if ep is None:
            continue
        for match in mol.GetSubstructMatches(ep):
            for i in range(len(match) - 1):
                excl_bonds.add(frozenset((match[i], match[i + 1])))

    dihedrals: list[tuple[int, int, int, int]] = []
    seen_bonds: set[frozenset[int]] = set()
    for b, c in bond_atom_pairs:
        bond_key = frozenset((b, c))
        if bond_key in seen_bonds or bond_key in excl_bonds:
            continue
        seen_bonds.add(bond_key)

        atom_b = mol.GetAtomWithIdx(b)
        atom_c = mol.GetAtomWithIdx(c)

        def best_neighbor(atom, exclude_idx):
            cands = [n for n in atom.GetNeighbors() if n.GetIdx() != exclude_idx]
            cands = [n for n in cands if n.GetAtomicNum() > 1] or cands
            if not cands:
                return None
            # prefer heavy + most connected
            cands.sort(key=lambda n: (n.GetAtomicNum() <= 1, -n.GetDegree()))
            return cands[0].GetIdx()

        a = best_neighbor(atom_b, c)
        d = best_neighbor(atom_c, b)
        if a is None or d is None:
            continue
        dihedrals.append((a, b, c, d))
    return dihedrals


def measure_dihedral(mol, a: int, b: int, c: int, d: int) -> float:
    from rdkit.Chem import rdMolTransforms
    conf = mol.GetConformer(0)
    return float(rdMolTransforms.GetDihedralDeg(conf, a, b, c, d))


def write_rot(
    out: Path,
    resname: str,
    dihedrals: list[tuple[int, int, int, int]],
    angles_deg: list[float],
    atom_names: list[str],
) -> None:
    """Write OSPREY .rot per Figure 3.

        ! comments...
        1                       (1 AA type)
        ANP 10 1                (resname num_dihedrals num_rotamers)
        OG1 PG N3B PB           (10 dihedral lines)
        ...
        -173 -81 -52 ... 107 -4 (one row of integer angles)
    """
    lines = [
        "! Auto-generated by detect_rotamers.py (RDKit)",
        "! Format: AA_name num_dihedrals num_rotamers",
        "!         <num_dihedrals> dihedral lines (atom-name quads)",
        "!         one rotamer-angle row (whole-degree integers)",
        "1",
        f"{resname} {len(dihedrals)} 1",
    ]
    for a, b, c, d in dihedrals:
        lines.append(f"{atom_names[a]} {atom_names[b]} {atom_names[c]} {atom_names[d]}")
    lines.append(" ".join(f"{int(round(x)):d}" for x in angles_deg))
    out.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=Path, required=True,
                    help="Protonated ligand PDB (output of protonate_ligand.py)")
    ap.add_argument("--resname", type=str, required=True,
                    help="3-letter PDB residue name of the ligand (e.g. ANP, 38Z)")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output .rot path")
    args = ap.parse_args()

    mol = load_mol(args.pdb)
    atom_names = get_pdb_atom_names(mol)
    dihedrals = find_rotatable_dihedrals(mol)

    if not dihedrals:
        print(f"[detect_rotamers] WARNING: no rotatable dihedrals found in {args.pdb}",
              file=sys.stderr)

    angles = [measure_dihedral(mol, a, b, c, d) for (a, b, c, d) in dihedrals]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_rot(args.out, args.resname, dihedrals, angles, atom_names)

    print(f"[detect_rotamers] {args.pdb} → {args.out}")
    print(f"                  {len(dihedrals)} rotatable dihedrals; "
          f"angles (deg, rounded) = "
          f"{[int(round(x)) for x in angles]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
