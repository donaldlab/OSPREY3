"""Open-source replacement for Maestro Step 1d (Epik), v2.

Pipeline (template-driven, NOT pKa-heuristic):
  1. Resolve a SMILES *template* for the ligand:
       - --smiles "..." (highest priority), else
       - cached/fetched RCSB CCD SMILES_CANONICAL for `--resname`.
     The template is the authoritative source of bond orders + formal charges.
  2. RDKit reads the input PDB for COORDINATES ONLY (sanitize=False, no H).
  3. AssignBondOrdersFromTemplate(template, pdb_mol) overwrites whatever bond
     orders/charges RDKit guessed from PDB geometry with the SMILES truth.
     -> osimertinib's acryloyl C=C survives, ATP/ANP polyphosphate stays -4,
        ANP's β–γ P–N bond is preserved.
  4. Chem.AddHs(mol, addCoords=True) places hydrogens by SMILES-implied
     valence. No pH/pKa adjustment is performed here; protonation state IS
     the SMILES template you fed in. (Use Dimorphite-DL upstream if you want
     pH-7.4 enumeration.)
  5. Write THREE outputs side-by-side:
       <out>           protonated PDB           (back-compat: existing readers)
       <out>.mol2      protonated mol2          (NEW: explicit bond table for
                                                 antechamber -fi mol2)
       <out>.charge.json  sidecar metadata
     antechamber must be invoked with `-fi mol2` to get the bond table; -fi
     pdb makes it re-guess bonds from coords and falsely flag ANP P–N
     ("multiple units" error).

CLI:
    python protonate_ligand.py \
        --in lig.pdb --out lig.h.pdb --resname ANP \
        [--smiles "<SMILES>"]            # override CCD lookup
        [--charge_override -4]           # bypass RDKit charge readout
        [--ccd_cache ~/.cache/evotraj/ccd]
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

DEFAULT_CCD_CACHE = Path(os.environ.get(
    "EVOTRAJ_CCD_CACHE",
    str(Path("~/.cache/evotraj/ccd").expanduser()),
))


# ---------------------------------------------------------------- CCD lookup
def fetch_ccd_smiles(resname: str, cache_dir: Path = DEFAULT_CCD_CACHE) -> str:
    """Return canonical SMILES for a 3-letter PDB resname from the RCSB CCD.

    Resolution order:
      1. cache_dir/<RESNAME>.smi  (one SMILES per file)
      2. https://files.rcsb.org/ligands/view/<RESNAME>.cif → parse
         _pdbx_chem_comp_descriptor for type=SMILES_CANONICAL preferring
         program=OpenEye (RCSB-canonical), falling back to RDKit then any.
    """
    resname = resname.upper().strip()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / f"{resname}.smi"
    if cached.exists():
        s = cached.read_text().strip()
        if s:
            return s

    url = f"https://files.rcsb.org/ligands/view/{resname}.cif"
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            cif = r.read().decode("utf-8", errors="replace")
    except Exception as e:
        raise RuntimeError(
            f"Could not fetch CCD for resname={resname!r} from {url}: {e}. "
            f"Pass --smiles explicitly or place SMILES at {cached}."
        )

    smiles = _parse_ccd_smiles(cif)
    if not smiles:
        raise RuntimeError(
            f"No SMILES_CANONICAL row found in CCD for {resname!r}. "
            f"Pass --smiles explicitly."
        )
    cached.write_text(smiles + "\n")
    return smiles


def _parse_ccd_smiles(cif_text: str) -> str:
    """Extract SMILES_CANONICAL from a CCD CIF. Handles the loop_ form used
    by RCSB ligand CIFs:
        loop_
        _pdbx_chem_comp_descriptor.comp_id
        _pdbx_chem_comp_descriptor.type
        _pdbx_chem_comp_descriptor.program
        _pdbx_chem_comp_descriptor.program_version
        _pdbx_chem_comp_descriptor.descriptor
        ANP SMILES           ACDLabs    10.04 "..."
        ANP SMILES_CANONICAL "OpenEye OEToolkits" 1.5.0 "..."
        ...
    """
    import shlex
    lines = cif_text.splitlines()
    i = 0
    candidates: list[tuple[str, str, str]] = []   # (type, program, descriptor)
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1; continue
        # collect column headers for this loop
        cols: list[str] = []
        i += 1
        while i < len(lines) and lines[i].lstrip().startswith("_"):
            cols.append(lines[i].strip())
            i += 1
        if "_pdbx_chem_comp_descriptor.descriptor" not in cols:
            continue
        idx_type = cols.index("_pdbx_chem_comp_descriptor.type")
        idx_prog = cols.index("_pdbx_chem_comp_descriptor.program")
        idx_desc = cols.index("_pdbx_chem_comp_descriptor.descriptor")
        # data rows until next '#' or 'loop_' or '_xxx' header
        while i < len(lines):
            row = lines[i]
            s = row.strip()
            if not s or s.startswith("#") or s == "loop_" or s.startswith("_"):
                break
            try:
                tok = shlex.split(row, posix=True)
            except ValueError:
                tok = row.split()
            if len(tok) >= len(cols):
                candidates.append((tok[idx_type], tok[idx_prog], tok[idx_desc]))
            i += 1

    # Prefer OpenEye SMILES_CANONICAL, then RDKit, then any SMILES_CANONICAL,
    # then any SMILES.
    def pick(pred):
        for t, p, d in candidates:
            if pred(t, p):
                return d
        return ""
    return (
        pick(lambda t, p: t == "SMILES_CANONICAL" and "OpenEye" in p)
        or pick(lambda t, p: t == "SMILES_CANONICAL" and "RDKit" in p)
        or pick(lambda t, p: t == "SMILES_CANONICAL")
        or pick(lambda t, p: t == "SMILES")
    )


# ---------------------------------------------------------------- pH 7.4 rules
# Minimal pH-7.4 deprotonation/protonation rules applied to the SMILES
# template BEFORE bond-order transfer. Heavy-atom count is preserved (we only
# flip H counts + formal charges on hetero atoms), so the
# AssignBondOrdersFromTemplate skeleton match still works. Each rule is a
# (SMARTS, atom_idx_in_match, new_charge, new_n_explicit_H_or_None_for_+1).
PH74_RULES = [
    # Phosphate / phosphonate OH (ATP α/β/γ, ANP β-O, AMP, etc.). pKa1<<2,
    # pKa2~7 ⇒ deprotonate every OH on a P=O.
    ("[OX2H1;$([OH]-[PX4]=O)]", 0, -1, 0),
    # Sulfate / sulfonate OH:
    ("[OX2H1;$([OH]-[SX4](=O)(=O))]", 0, -1, 0),
    # Carboxylic acid:
    ("[OX2H1;$([OH]-C(=O))]", 0, -1, 0),
    # Tetrazole NH (pKa ~ 4.9):
    ("[nH;$(n1nnnc1),$(n1cnnn1)]", 0, -1, 0),
    # Phosphoramidate N–H (ANP β–γ): pKa~8 ⇒ leave protonated, no rule.
    # Aliphatic primary/secondary amine: pKa ~ 9-10 ⇒ protonated. Skip
    # aromatic, amide, imine, nitrile, aniline-style, AND phosphoramidate
    # (N–P, e.g. ANP β–γ) / sulfonamide (N–S).
    ("[NX3;H2,H1;!$(NC=O);!$(N=*);!$(N#*);!$(N-a);!$(N-[PX4]);!$(N-[SX4])]",
     0, +1, None),
]


def apply_ph74_rules(template_mol):
    """Return a NEW Mol with pH-7.4 deprotonations/protonations applied.
    We mutate a copy, then SMILES-round-trip to get a clean RDKit Mol with
    consistent valence/aromaticity state. (Mutating in place breaks
    AssignBondOrdersFromTemplate's substructure search for some inputs —
    e.g. ANP — because partial sanitize leaves the cache inconsistent.)
    """
    from rdkit import Chem
    rw = Chem.RWMol(template_mol)
    n_changed = 0
    for smarts, idx, q, nH in PH74_RULES:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        for match in rw.GetSubstructMatches(patt):
            atom = rw.GetAtomWithIdx(match[idx])
            atom.SetFormalCharge(q)
            if nH is None:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
            else:
                atom.SetNumExplicitHs(nH)
            atom.SetNoImplicit(True)
            n_changed += 1
    if n_changed == 0:
        return template_mol, 0
    Chem.SanitizeMol(rw)
    smi2 = Chem.MolToSmiles(rw)
    fresh = Chem.MolFromSmiles(smi2)
    if fresh is None:
        # rules produced an invalid SMILES; back off to the original template
        return template_mol, 0
    return fresh, n_changed


# ---------------------------------------------------------------- core chem
def load_with_template(pdb_path: Path, smiles: str, apply_ph74: bool = True):
    """Read coords from PDB; transfer bond orders + formal charges from SMILES.

    Returns an RDKit Mol with NO hydrogens but correct heavy-atom topology.
    Raises if heavy-atom count or skeleton mismatch (catches CCD/SMILES
    version drift early — better than producing a silently-wrong ligand).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    pdb_mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=True, sanitize=False)
    if pdb_mol is None:
        raise RuntimeError(f"RDKit could not parse PDB: {pdb_path}")
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise RuntimeError(f"RDKit could not parse template SMILES: {smiles!r}")
    if apply_ph74:
        template, _n = apply_ph74_rules(template)

    n_pdb = pdb_mol.GetNumHeavyAtoms()
    n_tpl = template.GetNumHeavyAtoms()
    if n_pdb != n_tpl:
        raise RuntimeError(
            f"Heavy-atom count mismatch: PDB={n_pdb}, template SMILES={n_tpl}. "
            f"PDB may be missing atoms or template is wrong protomer/tautomer. "
            f"Pass a corrected --smiles."
        )
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
    except Exception as e:
        raise RuntimeError(
            f"AssignBondOrdersFromTemplate failed: {e}. "
            f"Heavy-atom skeleton in PDB likely doesn't match SMILES template."
        )
    return mol


def add_hydrogens(mol):
    """Add H by SMILES-implied valence with 3D coords. No pKa adjustment."""
    from rdkit import Chem
    return Chem.AddHs(mol, addCoords=True)


# ---------------------------------------------------------------- I/O
def write_pdb(mol, out_pdb: Path) -> None:
    from rdkit import Chem
    Chem.MolToPDBFile(mol, str(out_pdb))


def write_mol2_via_obabel(in_pdb: Path, out_mol2: Path) -> None:
    """RDKit doesn't ship a mol2 writer; round-trip the protonated PDB through
    obabel. Because the PDB now has correct bonds (CONECT records emitted by
    RDKit MolToPDBFile), obabel preserves them in the mol2 bond block."""
    if shutil.which("obabel") is None:
        raise RuntimeError("`obabel` not found on PATH "
                           "(conda install -c conda-forge openbabel)")
    proc = subprocess.run(
        ["obabel", "-ipdb", str(in_pdb), "-omol2", "-O", str(out_mol2)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0 or not out_mol2.exists():
        raise RuntimeError(
            f"obabel pdb→mol2 failed (exit {proc.returncode}):\n"
            f"  stdout: {proc.stdout}\n  stderr: {proc.stderr}"
        )


def normalize_atom_names(
    pdb: Path,
    resname: str,
    chain: str = "L",
    resnum: int = 1,
) -> None:
    """Rewrite atom names to <ELEMENT><counter> and force resname/chain/resnum.
    OSPREY/antechamber require unique atom names per residue.

    Also fixes a RDKit MolToPDBFile quirk: the element column (77-78) may
    contain a 3-char string like 'O1-' (element + formal-charge notation)
    when atoms have non-zero formal charge. OSPREY parses this as atom type
    "O1" and fails forcefield lookup. We split into a clean 2-char element
    (cols 77-78) and a separate 2-char charge field (cols 79-80, e.g. '1-').
    """
    import re
    counters: dict[str, int] = {}
    out_lines: list[str] = []
    chain_c = (chain or "L")[0]
    res_field = f"{int(resnum):>4d}"
    chg_re = re.compile(r"([A-Z][a-z]?)(\d?[+-]?)")
    for ln in pdb.read_text().splitlines():
        if not ln.startswith(("HETATM", "ATOM")):
            out_lines.append(ln)
            continue
        if len(ln) < 80:
            ln = ln.ljust(80)
        # Cols 77-80: parse "<ELEM><CHARGE>" — RDKit may pack "O1-" into 77-79.
        elem_field_raw = ln[76:80].strip().upper()
        m = chg_re.match(elem_field_raw)
        if m:
            element = m.group(1)
            charge = m.group(2) or ""
        else:
            element = elem_field_raw[:2].rstrip()
            charge = ""
        if not element:
            element = ln[12:14].strip()[:1].upper() or "X"
        counters[element] = counters.get(element, 0) + 1
        new_name = f"{element}{counters[element]}"
        if len(element) == 1 and len(new_name) <= 3:
            atom_field = f" {new_name:<3s}"
        else:
            atom_field = f"{new_name:<4s}"
        # Re-pad element (right-justified in 2 cols) + charge (right-justified
        # in 2 cols). Standard PDB format.
        elem_col = f"{element:>2s}"
        chg_col = f"{charge:>2s}" if charge else "  "
        new = (
            ln[:12]
            + atom_field          # 13-16: atom name
            + ln[16]              # 17: altLoc
            + f"{resname:>3s}"    # 18-20: resname
            + " "                 # 21
            + chain_c             # 22: chain
            + res_field           # 23-26: resnum
            + ln[26:76]           # 27-76: iCode + xyz + occ + bfactor
            + elem_col            # 77-78: element
            + chg_col             # 79-80: charge
        )
        out_lines.append(new.rstrip())
    pdb.write_text("\n".join(out_lines) + "\n")


def patch_mol2_resname(mol2: Path, resname: str) -> None:
    """obabel writes 'LIG1' or 'UNL1' as the SUBST_NAME. Rewrite to the
    requested resname so antechamber's `-rn` matches."""
    lines = mol2.read_text().splitlines()
    out: list[str] = []
    in_atom = False
    for ln in lines:
        if ln.startswith("@<TRIPOS>"):
            in_atom = (ln.strip() == "@<TRIPOS>ATOM")
            out.append(ln); continue
        if not in_atom or not ln.strip() or ln.startswith("#"):
            out.append(ln); continue
        # mol2 ATOM: id name x y z type subst_id subst_name charge ...
        parts = ln.split()
        if len(parts) >= 8:
            parts[7] = f"{resname}1"
            out.append(" ".join(parts))
        else:
            out.append(ln)
    mol2.write_text("\n".join(out) + "\n")


def sniff_chain_resnum(pdb: Path) -> tuple[str, int]:
    for ln in pdb.read_text().splitlines():
        if ln.startswith(("HETATM", "ATOM")) and len(ln) >= 26:
            chain = ln[21] if ln[21].strip() else "L"
            try:
                rnum = int(ln[22:26].strip())
            except ValueError:
                rnum = 1
            return chain, rnum
    return "L", 1


def sniff_resname(pdb: Path) -> str | None:
    for ln in pdb.read_text().splitlines():
        if ln.startswith(("HETATM", "ATOM")) and len(ln) >= 20:
            r = ln[17:20].strip()
            if r:
                return r
    return None


# ---------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", type=Path, required=True,
                    help="Input ligand PDB (heavy atoms; H optional)")
    ap.add_argument("--out", dest="out", type=Path, required=True,
                    help="Output protonated PDB. A sister .mol2 with explicit "
                         "bond block is written next to it for antechamber.")
    ap.add_argument("--resname", type=str, default=None,
                    help="3-letter PDB resname (e.g. ANP, ATP, AQ4). Inferred "
                         "from input if omitted.")
    ap.add_argument("--smiles", type=str, default=None,
                    help="Override SMILES template (skip CCD lookup). Use this "
                         "to pin a specific protomer/tautomer at pH 7.4.")
    ap.add_argument("--charge_override", type=int, default=None,
                    help="Force net charge instead of summing GetFormalCharge() "
                         "(only for exotic multimetal/radical edge cases).")
    ap.add_argument("--ccd_cache", type=Path, default=DEFAULT_CCD_CACHE,
                    help=f"CCD SMILES cache dir (default: {DEFAULT_CCD_CACHE})")
    ap.add_argument("--charge_json", type=Path, default=None,
                    help="Sidecar JSON path; defaults to <out>.charge.json")
    ap.add_argument("--no_pH74", action="store_true",
                    help="Disable the built-in pH-7.4 SMARTS deprotonation "
                         "rules (phosphate/sulfate/carboxylate→anion, "
                         "aliphatic amine→cation). Use when you pass a "
                         "--smiles already in the desired protonation state.")
    # Back-compat: accept (and ignore value of) --pH; v2 has its own pH model.
    ap.add_argument("--pH", type=float, default=7.4,
                    help="Compatibility flag. v2 only knows pH 7.4; pass "
                         "--no_pH74 to skip the rules entirely.")
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    if args.resname is None:
        args.resname = sniff_resname(args.inp)
        if args.resname is None:
            raise RuntimeError(f"Could not infer resname from {args.inp}; pass --resname.")

    chain, resnum = sniff_chain_resnum(args.inp)

    # 1. Resolve template SMILES
    smiles = args.smiles or fetch_ccd_smiles(args.resname, args.ccd_cache)

    # 2-4. Load + bond-order transfer + add Hs
    mol = load_with_template(args.inp, smiles, apply_ph74=not args.no_pH74)
    mol = add_hydrogens(mol)

    # 5. Write outputs
    out_pdb  = args.out
    out_mol2 = out_pdb.with_suffix(".mol2") if out_pdb.suffix else Path(str(out_pdb) + ".mol2")
    write_pdb(mol, out_pdb)
    normalize_atom_names(out_pdb, args.resname, chain=chain, resnum=resnum)
    write_mol2_via_obabel(out_pdb, out_mol2)
    patch_mol2_resname(out_mol2, args.resname)

    # Net charge: trust template unless overridden
    if args.charge_override is not None:
        net_q = int(args.charge_override)
    else:
        net_q = sum(a.GetFormalCharge() for a in mol.GetAtoms())

    n_atoms = mol.GetNumAtoms()
    n_heavy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
    from rdkit import Chem
    canon_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))

    sidecar = args.charge_json or Path(str(args.out) + ".charge.json")
    sidecar.write_text(json.dumps({
        "input": str(args.inp),
        "output_pdb": str(out_pdb),
        "output_mol2": str(out_mol2),
        "resname": args.resname,
        "template_smiles": smiles,
        "smiles_source": "user_override" if args.smiles else "rcsb_ccd",
        "ph74_rules_applied": (not args.no_pH74),
        "net_charge": int(net_q),
        "charge_source": "override" if args.charge_override is not None else "template",
        "n_atoms": n_atoms,
        "n_heavy": n_heavy,
        "smiles": canon_smiles,
    }, indent=2))

    print(f"[protonate v2] {args.inp.name}  resname={args.resname}  "
          f"charge={net_q:+d}  heavy={n_heavy}")
    print(f"               pdb : {out_pdb}")
    print(f"               mol2: {out_mol2}  (bond table for antechamber -fi mol2)")
    print(f"               json: {sidecar}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
