"""End-to-end open-source structure-prep pipeline for OSPREY/RESISTOR.

Replaces the STAR-Protocol Maestro-dependent recipe with a fully-scripted,
license-free path. One command takes a PDB ID + ligand resname and emits
all artifacts needed for `osprey affinity --design ...`:

    {tag}.protein.pdb   {tag}.ligand.h.pdb   {tag}.ligand.prepi
    {tag}.ligand.frcmod {tag}.ligand.tc      {tag}.ligand.rot
    {tag}.affinity.yaml

Pipeline (each step is a thin shell-out so failures are localized):

    [1] fetch_pdb       (RCSB or local cache)
    [2] pdb4amber       AmberTools — adds missing atoms, normalizes residue numbering
    [3] p4a_undo        restore canonical chain/resnum (script from RESISTOR deposit)
    [4] split           Python: separate ATOM (protein) vs HETATM-of-resname (ligand)
    [5] protonate_lig   OpenBabel/RDKit replacement for Maestro/Epik
    [6] antechamber     AmberTools — generate .prepi at correct net charge
    [7] parmchk2        AmberTools — generate .frcmod
    [8] gen_templ_coords   shell or Python — write .tc atom record block
    [9] detect_rotamers RDKit replacement for Maestro dihedral picker → .rot
    [10] build_yaml     assemble OSPREY affinity YAML
    [11] verify_design  osprey affinity --verify-design

Required tooling (one conda env):
    ambertools  rdkit  openbabel  python>=3.11

Maestro is *not* required at any step.

Usage:
    python prep_complex.py \
        --pdb_id 4qta --ligand_resname 38Z \
        --tag erk2_sch7 --out_dir data/erk2_sch7
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

# Sibling scripts
HERE = Path(__file__).resolve().parent
P4A_UNDO = HERE / "p4a-undo.py"


def _run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    print(f"$ {' '.join(map(str, cmd))}")
    return subprocess.run(cmd, check=True, **kw)


# --------------------------------------------------------------------- #
def fetch_pdb(pdb_id: str, out: Path) -> None:
    pdb_id = pdb_id.lower()
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        return
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"[fetch] {url} → {out}")
    urllib.request.urlretrieve(url, out)


def run_pdb4amber(in_pdb: Path, out_pdb: Path, renum_txt: Path) -> None:
    if shutil.which("pdb4amber") is None:
        raise RuntimeError("pdb4amber not on PATH — activate the AmberTools env")
    _run(["pdb4amber", "--add-missing-atoms", "-i", str(in_pdb), "-o", str(out_pdb)])
    p4a_renum = out_pdb.with_suffix(out_pdb.suffix + "_renum.txt")
    if p4a_renum.exists() and p4a_renum != renum_txt:
        shutil.move(str(p4a_renum), str(renum_txt))


def run_p4a_undo(p4a_pdb: Path, renum_txt: Path, out_pdb: Path) -> None:
    if not P4A_UNDO.exists():
        raise FileNotFoundError(f"{P4A_UNDO} missing (RESISTOR Zenodo deposit)")
    with out_pdb.open("w") as f:
        _run(["python", str(P4A_UNDO), str(p4a_pdb), str(renum_txt)], stdout=f)


def split_protein_ligand(
    in_pdb: Path, ligand_resname: str, prot_out: Path, lig_out: Path
) -> None:
    """Split an OSPREY-prep'd PDB into ATOM-only (protein) and HETATM-of-resn
    (ligand). Anything else is dropped (waters, ions).
    """
    prot_lines, lig_lines = [], []
    for ln in in_pdb.read_text().splitlines():
        if ln.startswith("ATOM"):
            prot_lines.append(ln)
        elif ln.startswith("HETATM"):
            resn = ln[17:20].strip().upper()
            if resn == ligand_resname.upper():
                lig_lines.append(ln)
        elif ln.startswith("TER"):
            prot_lines.append(ln)
            lig_lines.append(ln)
    if not lig_lines:
        raise RuntimeError(
            f"no HETATM with resname {ligand_resname!r} found in {in_pdb}"
        )
    prot_out.write_text("\n".join(prot_lines + ["END"]) + "\n")
    lig_out.write_text("\n".join(lig_lines + ["END"]) + "\n")


def protonate_ligand(
    in_pdb: Path, out_pdb: Path, resname: str, pH: float = 7.4,
) -> int:
    _run([
        "python", str(HERE / "protonate_ligand.py"),
        "--in", str(in_pdb), "--out", str(out_pdb),
        "--resname", resname, "--pH", str(pH),
    ])
    sidecar = out_pdb.with_suffix(out_pdb.suffix + ".charge.json")
    return int(json.loads(sidecar.read_text())["net_charge"])


def run_antechamber(
    lig_h_in: Path, prepi_out: Path, net_charge: int, resname: str = "LIG",
) -> None:
    """antechamber drops scratch files (ATOMTYPE.INF, ANTECHAMBER_*.AC, sqm.*)
    into cwd, so we cd into prepi_out.parent and pass *basenames* so all
    artifacts land beside the output.

    `-rn <resname>` forces antechamber to write the supplied 3-letter residue
    name into the prepi output instead of its default 'UNL'. OSPREY's YAML
    ligand block must use the same resname.

    Input format auto-detected from suffix: prefer .mol2 (explicit bond block,
    needed for ligands with non-standard bonds like ATP/ANP P–N where bond
    perception from PDB coords falsely reports 'multiple units'). Falls back
    to .pdb for back-compat.
    """
    if shutil.which("antechamber") is None:
        raise RuntimeError("antechamber not on PATH — activate AmberTools env")
    suf = lig_h_in.suffix.lower()
    fi = {".mol2": "mol2", ".pdb": "pdb", ".sdf": "mdl"}.get(suf)
    if fi is None:
        raise RuntimeError(
            f"unsupported antechamber input suffix {suf!r} for {lig_h_in}"
        )
    work = prepi_out.parent.resolve()
    in_arg = (lig_h_in.resolve().name
              if lig_h_in.resolve().parent == work
              else str(lig_h_in.resolve()))
    base_cmd = [
        "antechamber",
        "-i", in_arg, "-fi", fi,
        "-o", prepi_out.name, "-fo", "prepi",
        "-nc", str(net_charge),
        "-rn", resname,
        "-pf", "y",
    ]
    # Try AM1-BCC first (highest accuracy). If sqm fails to converge — common
    # for highly-charged polyphosphates like ATP/ANP at -4 — fall back to
    # Gasteiger ('-c gas'). Gasteiger is less accurate but always converges
    # and is acceptable for prototype-level OSPREY scoring.
    try:
        _run(base_cmd + ["-c", "bcc"], cwd=str(work))
        return
    except Exception as e:
        msg = str(e).lower()
        if "sqm" not in msg and "maxcyc" not in msg and "converge" not in msg \
           and not (work / "sqm.out").exists():
            # not a sqm-related failure; re-raise
            raise
        print(f"[antechamber] AM1-BCC failed ({e!s}); retrying with Gasteiger.")
        _run(base_cmd + ["-c", "gas"], cwd=str(work))


def run_parmchk2(prepi_in: Path, frcmod_out: Path) -> None:
    work = frcmod_out.parent.resolve()
    in_arg = (prepi_in.resolve().name
              if prepi_in.resolve().parent == work
              else str(prepi_in.resolve()))
    _run([
        "parmchk2",
        "-i", in_arg, "-f", "prepi",
        "-a", "Y",
        "-o", frcmod_out.name,
    ], cwd=str(work))


def gen_templ_coords(lig_h_pdb: Path, resname: str, out_tc: Path) -> None:
    """Write OSPREY template-coords block: ``RESNAME N`` header followed by
    one ``atomName x y z`` line per ATOM/HETATM record. Uses fixed-width PDB
    columns (cols 13-16 atom name, 31-38 x, 39-46 y, 47-54 z) so it is robust
    to header lines and irregular whitespace."""
    rows: list[str] = []
    for ln in lig_h_pdb.read_text().splitlines():
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        if len(ln) < 54:
            continue
        atom = ln[12:16].strip()
        x, y, z = ln[30:38].strip(), ln[38:46].strip(), ln[46:54].strip()
        rows.append(f"{atom} {x} {y} {z}")
    out_tc.write_text(f"{resname} {len(rows)}\n" + "\n".join(rows) + "\n")


def detect_rotamers(lig_h_pdb: Path, resname: str, out_rot: Path) -> None:
    _run([
        "python", str(HERE / "detect_rotamers.py"),
        "--pdb", str(lig_h_pdb), "--resname", resname, "--out", str(out_rot),
    ])


# --------------------------------------------------------------------- #
def _residue_block(
    chain: str, res_num: int, aa_type: str,
    is_flexible: bool = True,
    use_continuous: bool = False,
    include_structure_rotamer: bool = True,
    mutability: list[str] | None = None,
) -> str:
    """Build one OSPREY residue_configurations entry."""
    mut = mutability or []
    mut_block = (
        "      mutability: []"
        if not mut
        else "      mutability:\n" + "\n".join(f"        - {a}" for a in mut)
    )
    return (
        "    - identity:\n"
        f"        chain: {chain}\n"
        f"        res_num: {res_num}\n"
        f"        aa_type: {aa_type}\n"
        "      flexibility:\n"
        f"        is_flexible: {str(is_flexible).lower()}\n"
        f"        use_continuous: {str(use_continuous).lower()}\n"
        f"        include_structure_rotamer: {str(include_structure_rotamer).lower()}\n"
        f"{mut_block}"
    )


def parse_flex_spec(spec: str) -> tuple[str, int, str, list[str]]:
    """Parse `CHAIN:RESNUM:AA_TYPE[:MUT1,MUT2,...]` flex specifier.

    Examples:
        "A:105:GLN"               -> ("A", 105, "GLN", [])
        "A:105:GLN:ALA,VAL"       -> ("A", 105, "GLN", ["ALA","VAL"])
    """
    parts = spec.split(":")
    if len(parts) < 3:
        raise ValueError(f"flex spec needs CHAIN:RESNUM:AA_TYPE, got {spec!r}")
    chain, resnum, aa_type = parts[0], int(parts[1]), parts[2].upper()
    mut = parts[3].upper().split(",") if len(parts) > 3 and parts[3] else []
    return chain, resnum, aa_type, mut


def build_yaml(
    tag: str,
    protein_pdb: Path,
    lig_h_pdb: Path,
    lig_tc: Path,
    lig_prepi: Path,
    lig_rot: Path,
    out_yaml: Path,
    epsilon: float = 0.99,
    protein_flex: list[tuple[str, int, str, list[str]]] | None = None,
    ligand_flex: list[tuple[str, int, str, list[str]]] | None = None,
    scan_res: list[tuple[str, int, str, list[str]]] | None = None,
) -> None:
    def _yaml_block(text: str, n: int) -> str:
        """Indent `text` for embedding under a YAML `|` literal-block scalar.

        snakeyaml infers the block's indentation level from the FIRST non-empty
        line. AmberTools `prepi` files have inconsistent leading whitespace
        (line 1 has 4 leading spaces; other lines start at column 0), which
        breaks naive uniform indentation. Pre-normalize: left-pad every
        non-empty line so they all share the same maximum leading-space count,
        then apply the outer indent.
        """
        pad = " " * n
        lines = text.splitlines()
        max_lead = max(
            (len(ln) - len(ln.lstrip(" ")) for ln in lines if ln.strip()),
            default=0,
        )
        out: list[str] = []
        for ln in lines:
            if not ln.strip():
                out.append("")
            else:
                cur = len(ln) - len(ln.lstrip(" "))
                out.append(pad + " " * (max_lead - cur) + ln)
        return "\n".join(out)

    def _format_rcfgs(specs):
        if not specs:
            return "  residue_configurations: []\n"
        blocks = "\n".join(
            _residue_block(c, n, aa, mutability=mut) for (c, n, aa, mut) in specs
        )
        return f"  residue_configurations:\n{blocks}\n"

    def _format_scan(specs):
        if not specs:
            return "scan:\n  residues: []\n"
        lines = ["scan:", "  residues:"]
        for c, n, aa, _mut in specs:
            lines += [
                "    - identity:",
                f"        chain: {c}",
                f"        res_num: {n}",
                f"        aa_type: {aa}",
            ]
        return "\n".join(lines) + "\n"

    yaml = (
        f"osprey_version: '3.3'\n"
        f"design_name: {tag}\n"
        f"epsilon: {epsilon}\n"
        f"{_format_scan(scan_res)}"
        f"protein:\n"
        f"{_format_rcfgs(protein_flex)}"
        f"  coordinates: |\n{_yaml_block(protein_pdb.read_text(), 4)}\n"
        f"ligand:\n"
        f"{_format_rcfgs(ligand_flex)}"
        f"  coordinates: |\n{_yaml_block(lig_h_pdb.read_text(), 4)}\n"
        f"  extra_templates_coordinates: |\n{_yaml_block(lig_tc.read_text(), 4)}\n"
        f"  extra_templates: |\n{_yaml_block(lig_prepi.read_text(), 4)}\n"
        f"  extra_rotamers: |\n{_yaml_block(lig_rot.read_text(), 4)}\n"
        f"  translate_rotate: true\n"
    )
    out_yaml.parent.mkdir(parents=True, exist_ok=True)
    out_yaml.write_text(yaml)


def verify_design(yaml: Path) -> bool:
    import os
    osprey = os.environ.get("OSPREY_BIN") or shutil.which("osprey")
    if osprey is None:
        print(
            "[verify] `osprey` not found on PATH and $OSPREY_BIN not set; "
            "skipping --verify-design",
            file=sys.stderr,
        )
        return True
    proc = subprocess.run(
        [osprey, "affinity", "--design", str(yaml), "--verify-design"],
        capture_output=True, text=True,
    )
    out = (proc.stdout or "") + (proc.stderr or "")
    print(out)
    return "Design file validated" in out


# --------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_id", required=True)
    ap.add_argument("--ligand_resname", required=True,
                    help="3-letter HETATM resname of the ligand (e.g. 38Z, ANP)")
    ap.add_argument("--tag", required=True, help="Output filename stem")
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--pH", type=float, default=7.4)
    ap.add_argument("--epsilon", type=float, default=0.99)
    ap.add_argument("--skip_verify", action="store_true")
    ap.add_argument(
        "--protein_flex", action="append", default=[],
        metavar="CHAIN:RES:AA[:MUT1,MUT2]",
        help="Repeatable: mark a protein residue as flexible. Optional :MUT,..."
             " makes it mutable to those AAs.")
    ap.add_argument(
        "--ligand_flex", action="append", default=[],
        metavar="CHAIN:RES:AA",
        help="Repeatable: mark a ligand residue as flexible (rotamer search).")
    ap.add_argument(
        "--scan_res", action="append", default=[],
        metavar="CHAIN:RES:AA",
        help="Repeatable: hotspot residue to put into scan.residues for "
             "use with `osprey affinity --do-scan`.")
    ap.add_argument("--skip_prep", action="store_true",
                    help="Skip everything before build_yaml; reuses cached "
                         "antechamber/parmchk2 outputs in --out_dir.")
    args = ap.parse_args()

    od = args.out_dir
    od.mkdir(parents=True, exist_ok=True)
    tag = args.tag

    raw_pdb     = od / f"{args.pdb_id.lower()}.pdb"
    p4a_pdb     = od / f"{args.pdb_id.lower()}.p4a.pdb"
    p4a_renum   = od / f"{args.pdb_id.lower()}.p4a_renum.txt"
    renum_pdb   = od / f"{args.pdb_id.lower()}.renum.pdb"
    protein_pdb = od / f"{tag}.protein.pdb"
    lig_pdb     = od / f"{tag}.ligand.pdb"
    lig_h_pdb   = od / f"{tag}.ligand.h.pdb"
    lig_h_mol2  = od / f"{tag}.ligand.h.mol2"   # written by protonate_ligand v2
    prepi       = od / f"{tag}.ligand.prepi"
    frcmod      = od / f"{tag}.ligand.frcmod"
    tc          = od / f"{tag}.ligand.tc"
    rot         = od / f"{tag}.ligand.rot"
    yaml_out    = od / f"{tag}.affinity.yaml"

    if not args.skip_prep:
        fetch_pdb(args.pdb_id, raw_pdb)
        run_pdb4amber(raw_pdb, p4a_pdb, p4a_renum)
        run_p4a_undo(p4a_pdb, p4a_renum, renum_pdb)
        split_protein_ligand(renum_pdb, args.ligand_resname, protein_pdb, lig_pdb)
        net_q = protonate_ligand(lig_pdb, lig_h_pdb,
                                 resname=args.ligand_resname, pH=args.pH)
        # Feed mol2 (explicit bond table) to antechamber, NOT pdb. Required for
        # ligands where bond perception from coords fails (ATP/ANP P–N bonds,
        # osi acryloyl). protonate_ligand v2 writes the sister mol2 file.
        run_antechamber(lig_h_mol2, prepi, net_q, resname=args.ligand_resname)
        run_parmchk2(prepi, frcmod)
        gen_templ_coords(lig_h_pdb, args.ligand_resname, tc)
        detect_rotamers(lig_h_pdb, args.ligand_resname, rot)

    prot_flex = [parse_flex_spec(s) for s in args.protein_flex]
    lig_flex  = [parse_flex_spec(s) for s in args.ligand_flex]
    scan_res  = [parse_flex_spec(s) for s in args.scan_res]
    build_yaml(tag, protein_pdb, lig_h_pdb, tc, prepi, rot, yaml_out,
               epsilon=args.epsilon,
               protein_flex=prot_flex, ligand_flex=lig_flex,
               scan_res=scan_res)

    print(f"\n[done] all artifacts written to {od}/{tag}.*")
    print(f"       YAML: {yaml_out}  (frcmod: {frcmod})")

    if not args.skip_verify:
        ok = verify_design(yaml_out)
        return 0 if ok else 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
