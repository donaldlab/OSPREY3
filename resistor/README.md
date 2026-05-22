# RESISTOR open-source structure-prep pipeline

Drop-in, license-free replacement for the Maestro-dependent prep steps
in the OSPREY/RESISTOR STAR Protocol. Given a PDB ID + ligand resname,
emits every artifact needed by `osprey affinity --design ...`:

```
{tag}.protein.pdb   {tag}.ligand.h.pdb   {tag}.ligand.prepi
{tag}.ligand.frcmod {tag}.ligand.tc      {tag}.ligand.rot
{tag}.affinity.yaml
```

## Pipeline

`prep_complex.py` is a thin orchestrator; each step shells out so failures
are localised:

| Step | Tool | Replaces |
|------|------|----------|
| 1. fetch PDB        | `urllib` (RCSB)       | manual download |
| 2. normalize        | `pdb4amber`           | Maestro Protein Prep |
| 3. undo renumber    | `p4a-undo.py`         | (same as RESISTOR) |
| 4. split prot/lig   | Python                | (same) |
| 5. protonate ligand | `obabel` + RDKit (`protonate_ligand.py`) | Maestro / Epik |
| 6. build prepi      | `antechamber`         | (same) |
| 7. parameterize     | `parmchk2`            | (same) |
| 8. template coords  | Python (column-fixed) | shell `gen-templ-coords.sh` |
| 9. detect rotamers  | RDKit (`detect_rotamers.py`) | Maestro dihedral picker |
| 10. build YAML      | Python                | hand-written |
| 11. verify          | `osprey affinity --verify-design` | (same) |

Maestro is **not** required at any step.

## Install

```bash
conda env create -f environment.yml
conda activate resistor-prep
```

Requirements (all from conda-forge / bioconda):
`ambertools>=23` · `openbabel>=3.1` · `rdkit>=2023.09` · `python>=3.11`.

## Usage

```bash
python prep_complex.py \
    --pdb_id 4qta --ligand_resname 38Z \
    --tag erk2_sch7 --out_dir data/erk2_sch7 \
    --protein_flex A:105:GLN --protein_flex A:106:LEU:ALA,VAL \
    --scan_res     A:147:MET
```

Flags:
- `--protein_flex CHAIN:RES:AA[:MUT,...]` (repeatable) — flexible / mutable
  residue. With `:MUT,...` makes it mutable to those AAs.
- `--ligand_flex CHAIN:RES:AA` (repeatable) — flexible ligand residue.
- `--scan_res CHAIN:RES:AA` (repeatable) — hotspot for
  `osprey affinity --do-scan`.
- `--skip_prep` — reuse cached antechamber/parmchk2 outputs.
- `--skip_verify` — skip the final `osprey affinity --verify-design` check.

## Files

| File | Role |
|------|------|
| `prep_complex.py`     | end-to-end orchestrator |
| `protonate_ligand.py` | OpenBabel + RDKit ligand protonation, net-charge sidecar |
| `detect_rotamers.py`  | RDKit rotatable-bond → OSPREY `.rot` writer |
| `p4a-undo.py`         | restores canonical chain/resnum after `pdb4amber` |
| `environment.yml`     | minimal conda env |

`p4a-undo.py` is bundled here so the pipeline is fully self-contained;
historical copies live in the RESISTOR Zenodo deposit.
