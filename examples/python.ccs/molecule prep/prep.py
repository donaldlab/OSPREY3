
import osprey
osprey.start()

# import the prep module after starting Osprey
import osprey.prep


# Before we get started:
# Preparing a molecule for Osprey can be a very error-prone process.
# Automated scripts, by their very nature, don't provide any feedback
# to the user between processing steps, so these errors are very likely
# to go unnoticed. For that reason, we strongly recommend that molecule
# preparations be done interactively using the GUI rather than an automated
# script. Nevertheless, if you already know where the errors in the prep
# process will happen, you can automate the prep using a script like this
# one, after customizing it to your specific situation.


# let's load a PDB file
# this one is a crystal structure of the Bovine PTPase complexed with HEPES
pdb_path = '1dg9.pdb'
pdb = osprey.prep.loadPDB(open(pdb_path, 'r').read())

# what molecules did we get?
print('Loaded %d molecules:' % len(pdb))
for mol in pdb:
    print('\t%s: %s' % (mol, osprey.prep.molTypes(mol)))

# looks like the PDB file has a protein chain and a small molecule
# the protein chain must be PTPase, and the small molecule must be HEPES
# (ignore the solvent molecules, if any)
ptpase = pdb[0]
hepes = pdb[1]
mols = [ptpase, hepes]

# let's give our molecules better names
ptpase.setName('Bovine PTPase')
hepes.setName('HEPES')

# start the local service that calls AmberTools for us
# NOTE: this will only work on Linux machines
with osprey.prep.LocalService():

    # Molecule Preparation Step 1: remove duplicate atoms
    # Duplicate atoms don't usually appear in files from the PDB,
    # but these errors can happen sometimes in modified PDB files.
    for mol in mols:
        # remove all but the first duplicated atom from each group
        for group in osprey.prep.duplicateAtoms(mol):
            for atomi in range(1, len(group.getAtoms())):
                group.remove(atomi)
                print('removed duplicate atom %s' % group)

    # Molecule Preparation Step 2: add missing heavy atoms
    # Somtimes atom positions are not well-resolved in the electron density,
    # or protein chain end-caps are missing.
    # But we still want to include these atoms in the molecular models to
    # be able to infer bonds correctly in the next step.
    for mol in mols:
        for missing_atom in osprey.prep.inferMissingAtoms(mol):
            missing_atom.add()
            print('added missing atom: %s' % missing_atom)

    # Molecule Preparation Step 3: add bonds
    # PDB files contain no explicit information about bonds, so we have to
    # infer where they might be based on the atoms we can see.
    for mol in mols:
        bonds = osprey.prep.inferBonds(mol)
        for bond in bonds:
            mol.getBonds().add(bond)
        print('added %d bonds to %s' % (len(bonds), mol))

    # Molecule Preparation Step 4: add hydrogens
    # Most PDB files based on X-ray crystallography don't describe Hydrogen atoms,
    # usually because light atoms aren't well-resolved in the electron density.
    # Generally, automatic protonation inference works well on proteins,
    # so let's just protonate PTPase for now, we'll protonate HEPES later
    for mol in [ptpase]:
        protonated_atoms = osprey.prep.inferProtonation(mol)
        for protonated_atom in protonated_atoms:
            protonated_atom.add()
        print('added %d hydrogens to %s' % (len(protonated_atoms), mol))

    # Inferred protonation is often *very* wrong for small molecules,
    # so we'll specify the protonation state for HEPES manually here
    for atom in ['O3S', 'O8']:
        osprey.prep.protonate(hepes, hepes.getAtoms().findOrThrow(atom), 1, osprey.prep.Hybridization.Sp3)
    for atom in ['C10', 'C9', 'C6', 'C5', 'C3', 'C2', 'C7', 'C8']:
        osprey.prep.protonate(hepes, hepes.getAtoms().findOrThrow(atom), 2, osprey.prep.Hybridization.Sp3)

    # Molecule Preparation Step 5: formal net charges for small molecules
    # Formal net charges are used by Amber to calculate partial charges
    # which parameterize the electrostatic components of the forcefield for small molecules.
    # Let's also add a formal net charge for HEPES, but we don't need to add one for PTPase.
    hepes.setNetCharge(osprey.jvm.boxInt(0))

    # Molecule Preparation Step 6: minimize structures
    # Osprey picks sequences by looking for low-energy conformations.
    # Let's make sure our starting structure is low-energy
    # by doing an all-atom minimization against the Amber forcefield
    # but restrain the heavy atoms so they don't move around too much
    # print('minimizing ...')
    # def heavy_atoms(mol):
    #     return [atom for atom in mol.getAtoms() if atom.getElement().getSymbol() != 'H']
    # osprey.prep.minimize(
    #     [osprey.prep.minimizerInfo(mol, heavy_atoms(mol)) for mol in mols],
    #     numSteps = 100
    # )
    # print('minimization complete!')

    # Moleclue Preparation Step 7: save the results
    path = '1dg9.omol'
    open(path, 'w').write(osprey.prep.saveOMOL(mols))
    print('saved prepared molecules to %s' % path)


print('Moleclue preparation complete!')
