import os
from Bio.PDB import PDBParser, PDBIO

# this handles PDB prep including atomic labelling, bonds, protonation, minimization, and inversions
# this supposes a D-peptide in complex with an L-target, with the L-target listed 1st in the PDB (index 0)

# input_directory = where PDBs (from Scaffold_Generator) are located
# output_directory = where to save the processed matches


def DL_preprocess(input_directory: str, output_directory: str):

    for f in os.listdir(input_directory):
        input_file = os.path.join(input_directory, f)
        output_file = os.path.join(output_directory, f)

        # design preprocess step one: precheck PRO + change atomic labelling
        # commonly, PDB atom labels for carbon don't match ambertools templates. Let's change this manually to avoid
        # any issues with adding missing atoms

        # create the required BioPython modules so we can work with the PDB
        parser = PDBParser(PERMISSIVE=1)
        structure_id = 'complex'
        filename = input_file
        structure = parser.get_structure(structure_id, filename)

        # get 1st PDB in the ensemble
        model = structure[0]

        # OSPREY anchor coords prevent designs with proteins that contain N/C-term prolines
        # if the PDB has this N/C-term residue, skip the match file (we can't design with it, for now)
        need_skip = False
        for chain in model:
            clength = len(chain)
            counter = 0
            for res in chain:
                counter += 1
                if counter == 1 or counter == clength:
                    # check the term
                    atype = res.get_resname()
                    if atype == "PRO":
                        print("%s has an N or C-term proline, and will not be processed." % filename)
                        need_skip = True
                        break
        if need_skip:
            continue

        # MOST COMMON: change CD -> CD1 labelling
        # affected residues: ILE, LEU, PHE, TRP, TYR
        # loop over both chains + change CD labelling
        for chain in model:
            for res in chain:
                if res.get_resname() in ("ILE", "LEU", "PHE", "TRP", "TYR"):
                    for a in res:
                        if a.fullname == ' CD ':
                            a.fullname = ' CD1'
                            print("ATOMIC LABEL: changed residue " + str(res) + " to" + a.fullname)

        # design preprocess step two: flip the D-peptide to L-space, so we can work in ambertools w/ a homochiral system
        # get the chain names (note: assumes L-target is listed 1st in the PDB)
        ids = []
        for chain in model:
            ids.append(chain.id)
        peptide_id = ids[1]

        # get peptide chain
        peptide = model[peptide_id]

        # flip the D-peptide to an L-peptide
        # this separates the chains in space, so we don't need to worry about chiral-induced clashes when prepping
        for r in peptide:
            for a in r:
                a.coord[2] = a.coord[2] * -1

        # save the L:L complex with a different suffix
        io = PDBIO()
        io.set_structure(structure)
        L_filename = output_file[:-12] + "-L-complex.pdb"
        io.save(L_filename)

        # design preprocess step 4: prepare the L-ligand and L-target
        # ambertools has trouble with D-peptide protonation, so we'll prepare everything in L-space, then flip back

        # load in the L-L complex
        import osprey
        osprey.start()
        import osprey.prep
        pdb = osprey.prep.loadPDB(open(L_filename, 'r').read())
        target = pdb[0]
        peptide = pdb[1]
        mols = [target, peptide]

        with osprey.prep.LocalService():
            # remove all but the first duplicated atom from each group
            for mol in mols:
                for group in osprey.prep.duplicateAtoms(mol):
                    for atomi in range(1, len(group.getAtoms())):
                        group.remove(atomi)
                        print('removed duplicate atom %s' % group)

            # add missing heavy atoms
            for mol in mols:
                for missing_atom in osprey.prep.inferMissingAtoms(mol):
                    missing_atom.add()
                    print('added missing atom: %s' % missing_atom)

            # strip and add Hs
            # We don't know if the PDB depositor used reliable protonation software, so let's remove and use our own
            # it's also very possible no H are present, but it's a good idea to check
            for mol in mols:
                osprey.prep.deprotonate(mol)
                protonated_atoms = osprey.prep.inferProtonation(mol)
                for protonated_atom in protonated_atoms:
                    protonated_atom.add()
                print('added %d hydrogens to %s' % (len(protonated_atoms), mol))

            # save the L-target PDB
            L_target_file = output_file[:-12] + "-L-target.pdb"
            open(L_target_file, 'w').write(osprey.prep.savePDB(target))
            print('saved prepared PDB to %s' % L_target_file)

            # save the L-peptide PDB
            L_peptide_file = output_file[:-12] + "-L-peptide.pdb"
            open(L_peptide_file, 'w').write(osprey.prep.savePDB(peptide))
            print('saved prepared PDB to %s' % L_peptide_file)

        # get the L-target chain
        # OSPREY prep adds "END" to the PDB, which will throw a biopython warning. This can be safely ignored.
        target_structure = parser.get_structure('L_target', L_target_file)

        # flip the L-peptide back to D-space + add to the L-target
        peptide_structure = parser.get_structure('D_peptide', L_peptide_file)
        peptide_model = peptide_structure[0]
        peptide_chain = peptide_model[peptide_id]
        for r in peptide_chain:
            for a in r:
                a.coord[2] = a.coord[2] * -1
        peptide_chain.detach_parent()
        target_structure[0].add(peptide_chain)

        # save the D-L complex
        io = PDBIO()
        io.set_structure(target_structure)
        L_filename = output_file[:-12] + "-D-L-complex.pdb"
        io.save(L_filename)

        # design preprocess step 5: minimize the D-ligand wrt the L-target
        # SANDER requires the chemical context of the D-L context for minimization, so we must pass the complex

        # save the minimized complex

# example call
# DL_preprocess("scaffolds", "prepared-PDBs")