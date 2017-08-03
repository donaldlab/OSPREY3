
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb("2RL0.min.reduce.pdb")

# make sure all strands share the same template library (including wild-type rotamers)
templateLib = osprey.TemplateLibrary(ffparams.forcefld, moleculesForWildTypeRotamers=[mol])

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=[648, 654])
protein.flexibility[649].setLibraryRotamers(osprey.WILD_TYPE, 'TYR', 'ALA', 'VAL', 'ILE', 'LEU').addWildTypeRotamers().setContinuous()
protein.flexibility[650].setLibraryRotamers().addWildTypeRotamers().setContinuous()
protein.flexibility[651].setLibraryRotamers().addWildTypeRotamers().setContinuous()
protein.flexibility[654].setLibraryRotamers().addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=[155, 194])
ligand.flexibility[156].setLibraryRotamers().addWildTypeRotamers().setContinuous()
ligand.flexibility[172].setLibraryRotamers().addWildTypeRotamers().setContinuous()
ligand.flexibility[192].setLibraryRotamers().addWildTypeRotamers().setContinuous()
ligand.flexibility[193].setLibraryRotamers().addWildTypeRotamers().setContinuous()

# make the conf space for the protein
proteinConfSpace = osprey.ConfSpace(protein)

# make the conf space for the ligand
ligandConfSpace = osprey.ConfSpace(ligand)

# make the conf space for the protein+ligand complex
complexConfSpace = osprey.ConfSpace([protein, ligand])

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism)

# how should we define energies of conformations?
def confEcalcFactory(confSpace, ecalc):
	eref = osprey.ReferenceEnergies(confSpace, ecalc)
	return osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

# how should confs be ordered and searched?
def astarFactory(emat, rcs):
	return osprey.AStarTraditional(emat, rcs)
	# or
	# return osprey.AStarMPLP(emat, rcs, numIterations=5)

# run K*
kstar = osprey.KStar(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace,
	ecalc,
	confEcalcFactory,
	astarFactory,
	epsilon=0.99, # you proabably want something more precise in your real designs
	writeSequencesToConsole=True,
	writeSequencesToFile='kstar.results.tsv'
)
scoredSequences = kstar.run()

# use results
wildtype = osprey.Sequence.makeWildType(complexConfSpace)
for sequence in scoredSequences:
	print(sequence.toString(wildtype))
