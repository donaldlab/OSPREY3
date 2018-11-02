
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('../../test-resources/4npd.A_DomainA_noH_trim_his_clean.min.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=['A1', 'A37'])
protein.flexibility['A8'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['A38', 'A58'])
ligand.flexibility['A41'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A42'].setLibraryRotamers(osprey.WILD_TYPE,'THR','LYS').addWildTypeRotamers().setContinuous()
ligand.flexibility['A45'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A46'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# make the conf space for the protein+ligand complex
complexConfSpace = osprey.ConfSpace([protein, ligand])

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism)

# configure PAStE
paste = osprey.Paste(
	complexConfSpace,
	epsilon=0.99, # you proabably want something more precise in your real designs
	useWindowCriterion=True,
	writeSequencesToConsole=True,
	writeSequencesToFile='paste.results.tsv',
	mutFile='mut.txt'
)

# configure PAStE inputs for the conf space

# how should we define energies of conformations?
eref = osprey.ReferenceEnergies(paste.protein.confSpace, ecalc)
paste.protein.confEcalc = osprey.ConfEnergyCalculator(paste.protein.confSpace, ecalc, referenceEnergies=eref)

# compute the energy matrix
emat = osprey.EnergyMatrix(paste.protein.confEcalc, cacheFile='emat.%s.dat' % paste.protein.id)

# how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
def makeAStar(rcs, emat=emat):
	return osprey.AStarTraditional(emat, rcs, showProgress=False)
paste.protein.confSearchFactory = osprey.Paste.ConfSearchFactory(makeAStar)

# run PAStE
paste.run()
