
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# choose a molecule
mol = osprey.readPdb("2RL0.min.reduce.pdb")

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld, moleculesForWildTypeRotamers=[mol])

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=[648, 654])
protein.flexibility[649].setLibraryRotamers().addWildTypeRotamers().setContinuous()
protein.flexibility[650].setLibraryRotamers().addWildTypeRotamers().setContinuous()
protein.flexibility[651].setLibraryRotamers().addWildTypeRotamers().setContinuous()
protein.flexibility[654].setLibraryRotamers().addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=[155, 194])
ligand.flexibility[156].setLibraryRotamers().addWildTypeRotamers().setContinuous()
ligand.flexibility[172].setLibraryRotamers().addWildTypeRotamers().setContinuous()
ligand.flexibility[192].setLibraryRotamers().addWildTypeRotamers().setContinuous()
ligand.flexibility[193].setLibraryRotamers().addWildTypeRotamers().setContinuous()

# make the conf spaces
proteinConfSpace = osprey.ConfSpace(protein)
ligandConfSpace = osprey.ConfSpace(ligand)
complexConfSpace = osprey.ConfSpace([protein, ligand])

# how should we compute energies of molecules?
parallelism = osprey.Parallelism(cpuCores=4)
def ecalcFactory(confSpace):
	return osprey.EnergyCalculator(confSpace, ffparams)

# how should we define energies of conformations?
def confEcalcFactory(confSpace, ecalc):
	eref = osprey.ReferenceEnergies(confSpace, ecalc)
	return osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

# how should confs be ordered and searched?
def astarFactory(emat, pmat):
	return osprey.AStarTraditional(emat, pmat)
	# or
	# return osprey.AStarMPLP(emat, pmat, numIterations=5)

# run K*
osprey.KStar(proteinConfSpace, ligandConfSpace, complexConfSpace, ecalcFactory, confEcalcFactory, astarFactory)
