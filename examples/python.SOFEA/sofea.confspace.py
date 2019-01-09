
# this script can't be run on its own, it must be included into another script via execfile()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('2RL0.min.reduce.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the design strand
designStrand = osprey.Strand(mol, templateLib=templateLib, residues=['G648', 'G654'])
designStrand.flexibility['G649'].setLibraryRotamers(osprey.WILD_TYPE, 'TYR', 'ALA', 'VAL', 'ILE', 'LEU').addWildTypeRotamers().setContinuous()
designStrand.flexibility['G650'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
designStrand.flexibility['G651'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the ligand strand
ligandStrand = osprey.Strand(mol, templateLib=templateLib, residues=['A155', 'A194'])
ligandStrand.flexibility['A172'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligandStrand.flexibility['A192'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligandStrand.flexibility['A193'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# make a multi-state conf space
confSpace = osprey.MultiStateConfSpace([
	osprey.StateMutable('design', osprey.ConfSpace(designStrand)),
	osprey.StateUnmutable('ligand', osprey.ConfSpace(ligandStrand)),
	osprey.StateMutable('complex', osprey.ConfSpace([designStrand, ligandStrand]))
])
