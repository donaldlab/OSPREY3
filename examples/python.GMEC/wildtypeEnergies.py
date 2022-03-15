
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility['A2'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA', 'GLY').addWildTypeRotamers()

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# how should we compute energies of molecules?
ecalc = osprey.EnergyCalculator(confSpace, ffparams)

# calculate the full moleclue energy of the wild type
pmol = osprey.c.confspace.ParametricMolecule(strand.mol)
inters = osprey.c.energy.ResidueInteractions()
inters.addComplete(pmol.mol.residues)
energy = ecalc.calcEnergy(pmol, inters).energy
print('wild-type molecule energy: %.2f' % energy)

# how should we define energies of conformations?
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

# get the conformation of the wild-type
# no simple way to get it other than to filter down to the wild-type residue conformations for each design position
def wtIndex(pos):
	for rc in pos.resConfs:
		if rc.type == osprey.c.confspace.SimpleConfSpace.ResidueConf.Type.WildType:
			return rc.index
	raise Exception('no wild-type conformation found at design position: %s' % pos)
wildTypeConf = [wtIndex(pos) for pos in confSpace.positions]
print('wild-type conformation indices: %s' % wildTypeConf)

# some Osprey functions want the input in a slightly different form
wildTypeTuple = osprey.c.confspace.RCTuple(wildTypeConf)

# calculate the wild-type energy using the conformation system
energy = confEcalc.calcEnergy(wildTypeTuple).energy
print('wild-type conformation energy: %.2f' % energy)

# calculate the wild-type energy using the conformation system,
# but override to use all residue interations in the molecule
inters = osprey.c.energy.EnergyPartition.makeAll(confSpace, None, False, wildTypeTuple)
energy = confEcalc.calcEnergy(wildTypeTuple, inters).energy
print('wild-type total conformation energy: %.2f' % energy)

# calculate reference energies for the conformation space
eref = osprey.ReferenceEnergies(confSpace, ecalc)

# update the conformation energy calculator to use the reference energies
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

# calculate the wild-type energy using the conformation system, with reference energies
energy = confEcalc.calcEnergy(wildTypeTuple).energy
print('wild-type conformation energy with reference energies: %.2f' % energy)

# calculate the wild-type energy using the conformation system, with referenec energies,
# but override to use all residue interations in the molecule
inters = osprey.c.energy.EnergyPartition.makeAll(confSpace, eref, False, wildTypeTuple)
energy = confEcalc.calcEnergy(wildTypeTuple, inters).energy
print('wild-type total conformation energy with reference energies: %.2f' % energy)
