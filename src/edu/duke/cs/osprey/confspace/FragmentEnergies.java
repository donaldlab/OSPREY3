package edu.duke.cs.osprey.confspace;

public interface FragmentEnergies {

	default double getEnergy(SimpleConfSpace.Position pos, SimpleConfSpace.ResidueConf rc) {
		return getEnergy(pos.index, rc.index);
	}

	default double getEnergy(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf rc2) {
		return getEnergy(pos1.index, rc1.index, pos2.index, rc2.index);
	}

	double getEnergy(int pos, int rc);
	double getEnergy(int pos1, int rc1, int pos2, int rc2);
}
