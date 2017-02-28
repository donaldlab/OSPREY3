package edu.duke.cs.osprey.multistatekstar;

import java.util.Comparator;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;

public class EnergiedConfComparator implements Comparator<EnergiedConf> {

	@Override
	public int compare(EnergiedConf o1, EnergiedConf o2) {
		return o1.getEnergy() <= o2.getEnergy() ? 1 : -1;
	}

}
