package edu.duke.cs.osprey.multistatekstar;

import java.util.Comparator;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class ConfComparator implements Comparator<ScoredConf> {

	@Override
	public int compare(ScoredConf o1, ScoredConf o2) {
		double e1 = o1 instanceof EnergiedConf ? ((EnergiedConf)o1).getEnergy() : o1.getScore();
		double e2 = o2 instanceof EnergiedConf ? ((EnergiedConf)o2).getEnergy() : o2.getScore();
		return e1 <= e2 ? 1 : -1;
	}

}
