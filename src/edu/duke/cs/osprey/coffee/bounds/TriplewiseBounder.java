package edu.duke.cs.osprey.coffee.bounds;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.zmat.ClusterZMatrix;
import edu.duke.cs.osprey.tools.BigExp;


public class TriplewiseBounder implements Bounder {

	public final ClusterZMatrix zmat;

	private final PairwiseBounder pairwiseBounder;

	public TriplewiseBounder(ClusterZMatrix zmat) {

		this.zmat = zmat;

		pairwiseBounder = new PairwiseBounder(zmat);
	}

	@Override
	public BigExp g(ConfIndex index) {

		BigExp z = new BigExp(1.0);

		// start with the static-static energy
		z.mult(zmat.staticStatic());

		// multiply all the singles and pairs
		for (int i1=0; i1<index.numDefined; i1++) {
			int posi1 = index.definedPos[i1];
			int confi1 = index.definedRCs[i1];

			z.mult(zmat.single(posi1, confi1));

			for (int i2=0; i2<i1; i2++) {
				int posi2 = index.definedPos[i2];
				int confi2 = index.definedRCs[i2];

				z.mult(zmat.pair(posi1, confi1, posi2, confi2));
			}
		}

		return z;
	}

	@Override
	public BigExp h(ConfIndex index, RCs rcs) {
		return pairwiseBounder.h(index, rcs);
	}
}
