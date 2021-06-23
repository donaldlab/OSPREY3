package edu.duke.cs.osprey.coffee.bounds;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.coffee.nodedb.NodeTree;
import edu.duke.cs.osprey.coffee.zmat.ClusterZMatrix;
import edu.duke.cs.osprey.tools.BigExp;


public class PairwiseBounder implements Bounder {

	public final ClusterZMatrix zmat;

	protected final BigExp[][] optimizationCache;

	public PairwiseBounder(ClusterZMatrix zmat) {

		this.zmat = zmat;

		// pre-populate the optimization cache
		optimizationCache = new BigExp[zmat.numPairs()][];
		for (int posi1=0; posi1<zmat.confSpace.numPos(); posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {

				int pairi = zmat.pairIndex(posi1, posi2);
				optimizationCache[pairi] = new BigExp[zmat.confSpace.numConf(posi1)];

				for (int confi1=0; confi1<zmat.confSpace.numConf(posi1); confi1++) {

					BigExp max = zmat.pair(posi1, confi1, posi2, 0);
					for (int confi2=1; confi2<zmat.confSpace.numConf(posi2); confi2++) {
						BigExp z = zmat.pair(posi1, confi1, posi2, confi2);
						if (z.greaterThan(max)) {
							max = z;
						}
					}
					optimizationCache[pairi][confi1] = max;
				}
			}
		}
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
	public BigExp h(ConfIndex index, NodeTree tree) {

		// this is the classic A* heuristic

		// NOTE: applying higher-order corrections here isn't terribly useful
		// they're quite slow to multiply in, and don't improve zSumUpper that much
		// of course, they help get a better zPathTailUpper, but most of the zSumUpper looseness
		// comes from multiplying by the number of nodes rather than the looseness of zPathTailUpper

		BigExp z = new BigExp(1.0);

		// for each undefined position
		for (int i1=0; i1<index.numUndefined; i1++) {
			int posi1 = index.undefinedPos[i1];

			// optimize over possible assignments to pos1
			BigExp zpos1 = new BigExp(Double.NEGATIVE_INFINITY);
			for (int confi1 : tree.rcs.get(posi1)) {

				BigExp zrc1 = new BigExp(zmat.single(posi1, confi1));

				// interactions with defined residues
				for (int i2=0; i2<index.numDefined; i2++) {
					int posi2 = index.definedPos[i2];
					int confi2 = index.definedRCs[i2];

					zrc1.mult(zmat.pair(posi1, confi1, posi2, confi2));
				}

				// interactions with undefined residues
				for (int i2=0; i2<i1; i2++) {
					int posi2 = index.undefinedPos[i2];
					int pairi = zmat.pairIndex(posi1, posi2);

					// optimize over possible assignments to pos2
					zrc1.mult(optimizationCache[pairi][confi1]);
				}

				zpos1.max(zrc1);
			}

			assert (zpos1.isFinite());
			z.mult(zpos1);
		}

		return z;
	}
}
