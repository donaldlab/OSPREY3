package edu.duke.cs.osprey.coffee.bounds;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.nodedb.NodeTree;
import edu.duke.cs.osprey.coffee.zmat.ClusterZMatrix;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools.Optimizer;


public class TriplewiseBounder implements Bounder {

	public final ClusterZMatrix zmat;

	private final PairwiseBounder pairwiseBounder;
	private final Optimizer opt;

	public TriplewiseBounder(ClusterZMatrix zmat, Optimizer opt) {

		this.zmat = zmat;
		this.opt = opt;

		pairwiseBounder = new PairwiseBounder(zmat, opt);
	}
	public TriplewiseBounder(ClusterZMatrix zmat){
		this(zmat, Optimizer.Maximize);
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

				for (int i3=0; i3<i2; i3++) {
					int posi3 = index.definedPos[i3];
					int confi3 = index.definedRCs[i3];

					// TODO: accelerate this by pre-calculating all the posi3,confi3 values for p1c1p2c2 pairs
					//  then intersect against the defined positions in the conf
					// TODO: could also optimize the index calculation for the triples
					var triple = zmat.triple(posi1, confi1, posi2, confi2, posi3, confi3);
					if (triple != null) {
						z.mult(triple);
					}
				}
			}
		}

		return z;
	}

	@Override
	public BigExp h(ConfIndex index, NodeTree tree) {
		return pairwiseBounder.h(index, tree);
	}
}
