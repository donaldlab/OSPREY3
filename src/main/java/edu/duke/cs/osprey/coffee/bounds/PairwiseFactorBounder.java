package edu.duke.cs.osprey.coffee.bounds;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.coffee.nodedb.NodeTree;
import edu.duke.cs.osprey.coffee.zmat.ClusterZMatrix;
import edu.duke.cs.osprey.tools.BigExp;

public class PairwiseFactorBounder extends PairwiseBounder{
    public PairwiseFactorBounder(ClusterZMatrix zmat) {
        super(zmat);
    }

    @Override
    public BigExp h(ConfIndex index, NodeTree tree) {
        // New bound based on factoring partition function bounds

        // TODO: applying higher-order corrections here *may* be useful.
        // they're quite slow to multiply in, but we no longer multiply
        // by the number of nodes, so the looseness of zPathTailUpper is the only factor.

        BigExp z = new BigExp(1.0);

        // for each undefined position
        for (int i1=0; i1<index.numUndefined; i1++) {
            int posi1 = index.undefinedPos[i1];

            // sum over possible assignments to pos1
            BigExp zpos1 = new BigExp(0.0);
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

                zpos1.add(zrc1);
            }

            assert (zpos1.isFinite());
            z.mult(zpos1);
        }

        return z;
    }
}
