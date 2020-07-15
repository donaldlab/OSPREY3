package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.bounds.Bounder;
import edu.duke.cs.osprey.coffee.bounds.PairwiseBounder;
import edu.duke.cs.osprey.coffee.bounds.TriplewiseBounder;
import edu.duke.cs.osprey.coffee.zmat.ClusterZMatrix;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.BigExp;

import java.util.Arrays;
import java.util.stream.IntStream;


public class StateInfo {

	public final Coffee.StateConfig config;
	public final BoltzmannCalculator bcalc;

	public final ClusterZMatrix zmat;
	public final EnergyBoundStats energyBoundStats;
	public final int[] posPermutation;

	private final int[][] typesByConfByPos;
	private final int[] numTypesByPos;

	private Bounder bounder = null;

	public StateInfo(Coffee.StateConfig config, BoltzmannCalculator bcalc) {

		this.config = config;
		this.bcalc = bcalc;

		zmat = new ClusterZMatrix(config.confSpace, config.posInterGen, bcalc);
		energyBoundStats = new EnergyBoundStats();

		// sort positions so multi-sequence layers are first
		posPermutation = IntStream.range(0, config.confSpace.numPos())
			.boxed()
			.sorted((a, b) -> {

				boolean amut = config.confSpace.hasMutations(a);
				boolean bmut = config.confSpace.hasMutations(b);

				// prefer mutable positions first
				if (amut && !bmut) {
					return -1;
				} else if (!amut && bmut) {
					return +1;
				}

				// then, sort by ordering heuristic
				int order = Double.compare(calcOrderHeuristic(a), calcOrderHeuristic(b));
				if (order != 0) {
					// negate to sort in (weakly) descending order
					return -order;
				}

				// otherwise, sort by index
				return a - b;
			})
			.mapToInt(i -> i)
			.toArray();

		// calculate all the conf types by conf and pos
		typesByConfByPos = new int[config.confSpace.numPos()][];
		numTypesByPos = new int[config.confSpace.numPos()];
		for (int posi=0; posi<config.confSpace.numPos(); posi++) {

			typesByConfByPos[posi] = new int[config.confSpace.numConf(posi)];
			SeqSpace.Position seqPos = config.confSpace.seqSpace.getPosition(config.confSpace.name(posi));
			if (seqPos != null) {
				numTypesByPos[posi] = seqPos.resTypes.size();
				for (int confi=0; confi<config.confSpace.numConf(posi); confi++) {
					SeqSpace.ResType rt = seqPos.getResTypeOrThrow(config.confSpace.confType(posi, confi));
					typesByConfByPos[posi][confi] = rt.index;
				}
			} else {
				numTypesByPos[posi] = 0;
				Arrays.fill(typesByConfByPos[posi], Sequence.Unassigned);
			}
		}
	}

	private double calcOrderHeuristic(int posi) {

		// TODO: find a heuristic that works well here
		return 0.0;
	}

	public void initBounder() {

		// choose the best bounder for this state
		// TODO: find better bounders
		if (zmat.hasTriples()) {
			bounder = new TriplewiseBounder(zmat);
		} else {
			bounder = new PairwiseBounder(zmat);
		}
	}

	public ConfIndex makeConfIndex() {
		ConfIndex index = new ConfIndex(config.confSpace.numPos());
		index.updateUndefined();
		return index;
	}

	public BigExp leavesBySequenceUpper(ConfIndex index, RCs rcs) {

		BigExp count = new BigExp(1.0);

		for (int i=0; i<index.numUndefined; i++) {
			int posi = index.undefinedPos[i];

			int maxCount = 0;

			int numTypes = numTypesByPos[posi];
			if (numTypes > 0) {

				// count the confs by type
				int[] counts = new int[numTypes];
				for (int confi : rcs.get(posi)) {
					int confCount = ++counts[typesByConfByPos[posi][confi]];
					maxCount = Math.max(maxCount, confCount);
				}

			} else {

				// all confs are the same type
				maxCount = rcs.getNum(posi);
			}

			count.mult(maxCount);
		}

		return count;
	}

	public BigExp zSumUpper(ConfIndex index, RCs rcs) {
		BigExp out = zPathHead(index);
		out.mult(zPathTailUpper(index, rcs));
		out.mult(leavesBySequenceUpper(index, rcs));
		return out;
	}

	public BigExp zPathUpper(ConfIndex index, RCs rcs) {
		BigExp out = zPathHead(index);
		out.mult(zPathTailUpper(index, rcs));
		return out;
	}

	public BigExp zPathHead(ConfIndex index) {
		return bounder.g(index);
	}

	public BigExp zPathTailUpper(ConfIndex index, RCs rcs) {
		return bounder.h(index, rcs);
	}
}
