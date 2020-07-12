package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
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

				// also multiply triples, if any
				if (zmat.hasTriples()) {
					for (int i3=0; i3<i2; i3++) {
						int posi3 = index.definedPos[i3];
						int confi3 = index.definedRCs[i3];

						var triple = zmat.triple(posi1, confi1, posi2, confi2, posi3, confi3);
						if (triple != null) {
							z.mult(triple);
						}
					}
				}
			}
		}

		return z;
	}

	public BigExp zPathTailUpper(ConfIndex index, RCs rcs) {

		// this is the usual A* heuristic

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
			for (int confi1 : rcs.get(posi1)) {

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

					// optimize over possible assignments to pos2
					zrc1.mult(zmat.pairUpper(posi1, confi1, posi2));
				}

				zpos1.max(zrc1);
			}

			assert (zpos1.isFinite());
			z.mult(zpos1);
		}

		return z;
	}
}
