package edu.duke.cs.osprey.energy.compiled;


import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;

import java.util.Collections;
import java.util.List;


/**
 * A translation layer that allows the new conf energy calculators
 * to be somewhat compatible with the old ones.
 *
 * Ideally, all the consumers of the old ConfEnergyCalculator should be updated
 * to use the new ConfEnergyCalculator instead, and the new ConfSpace class.
 * Using the newer code will lead to better performance, but we all understand
 * that change is hard. This temporary (hah!) class should keep all the old code
 * happy for a time, but try not to write any new code that depends on this.
 */
@Deprecated
public class ConfEnergyCalculatorAdapter extends edu.duke.cs.osprey.energy.ConfEnergyCalculator {

	public static class Builder {

		public final ConfEnergyCalculator confEcalc;

		/** Defines what position interactions should be used for conformations and conformation fragments */
		private PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;

		/** Reference energies used to design against the unfolded state, useful for GMEC designs */
		private SimpleReferenceEnergies eref = null;

		private boolean minimize = true;

		/** True to include the static-static energies in conformation energies */
		private boolean includeStaticStatic = false;

		public Builder(ConfEnergyCalculator confEcalc) {
			this.confEcalc = confEcalc;
		}

		public Builder setPosInterDist(PosInterDist val) {
			posInterDist = val;
			return this;
		}

		public Builder setReferenceEnergies(SimpleReferenceEnergies val) {
			eref = val;
			return this;
		}

		public Builder setMinimize(boolean val) {
			minimize = val;
			return this;
		}

		public Builder setIncludeStaticStatic(boolean val) {
			includeStaticStatic = val;
			return this;
		}

		public ConfEnergyCalculatorAdapter build() {
			return new ConfEnergyCalculatorAdapter(
				confEcalc,
				new PosInterGen(posInterDist, eref),
				minimize,
				includeStaticStatic
			);
		}
	}

	public final ConfEnergyCalculator confEcalc;
	public final PosInterGen posInterGen;
	public final boolean minimize;
	public final boolean includeStaticStatic;

	// TODO: support approximator matrices?

	private ConfEnergyCalculatorAdapter(ConfEnergyCalculator confEcalc, PosInterGen posInterGen, boolean minimize, boolean includeStaticStatic) {
		super(confEcalc.tasks());

		this.confEcalc = confEcalc;
		this.posInterGen = posInterGen;
		this.minimize = minimize;
		this.includeStaticStatic = includeStaticStatic;
	}

	@Override
	public ConfSpaceIteration confSpaceIteration() {
		return confEcalc.confSpace();
	}

	@Override
	public ResidueInteractions makeFragInters(RCTuple frag) {
		// only used internally, so we don't need to implement it
		throw new UnsupportedOperationException();
	}

	@Override
	public ResidueInteractions makeSingleInters(int pos, int rc) {
		// only used by SimplerEnergyMatrixCalculator
		// so we only need to return a ResidueInteractions instance that has
		// a size useful enough to calibrate the progress bar
		ResidueInteractions inters = new ResidueInteractions();
		inters.addSingle("a"); // these resnums will never be seen
		return inters;
	}

	@Override
	public ResidueInteractions makePairInters(int pos1, int rc1, int pos2, int rc2) {
		// only used by SimplerEnergyMatrixCalculator
		// so we only need to return a ResidueInteractions instance that has
		// a size useful enough to calibrate the progress bar
		ResidueInteractions inters = new ResidueInteractions();
		inters.addPair("a", "b"); // these resnums will never be seen
		return inters;
	}

	@Override
	public ResidueInteractions makeTupleInters(RCTuple tuple) {
		throw new UnsupportedOperationException();
	}

	@Override
	public ResidueInteractions makeTripleCorrectionInters(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
		throw new UnsupportedOperationException();
	}

	@Override
	public ResidueInteractions makeQuadCorrectionInters(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, int pos4, int rc4) {
		throw new UnsupportedOperationException();
	}

	@Override
	public EnergyCalculator.EnergiedParametricMolecule calcSingleEnergy(int pos, int rc) {
		throw new UnsupportedOperationException();
	}

	private ParametricMolecule pmol(AssignedCoords coords) {

		// convert the coords to a pmol, without the dofs, since they can't be converted
		return new ParametricMolecule(coords.toMol(), null, null);
	}

	private EnergyCalculator.EnergiedParametricMolecule epmol(int[] conf, List<PosInter> inters) {

		if (minimize) {

			// minimize the conformation
			ConfEnergyCalculator.EnergiedCoords result = confEcalc.minimize(conf, inters);

			// wrap the energy and the dof values in an otherwise-empty epmol
			// and hope the caller doesn't want the pmol or the residue interactions
			return new EnergyCalculator.EnergiedParametricMolecule(pmol(result.coords), null, result.dofValues, result.energy);

		} else {

			ConfEnergyCalculator.EnergiedCoords result = confEcalc.calc(conf, inters);

			// wrap the energy in an otherwise-empty epmol
			// and hope the caller doesn't want the pmol or the residue interactions
			return new EnergyCalculator.EnergiedParametricMolecule(pmol(result.coords), null, result.energy);
		}
	}

	@Override
	public EnergyCalculator.EnergiedParametricMolecule calcIntraEnergy(int posi, int confi) {

		// convert the fragment to a conformation in the space
		int[] conf = confEcalc.confSpace().assign(posi, confi);

		// NOTE: don't need to add reference energies here, since this is used by the reference energy calculator

		// make the position interactions for just the intra energy at this position
		List<PosInter> inters = Collections.singletonList(new PosInter(posi, posi, 1.0, 0.0));

		return epmol(conf, inters);
	}

	@Override
	public EnergyCalculator.EnergiedParametricMolecule calcSingleEnergy(RCTuple frag) {

		// convert the fragment to a conformation in the space
		int posi = frag.pos.get(0);
		int confi = frag.RCs.get(0);
		int[] conf = confEcalc.confSpace().assign(posi, confi);

		// make the position interactions for a single tuple
		List<PosInter> inters = posInterGen.single(confEcalc.confSpace(), posi, confi);

		return epmol(conf, inters);
	}

	@Override
	public EnergyCalculator.EnergiedParametricMolecule calcPairEnergy(RCTuple frag) {

		// convert the fragment to a conformation in the space
		int posi1 = frag.pos.get(0);
		int confi1 = frag.RCs.get(0);
		int posi2 = frag.pos.get(1);
		int confi2 = frag.RCs.get(1);
		int[] conf = confEcalc.confSpace().assign(posi1, confi1, posi2, confi2);

		// make the position interactions for a single tuple
		List<PosInter> inters = posInterGen.pair(confEcalc.confSpace(), posi1, confi1, posi2, confi2);

		return epmol(conf, inters);
	}

	@Override
	public EnergyCalculator.EnergiedParametricMolecule calcTupleEnergy(RCTuple frag) {
		throw new UnsupportedOperationException();
	}

	@Override
	public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag) {

		int[] conf = confEcalc.confSpace().assign(frag);

		// make the position interactions for the whole conformation
		List<PosInter> inters;
		if (includeStaticStatic) {
			inters = posInterGen.all(confEcalc.confSpace(), conf);
		} else {
			inters = posInterGen.dynamic(confEcalc.confSpace(), conf);
		}

		return epmol(conf, inters);
	}

	@Override
	public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag, ResidueInteractions inters) {
		throw new UnsupportedOperationException();
	}
}
