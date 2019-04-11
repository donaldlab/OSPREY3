/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.pruning;


import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.structure.*;
import edu.duke.cs.osprey.tools.Progress;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.linear.*;

import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Function;


/**
 * Pruning of Local Unrealistic Geometries (PLUG)
 *
 * prunes RC tuples if probe-style clashes are unavoidable
 */
public class PLUG {

	public final SimpleConfSpace confSpace;

	/** max num iterations of generalized Newton iteration to find a "boundary" point where the violation function is zero */
	public int maxNumIterations = 30;

	/** distance threshold to claim that a violation value is close enough to zero */
	public double violationThreshold = 1e-2;

	/** factor of the voxel width used to approximate the gradient of the violation function */
	public double gradientDxFactor = 1e-4;

	private final Probe probe;
	private final AtomConnectivity connectivity;

	public PLUG(SimpleConfSpace confSpace) {

		this.confSpace = confSpace;

		// load probe
		this.probe = new Probe();
		this.probe.matchTemplates(this.confSpace);

		// init atom connectivity
		this.connectivity = new AtomConnectivity.Builder()
			.addTemplates(confSpace)
			.set15HasNonBonded(false) // follows probe convention
			.build();
	}

	public void pruneSingles(PruningMatrix pmat, double tolerance) {
		pruneSingles(pmat, tolerance, new TaskExecutor());
	}

	public void pruneSingles(PruningMatrix pmat, double tolerance, TaskExecutor tasks) {

		// count unpruned singles
		AtomicLong numSingles = new AtomicLong(0);
		pmat.forEachUnprunedSingle((pos1, rc1) -> {
			numSingles.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numSingles.get());

		// try to prune each single
		pmat.forEachUnprunedSingle((pos1, rc1) -> {
			tasks.submit(
				() -> shouldPruneTuple(new RCTuple(pos1, rc1), tolerance),
				(shouldPrune) -> {
					if (shouldPrune) {
						pmat.pruneSingle(pos1, rc1);
					}
					progress.incrementProgress();
				}
			);
			return PruningMatrix.IteratorCommand.Continue;
		});

		tasks.waitForFinish();
	}

	public void prunePairs(PruningMatrix pmat, double tolerance) {
		prunePairs(pmat, tolerance, new TaskExecutor());
	}

	public void prunePairs(PruningMatrix pmat, double tolerance, TaskExecutor tasks) {

		// count unpruned pairs
		AtomicLong numPairs = new AtomicLong(0);
		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			numPairs.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numPairs.get());

		// try to prune each pair
		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			tasks.submit(
				() -> shouldPruneTuple(new RCTuple(pos1, rc1, pos2, rc2), tolerance),
				(shouldPrune) -> {
					if (shouldPrune) {
						pmat.prunePair(pos1, rc1, pos2, rc2);
					}
					progress.incrementProgress();
				}
			);
			return PruningMatrix.IteratorCommand.Continue;
		});

		tasks.waitForFinish();
	}

	public void pruneTriples(PruningMatrix pmat, double tolerance) {
		pruneTriples(pmat, tolerance, new TaskExecutor());
	}

	public void pruneTriples(PruningMatrix pmat, double tolerance, TaskExecutor tasks) {

		// count unpruned triple
		AtomicLong numTriples = new AtomicLong(0);
		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			numTriples.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numTriples.get());

		// try to prune each triple
		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			tasks.submit(
				() -> shouldPruneTuple(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3), tolerance),
				(shouldPrune) -> {
					if (shouldPrune) {
						pmat.pruneTriple(pos1, rc1, pos2, rc2, pos3, rc3);
					}
					progress.incrementProgress();
				}
			);
			return PruningMatrix.IteratorCommand.Continue;
		});

		tasks.waitForFinish();
	}

	private class Voxel {

		final ParametricMolecule pmol;

		final int numDofs;
		final double[] width;
		final double[] width2;
		final double[] min;
		final double[] max;
		final double[] center;

		Voxel (ParametricMolecule pmol) {

			this.pmol = pmol;

			numDofs = pmol.dofs.size();
			width = new double[numDofs];
			width2 = new double[numDofs];
			min = new double[numDofs];
			max = new double[numDofs];
			center = new double[numDofs];

			for (int d=0; d<numDofs; d++) {
				width[d] = pmol.dofBounds.getWidth(d);
				width2[d] = width[d]*width[d];
				min[d] = pmol.dofBounds.getMin(d);
				max[d] = pmol.dofBounds.getMax(d);
				center[d] = pmol.dofBounds.getCenter(d);
			}
		}

		DegreeOfFreedom getDof(int d) {
			return pmol.dofs.get(d);
		}

		void applyDof(int d, double val) {
			getDof(d).apply(val);
		}
	}

	public boolean shouldPruneTuple(RCTuple tuple, double tolerance) {

		// make the molecule and get all the residue interactions for the tuple
		ParametricMolecule pmol = confSpace.makeMolecule(tuple);

		ResidueInteractions inters = ResInterGen.of(confSpace)
			.addIntras(tuple)
			.addInters(tuple)
			.addShell(tuple)
			.make();

		Voxel voxel = new Voxel(pmol);

		try {

			// get linear constraints for each atom pair
			List<LinearConstraint> constraints = getLinearConstraints(voxel, inters, tolerance);

			// no constraints? don't prune
			if (constraints.isEmpty()) {
				return false;
			}

			// use an LP solver (eg simplex) to determine if the constraints allow any feasible points
			new SimplexSolver().optimize(
				new SimpleBounds(voxel.min, voxel.max),
				new LinearConstraintSet(constraints),
				// dummy function: don't really need to minimize, but can't call simplex phase 1 solver directly
				new LinearObjectiveFunction(new double[pmol.dofs.size()], 0.0)
			);

			// if we got here, simplex didn't throw an exception
			// meaning at least one feasible point exists, so don't prune this tuple
			return false;

		} catch (NoFeasibleSolutionException ex) {

			// no feasible points, prune this tuple
			return true;
		}
	}

	private class AtomVoxel {

		final Atom atom;
		final Probe.AtomInfo probeInfo;
		final List<Integer> dofIndices = new ArrayList<>();

		AtomVoxel(Atom atom, Probe probe) {
			this.atom = atom;
			this.probeInfo = probe.getAtomInfo(atom);
		}

		AtomVoxel(Atom atom, Probe probe, Voxel voxel) {
			this(atom, probe);

			// determine which dofs affect this atom position
			for (int d=0; d<voxel.numDofs; d++) {

				// skip this dof if it's not continuously flexible
				if (voxel.width[d] <= 0.0) {
					continue;
				}

				// start at the center
				voxel.applyDof(d, voxel.center[d]);
				double[] start = atom.getCoords();

				// move a little bit away from the center
				voxel.applyDof(d, voxel.center[d] + gradientDxFactor*voxel.width[d]);
				double[] stop = atom.getCoords();

				// pick the dof if the positions are different by even a little
				if (!Arrays.equals(start, stop)) {
					dofIndices.add(d);
				}
			}
		}

		public boolean hasDofs() {
			return !dofIndices.isEmpty();
		}
	}

	private class AtomPairVoxel {

		final Voxel voxel;
		final Probe.AtomPair probePair;
		final List<Integer> dofIndices = new ArrayList<>();
		final int numDofs;

		AtomPairVoxel(Voxel voxel, AtomVoxel v1, AtomVoxel v2) {

			this.voxel = voxel;

			// make the probe atom pair
			this.probePair = probe.new AtomPair(v1.atom, v2.atom, v1.probeInfo, v2.probeInfo);

			// combine the dofs
			for (int d : v1.dofIndices) {
				if (!dofIndices.contains(d)) {
					dofIndices.add(d);
				}
			}
			for (int d : v2.dofIndices) {
				if (!dofIndices.contains(d)) {
					dofIndices.add(d);
				}
			}
			numDofs = dofIndices.size();
		}

		double min(int d) {
			d = dofIndices.get(d);
			return voxel.min[d];
		}

		double max(int d) {
			d = dofIndices.get(d);
			return voxel.max[d];
		}

		double center(int d) {
			d = dofIndices.get(d);
			return voxel.center[d];
		}

		double width(int d) {
			d = dofIndices.get(d);
			return voxel.width[d];
		}

		double width2(int d) {
			d = dofIndices.get(d);
			return voxel.width2[d];
		}

		void applyDof(int d, double val) {
			d = dofIndices.get(d);
			voxel.applyDof(d, val);
		}

		double getViolation(double[] x, double tolerance) {

			// apply the dofs
			for (int d=0; d<numDofs; d++) {
				applyDof(d, x[d]);
			}

			// get the violation
			return probePair.getViolation(tolerance);
		}

		double getViolationAlong(int d, double x, double dx, double tolerance) {

			// move along one dof
			applyDof(d, x + dx);

			// get the violation
			double violation = probePair.getViolation(tolerance);

			// put the dof back
			applyDof(d, x);

			return violation;
		}

		boolean outOfRange(double[] x) {
			for (int d=0; d<numDofs; d++) {
				if (x[d] < min(d) || x[d] > max(d)) {
					return true;
				}
			}
			return false;
		}
	}

	public List<LinearConstraint> getLinearConstraints(Voxel voxel, ResidueInteractions inters, double tolerance) {

		Map<Atom,AtomVoxel> atomVoxels = new HashMap<>();
		List<LinearConstraint> constraints = new ArrayList<>();

		// for each res pair
		for (ResidueInteractions.Pair resPair : inters) {
			Residue res1 = voxel.pmol.mol.residues.getOrThrow(resPair.resNum1);
			Residue res2 = voxel.pmol.mol.residues.getOrThrow(resPair.resNum2);

			// for each atom pair
			for (int[] atomPair : connectivity.getAtomPairs(res1, res2).getPairs(AtomNeighbors.Type.NONBONDED)) {
				Atom a1 = res1.atoms.get(atomPair[0]);
				Atom a2 = res2.atoms.get(atomPair[1]);

				// get voxel info for each atom, or skip the pair if no dofs
				AtomVoxel v1 = atomVoxels.computeIfAbsent(a1, (key) -> new AtomVoxel(a1, probe, voxel));
				AtomVoxel v2 = atomVoxels.computeIfAbsent(a2, (key) -> new AtomVoxel(a2, probe, voxel));
				if (!v1.hasDofs() && !v2.hasDofs()) {
					continue;
				}

				AtomPairVoxel pairVoxel = new AtomPairVoxel(voxel, v1, v2);
				LinearConstraint constraint = getLinearConstraint(pairVoxel, tolerance);
				if (constraint != null) {
					constraints.add(constraint);
				}
			}
		}

		return constraints;
	}

	public LinearConstraint getLinearConstraint(AtomPairVoxel voxel, double tolerance) {

		BoundaryPoint p = findBoundaryNewton(voxel, tolerance);

		if (p == null) {
			// dofs don't affect this atom pair, not useful for pruning, so don't make a constraint at all
			return null;
		}

		// did we not find a zero-point in the violation function?
		if (!p.atBoundary()) {

			// if we didn't find a zero-point, then assume all zero points lie outside the voxel
			// since all voxel points correspond to violations, this constraint is unsatisfiable
			if (p.violation > 0.0) {
				throw new NoFeasibleSolutionException();
			}

			// voxel always satisfiable for this atom pair, no constraint needed
			return null;
		}

		// use the boundary point to make a linear constraint on the dofs
		int n = p.dofValues.length;

		// make the linear constraint u.x >= w, where:
		//    u = -g
		//    w = -g.x*
		//    x* is the boundary point where the atom pair overlap is approx 0
		//    g is the gradient at x*
		// ie, the tangent hyperplane (d-1 linear subspace) to the isosurface at this point in the violation function
		RealVector u = new ArrayRealVector(n);
		double w = 0.0;
		for (int d=0; d<n; d++) {
			double g = -p.gradient[d];
			u.setEntry(d, g);
			w += p.dofValues[d]*g;
		}

		return new LinearConstraint(u, Relationship.GEQ, w);
	}

	public static class BoundaryPoint {

		double[] dofValues;
		double violation;
		double[] gradient;

		public BoundaryPoint(double[] dofValues, double violation, double[] gradient) {
			this.dofValues = dofValues.clone();
			this.violation = violation;
			this.gradient = gradient.clone();
		}

		public BoundaryPoint(double violation) {
			this.dofValues = null;
			this.violation = violation;
			this.gradient = null;
		}

		public boolean atBoundary() {
			return gradient != null;
		}

		@Override
		public String toString() {
			StringBuilder buf = new StringBuilder();
			buf.append(String.format("violation: %.3f", violation));
			if (dofValues != null) {
				for (int d=0; d<dofValues.length; d++) {
					buf.append(String.format("\n\tx[%d]: %.3f", d, dofValues[d]));
				}
			}
			if (gradient != null) {
				for (int d=0; d<gradient.length; d++) {
					buf.append(String.format("\n\tg[%d]: %.3f", d, gradient[d]));
				}
			}
			return buf.toString();
		}
	}

	/**
	 * find a point in the voxel where the atom pair overlap is close to 0
	 * using a generalization of Newton's method starting at the voxel center
	 */
	public BoundaryPoint findBoundaryNewton(AtomPairVoxel pairVoxel, double tolerance) {

		// objective function
		Function<double[],Double> f = (x) ->
			pairVoxel.getViolation(x, tolerance);

		// gradient function (approximate)
		double[] gout = new double[pairVoxel.numDofs];
		Function<double[],double[]> g = (x) -> {
			double baseViolation = f.apply(x);
			for (int d=0; d<pairVoxel.numDofs; d++) {
				double dx = gradientDxFactor*pairVoxel.width(d);
				gout[d] = (pairVoxel.getViolationAlong(d, x[d], dx, tolerance) - baseViolation)/dx;
			}
			return gout;
		};

		// start at the center of the voxel
		double[] x = new double[pairVoxel.numDofs];
		for (int d=0; d<pairVoxel.numDofs; d++) {
			x[d] = pairVoxel.center(d);
		}

		// get the initial violation
		double violation = f.apply(x);

		if (violation == 0.0) {
			// well, that was easy
			return new BoundaryPoint(x, violation, g.apply(x));
		}

		// iterate a newton-like method until convergence, or we hit a voxel boundary
		for (int i=0; i<maxNumIterations; i++) {

			// calculate the gradient at x
			double[] grad = g.apply(x);

			double s = 0.0;
			for (int d=0; d<pairVoxel.numDofs; d++) {
				s += pairVoxel.width2(d)*grad[d]*grad[d];
			}

			if (s == 0.0) {
				// gradient was zero (assuming no voxel width was zero), so can't find f=0 point
				// this usually means the DoFs have no effect on distance for this atom pair, and hence no boundary point exists
				return null;
			}

			// take a step along the gradient direction
			for (int d=0; d<pairVoxel.numDofs; d++) {
				x[d] -= violation*grad[d]*pairVoxel.width2(d)/s;
			}

			// did we step out of the voxel?
			if (pairVoxel.outOfRange(x)) {

				// the gradient suggests all boundary points are outside of the voxel
				// so just return the violation of the initial point, and assume that represents the entire voxel
				return new BoundaryPoint(violation);
			}

			// is the current violation close enough to 0?
			violation = f.apply(x);
			double diff = Math.abs(violation);
			if (diff <= violationThreshold) {
				return new BoundaryPoint(x, violation, grad);
			}
		}

		// ran out of iterations and didn't find a boundary
		final boolean developerIsInvestigatingFrequencyOfThisHappening = false;
		if (developerIsInvestigatingFrequencyOfThisHappening) {
			throw new Error("can't find boundary point after " + maxNumIterations + " iterations for " + pairVoxel.probePair);
		}

		// don't model this atom pair with a constraint at all, just to be conservative
		return null;
	}
}
