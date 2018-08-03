package edu.duke.cs.osprey.pruning;


import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.structure.*;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.linear.*;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;


public class PLUG {

	public final SimpleConfSpace confSpace;

	/** max num iterations of generalized Newton iteration to find a "boundary" point where the violation function is zero */
	public int maxNumIterations = 30;

	/** distance threshold to claim that a violation value is close enough to zero */
	public double violationThreshold = 1e-2;

	/** percentage of the voxel width used to approximate the gradient of the violation function */
	public double gradientDxPercent = 1e-4;

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
		pmat.forEachUnprunedSingle((pos1, rc1) -> {
			if (shouldPruneTuple(new RCTuple(pos1, rc1), tolerance)) {
				pmat.pruneSingle(pos1, rc1);
			}
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	public void prunePairs(PruningMatrix pmat, double tolerance) {
		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			if (shouldPruneTuple(new RCTuple(pos1, rc1, pos2, rc2), tolerance)) {
				pmat.prunePair(pos1, rc1, pos2, rc2);
			}
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	public void pruneTriples(PruningMatrix pmat, double tolerance) {
		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			if (shouldPruneTuple(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3), tolerance)) {
				pmat.pruneTriple(pos1, rc1, pos2, rc2, pos3, rc3);
			}
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	public boolean shouldPruneTuple(RCTuple tuple, double tolerance) {

		// make the molecule and get all the residue interactions for the tuple
		ParametricMolecule pmol = confSpace.makeMolecule(tuple);
		ResidueInteractions inters = ResInterGen.of(confSpace)
			.addIntras(tuple)
			.addInters(tuple)
			.addShell(tuple)
			.make();

		try {

			// get linear constraints for each atom pair
			List<LinearConstraint> constraints = getLinearConstraints(pmol, inters, tolerance);

			// no constraints? don't prune
			if (constraints.isEmpty()) {
				return false;
			}

			// use an LP solver (eg simplex) to determine if the constraints allow any feasible points
			// TODO: get the simple bounds from the objective function obj?
			double[] vmin = new double[pmol.dofs.size()];
			double[] vmax = new double[pmol.dofs.size()];
			for (int d=0; d<pmol.dofs.size(); d++) {
				vmin[d] = pmol.dofBounds.getMin(d);
				vmax[d] = pmol.dofBounds.getMax(d);
			}
			new SimplexSolver().optimize(
				new SimpleBounds(vmin, vmax), // voxel bounds
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

	public List<LinearConstraint> getLinearConstraints(ParametricMolecule pmol, ResidueInteractions inters, double tolerance) {

		List<LinearConstraint> constraints = new ArrayList<>();
		for (ResidueInteractions.Pair resPair : inters) {

			Residue res1 = pmol.mol.residues.getOrThrow(resPair.resNum1);
			Residue res2 = pmol.mol.residues.getOrThrow(resPair.resNum2);

			for (int[] atomPair : connectivity.getAtomPairs(res1, res2).getPairs(AtomNeighbors.Type.NONBONDED)) {

				Probe.AtomPair probePair = probe.new AtomPair(
					res1.atoms.get(atomPair[0]),
					res2.atoms.get(atomPair[1])
				);

				LinearConstraint constraint = getLinearConstraint(pmol, probePair, tolerance);
				if (constraint != null) {
					constraints.add(constraint);
				}
			}
		}
		return constraints;
	}

	public LinearConstraint getLinearConstraint(ParametricMolecule pmol, Probe.AtomPair pair, double tolerance) {


		BoundaryPoint p = findBoundaryNewton(pmol, pair, tolerance);
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
	public BoundaryPoint findBoundaryNewton(ParametricMolecule pmol, Probe.AtomPair pair, double tolerance) {

		final int numDofs = pmol.dofs.size();

		// make voxel properties easy to access
		// TODO: cache this somehow? make objective function object, with swappable pairs?
		double[] vw = new double[numDofs];
		double[] vw2 = new double[numDofs];
		double[] vmin = new double[numDofs];
		double[] vmax = new double[numDofs];
		for (int d=0; d<numDofs; d++) {
			vw[d] = pmol.dofBounds.getMax(d) - pmol.dofBounds.getMin(d);
			vw2[d] = vw[d]*vw[d];
			vmin[d] = pmol.dofBounds.getMin(d);
			vmax[d] = pmol.dofBounds.getMax(d);
		}

		// objective function
		Function<double[],Double> f = (x) -> {
			for (int d=0; d<numDofs; d++) {
				pmol.dofs.get(d).apply(x[d]);
			}
			return pair.getViolation(tolerance);
		};

		// gradient
		Function<double[],double[]> g = (x) -> {
			double[] y = x.clone();
			double[] out = new double[numDofs];
			for (int d=0; d<numDofs; d++) {
				double dx = gradientDxPercent*vw[d];
				y[d] = x[d] - dx;
				double fdm = f.apply(y);
				y[d] = x[d] + dx;
				double fdp = f.apply(y);
				y[d] = x[d];
				out[d] = (fdp - fdm)/dx/2;
			}
			return out;
		};

		// start at the center of the voxel
		double[] x = new double[numDofs];
		for (int d=0; d<numDofs; d++) {
			x[d] = (vmin[d] + vmax[d])/2;
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
			for (int d=0; d<numDofs; d++) {
				s += vw2[d]*grad[d]*grad[d];
			}

			if (s == 0.0) {
				// gradient was zero (assuming no voxel width was zero), so can't find f=0 point
				// this usually means the DoFs have no effect on distance for this atom pair, and hence no boundary point exists
				return null;
			}

			// take a step along the gradient direction
			for (int d=0; d<numDofs; d++) {
				x[d] -= violation*grad[d]*vw2[d]/s;
			}

			// did we step out of the voxel?
			boolean outOfVoxel = false;
			for (int d=0; d<numDofs; d++) {
				if (x[d] < vmin[d] || x[d] > vmax[d]) {
					outOfVoxel = true;
					break;
				}
			}
			if (outOfVoxel) {

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

		// ran out of iterations
		// TODO: what to do here?
		throw new Error("ran out of iterations, didn't find f=0");
	}
}
