package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.math.BigDecimal;

public class SHARKStarBound extends MARKStarBoundFastQueues {

	private MARKStarNode precomputedRootNode;
	private SimpleConfSpace confSpace;

	/**
	 * Constructor to make a default MARKStarBound Class
	 * @param confSpace the partition function conformation space
	 * @param rigidEmat the rigid pairwise energy matrix
	 * @param minimizingEmat the parwise-minimized energy matrix
	 * @param minimizingConfEcalc the energy calculator to calculate minimized conf energies
	 * @param rcs information on possible rotamers at all design positions
	 * @param parallelism information for threading
	 */
	public SHARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
						 ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
	    super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
		this.confSpace = confSpace;
    }

	public SHARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
						  ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism,
						  SHARKStarBound precomputedFlex) {
		super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);

		/*
		Now we need to do a couple things.
		For now, let's assume we are working with a single sequence

		 */

		this.confSpace = confSpace;
		precomputedRootNode = precomputedFlex.rootNode;
		//updateConfTree(precomputedFlex.rootNode);
		updateBound();


		/*
		TODO: Go through the tree to make sure that the assignments are compatible with the new confspace

		TODO: Make sure all of the energy matrices (including the correction emats) are compatible with the new confspace

		TODO: Do something with score / object contexts?

		TODO: Go through and update bounds

		TODO: Populate queue
		 */
	}

	/**
	 * Returns the partition function lower bound for a particular sequence
	 *
	 * Note that SHARKStarBound will eventually contain a multi-sequence confTree, although this isn't currently the case
	 *
	 * @param seq Sequence for which to get pfunc lower bound
	 * @return BigDecimal pfunc lower bound
	 */
	public BigDecimal getLowerBound (Sequence seq){
		throw new UnsupportedOperationException("getLowerBound(seq) is not yet implemented");
	}

	/**
	 * Returns the partition function upper bound for a particular sequence
	 *
	 * Note that SHARKStarBound will eventually contain a multi-sequence confTree, although this isn't currently the case
	 *
	 * @param seq Sequence for which to get pfunc upper bound
	 * @return BigDecimal pfunc upper bound
	 */
	public BigDecimal getUpperBound (Sequence seq){
		throw new UnsupportedOperationException("getUpperBound(seq) is not yet implemented");
	}

	/**
	 * Returns the partition function lower bound for the whole confTree
	 *
	 * Note that SHARKStarBound will eventually contain a multi-sequence confTree, although this isn't currently the case
	 *
	 * @return BigDecimal pfunc lower bound
	 */
	public BigDecimal getLowerBound (){
		return rootNode.getLowerBound();
	}

	/**
	 * Returns the partition function upper bound for the whole confTree
	 *
	 * Note that SHARKStarBound will eventually contain a multi-sequence confTree, although this isn't currently the case
	 *
	 * @return BigDecimal pfunc upper bound
	 */
	public BigDecimal getUpperBound (){
		return rootNode.getUpperBound();
	}

	/**
	 * Returns the partition function lower bound for the precomputed confspace
	 *
	 * @return BigDecimal precomputed pfunc lower bound
	 */
	public BigDecimal getPrecomputedLowerBound (){
		return precomputedRootNode.getLowerBound();
	}

	/**
	 * Returns the partition function upper bound for the precomputed confTree
	 *
	 * @return BigDecimal precomputed pfunc upper bound
	 */
	public BigDecimal getPrecomputedUpperBound (){
		return precomputedRootNode.getUpperBound();
	}

	/**
	 * Makes the current confTree consistent with the current confSpace
	 *
	 * When we precompute flexible residues, we will have a tree that is for a flexible confspace.
	 * However, when we want to compute for mutable residues, we need to extend the length of assignments in our tree
	 */
	public void updateConfTree(){
		//System.out.println("The precomputed root node is " + precomputedRootNode.toTuple());
		precomputedRootNode.printTree();
		System.out.println("\n###############\nFull root upper: "+rootNode.getUpperBound()+" lower: "+rootNode.getLowerBound());

		/*
		Here's the plan: use the permutation matrix to map assignements onto the new tree.
		The g scores should map fine
		the h scores may need to be updated, which is kind of awkward I suppose.
		But, any full conformations should be minimized and should be added to the MAE energy matrix
		likewise for any partial minimizations
		 */
	}

	/**
	 * Generate a permutation matrix that lets us map positions from the precomputed confspace to the new confspace
	 */
	public int[] genConfSpaceMapping(){
		return null;
	}

}
