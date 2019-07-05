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

import java.math.BigInteger;
import java.util.*;

import edu.duke.cs.osprey.confspace.*;

/**
 *
 * @author mhall44
 */
public class PruningMatrix extends TupleMatrixBoolean {
	
	private static final long serialVersionUID = 1212622649775905467L;
	
    //similar to energy matrix, but indicates what RCs and tuples of RCs are pruned
    //pruning indicated by true boolean
    //a conformation is pruned if it contains any pruned RC or tuple
    
    //private HigherTupleFinder[][][][] higherTerms;//look up higher terms by pair
    
    //maybe separate intra too?
    
	protected PruningMatrix(){//no allocation (for overriding by UpdatedPruningMatrix)
		super();
    }
	
	public PruningMatrix(PruningMatrix other) {
		super(other);
	}
        
    public PruningMatrix(ConfSpace cSpace, double pruningInterval){
        super(cSpace, pruningInterval, false);//higher tuples are unpruned unless otherwise indicated
        
        //We'll want to initialize everything to be unpruned, because this will be looked up during pruning
        //currently all entries in oneBody and pairwise are null
        fill(false);
    }

	public PruningMatrix(SimpleConfSpace confSpace) {
    	super(confSpace, 0, false);

    	// start with everything unpruned
		fill(false);
	}
    
    public PruningMatrix(int numPos, int[] numAllowedAtPos, double pruningInterval) {
    	super(numPos, numAllowedAtPos, pruningInterval, false);
    }

    @Override
	public Boolean getTuple(RCTuple tuple) {
		Boolean val = super.getTuple(tuple);
		if (val != null) {
			return val;
		}
		return false;
	}

	public void unprunedRCsAtPos(ArrayList<Integer> out, int pos) {
    	out.clear();
    	int numRCs = getNumConfAtPos(pos);
		for (int index=0; index<numRCs; index++) {
			if(!getOneBody(pos,index))
				out.add(index);
		}
    }
    
    
    public ArrayList<Integer> unprunedRCsAtPos(int pos){
        //which RCs at the given position are unpruned?
        //Return index of the RCs within the position
        ArrayList<Integer> out = new ArrayList<>();
        unprunedRCsAtPos(out, pos);
        return out;
    }
    

    public void prunedRCsAtPos(ArrayList<Integer> out, int pos) {
    	out.clear();
    	int numRCs = getNumConfAtPos(pos);
		for (int index=0; index<numRCs; index++) {
			if(getOneBody(pos,index))
				out.add(index);
		}
    }
    
    
    public ArrayList<Integer> prunedRCsAtPos(int pos){
        //which RCs at the given position are unpruned?
        //Return index of the RCs within the position
        ArrayList<Integer> out = new ArrayList<>();
        prunedRCsAtPos(out, pos);
        return out;
    }
    
    
    public ArrayList<RCTuple> unprunedRCTuplesAtPos(ArrayList<Integer> pos){
        //get a list of unpruned RCTuples with the given positions
        //this method tests a few things more than once, so it could be sped up if needed, but it is convenient
        
        int numPos = pos.size();
        ArrayList<RCTuple> unpruned = new ArrayList<>();
        
        if(numPos==1){
            int posNum = pos.get(0);
            for(int rc=0; rc<getNumConfAtPos(posNum); rc++){
                if(!getOneBody(posNum,rc))
                    unpruned.add(new RCTuple(posNum,rc));
            }
        }
        else {
            //get unpruned tuples of RCs at all but the last position
            //then see what RCs at the last position we can add
            ArrayList<Integer> posReduced = (ArrayList<Integer>)pos.clone();
            posReduced.remove(numPos-1);
            
            ArrayList<RCTuple> tupsReduced = unprunedRCTuplesAtPos(posReduced);
            
            int lastPos = pos.get(numPos-1);
            
            for(int rc=0; rc<getNumConfAtPos(lastPos); rc++){
                if(!getOneBody(lastPos,rc)){
                    for(RCTuple reducedTup : tupsReduced){//try to combine into an unpruned RC
                        
                        ArrayList<Integer> fullRCList = (ArrayList<Integer>)reducedTup.RCs.clone();
                        fullRCList.add(rc);
                        
                        RCTuple fullTup = new RCTuple(pos,fullRCList);
                        if(!isPruned(fullTup))
                            unpruned.add(fullTup);
                    }
                }
            }
        }
        
        return unpruned;     
    }
    
    
    public boolean isPruned(RCTuple tup) {
    	
    	// OPTIMIZATION: this function gets hit a lot
    	// so even pedantic optimizations can have a noticeable impact
    	
    	// copy some references to stack
    	ArrayList<Integer> tuppos = tup.pos;
    	ArrayList<Integer> tupRCs = tup.RCs;
    	
    	// OPTIMIZATION: skipping even the check for higher order terms
    	// improves CPU cache performance a lot when we don't actually have any terms to use
    	boolean hasHigherOrderTerms = hasHigherOrderTerms();
    	
        //can be prune per se, or check if some singles in it are pruned, or pairs, etc.
    	
    	// OPTIMIZATION: checking singles separately from pairs improves CPU cache performance
    	// since we're memory-bound for this workload, trading the extra CPU instructions for
    	// fewer cache misses has a noticeable impact on performance
    	
    	// check singles
    	int numTupPos = tuppos.size();
        for (int i1=0; i1<numTupPos; i1++) {
            int pos1 = tuppos.get(i1);
            int rc1 = tupRCs.get(i1);
            
            if (getOneBody(pos1, rc1)) {
                return true;
            }
        }
            
        // check pairs
        for (int i1=0; i1<numTupPos; i1++) {
            int pos1 = tuppos.get(i1);
            int rc1 = tupRCs.get(i1);
            
            for (int i2=0; i2<i1; i2++) {
                int pos2 = tuppos.get(i2);
                int rc2 = tupRCs.get(i2);
            
                if (getPairwise(pos1, rc1, pos2, rc2)) {
                    return true;
                }
                
                if (hasHigherOrderTerms) {
					HigherTupleFinder<Boolean> htf = getHigherOrderTerms(pos1, rc1, pos2, rc2);
					if (htf != null) {
						if (isPrunedHigherOrder(tup, i2, htf)) {
							return true;
						}
					}
                }
            }
        }

        // check triples if needed using tuple trees
		if (hasHigherOrderTuples()) {

        	// allocate just one tuple, but update it inside the loops
        	RCTuple tuple = new RCTuple(0, 0, 0, 0, 0, 0);

			for (int i1=2; i1<numTupPos; i1++) {
				int pos1 = tuppos.get(i1);
				int rc1 = tupRCs.get(i1);

				// update tuple in reverse order, so positions are sorted
				tuple.pos.set(2, pos1);
				tuple.RCs.set(2, rc1);

				for (int i2=1; i2<i1; i2++) {
					int pos2 = tuppos.get(i2);
					int rc2 = tupRCs.get(i2);

					tuple.pos.set(1, pos2);
					tuple.RCs.set(1, rc2);

					for (int i3=0; i3<i2; i3++) {
						int pos3 = tuppos.get(i3);
						int rc3 = tupRCs.get(i3);

						tuple.pos.set(0, pos3);
						tuple.RCs.set(0, rc3);

						if (getTuple(tuple)) {
							return true;
						}
					}
				}
			}
		}
        
        return false;
    }
    
    
    public boolean isPrunedHigherOrder(RCTuple tup, int curIndex, HigherTupleFinder<Boolean> htf){
        //Checks if tup is pruned based on interactions in htf (corresponds to some sub-tuple of tup)
        //with RCs whose indices in tup are < curIndex
        ArrayList<Integer> interactingPos = htf.getInteractingPos();
        
        for(int ipos : interactingPos){
            
            //see if ipos is in tup with index < curIndex
            int iposIndex = -1;
            for(int ind=0; ind<curIndex; ind++){
                if(tup.pos.get(ind)==ipos){
                    iposIndex = ind;
                    break;
                }
            }

            if(iposIndex > -1){//ipos interactions need to be counted
                int iposRC = tup.RCs.get(iposIndex);
                if( htf.getInteraction(ipos, iposRC) )//sub-tuple plus (ipos,iposRC) is pruned
                    return true;
                
                //see if need to go up to highers order again...
                HigherTupleFinder htf2 = htf.getHigherInteractions(ipos,iposRC);
                if(htf2!=null){
                    if( isPrunedHigherOrder(tup,iposIndex,htf2) )
                        return true;
                }
            }
        }
        
        //if we get here, not pruned
        return false;
    }

    public void pruneSingle(int pos1, int rc1) {
    	setOneBody(pos1, rc1, true);
    	prunePairsFromSingle(pos1, rc1);
	}

	public void prunePair(int pos1, int rc1, int pos2, int rc2) {
    	setPairwise(pos1, rc1, pos2, rc2, true);
	}

	public void pruneTriple(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
		this.setTuple(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3).sorted(), true);
	}

	public void prunePairsFromSingles() {
		int n = getNumPos();
		for (int pos1=0; pos1<n; pos1++) {
			int n1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<n1; rc1++) {

				if (getOneBody(pos1, rc1)) {
					prunePairsFromSingle(pos1, rc1);
				}
			}
		}
	}

	public void prunePairsFromSingle(int pos1, int rc1) {
		int n = getNumPos();
		for (int pos2=0; pos2<n; pos2++) {

			if (pos1 == pos2) {
				continue;
			}
			int n2 = getNumConfAtPos(pos2);
			for (int rc2=0; rc2<n2; rc2++) {
				setPairwise(pos1, rc1, pos2, rc2, true);
			}
		}
	}

	public boolean isSinglePruned(int pos, int rc) {
    	return getOneBody(pos, rc);
	}

	public boolean isPairPruned(int pos1, int rc1, int pos2, int rc2) {
    	return getPairwise(pos1, rc1, pos2, rc2)
			|| isSinglePruned(pos1, rc1)
			|| isSinglePruned(pos2, rc2);
	}

	public boolean isTriplePruned(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
		return isSinglePruned(pos1, rc1)
			|| isSinglePruned(pos2, rc2)
			|| isSinglePruned(pos3, rc3)
			|| getPairwise(pos1, rc1, pos2, rc2)
			|| getPairwise(pos1, rc1, pos3, rc3)
			|| getPairwise(pos2, rc2, pos3, rc3)
			|| getTuple(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3).sorted());
	}

	public boolean isQuadruplePruned(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, int pos4, int rc4) {
		RCTuple tuple = new RCTuple(0, 0, 0, 0, 0, 0);
		return isSinglePruned(pos1, rc1)
			|| isSinglePruned(pos2, rc2)
			|| isSinglePruned(pos3, rc3)
			|| isSinglePruned(pos4, rc4)
			|| getPairwise(pos1, rc1, pos2, rc2)
			|| getPairwise(pos1, rc1, pos3, rc3)
			|| getPairwise(pos1, rc1, pos4, rc4)
			|| getPairwise(pos2, rc2, pos3, rc3)
			|| getPairwise(pos2, rc2, pos4, rc4)
			|| getPairwise(pos3, rc3, pos4, rc4)
			|| getTuple(tuple.set(pos1, rc1, pos2, rc2, pos3, rc3).sorted())
			|| getTuple(tuple.set(pos1, rc1, pos2, rc2, pos4, rc4).sorted())
			|| getTuple(tuple.set(pos1, rc1, pos3, rc3, pos4, rc4).sorted())
			|| getTuple(tuple.set(pos2, rc2, pos3, rc3, pos4, rc4).sorted());
	}

	public boolean isQuintuplePruned(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, int pos4, int rc4, int pos5, int rc5) {
		RCTuple tuple = new RCTuple(0, 0, 0, 0, 0, 0);
		return isSinglePruned(pos1, rc1)
			|| isSinglePruned(pos2, rc2)
			|| isSinglePruned(pos3, rc3)
			|| isSinglePruned(pos4, rc4)
			|| isSinglePruned(pos5, rc5)
			|| getPairwise(pos1, rc1, pos2, rc2)
			|| getPairwise(pos1, rc1, pos3, rc3)
			|| getPairwise(pos1, rc1, pos4, rc4)
			|| getPairwise(pos1, rc1, pos5, rc5)
			|| getPairwise(pos2, rc2, pos3, rc3)
			|| getPairwise(pos2, rc2, pos4, rc4)
			|| getPairwise(pos2, rc2, pos5, rc5)
			|| getPairwise(pos3, rc3, pos4, rc4)
			|| getPairwise(pos3, rc3, pos5, rc5)
			|| getPairwise(pos4, rc4, pos5, rc5)
			|| getTuple(tuple.set(pos1, rc1, pos2, rc2, pos3, rc3).sorted())
			|| getTuple(tuple.set(pos1, rc1, pos2, rc2, pos4, rc4).sorted())
			|| getTuple(tuple.set(pos1, rc1, pos2, rc2, pos5, rc5).sorted())
			|| getTuple(tuple.set(pos1, rc1, pos3, rc3, pos4, rc4).sorted())
			|| getTuple(tuple.set(pos1, rc1, pos3, rc3, pos5, rc5).sorted())
			|| getTuple(tuple.set(pos1, rc1, pos4, rc4, pos5, rc5).sorted())
			|| getTuple(tuple.set(pos2, rc2, pos3, rc3, pos4, rc4).sorted())
			|| getTuple(tuple.set(pos2, rc2, pos3, rc3, pos5, rc5).sorted())
			|| getTuple(tuple.set(pos2, rc2, pos4, rc4, pos5, rc5).sorted())
			|| getTuple(tuple.set(pos3, rc3, pos4, rc4, pos5, rc5).sorted());
	}

	public void markAsPruned(RCTuple tup){
        setTupleValue(tup, true);
        /*
        int tupSize = tup.pos.size();
        if(tupSize==1)
            setOneBody(tup.pos.get(0), tup.RCs.get(0), true);
        else if(tupSize==2)
            setPairwise(tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1), true);
        else
        */
    }
    
    
    public int countPrunedRCs(){
        //how many RCs are pruned overall?
        int count = 0;
        int numPos = getNumPos();
        for (int res1=0; res1<numPos; res1++) {
        	int m1 = getNumConfAtPos(res1);
        	for (int i1=0; i1<m1; i1++) {
        		if (getOneBody(res1, i1) == true) {
        			count++;
        		}
        	}
        }
        return count;
    }
    
    
    public int countPrunedPairs(){
        //how many pairs are pruned overall?
        int count = 0;
        int numPos = getNumPos();
        for (int res1=0; res1<numPos; res1++) {
        	int m1 = getNumConfAtPos(res1);
        	for (int i1=0; i1<m1; i1++) {
        		for (int res2=0; res2<res1; res2++) {
        			int m2 = getNumConfAtPos(res2);
        			for (int i2=0; i2<m2; i2++) {
        				if (getPairwise(res1, i1, res2, i2) == true) {
        					count++;
        				}
        			}
        		}
        	}
        }
        return count;
    }

	public int countPrunedTriples() {

		int numPos = getNumPos();
		int count = 0;
		RCTuple tuple = new RCTuple(0, 0, 0, 0, 0, 0);

		for (int pos1=0; pos1<getNumPos(); pos1++) {
			tuple.pos.set(0, pos1);
			for (int rc1=0; rc1<getNumConfAtPos(pos1); rc1++) {
				tuple.RCs.set(0, rc1);

				for (int pos2=pos1+1; pos2<numPos; pos2++) {
					tuple.pos.set(1, pos2);
					for (int rc2=0; rc2<getNumConfAtPos(pos2); rc2++) {
						tuple.RCs.set(1, rc2);

						TupleTree<Boolean> tree = getHigherOrderTuples(pos1, rc1, pos2, rc2);
						if (tree == null) {
							continue;
						}

						for (int pos3=pos2+1; pos3<numPos; pos3++) {
							tuple.pos.set(2, pos3);
							for (int rc3=0; rc3<getNumConfAtPos(pos3); rc3++) {
								tuple.RCs.set(2, rc3);

								Boolean val = tree.get(tuple);
								if (val != null && val == true) {
									count++;
								}
							}
						}
					}
				}
			}
		}

		return count;
	}

    /*boolean isPruned(RC rc){
        //look up 1-body
        return getOneBody(pos,rcNum);
    }*/

    public int countUnprunedSingles(int pos) {
    	int count = 0;
    	for (int rc=0; rc<getNumConfAtPos(pos); rc++) {
    		if (!isSinglePruned(pos, rc)) {
    			count++;
			}
		}
    	return count;
	}

	/**
	 * Calculates an upper bound on the number of conformations remaining in
	 * the conformation space after all singles and pairs pruning.
	 *
	 * Uses a dynamic programming algorithm that runs in O(nm^2) time and O(nm) space
	 * where n is the number of design positions, and m is the max number of RCs at a position.
	 *
	 * This algorithm essentially under-estimates pruning by ignoring pruned pairs
	 * that are not between consecutive positions. The then the leftover pruned pairs
	 * have a nice structure that can be analyzed in polynomial time using dynamic programming.
	 *
	 * NOTE: Counting the exact size of the conformation space after pairs pruning is
	 * apparently #P-complete, and hence very hard to do. =P
	 */
	public BigInteger calcUnprunedConfsUpperBound() {

		// start with the canonical permutation
		int n = getNumPos();
		List<Integer> posPermutation = new ArrayList<>(n);
		for (int i=0; i<n; i++) {
			posPermutation.add(i);
		}

		// since we're ignoring some pruned pairs to get a bound,
		// we have some freedom to choose which pairs to ignore

		// by enforcing an order for positions,
		// we're choosing to ignore pairs between non-consecutive positions
		// the consecutive positions then form the "diagonal" of a position pair matrix

		// ideally, we want to choose the position ordering that minimizes the upper bound
		// but I don't know how to solve that problem optimally yet
		// for now, let's just try a new greedy heuristic:

		// let's sort positions by number of RCs
		posPermutation.sort(Comparator.comparing(pos -> getNumConfAtPos(pos)));

		return countConfsAfterOnlyDiagonalPairsPruning(this, posPermutation);
	}

	private static BigInteger countConfsAfterOnlyDiagonalPairsPruning(PruningMatrix pmat, List<Integer> posPermutation) {

		// allocate space for the intermediate counts
		// (columns are design positions, rows are RCs at that positions)
		int n = pmat.getNumPos();
		BigInteger[][] counts = new BigInteger[n][];
		for (int i=0; i<n; i++) {
			int pos = posPermutation.get(i);
			counts[i] = new BigInteger[pmat.getNumConfAtPos(pos)];
		}

		// initialize all counts to zero, except for RCs at the last position, which are one
		for (int i=0; i<n; i++) {
			int pos = posPermutation.get(i);
			for (int rc=0; rc<pmat.getNumConfAtPos(pos); rc++) {
				if (i < n - 1) {
					counts[i][rc] = BigInteger.ZERO;
				} else {
					counts[i][rc] = BigInteger.ONE;
				}
			}
		}

		// staring at the second-to-last position and moving backwards...
		for (int i1 = n - 2; i1 >= 0; i1--) {
			int i2 = i1 + 1;
			int pos1 = posPermutation.get(i1);
			int pos2 = posPermutation.get(i2);

			// for each RC at this position...
			for (int rc1=0; rc1<pmat.getNumConfAtPos(pos1); rc1++) {

				// skip pruned singles
				if (pmat.isSinglePruned(pos1, rc1)) {
					continue;
				}

				// sum the counts of the unpruned RCs at the next position
				// and store the sum at this pos,rc

				// for each RC at the next position...
				for (int rc2=0; rc2<pmat.getNumConfAtPos(pos2); rc2++) {

					// skip pruned singles and pairs
					if (pmat.isSinglePruned(pos2, rc2) || pmat.isPairPruned(pos1, rc1, pos2, rc2)) {
						continue;
					}

					// update the intermediate counts at this pos
					counts[i1][rc1] = counts[i1][rc1].add(counts[i2][rc2]);
				}
			}
		}

		// add counts from all RCs at the first position to get the total count
		BigInteger sum = BigInteger.ZERO;
		for (int rc=0; rc<pmat.getNumConfAtPos(posPermutation.get(0)); rc++) {
			sum = sum.add(counts[0][rc]);
		}

		return sum;
	}

	/**
	 * calculates an lower bound on the number of conformations remaining in
	 * the conformation space after all singles and pairs pruning.
	 *
	 * Uses the same dynamic programming algorithm as the upper bound, but we
	 * perform a preprocessing step on the pruning matrix first.
	 *
	 * This algorithm essentially over-estimates pruning, by "upgrading" pruned pairs
	 * not between consecutive positions to pruned singles. The then the leftover pruned pairs
	 * have a nice structure that can be analyzed in polynomial time.
	 *
	 * NOTE: Counting the exact size of the conformation space after pairs pruning is
	 * apparently #P-complete, and hence very hard to do. =P
	 */
	public BigInteger calcUnprunedConfsLowerBound() {

		// first, transform the pruning matrix
		PruningMatrix expandedPmat = new PruningMatrix(this);

		// for each "off-diagonal" position pair (e.g. "diagonal" pairs are 01, 12, 23, etc),
		// create a pruned single for each pruned pair that "covers" the tuples pruned by the pair
		// (inevitably some extra tuples get pruned this way, but that's why we're over-estimating pruning)
		int n = getNumPos();
		for (int pos1=2; pos1<n; pos1++) {

			for (int pos2=0; pos2<pos1 - 1; pos2++) {

				for (int rc1=0; rc1<getNumConfAtPos(pos1); rc1++) {
					for (int rc2=0; rc2<getNumConfAtPos(pos2); rc2++) {

						if (getPairwise(pos1, rc1, pos2, rc2)) {

							// replace this pruned pair with the pruned single that
							// minimizes the number of extra pruned tuples

							// but which is better, pos1, or pos2?

							// always prefer a single that has already been pruned
							if (expandedPmat.isSinglePruned(pos1, rc1)) {
								// already pruned, nothing to do
							} else if (expandedPmat.isSinglePruned(pos2, rc2)) {
								// already pruned, nothing to do
							} else {

								// otherwise, pick the pos with smaller number of unpruned RCs
								if (expandedPmat.countUnprunedSingles(pos1) < expandedPmat.countUnprunedSingles(pos2)) {
									expandedPmat.pruneSingle(pos1, rc1);
								} else {
									expandedPmat.pruneSingle(pos2, rc2);
								}
							}
						}
					}
				}
			}
		}

		// use the canonical permutation
		List<Integer> posPermutation = new ArrayList<>(n);
		for (int i=0; i<n; i++) {
			posPermutation.add(i);
		}

		return countConfsAfterOnlyDiagonalPairsPruning(expandedPmat, posPermutation);
	}


	/*
	 * Iteration that is simple to express in a for-loop, is often devilishly difficult to express as a
	 * java.util.Iterator subclass, since Java doesn't have continuations or suspendable functions like other languages.
	 * These functions try to keep the simple for-loop expressions for iteration over tuples, but still allow convenient
	 * iteration using functional-style programming (ie, passing in functions to the iterator).
	 *
	 * sure, we could just use a for-loop to buffer tuples into a list which has easy access to an Iterator,
	 * but that's horribly inefficient. let's do better! =)
	 *
	 * The argument function must return a command to the iterator at each application: Either Stop to terminate
	 * the iterator, or Continue to ask to be called again with the next tuple.
	 */

	public static enum IteratorCommand {
		Break,
		Continue
	}


	// FOR SINGLES

	public static interface SingleConsumer {
		public IteratorCommand apply(int pos, int rc);
	}

	public IteratorCommand forEachUnprunedSingle(SingleConsumer consumer) {
		int n = getNumPos();
		for (int pos1=0; pos1<n; pos1++) {
			switch (forEachUnprunedSingleAt(pos1, consumer)) {
				case Break: return IteratorCommand.Break;
			}
		}
		return IteratorCommand.Continue;
	}

	public IteratorCommand forEachUnprunedSingleAt(int pos1, SingleConsumer consumer) {
		for (int rc1=0; rc1<getNumConfAtPos(pos1); rc1++) {

			// skip pruned stuff
			if (isSinglePruned(pos1, rc1)) {
				continue;
			}

			switch (consumer.apply(pos1, rc1)) {
				case Break: return IteratorCommand.Break;
			}
		}
		return IteratorCommand.Continue;
	}


	// FOR PAIRS

	public static interface PairConsumer {
		public IteratorCommand apply(int pos1, int rc1, int pos2, int rc2);
	}

	public IteratorCommand forEachUnprunedPair(PairConsumer consumer) {
		for (int pos1=1; pos1<getNumPos(); pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {
				switch (forEachUnprunedPairAt(pos1, pos2, consumer)) {
					case Break: return IteratorCommand.Break;
				}
			}
		}
		return IteratorCommand.Continue;
	}

	public IteratorCommand forEachUnprunedPairAt(int pos1, int pos2, PairConsumer consumer) {
		for (int rc1=0; rc1<getNumConfAtPos(pos1); rc1++) {

			// skip pruned stuff
			if (getOneBody(pos1, rc1)) {
				continue;
			}

			for (int rc2=0; rc2<getNumConfAtPos(pos2); rc2++) {

				// skip pruned stuff
				if (getOneBody(pos2, rc2) || getPairwise(pos1, rc1, pos2, rc2)) {
					continue;
				}

				switch (consumer.apply(pos1, rc1, pos2, rc2)) {
					case Break: return IteratorCommand.Break;
				}
			}
		}
		return IteratorCommand.Continue;
	}

	public static interface TripleConsumer {
		public IteratorCommand apply(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3);
	}

	public IteratorCommand forEachUnprunedTriple(TripleConsumer consumer) {
		for (int pos1=2; pos1<getNumPos(); pos1++) {
			for (int pos2=1; pos2<pos1; pos2++) {
				for (int pos3=0; pos3<pos2; pos3++) {
					switch (forEachUnprunedTripleAt(pos1, pos2, pos3, consumer)) {
						case Break: return IteratorCommand.Break;
					}
				}
			}
		}
		return IteratorCommand.Continue;
	}

	public IteratorCommand forEachUnprunedTripleAt(int pos1, int pos2, int pos3, TripleConsumer consumer) {

		// allocate one tuple for the whole function, but update it at each iteration
		RCTuple tuple = new RCTuple(0, 0, 0, 0, 0, 0);

		for (int rc1=0; rc1<getNumConfAtPos(pos1); rc1++) {

			// skip pruned stuff
			if (getOneBody(pos1, rc1)) {
				continue;
			}

			// update the tuple in reverse order, so tuple positions are sorted
			tuple.pos.set(2, pos1);
			tuple.RCs.set(2, rc1);

			for (int rc2=0; rc2<getNumConfAtPos(pos2); rc2++) {

				// skip pruned stuff
				if (getOneBody(pos2, rc2) || getPairwise(pos1, rc1, pos2, rc2)) {
					continue;
				}

				tuple.pos.set(1, pos2);
				tuple.RCs.set(1, rc2);

				for (int rc3=0; rc3<getNumConfAtPos(pos3); rc3++) {

					tuple.pos.set(0, pos3);
					tuple.RCs.set(0, rc3);

					// skip pruned stuff
					if (getOneBody(pos3, rc3)
						|| getPairwise(pos1, rc1, pos3, rc3)
						|| getPairwise(pos2, rc2, pos3, rc3)
						|| getTuple(tuple)
					) {
						continue;
					}

					switch (consumer.apply(pos1, rc1, pos2, rc2, pos3, rc3)) {
						case Break: return IteratorCommand.Break;
					}
				}
			}
		}
		return IteratorCommand.Continue;
	}
}
