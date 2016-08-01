package edu.duke.cs.osprey.kstar;

import java.util.ArrayList;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.FullAStarNode;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */


@SuppressWarnings("serial")
public class KAStarConfTree extends ConfTree<FullAStarNode> {

	
	protected KSSearchProblem sp;
	protected PruningMatrix reducedPmat;
	protected PruningMatrix panPmat;
	protected ArrayList<ArrayList<Integer>> panUprunedRCsAtPos;
	protected boolean energyLB = true; // compute lb xor ub
	
	
	public KAStarConfTree(KSSearchProblem reducedSP, PruningMatrix reducedPmat, PruningMatrix panPmat) {
		super(new FullAStarNode.Factory(reducedPmat.getNumPos()), reducedSP, reducedPmat);
		
		this.sp = reducedSP;
		this.reducedPmat = reducedPmat;
		this.panPmat = panPmat;
		this.energyLB = reducedSP.contSCFlex ? true : false;
		initVars();
	}


	private void initVars() {
        // reduce the energy matrix, too
        emat = sp.getReducedEnergyMatrix();
        
        // cache this rather than recompute it
        panUprunedRCsAtPos = new ArrayList<>(panPmat.getNumPos());
        for(int index = 0; index < panPmat.getNumPos(); ++index) {
        	ArrayList<Integer> tmp = panPmat.unprunedRCsAtPos(index); tmp.trimToSize();
        	panUprunedRCsAtPos.add(tmp);
        }
    }
	
	
	protected ArrayList<Integer> getUndefinedPos( RCTuple definedTuple ) {
		// defined tuple positions are wrt reduced matrix positions
		ArrayList<Integer> ans = sp.getMaxPosNums();
		
		ArrayList<Integer> definedPos = new ArrayList<>(definedTuple.pos.size());
		for(int pos : definedTuple.pos) {
			// convert to absolute positions
			definedPos.add(sp.posNums.get(pos));
		}

		// remove defined positions
		ans.removeAll(definedPos);

		ans.trimToSize();
		return ans;
	}


	protected ArrayList<Integer> allowedRCsAtLevel(int level, int[] partialConf, ArrayList<Integer> undefinedPos) {
		// level is the absolute pos, defined WRT the all sequences pruning matrix
		ArrayList<Integer> allowedRCs;

		if(undefinedPos.contains(level)) { // level is undefined
			// allowedRCs = panPmat.unprunedRCsAtPos(level);
			allowedRCs = panUprunedRCsAtPos.get(level);
		}

		else {
			// level is defined (i.e. in posNums). its index is the level in the reducedSP
			allowedRCs = new ArrayList<>();
			
			// if the position in the partial conf is assigned, then return the solitary rc
			// else return the allowed rcs at that position in the reduced energy matrix
			int level2 = sp.posNums.indexOf(level);
			
			if(partialConf[level2] != -1) {
				// defined and assigned
				allowedRCs.add(partialConf[level2]);
			}
			
			else {
				// defined but unassigned
				//allowedRCs.addAll(unprunedRCsAtPos.get(level2));
				for(int i : unprunedRCsAtPos[level2]) allowedRCs.add(i);
			}
			
			// allowedRCs.trimToSize();
		}

		return allowedRCs;
	}


	protected double RCContribution(int level, int rc, RCTuple definedTuple, int[] partialConf, ArrayList<Integer> undefinedPos) {
		//Provide a lower bound on what the given rc at the given level can contribute to the energy
		//assume partialConf and definedTuple

		double rcContrib = sp.getEnergyMatrix().getOneBody(level, rc);

		//for this kind of lower bound, we need to split up the energy into the defined-tuple energy
		//plus "contributions" for each undefined residue
		//so we'll say the "contribution" consists of any interactions that include that residue
		//but do not include higher-numbered undefined residues
		for(int level2 = 0; level2 < sp.confSpace.numPos; level2++){

			if(!undefinedPos.contains(level2) || level2 < level) { //defined or lower numbered residues

				double levelBestE = energyLB ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;//best pairwise energy

				ArrayList<Integer> allowedRCs = allowedRCsAtLevel(level2, partialConf, undefinedPos);

				for( int rc2 : allowedRCs ) {

					double interactionE = sp.getEnergyMatrix().getPairwise(level, rc, level2, rc2);

					double higherOrderE = higherOrderContrib(level, rc, level2, rc2, partialConf, undefinedPos);
					//add higher-order terms that involve rc, rc2, and parts of partialConf

					interactionE += higherOrderE;

					//besides that only residues in definedTuple or levels below level2
					levelBestE = energyLB ? Math.min(levelBestE, interactionE) : Math.max(levelBestE, interactionE);
				}

				rcContrib += levelBestE;
			}
		}

		return rcContrib;
	}


	protected double higherOrderContrib(int pos1, int rc1, int pos2, int rc2, 
			int[] partialConf, ArrayList<Integer> undefinedPos) {
		//higher-order contribution for a given RC pair, when scoring a partial conf

		HigherTupleFinder<Double> htf = sp.getEnergyMatrix().getHigherOrderTerms(pos1, rc1, pos2, rc2);

		if(htf==null)
			return 0;//no higher-order interactions
		else {
			RCTuple curPair = new RCTuple(pos1, rc1, pos2, rc2);
			return higherOrderContrib(htf, curPair, partialConf, undefinedPos);
		}
	}


	@SuppressWarnings("unchecked")
	double higherOrderContrib(HigherTupleFinder<Double> htf, RCTuple startingTuple, int[] partialConf, ArrayList<Integer> undefinedPos) {
		//recursive function to get bound on higher-than-pairwise terms
		//this is the contribution to the bound due to higher-order interactions
		//of the RC tuple startingTuple (corresponding to htf)

		double contrib = 0;

		//to avoid double-counting, we are just counting interactions of starting tuple
		//with residues before the "earliest" one (startingLevel) in startingTuple
		//"earliest" means lowest-numbered, except non-mutating res come before mutating
		int startingLevel = startingTuple.pos.get( startingTuple.pos.size()-1 );

		for(int iPos : htf.getInteractingPos()){//position has higher-order interaction with tup
			if(posComesBefore(iPos,startingLevel,undefinedPos)) {//interaction in right order
				//(want to avoid double-counting)

				double levelBestE = energyLB ? Double.POSITIVE_INFINITY : 
					Double.NEGATIVE_INFINITY;//best value of contribution from tup-iPos interaction

				ArrayList<Integer> allowedRCs = allowedRCsAtLevel(iPos, partialConf, undefinedPos);

				for( int rc : allowedRCs ){

					RCTuple augTuple = startingTuple.addRC(iPos, rc);

					double interactionE = htf.getInteraction(iPos, rc);

					//see if need to go up to highers order again...
					@SuppressWarnings("rawtypes")
					HigherTupleFinder htf2 = htf.getHigherInteractions(iPos, rc);
					if(htf2!=null){
						interactionE += higherOrderContrib(htf2, augTuple, partialConf, undefinedPos);
					}

					levelBestE = energyLB ? Math.min(levelBestE, interactionE) : Math.max(levelBestE, interactionE);
				}

				contrib += levelBestE;//add up contributions from different interacting positions iPos
			}
		}

		return contrib;
	}


	protected boolean posComesBefore(int pos1, int pos2, ArrayList<Integer> undefinedPos){
		//for purposes of contributions to traditional conf score, 
		//we go through defined and then through undefined positions (in partialConf);
		//within each of these groups we go in order of position number
		if(!undefinedPos.contains(pos2)){//pos2 defined
			return (pos1<pos2 && !undefinedPos.contains(pos1));//pos1 must be defined to come before pos2
		}
		else//pos1 comes before pos2 if it's defined, or if pos1<pos2
			return (pos1<pos2 || !undefinedPos.contains(pos1));
	}
	
	
    protected double scoreConfDifferential(FullAStarNode parentNode, int childPos, int childRc) {
    	assertSplitPositions();
    	
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations can make an impact
		
    	// get the full conf, start with the parent first
    	int[] conf = parentNode.getNodeAssignments();
    	
    	// but if this is actually a child node, switch to the child conf
    	if (childPos >= 0) {
    		
    		// parent shouldn't be assigned here
    		assert (conf[childPos] < 0);
    		
    		// make the child conf
    		System.arraycopy(conf, 0, childConf, 0, numPos);
    		childConf[childPos] = childRc;
    		conf = childConf;
    	}
		
    	double ans = scoreNode(conf);
    	
		return ans;
    }
	
    
	protected double scoreNode(int[] partialConf){

		if(traditionalScore) {
			// RCTuple definedTuple = new RCTuple(partialConf);
			rcTuple.set(partialConf);

			//"g-score"
			double score = emat.getConstTerm() + emat.getInternalEnergy(rcTuple);

			//"h-score"
			//score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
			//plus contributions associated with each of the undefined res ("h-score")
			ArrayList<Integer> undefinedLevels = getUndefinedPos(rcTuple);

			for(int level : undefinedLevels) {

				double bestInteractionE = energyLB ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;

				// ArrayList<Integer> rcsAtUndefinedLevel = panPmat.unprunedRCsAtPos(level);
				ArrayList<Integer> rcsAtUndefinedLevel = panUprunedRCsAtPos.get(level);
				for(int rc : rcsAtUndefinedLevel) {
					double rcContribution = RCContribution(level, rc, rcTuple, partialConf, undefinedLevels);

					bestInteractionE = energyLB ? Math.min(bestInteractionE, rcContribution) : 
						Math.max(bestInteractionE, rcContribution);
				}

				score += bestInteractionE;
			}

			return score;
		}

		else {
			//other possibilities include MPLP, etc.
			//But I think these are better used as refinements
			//we may even want multiple-level refinement
			throw new RuntimeException("Advanced A* scoring methods not implemented yet!");
		}
	}


	public double confEnergyBound(int[] partialConf) {
		return scoreNode(partialConf);
	}
}
