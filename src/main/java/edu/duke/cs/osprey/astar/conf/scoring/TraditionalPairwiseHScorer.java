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

package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.tools.MathTools;

public class TraditionalPairwiseHScorer implements AStarScorer {
	
	public final EnergyMatrix emat;
	public final RCs rcs;
	public final MathTools.Optimizer optimizer;
	
	private double[][][] undefinedEnergies; // indexed by (pos1,pos2), rc at pos1
	private ConfAStarNode cachedNode;
	private double[][] cachedEnergies;
	
	public TraditionalPairwiseHScorer(EnergyMatrix emat, RCs rcs) {
		this(emat, rcs, MathTools.Optimizer.Minimize);
	}

	public TraditionalPairwiseHScorer(EnergyMatrix emat, RCs rcs, MathTools.Optimizer optimizer) {
		this.emat = emat;
		this.rcs = rcs;
		this.optimizer = optimizer;
		
		int numPos = emat.getNumPos();
		
		// pre-compute all undefined energy terms
		undefinedEnergies = new double[numPos][][];
		for (int pos1=0; pos1<numPos; pos1++) {
			
			int numRCs = rcs.get(pos1).length;
			undefinedEnergies[pos1] = new double[numRCs][];
			
			for (int i=0; i<numRCs; i++) {
				int rc1 = rcs.get(pos1)[i];
				
				undefinedEnergies[pos1][i] = new double[numPos];
			
				for (int pos2=0; pos2<pos1; pos2++) {
					
					// optimize over rc2
					double optEnergy = optimizer.initDouble();
					for (int rc2 : rcs.get(pos2)) {
						optEnergy = optimizer.opt(optEnergy, emat.getPairwise(pos1, rc1, pos2, rc2));
					}
					
					undefinedEnergies[pos1][i][pos2] = optEnergy;
				}
			}
		}
		
		// allocate space for the cache
		cachedNode = null;
		cachedEnergies = new double[numPos][];
		for (int pos=0; pos<numPos; pos++) {
			cachedEnergies[pos] = new double[rcs.get(pos).length];
		}
	}
	
	public TraditionalPairwiseHScorer make() {
		return new TraditionalPairwiseHScorer(emat, rcs, optimizer);
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {
		
    	// bound energy of undefined conf
    	double hscore = 0;
    	
    	calcCachedEnergies(confIndex, rcs);
    	
		// for each undefined pos...
		for (int i=0; i<confIndex.numUndefined; i++) {
			int pos = confIndex.undefinedPos[i];
			
			// optimize over rcs at this pos
			double optRCEnergy = optimizer.initDouble();
			for (int j=0; j<rcs.get(pos).length; j++) {
				optRCEnergy = optimizer.opt(optRCEnergy, cachedEnergies[pos][j]);
			}
			
			hscore += optRCEnergy;
		}
		
		return hscore;
	}

	@Override
	public double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations can make an impact
		
		// if the intermediate energies aren't cached, calculate them now
		if (cachedNode != confIndex.node) {
			calcCachedEnergies(confIndex, rcs);
			cachedNode = confIndex.node;
		}
		
    	// compute the h-score
    	double hscore = 0;
    	for (int i=0; i<confIndex.numUndefined; i++) {
    		int pos = confIndex.undefinedPos[i];
    		
    		// don't score at nextPos, it's defined now
    		if (pos == nextPos) {
    			continue;
    		}
    		
    		// optimize energy over all rcs
    		double optRCEnergy = optimizer.initDouble();
    		
    		double[] cachedEnergiesAtPos = cachedEnergies[pos];
    		double[][] undefinedEnergiesAtPos = undefinedEnergies[pos];
    		
			// for each rc at this pos...
			int[] rcsAtPos = rcs.get(pos);
			int n = rcsAtPos.length;
			for (int j=0; j<n; j++) {
				int rc = rcsAtPos[j];
				
				double rcEnergy = cachedEnergiesAtPos[j];
				
				// subtract undefined contribution
				if (pos > nextPos) {
					rcEnergy -= undefinedEnergiesAtPos[j][nextPos];
				}
				
				// add defined contribution
				rcEnergy += emat.getPairwise(pos, rc, nextPos, nextRc);
				
				optRCEnergy = optimizer.opt(optRCEnergy, rcEnergy);
			}
			
			hscore += optRCEnergy;
    	}
    	
    	return hscore;
	}

	private void calcCachedEnergies(ConfIndex confIndex, RCs rcs) {
		
		// for each undefined pos...
		for (int i=0; i<confIndex.numUndefined; i++) {
			int pos1 = confIndex.undefinedPos[i];
			
			// for each rc...
			int[] rcs1 = rcs.get(pos1);
			int n1 = rcs1.length;
			for (int j=0; j<n1; j++) {
				int rc1 = rcs1[j];
				
				// start with the one-body energy
				double energy = emat.getOneBody(pos1, rc1);
				
				// add defined energies
				for (int k=0; k<confIndex.numDefined; k++) {
					int pos2 = confIndex.definedPos[k];
					int rc2 = confIndex.definedRCs[k];
					
					energy += emat.getPairwise(pos1, rc1, pos2, rc2);
				}
				
				// add undefined energies
				double[] energies = undefinedEnergies[pos1][j];
				for (int k=0; k<confIndex.numUndefined; k++) {
					int pos2 = confIndex.undefinedPos[k];
					if (pos2 < pos1) {
						energy += energies[pos2];
					}
				}
				
				cachedEnergies[pos1][j] = energy;
			}
		}
	}
}
