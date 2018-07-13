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

package edu.duke.cs.osprey.kstar.pfunc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.GMECConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.gmec.MinimizingConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFAdapter;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFParallel0;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFParallel1;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFParallel2;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFTraditional;
import edu.duke.cs.osprey.kstar.pfunc.impl.PFUB;
import edu.duke.cs.osprey.parallelism.Parallelism;


/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 * Chooses the type of partition function implementation, depending upon user
 * supplied value
 *
 * Deprecated, use GradientDescentPfunc instead
 *
 * These implementations are no longer tested or regularly used, so they
 * may no longer work correctly. Use only with extreme skepticism.
 */
@Deprecated
public class PFFactory {

	public static PFAbstract getPartitionFunction( String implementation, int strand, 
			ArrayList<String> sequence, ArrayList<Integer> absolutePos, String checkPointPath, 
			String searchProblemName, KSConfigFileParser cfp, KSSearchProblem sp ) {

		switch( implementation.toLowerCase() ) {

		case "traditional":
			return new PFTraditional( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );
		
		case "ub":
			return new PFUB( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );

		case "parallel0":
			return new PFParallel0( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );
			
		case "parallel1":
			return new PFParallel1( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );

		case "parallel2":
			return new PFParallel2( strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp );
			
		case "parallelconf":
			return new PFAdapter(implementation, strand, sequence, absolutePos, checkPointPath, searchProblemName, cfp, sp) {
				
				private static final long serialVersionUID = 8718298627856880018L;
				
				private StrandInfo strandInfo;
				
				{
					// get the two search problems
					// TODO: there really should be a simpler way to do this...
					// all that's different is the allowed AAs and the pruning matrix, right?
					// except it's actually reducedMat... so it doesn't actually replace the pruning matrix
					KSSearchProblem multiSeqSearch = sp;
					KSSearchProblem singleSeqSearch = getReducedSearchProblem();
					
					// just in case...
					assert (singleSeqSearch.emat == multiSeqSearch.emat);
					assert (singleSeqSearch.confSpace == multiSeqSearch.confSpace);
					assert (singleSeqSearch.shellResidues == multiSeqSearch.shellResidues);
					
					// get once-per-stand things
					StrandInfo strandInfo;
					if (strandInfoCache != null) {
						strandInfo = strandInfoCache.get(multiSeqSearch.confSpace);
						if (strandInfo == null) {
							strandInfo = new StrandInfo(cfp, multiSeqSearch);
							strandInfoCache.put(multiSeqSearch.confSpace, strandInfo);
						}
					} else {
						strandInfo = new StrandInfo(cfp, multiSeqSearch);
						this.strandInfo = strandInfo;
					}
					
					// make the conf search (eg A*)
					ConfSearchFactory confSearchFactory = ConfSearchFactory.Tools.makeFromConfig(singleSeqSearch, cfp);
					
					// make the partition function
					ParallelConfPartitionFunction pfunc = new ParallelConfPartitionFunction(
						singleSeqSearch.emat,
						singleSeqSearch.reducedMat,
						confSearchFactory,
						strandInfo.ecalc
					);
					pfunc.setReportProgress(!PFAbstract.suppressOutput);
					setPartitionFunction(pfunc);
				}

				@Override
				public void cleanup() {
					super.cleanup();
					
					setPartitionFunction(null);
					
					if (strandInfo != null) {
						strandInfo.cleanup();
						strandInfo = null;
					}
				}
			};
		
		default:
			throw new RuntimeException("ERROR: specified value of parameter kStarPFuncMethod is invalid: " + implementation);

		}
	}
	
	// we can re-use minimizers for the same strand to save overhead on each sequence
	private static Map<ConfSpace,StrandInfo> strandInfoCache = null;
	
	private static class StrandInfo {
		
		public final MinimizingConfEnergyCalculator ecalc;
		
		public StrandInfo(KSConfigFileParser cfp, KSSearchProblem multiSeqSearch) {
			Parallelism parallelism = Parallelism.makeFromConfig(cfp);
			ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
			ecalc = MinimizingConfEnergyCalculator.make(ffparams, multiSeqSearch, parallelism);
		}
		
		public void cleanup() {
			ecalc.clean();
		}
	}
	
	public static void initStrandInfoCache() {
		if (strandInfoCache != null) {
			cleanupStrandInfoCache();
		}
		strandInfoCache = new HashMap<>();
	}
	
	public static void cleanupStrandInfoCache() {
		if (strandInfoCache != null) {
			for (StrandInfo strandInfo : strandInfoCache.values()) {
				strandInfo.cleanup();
			}
			strandInfoCache = null;
		}
	}
}
