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

package edu.duke.cs.osprey.kstar.pfunc.impl;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.kstar.KSConf;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.kstar.pfunc.PFAbstract;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

public class PFAdapter extends PFAbstract {
	
	private static final long serialVersionUID = -4339776095146902581L;
	
	private String pfImpl;
	private PartitionFunction pfunc;

	public PFAdapter(String pfImpl, int strand, ArrayList<String> sequence, ArrayList<Integer> absolutePos, String checkPointPath, String reducedSPName, KSConfigFileParser cfp, KSSearchProblem panSP) {
		super(strand, sequence, absolutePos, checkPointPath, reducedSPName, cfp, panSP);
		this.pfImpl = pfImpl;
	}
	
	public PartitionFunction getPartitionFunction() {
		return pfunc;
	}
	
	public void setPartitionFunction(PartitionFunction pfunc) {
		this.pfunc = pfunc;
	}

	@Override
	protected void printHeader() {
		// nothing to do
	}

	@Override
	public void start() {
		
		setRunState(RunState.STARTED);
		pfunc.init(null, null, PFAbstract.targetEpsilon);
		
		// report top confs if needed
		if (isFullyDefined() && saveTopConfsAsPDB) {
			pfunc.setConfListener((ScoredConf sconf) -> {	
				double energy = sconf instanceof EnergiedConf ? ((EnergiedConf)sconf).getEnergy() : sconf.getScore();
				saveTopConf(new KSConf(
					((EnergiedConf)sconf).getAssignments(),
					((EnergiedConf)sconf).getScore(),
					energy
				));
			});
		}
	}

	@Override
	protected void computeSlice() {
		iterate();
	}

	@Override
	protected void compute() {
		pfunc.compute();
		update();
	}

	@Override
	protected void iterate() {
		pfunc.compute(pfunc.getParallelism());
		update();
	}

	@Override
	public String getImpl() {
		return pfImpl;
	}
	
	private void update() {
		
		// convert the PartitionFunction status to a PFAbstract status
		this.eAppx = mapStatus(pfunc.getStatus());
		
		// copy the pfunc values
		PartitionFunction.Values values = pfunc.getValues();
		this.qStar = values.qstar;
		this.qPrime = values.qprime;
		this.pStar = values.pstar;
		this.effectiveEpsilon = values.getEffectiveEpsilon();
	}

	private PFAbstract.EApproxReached mapStatus(PartitionFunction.Status status) {
		switch (status) {
			case Estimating: return EApproxReached.FALSE;
			case Estimated: return EApproxReached.TRUE;
			default: return EApproxReached.NOT_POSSIBLE;
		}
	}
}
