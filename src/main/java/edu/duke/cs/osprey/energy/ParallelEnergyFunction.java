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

package edu.duke.cs.osprey.energy;

import java.util.ArrayList;

import edu.duke.cs.osprey.parallelism.WorkCrew;
import edu.duke.cs.osprey.parallelism.Worker;

public class ParallelEnergyFunction implements EnergyFunction {
	
	private static class EnergyWorker extends Worker {
		
		private int startIndex;
		private int stopIndex;
		private double energy;

		public EnergyWorker(WorkCrew<EnergyWorker> crew) {
			super(crew);
		}

		@Override
		protected void workIt() {
			
			// copy some references/values to the stack for that little extra bit of speed
			ArrayList<EnergyFunction> terms = efunc.terms;
			ArrayList<Double> coeffs = efunc.coeffs;
			int startIndex = this.startIndex;
			int stopIndex = this.stopIndex;
			
			// do it!
			double energy = 0;
			for (int i=startIndex; i<=stopIndex; i++) {
				energy += terms.get(i).getEnergy()*coeffs.get(i);
			}
			
			// save results to 'this' memory
			this.energy = energy;
		}
	}

	private static WorkCrew<EnergyWorker> crew;
	private static ParallelEnergyFunction efunc;
	
	static {
		crew = null;
		efunc = null;
	}
	
	public static boolean isCrewStarted() {
		return crew != null;
	}
	
	public static void startCrew(int numThreads) {
		crew = new WorkCrew<>("Energy");
		for (int i=0; i<numThreads; i++) {
			new EnergyWorker(crew);
		}
		crew.start();
	}
	
	public static void stopCrew() {
		if (crew == null) {
			return;
		}
		crew.askToStop();
		crew = null;
	}
	
	public static void setEFunc(ParallelEnergyFunction val) {
		
		if (efunc == val) {
			return;
		}
		
		efunc = val;
		
		// set up partition
		// partition terms among workers
		int numTerms = efunc.terms.size();
		int numWorkers = crew.getWorkers().size();
		int width = (numTerms + numWorkers - 1)/numWorkers;
		int startIndex = 0;
		int stopIndex = startIndex + width - 1;
		for (EnergyWorker worker : crew.getWorkers()) {
			worker.startIndex = startIndex;
			worker.stopIndex = stopIndex;
			startIndex += width;
			stopIndex = Math.min(stopIndex + width, numTerms - 1);
		}
	}

	private static final long serialVersionUID = -2789380428939629566L;
	
	private ArrayList<EnergyFunction> terms;
	private ArrayList<Double> coeffs;
	
	public ParallelEnergyFunction(ArrayList<EnergyFunction> terms, ArrayList<Double> coeffs) {
		this.terms = terms;
		this.coeffs = coeffs;
	}
	
	public ArrayList<EnergyFunction> getTerms() {
		return terms;
	}

	public ArrayList<Double> getCoeffs() {
		return coeffs;
	}
	
	@Override
	public double getEnergy() {
		
		setEFunc(this);
		
		// start all the processors
		crew.sendWork();
		
		// wait for the work to finish
		try {
			// wait 10 seconds at first
			boolean finished = crew.waitForResults(10000);
			if (!finished) {
				System.err.println("WARNING: ParallelEnergyFunction is taking more than 10 seconds to evaluate. Maybe something is wrong?");
				
				// wait 10000 more seconds
				//sometimes there are temporary slowdowns on the cluster and runs crash unnecessarily if we have a hard cap around 10 s
				finished = crew.waitForResults(10000000);
				if (!finished) {
					throw new Error("Timed out waiting 10000 seconds for energy calculations to finish!"
							+ "\nEnergy calculation shouldn't take more than 10000 seconds, right?");
				}
			}
		} catch (InterruptedException ex) {
			// something wanted us to stop, so stop, then forward the exception
			throw new Error(ex);
		}
		
		// collect the results
		double energy = 0;
		for (EnergyWorker worker : crew) {
			energy += worker.energy;
		}
		return energy;
	}
}
