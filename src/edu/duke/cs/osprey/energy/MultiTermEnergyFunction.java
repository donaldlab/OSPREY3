/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy;

import java.util.ArrayList;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.structure.Residue;

/**
 *
 * @author mhall44
 */
//This an energy function consisting of terms that are other energy functions
//For example, it can be AMBER energy + dihedral energy
//or ( energy of a pair of residues ) - (energy of residue 1) - (energy of residue 2)
//Total energy = sum_i (coeff[i] * energy i)

public class MultiTermEnergyFunction implements EnergyFunction {

	private static final long serialVersionUID = -3516267414126293331L;

	private static int NUM_THREADS = 1;
	public static boolean useParallelEFunc = true;

	ArrayList<EnergyFunction> terms = new ArrayList<>();
	ArrayList<Double> coeffs = new ArrayList<>();
	ArrayList<Double> partialE = new ArrayList<>();
	ArrayList<Integer> indexes = new ArrayList<>();
	double preCompE = 0.0;

	private ParallelEnergyFunction parallelEfunc = null;


	//Constructor with no terms
	public MultiTermEnergyFunction(){

	}

	public void addTerm(EnergyFunction ef){
		terms.add(ef);
		coeffs.add(1.);
		partialE.add(0.0);
		indexes.add(indexes.size());
	}

	//add a term to this with a given coefficient
	public void addTermWithCoeff(EnergyFunction ef, double coeff){
		terms.add(ef);
		coeffs.add(coeff);
		partialE.add(0.0);
		indexes.add(indexes.size());
	}


	public static void setNumThreads( int threads ) { 
		NUM_THREADS = threads;

		if(NUM_THREADS < 1) NUM_THREADS = 1;

		else if(NUM_THREADS > Runtime.getRuntime().availableProcessors()) 
			NUM_THREADS = Runtime.getRuntime().availableProcessors();
	}


	public static int getNumThreads() {
		return NUM_THREADS;
	}


	@Override
	public double getEnergy(){

		double E = 0;

		if(terms.size()!=coeffs.size()){
			throw new RuntimeException("ERROR: MultiTermEnergyFunction has "+terms.size()
			+" terms but "+coeffs.size()+" coefficients");
		}

		/* OPTIMIZATION: It might seem like we could optimize more here by only
        synchronizing with the thread pool if there are enough terms to make paying
        the syncronization cost worth it. ie, only call the threads if numTerms > c.
        Counterintuitively, this actually performs worse in practice! My tests indicated
        we paid about a 20% performance penalty for forcing the main thread to do the
        energy calculations for even just one energy term. I suspect this has to do
        with cpu cache performance. Ignoring inter-core thread migration issues and thread
        pinning, I suspect that the main thread spends most of its time on one core, and
        the processor threads on other cores in numTerms<numProcessors cases. That means
        the cpu cache on the processing cores is already "warmed up" for energy calculations.
        Forcing the main thread to do that work (and alternate between work in CCD) is probably
        causing many cpu cache misses and slowing down overall performance.
        But that's just a theory. Either way, just use parallelism all the time even for 1 or
        2 energy term because emperically it works much better! =)
		 */        
		if(NUM_THREADS == 1) {
			for(int termNum=0; termNum<terms.size(); termNum++){
				double termE = terms.get(termNum).getEnergy();
				E += coeffs.get(termNum)*termE;
			}
		} else {

			if(useParallelEFunc) {
				if (!ParallelEnergyFunction.isCrewStarted()) {
					ParallelEnergyFunction.startCrew(NUM_THREADS);
				}
				if (parallelEfunc == null) {
					parallelEfunc = new ParallelEnergyFunction(terms, coeffs);
				}
				E = parallelEfunc.getEnergy();
			}
			
			// AAO 2016: can allow concurrent minimizers, which is 25x faster for systems 
			// that cannot saturate the threadpool with a single minimizer object. only my code actually enters this block
                        //MH: setting up a config option in Main to choose this
			else {
				indexes.parallelStream().forEach((term) -> partialE.set(term, terms.get(term).getEnergy()*coeffs.get(term)));
				for(int term = 0; term < indexes.size(); ++term) E += partialE.get(term);
			}
		}

		if(Double.isNaN(E) || Double.isInfinite(E))//This can happen if there are positive and negative terms
			//with infinite energy...we assume this to be an impossible conformation
			//and thus return inifinity
			return Double.POSITIVE_INFINITY;

		return preCompE = E;
	}

	public ArrayList<EnergyFunction> getTerms() {
		return terms;
	}

	public ArrayList<Double> getCoeffs() {
		return coeffs;
	}

	public ArrayList<EnergyFunction> makeDOFPartialEFuncs(ArrayList<DegreeOfFreedom> dofs) {
		ArrayList<EnergyFunction> dofEfuncs = new ArrayList<>();
		for (DegreeOfFreedom dof : dofs) {

			if(dof.getResidue() == null){//degree of freedom moves multiple residues
				dofEfuncs.add(this);//use the full energy
			}
			else {//degree of freedom moves just one residue
				dofEfuncs.add(makeResidueEfunc(dof.getResidue()));
			}
		}
		return dofEfuncs;
	}

	private EnergyFunction makeResidueEfunc(Residue residue) {

		// find all the terms that involve this residue
		MultiTermEnergyFunction resEfunc = new MultiTermEnergyFunction();
		for (EnergyFunction term : terms) {

			// TODO: could make the energy term worry about residue association queries
			// instead of switching on type here
			if (term instanceof SingleResEnergy) {
				SingleResEnergy singleResTerm = (SingleResEnergy)term;
				if (singleResTerm.getRes() == residue) {
					resEfunc.addTerm(singleResTerm);
				}
			} else if (term instanceof ResPairEnergy) {
				ResPairEnergy resPairTerm = (ResPairEnergy)term;
				if (resPairTerm.getRes1() == residue || resPairTerm.getRes2() == residue) {
					resEfunc.addTerm(resPairTerm);
				}
			} else {//not a one-body or pairwise forcefield term--may involve all residues
				resEfunc.addTerm(term);
				//throw new Error("Unsupported energy function term: " + term.getClass().getName());
			}
		}
		return resEfunc;
	}

	public double getPreCompE() {
		return preCompE;
	}

	public void setPreCompE(double in) {
		preCompE = in;
	}
}

