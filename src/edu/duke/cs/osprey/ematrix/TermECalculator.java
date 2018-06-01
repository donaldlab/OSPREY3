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

package edu.duke.cs.osprey.ematrix;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.ematrix.epic.EPICFitter;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.ematrix.epic.EPoly;
import edu.duke.cs.osprey.ematrix.epic.FitParams;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.handlempi.MPISlaveTask;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 *
 * @author mhall44
 */

//This class is meant to precalculate the minimized energy of an energy term
//(e.g. intra+shell or pairwise)
//it should be suitable for calculating at a slave node

public class TermECalculator implements MPISlaveTask {
    
    ConfSpace confSpace;
    boolean doingEPIC;//doing EPIC fit instead of just minimum computation
    boolean doingIntra;//doing just intra energy (only relevant for one-body energies)
    //if false, one-body energies will be intra+shell
    
    int res[];//the residues to handle (can be either 1 or 2)
    
    EnergyFunction termE;//evaluates the energy for this term
    
    PruningMatrix pruneMat;//pruning matrix so EPIC doesn't handle pruned RCs or pairs; this can be null 
    //for scalar (non-EPIC) energy calculations
    
    EPICSettings epicSettings = null;//needed if doing EPIC
    
    boolean addResEntropy;//add residue entropy to one-body energies
    
    //We'll be calculating either one-body or pairwise energies
    //as either scalar or EPIC energies
    //So we're trying to fill in one of these four lists
    ArrayList<Double> oneBodyE = new ArrayList<>();
    ArrayList<ArrayList<Double>> pairwiseE = new ArrayList<>();
    ArrayList<EPoly> oneBodyPoly = new ArrayList<>();
    ArrayList<ArrayList<EPoly>> pairwisePoly = new ArrayList<>();
    HashMap<ArrayList<Integer>, Double> nBodyE = new HashMap<>();
    ArrayList<MoleculeModifierAndScorer> mofs = null;
    
    
    private static final boolean SAPEKeepStandalone = false;
    //Make SAPE terms in the EPIC matrix retain ability to be evaluated standalone
    //getting rid of this is fine for normal operation as of Aug '16 and could save memory
    
    
    public TermECalculator(ConfSpace s, ArrayList<Residue> shellResidues, 
            boolean doEPIC, boolean doIntra, PruningMatrix prm, EPICSettings es, 
            boolean addResEnt, int... resToCalc){
        
        confSpace = s;
        doingEPIC = doEPIC;
        doingIntra = doIntra;
        res = resToCalc;
        pruneMat = prm;
        epicSettings = es;
        addResEntropy = addResEnt;
        
        Residue firstRes = confSpace.posFlex.get(res[0]).res;
        
        if(doingIntra)//intra energy only
            termE = EnvironmentVars.curEFcnGenerator.singleResEnergy(firstRes);
        else{
            
            if(res.length==1){//intra+shell
                termE =EnvironmentVars.curEFcnGenerator.intraAndShellEnergy(firstRes, shellResidues);
            }
            else if(res.length==2){//pairwise
                Residue secondRes = confSpace.posFlex.get(res[1]).res;
                termE = EnvironmentVars.curEFcnGenerator.resPairEnergy(firstRes, secondRes);
            }
            else if(res.length > 2){
				// add all respair energies
				termE = new MultiTermEnergyFunction();
				
				for (int i = 0; i < res.length; i++) {
					Residue l = confSpace.posFlex.get(res[i]).res;
					for (int j = i+1; j < res.length; j++) {
						Residue r = confSpace.posFlex.get(res[j]).res;
						((MultiTermEnergyFunction)termE).addTerm(EnvironmentVars.curEFcnGenerator.resPairEnergy(l, r));
					}
				}
			}
            else{
                throw new UnsupportedOperationException("ERROR: Can only precompute energy for >= 1 body terms");
                //we are excluding shell-shell interactions throughout this version of OSPREY,
                //since they don't change and are time-consuming
            }
        }
    }
    

    @Override
    public Object doCalculation() {
        //This operation, which can be handled on a slave node once this TermMinECalculator is sent there,
        //calculates the specified energies
        
        if(res.length==1){//1-body calculation
            
            oneBodyCalc();
            
            if(doingEPIC)
                return oneBodyPoly;
            else
                return oneBodyE;
        }   
        else if(res.length==2){//pairwise calculation
            
            pairwiseCalc();
            
            if(doingEPIC)
                return pairwisePoly;
            else
                return pairwiseE;
        }
        else if(res.length > 2) {
			
			nBodyCalc();

			if(doingEPIC)
				throw new UnsupportedOperationException("ERROR: have not implemented >2-body energies with EPIC");
			else
				return nBodyE;
		}
        else
            throw new UnsupportedOperationException("ERROR: Trying to precompute term for "+res.length+" bodies");
    }
    
    
    public void oneBodyCalc(){
        //list minimized one-body energies for all the RCs in res
        //the energy minimized can be just intra or intra+shell (decided in constructor)
        
        ArrayList<RC> RCList = confSpace.posFlex.get(res[0]).RCs;
        
        for(int RCNum=0; RCNum<RCList.size(); RCNum++){
            RCTuple RCTup = new RCTuple(res[0],RCNum);
            calcTupleEnergy(RCTup);
        }
    }
    
    
    public void pairwiseCalc(){
        //list minimized one-body energies for all the RCs in res
        
        ArrayList<RC> RCList1 = confSpace.posFlex.get(res[0]).RCs;
        ArrayList<RC> RCList2 = confSpace.posFlex.get(res[1]).RCs;
        
        for(int firstRCNum=0; firstRCNum<RCList1.size(); firstRCNum++){
            
            for(int secondRCNum=0; secondRCNum<RCList2.size(); secondRCNum++){
                RCTuple RCTup = new RCTuple(res[0],firstRCNum,res[1],secondRCNum);
                calcTupleEnergy(RCTup);
            }
        }
    }
    
    
	public void nBodyCalc() {
		if(pruneMat == null)
			throw new RuntimeException("ERROR: n-body calc requires a pruning matrix");
		
		ArrayList<ArrayList<RC>> RCLists = new ArrayList<>(res.length);
		for(int index = 0; index < res.length; ++index) RCLists.add(new ArrayList<>());
		
		for(int i = 0; i < res.length; ++i) {
			for(Integer rcNum : pruneMat.unprunedRCsAtPos(res[i]))
				RCLists.get(i).add(confSpace.posFlex.get(res[i]).RCs.get(rcNum));
		}

		Integer[] resElements = new Integer[res.length]; Arrays.fill(resElements, -1);
		Integer[] rcNums = new Integer[res.length]; Arrays.fill(rcNums, -1);
		createNBodyTuples(RCLists, resElements, rcNums, 0);
	}


	private void createNBodyTuples( ArrayList<ArrayList<RC>> RCLists, Integer[] resElements, Integer[] rcNums, int depth ) {

		if(depth == resElements.length) {
			// resIndex and rcNums are fully populated
			RCTuple nTuple = new RCTuple(new ArrayList<>(Arrays.asList(resElements)), new ArrayList<>(Arrays.asList(rcNums)));
			calcTupleEnergy(nTuple);
			return;
		}

		if(resElements[depth] == -1) 
			resElements[depth] = new Integer(res[depth]);

		for(int rcNum = 0; rcNum < RCLists.get(depth).size(); ++rcNum) {
			rcNums[depth] = rcNum;
			createNBodyTuples(RCLists, resElements, rcNums, depth+1);
		}
	}
	
	
    public void calcTupleEnergyLazy(RCTuple RCs) {
        // supports parallel version of the function
    	mofs = new ArrayList<>();
    	
        boolean skipTuple = false;
        
        if(RCs.pos.size()==2){//pair: need to check for parametric incompatibility
            //If there are DOFs spanning multiple residues, then parametric incompatibility
            //is whether the pair is mathematically possible (i.e. has a well-defined voxel)
            RC rc1 = confSpace.posFlex.get( RCs.pos.get(0) ).RCs.get( RCs.RCs.get(0) );
            RC rc2 = confSpace.posFlex.get( RCs.pos.get(1) ).RCs.get( RCs.RCs.get(1) );
            if(rc1.isParametricallyIncompatibleWith(rc2)){
                skipTuple = true;
            }
        }
        
		else if(RCs.pos.size() > 2){
			ArrayList<RC> indivRCs = new ArrayList<>(RCs.pos.size());
			for(int i = 0; i < RCs.pos.size(); ++i) indivRCs.add(confSpace.posFlex.get(RCs.pos.get(i)).RCs.get(RCs.RCs.get(i)));

			for(int i = 0; i < indivRCs.size(); ++i) {
				RC rc1 = indivRCs.get(i);

				for(int j = i+1; j < indivRCs.size(); ++j) {
					RC rc2 = indivRCs.get(j);

					if(rc1.isParametricallyIncompatibleWith(rc2)) {
						skipTuple = true;
						break;
					}
				}

				if(skipTuple) break;
			}
		}
        
        if( (pruneMat!=null) && (!skipTuple) ){
            if(pruneMat.isPruned(RCs))
                skipTuple = true;
        }
        
        if(!skipTuple) {
            MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(termE,confSpace,RCs);
            mofs.add(mof);

            DoubleMatrix1D bestDOFVals;

            if(mof.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
                CCDMinimizer ccdMin = new CCDMinimizer(mof,true);
                bestDOFVals = ccdMin.minimize().dofValues;
            }
            else//molecule is already in the right, rigid conformation
                bestDOFVals = DoubleFactory1D.dense.make(0);

            double minEnergy = mof.getValue(bestDOFVals);

            if(doingEPIC){
               throw new UnsupportedOperationException("ERROR: EPIC is not supprted in this function");
            }
        }
    }
    
    
    public void calcTupleEnergy(RCTuple RCs){
        //calculate the rigid energy, minimum energy or EPIC fit for an RC tuple
        //(depending on what type of matrix is being calculated)
        //store it in our list of results

        boolean skipTuple = false;
        double minEnergy = Double.POSITIVE_INFINITY;//rigid or voxel-minimum energy
        EPoly EPICFit = null;//can use null poly for pruned term
        
        if(RCs.pos.size()==2){//pair: need to check for parametric incompatibility
            //If there are DOFs spanning multiple residues, then parametric incompatibility
            //is whether the pair is mathematically possible (i.e. has a well-defined voxel)
            RC rc1 = confSpace.posFlex.get( RCs.pos.get(0) ).RCs.get( RCs.RCs.get(0) );
            RC rc2 = confSpace.posFlex.get( RCs.pos.get(1) ).RCs.get( RCs.RCs.get(1) );
            if(rc1.isParametricallyIncompatibleWith(rc2)){
                skipTuple = true;
            }
        }
        
		else if(RCs.pos.size() > 2){
			ArrayList<RC> indivRCs = new ArrayList<>(RCs.pos.size());
			for(int i = 0; i < RCs.pos.size(); ++i) indivRCs.add(confSpace.posFlex.get(RCs.pos.get(i)).RCs.get(RCs.RCs.get(i)));

			for(int i = 0; i < indivRCs.size(); ++i) {
				RC rc1 = indivRCs.get(i);

				for(int j = i+1; j < indivRCs.size(); ++j) {
					RC rc2 = indivRCs.get(j);

					if(rc1.isParametricallyIncompatibleWith(rc2)) {
						skipTuple = true;
						break;
					}
				}

				if(skipTuple) break;
			}
		}
        
        if( (pruneMat!=null) && (!skipTuple) ){
            if(pruneMat.isPruned(RCs))
                skipTuple = true;
        }
        
        if(!skipTuple){
            MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(termE,confSpace,RCs);

            DoubleMatrix1D bestDOFVals;

            if(mof.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
                CCDMinimizer ccdMin = new CCDMinimizer(mof,true);
                bestDOFVals = ccdMin.minimize().dofValues;
            }
            else//molecule is already in the right, rigid conformation
                bestDOFVals = DoubleFactory1D.dense.make(0);

            minEnergy = mof.getValue(bestDOFVals);

            if(doingEPIC){
                EPICFit = compEPICFit(mof,minEnergy,bestDOFVals,RCs);
            }
        }
        
        //now store the energy we calculated in the appropriate list of results
        int numBodies = RCs.pos.size();
        
        if ( numBodies == 1 ) {//one-body term
            
            if(doingEPIC)
                oneBodyPoly.add(EPICFit);
            else{
                if(addResEntropy)
                    minEnergy += getResEntropy(RCs);
                
                oneBodyE.add(minEnergy);
            }
        }
        else if ( numBodies == 2 ) {//pairwise term
            
            int firstRCNum = RCs.RCs.get(0);
            //we have a list of lists, so we use firstRCNum to see
            //which list to add to
            
            if(doingEPIC){
                if(pairwisePoly.size()<=firstRCNum)//this is the first term for RC #firstRCNum
                    pairwisePoly.add(new ArrayList<EPoly>());
                
                pairwisePoly.get(firstRCNum).add(EPICFit);
            }
            else {
                if(pairwiseE.size()<=firstRCNum)
                    pairwiseE.add(new ArrayList<Double>());
                
                pairwiseE.get(firstRCNum).add(minEnergy);
            }
        }
        else if( numBodies > 2 ) {
			if(doingEPIC)
				throw new RuntimeException("ERROR: have not implemented > 2-body terms with EPIC");
			
			// store rc and its energy
			if(!skipTuple)
				nBodyE.put(RCs.RCs, minEnergy);
		}
        else
            throw new UnsupportedOperationException("ERROR: Trying to precompute term for "+numBodies+" bodies");
    }
    
    
    double getResEntropy(RCTuple RCs){
        //Given a singleton RC tuple, return the residue entropy for its amino acid type
        if( RCs.pos.size() != 1 ){
            throw new RuntimeException("ERROR: Trying to get res entropy for non-singleton RC"
                    + " tuple "+RCs.stringListing());
        }
        
        double resEntropy = confSpace.getRCResEntropy( RCs.pos.get(0), RCs.RCs.get(0) );
        return resEntropy;
    }
    
    
    private EPoly compEPICFit ( MoleculeModifierAndScorer mof, double minEnergy, 
            DoubleMatrix1D bestDOFVals, RCTuple RCList ){
        
            //Do the EPIC series fits for this rotamer pair
            //we are given an energy objective function valid on the voxel,
            //plus the energy and DOF values for the minimum-energy point

            EPICFitter fitter = new EPICFitter(mof,epicSettings,bestDOFVals,minEnergy);
            
            EPoly bestFit = null;
            
            if(mof.getNumDOFs()>0){//need to actually compute polynomial fits
                
                FitParams curFitParams = FitParams.quadratic(mof.getNumDOFs(),false);
                double bestResid = Double.POSITIVE_INFINITY;
                int bestBound = -1;
                ArrayList<EPoly> series = new ArrayList<>();//list of fit series
                  
                for( int boundCount=0; bestResid>epicSettings.EPICGoalResid && curFitParams!=null; boundCount++ ){
                    
                    System.out.println("FIT NUMBER "+boundCount+": "+curFitParams.getDescription());
                    
                    EPoly curSeries;
                    try{
                        curSeries = fitter.doFit(curFitParams);
                    }
                    catch(Exception e){//sometimes singular matrices will arise during fitting
                        //if so we skip that order.  More SAPE etc. might help
                        System.err.println("Fit failed: "+e.getMessage());
                        e.printStackTrace();
                        series.add(null);
                        
                        //go up to next order
                        curFitParams = fitter.raiseFitOrder(curFitParams);
                        continue;
                    }
                    
                    double meanResid = fitter.crossValidateSeries(curSeries,curFitParams);
                    series.add(curSeries);
                    
                    if( curFitParams.order==2 && curFitParams.PCOrder<=2 && 
                            curFitParams.SAPECutoff==0 )//pure quadratic--to use as template for PCs
                        fitter.PCTemplate = curSeries;
                    
                    if(meanResid<bestResid){
                        if(checkStaysPositive(curSeries, mof)){
                            bestResid = meanResid;
                            bestBound = boundCount;
                        }
                    }
                    
                    curFitParams = fitter.raiseFitOrder(curFitParams);
                }
                
                if(bestBound==-1){
                    throw new RuntimeException("ERROR: No EPIC fit found without serious errors"
                            + " (e.g. significant error at minimum)");
                }
                
                if(bestResid > epicSettings.EPICGoalResid){
                    System.out.println("Failed to reach goal residual.");
                }
                System.out.println("Best residual: "+bestResid+" for bound number "+bestBound);

                printFitTests(fitter, RCList, minEnergy, mof, bestDOFVals, series);

                bestFit = series.get(bestBound);
            }
            else//no DOFs
                bestFit = fitter.blank();

            bestFit.setMinE(minEnergy);
            
            if(!SAPEKeepStandalone)
                bestFit.deleteMOFStandalone();
            
            return bestFit;
    }
    
    
    private void printFitTests(EPICFitter fitter, RCTuple RCList, double minEnergy,
            MoleculeModifierAndScorer mof, DoubleMatrix1D bestDOFVals, ArrayList<EPoly> series){
        //Do some tests on fit performance, and print the results
        int numDOFs = fitter.numDOFs;

        //TESTING FITS
        System.out.println("RCs: "+RCList.stringListing());
        System.out.println("Minimum energy: "+minEnergy);

        double testScales[] = new double[] { 0.01, 0.5, 5, 100 };//100
        int samplesPerScale = 3;



        double relMax[] = new double[numDOFs];//maximum shifts of degrees of freedom relative to minimum point (startVec)
        double relMin[] = new double[numDOFs];
        DoubleMatrix1D constr[] = mof.getConstraints();
        for(int dof=0; dof<numDOFs; dof++){
            relMax[dof] = constr[1].get(dof) - bestDOFVals.get(dof);
            relMin[dof] = constr[0].get(dof) - bestDOFVals.get(dof);
        }


        for(double scale : testScales){
            for(int s=0; s<samplesPerScale; s++){

                //Generate vector relative to minimum
                double dx[] = new double[numDOFs];
                //and absolute
                DoubleMatrix1D sampAbs = DoubleFactory1D.dense.make(numDOFs);
                for(int dof=0; dof<numDOFs; dof++){
                    double top = Math.min(relMax[dof], scale);
                    double bottom = Math.max(relMin[dof], -scale);
                    dx[dof] = bottom + Math.random()*(top-bottom);
                    sampAbs.set(dof, bestDOFVals.get(dof)+dx[dof]);
                }

                double trueVal = mof.getValue(sampAbs) - minEnergy;

                System.out.print("TEST: scale="+scale+" dx=");
                for(int dof=0; dof<numDOFs; dof++)
                    System.out.print(dx[dof]+" ");

                System.out.print("TRUE="+trueVal+" FIT=");

                for(EPoly b : series){
                    if(b!=null)
                        System.out.print(b.evaluate(sampAbs,false,false)+" ");
                }

                System.out.println();
            }
        }
    }
    
    
    
    private boolean checkStaysPositive(EPoly term, MoleculeModifierAndScorer mof){
        //Check, by minimization, that this EPIC term doesn't go too far negative
        //mof defines the voxel for the term conveniently
        
        ArrayList<EPoly> termList = new ArrayList<>();
        termList.add(term);
        EPICEnergyFunction epicEF = new EPICEnergyFunction(termList, false);
        
        MoleculeModifierAndScorer ofEPIC = new MoleculeModifierAndScorer( epicEF, mof.getConstraints(),
                mof.getMolec(), mof.getDOFs() );
        
        CCDMinimizer emin = new CCDMinimizer(ofEPIC, true);
        DoubleMatrix1D lowPoint = emin.minimize().dofValues;
        
        //energies will be relative to the voxel center energy (from the initial minimization)
        double lowPointTrueE = mof.getValue(lowPoint) - term.getMinE();
        double lowPointEPICE = ofEPIC.getValue(lowPoint);
        
        
        double tol = 0.1;//we'll only consider this a problem if the minimum
        //(where EPIC should be pretty accurate) is well below the actual minimum
        if(lowPointEPICE < lowPointTrueE - tol){
            System.out.println("Rejecting fit because of low minimum for EPIC term, "
                    +lowPointEPICE+".  Corresponding true E: "+lowPointTrueE);
            return false;
        }
        
        //But let's warn if the minimum is low in absolute terms
        if(lowPointEPICE < -tol){
            System.out.println("Warning: low minimum for EPIC term, "
                    +lowPointEPICE+".  Corresponding true E: "+lowPointTrueE);
        }
        
        return true;
    }
    
}
