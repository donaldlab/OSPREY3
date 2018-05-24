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
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.ematrix.epic.EPICFitter;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.ematrix.epic.EPoly;
import edu.duke.cs.osprey.ematrix.epic.FitParams;
import edu.duke.cs.osprey.ematrix.epic.NewEPICMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;

/**
 *
 * 
 * @author mhall44
 */
public class NewEPICMatrixCalculator {
    
    SimpleConfSpace searchSpace;
    //ArrayList<Residue> shellResidues;
    ConfEnergyCalculator confECalc;
    PruningMatrix pruneMat;
    EPICSettings epicSettings;
    
    //we'll allocate this and fill it in
    private NewEPICMatrix epicMat = null;
    
    static boolean SAPEKeepStandalone = false;
    
    
    public NewEPICMatrixCalculator(SimpleConfSpace confSpace, ConfEnergyCalculator confECalc, 
            PruningMatrix pruneMat, EPICSettings epicSettings){
        searchSpace = confSpace;
        this.confECalc = confECalc;
        this.pruneMat = pruneMat;
        this.epicSettings = epicSettings;
    }
    
    
    
    
    
   //Calculate a pairwise EPIC matrix based on a pairwise energy function
   //since these are not bounds no need to do any fancy partitioning, just do
   //intra+shell and pairwise terms
   public void calcPEM(){
       
       System.out.println();
       System.out.println("BEGINNING EPIC MATRIX PRECOMPUTATION");
       System.out.println();
       
       initMatrix();
       
       for(int pos=0; pos<searchSpace.getNumPos(); pos++){
            
            System.out.println("Starting intra+shell energy calculations for residue "+pos);
            
            for(int rc=0; rc<searchSpace.getNumResConfs(pos); rc++){
                EPoly singlePoly = makeEPoly(new RCTuple(pos,rc));
                epicMat.setOneBody(pos, rc, singlePoly);
            }

            for(int pos2=0; pos2<pos; pos2++){
                System.out.println("Starting pairwise energy calculations for residues "+pos+", "+pos2);

                for(int rc=0; rc<searchSpace.getNumResConfs(pos); rc++){
                    for(int rc2=0; rc2<searchSpace.getNumResConfs(pos2); rc2++){
                        
                        EPoly pairPoly = makeEPoly(new RCTuple(pos, rc, pos2, rc2));
                        epicMat.setPairwise(pos, rc, pos2, rc2, pairPoly);
                    }
                }
            }
        }
       
       System.out.println("EPIC MATRIX CALCULATION DONE");
   }
    
    
    
    
    private void initMatrix(){
        //initialize the matrix we're calculating
        epicMat = new NewEPICMatrix(searchSpace, pruneMat.getPruningInterval());
    }
    
    
    public NewEPICMatrix getEPICMatrix(){
        if(epicMat==null)
            throw new RuntimeException("ERROR: EPIC matrix is null after calculation");
        
        return epicMat;
    }
    
    
    public EPoly makeEPoly(RCTuple RCs){
        //calculate the EPIC fit for an RC tuple

        if(RCs.pos.size()==2){//pair: need to check for parametric incompatibility
            //If there are DOFs spanning multiple residues, then parametric incompatibility
            //is whether the pair is mathematically possible (i.e. has a well-defined voxel)
            ResidueConf rc1 = searchSpace.positions.get( RCs.pos.get(0) ).resConfs.get( RCs.RCs.get(0) );
            ResidueConf rc2 = searchSpace.positions.get( RCs.pos.get(1) ).resConfs.get( RCs.RCs.get(1) );
            if(rc1.isParametricallyIncompatibleWith(rc2)){
                return null;//impossible combination of RC's
            }
        }
        
        if(pruneMat.isPruned(RCs))
            return null;//pruned tuple...effectively infinite energy
        
        MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(makeObjectiveFunction(RCs));
        //new MoleculeModifierAndScorer(termE,confSpace,RCs);
         
        DoubleMatrix1D bestDOFVals;

        if(mof.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
            CCDMinimizer ccdMin = new CCDMinimizer(mof,true);
            bestDOFVals = ccdMin.minimize().dofValues;
        }
        else//molecule is already in the right, rigid conformation
            bestDOFVals = DoubleFactory1D.dense.make(0);

        double minEnergy = mof.getValue(bestDOFVals);
        EPoly epicFit = compEPICFit(mof,minEnergy,bestDOFVals,RCs);
        
        EnergyFunction.Tools.cleanIfNeeded(mof.getEfunc());//apparently some of these efuncs can get messy
        return epicFit;
    }
    
    
    private MoleculeObjectiveFunction makeObjectiveFunction(RCTuple RCs){
        switch(RCs.size()){
            case 1:
                return confECalc.makeIntraShellObjFcn(RCs.pos.get(0), RCs.RCs.get(0));
            case 2:
                return confECalc.makePairwiseObjFcn(RCs.pos.get(0), RCs.RCs.get(0), RCs.pos.get(1), RCs.RCs.get(1));
            default:
                throw new RuntimeException("ERROR: Can't calculate EPIC term for RCTuple "+RCs.stringListing());
        }
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
        
        ObjectiveFunction ofEPIC = new MoleculeModifierAndScorer( epicEF, mof.getConstraints(),
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
        
        epicEF.unassignSharedMolec();
        
        return true;
    }
    
    
    
}
