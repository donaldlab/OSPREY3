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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.ematrix.epic.EPoly;
import edu.duke.cs.osprey.handlempi.MPIMaster;
import edu.duke.cs.osprey.handlempi.MPISlaveTask;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;

/**
 *
 * @author mhall44
 */
public class EnergyMatrixCalculator {
    
    ConfSpace searchSpace;
    ArrayList<Residue> shellResidues;
    
    
    boolean doEPIC;
    //if we're doing EPIC, then we'll be wanting to avoid pruned RCs, and we'll need EPICSettings
    PruningMatrix pruneMat = null;
    EPICSettings epicSettings = null;
    
    
    boolean useERef = false;
    //If using E ref, will compute a reference energy for each amino-acid type
    //and correct the intra+shell energies based on the reference energies
    boolean addResEntropy = false;//add residue entropy to one-body energies
    
    //We are calculating either a scalar or EPIC matrix, so we'll allocate and fill in just one of these
    private EnergyMatrix emat = null;
    private EPICMatrix epicMat = null;
    
    
    //constructor for calculating a scalar energy matrix (rigid or pairwise lower bounds)
    public EnergyMatrixCalculator(ConfSpace s, ArrayList<Residue> sr, boolean useERef, 
            boolean addResEntropy) {
        
        searchSpace = s;
        shellResidues = sr;
        doEPIC = false;
        this.useERef = useERef;
        this.addResEntropy = addResEntropy;
    }
    
    
    //Constructor for calculating an EPIC matrix
    public EnergyMatrixCalculator(ConfSpace s, ArrayList<Residue> sr, PruningMatrix pr, EPICSettings es){
        searchSpace = s;
        shellResidues = sr;
        doEPIC = true;
        pruneMat = pr;
        epicSettings = es;
    }
    
    
	public EnergyMatrixCalculator(ConfSpace confSpace, ArrayList<Residue> shellResidues, 
			boolean doEPIC, PruningMatrix prm, EPICSettings es, 
			boolean addResEntropy, EnergyMatrix emat) {

		this.searchSpace = confSpace;
		this.shellResidues = shellResidues;
		this.doEPIC = doEPIC;
		this.pruneMat = prm;
		this.epicSettings = es;
		this.addResEntropy = addResEntropy;

		this.emat = emat;
	}
	
	
	public void addEnergyTerms(boolean doIntra, int... posNums) {
		// merge residues at the specified positions; this updates the energy matrix.
		// the energy matrix currently only support pairs and triples
		TermECalculator hotECalc = new TermECalculator(searchSpace, shellResidues, 
				doEPIC, doIntra, pruneMat, epicSettings, addResEntropy, posNums);

		Object hotEnergies = hotECalc.doCalculation();
		storeEnergy(hotEnergies, posNums);
	}
    
    
   //Calculate a pairwise energy matrix based on a pairwise energy function
   public void calcPEM(){
       
       System.out.println();
       if(doEPIC)
           System.out.println("BEGINNING EPIC MATRIX PRECOMPUTATION");
       else
           System.out.println("BEGINNING ENERGY MATRIX PRECOMPUTATION");
       System.out.println();
       
       initMatrix();
       
       if(EnvironmentVars.useMPI)
           calcPEMDistributed();
       else
           calcPEMLocally();
       
       if(useERef){
           System.out.println("COMPUTING REFERENCE ENERGIES");
           emat.seteRefMat(new ReferenceEnergies(searchSpace));
       }
       
       System.out.println("ENERGY MATRIX CALCULATION DONE");
   }
    
    
    public void calcPEMLocally(){
        //do the energy calculation here
        
        for(int res=0; res<searchSpace.numPos; res++){
            
            System.out.println("Starting intra+shell energy calculations for residue "+res);
            
            TermECalculator oneBodyECalc = new TermECalculator(searchSpace,shellResidues,doEPIC,
                    false,pruneMat,epicSettings,addResEntropy,res);
            
            Object oneBodyE = oneBodyECalc.doCalculation();
            storeEnergy(oneBodyE, res);

            for(int res2=0; res2<res; res2++){
                
                System.out.println("Starting pairwise energy calculations for residues "+res+", "+res2);
                
                TermECalculator pairECalc = new TermECalculator(searchSpace,shellResidues,doEPIC,
                        false,pruneMat,epicSettings,false,res,res2);
                Object pairE = pairECalc.doCalculation();
                storeEnergy(pairE, res, res2);
            }
        }
    }
    
    
    public void calcPEMDistributed(){
        //do energy calculation on slave nodes via MPI
        
        MPIMaster mm = MPIMaster.getInstance();//we'll only be running one MPI at once
        ArrayList<MPISlaveTask> tasks = new ArrayList<>();
        
        //generate TermMinECalc objects, in the same order as for local calculation,
        //but this time pass them off to MPI
        for(int res=0; res<searchSpace.numPos; res++){
            
            tasks.add( new TermECalculator(searchSpace,shellResidues,doEPIC,false,
                    pruneMat,epicSettings,addResEntropy,res) );

            for(int res2=0; res2<res; res2++)
                tasks.add( new TermECalculator(searchSpace,shellResidues,doEPIC,false,
                        pruneMat,epicSettings,false,res,res2) );
        }
        
        ArrayList<Object> calcResults = mm.handleTasks(tasks);
        
        //Now go through our task results in the same order and put the energies in our matrix
        int resultCount = 0;
        
        for(int res=0; res<searchSpace.numPos; res++){
            
            storeEnergy( calcResults.get(resultCount), res );
            resultCount++;

            for(int res2=0; res2<res; res2++){
                storeEnergy( calcResults.get(resultCount), res, res2 );
                resultCount++;
            }
        }
    }
    
    
    private void initMatrix(){
        //initialize the matrix we're calculating
        if(doEPIC)
            epicMat = new EPICMatrix(searchSpace, pruneMat.getPruningInterval());
        else
            emat = new EnergyMatrix(searchSpace, Double.POSITIVE_INFINITY);
            //all RCs included (infinite pruning interval)
    }
    
    private void storeEnergy(Object calcResult, int... res){
        //eCalc has performed its calculations, for the residue or pair denoted by res.
        //store the results of this calculation in our matrix.  
        
        if (doEPIC) {
            if (res.length==1)
                epicMat.setOneBody(res[0], (ArrayList<EPoly>)calcResult);
            else
                epicMat.setPairwise(res[0], res[1], (ArrayList<ArrayList<EPoly>>)calcResult);
        }
        else {
            if (res.length==1)
                emat.setOneBody(res[0], (ArrayList<Double>)calcResult);
            
            else if(res.length==2)
                emat.setPairwise(res[0], res[1], (ArrayList<ArrayList<Double>>)calcResult);
            
			else if(res.length > 2) {
				Integer[] pos = new Integer[res.length]; for(int i = 0; i < pos.length; ++i) pos[i] = new Integer(res[i]);

				HashMap<ArrayList<Integer>, Double> nBody2E = (HashMap<ArrayList<Integer>, Double>)calcResult;
				for(ArrayList<Integer> rc : nBody2E.keySet()) {

					double nbE = nBody2E.get(rc);

					// subtract out pairwise terms
					double pwE = 0;
					for(int i = 0; i < pos.length; ++i) {
						for(int j = i+1; j < pos.length; ++j) {
							pwE += emat.getPairwise(pos[i], rc.get(i), pos[j], rc.get(j));
						}
					}
					
					// AAO 2016: 
					// TODO: subtract out 3 to n-1 order terms
					
					RCTuple nBody = new RCTuple(new ArrayList<>(Arrays.asList(pos)), new ArrayList<>(rc));
					emat.setHigherOrder(nBody, nbE-pwE);
				}
			}
        }
    }
    
    public EnergyMatrix getEMatrix(){
        if(emat==null)
            throw new RuntimeException("ERROR: Energy matrix is null after calculation");
        
        return emat;
    }
    
    public EPICMatrix getEPICMatrix(){
        if(epicMat==null)
            throw new RuntimeException("ERROR: EPIC matrix is null after calculation");
        
        return epicMat;
    }
    
    
}
