/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

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
import java.util.ArrayList;
import java.util.Arrays;

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


	public void addEnergyTerms(boolean doIntra, int... resToCalc) {
		// merge residues at the specified positions; this updates the energy matrix.
		// the energy matrix currently only support pairs and triples
		TermECalculator hotECalc = new TermECalculator(searchSpace, shellResidues, 
				doEPIC, doIntra, pruneMat, epicSettings, addResEntropy, resToCalc);

		Object hotE = hotECalc.doCalculation();
		storeEnergy(hotE, resToCalc);
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
			emat.eRefMat = new ReferenceEnergies(searchSpace);
			System.out.println("CORRECTING ENERGY MATRIX BASED ON REFERENCE ENERGIES");
			emat.eRefMat.correctEnergyMatrix(emat);
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

		if(doEPIC){
			if(res.length==1)//intra+shell energy
				epicMat.oneBody.set( res[0], (ArrayList<EPoly>) calcResult );
			else if(res.length==2)// pairwise
				epicMat.pairwise.get(res[0]).set( res[1], (ArrayList<ArrayList<EPoly>>) calcResult );
			else
				throw new RuntimeException("ERROR: have not implemented three body terms with EPIC");
		}
		else {

			if(res.length==1)//intra+shell energy
				emat.oneBody.set( res[0], (ArrayList<Double>) calcResult );

			else if(res.length==2)//pairwise
				emat.pairwise.get(res[0]).set( res[1], (ArrayList<ArrayList<Double>>) calcResult );

			else if(res.length==3){// three-body	
				ArrayList<ArrayList<ArrayList<Double>>> threeBody = (ArrayList<ArrayList<ArrayList<Double>>>)calcResult;

				for(int rc0 = 0; rc0 < threeBody.size(); ++rc0 ) {
					ArrayList<ArrayList<Double>> indexPos1 = threeBody.get(rc0);

					for(int rc1 = 0; rc1 < indexPos1.size(); ++rc1 ) {
						ArrayList<Double> indexPos2 = indexPos1.get(rc1);

						for(int rc2 = 0; rc2 < indexPos2.size(); ++rc2) {
							double threeBodyE = threeBody.get(rc0).get(rc1).get(rc2);

							// subtract out pairwise energies, which have already been set, from this RC
							double pairWiseE = emat.getPairwise(res[0], rc0, res[1], rc1);
							pairWiseE += emat.getPairwise(res[0], rc0, res[2], rc2);
							pairWiseE += emat.getPairwise(res[1], rc1, res[2], rc2);

							RCTuple triple = new RCTuple( new ArrayList<>(Arrays.asList(res[0], res[1], res[2])), 
									new ArrayList<>(Arrays.asList(rc0, rc1, rc2)) );

							emat.setHigherOrder(triple, threeBodyE-pairWiseE);
						}
					}
				}
			}

			else if(res.length==4){// four-body

				ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> fourBody = (ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>)calcResult;
				Integer[] ires = new Integer[res.length]; for(int i = 0; i < ires.length; ++i) ires[i] = new Integer(res[i]);
				Integer[] rcs = new Integer[ires.length];

				for(rcs[0] = 0; rcs[0] < fourBody.size(); ++rcs[0]) {
					ArrayList<ArrayList<ArrayList<Double>>> indexPos1 = fourBody.get(rcs[0]);

					for(rcs[1] = 0; rcs[1] < indexPos1.size(); ++rcs[1] ) {
						ArrayList<ArrayList<Double>> indexPos2 = indexPos1.get(rcs[1]);

						for(rcs[2] = 0; rcs[2] < indexPos2.size(); ++rcs[2]) {
							ArrayList<Double> indexPos3 = indexPos2.get(rcs[2]);

							for(rcs[3] = 0; rcs[3] < indexPos3.size(); ++rcs[3]) {
								double fourBodyE = fourBody.get(rcs[0]).get(rcs[1]).get(rcs[2]).get(rcs[3]);

								// subtract out pairwise energies
								double pairWiseE = 0;
								for(int i = 0; i < ires.length; ++i) {
									for(int j = i+1; j < ires.length; ++j) {
										pairWiseE += emat.getPairwise(ires[i], rcs[i], ires[j], rcs[j]);
									}
								}

								RCTuple quadruple = new RCTuple(new ArrayList<>(Arrays.asList(ires)), new ArrayList<>(Arrays.asList(rcs)));

								emat.setHigherOrder(quadruple, fourBodyE-pairWiseE);
							}
						}
					}
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
