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

package edu.duke.cs.osprey.dof.deeper.perts;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DOFBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.deeper.GenChi1Calc;
import edu.duke.cs.osprey.dof.deeper.ResBBState;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.TreeSet;

/**
 *
 * In DEEPer, perturbations don't necessarily commute
 * We'll put them into "perturbation blocks," which do
 * this way, when we want to apply a parameter value for a perturbation,
 * we can appropriately handle any perturbations that come after it
 * 
 * 
 * Note: All sidechains within a perturbation block will be idealized (in the C-beta sense)
 * since the original sidechain conformations are no longer valid for perturbed backbones
 * 
 * @author mhall44
 */
public class PerturbationBlock implements Serializable, DOFBlock {
    
    
    ArrayList<Perturbation> perts;
    //perturbations, in order of application
    
    
    ArrayList<LinkedHashMap<Residue,ResBBState>> prePertBBStates;
    //Backbone states immediately prior to application of each perturbation
    //same order of perturbation as in perts
    
    
    ArrayList<ArrayList<Perturbation>> successors;
    //a successor s of a perturbation p is a perturbation whose starting point 
    //depends on the parameter value of p
    //(i.e. s must come after p in perts, 
    //and there must be a chain p,q,...,r,s (in the same order as perts)
    //such that p and q have overlapping res, ... , r and s have overlapping res)
    
    
    ArrayList<ArrayList<Residue>> dependentResidues;
    //for each perturbation p, a list of all the residues whose conformations depend on the parameter value of p
    
    ArrayList<Residue> allResidues;//all residues in the block
    
    
    
    
    public PerturbationBlock(){
        
    }
    
    
    //initialize block with all sidechains idealized (preserving gen chi1),
    //but all perturbations at 0 parameter (so no backbone motion
    //and thus all pre-pert states are just the initial BB state)
    public PerturbationBlock(ArrayList<Perturbation> perts){
        
        this.perts = perts;
        
        prePertBBStates = new ArrayList<>();
        
        //assign the perturbations to this block
        for(int pertNum=0; pertNum<perts.size(); pertNum++){
            Perturbation pert = perts.get(pertNum);
            pert.block = this;
            pert.indexInBlock = pertNum;
            
            prePertBBStates.add(new LinkedHashMap<>());
        }
                
        computeSuccessors();
        computeDependentResidues();
        
        //Before any perturbing begins, we have the current BB conformation.  Record it.
        initPrePertStates();
        
        //we also need to idealize all the sidechains
        for(Residue res : allResidues){
            double chi1 = GenChi1Calc.getGenChi1(res);
            SidechainIdealizer.idealizeSidechain(EnvironmentVars.resTemplates, res);
            GenChi1Calc.setGenChi1(res, chi1);
        }
    }
    
    
    
    void updateSuccessorPrePertStates(int indexInBlock){
        //we've just applied the perturbation with the specified index in block
        //so the current molecular geometry defines the starting BB states for
        //ensuing perturbations
        for( Residue res : dependentResidues.get(indexInBlock) ){
            ResBBState state = new ResBBState(res);
            
            //set as pre-pert state for all perturbations after current one
            for(int index=indexInBlock+1; index<perts.size(); index++){
                prePertBBStates.get(index).put(res, state);
            }
        }
    }
    
    private void initPrePertStates(){
        //No perturbations have been performed yet
        //so initialize the current conformation to be the pre-pert state
        //for all perturbations
        for(Residue res : allResidues){
            ResBBState state = new ResBBState(res);
            
            //set as pre-pert state for all perturbations
            for(int index=0; index<perts.size(); index++){
                prePertBBStates.get(index).put(res, state);
            }
        }
    }
    
    
    private void computeSuccessors(){
        //compute the successors for all the perturbations
        
        successors = new ArrayList<>();//initialize list with null so can add in reverse order
        for(int pertNum=0; pertNum<perts.size(); pertNum++)
            successors.add(new ArrayList<>());
        
        for(int pertNum=perts.size()-1; pertNum>=0; pertNum--){
            TreeSet<Integer> successorIndices = new TreeSet<>();
            
            for(int pert2=pertNum+1; pert2<perts.size(); pert2++){
                if(perturbationsOverlap(pertNum,pert2)){//if pert2 overlaps pertNum, then it's a successor, and so are (indirectly) all its successors
                    successorIndices.add(pert2);
                    for(Perturbation pert2Successor : successors.get(pert2))
                        successorIndices.add(pert2Successor.indexInBlock);
                }
            }
            
            //now convert successorIndices to list of perturbations (in same order as perts, because using TreeSet)
            for(int ind : successorIndices)
                successors.get(pertNum).add(perts.get(ind));
        }
    }
    
    
    private void computeDependentResidues(){
        //compute the dependent residues for all perturbations
        
        dependentResidues = new ArrayList<>();
        LinkedHashSet<Residue> allRes = new LinkedHashSet<>();
        
        for(int pertNum=0; pertNum<perts.size(); pertNum++){
            
            Perturbation pert = perts.get(pertNum);
            LinkedHashSet<Residue> pertDep = new LinkedHashSet<>();
            
            for(Residue res : pert.resDirectlyAffected)
                pertDep.add(res);
            
            for(Perturbation successor : successors.get(pertNum)){
                for(Residue res : successor.resDirectlyAffected)
                    pertDep.add(res);
            }
            
            dependentResidues.add(new ArrayList<>(pertDep));
            allRes.addAll(pertDep);
        }
        
        allResidues = new ArrayList<>(allRes);
    }
    
    
    private boolean perturbationsOverlap(int index1, int index2){
        //Do the perturbations with the specified indices overlap in terms of their directly affected residues?
        for(Residue res1 : perts.get(index1).resDirectlyAffected){
            for(Residue res2 : perts.get(index2).resDirectlyAffected){
                if(res1==res2)
                    return true;
            }
        }
        
        return false;
    }

    @Override
    public DOFBlock copyForNewMolecule(Molecule mol, LinkedHashMap<DegreeOfFreedom, DegreeOfFreedom> copiedDOFMap) {
        
        PerturbationBlock copiedBlock = new PerturbationBlock();
        
        //Make copies of all the perturbations
        copiedBlock.perts = new ArrayList<>();
        for(Perturbation pert : perts){
            Perturbation copiedPert = pert.copyForNewMolecule(mol,copiedBlock);
            copiedDOFMap.put(pert, copiedPert);
            copiedBlock.perts.add(copiedPert);
        }
        
        copiedBlock.prePertBBStates = new ArrayList<>();
        for(LinkedHashMap<Residue,ResBBState> bbMap : prePertBBStates){
            LinkedHashMap<Residue,ResBBState> copiedBBMap = new LinkedHashMap<>();
            for(Residue res : bbMap.keySet()){
                Residue otherRes = res.equivalentInMolec(mol);
                copiedBBMap.put( otherRes, new ResBBState(bbMap.get(res)) );
            }
            copiedBlock.prePertBBStates.add(copiedBBMap);
        }
        
        copiedBlock.successors = new ArrayList<>();
        for(ArrayList<Perturbation> succ : successors){
            ArrayList<Perturbation> copiedSucc = new ArrayList<>();
            for(Perturbation pert : succ){
                copiedSucc.add( (Perturbation) copiedDOFMap.get(pert) );
            }
            copiedBlock.successors.add(copiedSucc);
        }
       
        copiedBlock.dependentResidues = new ArrayList<>();
        for(ArrayList<Residue> dep : dependentResidues){
            copiedBlock.dependentResidues.add( Residue.equivalentInMolec(dep,mol) );
        }
        
        copiedBlock.allResidues = Residue.equivalentInMolec(allResidues, mol);
    
        return copiedBlock;
    }

    public ArrayList<Perturbation> getPerts() {
        return perts;
    }
    
    
    
}
