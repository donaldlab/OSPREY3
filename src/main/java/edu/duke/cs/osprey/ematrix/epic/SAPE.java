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

package edu.duke.cs.osprey.ematrix.epic;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache.ResPair;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.ResidueForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.energy.forcefield.SparseFFEnergy;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * Sparse Atom-Pair Energies
 * This is a listing of atom-pair energies for just a small number of the atom pairs
 * involved in an interactions
 * 
 * @author mhall44
 */
public class SAPE implements Serializable {
    
    //We want to be able to evaluate these SAPE in two cases:
    //as a stand-alone energy term, or 
    //as a term in an MolecEObjFunction where the molecule DOFs are set once
    //and then used for all SAPE and other terms directly calculated from conformation
    
    //stuff for standalone evaluation
    //Molecule mStandalone;//PROBABLY JUST HANDLE W/I OBJ FCN
    MoleculeModifierAndScorer mofStandalone;//objective function that sets the molecule DOFs and 
    //returns the energy
    
    //stuff for shared-molecule evaluation
    //we assume the shared molecule has the same residue numbers as the original molecule!
    
    //information sufficient to regenerate a force-field energy, given a molecule
    ForcefieldParams ffParams;
    ArrayList<int[]> interactingRes;//list of pairs of interacting residues
    ArrayList<ArrayList<int[]>> atomPairList;//for each pair of interacting residues, list of interacting atom pairs
    EnergyFunction sharedMolecEnergyFunction = null;//energy function linked to the shared molecule
    
    
    
    public SAPE (MoleculeModifierAndScorer objFcn, double distCutoff, DoubleMatrix1D[] sampAbs){
        //We'll initialize the SAPE for standalone work, and also
        //initialize information sufficient to reinitialize with a shared molecule
        
        
        //The standalone molecule and energy function for SAPE will belong only to the SAPE object
        //(we'll deep-copy, deleting fields not needed for SAPE)
        EnergyFunction mainEF = objFcn.getEfunc();
        objFcn.setEfunc(null);
        mofStandalone = (MoleculeModifierAndScorer) ObjectIO.deepCopy(objFcn);
        objFcn.setEfunc(mainEF);
        
        
        Molecule standaloneMolec = mofStandalone.getMolec();
        //compressMolecule(standaloneMolec);//if needed...see comments below (compressSVE)
        
        ffParams = null;
        findFFParams(mainEF);//get forcefield params from some component of mainEF
        if(ffParams==null){
            throw new RuntimeException("ERROR: Trying to approximate energy function with SAPE but can't find any forcefield-type energies");
        }

        pickAtomPairs(mainEF,distCutoff,sampAbs);//select which atom-pairs to include in SAPE
        
        //now make a standalone energy function including only these pairs
        EnergyFunction standaloneEFcn = makeSparseEFcn( standaloneMolec );
        mofStandalone.setEfunc(standaloneEFcn);
    }
    
    

    
    /*
    static ContSCObjFunction compressSVE(ContSCObjFunction of, SparseVDWEnergy sve){
        //store of and m as copies in their current state!
        //we temporarily remove the biggest things to avoid running out of memory during deepCopy
        //return copied objective function
            
        EnergyFunction storeEF = of.efunc;
        of.efunc = null;
        
        ContSCObjFunction retObjFcn = null;
        
        try{retObjFcn = (ContSCObjFunction)KSParser.deepCopy(of);}
        catch(Exception e){
            System.err.println("woah sveOF deepCopy failed");
            System.err.println(e.getMessage());
        }
        //ans.sveOF.efunc = null;
        of.m.bondedMatrix = storeBondedMatrix;
        of.m.nonBonded = storeNonBonded;
        of.efunc = storeEF;

        //won't need these big objects that aren't needed to setDOFs...
        retObjFcn.de = null;
        retObjFcn.strandRot = null;//won't need this

        sve.m = retObjFcn.m;

        //now clear out the atoms we won't need for more compact storage...
        //we need atoms in residues that are either (1) flexible, or (2) 
        //involved in perturbations that affect flexible residues
        HashSet<Residue> resNeeded = new HashSet<>();
        for(Residue res : sve.m.residue){
            if(res.flexible){
                resNeeded.add(res);
                for(int p : res.perts){
                    for(int molResNum : sve.m.perts[p].resDirectlyAffected){
                        //the resDirectlyAffected need to be complete
                        //because they are used to apply the perturbation
                        resNeeded.add(sve.m.residue[molResNum]);
                    }
                }
            }
        }

        for(Residue res : sve.m.residue){

            if(!resNeeded.contains(res)){

                for(int atNum=0; atNum<res.numberOfAtoms; atNum++){
                    int molAtNum = res.atom[atNum].moleculeAtomNumber;
                    res.atom[atNum] = null;
                    sve.m.atom[molAtNum] = null;
                }
            }
        }

        //since we are only using sve.m for sveOF.setDOFs and sve.getEnergy()
        //we can delete unnecessary bulky parts of the molecule
        //sve.m.bondedMatrix = null;//this is now removed before deepCopy
        //sve.m.nonBonded = null;
        sve.m.gradient = null;
        sve.m.connected = null;
        sve.m.connected12 = null;
        sve.m.connected13 = null;
        sve.m.connected14 = null;


        sve.DOFVals = retObjFcn.curDOFVals;
        
        return retObjFcn;
    }
    */
    

    
    
    
    private void pickAtomPairs(EnergyFunction EFToApprox, double distCutoff, DoubleMatrix1D[] sampAbs){
        //given the interacting residues in forcefield-type energy-function we're approximating,
        //pick the atom pairs to include
        
        //first figure out the interacting residues
        interactingRes = new ArrayList<>();
        addInteractingRes(EFToApprox);
        
        if(interactingRes.isEmpty()){
            throw new RuntimeException("ERROR: Trying to get SAPE approximation of energy function with no electrostatic/VDW interactions");
        }

        
        Molecule molec = mofStandalone.getMolec();
        
        //now pick the interacting atoms for each pairs of residues
        //based on the distances between them at each sample
        ArrayList<boolean[][]> atomPairsCloseEnough = new ArrayList<>();
        for(int resPair[] : interactingRes){
            Residue res1 = molec.residues.get(resPair[0]);
            Residue res2 = molec.residues.get(resPair[1]);
            atomPairsCloseEnough.add(new boolean[res1.atoms.size()][res2.atoms.size()]);
        }
        
        //For each pair of atoms at each residue pair whose interactions we're considering,
        //figure out if that atom pair is close enough to count in any of the samples
        for(DoubleMatrix1D samp : sampAbs){
                mofStandalone.setDOFs(samp);//apply the sample geometry to the molecule
                for(int resPairNum=0; resPairNum<interactingRes.size(); resPairNum++){
                    
                    Residue res1 = molec.residues.get(interactingRes.get(resPairNum)[0]);
                    Residue res2 = molec.residues.get(interactingRes.get(resPairNum)[1]);
                    
                    for(int at1=0; at1<res1.atoms.size(); at1++){
                        for(int at2=0; at2<res2.atoms.size(); at2++){
                            
                            if(res1==res2 && at2>=at1)//avoid double-counting in intra interactions
                                break;//(and also no need to include the interaction of an atom with itself)
                            
                            double dist = VectorAlgebra.distance(res1.coords, at1, res2.coords, at2);
                            
                            if(dist<distCutoff)
                                atomPairsCloseEnough.get(resPairNum)[at1][at2] = true;
                        }
                    }
                }
        }
    
        //now convert the boolean matrices to lists of atom pairs
        atomPairList = new ArrayList<>();
        for(boolean[][] closenessMatrix : atomPairsCloseEnough){
            ArrayList<int[]> atomPairsAtResPair = new ArrayList<>();
            
            for(int at1=0; at1<closenessMatrix.length; at1++){
                for(int at2=0; at2<closenessMatrix[at1].length; at2++){
                    
                    if(closenessMatrix[at1][at2]){
                        atomPairsAtResPair.add(new int[] {at1,at2});
                    }
                }
            }
            
            atomPairList.add(atomPairsAtResPair);
        }
        
        
        //finally, remove residue pairs that ended up not interacting
        for(int resPairNum = interactingRes.size()-1; resPairNum>=0; resPairNum--){
            //traverse backwards to avoid messing up numbering
            if(atomPairList.get(resPairNum).isEmpty()){
                atomPairList.remove(resPairNum);
                interactingRes.remove(resPairNum);
            }
        }
    }
    
    
    private void addInteractingRes(EnergyFunction EFToApprox) {
        //pick out the interacting residues in the energy function, and add them to interactingRes
                
        if(EFToApprox instanceof MultiTermEnergyFunction){//need to recurse to terms
            for( EnergyFunction term : ((MultiTermEnergyFunction)EFToApprox).getTerms() ){
                addInteractingRes(term);
            }
        }
        
        //we could be approximating either a single-residue or residue-pair energy
        else if(EFToApprox instanceof SingleResEnergy){
            Residue res = ((SingleResEnergy)EFToApprox).getRes();
            
            interactingRes.add( new int[] {res.indexInMolecule, res.indexInMolecule} );
        }
        
        else if(EFToApprox instanceof ResPairEnergy){
            Residue res1 = ((ResPairEnergy)EFToApprox).getRes1();
            Residue res2 = ((ResPairEnergy)EFToApprox).getRes2();
            
            interactingRes.add( new int[] {res1.indexInMolecule, res2.indexInMolecule} );
        }
        
        else if(EFToApprox instanceof ResidueForcefieldEnergy){
            for(ResPair rp : ((ResidueForcefieldEnergy)EFToApprox).resPairs){
                interactingRes.add( new int[] {rp.resIndex1,rp.resIndex2} );
            }
        }
        
        //other types of energies not approximated in SAPE.  
    }
    
    
    private void findFFParams(EnergyFunction EFToApprox){
        //Given the energy function we're approximating, choose forcefield parameters for this SAPE object
        //we expect that EFToApprox has at least some terms containing force-field type parameters
        //we'll look recursively for a forcefield type energy as in addInteractingRes
        //we call this function when ffParams is null and needs to be filled in
        
        if(EFToApprox instanceof ResidueForcefieldEnergy){
            ffParams = ((ResidueForcefieldEnergy)EFToApprox).resPairCache.ffparams;
            return;
        }
        
        if(EFToApprox instanceof MultiTermEnergyFunction){//need to recurse to terms
            for( EnergyFunction term : ((MultiTermEnergyFunction)EFToApprox).getTerms() ){
                findFFParams(term);
                if(ffParams!=null)//found params
                    return;
            }
        }
        
        //we could be approximating either a single-residue or residue-pair energy
        else if(EFToApprox instanceof SingleResEnergy)
            ffParams = ((SingleResEnergy)EFToApprox).getFFParams();
        else if(EFToApprox instanceof ResPairEnergy)
            ffParams = ((ResPairEnergy)EFToApprox).getFFParams();
    }
        
    
    
    EnergyFunction makeSparseEFcn(Molecule m){
        //generate a forcefield energy function for the atom pairs of interest
        //that evaluates the energy based on the conformation of m
        
        MultiTermEnergyFunction ans = new MultiTermEnergyFunction();//we'll make a term for each interacting residue pair
        
        for(int resPairNum=0; resPairNum<interactingRes.size(); resPairNum++){
            
            int[] resPair = interactingRes.get(resPairNum);
            
            Residue res1 = m.residues.get(resPair[0]);
            Residue res2 = m.residues.get(resPair[1]);

            ArrayList<Atom[]> interactingAtomPairs = new ArrayList<>();
            for(int atNumPair[] : atomPairList.get(resPairNum)){
                Atom atom1 = res1.atoms.get(atNumPair[0]);
                Atom atom2 = res2.atoms.get(atNumPair[1]);
                interactingAtomPairs.add( new Atom[] {atom1,atom2} );
            }

            SparseFFEnergy sparseE = new SparseFFEnergy( interactingAtomPairs, ffParams );
            ans.addTerm(sparseE);
        }
        
        return ans;
    }
    
    
    
    
    
    public void assignSharedMolecule(Molecule m){
        sharedMolecEnergyFunction = makeSparseEFcn(m);
    }
    
    
    
    public double getEnergyStandalone(DoubleMatrix1D coordVals){
        return mofStandalone.getValue(coordVals);
    }
    
    
    public double getEnergySharedMolec(){
        //in this case, we assume the shared molecule is already in the right conformation
        if(sharedMolecEnergyFunction==null)
            throw new RuntimeException("ERROR: Haven't initialized shared-molecule energy function for SAPE");
        
        return sharedMolecEnergyFunction.getEnergy();
    }
    
}
