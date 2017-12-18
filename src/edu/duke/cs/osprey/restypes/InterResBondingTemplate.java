/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.io.Serializable;
import java.util.Set;

/**
 *
 * This class is meant to be part of a residue template; it specifies
 * what inter-residue bonds the template can form
 * Since residues may be disordered, we do not always expect these bonds to be present though
 *
 * @author mhall44
 */
public abstract class InterResBondingTemplate implements Serializable {


    public abstract void connectInterResBonds(Residue res, boolean lowerNumResOnly);
    //connect the inter-residue bonds involving res.  If lowerNumResOnly
    //only connect to lower-numbered residue (this is for when we're connecting up the whole molecule)
    //This function does not mark a residue as having inter-res bonds connected...that's up
    //to the residue or molecule calling it (would be confusing because of lowerNumResOnly)

    public abstract boolean atomCanBondOtherRes(Atom atom);
    //atom is from this residue; can it bond other residues?

    public static class NoBondingTemplate extends InterResBondingTemplate {
        //for residues that can't bond to anything else
        @Override
        public void connectInterResBonds(Residue res, boolean lowerNumResOnly) {}

        @Override
        public boolean atomCanBondOtherRes(Atom atom) {
            return false;
        }
    }

    public static class PeptideBondingTemplate extends InterResBondingTemplate {
        @Override
        public void connectInterResBonds(Residue res, boolean lowerNumResOnly) {
            if(res.indexInMolecule>0){//res is not the first residue
                Residue prevRes = res.molec.residues.get(res.indexInMolecule-1);
                tryToMakePeptideBond(prevRes,res);
            }
            if(res.indexInMolecule<res.molec.residues.size()-1
                    && !lowerNumResOnly){//res is not the last residue, and we need to connect to higher-numbered res
                Residue nextRes = res.molec.residues.get(res.indexInMolecule+1);
                tryToMakePeptideBond(res,nextRes);
            }
        }

        @Override
        public boolean atomCanBondOtherRes(Atom atom) {
            return atom.name.equalsIgnoreCase("N") || atom.name.equalsIgnoreCase("C");
        }
    }


    public static class CysteineBondingTemplate extends PeptideBondingTemplate {
        @Override
        public void connectInterResBonds(Residue res, boolean lowerNumResOnly) {
            super.connectInterResBonds(res, lowerNumResOnly);//connect peptide bonds

            for(Residue res2 : res.molec.residues){
                if(res2!=res){
                    if(res2.indexInMolecule<res.indexInMolecule || !lowerNumResOnly){
                        if(res2.template.interResBonding instanceof CysteineBondingTemplate)
                            tryToMakeDisulfideBond(res,res2);
                    }
                }
            }
        }

        @Override
        public boolean atomCanBondOtherRes(Atom atom) {
            return atom.name.equalsIgnoreCase("N") || atom.name.equalsIgnoreCase("C")
                    || atom.name.equalsIgnoreCase("SG");
        }
    }


    public static class SpecifiedBondingAtomsTemplate extends InterResBondingTemplate {
        //A list of atom names in this residue is specified as being available for bonding
        //other residues; it will bond anything within range that is available for bonding
        Set<String> bondingAtomNames;

        public SpecifiedBondingAtomsTemplate(Set<String> bondingAtNames){
            bondingAtomNames = bondingAtNames;
        }

        @Override
        public void connectInterResBonds(Residue res, boolean lowerNumResOnly) {
            for(Residue res2 : res.molec.residues){
                if(res2!=res){
                    if(res2.indexInMolecule<res.indexInMolecule || !lowerNumResOnly){
                        for(int atnum2=0; atnum2<res2.atoms.size(); atnum2++){
                            Atom at2 = res2.atoms.get(atnum2);
                            if(res2.template.interResBonding.atomCanBondOtherRes(at2)){
                                for(Atom at : res.atoms){
                                    if(atomCanBondOtherRes(at)){
                                        double dist = VectorAlgebra.distance(res.coords,
                                                at.indexInRes, res2.coords, atnum2);
                                        if(dist<HardCodedResidueInfo.bondLengthUpperBound(at.elementNumber,at2.elementNumber)){
                                            at.addBond(res2.atoms.get(atnum2));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        @Override
        public boolean atomCanBondOtherRes(Atom atom) {
            return bondingAtomNames.contains(atom.name);
        }

        public Set<String> getBondingAtomNames(){return bondingAtomNames;}
    }



    public static void tryToMakePeptideBond(Residue res1, Residue res2){
        //Given consecutive residues res1 and res2, make a peptide bond between them if appropriate
        /*if(res1.template==null || res2.template==null) {
            throw new RuntimeException("ERROR: Trying to peptide-bond residue without template");
        }*/ // sometimes this is useful even without template
        if( (HardCodedResidueInfo.hasAminoAcidBB(res1) && HardCodedResidueInfo.hasAminoAcidBB(res2)) ){//can only peptide-bond amino acids

            int CIndex = res1.getAtomIndexByName("C");
            int NIndex = res2.getAtomIndexByName("N");

            //Get distance between these atoms
            double CNDist = VectorAlgebra.distance(res1.coords, CIndex, res2.coords, NIndex);
            if(CNDist<HardCodedResidueInfo.bondLengthUpperBound(6,7)){
                Atom C = res1.atoms.get(CIndex);
                Atom N = res2.atoms.get(NIndex);
                C.addBond(N);
            }
        }
    }


    public static void tryToMakeDisulfideBond(Residue res1, Residue res2){
        //Given CYX residues res1 and res2, make a disulfide bond between them if appropriate
        if(res1.template==null || res2.template==null)
            throw new RuntimeException("ERROR: Trying to disulfide-bond residue without template");

        int SIndex1 = res1.getAtomIndexByName("SG");
        int SIndex2 = res2.getAtomIndexByName("SG");

        if(SIndex1==-1 || SIndex2==-1)
            throw new RuntimeException("ERROR: Trying to disulfide-bond residue without SG");

        //Get distance between these atoms
        double SSDist = VectorAlgebra.distance(res1.coords, SIndex1, res2.coords, SIndex2);
        if(SSDist<HardCodedResidueInfo.bondLengthUpperBound(16,16)){
            Atom S1 = res1.atoms.get(SIndex1);
            Atom S2 = res2.atoms.get(SIndex2);
            S1.addBond(S2);
        }
    }


}
