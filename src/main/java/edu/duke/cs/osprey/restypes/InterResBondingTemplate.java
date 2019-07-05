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

package edu.duke.cs.osprey.restypes;

import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.io.Serializable;
import java.util.List;

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

	public abstract boolean isInterResBondedForward(Residue res1, Residue res2);
	public abstract boolean makeInterResBondForward(Residue res1, Residue res2);

    public static class NoBondingTemplate extends InterResBondingTemplate {

        //for residues that can't bond to anything else
        @Override
        public void connectInterResBonds(Residue res, boolean lowerNumResOnly) {}

        @Override
        public boolean atomCanBondOtherRes(Atom atom) {
            return false;
        }

        @Override
		public boolean isInterResBondedForward(Residue res1, Residue res2) {
        	return false;
		}

		@Override
		public boolean makeInterResBondForward(Residue res1, Residue res2) {
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

		@Override
		public boolean isInterResBondedForward(Residue res1, Residue res2) {

			// make res1:C to res2:N bond
			Atom C = res1.getAtomByName("C");
			Atom N = res2.getAtomByName("N");

			if (C == null || N == null) {
				return false;
			}

			return C.bonds.contains(N);
		}

		@Override
		public boolean makeInterResBondForward(Residue res1, Residue res2) {

			// check for res1:C to res2:N bond
			Atom C = res1.getAtomByName("C");
			Atom N = res2.getAtomByName("N");

			// no way to make a peptide bond? then fail
			if (C == null || N == null) {
				return false;
			}

			// don't worry if the atoms are too far apart, just make the bond
			C.addBond(N);
			return true;
		}
	}

	public static class NucleotideBondingTemplate extends InterResBondingTemplate {

		@Override
		public void connectInterResBonds(Residue res, boolean lowerNumResOnly) {
			if(res.indexInMolecule>0){//res is not the first residue
				Residue prevRes = res.molec.residues.get(res.indexInMolecule-1);
				tryToMakeNucleotideBond(prevRes,res);
			}
			if(res.indexInMolecule<res.molec.residues.size()-1
				&& !lowerNumResOnly){//res is not the last residue, and we need to connect to higher-numbered res
				Residue nextRes = res.molec.residues.get(res.indexInMolecule+1);
				tryToMakeNucleotideBond(res,nextRes);
			}
		}

		@Override
		public boolean atomCanBondOtherRes(Atom atom) {
			return atom.name.equalsIgnoreCase("P") || atom.name.equalsIgnoreCase("O3'");
		}

		@Override
		public boolean isInterResBondedForward(Residue res1, Residue res2) {

			// make res1:O3' to res2:P bond
			Atom O = res1.getAtomByName("O3'");
			Atom P = res2.getAtomByName("P");

			if (O == null || P == null) {
				return false;
			}

			return O.bonds.contains(P);
		}

		@Override
		public boolean makeInterResBondForward(Residue res1, Residue res2) {

			// check for res1:O3' to res2:P bond
			Atom O = res1.getAtomByName("O3'");
			Atom P = res2.getAtomByName("P");

			// no way to make a nucleotide bond? then fail
			if (O == null || P == null) {
				return false;
			}

			// don't worry if the atoms are too far apart, just make the bond
			O.addBond(P);
			return true;
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

		@Override
		public boolean isInterResBondedForward(Residue res1, Residue res2) {
			throw new UnsupportedOperationException("atom connectivity queries for cysteine bonds are not yet implemented");
		}

		@Override
		public boolean makeInterResBondForward(Residue res1, Residue res2) {
			throw new UnsupportedOperationException("atom connectivity queries for cysteine bonds are not yet implemented");
		}
	}


    public static class SpecifiedBondingAtomsTemplate extends InterResBondingTemplate {

        //A list of atom names in this residue is specified as being available for bonding
        //other residues; it will bond anything within range that is available for bonding
        List<String> bondingAtomNames;

        public SpecifiedBondingAtomsTemplate(List<String> bondingAtNames){
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

        public List<String> getBondingAtomNames(){return bondingAtomNames;}

        private static String notSupportedMsg = "atom connectivity queries for generic inter-residue bonds are not supported."
			+ " Generic inter-residue bonds do not define a \"forward\" direction";

		@Override
		public boolean isInterResBondedForward(Residue res1, Residue res2) {
			throw new UnsupportedOperationException(notSupportedMsg);
		}

		@Override
		public boolean makeInterResBondForward(Residue res1, Residue res2) {
			throw new UnsupportedOperationException(notSupportedMsg);
		}
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

	public static void tryToMakeNucleotideBond(Residue res1, Residue res2) {

		// Given consecutive residues res1 and res2, make a nucleotide bond between them if appropriate

		int Oi = res1.getAtomIndexByName("O3'");
		if (Oi < 0) {
			return;
		}
		int Pi = res2.getAtomIndexByName("P");
		if (Pi < 0) {
			return;
		}
		Atom O = res1.atoms.get(Oi);
		Atom P = res2.atoms.get(Pi);

		// Get distance between these atoms
		double dist = VectorAlgebra.distance(res1.coords, Oi, res2.coords, Pi);
		if (dist < HardCodedResidueInfo.bondLengthUpperBound(O.elementNumber, P.elementNumber)) {
			O.addBond(P);
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
