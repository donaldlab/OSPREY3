/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.structure;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResTemplateMatching;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.Protractor;
import edu.duke.cs.osprey.tools.StringParsing;
import edu.duke.cs.osprey.tools.VectorAlgebra;

import java.io.Serializable;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class Residue implements Serializable {
    //This is a residue
    //in the general sense: it's a piece of a molecular system that can mutate, take on RCs, etc.
    //and is the basic unit for energy calculations: energy is broken down to intra-residue energies,
    //pairwise interactions between residues, and maybe a few triples etc.
    
    //residue information
    public String fullName;//short name is the first three characters of this
    
    public int indexInMolecule = -1;//index of this residue in the molecule it's in
    public Molecule molec;//the molecule it's in
    public ResidueTemplate template = null;
    //this gives the name of the residue like "GLU" or "H2O" or whatever
    //and any information about bonds, force-field, rotatable dihedrals, and rotamers
    //ONLY ASSIGN USING ASSIGNTEMPLATE (will reorder atoms to match template)
    
    
    //Residue[] neighbors;//The ResidueTemplate specifies a certain number of residues that should be bonded to this one
    //usually two, through backbone bonds, but there could be exceptions (disulfides, unbonded ligands without neighbors,
    //weird non-protein stuff)
    //if neighbors[i] is null that means the i'th neighbor is missing from the structure
    
    //atom information
    public ArrayList<Atom> atoms;//information on individual atoms
    //(when assigning the template we will copy this from the template)
    
    
    public double coords[];//coordinates of all the atoms: x1, y1, z1, x2, ...
    //we want the coordinates to be all in one place and fairly cache-friendy
    //(e.g. when doing pairwise minimizations, polynomial fits, etc. we can easily cache these coordinates
    //for our residues of interest and access them quickly)
    
    
    //flags indicating that we've marked the bonds between the atoms
    //intra- and inter- may be marked at different times
    //forcefield energies can only be meaningfully evaluated if all the bonds are marked
    public boolean intraResBondsMarked = false;
    public boolean interResBondsMarked = false;//for bonds between a pair of residues,
    //if interResBondsMarked is true for both, then we expect bonds between the two to be in place
    
    
    public ArrayList<ConfProblem> confProblems = new ArrayList<>();
    //If the residue has a ring that might need idealization, we need to store its pucker DOF
    public ProlinePucker pucker = null;
    
    
    
    public int secondaryStruct;
    //Types of secondary structure
    public final static int HELIX = 0;
    public final static int SHEET = 1;
    public final static int LOOP = 2;
    
    
    
    
    public Residue(ArrayList<Atom> atomList, ArrayList<double[]> coordList, String nameFull, Molecule m){
        //generate a residue.  Don't put in template and related stuff like bonds (that can be assigned later)
        
        atoms = atomList;
        molec = m;
        
        int numAtoms = atoms.size();
        
        for(int a=0; a<numAtoms; a++){
            atoms.get(a).res = this;
            atoms.get(a).indexInRes = a;
        }
        
        
        //record coordinates, unless coordList is null (then leave coords null too)
        if(coordList==null)
            coords = null;
        else {
        
            if(numAtoms != coordList.size()){
                throw new RuntimeException("ERROR: Trying to instantiate residue with "+numAtoms+
                        " atoms but "+coordList.size()+" atom coordinates");
            }

            coords = new double[3*numAtoms];
            for(int a=0; a<numAtoms; a++){
                System.arraycopy(coordList.get(a), 0, coords, 3*a, 3);
            }
        }
        
        fullName = nameFull;
    }
    
    public String getPDBResNumber() {
            if (fullName.length() > 5)
                    return( (StringParsing.getToken(fullName.substring(5),1)) );
            return Integer.toString(indexInMolecule+1);
    }
    
    boolean assignTemplate (){
        //assign a template to this residue if possible, using the ResidueTemplateLibrary
        //return if successful or not
        
        //we'll use the default ResidueTemplateLibrary (from EnvironmentVars)
        ResidueTemplateLibrary templateLib = EnvironmentVars.resTemplates;
        
        //first see if there are any templates matching this name
        ArrayList<ResidueTemplate> templCandidates = new ArrayList<>();
        
        String templateName = HardCodedResidueInfo.getTemplateName(this);
        //Residues will usually have a template name that's the first three residue of their full name,
        //but there are a few exceptions
        
        for(ResidueTemplate templ : templateLib.templates){
            if(templ.name.equalsIgnoreCase(templateName))
                templCandidates.add(templ);
        }
        
        if(templCandidates.isEmpty())
            return false;
        
        
        //now try to match up atoms
        ResTemplateMatching bestMatching = null;
        double bestScore = Double.POSITIVE_INFINITY;
        
        for(ResidueTemplate templ : templCandidates){
            ResTemplateMatching templMatching = new ResTemplateMatching(this,templ);
            if(templMatching.score < bestScore){
                bestScore = templMatching.score;
                bestMatching = templMatching;
            }
        }
        
        if(bestMatching==null)//no good matching found
            return false;
        
        //now rename/order atoms in this residue to match template!
        bestMatching.assign();
        return true;//matched successfully!
    }
    
    /**
     * Computes the dihedral values for the template that is currently assigned to this type.
     * @return An array of dihedral values for the wildtype rotamer..
     */
    public double[] getCurrentRotamerDihedrals(){
    	int numDih = this.template.numDihedrals;
    	double dihedrals[] = new double[this.template.numDihedrals];
    	for(int dih = 0; dih < numDih; dih++){
    		// Get the name of the four atoms that move for each dihedral
    		String atomNamesDihedral[] = new String[4];
    		for(int i = 0; i < 4; i++){
    			atomNamesDihedral[i] = this.template.templateRes.atoms.get(template.dihedral4Atoms[dih][i]).name;
    		}
    		// Find the four atoms in the arraylist of atoms for this residue.
    		double curDihCoordinates[][] = new double[4][];
    		for(Atom at: this.atoms){
    			for (int i = 0; i < 4; i++){
    				if(at.name.equals(atomNamesDihedral[i])){
    					curDihCoordinates[i] = at.getCoords();
    				}
    			}
    		}
    		// Compute the dihedral of the four atoms.
    		dihedrals[dih] = Protractor.measureDihedral(curDihCoordinates);
    		
    	}
    	return dihedrals;
    }
    
    public void markIntraResBondsByTemplate(){
        //assign all the bonds between atoms in this residue, based on the template

        
        int numAtoms = atoms.size();
        ArrayList<Atom> templateAtoms = template.templateRes.atoms;
        
        if(templateAtoms.size()!=numAtoms)
            throw new RuntimeException("ERROR: Template for "+fullName+" has the wrong number of atoms");
        
        //Within-residue bonds are exactly as in the template
        boolean[][] intraResBondMatrix = templateIntraResBondMatrix();
     
        //now apply this bond matrix
        for(int atNum1=0; atNum1<numAtoms; atNum1++){
            
            Atom atom1 = atoms.get(atNum1);

            //check that all atoms have the same names as
            //corresponding template atoms, in the same order
            if( ! atom1.name.equalsIgnoreCase( templateAtoms.get(atNum1).name ) ){
                throw new RuntimeException("ERROR: Atom name mismatch with template found in "+fullName);
            }
            
            for(int atNum2=0; atNum2<numAtoms; atNum2++){
                if(intraResBondMatrix[atNum1][atNum2]){
                    Atom atom2 = atoms.get(atNum2);
                    atom1.bonds.add(atom2);
                }
            }
            
            /*
            if(atom1.bonds[bondedAtomCount]==null){
                //atom1 has bonds that are null in the template
                //need to look outside the residue
                
                //currently expecting just one interaction, which we expect will be to the closest atom
                //from another residue that is also expecting an out-of-residue bond
                if(bondedAtomCount<atom1.bonds.length-1)//more than one missing interaction
                    throw new UnsupportedOperationException("ERROR: Not yet supporting atoms"
                            + " with multiple bonds outside the residue");
                
                
                Atom bondedAtom = null;//the atom in another residue bonded to atom1...we'll search for this
                double bondedAtomDist = Double.POSITIVE_INFINITY;//atom1 -- bondedAtom distance
                
                for(Residue otherRes : molec.residues){
                    if(otherRes!=this){
                        
                        for(int otherAtNum=0; otherAtNum<otherRes.atoms.size(); otherAtNum++){
                            Atom otherAtom = otherRes.atoms.get(otherAtNum);
                            
                            for(Atom otherAtomBond : otherAtom.bonds){
                                if(otherAtomBond!=null){
                                    if(otherAtomBond.res==otherRes)
                                        continue;//no opening here to bond to atom1
                                }
                                    
                                //if bond is null or is to other residue, there is an opening
                                //compute 
                                double dist = VectorAlgebra.distance(coords,atNum1,otherRes.coords,otherAtNum);
                                
                                //if dist is shorter than distances to other bondedAtoms,
                                //then we'll want to make this bond
                                if(dist<bondedAtomDist){
                                    bondedAtomDist = dist;
                                    bondedAtom = otherAtom;
                                }
                            }
                        }
                    }
                }
                
                
                //if template not assigned to other residue, then leave null
                //will connect up when assigning
                if(bondedAtom.res.template!=null){
                    atom1.bonds[bondedAtomCount] = bondedAtom;
                }
            }*/
        }
        
        
        checkBonds();
        intraResBondsMarked = true;
    }
    
    
    public boolean[][] templateIntraResBondMatrix(){
        //matrix(i,j) indicates whether atoms i and j are bonded to each other in the template

        int numAtoms = atoms.size();
        ArrayList<Atom> templateAtoms = template.templateRes.atoms;
        
        boolean intraResBondMatrix[][] = new boolean[numAtoms][numAtoms];
        for(int atNum1=0; atNum1<numAtoms; atNum1++){
            for(int atNum2=0; atNum2<numAtoms; atNum2++){
                
                Atom atom1 = templateAtoms.get(atNum1);
                Atom atom2 = templateAtoms.get(atNum2);
                
                //check if atom2 is in bond list for atom1
                //if so set bond matrix boolean to true
                for(Atom bondedAtom : atom1.bonds){
                    if(bondedAtom==atom2){
                        intraResBondMatrix[atNum1][atNum2] = true;
                        break;
                    }
                }
            }
        }
        
        return intraResBondMatrix;
    }
    
    
    public void checkTemplateAtomNames(){
        //check that atom names match between the residue and the template
        int numAtoms = atoms.size();
        ArrayList<Atom> templateAtoms = template.templateRes.atoms;
        
        if(templateAtoms.size()!=numAtoms)
            throw new RuntimeException("ERROR: Template has wrong number of atoms for "+fullName);
        
        for(int atNum=0; atNum<numAtoms; atNum++){
            String resAtomName = atoms.get(atNum).name;
            String templateAtomName = templateAtoms.get(atNum).name;
            
            if(!resAtomName.equalsIgnoreCase(templateAtomName)){
                throw new RuntimeException("ERROR: Atom name mismatch between template and atom list in "+
                        fullName+": "+resAtomName+" in atom list, "+templateAtomName+" in template");
            }
        }
    }
    
    
    public void checkBonds(){
        //ensure no bonds are null, and all bonds are listed from both sides
        for(Atom atom1 : atoms){
            
            for(Atom atom2 : atom1.bonds){
                if(atom2==null)
                    throw new RuntimeException("ERROR: checkBonds found a null bond");
                
                //make sure atom2 lists the bond to atom1 too
                boolean bothSides = false;
                for(Atom bondedAtom : atom2.bonds){
                    if(bondedAtom==atom1){
                        bothSides = true;
                        break;
                    }
                }
                
                if(!bothSides){
                    throw new RuntimeException("ERROR: checkBonds found atom1 bonded to atom2"
                            + " but atom2 not bonded to atom1");
                }
            }
        }
    }
    
    
    
    public int getAtomIndexByName(String n){
        //index in atoms of atom with name n
        //-1 if there is no atom here by that name
        
        for(int atNum=0; atNum<atoms.size(); atNum++){
            if(atoms.get(atNum).name.equalsIgnoreCase(n))
                return atNum;
        }
        
        return -1;
    }
    
    public double[] getCoordsByAtomName(String n){
        //return coordinates of atom named n
        //return null if there is none
        int index = getAtomIndexByName(n);
        if(index==-1)
            return null;
        
        return atoms.get(index).getCoords();
    }
    
    public void setCoordsByAtomName(String name, double atomCoords[]){
        //Set the atom with this name to have the specified coordinates
        //Use this function only if we know this atom exists in this residue!
        int index = getAtomIndexByName(name);
        System.arraycopy(atomCoords, 0, coords, 3*index, 3);
    }
    
    
    public double distanceTo(Residue res2){
        //distance between the residues (measured at the closest atoms)
        double minDist = Double.POSITIVE_INFINITY;//minimum distance
        
        //for both residues, coordinates are concatenated 3-vectors
        for(int atNum1=0; atNum1<atoms.size(); atNum1++){
            for(int atNum2=0; atNum2<res2.atoms.size(); atNum2++){
                
                double dist = VectorAlgebra.distance(coords, atNum1, res2.coords, atNum2);
                minDist = Math.min(minDist,dist);
            }
        }
        
        return minDist;
    }
    
    
    
    
    public double[][] atomDistanceMatrix(){
        //matrix of distances between all atoms
        int numAtoms = atoms.size();
        double[][] ans = new double[numAtoms][numAtoms];
        
        for(int atNum1=0; atNum1<numAtoms; atNum1++){
            for(int atNum2=0; atNum2<numAtoms; atNum2++){
                
                double dist = VectorAlgebra.distance(coords,atNum1,coords,atNum2);
                ans[atNum1][atNum2] = dist;
            }
        }
        
        return ans;
    }
    
    public double[][] estBondDistanceMatrix(){
        //similar to atomDistanceMatrix, but only has nonzero values for bonded atoms,
        //and these are estimated based on forcefield equilibrium bond lengths
        //useful when we don't have coordinate available (e.g. in some template assignments)
        int numAtoms = atoms.size();
        double[][] ans = new double[numAtoms][numAtoms];
        
        ForcefieldParams ffParams = EnvironmentVars.resTemplates.ffParams;
        
        for(int atNum1=0; atNum1<numAtoms; atNum1++){
            
            Atom atom1 = atoms.get(atNum1);
            int atomType1 = atom1.type;
            
            for(Atom atom2 : atom1.bonds){
                int atNum2 = getAtomIndexByName(atom2.name);
                int atomType2 = atom2.type;
                
                ans[atNum1][atNum2] = ffParams.getBondEBL(atomType1, atomType2);
                
                if(Double.isNaN(ans[atNum1][atNum2])){//No EBL for these atom types
                    //so estimate based on element types
                    ans[atNum1][atNum2] = ffParams.estBondEBL(atomType1, atomType2);
                }
            }
        }
        
        return ans;
    }
    
    
    
    public void removeInterResBonds(){
        //disconnect all bonds between this and other residues
        interResBondsMarked = false;
        
        for(Atom atom : atoms){
            
            for(int bondNum=atom.bonds.size()-1; bondNum>=0; bondNum--){
                //iterate backwards so can remove bonds without messing up the index
                
                Atom bondedAtom = atom.bonds.get(bondNum);
                
                if(bondedAtom.res!=this){//inter-residue bond
                    
                    atom.bonds.remove(bondedAtom);
                    bondedAtom.bonds.remove(atom);
                }
            }
        }
    }
    
}
