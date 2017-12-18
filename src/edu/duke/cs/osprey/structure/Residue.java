/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResTemplateMatching;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.Protractor;
import edu.duke.cs.osprey.tools.StringParsing;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.BitSet;
import java.util.List;

/**
 *
 * @author mhall44
 */
public class Residue implements Serializable {
    //This is a residue
    //in the general sense: it's a piece of a molecular system that can mutate, take on RCs, etc.
    //and is the basic unit for energy calculations: energy is broken down to intra-residue energies,
    //pairwise interactions between residues, and maybe a few triples etc.
    
    private static final long serialVersionUID = -6188262155881479376L;
    
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
    
    
    public static enum SecondaryStructure {
        HELIX,
        SHEET,
        LOOP;
    }
    public SecondaryStructure secondaryStruct = SecondaryStructure.LOOP;
    
    
    public Residue(Residue other) {
        this(copyAtoms(other.atoms), Arrays.copyOf(other.coords, other.coords.length), other.fullName, other.molec);
        
        this.indexInMolecule = other.indexInMolecule;
        this.template = other.template;
        this.confProblems = new ArrayList<>(other.confProblems);
        this.pucker = other.pucker;
        this.secondaryStruct = other.secondaryStruct;
        
        // init intra-res atom bonds
        if (this.template != null) {
            markIntraResBondsByTemplate();
        }
        
        // NOTE: we don't copy inter-res atom bonds, so the molecule will have to re-bond everything
    }
    
    public static ArrayList<Atom> copyAtoms(ArrayList<Atom> atoms) {
        ArrayList<Atom> out = new ArrayList<>();
        for (Atom atom : atoms) {
            out.add(atom.copy());
        }
        return out;
    }
    
    public Residue(ArrayList<Atom> atoms, ArrayList<double[]> coords, String fullName, Molecule molec) {
        this(atoms, convertCoords(coords), fullName, molec);
    }
    
    private static double[] convertCoords(ArrayList<double[]> coords) {
        
        //record coordinates, unless coordList is null (then leave coords null too)
        if (coords == null) {
            return null;
        }
        
        int numAtoms = coords.size();
        double[] out = new double[3*numAtoms];
        for (int i=0; i<numAtoms; i++) {
            System.arraycopy(coords.get(i), 0, out, 3*i, 3);
        }
        return out;
    }
    
    public Residue(ArrayList<Atom> atoms, double[] coords, String fullName, Molecule molec) {
        //generate a residue.  Don't put in template and related stuff like bonds (that can be assigned later)
        
        this.atoms = atoms;
        this.fullName = fullName;
        this.molec = molec;
        
        int numAtoms = atoms.size();
        this.coords = coords;
        if (coords != null) {
            if(numAtoms*3 != coords.length){
                throw new RuntimeException("ERROR: Trying to instantiate residue with "+numAtoms+
                        " atoms but "+coords.length+"/3 atom coordinates");
            }
        }
        for(int a=0; a<numAtoms; a++){
            atoms.get(a).res = this;
            atoms.get(a).indexInRes = a;
        }
    }
    
    // cache res numbers for performance
    // they're used a LOT and profiling shows this is actually a performance bottleneck!
    private String resNum = null;
    
    public String getPDBResNumber() {
        
        // populate the cache if needed
        if (resNum == null) {
            if (fullName.length() > 5) {
            		//let's include the chain ID in case there are multiple chains
            		//with overlapping numbers
            		resNum = StringParsing.getToken(fullName.substring(3,5),1)
            					+ StringParsing.getToken(fullName.substring(5),1);
            } else {
                resNum = Integer.toString(indexInMolecule+1);
            }
        }
        
        return resNum;
    }
    
    public char getChainId() {
        return fullName.charAt(4);
    }
    
    public boolean assignTemplate(ResidueTemplateLibrary templateLib) {
        //assign a template to this residue if possible, using the ResidueTemplateLibrary
        //return if successful or not
        
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
        
        return assignTemplate(templCandidates, templateLib.ffparams);
    }
            
    public boolean assignTemplate(List<ResidueTemplate> templCandidates, ForcefieldParams ffParams) {    
        //now try to match up atoms
        ResTemplateMatching bestMatching = null;
        double bestScore = Double.POSITIVE_INFINITY;
        
        for(ResidueTemplate templ : templCandidates){
            ResTemplateMatching templMatching = new ResTemplateMatching(this, templ, ffParams);
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
    
    
    
    public void markIntraResBondsByTemplate(){
        //assign all the bonds between atoms in this residue, based on the template

        
        int numAtoms = atoms.size();
        ArrayList<Atom> templateAtoms = template.templateRes.atoms;
        
        if(templateAtoms.size()!=numAtoms)
            throw new RuntimeException("ERROR: Template for "+fullName+" has the wrong number of atoms");
        
        //Within-residue bonds are exactly as in the template
        copyIntraBondsFrom(template.templateRes);
        intraResBondsMarked = true;
    }

    public void reconnectInterResBonds(){
        template.interResBonding.connectInterResBonds(this, false);
        interResBondsMarked = true;
    }
    
    
    public boolean[][] getIntraResBondMatrix(){
        //matrix(i,j) indicates whether atoms i and j are bonded to each other

        int numAtoms = atoms.size();
        boolean intraResBondMatrix[][] = new boolean[numAtoms][numAtoms];
        for(int atNum1=0; atNum1<numAtoms; atNum1++){
            for(int atNum2=0; atNum2<numAtoms; atNum2++){
                
                Atom atom1 = atoms.get(atNum1);
                Atom atom2 = atoms.get(atNum2);
                
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
    
    public void copyIntraBondsFrom(Residue other) {
        
        boolean[][] isBonded = other.getIntraResBondMatrix();
        int numAtoms = atoms.size();
        
        // double check the residues match atoms
        for (int i=0; i<numAtoms; i++) {
            if (!atoms.get(i).name.equalsIgnoreCase(other.atoms.get(i).name)) {
                throw new Error("ERROR: Atom names don't match aross residues, can't copy bonds");
            }
        }

        // copy the bonds
        for(int i=0; i<numAtoms; i++){
            Atom atomi = atoms.get(i);
            for(int j=0; j<numAtoms; j++){
                if(isBonded[i][j]){
                    Atom atomj = atoms.get(j);
                    atomi.bonds.add(atomj);
                }
            }
        }
        
        checkBonds();
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
    
    public Atom getAtomByName(String name) {
        int index = getAtomIndexByName(name);
        if (index >= 0) {
            return atoms.get(index);
        }
        return null;
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

    public boolean isBondedTo(Residue res2){
        for(Atom at : atoms){
            for(Atom at2 : res2.atoms){
                if(at.bonds.contains(at2))
                    return true;
            }
        }
        return false;
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
    
    public double[][] estBondDistanceMatrix(ForcefieldParams ffParams) {
        //similar to atomDistanceMatrix, but only has nonzero values for bonded atoms,
        //and these are estimated based on forcefield equilibrium bond lengths
        //useful when we don't have coordinate available (e.g. in some template assignments)
        int numAtoms = atoms.size();
        double[][] ans = new double[numAtoms][numAtoms];
        
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
    
    public int getNumDihedrals() {
        return template.numDihedrals;
    }
    
    public double getDihedralAngle(int i) {
        return Protractor.measureDihedral(coords, template.getDihedralDefiningAtoms(i));
    }
    
    public double[] getDihedralAngles() {
        double[] dihedrals = new double[getNumDihedrals()];
        for(int i=0; i<dihedrals.length; i++) {
            dihedrals[i] = getDihedralAngle(i);
        }
        return dihedrals;
    }
    
    
    public Residue equivalentInMolec(Molecule mol){
        //get the equivalent of this residue in molecule mol
        return mol.getResByPDBResNumber(getPDBResNumber());
    }
    
    public static ArrayList<Residue> equivalentInMolec(ArrayList<Residue> resList, Molecule mol){
        //list the equivalent in mol of each residue in resList
        ArrayList<Residue> ans = new ArrayList<>();
        for(Residue res : resList)
            ans.add(res.equivalentInMolec(mol));
        return ans;
    }
    
    public void deleteHydrogens(){
        //remove any hydrogens
        BitSet isHeavy = new BitSet();
        ArrayList<Atom> newAtoms = new ArrayList<>();
        for(int at=0; at<atoms.size(); at++){
            if(!atoms.get(at).isHydrogen()){
                isHeavy.set(at);
                newAtoms.add(atoms.get(at));
            }
        }
        
        double newCoords[] = new double[3*isHeavy.cardinality()];
        int newAtCounter = 0;
        for(int at=0; at<atoms.size(); at++){
            if(isHeavy.get(at)){
                System.arraycopy(coords, 3*at, newCoords, 3*newAtCounter, 3);
                newAtCounter++;
            }
        }
        
        atoms = newAtoms;
        coords = newCoords;
    }

    public void addAtom(Atom atom, double newAtomCoords[]){
        atoms.add(atom);
        double[] newCoords = new double[coords.length+3];
        System.arraycopy(coords, 0, newCoords, 0, coords.length);
        System.arraycopy(newAtomCoords, 0, newCoords, coords.length, 3);
        coords = newCoords;
        atom.res = this;
        atom.indexInRes = atoms.size()-1;
    }
}
